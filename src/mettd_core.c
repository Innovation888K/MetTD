/* MetTD core kernel: enzyme scoring + matrix-permutation test, reimplemented from MetTD.R
 *
 * Modes:
 *   --zscore var|sd      standardization denominator (var = legacy R behaviour, sd = revised)
 *   --level meandiff|fc  which level column of ref.csv to use (both exported from R)
 *   --pcorr old|add1     old = max(cnt>obs,1)/B (legacy R); add1 = (cnt>=obs+1)/(B+1) finite-permutation correction
 *   --gamma G            distance exponent in weight w = 1/D^G (legacy 2)
 *   --radius R           0 = unlimited; else drop D > R
 *   --perms B --seed S
 * Null model replicates shuffle_matrix(): the values occupying the (ref x enzyme) cells
 * are a uniform subset of ALL N*N matrix entries without replacement.
 *
 * dist binary layout: column-major doubles (R native), element (i,j) at offset i + j*N, 0-based.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static char *dist_path, *enzyme_path, *ref_path, *out_path;
static long N = 0;
static int zmode = 0;        /* 0=var, 1=sd */
static int lmode = 0;        /* 0=meandiff, 1=fc */
static int pmode = 1;        /* 0=old, 1=add1 */
static double gamma_exp = 2.0;
static double radius = 0.0;
static long perms = 1000;
static uint64_t seed = 42;
static int obs_only = 0;
static int formula = 0;         /* 0 = (Σw·z/Σw)^2 [legacy mean^2]; 1 = Σ z_i^2 / D^2  (Σ w·z_i^2) */

/* ---- xoshiro256++ ---- */
static uint64_t s_rng[4];
static uint64_t rotl(const uint64_t x, int k){ return (x<<k)|(x>>(64-k)); }
static uint64_t next_u64(void){
  const uint64_t result = rotl(s_rng[0]+s_rng[3],23)+s_rng[0];
  const uint64_t t = s_rng[1]<<17;
  s_rng[2]^=s_rng[0]; s_rng[3]^=s_rng[1]; s_rng[1]^=s_rng[2]; s_rng[0]^=s_rng[3];
  s_rng[2]^=t; s_rng[3]=rotl(s_rng[3],45);
  return result;
}
static void rng_init(uint64_t sd){
  uint64_t z = sd + 0x9E3779B97F4A7C15ULL;
  for(int i=0;i<4;i++){ z += 0x9E3779B97F4A7C15ULL; uint64_t r = z;
    r = (r ^ (r>>30)) * 0xBF58476D1CE4E5B9ULL; r = (r ^ (r>>27)) * 0x94D049BB133111EBULL;
    s_rng[i] = r ^ (r>>31); }
  for(int i=0;i<64;i++) next_u64();
}

typedef struct { long idx; char name[64]; } EnzRow;
typedef struct { long idx; double lv_md, lv_fc; } RefRow;

static EnzRow *enz; static long n_enz = 0;
static RefRow *ref; static long n_ref = 0;
static double *dist = NULL;      /* N*N */
static uint32_t *pool = NULL;    /* N*N */

static double *zlev = NULL;
static double *obs_score = NULL;

static void die(const char *m){ fprintf(stderr,"ERROR: %s\n",m); exit(1); }

static void read_dist(void){
  size_t total = (size_t)N*(size_t)N;
  dist = malloc(total*sizeof(double));
  if(!dist) die("alloc dist");
  FILE *f = fopen(dist_path,"rb"); if(!f) die("open dist bin");
  size_t got = 0, chunk = 1<<22;
  char *p = (char*)dist;
  size_t need = total*sizeof(double);
  while(got < need){
    size_t want = need - got;
    if(chunk < want) want = chunk;      /* never request more than remains: some C
                                          runtimes return 0 rather than a short count */
    size_t r = fread(p+got,1,want,f);
    if(r==0) die("short read dist bin");
    got += r;
  }
  fclose(f);
}

static void read_enzyme(void){
  FILE *f=fopen(enzyme_path,"r"); if(!f) die("open enzyme.csv");
  char line[512]; long cap=2048; n_enz=0;
  enz=malloc(cap*sizeof(EnzRow));
  if(!fgets(line,sizeof line,f)) die("enzyme header");
  while(fgets(line,sizeof line,f)){
    if(n_enz>=cap){cap*=2;enz=realloc(enz,cap*sizeof(EnzRow));}
    char name[64];
    if(sscanf(line,"%ld,%63[^\n]",&enz[n_enz].idx,name)!=2){
      if(sscanf(line,"%ld,%63s",&enz[n_enz].idx,name)!=2) continue;
    }
    strncpy(enz[n_enz].name,name,63); enz[n_enz].name[63]=0;
    n_enz++;
  }
  fclose(f);
}

static void read_ref(void){
  FILE *f=fopen(ref_path,"r"); if(!f) die("open ref.csv");
  char line[512]; long cap=4096; n_ref=0; long line_no=1;
  ref=malloc(cap*sizeof(RefRow));
  if(!fgets(line,sizeof line,f)) die("ref header");
  while(fgets(line,sizeof line,f)){
    line_no++;
    if(line[0]=='\n'||line[0]=='\r') continue;
    long idx; double a,b;
    if(sscanf(line,"%ld,%lf,%lf",&idx,&a,&b)!=3)
      die("ref.csv: unparsable row (node index must be an integer, three columns "
          "idx,level_meandiff,level_fc)");
    char extra;
    if(sscanf(line,"%ld,%lf,%lf %c",&idx,&a,&b,&extra)==4)
      die("ref.csv: trailing characters after the third column");
    if(n_ref>=cap){cap*=2;ref=realloc(ref,cap*sizeof(RefRow));}
    ref[n_ref].idx=idx; ref[n_ref].lv_md=a; ref[n_ref].lv_fc=b;
    n_ref++;
  }
  fclose(f);
}

static void standardize(void){
  zlev=malloc(n_ref*sizeof(double));
  double sum=0; for(long i=0;i<n_ref;i++) sum += lmode? ref[i].lv_fc : ref[i].lv_md;
  double mean=sum/n_ref;
  double ss=0; for(long i=0;i<n_ref;i++){double d=(lmode?ref[i].lv_fc:ref[i].lv_md)-mean; ss+=d*d;}
  double var=ss/(n_ref-1);
  double den = zmode ? sqrt(var) : var;   /* zmode 1 = sd */
  for(long i=0;i<n_ref;i++) zlev[i]=((lmode?ref[i].lv_fc:ref[i].lv_md)-mean)/den;
}

/* score one enzyme given a function to fetch D for (ref r, enzyme e) */
static double score_enzyme_from_array(const double *dcol /* n_ref values already in ref order */){
  double num=0, wsum=0;
  if(formula==1){
    /* S = Σ z_i^2 / D_i^2 = Σ w * z_i^2  (no normalization, no squaring the sum) */
    for(long r=0;r<n_ref;r++){
      double D=dcol[r];
      if(D==0.0) continue;
      if(radius>0.0 && D>radius) continue;
      double w = pow(D, -gamma_exp);
      if(!(w>0)) continue;
      num += w*zlev[r]*zlev[r];
    }
    if(num>0) return num;
    return -10000.0;
  }
  for(long r=0;r<n_ref;r++){
    double D=dcol[r];
    if(D==0.0) continue;                       /* legacy: drop D==0 */
    if(radius>0.0 && D>radius) continue;       /* optional radius cut */
    double w = pow(D, -gamma_exp);             /* Inf -> 0 */
    if(!(w>0)) continue;                       /* Inf or NaN distance contribute nothing */
    wsum += w; num += w*zlev[r];
  }
  if(wsum>0) return (num/wsum)*(num/wsum);
  return -10000.0;
}

static double *dcol_buf=NULL;

static void compute_observed(void){
  obs_score=malloc(n_enz*sizeof(double));
  dcol_buf=malloc(n_ref*sizeof(double));
  for(long e=0;e<n_enz;e++){
    long c = enz[e].idx-1;                     /* 1-based -> 0-based column */
    for(long r=0;r<n_ref;r++) dcol_buf[r]=dist[(ref[r].idx-1)+(size_t)c*N];
    obs_score[e]=score_enzyme_from_array(dcol_buf);
  }
}

static void permute(void){
  size_t total=(size_t)N*(size_t)N;
  pool=malloc(total*sizeof(uint32_t));
  if(!pool) die("alloc pool");
  for(size_t i=0;i<total;i++) pool[i]=(uint32_t)i;
  long K=n_ref*n_enz;
  if((size_t)K>total) die("K exceeds pool");
  long *cnt=malloc(n_enz*sizeof(long));       /* counts (old: >, add1: >=) */
  memset(cnt,0,n_enz*sizeof(long));

  double *num=malloc(n_enz*sizeof(double));
  double *wsum=malloc(n_enz*sizeof(double));
  long *nvalid=malloc(n_enz*sizeof(long));

  for(long p=0;p<perms;p++){
    memset(num,0,n_enz*sizeof(double));
    memset(wsum,0,n_enz*sizeof(double));
    memset(nvalid,0,n_enz*sizeof(long));
    /* partial Fisher-Yates over pool, enzyme-major placement */
    long j;
    for(j=0;j<K;j++){
      uint32_t r = (uint32_t)(j + (next_u64() % (uint64_t)(total-(size_t)j)));
      uint32_t tmp=pool[j]; pool[j]=pool[r]; pool[r]=tmp;
      double D = dist[pool[j]];
      long e = j / n_ref, rr = j % n_ref;
      if(D==0.0) continue;
      if(radius>0.0 && D>radius) continue;
      double w = pow(D,-gamma_exp);
      if(!(w>0)) continue;
      wsum[e]+=w;
      num[e]+= (formula==1) ? (w*zlev[rr]*zlev[rr]) : (w*zlev[rr]);
      nvalid[e]++;
    }
    for(long e=0;e<n_enz;e++){
      double sc;
      if(formula==1) sc = (nvalid[e]>0)? num[e] : -10000.0;
      else           sc = (wsum[e]>0)? (num[e]/wsum[e])*(num[e]/wsum[e]) : -10000.0;
      int hit = pmode ? (sc >= obs_score[e]) : (sc > obs_score[e]);
      if(hit) cnt[e]++;
    }
  }

  /* p-values */
  double *pv=malloc(n_enz*sizeof(double));
  for(long e=0;e<n_enz;e++)
    pv[e] = pmode ? ((double)cnt[e]+1.0)/((double)perms+1.0)
                  : (double)(cnt[e]>0?cnt[e]:1)/(double)perms;

  /* BH q-values */
  long *ord=malloc(n_enz*sizeof(long));
  for(long i=0;i<n_enz;i++) ord[i]=i;
  /* simple insertion-free sort: use qsort */
  int cmp(const void *a,const void *b){
    double pa=pv[*(const long*)a], pb=pv[*(const long*)b];
    return pa<pb?-1:(pa>pb?1:0);
  }
  qsort(ord,n_enz,sizeof(long),cmp);
  double *q=malloc(n_enz*sizeof(double));
  double prev=1.0;
  for(long k=n_enz-1;k>=0;k--){
    long i=ord[k]; double r=(double)n_enz/(double)(k+1);
    double val=pv[i]*r; if(val>prev) val=prev; prev=val; q[i]=val;
  }

  FILE *f=fopen(out_path,"w");
  fprintf(f,"enzyme,score,p,q\n");
  for(long e=0;e<n_enz;e++)
    fprintf(f,"%s,%.10g,%.10g,%.10g\n",enz[e].name,obs_score[e],pv[e],q[e]);
  fclose(f);
}

int main(int argc,char**argv){
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"--dist")) dist_path=argv[++i];
    else if(!strcmp(argv[i],"--enzyme")) enzyme_path=argv[++i];
    else if(!strcmp(argv[i],"--ref")) ref_path=argv[++i];
    else if(!strcmp(argv[i],"--out")) out_path=argv[++i];
    else if(!strcmp(argv[i],"--n")) N=atol(argv[++i]);
    else if(!strcmp(argv[i],"--zscore")) zmode=!strcmp(argv[++i],"sd");
    else if(!strcmp(argv[i],"--level")) lmode=!strcmp(argv[++i],"fc");
    else if(!strcmp(argv[i],"--pcorr")) pmode=!strcmp(argv[++i],"add1");
    else if(!strcmp(argv[i],"--gamma")) gamma_exp=atof(argv[++i]);
    else if(!strcmp(argv[i],"--formula")) formula=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--radius")) radius=atof(argv[++i]);
    else if(!strcmp(argv[i],"--perms")) perms=atol(argv[++i]);
    else if(!strcmp(argv[i],"--seed")) seed=strtoull(argv[++i],NULL,10);
    else if(!strcmp(argv[i],"--obs-only")) obs_only=1;
    else { fprintf(stderr,"unknown arg %s\n",argv[i]); return 1; }
  }
  if(!dist_path||!enzyme_path||!ref_path||!out_path||N<=0) die("missing args");
  fprintf(stderr,"loading dist (%ld x %ld)...\n",N,N);
  read_dist(); read_enzyme(); read_ref();
  fprintf(stderr,"n_enz=%ld n_ref=%ld\n",n_enz,n_ref);
  standardize(); compute_observed();
  fprintf(stderr,"observed scores done\n");
  if(obs_only){
    FILE *f=fopen(out_path,"w");
    fprintf(f,"enzyme,score\n");
    for(long e=0;e<n_enz;e++) fprintf(f,"%s,%.10g\n",enz[e].name,obs_score[e]);
    fclose(f); return 0;
  }
  permute();
  fprintf(stderr,"done -> %s\n",out_path);
  return 0;
}
