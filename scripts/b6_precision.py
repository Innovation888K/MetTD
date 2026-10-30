"""B6: ranking-quality metrics for the 13-drug benchmark (MetTD vs SEA vs SwissTargetPrediction).
Outputs first-hit rank, hits@N, precision@N, hypergeometric p, enrichment vs random.
"""
import csv, sys, os, zipfile, re
from math import comb
from xml.etree import ElementTree as ET

os.chdir(r'D:/ac结果/fig3干实验验证')
sys.stdout.reconfigure(encoding='utf-8')

DRUGS = ['5-Fluorouracil','Bleomycin','Bortezomib','Ceftriaxone','Enasidenib','Pemetrexed',
         'Pimozide','acyclovir','cannabidiol','dihydroartemisinin','metformin',
         'phosphoric acid','rosiglitazone']
NS = [5,10,20,30,50]

def read_xlsx_rows(path):
    z=zipfile.ZipFile(path); ss=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')):
            ss.append(''.join(t.text or '' for t in si.iter() if t.tag.endswith('}t')))
    sh=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    ns='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    rows=[]
    for row in sh.iter(ns+'row'):
        vals=[]
        for c in row:
            v=c.find(ns+'v')
            if v is None: vals.append(''); continue
            vals.append(ss[int(v.text)] if c.get('t')=='s' else v.text)
        rows.append(vals)
    return rows

# UniProt accession -> gene map
acc2gene = {}
for r in read_xlsx_rows(r'D:/RSpace/MetTD/uniprot_human.xlsx')[1:]:
    if len(r) > 3 and r[0]:
        acc2gene[r[0]] = r[3].split(' ')[0]

def targets(drug):
    p = f'targets/{drug}.csv'
    with open(p, encoding='utf-8-sig', errors='replace') as f:
        return {r['genename'].strip() for r in csv.DictReader(f) if r.get('genename')}

def mettd_ranked(drug):
    p = f'Topn_hitrate/MetTD/{drug}/result.xlsx'
    if not os.path.exists(p):
        p = f'D:/wyc/MetTD代谢组_TopN_hit_rate/{drug}/result.xlsx'
    rows = read_xlsx_rows(p)[1:]
    data = [(float(r[0]), r[2]) for r in rows if r and r[0] and r[2]]
    data.sort(key=lambda x: x[0])
    ecs = [e for _, e in data]
    # EC -> genes (EC2genename: rows in same enzyme order as the p-ordered significant list)
    gpath = f'Topn_hitrate/MetTD/EC2genename/{drug}.csv'
    gene_lists = []
    with open(gpath, encoding='utf-8-sig', errors='replace') as f:
        for r in csv.reader(f):
            if not r or r[0] in ('', 't.x.'): continue
            genes = [x for x in r[1:] if x and x != 'NA']
            gene_lists.append(genes)
    ranked = []
    for i, ec in enumerate(ecs):
        gl = gene_lists[i] if i < len(gene_lists) else []
        ranked.extend(gl)
    return ranked

def sea_ranked(drug):
    import io
    p = f'Topn_hitrate/sea/{drug}.xls'
    text = None
    try:
        with open(p, encoding='utf-8-sig', errors='replace') as f:
            text = f.read()
        if len(text.strip().splitlines()) < 3:
            text = None
    except FileNotFoundError:
        text = None
    if text is None:
        with zipfile.ZipFile(f'sea/{drug}.zip') as z:
            text = z.read('sea-results.xls').decode('utf-8', errors='replace')
    rows = []
    lines = [l for l in text.splitlines() if l.strip()]
    hdr = lines[0].split('	')
    ti = hdr.index('Target ID') if 'Target ID' in hdr else 1
    pi = hdr.index('P-Value') if 'P-Value' in hdr else 3
    for l in lines[1:]:
        parts = l.split('	')
        if len(parts) <= max(ti, pi): continue
        t = parts[ti].strip().split('_')[0]
        try: pv = float(parts[pi])
        except: pv = 1.0
        g = acc2gene.get(t, t)
        rows.append((pv, g))
    rows.sort(key=lambda x: x[0])
    seen=set(); out=[]
    for _, g in rows:
        if g and g not in seen: seen.add(g); out.append(g)
    return out

def swiss_ranked(drug):
    p = f'Topn_hitrate/swiss/{drug}.csv'
    rows = []
    with open(p, encoding='utf-8-sig', errors='ignore') as f:
        for r in csv.DictReader(f):
            acc = (r.get('Uniprot ID') or '').strip()
            g = acc2gene.get(acc) or (r.get('Common name') or '').strip()
            try: pr = float(r.get('Probability*', 0) or 0)
            except: pr = 0.0
            rows.append((pr, g))
    rows.sort(key=lambda x: -x[0])
    seen=set(); out=[]
    for _, g in rows:
        if g and g not in seen: seen.add(g); out.append(g)
    return out

def hyper_p(k, N, T, G):
    if G < N: N2 = G
    else: N2 = N
    if T > G: return 1.0
    return sum(comb(T,i)*comb(G-T,N2-i) for i in range(k, min(T,N2)+1)) / comb(G, N2)

enz_genes = set()
for r in csv.DictReader(open(r'D:/MetTD_recompute/out/EC_gene_map.csv', encoding='utf-8-sig')):
    enz_genes.add(r['gene_symbol'])

out = []
for drug in DRUGS:
    T = targets(drug)
    for method, fn in [('MetTD', mettd_ranked), ('SEA', sea_ranked), ('SwissTP', swiss_ranked)]:
        try:
            ranked = fn(drug)
        except Exception as e:
            print(f'{drug}/{method} ERR {e}'); continue
        G = len(ranked); hits_all = 0; first = None
        hset = set(T)
        for i, g in enumerate(ranked):
            if g in hset:
                hits_all += 1
                if first is None: first = i + 1
        row = dict(drug=drug, method=method, T=len(T), G=G, hits_total=hits_all, first_hit=first or 0)
        for N in NS:
            hitsN = sum(1 for g in ranked[:N] if g in hset)
            row[f'hits@{N}'] = hitsN
            row[f'prec@{N}'] = round(hitsN / N, 3)
            row[f'hyperp@{N}'] = f'{hyper_p(hitsN, N, len(hset), G):.2e}' if G else ''
        row['enrich@30'] = round((row['hits@30']/30) / (len(hset)/G), 2) if G else ''
        eset = hset & enz_genes
        row['T_enz'] = len(eset)
        hN = sum(1 for g in ranked[:30] if g in eset)
        row['hits@30_enz'] = hN
        row['prec@30_enz'] = round(hN/30, 3)
        fe = next((i+1 for i,g in enumerate(ranked) if g in eset), 0)
        row['first_hit_enz'] = fe
        out.append(row)

w = csv.writer(open(r'D:/MetTD_recompute/out/b6_precision.csv','w',newline='',encoding='utf-8-sig'))
hdr = ['drug','method','T','G','hits_total','first_hit'] + [f'hits@{n}' for n in NS] + [f'prec@{n}' for n in NS] + [f'hyperp@{n}' for n in NS] + ['enrich@30','T_enz','hits@30_enz','prec@30_enz','first_hit_enz']
w.writerow(hdr)
for r in out: w.writerow([r.get(h,'') for h in hdr])

import statistics
for m in ['MetTD','SEA','SwissTP']:
    rs = [r for r in out if r['method']==m]
    print(m, 'drugs:', len(rs),
          'median first-hit:', (statistics.median([r['first_hit'] for r in rs if r['first_hit']]) if any(r['first_hit'] for r in rs) else 'none'),
          'mean prec@30: %.3f' % statistics.mean([r['prec@30'] for r in rs]),
          'median enrich@30:', statistics.median([float(r['enrich@30'] or 0) for r in rs]),
          'hits@30 total:', sum(r['hits@30'] for r in rs))
    rs0 = rs
    fhe = [r['first_hit_enz'] for r in rs0 if r['first_hit_enz']]
    print(m, 'ENZYME-split: median first-hit-enzyme:', statistics.median(fhe) if fhe else 'none',
          ' mean prec@30_enz: %.3f' % statistics.mean([float(r['prec@30_enz'] or 0) for r in rs0]),
          ' hits@30_enz total:', sum(int(r['hits@30_enz']) for r in rs0))
