#!/usr/bin/env python3
"""Fill time/dose and the six new studies' descriptors from isa_cache into the v3 table."""
import sys, os, csv, json, glob, re
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd

CACHE = 'D:/MetTD_recompute/isa_cache'
TAB = r'D:/ac结果/第一次返修/返修新材料/基准数据集元数据表_v3_100条.csv'
D = pd.read_csv(TAB, encoding='utf-8-sig', dtype=str).fillna('')
NEW6 = {'MTBLS1967', 'MTBLS1974', 'MTBLS2085', 'MTBLS11858', 'MTBLS3129', 'MTBLS8062'}

# ---------- MTBLS ISA ----------
def read_tsv(path):
    raw = open(path, encoding='utf-8', errors='replace').read().splitlines()
    rows = [l.split('\t') for l in raw if l.strip()]
    return rows

inv_info, sample_factors = {}, {}
for f in sorted(glob.glob(f'{CACHE}/*_i_Investigation.txt')):
    acc = os.path.basename(f).split('_')[0]
    txt = open(f, encoding='utf-8', errors='replace').read()
    def grab(label):
        def ok(v):
            v = v.strip()
            if not v or '\t' in v:
                return ''
            if re.match(r'^(PMID|DOI|Study|Investigation|Comment|Protocol|Name)', v):
                return ''
            return v if len(v) < 200 else ''
        m = re.search(rf'^{label}\t"?([^"\r\n]*)"?[ \t]*$', txt, re.M)
        if m and ok(m.group(1)):
            return m.group(1).strip()
        m = re.search(rf'^{label}\t?[ \t]*\r?\n\t?"?([^"\r\n\t]+)"?', txt, re.M)
        return ok(m.group(1)) if m else ''
    inv_info[acc] = dict(title=grab('Study Title'), doi=grab('Study Publication DOI'),
                         pmid=grab('Study PubMed ID'))
    inv_info[acc]['pubtitle'] = grab('Study Publication Title')

for f in sorted(glob.glob(f'{CACHE}/*_s_*.txt')):
    acc = os.path.basename(f).split('_')[0]
    rows = read_tsv(f)
    hdr_i = None
    for i, r in enumerate(rows[:5]):
        if any('Sample Name' in c for c in r):
            hdr_i = i; break
    if hdr_i is None:
        continue
    hdr = rows[hdr_i]
    idx = {c: j for j, c in enumerate(hdr)}
    sname = [c for c in hdr if c.strip() == 'Sample Name']
    fcols = [c for c in hdr if c.startswith('Factor Value') or c.startswith('Characteristics[')]
    if not sname:
        continue
    m = {}
    for r in rows[hdr_i + 1:]:
        if len(r) < len(hdr):
            continue
        key = r[idx[sname[0]]].strip()
        if key:
            m.setdefault(key, {c: r[idx[c]].strip() for c in fcols})
    sample_factors[acc] = m

# ---------- Metabolomics Workbench ----------
mw = {}
for f in sorted(glob.glob(f'{CACHE}/*_mw.json')):
    acc = os.path.basename(f).split('_')[0]
    try:
        d = json.loads(open(f, encoding='utf-8', errors='replace').read())
        mw[acc] = d
    except Exception:
        pass

# organism / tissue from sample-table Characteristics; instrument from a_ assay files
def col_mode(sf, col):
    vals = [v.get(col, '') for v in sf.values() if v.get(col, '').strip()]
    if not vals:
        return ''
    return max(set(vals), key=vals.count)

org_tissue, instruments = {}, {}
for acc, sf in sample_factors.items():
    org_tissue[acc] = (col_mode(sf, 'Characteristics[Organism]'),
                       col_mode(sf, 'Characteristics[Organism part]') or
                       col_mode(sf, 'Characteristics[Tissue]') or
                       col_mode(sf, 'Characteristics[tissue type]'))
for f in sorted(glob.glob(f'{CACHE}/*_a_*.txt')):
    acc = os.path.basename(f).split('_')[0]
    rows = read_tsv(f)
    for r in rows[:5]:
        if any('Instrument' in c for c in r):
            hdr = r
            icol = [c for c in hdr if 'Instrument' in c]
            if not icol:
                continue
            j = hdr.index(icol[0])
            vals = [x[j].strip() for x in rows[rows.index(r) + 1:] if len(x) > j and x[j].strip()]
            if vals:
                instruments.setdefault(acc, set()).add(max(set(vals), key=vals.count))
            break

print('ISA investigation files:', sorted(inv_info))
print('ISA sample tables:', {k: len(v) for k, v in sorted(sample_factors.items())})
print('MW JSON:', sorted(mw))
print()
for acc, info in sorted(inv_info.items()):
    print(f'{acc:12s} {info["title"][:70]}')
    if info['pmid'] or info['doi']:
        print(f'             PMID {info["pmid"]} DOI {info["doi"]}')

# ---------- fill ----------
fixed = {'organism': {}, 'tissue': {}, 'matrix': {}, 'platform': {},
         'study_title': {}, 'study_reference': {}, 'time': {}, 'dose': {}}
for acc, info in inv_info.items():
    fixed['study_title'][acc] = info['title']
    fixed['study_reference'][acc] = (f'PMID {info["pmid"]}' if info['pmid'] else '') or \
                                    (f'DOI {info["doi"]}' if info['doi'] else '')
for acc, d in mw.items():
    t = (d.get('study_title') or '').strip()
    if t:
        fixed['study_title'][acc] = t
    if d.get('pubmed_id'):
        fixed['study_reference'][acc] = f'PMID {d.get("pubmed_id")}'
    elif d.get('doi'):
        fixed['study_reference'][acc] = f'DOI {d.get("doi")}'
    fixed['organism'][acc] = (d.get('organism') or d.get('species') or '')
    fixed['tissue'][acc] = ''
    fixed['platform'][acc] = (d.get('instrumentType') or d.get('instrument_type') or
                              d.get('analysis_type') or '')

def guess_matrix(factors):
    s = ' '.join(factors.values()).lower()
    for k, v in [('plasma', 'plasma'), ('serum', 'serum'), ('urine', 'urine'),
                 ('tissue', 'tissue'), ('brain', 'brain'), ('kidney', 'kidney'),
                 ('liver', 'liver'), ('fibroblast', 'fibroblasts'), ('cell', 'cells'),
                 ('pbmc', 'PBMC'), ('blood', 'whole blood')]:
        if k in s:
            return v
    return ''

def factor_td(factors):
    t = d = ''
    for k, v in factors.items():
        kl = k.lower()
        if 'time' in kl and not t and v.strip():
            t = v.strip()
        if ('concentration' in kl or 'dose' in kl) and not d and v.strip():
            d = v.strip()
    return t, d

def match_sample(sf, n):
    if n in sf:
        return sf[n]
    for k, v in sf.items():
        if k.startswith(n) or n.startswith(k) or k.rstrip('0123456789') == n.rstrip('0123456789'):
            return v
    return None

D['match_rate_s_table'] = ''
for i, r in D.iterrows():
    acc = r.accession
    if acc in NEW6:
        org, tis = org_tissue.get(acc, ('', ''))
        plat = ';'.join(sorted(instruments.get(acc, set())))[:80]
        vals = {'organism': org, 'tissue': tis, 'platform': plat,
                'study_title': fixed['study_title'].get(acc, ''),
                'study_reference': fixed['study_reference'].get(acc, '')}
        for col, v in vals.items():
            D.at[i, col] = v        # always overwrite: clears stale/MISSING cells
        if not D.at[i, 'matrix'].strip():
            sf = sample_factors.get(acc)
            if sf:
                D.at[i, 'matrix'] = guess_matrix(list(sf.values())[0])
    # time/dose from sample factors
    sf = sample_factors.get(acc)
    if sf:
        names = [x.strip() for x in f'{r.control_group};{r.model_group}'.split(';') if x.strip()]
        tvals, dvals, hit = [], [], 0
        for n in names:
            fv = match_sample(sf, n)
            if fv is None:
                continue
            hit += 1
            t, dd = factor_td(fv)
            if t and not re.match(r'^\d{1,2}:\d{2}$', t):   # drop clock-style run times
                tvals.append(t)
            if dd:
                dvals.append(dd)
        if tvals and any(re.search(r'\d', v) for v in tvals):
            D.at[i, 'time_point'] = ';'.join(sorted(set(tvals)))
        if dvals and any(re.search(r'\d', v) for v in dvals):
            D.at[i, 'dose'] = ';'.join(sorted(set(dvals)))
        D.at[i, 'match_rate_s_table'] = f'{hit}/{len(names)}'

D.to_csv(TAB, index=False, encoding='utf-8-sig')
print('\nrewrote', TAB)
n6 = D[D.accession.isin(NEW6)]
print('\nsix new studies now:')
print(n6[['accession', 'study_title', 'organism', 'platform', 'study_reference']]
      .drop_duplicates().to_string(index=False, max_colwidth=46))
print('\nrecords with a real time/dose value:',
      int((~D.time_point.str.startswith('not applicable')).sum()), '/100')
