#!/usr/bin/env python3
"""Build the revised benchmark metadata table for benchmark_v3 (100 records).

New file only; the original 基准数据集元数据表_完成版.csv is not modified.
- PMM2-CDG label -> PMM2 (the deficient enzyme)
- adds evidence_tier and recodes perturbation_type to genetic/pharmacological/inferred
- restores full sample-group strings from benchmark_v3 (the old table truncated at 80 chars)
- fills time point and dose for every record (parsed from sample names, else 'not applicable')
- adds the six studies introduced in benchmark_v3
"""
import sys, re
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np, pandas as pd

RD = 'D:/MetTD_recompute'
META = r'D:/ac结果/第一次返修/返修新材料/基准数据集元数据表_完成版.csv'
OUT = r'D:/ac结果/第一次返修/返修新材料/基准数据集元数据表_v3_100条.csv'

v3 = pd.read_csv(f'{RD}/benchmark_v3.csv', encoding='utf-8-sig', dtype=str).fillna('')
ev = pd.read_csv(f'{RD}/benchmark_v3_mettd_eval.csv', encoding='utf-8-sig')
assert (ev.Program.values == v3.Program.values).all()
T = pd.read_csv(f'{RD}/benchmark_v3_groundtruth_tiers.csv', encoding='utf-8-sig')
tier_of = dict(zip(T.accession, T.tier))
m = pd.read_csv(META, encoding='utf-8-sig', dtype=str).fillna('')
m = m[m.accession != '']
first = m.groupby('accession').first()

g = pd.read_csv(f'{RD}/out/EC_gene_map.csv', dtype=str).fillna('')
gene2ec = {}
for _, r in g.iterrows():
    if r['gene_symbol']:
        gene2ec.setdefault(r['gene_symbol'].strip().upper(), set()).add(r['EC_number'].strip())

PT = {'T1 genetic': 'genetic', 'T2 pharmacological': 'pharmacological',
      'T3 inferred': 'inferred', 'T? new study': ''}
RATIONALE = {
 'MTBLS549': '58-h fasting; IDH1 assigned from the fasting-induced remodelling of TCA-cycle '
             'and NADP(+)-dependent metabolism reported in the source study',
 'ST000604': 'high-sugar diet in mice; IDH1 assigned from the reported shift in arginine/TCA '
             'metabolism, not from manipulation of IDH1',
 'MTBLS3128': 'rosuvastatin pre-dose serum profiling; IDH1 is not the drug target (HMGCR is); '
              'label assigned from the reported metabolic correlate',
 'ST002532': 'cisplatin-treated kidney tissue; IDH1 assigned from reported IDH1 involvement in '
             'cisplatin nephrotoxicity',
 'MTBLS3020': 'arginine challenge in Streptococcus mutans culture; IDH1 assigned by metabolic '
              'association',
 'ST000583': 'arginine/fluoride exposure of oral biofilms; IDH1 assigned by metabolic association',
 'ST002527': 'arginase inhibition in PDAC interstitial fluid; IDH1 assigned by metabolic association',
 'ST000302': 'study was designed around IDH1 in glioma (fluxomics); label is the study subject',
 'ST002221': 'ACLY-dependent glutaminolysis in ccRCC is the study subject; ACLY label direct',
 'ST002963': 'GSH/GSSG in KRAS-driven lung tumours; G6PD assigned as the NADPH source for '
             'glutathione regeneration',
 'ST001451-2': 'eleostearic acid treatment; ACSL1 assigned from acyl-CoA/lipid remodelling',
 'MTBLS1973': 'IDH1-mutant cells; target genetically altered',
 'ST000548': 'leukaemia-associated IDH1/IDH2 mutants; target genetically altered',
 'ST000276': 'IDH1 knockdown; target genetically altered',
 'ST000714': 'glioblastoma lines differing at IDH1; target genetically altered',
 'ST002231': 'IDH1-R132H versus EV cells with/without an IDH1 inhibitor; target genetically and '
             'pharmacologically altered',
 'ST000694': 'ACSL5 knockout versus floxed mice; target genetically altered',
 'ST001441': 'SUCLA2-deficient patient fibroblasts; target genetically altered',
 'ST002423': 'ATP7B (Wilson disease) model; target genetically altered',
 'ST002696': 'ATP10A-deficient mice; target genetically altered',
 'ST003347': 'ADSL-deficient myoblasts; target genetically altered',
 'ST002563': 'PMM2-CDG fibroblasts with AKR1B1 siRNA; label set to PMM2, the deficient enzyme',
 'MTBLS2159': 'macrophage ACLY deficiency; target genetically altered',
 'ST003423': 'ACLY/ACSS2-directed intervention in MASLD; label is the agent-directed enzyme',
 'MTBLS2130': 'drug-response phenotyping in colorectal cancer organoids; IDH1 label taken from '
              'the reported responder signature',
}
TIME_RE = re.compile(r'(\d+(?:\.\d+)?)\s*(hr|h|min|month|mo|day|d)(?![A-Za-z])', re.I)
DOSE_RE = re.compile(r'(\d+(?:\.\d+)?)\s*(nM|uM|µM|mM|mg/kg|ug/mL|mg)(?![A-Za-z])', re.I)

def parse(s, rx):
    out = set()
    for x in s.split(';'):
        mm = rx.search(x)
        if mm:
            out.add(f'{mm.group(1)}{mm.group(2)}'.lower())
    return ';'.join(sorted(out))

rows = []
for i, e in v3.iterrows():
    acc, tgt = e.Program, e.Target
    lab = 'PMM2' if tgt == 'PMM2-CDG' else tgt
    tier = tier_of.get(acc, 'T? new study')
    src = first.loc[acc] if acc in first.index else None
    ec = ';'.join(sorted(gene2ec.get(lab.upper(), set())))
    tp = parse(e.Ctrl_group, TIME_RE) or parse(e.Model_group, TIME_RE)
    dz = parse(e.Ctrl_group, DOSE_RE) or parse(e.Model_group, DOSE_RE)
    rows.append(dict(
        entry_id=e.ID, accession=acc, status=e.status, target_label=lab, target_EC=ec,
        evidence_tier=tier.replace('T1 ', 'T1 ').replace('T? new study', 'T? (new study, tier to confirm)'),
        perturbation_type=PT.get(tier, ''),
        label_rationale=RATIONALE.get(acc, ''),
        label_supporting_reference=(src.reference if src is not None else ''),
        control_group=e.Ctrl_group, model_group=e.Model_group,
        n_ctrl=len([x for x in e.Ctrl_group.split(';') if x]), n_model=len([x for x in e.Model_group.split(';') if x]),
        time_point=tp or 'not applicable (between-group design)',
        dose=dz or 'not applicable (between-group design)',
        organism=(src.organism if src is not None else 'MISSING'),
        tissue=(src.tissue if src is not None else 'MISSING'),
        matrix=(src.matrix_inferred if src is not None else 'MISSING'),
        platform=(src.platform if src is not None else 'MISSING'),
        study_title=(src.study_title if src is not None else 'MISSING'),
        study_reference=(src.reference if src is not None else 'MISSING'),
        mettd_target_rank=pd.to_numeric(ev.Target_rank.iloc[i], errors='coerce')))
D = pd.DataFrame(rows)
D.to_csv(OUT, index=False, encoding='utf-8-sig')

print('wrote', OUT, D.shape)
print('\nrows per tier:'); print(D.evidence_tier.value_counts().to_string())
print('\nrows per perturbation_type:'); print(D.perturbation_type.value_counts().to_string())
print('\nblank/MISSING cells per column:')
for c in D.columns:
    b = int((D[c].astype(str).str.strip().isin(['', 'MISSING', 'nan'])).sum())
    if b: print(f'  {c:28s} {b}')
print('\ntime_point filled with a real value:', int((~D.time_point.str.startswith('not applicable')).sum()), '/100')
print('dose filled with a real value:', int((~D.dose.str.startswith('not applicable')).sum()), '/100')
print('\nlabel_rationale missing for:'); print(sorted(set(D[D.label_rationale == ''].accession)))
