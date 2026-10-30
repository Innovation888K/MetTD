"""B2 finalize: full Fig4d/ST1 regen + drug q-count summary + v2 backfill text."""
import csv, os, subprocess, sys
sys.stdout.reconfigure(encoding='utf-8')
os.chdir(r'D:/MetTD_recompute')

subprocess.run([sys.executable, 'fig4d_trend.py'], check=True)

print('\n=== drug rerun q counts (corrected null) ===')
rows_out = []
for f in sorted(os.listdir('out/b2')):
    if f.startswith('drug_'):
        d = list(csv.DictReader(open(f'out/b2/{f}')))
        q = sum(1 for r in d if float(r['q']) < 0.05)
        p = sum(1 for r in d if float(r['p']) < 0.05)
        top = [r['enzyme'] for r in sorted(d, key=lambda r: -float(r['score']))[:3]]
        rows_out.append([f.replace('drug_', '').replace('.csv', ''), p, q, ';'.join(top)])
        print(f'{f}: p<0.05={p} q<0.05={q} top3={top}')
w = csv.writer(open(r'D:/ac结果/第一次返修/返修新材料/13药重跑_q计数_自产6药.csv', 'w', newline='', encoding='utf-8-sig'))
w.writerow(['drug', 'p<0.05', 'q<0.05', 'top3_by_score'])
w.writerows(rows_out)

print('\n=== v2 backfill lines ===')
tq = {}
for f in os.listdir('out/b2'):
    if f.startswith('t') and f.endswith('.csv'):
        d = list(csv.DictReader(open(f'out/b2/{f}')))
        tq[f.replace('.csv', '')] = sum(1 for r in d if float(r['q']) < 0.05)
met = {}
for sx in ['5-2', '5-3', '5-4']:
    f = 'out/fig3c_sd_meandiff_matrix.csv' if sx == '5-2' else f'out/b2/met{sx}.csv'
    if os.path.exists(f):
        d = list(csv.DictReader(open(f)))
        order = sorted(d, key=lambda r: -float(r['score']))
        idx = {r['enzyme']: i + 1 for i, r in enumerate(order)}
        met[sx] = (sum(1 for r in d if float(r['q']) < 0.05), idx.get('4.2.1.3'), idx.get('1.1.1.42'), idx.get('2.3.3.1'))
        print(f'met{sx}: q<0.05={met[sx][0]} ACO2 rank={met[sx][1]} IDH1 rank={met[sx][2]} CS rank={met[sx][3]}')
for k in sorted(tq, key=lambda x: (float(x[1:-1]), x[-1])):
    print(f'{k}: q<0.05 = {tq[k]}')
