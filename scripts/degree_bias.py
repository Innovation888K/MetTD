import csv, sys, glob, math
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
deg = {}
M = np.memmap('model_dist.bin', dtype=np.float64, mode='r').reshape(9492, 9492)
enz = {r['enzyme']: int(r['idx']) for r in csv.DictReader(open('enzyme.csv'))}
A = (M == 2.0); np.fill_diagonal(A, False); dcol = A.sum(axis=0)
for e, i in enz.items(): deg[e] = int(dcol[i - 1])
def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float); ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    return float((rx * ry).sum() / math.sqrt((rx ** 2).sum() * (ry ** 2).sum()))
rows_out = [['dataset', 'spearman_degree_score', 'top5_enzymes', 'top5_degree_percentile']]
for f in sorted(glob.glob('out/b2/*.csv')) + ['out/fig3c_sd_meandiff_matrix.csv']:
    rows = list(csv.DictReader(open(f)))
    sc = [float(r['score']) for r in rows]; dg = [deg.get(r['enzyme'], 0) for r in rows]
    top5 = sorted(rows, key=lambda r: -float(r['score']))[:5]
    tp = [round(100 * sum(1 for v in deg.values() if v <= deg.get(r['enzyme'], 0)) / len(deg)) for r in top5]
    name = f.replace('\\', '/').split('/')[-1].replace('.csv', '')
    rows_out.append([name, f'{spearman(dg, sc):+.3f}', ';'.join(r['enzyme'] for r in top5), str(tp)])
w = csv.writer(open(r'D:/ac结果/第一次返修/返修新材料/SI_度偏倚检验.csv', 'w', newline='', encoding='utf-8-sig'))
w.writerows(rows_out)
for r in rows_out: print(r)
