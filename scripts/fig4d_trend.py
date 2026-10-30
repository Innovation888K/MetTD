"""Fig4d rerun: significant-enzyme counts across time x dose + key-enzyme trajectories.
Reads out/b2/t*.csv (production kernel outputs)."""
import csv, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding='utf-8')
os.chdir(r'D:/MetTD_recompute')
times = ['0.5', '2', '4', '8', '24']
doses = {'l': '5 uM', 'm': '25 uM', 'h': '50 uM'}
KEYS = [('1.1.1.42', 'IDH1'), ('4.2.1.3', 'ACO2'), ('2.3.3.1', 'CS')]

conds = {}
for t in times:
    for d in 'lmh':
        f = f'out/b2/t{t}{d}.csv'
        if os.path.exists(f):
            rows = list(csv.DictReader(open(f)))
            conds[(t, d)] = rows

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
for d, lab in doses.items():
    xs, ys = [], []
    for t in times:
        if (t, d) in conds:
            xs.append(float(t) if t != '0.5' else 0.5)
            ys.append(sum(1 for r in conds[(t, d)] if float(r['q']) < 0.05))
    ax.plot(xs, ys, 'o-', label=lab)
ax.set_xscale('log'); ax.set_xticks([0.5, 2, 4, 8, 24]); ax.set_xticklabels(['0.5', '2', '4', '8', '24'])
ax.set_xlabel('Time (h)'); ax.set_ylabel('Enzymes with q < 0.05'); ax.legend(title='Metformin')
ax.set_title('Significant enzymes over time x dose (rerun)')

ax = axes[1]
for ec, nm in KEYS:
    for d, lab in doses.items():
        xs, ys = [], []
        for t in times:
            if (t, d) in conds:
                r0 = [r for r in conds[(t, d)] if r['enzyme'] == ec]
                if r0:
                    order = sorted(range(len(conds[(t, d)])),
                                   key=lambda i: -float(conds[(t, d)][i]['score']))
                    rank = order.index([i for i, r in enumerate(conds[(t, d)]) if r['enzyme'] == ec][0]) + 1
                    xs.append(float(t) if t != '0.5' else 0.5); ys.append(rank)
        if xs:
            ax.plot(xs, ys, 'o--', label=f'{nm} {lab}', alpha=0.8)
ax.set_yscale('log'); ax.invert_yaxis()
ax.set_xscale('log'); ax.set_xticks([0.5, 2, 4, 8, 24]); ax.set_xticklabels(['0.5', '2', '4', '8', '24'])
ax.set_xlabel('Time (h)'); ax.set_ylabel('Score rank (1126 enzymes)')
ax.legend(fontsize=7, ncol=2); ax.set_title('Key TCA enzyme trajectories')
fig.tight_layout()
fig.savefig(r'D:/ac结果/第一次返修/返修新材料/Fig4d_rerun_trend.png', dpi=300)
fig.savefig(r'D:/ac结果/第一次返修/返修新材料/Fig4d_rerun_trend.pdf')

# ST1 sheet-7 style table
w = csv.writer(open(r'D:/ac结果/第一次返修/返修新材料/ST1_伪时序重跑_q计数.csv', 'w', newline='', encoding='utf-8-sig'))
w.writerow(['condition', 'q<0.05', 'p<0.05', 'top1_enzyme'])
for (t, d), rows in sorted(conds.items(), key=lambda x: (float(x[0][0]), x[0][1])):
    q = sum(1 for r in rows if float(r['q']) < 0.05)
    p = sum(1 for r in rows if float(r['p']) < 0.05)
    top = max(rows, key=lambda r: float(r['score']))['enzyme']
    w.writerow([f'{t}h-{doses[d]}', q, p, top])
print('conditions plotted:', len(conds))
for (t, d) in sorted(conds, key=lambda x: (float(x[0]), x[1])):
    q = sum(1 for r in conds[(t, d)] if float(r['q']) < 0.05)
    print(f'  {t}h {doses[d]}: q<0.05 = {q}')
