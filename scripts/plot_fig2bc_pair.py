#!/usr/bin/env python3
"""Two-panel Fig2b+Fig2c for benchmark_v3 under the agreed protocol:
MetTD scored on all contrasts (no split); ML methods under leave-one-study-out CV."""
import sys, json
sys.stdout.reconfigure(encoding='utf-8')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9,
    'axes.labelsize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 7.5, 'axes.linewidth': 0.8, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'xtick.direction': 'out', 'ytick.direction': 'out'})

R = json.load(open('D:/MetTD_recompute/fig2bc_v3_both.json'))
ORDER = ['MetTD', 'XGBoost', 'SVM', 'BP', 'RF']
COL = {'MetTD': '#4E6A8F', 'XGBoost': '#C0726A', 'SVM': '#7FA07A',
       'BP': '#D8A657', 'RF': '#9B87B5'}
KEYS = ['acc', 'prec', 'rec', 'f1', 'spec']
LAB = ['Accuracy', 'Precision', 'Recall', 'F1-score', 'Specificity']
PROTO = sys.argv[1] if len(sys.argv) > 1 else 'LOSO'
SUF = sys.argv[2] if len(sys.argv) > 2 else ''
blk = R[PROTO]

fig = plt.figure(figsize=(17.4 / 2.54, 7.6 / 2.54), facecolor='white')

axA = fig.add_axes([0.070, 0.125, 0.380, 0.800])
axA.plot([0, 1], [0, 1], ls=(0, (4, 3)), lw=0.8, c='#AEAEAE', zorder=1)
for k in ORDER:
    r = blk['roc'][k]
    axA.plot(r['fpr'], r['tpr'], lw=1.6, c=COL[k], zorder=3,
             label=f"{k} ({blk['metrics'][k]['auc']:.2f})")
axA.set_xlabel('1 - Specificity')
axA.set_ylabel('Sensitivity')
axA.set_xlim(0, 1); axA.set_ylim(0, 1)
axA.set_xticks(np.arange(0, 1.01, 0.2)); axA.set_yticks(np.arange(0, 1.01, 0.2))
axA.legend(frameon=False, loc='lower right', handlelength=1.7)
axA.text(-0.115, 1.020, 'b', transform=axA.transAxes, fontsize=12)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)

axB = fig.add_axes([0.555, 0.125, 0.425, 0.800])
M = np.array([[blk['metrics'][k][x] for x in KEYS] for k in ORDER])
w = 0.15
x = np.arange(len(KEYS))
for i, k in enumerate(ORDER):
    xs = x + (i - 2) * w
    axB.bar(xs, M[i], width=w, color=COL[k], ec='white', lw=0.4, label=k)
    for xi, v in zip(xs, M[i]):
        axB.text(xi, v + 0.012, f'{v:.2f}', ha='center', va='bottom', fontsize=5.4)
axB.set_xticks(x); axB.set_xticklabels(LAB)
axB.set_ylabel('Value')
axB.set_ylim(0, 1.30)
axB.set_yticks(np.arange(0, 1.01, 0.2))
axB.legend(frameon=False, ncol=3, loc='upper center', handlelength=1.0,
           columnspacing=1.0, borderaxespad=0.2, labelspacing=0.35)
axB.text(-0.075, 1.020, 'c', transform=axB.transAxes, fontsize=12)
for s in ('top', 'right'):
    axB.spines[s].set_visible(False)

out = f'D:/MetTD_recompute/out/Fig2bc_v3_{PROTO}{SUF}'
fig.savefig(out + '.pdf')
fig.savefig(out + '.png', dpi=600, facecolor='white')
print('wrote', out + '.png / .pdf')
for k in ORDER:
    m = blk['metrics'][k]
    print(f"  {k:8s} AUC={m['auc']:.3f} Acc={m['acc']:.2f} Prec={m['prec']:.2f} "
          f"Rec={m['rec']:.2f} F1={m['f1']:.2f} Spec={m['spec']:.2f}")
