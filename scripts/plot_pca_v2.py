#!/usr/bin/env python3
"""Redraw the Fig4 metabolite PCA score plot from the exported coordinate source data."""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

CSV = r'D:/ac结果/fig4结果/代谢组PCA_坐标源数据.csv'
VAR = r'D:/ac结果/fig4结果/代谢组PCA_方差贡献.csv'
OUT = r'D:/ac结果/fig4结果/代谢组PCA_v2'

plt.rcParams.update({'font.family': 'Arial', 'font.size': 11, 'axes.titlesize': 15,
                     'axes.labelsize': 13, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
                     'legend.fontsize': 11, 'axes.linewidth': 0.9,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})

d = pd.read_csv(CSV)
v = pd.read_csv(VAR)
lab = {r['PC']: f"{r['PC']}({r['proportion_of_variance'] * 100:.2f}%)" for _, r in v.iterrows()}
ORDER = ['Control', '5nM', '25nM', '50nM', 'QC']
COL = {'Control': '#C77CFF', '5nM': '#00BFC4', '25nM': '#F8766D', '50nM': '#7CAE00',
       'QC': '#7A7A7A'}

fig = plt.figure(figsize=(6.0, 5.0), facecolor='white')
ax = fig.add_axes([0.13, 0.14, 0.62, 0.80])
ax.set_facecolor('white')
for g in ORDER:
    s = d[d.group == g]
    ax.scatter(s.PC1, s.PC2, s=58, c=COL[g], label=g, zorder=3, linewidths=0)
ax.set_xlabel(lab['PC1'])
ax.set_ylabel(lab['PC2'])
ax.set_title('PCA Scores Plot', pad=8)
ax.tick_params(direction='out', length=3.5, width=0.9)
ax.set_xlim(d.PC1.min() - 12, d.PC1.max() + 12)
ax.set_ylim(d.PC2.min() - 12, d.PC2.max() + 12)
leg = ax.legend(title='Condition', frameon=False, loc='center left',
                bbox_to_anchor=(1.02, 0.5))
leg.get_title().set_fontsize(13)
fig.savefig(OUT + '.pdf')
fig.savefig(OUT + '.png', dpi=300, facecolor='white')
print('wrote', OUT + '.png / .pdf | points plotted:', len(d),
      {g: int((d.group == g).sum()) for g in ORDER})
