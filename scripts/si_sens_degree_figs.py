# -*- coding: utf-8 -*-
"""SI sensitivity + degree-bias figures for MetTD revision (R2.M2/M3)."""
import csv, re, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "font.size": 7.5,
    "axes.linewidth": 0.7,
    "axes.unicode_minus": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
})

# Morandi palette
M_GREEN = "#8EA69B"; M_ROSE = "#C7A9A3"; M_MAUVE = "#9A8FA8"
M_SAND  = "#C9B98F"; M_GRAY = "#7D7D7D"; M_BLUE = "#93A8B8"
INK = "#1A1A1A"

OUT = "D:/ac结果/第一次返修/返修新材料/"
MM = 1 / 25.4

# ---------------- data: sensitivity table ----------------
settings = ["Reference (γ=2, B=5000)", "B=100", "B=500", "B=1000", "γ=1", "γ=3",
            "radius ≤ 6 (3 hops)", "radius ≤ 10 (5 hops)", "− currency (131 nodes)",
            "High-conf subset (64% feats)", "Dedupe (unique nodes)"]
aco2 = [1, 1, 1, 1, 1, 1, 3, 1, 1, 205, 1]
idh1 = [3, 3, 3, 3, 77, 4, 23, 3, 3, 989, 2]
cs   = [71, 71, 71, 71, 1105, 9, 168, 72, 71, 473, 22]
qcnt = [535, 535, 535, 535, 531, 510, 399, 497, 316, 716, 596]

fig = plt.figure(figsize=(174 * MM, 105 * MM), dpi=300)
axA = fig.add_axes([0.235, 0.155, 0.310, 0.745])
axB = fig.add_axes([0.630, 0.155, 0.330, 0.745])

ypos = list(range(len(settings)))
# panel A: rank stability (log x, rank 1 at left, reference at top)
axA.plot(aco2, ypos, "o-", color=M_GREEN, lw=1.4, ms=4.2, label="ACO2 (4.2.1.3)", zorder=5)
axA.plot(idh1, ypos, "s-", color=M_ROSE, lw=1.4, ms=3.8, label="IDH1 (1.1.1.42)", zorder=4)
axA.plot(cs, ypos, "^-", color=M_MAUVE, lw=1.4, ms=4.0, label="CS (2.3.3.1)", zorder=3)
axA.set_xscale("log")
axA.set_xlim(2400, 0.8)  # inverted: rank 1 at left
axA.set_xticks([1, 3, 10, 30, 100, 300, 1000])
axA.xaxis.set_minor_locator(FixedLocator([]))
axA.set_yticks(ypos)
axA.set_yticklabels(settings, fontsize=6.5, color=INK)
axA.set_ylim(10.6, -0.6)
axA.axhline(0.5, color="#BBBBBB", lw=0.6, ls=(0, (3, 2)))
axA.annotate("205", (205, 9), xytext=(5, 2), textcoords="offset points",
             ha="left", fontsize=6.2, color=M_GREEN)
axA.annotate("989", (989, 9), xytext=(5, 2), textcoords="offset points",
             ha="left", fontsize=6.2, color=M_ROSE)
axA.annotate("77", (77, 4), xytext=(5, 3), textcoords="offset points",
             ha="left", fontsize=6.2, color=M_ROSE)
axA.annotate("1105", (1105, 4), xytext=(5, -5), textcoords="offset points",
             ha="left", fontsize=6.2, color=M_MAUVE)
axA.set_xlabel("MetTD rank (log scale)", color=INK)
axA.set_title("A  Target-rank stability across 11 settings", loc="left",
              fontsize=8.5, color=INK, pad=6)
axA.legend(fontsize=6.5, frameon=False, loc="lower right", handlelength=1.4)
axA.grid(axis="x", color="#DDDDDD", lw=0.5, zorder=0)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)

# panel B: significant-enzyme counts per setting (horizontal bars)
cols = [M_GRAY] + [M_BLUE] * 10
cols[9] = M_SAND
axB.barh(ypos, qcnt, height=0.62, color=cols, edgecolor="none", zorder=3)
for yi, v in zip(ypos, qcnt):
    axB.annotate(str(v), (v, yi), xytext=(3, 0), textcoords="offset points",
                 va="center", ha="left", fontsize=5.8, color=INK)
axB.set_xlim(0, 840)
axB.set_xticks([0, 200, 400, 600, 800])
axB.set_yticks(ypos)
axB.set_yticklabels([])
axB.set_ylim(10.6, -0.6)
axB.axhline(0.5, color="#BBBBBB", lw=0.6, ls=(0, (3, 2)))
axB.set_xlabel("Enzymes with q < 0.05", color=INK)
axB.set_title("B  FDR-significant enzyme count per setting", loc="left",
              fontsize=8.5, color=INK, pad=6)
axB.grid(axis="x", color="#DDDDDD", lw=0.5, zorder=0)
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)

fig.savefig(OUT + "SI_敏感性分析图.pdf")
fig.savefig(OUT + "SI_敏感性分析图.png")
plt.close(fig)

# ---------------- data: degree-bias table ----------------
rows = list(csv.DictReader(open(OUT + "SI_度偏倚检验.csv", encoding="utf-8-sig")))
names = [r["dataset"] for r in rows]
rho = [float(r["spearman_degree_score"]) for r in rows]

def grp(n):
    if n.startswith("drug_"): return M_GREEN
    if n.startswith("met"):   return M_ROSE
    if n.startswith("t"):     return M_MAUVE
    return M_GRAY
colors = [grp(n) for n in names]
order = sorted(range(len(names)), key=lambda i: rho[i])
pcts = []
for r in rows:
    pcts += [int(v) for v in re.findall(r"\d+", r["top5_degree_percentile"])]

fig = plt.figure(figsize=(174 * MM, 118 * MM), dpi=300)
axC = fig.add_axes([0.145, 0.170, 0.520, 0.755])
axD = fig.add_axes([0.720, 0.170, 0.225, 0.755])

def short_name(n):
    if n.startswith("drug_"): return n[5:]
    if n.startswith("fig3c"): return "fig3c"
    return n
short = [short_name(n) for n in names]

ys = range(len(order))
for y0, i in zip(ys, order):
    axC.plot([0, rho[i]], [y0, y0], color=colors[i], lw=1.2, zorder=3)
    axC.plot(rho[i], y0, "o", color=colors[i], ms=4.0, zorder=4)
axC.axvspan(-0.1, 0.1, color="#EFEFEF", lw=0, zorder=1)
axC.axvline(0, color=INK, lw=0.7, zorder=2)
axC.set_yticks(list(ys))
axC.set_yticklabels([short[i] for i in order], fontsize=6.2, color=INK)
axC.set_xlim(-0.15, 0.30)
axC.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
axC.set_xlabel("Spearman ρ (enzyme degree vs MetTD score)", color=INK)
axC.set_title("C  Degree–score correlation per dataset (24 datasets)", loc="left",
              fontsize=8.5, color=INK, pad=6)
axC.annotate("median +0.020\nmax +0.239", (0.239, order[-1]),
             xytext=(-2, -14), textcoords="offset points", ha="right",
             fontsize=6.2, color=INK)
axC.grid(axis="x", color="#DDDDDD", lw=0.5, zorder=0)
for s in ("top", "right", "left"):
    axC.spines[s].set_visible(False)

# panel D: top-5 degree percentile histogram
bins = list(range(10, 111, 10))
counts, _, _ = axD.hist(pcts, bins=bins, color=M_BLUE, edgecolor="white", lw=0.6, zorder=3)
med = st.median(pcts)
ymax = max(counts) * 1.18
axD.set_ylim(0, ymax)
axD.set_yticks([0, 20, 40, 60])
axD.axvline(med, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=4)
axD.annotate("median 60", (med, max(counts) * 1.12), xytext=(3, 0),
             textcoords="offset points", va="top", fontsize=6.0, color=INK)
axD.set_xlim(10, 110)
axD.set_xticks([10, 30, 50, 70, 90, 110])
axD.set_xlabel("Degree percentile", color=INK, fontsize=6.8)
axD.set_ylabel("Top-5 enzymes (n = 120)", color=INK)
axD.set_title("D", loc="left", fontsize=8.5, color=INK, pad=6)
axD.grid(axis="y", color="#DDDDDD", lw=0.5, zorder=0)
for s in ("top", "right"):
    axD.spines[s].set_visible(False)

fig.savefig(OUT + "SI_度偏倚检验图.pdf")
fig.savefig(OUT + "SI_度偏倚检验图.png")
plt.close(fig)
print("done")
