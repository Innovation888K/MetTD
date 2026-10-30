"""B4v2 pathway figure COMPLETE REWRITE.
No path_effects, no tight bbox, explicit axes positions, compact layout.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np, math

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'pdf.fonttype': 42,
})

BLACK = '#000000'
C = {
    'ds': '#C4A0A0', 'dm': '#C4B0A4', 'dl': '#D0C4B8',
    'n': '#D4D4D4', 'enz': '#A0AEC0', 'tgt': '#C4766A', 'arr': '#9A9A9A',
}

BASE = 10000
# Compact pathway: just the essential metabolites, clean flow
# (abbrev, baseline, perturbed, is_d2)
# Scenario C values (ACO2 activated + IDH1 inhibited)
METS = [
    # key          abbrev    base   pert   d2    x      y
    ('OAA',        'OAA',    BASE,  BASE,  False, 0.0,   1.0),
    ('AcCoA',      'AcCoA',  BASE,  BASE,  False, 0.0,  -1.0),
    ('Cit',        'Cit',    BASE,  7500,  True,  2.0,   0.0),
    ('cAcon',      'cAcon',  BASE,  12500, True,  4.0,   0.0),
    ('IsoCit',     'IsoCit', BASE,  12000, True,  6.0,   0.0),
    ('aKG',        'aKG',    BASE,  6500,  False, 8.0,   0.0),
    ('SucCoA',     'SucCoA', BASE,  BASE,  False, 10.0,  0.0),
    ('Suc',        'Suc',    BASE,  9500,  False, 12.0,  0.0),
    ('OxSuc',      'OxSuc',  BASE,  8000,  False, 7.0,   1.5),
    ('Glu',        'Glu',    BASE,  9500,  False, 8.0,  -1.5),
]

# Local pathway scores + permutation p-values (10 metabolites × 5 enzymes, B=10000)
# (name, EC, score, rank, p_value, x, y, is_target)
ENZS = [
    ('CS',    '2.3.3.1',  0.45,  5,  0.768,  1.0,  0.0,  False),
    ('ACO2',  '4.2.1.3',  1.46,  1,  0.037,  3.0,  0.0,  True),
    ('2-MCD', '4.2.1.56', 1.19,  3,  0.035,  5.0,  0.0,  False),
    ('IDH1',  '1.1.1.42', 1.38,  2,  0.071,  7.0,  0.0,  False),
    ('aKGDH', '1.2.4.1',  0.75,  4,  0.288,  9.0,  0.0,  False),
]

# Reactions: (from_x, from_y, to_x, to_y, bidirectional)
ARROWS = [
    (0.0, 1.0, 2.0, 0.0, False),   # OAA → Cit
    (0.0, -1.0, 2.0, 0.0, False),  # AcCoA → Cit
    (2.0, 0.0, 4.0, 0.0, True),    # Cit ↔ cAcon
    (4.0, 0.0, 6.0, 0.0, False),   # cAcon → IsoCit
    (6.0, 0.0, 8.0, 0.0, False),   # IsoCit → aKG
    (6.0, 0.0, 7.0, 1.5, False),   # IsoCit → OxSuc
    (7.0, 1.5, 8.0, 0.0, False),   # OxSuc → aKG
    (8.0, 0.0, 8.0, -1.5, True),   # aKG ↔ Glu
    (8.0, 0.0, 10.0, 0.0, False),  # aKG → SucCoA
    (10.0, 0.0, 12.0, 0.0, False), # SucCoA → Suc
]

def get_color(fc):
    if fc <= 0.80: return C['ds'], 0.75
    elif fc <= 0.88: return C['dm'], 0.65
    elif fc <= 0.96: return C['dl'], 0.55
    else: return C['n'], 0.35

def draw_panel(ax, show_scores, label, subtitle):
    """Draw one pathway panel. All coordinates in data space."""
    # Draw arrows
    for x1, y1, x2, y2, bidir in ARROWS:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                   arrowprops=dict(arrowstyle='<->' if bidir else '->',
                                 color=C['arr'], lw=0.7, alpha=0.4,
                                 shrinkA=18, shrinkB=18))

    # Draw enzymes (rectangles)
    for name, ec, score, rank, p_val, ex, ey, is_t in ENZS:
        bw, bh = 1.3, 0.45
        rect = mpatches.FancyBboxPatch(
            (ex-bw/2, ey-bh/2), bw, bh,
            boxstyle='round,pad=0.05',
            facecolor=C['tgt'] if is_t else C['enz'],
            edgecolor='#6A6A6A',
            linewidth=1.2 if is_t else 0.5,
            alpha=0.85 if is_t else 0.4, zorder=3)
        ax.add_patch(rect)
        ax.text(ex, ey+0.02, name, ha='center', va='center', fontsize=7,
                fontweight='bold', color=BLACK, zorder=4)
        ax.text(ex, ey-0.15, ec, ha='center', va='center', fontsize=5,
                color=BLACK, zorder=4)
        if show_scores:
            # Format: score + rank + p-value with significance marker
            if p_val < 0.01:
                sig = '**'
            elif p_val < 0.05:
                sig = '*'
            elif p_val < 0.10:
                sig = '†'
            else:
                sig = ''
            ax.text(ex, ey-bh/2-0.08, f'S={score:.2f} #{rank} P={p_val:.3f}{sig}',
                    ha='center', va='top', fontsize=6,
                    color=BLACK, fontstyle='italic', zorder=4)

    # Draw metabolites (circles + labels)
    for key, abbrev, base, pert, is_d2, mx, my in METS:
        fc = pert / base
        if not show_scores:
            color, alpha = C['n'], 0.35; ab = base
        else:
            color, alpha = get_color(fc); ab = pert

        r = 0.30
        circ = plt.Circle((mx, my), r, facecolor=color,
                          edgecolor='#8A6A6A' if is_d2 else '#B0B0B0',
                          linewidth=1.5 if is_d2 else 0.7,
                          alpha=alpha, zorder=4)
        ax.add_patch(circ)

        if not show_scores:
            lbl = f'{abbrev}\n{ab:,.0f}'
        else:
            lbl = f'{abbrev}\n{ab:,.0f}' + (f'\nFC={fc:.2f}' if fc != 1.0 else '')

        # Label position: above for y>0, below for y<0, right for y=0 side-branches
        if my > 0.5:
            ly = my + r + 0.05; va = 'bottom'; ha = 'center'
        elif my < -0.5:
            ly = my - r - 0.05; va = 'top'; ha = 'center'
        else:
            # Main pathway (y=0): label above unless score label is below
            ly = my + r + 0.05; va = 'bottom'; ha = 'center'

        ax.text(mx, ly, lbl, ha=ha, va=va, fontsize=7,
                fontweight='bold' if (show_scores and fc != 1.0) else 'normal',
                fontstyle='italic', color=BLACK, zorder=5,
                linespacing=0.9)

    # Panel label
    ax.text(0.02, 0.98, label, transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=BLACK, va='top')
    ax.text(0.06, 0.98, subtitle, transform=ax.transAxes,
            fontsize=8, color=BLACK, va='top')

# Create figure: A4 width, moderate height
FW = 170 / 25.4  # 6.69 inches
FH = 160 / 25.4  # 6.30 inches
fig = plt.figure(figsize=(FW, FH))

# Explicit axes positions: [left, bottom, width, height] in figure fractions
# Top panel (A)
ax1 = fig.add_axes([0.05, 0.52, 0.93, 0.42])
# Bottom panel (B)
ax2 = fig.add_axes([0.05, 0.08, 0.93, 0.42])

# Set same data limits for both
for ax in [ax1, ax2]:
    ax.set_xlim(-1.0, 13.0)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect('auto')  # NOT equal — let matplotlib fit
    ax.axis('off')

draw_panel(ax1, False, 'A', '  Baseline: all metabolites at 10,000')
draw_panel(ax2, True,  'B', '  After ACO2 perturbation (pathway-local scores shown)')

# Figure-level title and legend
fig.text(0.5, 0.97, 'MetTD simulation: ACO2 perturbation on the TCA cycle',
         ha='center', fontsize=12, fontweight='bold', color=BLACK)
fig.text(0.5, 0.94, 'ACO2 ranks #1 (P=0.037*), IDH1 ranks #2 (P=0.071\u2020) among 5 pathway enzymes',
         ha='center', fontsize=9, color=BLACK, fontweight='bold')

# Legend
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C['ds'],
           markersize=8, alpha=0.75, label='FC\u22640.80'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C['dm'],
           markersize=8, alpha=0.65, label='FC 0.80\u20130.88'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C['dl'],
           markersize=8, alpha=0.55, label='FC 0.88\u20130.96'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C['n'],
           markersize=8, alpha=0.35, label='Unchanged'),
    mpatches.Patch(facecolor=C['tgt'], alpha=0.85, label='ACO2'),
    mpatches.Patch(facecolor=C['enz'], alpha=0.4, label='Enzyme'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='none',
           markeredgecolor='#8A6A6A', markersize=8, markeredgewidth=1.5,
           label='D=2 metabolite'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=7,
           fontsize=6, frameon=True, fancybox=True,
           bbox_to_anchor=(0.5, 0.005), columnspacing=0.8)

# Save to ASCII path only (Chinese path causes encoding issues in Python/PIL)
out_ascii = 'D:/MetTD_recompute/out/pathway_fig'
fig.savefig(out_ascii + '.pdf', facecolor='white', dpi=300)
fig.savefig(out_ascii + '.png', facecolor='white', dpi=300)
print(f"Saved: {out_ascii}.pdf / .png")
print(f"Figure: {FW*25.4:.0f} x {FH*25.4:.0f} mm")

# === PIXEL-LEVEL VERIFICATION ===
from PIL import Image
img = Image.open(out_ascii + '.png')
W, H = img.size
arr = np.array(img)
# Check margins (no content in outer 5 pixels)
margin = 5
top_strip = arr[:margin, :, :3]
bottom_strip = arr[-margin:, :, :3]
left_strip = arr[:, :margin, :3]
right_strip = arr[:, -margin:, :3]
white = np.array([255, 255, 255])
for name, strip in [('top', top_strip), ('bottom', bottom_strip),
                     ('left', left_strip), ('right', right_strip)]:
    non_white = np.any(strip != white, axis=-1).sum()
    if non_white > 0:
        print(f"WARNING: {non_white} non-white pixels in {name} margin")
    else:
        print(f"  {name} margin: clean")
print(f"Image: {W}x{H} pixels")
