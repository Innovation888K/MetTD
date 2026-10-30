"""Local pathway permutation test: ONLY the 10 metabolites and 5 enzymes shown in the pathway figure.
No extra metabolites from the network. Clean, minimal, exactly what's in the diagram.
"""
import sys, math, random
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

# === EXACTLY the 10 metabolites from the pathway figure ===
# (name, baseline, perturbed)
METS = [
    ('OAA',    10000, 10000),
    ('AcCoA',  10000, 10000),
    ('Cit',    10000,  7500),
    ('cAcon',  10000, 12500),
    ('IsoCit', 10000, 12000),
    ('aKG',    10000,  6500),
    ('SucCoA', 10000, 10000),
    ('Suc',    10000,  9500),
    ('OxSuc',  10000,  8000),
    ('Glu',    10000,  9500),
]

# === EXACTLY the 5 enzymes from the pathway figure ===
# (name, EC, which metabolites at D=2)
# D=2 means directly connected (substrate/product)
ENZYMES = [
    ('CS',    '2.3.3.1',  ['OAA', 'AcCoA', 'Cit']),           # CS: OAA+AcCoA→Cit
    ('ACO2',  '4.2.1.3',  ['Cit', 'cAcon', 'IsoCit']),         # ACO2: Cit→cAcon→IsoCit
    ('2-MCD', '4.2.1.56', ['cAcon', 'IsoCit']),                # 2-MCD: cAcon→IsoCit (parallel)
    ('IDH1',  '1.1.1.42', ['IsoCit', 'OxSuc', 'aKG']),         # IDH1: IsoCit(+OxSuc)→aKG
    ('aKGDH', '1.2.4.1',  ['aKG', 'SucCoA']),                  # aKGDH: aKG→SucCoA
]

# D=6 metabolites (indirect, one enzyme hop away)
# In the real bipartite network, these are 3 edges away
D6_MAP = {
    'CS':    ['IsoCit', 'aKG', 'Suc'],          # 2 hops through ACO2/IDH1
    'ACO2':  ['OAA', 'AcCoA', 'aKG', 'OxSuc', 'Suc', 'Glu', 'SucCoA'],
    '2-MCD': ['Cit', 'OAA', 'AcCoA', 'aKG', 'OxSuc'],
    'IDH1':  ['Cit', 'cAcon', 'SucCoA', 'Suc', 'Glu'],
    'aKGDH': ['IsoCit', 'OxSuc', 'Suc', 'Glu'],
}

# Compute z-scores from the 10 metabolites
levels = np.array([m[2] - m[1] for m in METS])
met_names = [m[0] for m in METS]
mean_l = levels.mean()
var_l = levels.var(ddof=1)
sd_l = math.sqrt(var_l)
z = (levels - mean_l) / sd_l
z_dict = dict(zip(met_names, z))

print("10 metabolites, z-scores:")
for m, zv in zip(met_names, z):
    print(f"  {m:8s} Δ={dict(zip(met_names,levels))[m]:+.0f}  z={zv:+.3f}")

# Distance weights: D=2 → w=1/4, D=6 → w=1/36
def weight(dist):
    return 1.0 / (dist ** 2)

def score_enzyme(enzy_name, z_map):
    """Formula 1: S = Σ(z²/D²), summing over D=2 and D=6 metabolites."""
    enz = [e for e in ENZYMES if e[0] == enzy_name][0]
    s = 0.0
    # D=2 metabolites
    for met in enz[2]:
        if met in z_map:
            s += weight(2) * z_map[met] ** 2
    # D=6 metabolites
    for met in D6_MAP.get(enzy_name, []):
        if met in z_map:
            s += weight(6) * z_map[met] ** 2
    return s

# Observed scores
print(f"\nObserved scores (S = Σ z²/D², 10 metabolites only):")
obs = {}
for name, ec, _ in ENZYMES:
    obs[name] = score_enzyme(name, z_dict)
    print(f"  {name:6s} score = {obs[name]:.4f}")

obs_order = sorted(obs, key=lambda e: -obs[e])
print(f"\nRanking: {' > '.join(obs_order)}")

# === PERMUTATION TEST ===
# Shuffle z-scores among the 10 metabolites, recompute 5 enzyme scores
B = 10000
rng = random.Random(42)

rank_count = {e: {r: 0 for r in range(1, 6)} for e in obs}
perm_scores = {e: [] for e in obs}

for b in range(B):
    z_vals = list(z)
    rng.shuffle(z_vals)
    z_perm = dict(zip(met_names, z_vals))

    sc = {}
    for name, _, _ in ENZYMES:
        sc[name] = score_enzyme(name, z_perm)

    order = sorted(sc, key=lambda e: -sc[e])
    for rank, e in enumerate(order, 1):
        rank_count[e][rank] += 1
    for e in sc:
        perm_scores[e].append(sc[e])

print(f"\n{'='*65}")
print(f"PERMUTATION TEST ({B} shuffles, 10 metabolites × 5 enzymes)")
print(f"{'='*65}")
print(f"\n{'Enzyme':8s} {'Obs':>8s} {'PermMean':>9s} {'P(≥obs)':>8s} {'P(#1)':>7s} {'MeanRank':>9s}")
print("-" * 55)

for e in obs_order:
    ps = perm_scores[e]
    p_ge = sum(1 for s in ps if s >= obs[e]) / B
    p_1 = rank_count[e][1] / B
    mr = sum(r * c for r, c in rank_count[e].items()) / B
    print(f"{e:8s} {obs[e]:>8.4f} {np.mean(ps):>9.4f} {p_ge:>8.4f} {p_1:>7.4f} {mr:>9.2f}")

print(f"\n{'Enzyme':8s} {'#1':>7s} {'#2':>7s} {'#3':>7s} {'#4':>7s} {'#5':>7s}")
print("-" * 45)
for e in obs_order:
    fr = [rank_count[e][r]/B for r in range(1,6)]
    print(f"{e:8s} " + " ".join(f"{f:>7.3f}" for f in fr))
