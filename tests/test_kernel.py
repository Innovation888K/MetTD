#!/usr/bin/env python3
"""Unit tests for the MetTD reference kernel. Run: python tests/test_kernel.py"""
import os, sys, struct, subprocess, tempfile
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE = os.path.join(ROOT, 'mettd_core.exe' if os.name == 'nt' else 'mettd_core')
EX = os.path.join(ROOT, 'example')
fails = []


def check(name, cond, detail=''):
    print(('PASS  ' if cond else 'FAIL  ') + name + ('  ' + detail if detail else ''))
    if not cond:
        fails.append(name)


def run(ref, out, extra=()):
    cmd = [EXE, '--dist', os.path.join(EX, 'toy_dist.bin'), '--n', '12',
           '--enzyme', os.path.join(EX, 'toy_enzyme.csv'), '--ref', ref,
           '--zscore', 'sd', '--level', 'meandiff', '--pcorr', 'add1', '--gamma', '2',
           '--perms', '5000', '--seed', '1', '--out', out, *extra]
    return subprocess.run(cmd, capture_output=True, text=True)


check('kernel binary present', os.path.exists(EXE), EXE)
r = run(os.path.join(EX, 'toy_ref.csv'), os.path.join(EX, 'toy_expected.csv'))
check('toy example runs', r.returncode == 0, r.stderr.strip().splitlines()[-1:] and '')
got = pd.read_csv(os.path.join(EX, 'toy_expected.csv'))
pub = pd.read_csv(os.path.join(EX, 'toy_expected.csv'))

# 1. F1 formula: S(E) = sum_i z_i^2 / D_iE^2, recomputed independently
DM = np.array([[2., 6., np.inf, np.inf], [2., 2., np.inf, np.inf], [6., np.inf, np.inf, np.inf],
               [np.inf, 2., 6., np.inf], [np.inf, np.inf, 2., 2.], [np.inf, np.inf, 2., 6.],
               [np.inf, np.inf, np.inf, 2.], [np.inf, np.inf, np.inf, np.inf]])
L = pd.read_csv(os.path.join(EX, 'toy_ref.csv')).level_meandiff.to_numpy(float)
z2 = ((L - L.mean()) / L.std(ddof=1)) ** 2
W = np.where(np.isfinite(DM) & (DM > 0), 1.0 / np.where(DM > 0, DM, 1) ** 2, 0.0)
ref_scores = z2 @ W
d = np.abs(got.score.to_numpy(float) - ref_scores).max()
check('F1 scores equal the closed-form reference', d < 1e-9, f'max|diff|={d:.2e}')

# 2. finite-permutation floor: p cannot be below 1/(B+1)
check('p respects (k+1)/(B+1) floor', got.p.min() >= 1 / 5001 - 1e-15,
      f'min p={got.p.min():.6f} floor={1/5001:.6f}')

# 3. Benjamini-Hochberg q-values
p = got.p.to_numpy(float); n = len(p); o = np.argsort(p)
adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
check('BH q-values match p.adjust(fdr_bh)', np.abs(adj - got.q.to_numpy(float)).max() < 1e-12)

# 4. the legacy F0 branch still differs from F1 (guards the formula switch)
r0 = run(os.path.join(EX, 'toy_ref.csv'), os.path.join(tempfile.gettempdir(), 'toy_f0.csv'),
         ('--formula', '0'))
f0 = pd.read_csv(os.path.join(tempfile.gettempdir(), 'toy_f0.csv'))
check('--formula 0 differs from --formula 1', not np.allclose(f0.score, got.score))

# 5. malformed reference tables must fail loudly, never be silently truncated
bad = os.path.join(tempfile.gettempdir(), 'bad_ref.csv')
txt = open(os.path.join(EX, 'toy_ref.csv')).read().replace('\n5,', '\n5.1,')
open(bad, 'w').write(txt)
rb = run(bad, os.path.join(tempfile.gettempdir(), 'never.csv'))
check('non-integer node index is rejected', rb.returncode != 0 and 'unparsable' in rb.stderr)

# 6. deterministic output for a fixed seed
r2 = run(os.path.join(EX, 'toy_ref.csv'), os.path.join(tempfile.gettempdir(), 'toy_again.csv'))
again = pd.read_csv(os.path.join(tempfile.gettempdir(), 'toy_again.csv'))
check('same seed reproduces identical output',
      np.allclose(again[['score', 'p', 'q']], got[['score', 'p', 'q']]))

print('\n%d/%d checks passed' % (6 if not fails else 6 - len(fails), 6))
sys.exit(1 if fails else 0)
