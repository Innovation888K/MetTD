#!/usr/bin/env python3
"""Generate the toy example (distance matrix, enzyme table, reference table) and,
with --run, the expected output produced by the compiled C kernel.

Run:  python example/make_example.py --run --exe ./mettd_core
"""
import sys, os, struct, subprocess, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
N_ENZ, N_MET = 4, 8
N = N_ENZ + N_MET                      # nodes 0..3 are enzymes, 4..11 are metabolites
# distance from metabolite j to enzyme i; np.inf = unreachable in the bipartite graph
DM = np.array([[2., 6., np.inf, np.inf],
               [2., 2., np.inf, np.inf],
               [6., np.inf, np.inf, np.inf],
               [np.inf, 2., 6., np.inf],
               [np.inf, np.inf, 2., 2.],
               [np.inf, np.inf, 2., 6.],
               [np.inf, np.inf, np.inf, 2.],
               [np.inf, np.inf, np.inf, np.inf]])
LEVELS = np.array([3.1, -2.4, 1.8, 2.9, -1.1, 0.7, -2.8, 1.4])

ap = argparse.ArgumentParser()
ap.add_argument('--run', action='store_true')
ap.add_argument('--exe', default='./mettd_core')
a = ap.parse_args()

DENSE = np.full((N, N), np.inf)
np.fill_diagonal(DENSE, 0.0)
for j in range(N_MET):
    for i in range(N_ENZ):
        d = DM[j, i]
        DENSE[N_ENZ + j, i] = d
        DENSE[i, N_ENZ + j] = d
with open(os.path.join(HERE, 'toy_dist.bin'), 'wb') as f:
    f.write(struct.pack('<%dd' % (N * N), *DENSE.T.reshape(-1)))   # column-major, as read
pd.DataFrame({'idx': np.arange(1, N_ENZ + 1),
              'enzyme': ['1.1.1.1', '2.7.1.1', '4.2.1.3', '1.1.1.42']}).to_csv(
    os.path.join(HERE, 'toy_enzyme.csv'), index=False)
pd.DataFrame({'idx': N_ENZ + 1 + np.arange(N_MET), 'level_meandiff': LEVELS,
              'level_fc': LEVELS / 10.0}).to_csv(os.path.join(HERE, 'toy_ref.csv'),
                                                 index=False)
print('wrote toy_dist.bin (%d x %d doubles, column-major), toy_enzyme.csv (%d enzymes), '
      'toy_ref.csv (%d metabolites)' % (N, N, N_ENZ, N_MET))

if a.run:
    out = os.path.join(HERE, 'toy_expected.csv')
    cmd = [a.exe, '--dist', os.path.join(HERE, 'toy_dist.bin'), '--n', str(N),
           '--enzyme', os.path.join(HERE, 'toy_enzyme.csv'),
           '--ref', os.path.join(HERE, 'toy_ref.csv'), '--zscore', 'sd',
           '--level', 'meandiff', '--pcorr', 'add1', '--gamma', '2', '--formula', '1',
           '--perms', '5000', '--seed', '1', '--out', out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    sys.stderr.write(r.stderr[-300:])
    x = pd.read_csv(out)
    print('\ntoy_expected.csv:')
    print(x.to_string(index=False))
    # independent re-derivation of the F1 scores, no shared code with the kernel
    z = (LEVELS - LEVELS.mean()) / LEVELS.std(ddof=1)
    W = np.where(np.isfinite(DM) & (DM > 0), 1.0 / np.where(DM > 0, DM, 1) ** 2, 0.0)
    ref = z ** 2 @ W
    got = x['score'].to_numpy(float)
    print('\nindependent numpy check: max |kernel - numpy| = %.3e' % np.abs(got - ref).max())
    assert np.abs(got - ref).max() < 1e-9, 'kernel disagrees with the reference formula'
    print('p floor with B=5000 is 1/5001 = %.6f; observed min p = %.6f'
          % (1 / 5001, x['p'].min()))
    assert x['p'].min() >= 1 / 5001 - 1e-12, 'p below the (k+1)/(B+1) floor'
    print('OK: equations, permutation floor and example output are consistent.')
