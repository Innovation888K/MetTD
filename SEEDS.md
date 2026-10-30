# Random seeds and stochastic settings

Every stochastic quantity in the manuscript is fixed and listed here.

| Use | Setting | Where |
|---|---|---|
| Network permutation null (all published p/q values) | `--seed 1`, B = 5000 | `src/mettd_core.c`, all runs |
| Permutation estimator | (k+1)/(B+1), so the smallest attainable p is 1/5001 = 2.000e-04 | `--pcorr add1` |
| Distance matrix | deterministic function of the BioModels SBML; no randomness | `data/network/` |
| XGBoost baseline | `random_state=123` | `scripts/fig2bc_v3_losoc.py` |
| Random forest baseline | `random_state=123` | same |
| SVM baseline | `random_state=1` | same |
| BP (MLP) baseline | `random_state=1` | same |
| Cross-validation splits | leave-one-study-out (deterministic given the study labels); the contrast-wise comparison uses `StratifiedKFold(5, shuffle=True, random_state=123)` | same |
| Bootstrap CIs (Fig. 2 metrics) | `numpy.random.default_rng(1)`, 10,000 resamples over drugs | `scripts/b6_precision.py` |
| Pathway-local simulation (Fig. R1) | `random.Random(42)`, 10,000 shuffles | `scripts/pathway_perm_figR1.py` |

Re-running `python tests/test_kernel.py` after a rebuild must reproduce
`example/toy_expected.csv` byte-for-byte, because the kernel seeds its Mersenne
stream from `--seed`.
