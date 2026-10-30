# MetTD 1.1.0

Version used for the revised manuscript.

## Behaviour changes relative to 1.0.x
- Score is now `S(E) = sum_i z_i^2 / D_iE^2` (each standardized change squared
  individually, then weighted by inverse squared distance). The previous
  `(sum_i w_i z_i / sum_i w_i)^2` remains available as `--formula 0` /
  `formula = "F0"` for regenerating older result files.
- Permutation p-values use `(k+1)/(B+1)`; the smallest attainable p with B = 5,000
  is 2.0e-04. `--pcorr old` reproduces the previous `max(k,1)/B`.
- Benjamini-Hochberg q-values across all 1,126 enzymes are reported with every
  result table.

## Defects fixed
- `read_dist()` requested a 4 MiB block and treated a short return as an error, so
  the kernel aborted on any distance file smaller than 4 MiB (it only worked on the
  full 720 MB matrix). Small test networks now load correctly.
- `read_ref()` silently skipped rows whose node index was not an integer. Reference
  tables produced by the original `export_ref.R` wrote indices as `5211.1`, so those
  metabolites were dropped without warning - 1,508 of 2,783 rows in one dataset. The
  parser now aborts with an explicit message.

## Added
- `src/mettd_core.c` reference kernel; `example/` toy network with expected output;
  `tests/` covering the equation, the p floor, BH correction, seed determinism and
  input validation.
- `data/network/` (nodes, enzymes, node-id map, EC-gene map, distance-matrix checksum
  and regeneration script) and `data/benchmark/` (labels, evidence tiers, ranks,
  per-record metadata).
- `scripts/` for every figure and SI table, `REPRODUCIBILITY.md`, `SEEDS.md`,
  `requirements.txt`, `environment.yml`, `R_sessionInfo.txt`, MIT `LICENSE`,
  `CITATION.cff`.
