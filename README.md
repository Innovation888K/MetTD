# MetTD

Metabolomics-based Target iDentification: prioritises the enzymes whose local
metabolic neighbourhood is perturbed, by combining standardized metabolite abundance
changes with metabolite-enzyme topological distance in a genome-scale reaction network.



## Equation implemented

For an enzyme *E* and the *n* metabolites mapped in one comparison,

```
S(E) = sum_{i=1..n}  z_i^2 / D_iE^2
```

where `z_i = (Level_i - mean Level) / sd(Level)` (sample sd, ddof = 1) and
`Level_i` is the model-group mean minus the baseline-group mean. `D_iE` is the network
distance, twice the shortest-path hop count in the enzyme-metabolite bipartite graph,
so a directly connected pair has `D = 2`. Each term is squared **individually before**
weighting, so the score is a non-negative, direction-agnostic magnitude. Unreachable
pairs contribute weight zero; zero-distance entries are excluded.

Significance uses a network permutation test: in each of B = 5,000 permutations the
`n_ref x n_enz` distances actually used are resampled without replacement from the full
9,492 x 9,492 distance matrix, and p = (k+1)/(B+1). p-values are corrected across the
1,126 enzymes by Benjamini-Hochberg.

## What changed in 1.1.0

| Item | 1.0.x (original submission) | 1.1.0 (this release) |
|---|---|---|
| Score | `(sum_i w_i z_i / sum_i w_i)^2` | `sum_i z_i^2 / D_iE^2` |
| Permutation p | `max(k,1)/B` | `(k+1)/(B+1)` |
| Multiple testing | not reported | BH q across all 1,126 enzymes |
| Reference tables | node indices written as `5211.1`, silently dropped by the parser | integer indices; the parser now aborts on a non-integer index instead of dropping the row |
| Distance-matrix reader | requested 4 MiB blocks and aborted on smaller files | reads exactly the remaining bytes |

The legacy behaviours remain reachable (`--formula 0`, `--pcorr old` in the C kernel;
`formula = "F0"`, `correction = "old"` in R) so that earlier result files can still be
regenerated, but they are no longer the defaults and are not used for any number in the
revised manuscript.

## Layout

```
main/MetTD.R              R implementation (scoring, permutation, plotting, mapping)
main/Single_work_flow.R   end-to-end worked example in R
src/mettd_core.c          reference C kernel used for every published number
data/network/             nodes.csv, enzyme.csv, nodeid_map.csv, EC_gene_map.csv,
                          DISTANCE_MATRIX.md (checksum + how to rebuild model_dist.bin)
data/benchmark/           benchmark_v3.csv, evidence tiers, per-record MetTD ranks,
                          benchmark_metadata_v3.csv (subject groups, tiers, rationale)
scripts/                  figure- and table-producing scripts (see REPRODUCIBILITY.md)
example/                  toy 12-node network: input, expected output, generator
tests/                    unit tests for equation, p floor, BH, seeds, input validation
data_anal/                machine-learning comparators and Top-N hit-rate analysis
```

## Build and quickstart

```sh
gcc -O2 -std=c11 -o mettd_core src/mettd_core.c -lm
python example/make_example.py --run --exe ./mettd_core   # regenerates and checks the
                                                          # expected output
python tests/test_kernel.py
```

Full-scale run (requires `model_dist.bin`, see `data/network/DISTANCE_MATRIX.md`):

```sh
./mettd_core --dist model_dist.bin --n 9492 \
  --enzyme data/network/enzyme.csv --ref my_comparison/ref.csv \
  --zscore sd --level meandiff --pcorr add1 --gamma 2 --formula 1 \
  --perms 5000 --seed 1 --out result.csv
```

`ref.csv` has three columns, `idx,level_meandiff,level_fc`, where `idx` is the
1-based node index from `data/network/nodes.csv`.

## Data availability

* In-house HepG2 LC-MS data, sample metadata, annotations and target-ranking outputs:
  deposited under MetaboLights accession **[to be filled on submission]**; the raw and
  processed files are mirrored in the archive linked from [CITATION.cff](CITATION.cff).
* Public benchmark: every record in `data/benchmark/benchmark_v3.csv` names its
  MetaboLights / Metabolomics Workbench accession and its two sample groups, so the
  analysis can be re-run from the listed accessions with `scripts/` plus `main/MetTD.R`.
* The 9,492 x 9,492 distance matrix (720,784,512 bytes) is too large for Git; its
  SHA-256 and a regeneration script are in the repository, and the file itself is in the
  archived release.

## License and citation

MIT - see [LICENSE](LICENSE). Please cite the software release and the paper; see
[CITATION.cff](CITATION.cff). Archived release and DOI: **[reserved - Zenodo DOI to be
issued at acceptance]**.
