# Reproducing every figure and table

| Manuscript item | Script | Inputs | Output |
|---|---|---|---|
| Fig. 2b ROC, Fig. 2c metrics (revised protocol) | `scripts/fig2bc_v3_losoc.py` then `scripts/fig2bc_v3_from_eval.py`, `scripts/plot_fig2bc_pair.py` | `data/benchmark/benchmark_v3.csv`, `data/network/*`, `model_dist.bin` | `fig2bc_v3_evalcsv.json`, `Fig2bc_v3_*.png/pdf` |
| Fig. 2d/e target ranking and precision@N | `scripts/b6_precision.py` | 13-drug ranking tables | `Fig2_precision_at_N.csv` |
| Fig. 3c score-p scatter (metformin 2 h) | `src/mettd_core ... --ref t2h/ref.csv` then `scripts/finalize_b2.py` | `model_dist.bin`, `t2h/ref.csv` | `out/formula1/t2h.csv`, `Fig3c_rerun_scatter.*` |
| Fig. 4d enzyme trend over time | `scripts/fig4d_trend.py` | time-course score tables | `Fig4d_rerun_trend.*` |
| Fig. 4 PCA source coordinates | `scripts/plot_pca_v2.py` | exported PC1/PC2 table | `代谢组PCA_v2.*` |
| Fig. 5 network statistics | `scripts/net_stats.py` | `model_dist.bin`, `nodes.csv`, `enzyme.csv` | printed statistics (9,492 nodes / 17,253 edges / 2,625 components / 2,621 isolated) |
| SI sensitivity figure (Fig. R2) | `scripts/si_sens_degree_figs.py` | `out/formula1/sens_*.csv`, `out/sens_*.csv` | `SI_敏感性分析图.*` |
| SI degree-bias figure (Fig. R3) | `scripts/degree_bias.py` | score tables | `SI_度偏倚检验.*` |
| Fig. R1 pathway-local simulation | `scripts/pathway_perm_figR1.py`, `scripts/b4v2_pathway_final.py` | 10 metabolites x 5 enzymes | ACO2 rank 1 (S=1.46, P=0.037), IDH1 rank 2 (S=1.38, P=0.068) |
| Benchmark composition and evidence tiers | `scripts/build_metadata_v3.py`, `scripts/fill_isa.py` | `data/benchmark/*`, ISA-Tab cache | `benchmark_metadata_v3.csv`, tier counts |
| TCGA survival (HR, CI, PH, adjusted, BH) | `scripts/tcga_adjusted.R` | `tcga_cox_input.csv` | `tcga_survival_adjusted.csv` |

Each script is self-contained, takes its inputs from the paths above, and prints the
numbers that appear in the manuscript, so a reader can diff the printed values against
the article text.
