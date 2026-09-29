# Reproduction materials: three-group variance-ratio tests for heteroscedasticity

Code and per-replication results for the paper

> El-Kotory, A. and Ebrahim, E. K. *Three-group variance-ratio tests for heteroscedasticity in linear regression: an exact null distribution and an outlier-resistant version.*

The tests themselves are in the R package **KOTORY** (on CRAN; source included here as a Git submodule pinned to the version used for the paper).

## Layout

```
simulations/          every script, in the order it was run, with its stored results
KOTORY/               the R package (Git submodule, version 0.1.0)
paper_stat_papers/    figures/ and tables/: the scripts write the paper's figures and tables here
```

The scripts use paths relative to `simulations/` (`../KOTORY`, `../paper_stat_papers/figures`), so run them from inside that folder:

```bash
git clone --recursive https://github.com/ebrahimkhaled/KOTORY-reproduction.git
cd KOTORY-reproduction/simulations
Rscript 29_factorial_summary.R
```

R packages used: `robustbase`, `lmtest`, `skedastic`, `SuppDists` (one check only), `carData`, `pkgload` and `devtools` (to load KOTORY from source), and `parallel`.

## Where each result comes from

Every simulation stores the p-value of every test in every replication, so each number can be recomputed from the stored files without rerunning the simulation. Summary scripts (`17`, `29`, `34`) read only stored results and run in seconds.

| Paper item | Simulation | Summary / output |
|---|---|---|
| Theorem 1 (exact law of KaH-III): Kolmogorov–Smirnov checks and Hartley's tabulated points | `01_kah3_vs_hartley.R`, `02_kah3_exact_check.R` | the csv files they write |
| Hartley's maximum F distribution | `fmax_exact.R` (sourced by the other scripts) | |
| Table 1, Fig. 1: effective degrees of freedom of the LTS scale | `03_lts_effective_df.R` (the table shipped in KOTORY is built by `KOTORY/data-raw/nu_table.R`) | `14_nu_star_theory.R` (the limit r(α)), `16_paper_figures.R`, `17_paper_tables.R` |
| Fig. 2: null distributions | `04_validate_robust.R`, `07_figures.R` (`robust_full/`) | `16_paper_figures.R` |
| Fig. 3: size of all tests | `08_alpha_check.R`, `11_all_tests.R` (`all_tests/`) | `22_size_heatmap.R` |
| Figs. 4 and 5: factorial comparison | `28_factorial.R` (`factorial/`) | `29_factorial_summary.R`, `30_factorial_fig.R`, `31_scenario_grid.R`, `32_winning_cells.R` |
| Fig. 6: the usage rule against the best other test in each setting | `28_factorial.R` | `36_rule_vs_best_fig.R` |
| Table 2: heavy-tailed errors | `33_heavy_tails.R` (`heavy_tails/`) | `34_heavy_tails_summary.R` |
| Table 3: power | `11_all_tests.R` | `17_paper_tables.R` |
| Fig. 7: power curves | `21_power_curves.R` (`power_curves/`) | `23_power_curves_fig.R` |
| Table 4: size-adjusted power | `11_all_tests.R` | `20_size_adjusted_power.R`, `17_paper_tables.R` |
| Table 5: the robust Goldfeld–Quandt test | `rgq_alih_ong.R`, `13_rgq.R` (`rgq/`) | `17_paper_tables.R` |
| Table 6: benchmark data sets | `15_real_data.R` | `17_paper_tables.R` |
| Section 4, the wages example (SLID data from `carData`) | `35_wages_example.R` | printed |
| BAMSET computed from its definition (Appendix) | `18_bamset_check.R`, `18b_bamset_own.R`, `19_bamset_rerun.R` | |

Scripts `05`, `06`, `10`, `12` and `24`–`27` are exploratory steps (checks, variance shapes, further scenarios) kept for completeness; the paper's results do not depend on them.

## Random-number safety

Two functions of the `skedastic` package reset the random-number generator when called, which would make every later replication of a simulation loop reuse the same data. Every simulation here calls foreign code through a wrapper that saves and restores the generator state, and stores a fingerprint of each data set; a script stops if two replications share a fingerprint. Runs made before this guard existed were discarded and are not included. From step 09 only the least trimmed squares table (`competitors/k0.rds`) is used, which does not call `skedastic`.

## Licence

Code: MIT. Stored results: CC BY 4.0. See `LICENSE`.
