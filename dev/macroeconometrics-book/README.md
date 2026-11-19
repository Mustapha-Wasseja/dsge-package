# Paccagnini, *Macroeconometrics*, Chapter 6 (DSGE models) with dsge

Companion code: <https://github.com/AlessiaPaccagnini/Macroeconometrics_book>
(De Gruyter, forthcoming). Chapter 6 is the book's only DSGE chapter. Its
model-based examples are run here with dsge, in R, without Dynare or MATLAB.

| Book item | Book code | What it computes | Here |
|---|---|---|---|
| Table 6.1 | `Codes_Figures/chapter6_dsge_companion.R`, Section 1 | HP-filtered business-cycle moments of US data | not a DSGE model; nothing for dsge to do |
| Figure 6.1 | same file, Section 2 | RBC impulse responses, hand-derived log-linearisation solved by Blanchard-Kahn | `ch6_figure_6_1_rbc_irf.R` |
| Figure 6.2 | same file, Section 3 | prior vs posterior for six SW (2007) parameters, drawn as normal approximations around the published mode and 90% interval | `ch6_smets_wouters_estimation.R` (from an actual estimation) |
| Empirical example | `Codes_Empirical_Examples/sw_table1.mod` (+ 3 included `.mod` files), `usmodel_data.mat`, `table1.m` | Bayesian estimation of Smets-Wouters (2007) in Dynare; Tables 1A and 1B | `ch6_smets_wouters_estimation.R` |

The book's "Extra Code" folder only links to an external MATLAB toolbox
(DSGE-SVt), so there is nothing there to run.

## Running

The book's files are not copied here (the book repository has no licence
file). Point `MACRO_BOOK_DIR` at a clone, or leave it unset and the scripts
clone the repository into a temporary directory.

```bash
export MACRO_BOOK_DIR=/path/to/Macroeconometrics_book
Rscript dev/macroeconometrics-book/ch6_figure_6_1_rbc_irf.R
Rscript dev/macroeconometrics-book/ch6_smets_wouters_estimation.R
```

The estimation needs `R.matlab` (to read `usmodel_data.mat`). Its posterior
mode and draws are cached in `output/*.rds` (not committed), so a re-run only
redraws the tables and figures; set `SW_REFIT_MODE=1` / `SW_REFIT_MCMC=1` to
recompute. The chain length is set with `SW_CHAINS`, `SW_ITER`, `SW_WARMUP`,
`SW_CORES` and `SW_SEED`.

<!-- RESULTS -->
