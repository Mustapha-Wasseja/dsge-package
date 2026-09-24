# Validating `read_dynare()` against Dynare

`validate.R` runs each `.mod` file through Dynare (in Octave) and through
`read_dynare()` + `solve_dsge()`, then compares first-order impulse
responses of every endogenous variable to every shock over 20 periods.
`extra/validate_extra.R` compares log-likelihoods, an optimal simple
rule and OccBin simulations.

```sh
apt-get install octave dynare        # Dynare 6.0 on Ubuntu 24.04
Rscript dev/dynare-validation/validate.R                 # bundled models
Rscript dev/dynare-validation/validate.R path/to/*.mod   # any other files
Rscript dev/dynare-validation/extra/validate_extra.R
Rscript dev/dynare-validation/extra/validate_smets_wouters.R path/to/DSGE_mod/Smets_Wouters_2007
Rscript dev/dynare-validation/extra/batch_dsge_mod.R path/to/DSGE_mod results.csv
Rscript dev/dynare-validation/extra/validate_higher_order.R [ORDER] path/to/*.mod
Rscript dev/dynare-validation/extra/validate_estimation.R path/to/DSGE_mod/RBC_baseline
```

## Impulse responses (2026-09-24, Dynare 6.0, Octave 8.4)

| Model | Features exercised | IRFs | Max abs. difference |
|---|---|---|---|
| `inst/examples/rbc.mod` | `k(-1)` timing, `steady_state_model` | 5 | 4.8e-11 |
| `long_leads_lags.mod` | leads/lags of 2 periods, lagged shock | 6 | 2.7e-13 |
| `nk_linear.mod` | `model(linear)`, local variable, tags, 3 shocks | 18 | 2.5e-14 |
| `rbc_predetermined.mod` | `predetermined_variables`, `initval` only | 5 | 5.7e-10 |
| `ss_leads_det.mod` | `STEADY_STATE()`, shock lead, `varexo_det` | 3 | 5.6e-12 |
| `nk_ramsey.mod` | `ramsey_model`, linear-quadratic | 10 | 9.1e-16 |
| `ramsey_growth.mod` | nonlinear Ramsey problem | 3 | 3.0e-09 |
| `nk_discretion.mod` | `discretionary_policy` | 10 | 1.5e-09 |
| Dynare `example1.mod` | correlated shocks (`var e, u = ...`) | 12 | 2.7e-11 |
| Dynare `example2.mod` | variance statements | 12 | 4.6e-08 |
| Dynare `agtrend.mod` | macro processor (imported unexpanded) | 30 | 2.2e-08 |
| Dynare `bkk.mod` | macro processor, 4-period time to build, `corr` | 46 | 1.6e-06 |

## Other checks (`extra/`)

| Check | Dynare | dsge |
|---|---|---|
| Log-likelihood, NK with measurement error on `r` (`nk_me.mod`) | 539.5442 | 539.544231 |
| Log-likelihood, 2 observables / 3 shocks (`nk_fewobs.mod`) | 332.5837 | 332.583721 |
| OSR loss at optimum (`nk_osr.mod`) | 0.0005290500965 | 0.0005290500964 |
| OSR `phi_pi` (`phi_y` at its bound 2) | 2.336216 | 2.336203 |
| OccBin ZLB, one surprise shock (`nk_zlb_occbin.mod`) | 10 periods at the bound | same; paths within 2.1e-13 |
| OccBin ZLB, two surprise shocks | 15 periods at the bound | same; paths within 2.2e-13 |

## Johannes Pfeifer's DSGE_mod collection (`extra/batch_dsge_mod.R`)

All 68 `.mod` files of https://github.com/JohannesPfeifer/DSGE_mod (GPL,
not included). Each file is cut at its first computing command
(`stoch_simul`, `estimation`, ...), so that Dynare and dsge see the same
model and calibration, and a first-order `stoch_simul` is appended (for
`ramsey_policy` / `discretionary_policy` files, the policy command
itself). Shocks without a declared variance get a standard deviation of
0.01 in both programs. dsge imports the cut file with `read_dynare()`
(MATLAB statements, `verbatim` blocks and `_steadystate.m` files run
through its MATLAB interpreter) and compares all first-order IRFs over 20
periods.

Results (2026-09-24, Dynare 6.0, Octave 8.4):

| Outcome | Files |
|---|---|
| IRFs agree to within 1e-6 (relative to the largest response) | 52 |
| IRFs agree to 7e-5 relative (Kiyotaki-Moore 1997: `k + m*kp` is constant, so the states are linearly dependent and the solution is ill-conditioned) | 1 |
| Perfect-foresight models (`simul`, `perfect_foresight_*`), not compared | 9 |
| Dynare fails in Octave (`hpfilter`, `fmincon`, `ksdensity` or `verLessThan` missing; a `.mat` data file not in the repository; Smets-Wouters parameters that only `estimated_params` initialises) | 6 |

`read_dynare()` imports 67 of the 68 files unchanged; the exception loads
a `.mat` file that is not in the repository (Dynare fails on it too). Most
matching models agree to 1e-9 or better; the models the importer
previously could not handle include files whose calibration or steady
state is computed in MATLAB (e.g. `Gali_2010`, `Basu_Bundick_2017`,
`Born_Pfeifer_2020`, `Chari_et_al_2007`, `Ghironi_Melitz_2005`,
`Jermann_Quadrini_2012_NK`), unit-root models (`McCandless_2008_Chapter_9`),
discretionary policy in purely forward-looking models (`Gali_2008/2015`
chapter 5) and Ramsey policy with `STEADY_STATE()`
(`Gali_2015_chapter_6_4`).

## Second- and third-order solutions (`extra/validate_higher_order.R`)

Dynare solves each model at order 2 and 3 and exports its decision rules
(`ghx`, `ghu`, `ghxx`, `ghxu`, `ghuu`, `ghs2` and, at order 3, `ghxxx`,
`ghxxu`, `ghxuu`, `ghuuu`, `ghxss`, `ghuss`). Both decision rules are
evaluated at the same 20 random points (states drawn from the first-order
ergodic distribution, shocks from their distribution); the table gives
the largest difference over all variables and points relative to the
largest deviation from the steady state (models whose shocks have unit
standard deviations have very large deviations).

| Model | Variables | Order 2 | Order 3 |
|---|---|---|---|
| rbc.mod (inst/examples) | 5 | 1.9e-11 | 1.9e-11 |
| Collard_2001_example1 (correlated shocks) | 6 | 1.5e-11 | 1.5e-11 |
| SGU_2004 | 3 | 1.3e-11 | 3.2e-11 |
| SGU_2003 | 12 | 5.7e-10 | 8.9e-10 |
| Hansen_1985 | 9 | 2.7e-11 | 2.7e-11 |
| Sims_2012_RBC | 13 | 1.2e-10 | 1.2e-10 |
| RBC_baseline | 15 | 3.9e-10 | 5.6e-10 |
| RBC_state_dependent_GIRF | 9 | 1.9e-10 | 1.9e-10 |
| RBC_capitalstock_shock | 6 | 4.4e-09 | 5.3e-09 |
| RBC_news_shock_model (news shocks) | 8 | 1.1e-09 | 1.9e-09 |
| McCandless_2008_Chapter_13 (leads of 2 in nonlinear terms) | 14 | 2.3e-10 | 4.0e-09 |
| GarciaCicco_et_al_2010 | 18 | 5.4e-11 | 6.1e-10 |
| Aguiar_Gopinath_2007 | 21 | 3.1e-08 | 7.1e-08 |
| Caldara_et_al_2012 (Epstein-Zin) | 12 | 2.3e-08 | 2.6e-08 |
| Chari_et_al_2007 (`_steadystate.m`) | 13 | 3.6e-08 | 3.3e-08 |
| Ghironi_Melitz_2005 (`_steadystate.m`) | 35 | 2.3e-07 | 2.3e-07 |
| BP2020_CES (`_steadystate.m`, 43 variables) | 43 | 1.1e-09 | 1.1e-09 |
| Andreasen_2012_rare_disasters (134 variables) | 134 | 7.5e-07 | 1.1e-06 |

In three further models (Jermann 1998, Gali 2015 ch. 3 nonlinear, Basu and
Bundick 2017) the shocks have no variance before the first computing
command, so only the steady state and the zero risk correction are
compared; they agree to 2.4e-9 or better. The remaining small differences
come from dsge's numerical first-order Jacobian (Dynare differentiates
analytically); second and third derivatives are symbolic in both.

## Smets & Wouters (2007), end to end (`extra/validate_smets_wouters.R`)

Johannes Pfeifer's replication file (`DSGE_mod/Smets_Wouters_2007`, GPL,
not included) imported unchanged; comparison at the published posterior
mode on the US data (`presample = 4`; Dynare with `lik_init = 1`, i.e.
`LIK_INIT=1` in the environment, unless stated otherwise):

| Check | Dynare | dsge |
|---|---|---|
| IRFs, 40 variables x 7 shocks (280), max abs diff | | 2.5e-12 |
| Log-likelihood | -1714.061158377 | -1714.061158377 |
| Log-prior (36 priors, incl. 7 `inv_gamma_pdf`) | -23.994069948 | -23.994069948 |
| Log-posterior kernel | -1738.055228325 | -1738.055228325 |
| Log-likelihood with the file's own `lik_init = 2` | -1738.513893160 | -1738.513893160 |
| At dsge's posterior median (600-draw `bayes_dsge()` run): log-likelihood | -1471.497498370 | -1471.497498370 |
| At dsge's posterior median: log-prior | -31.952176636 | -31.952176636 |

The short `bayes_dsge()` run (1 chain, 600 draws, 38 minutes) moved to a
region whose posterior kernel is about 235 log points above the published
mode under these settings (`first_obs = 1` on all 230 quarters,
`lik_init = 1`); Dynare evaluates exactly the same values there, so the
two programs define the same posterior.

Dynare prints the log-likelihood to four decimals. The Dynare example
files are GPL-licensed and are not included in this repository; they were
taken from the DynareR package's demos. The small IRF differences come
from dsge's finite-difference Jacobians (Dynare uses analytic
derivatives); `bkk` has the most leads/lags and a root of 0.995, so the
error accumulates the most there (about 2e-5 relative to the response).

## Full Bayesian estimation (`extra/validate_estimation.R`)

RESULTS_EST
