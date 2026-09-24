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

## Smets & Wouters (2007), end to end (`extra/validate_smets_wouters.R`)

Johannes Pfeifer's replication file (`DSGE_mod/Smets_Wouters_2007`, GPL,
not included) imported unchanged; comparison at the published posterior
mode on the US data (Dynare with `lik_init = 1`, `presample = 4`):

| Check | Dynare | dsge |
|---|---|---|
| IRFs, 40 variables x 7 shocks (280), max abs diff | | 2.5e-12 |
| Log-likelihood | -1714.061158377 | -1714.061158377 |
| Log-prior (36 priors, incl. 7 `inv_gamma_pdf`) | -23.994069948 | -23.994069948 |
| Log-posterior kernel | -1738.055228325 | -1738.055228325 |
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
