# Validating `read_dynare()` against Dynare

`validate.R` runs each `.mod` file through Dynare (in Octave) and through
`read_dynare()` + `solve_dsge()`, then compares first-order impulse
responses of every endogenous variable to every shock over 20 periods.

```sh
apt-get install octave dynare        # Dynare 6.0 on Ubuntu 24.04
Rscript dev/dynare-validation/validate.R                 # bundled models
Rscript dev/dynare-validation/validate.R path/to/*.mod   # any other files
```

Files that use the macro processor must be expanded first with
`dynare model.mod savemacro onlymacro`.

## Results (2026-09-24, Dynare 6.0, Octave 8.4)

| Model | Features exercised | IRFs | Max abs. difference |
|---|---|---|---|
| `inst/examples/rbc.mod` | `k(-1)` timing, `steady_state_model` | 5 | 4.8e-11 |
| `long_leads_lags.mod` | leads/lags of 2 periods, lagged shock | 6 | 2.7e-13 |
| `nk_linear.mod` | `model(linear)`, local variable, tags, 3 shocks | 18 | 2.5e-14 |
| `rbc_predetermined.mod` | `predetermined_variables`, `initval` only | 5 | 5.7e-10 |
| Dynare `example1.mod` | correlated shocks (`var e, u = ...`) | 12 | 2.7e-11 |
| Dynare `example2.mod` | variance statements | 12 | 4.6e-08 |
| Dynare `agtrend.mod` (macro-expanded) | 15 variables | 30 | 2.2e-08 |
| Dynare `bkk.mod` (macro-expanded) | 4-period time to build, `corr` shocks | 46 | 1.6e-06 |

The Dynare example files are GPL-licensed and are not included in this
repository; they were taken from the DynareR package's demos. The small
differences come from dsge's finite-difference Jacobians (Dynare uses
analytic derivatives); `bkk` has the most leads/lags and a root of 0.995,
so the error accumulates the most there (about 2e-5 relative to the
response).
