## Submission summary

This is an update of the dsge package from 1.0.0 (the version on CRAN)
to 1.2.0; version 1.1.0 was not submitted. It adds new functionality and
fixes bugs; no function exported in 1.0.0 has been removed. See NEWS.md
for details.
Highlights:

- `read_dynare()` imports Dynare `.mod` model files into R, including
  Dynare's macro language, calibration, steady state, shocks,
  measurement errors, priors and estimation settings, and Ramsey,
  discretionary and optimal-simple-rule policy and occasionally binding
  constraints (OccBin). It is written entirely in R: it does not call
  Dynare, MATLAB or Octave and adds no new package dependencies. Results
  were checked against Dynare 6.0 during development; the comparison
  scripts are in the GitHub repository and excluded from the package
  build.
- From 1.1.0: Bayes factor model comparison, parallel MCMC chains,
  third-order perturbation, a bootstrap particle filter with particle
  marginal Metropolis-Hastings, and Ramsey optimal policy.
- Nonlinear perfect foresight via a stacked-time Newton solver.
- Variance decomposition (unconditional and forecast-error).
- Optimal simple rules, discretionary policy, conditional forecasts,
  impulse-response matching, GMM/SMM estimation, DSGE-VAR (with joint
  Metropolis-Hastings and forecasting), a tempered SMC sampler, extended
  path, endogenous priors and global sensitivity analysis.
- A Kalman filter for skew-normal shocks, a Markov-switching volatility
  filter, polynomial adjustment cost equations, and LaTeX export of
  model equations.
- Bug fix: the unconditional state covariance now falls back to the
  doubling algorithm when the Kronecker-form Lyapunov system is
  near-singular, instead of returning a placeholder matrix.

## Test environments

- local: Ubuntu 24.04, R 4.3.3
- GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)

## R CMD check results

`R CMD check --as-cran` with the CRAN incoming checks enabled
(`_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`, spelling via aspell):

0 errors | 0 warnings | 2 notes

- "checking CRAN incoming feasibility ... NOTE: Possibly misspelled words
  in DESCRIPTION: Andrieu, Grohe, Juillard, Kass, Raftery, Schmitt,
  Uribe, al, et". These are author names in the cited references and
  "et al.", and are spelled correctly.
- "checking for future file timestamps ... NOTE: unable to verify
  current time". This comes from the check machine being unable to reach
  an external time server, not from the package.

The DOI for Schmitt-Grohe and Uribe (2004) in DESCRIPTION has been
corrected to 10.1016/S0165-1889(03)00043-5 (the previous one did not
resolve).

## Reverse dependencies

There are currently no reverse dependencies: no CRAN package depends on,
imports, links to or suggests dsge
(`tools::package_dependencies("dsge", reverse = TRUE, which = "all")`,
checked 2026-09-24).
