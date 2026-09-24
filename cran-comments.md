## Submission summary

This is an update of the dsge package (version 1.2.0).  It adds new
functionality and fixes bugs; no function exported in 1.0.0 has been
removed.  See NEWS.md for details.
Highlights:

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

0 errors | 0 warnings | 1 note

- "checking for future file timestamps ... NOTE: unable to verify
  current time". This comes from the check machine being unable to reach
  an external time server, not from the package.

## Reverse dependencies

TODO before submitting: run `revdepcheck::revdep_check()` (or
`tools::package_dependencies("dsge", reverse = TRUE)`) and record the
result here.
