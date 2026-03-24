# CRAN submission comments — dsge 1.0.0

## R CMD check results

0 errors | 0 warnings | 1 note

The single NOTE ("unable to verify current time") is an environment-specific
issue with the local system clock and does not reflect a package problem.
It does not appear on standard CRAN check infrastructure.

Tested on:

* Windows 11, R 4.5.2 (local, R CMD check)

## Package purpose

The dsge package provides tools for specifying, solving, and estimating
dynamic stochastic general equilibrium (DSGE) models by maximum likelihood
and Bayesian methods. It supports linear models (formula interface) and
nonlinear models (first-order perturbation). This is the first major
release (v1.0.0) consolidating the full estimation and postestimation
toolkit.

## Major features in v1.0.0

* Linear DSGE: specification, solution (Klein 2000), ML estimation
* Nonlinear DSGE: string-based equations, Newton-Raphson steady state,
  automatic linearization, ML estimation
* Bayesian estimation: adaptive RWMH for linear and nonlinear models
* Kalman smoother and smoothed structural shocks
* Historical shock decomposition of observed variables
* Local identification diagnostics (Jacobian rank/SVD)
* Parameter sensitivity analysis (likelihood, IRFs, steady state)
* Prior-vs-posterior informativeness comparison
* Perfect foresight / deterministic transition paths
* Second-order perturbation with pruned simulation
* Occasionally binding constraints (piecewise-linear simulation)
* Model-implied covariance/correlation matrices
* Prediction tools: fitted values, intervals, accuracy (RMSE/MAE)
* Robust (sandwich) standard errors for ML estimation
* Posterior predictive checks, marginal likelihood, Geweke diagnostics
* Impulse-response functions with delta-method SEs and Bayesian bands

## Dependencies

Hard dependencies are minimal: stats, graphics, grDevices, numDeriv.
No compiled code. Suggested packages (coda, testthat, knitr, rmarkdown)
are used only for optional MCMC conversion, testing, and vignettes.

## Test suite

570 tests, all passing. Tests cover model specification, solution,
estimation, Kalman filter, Bayesian estimation, smoother, identification,
sensitivity, perfect foresight, second-order perturbation, OBC simulation,
covariance reporting, prediction tools, robust SEs, and Bayesian
diagnostics.

## Downstream reverse dependencies

None. This is a new package.
