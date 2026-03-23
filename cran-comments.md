# CRAN submission comments — dsge 0.4.0

## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on:

* Windows 11, R 4.5.2 (local, R CMD check --as-cran)

## Changes since v0.3.0

This release adds Bayesian estimation for nonlinear DSGE models:

* `bayes_dsge()` now accepts `dsgenl_model` objects. For each MCMC draw,
  the steady state is re-solved and the model re-linearized at the
  candidate parameters.
* Failed steady-state solves and linearization failures are handled
  gracefully (proposals are rejected without crashing the sampler).
* Validated on nonlinear RBC and NK models with successful parameter
  recovery.

## Dependencies

Hard dependencies remain minimal: `stats`, `graphics`, `grDevices`,
and `numDeriv`. No compiled code. `coda` is in Suggests only (optional,
for `as.mcmc()` conversion).

## Scope

* Linear DSGE: ML and Bayesian estimation
* Nonlinear DSGE: ML and Bayesian estimation (first-order perturbation)
* Bayesian nonlinear uses first-order approximation around parameter-
  specific steady states; this is not a fully nonlinear filtering solution

## Downstream reverse dependencies

None.
