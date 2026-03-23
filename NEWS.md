# dsge 0.4.0

## Bayesian nonlinear DSGE estimation

* `bayes_dsge()` now accepts `dsgenl_model` objects for Bayesian estimation
  of nonlinear DSGE models via first-order perturbation.
* For each MCMC draw, the model's steady state is re-solved and the system
  re-linearized at the candidate parameters, ensuring correct parameter-
  dependent dynamics.
* Failed steady-state solves, linearization failures, and Blanchard-Kahn
  violations are handled gracefully — proposals are rejected without
  crashing the sampler. Failure counts are reported in the output.
* Data is automatically centered around the parameter-specific steady
  state for nonlinear models.
* Posterior IRFs work for nonlinear models via the same `irf()` interface.
* Example scripts: `bayesian_rbc_nonlinear.R` and
  `bayesian_nk_nonlinear.R` in `inst/examples/`.

# dsge 0.3.0

## Bayesian estimation (linear)

* New `bayes_dsge()` function for Bayesian estimation of linear DSGE models
  via adaptive Random-Walk Metropolis-Hastings (RWMH).
* New `prior()` constructor for specifying prior distributions: `normal`,
  `beta`, `gamma`, `uniform`, and `inv_gamma`.
* Multiple chains with dispersed starting points from prior draws.
* Adaptive proposal covariance during warmup targeting ~25% acceptance rate.
* Parameter transformations (log, logit) with Jacobian corrections for
  constrained parameters.
* MCMC diagnostics: effective sample size (ESS), R-hat (split), MCSE,
  and acceptance rates.
* Posterior summary: mean, SD, median, 95% credible intervals.
* `print()`, `summary()`, `coef()`, and `plot()` methods for `dsge_bayes`
  objects. Plot supports trace plots and posterior density with prior overlay.
* Posterior impulse-response functions via `irf()` with pointwise credible
  bands computed from posterior draws.
* Default `inv_gamma(0.01, 0.01)` priors for shock standard deviations.
* `irf()` is now an S3 generic, dispatching to `irf.default` (ML fits and
  solutions) and `irf.dsge_bayes` (Bayesian fits).

## Real-data examples

* Bayesian NK model estimated on real FRED data (GDPDEF inflation and
  federal funds rate, 1955Q1–2015Q4, 244 obs), producing economically
  plausible posteriors and impulse-response functions.
* Extended 3-observable NK model estimated on FRED data (inflation, FFR,
  HP-filtered output gap) with cost-push shock.
* Example scripts: `bayesian_nk_fred.R` and `bayesian_nk_extended.R`
  in `inst/examples/`.

## Bug fixes

* Fixed data-informed initialization of shock standard deviations in
  `estimate()` to avoid convergence to local optima.

# dsge 0.2.0

## Nonlinear DSGE support

* New `dsgenl_model()` constructor for nonlinear DSGE models defined via
  string-based equations with `VAR(+1)` lead notation.
* `steady_state()` generic for computing the deterministic steady state
  via Newton-Raphson with damped line search. Supports user-supplied
  analytical steady-state functions.
* `linearize()` computes first-order Taylor expansion around steady state,
  mapping the Jacobian blocks into the existing canonical structural form.
* `solve_dsge()` now accepts `dsgenl_model` objects: automatically computes
  steady state, linearizes, and solves via the existing Klein solver.
* `estimate()` now accepts `dsgenl_model` objects: re-linearizes at each
  parameter evaluation, subtracts model-implied steady state from data.
* All existing postestimation tools (policy matrix, transition matrix,
  stability, IRFs, forecasts, predict, residuals) work with nonlinear
  models after linearization.
* Extended Klein solver to handle A4 matrix (lead-state terms in control
  equations) for broader model classes.

## Limitations

* Nonlinear support uses first-order perturbation only; second-order and
  higher approximations are not yet implemented.
* Steady-state solving is numerical; convergence depends on initial guesses.

# dsge 0.1.0

Initial CRAN release.

## Model specification

* Formula-based model specification via `obs()`, `unobs()`, and `state()`
  equation wrappers.
* Forward expectations via `lead(x)` (core primitive) and `E(x)` (alias).
* Fixed and free parameter separation with arbitrary nonlinear
  parameter expressions in coefficients.
* Endogenous state variables via `state(, shock = FALSE)`.

## Solution

* Klein (2000) solution via the method of undetermined coefficients.
* Blanchard-Kahn saddle-path stability diagnostics.
* Policy matrix G and transition matrix H extraction.

## Estimation

* Maximum likelihood estimation via Kalman filter log-likelihood.
* Log-parameterization of shock standard deviations for positivity.
* Delta-method standard errors for all estimated quantities.
* Standard R methods: `coef()`, `vcov()`, `logLik()`, `nobs()`,
  `predict()`, `residuals()`, `summary()`.

## Postestimation

* Policy and transition matrix extraction with delta-method standard errors.
* Stability check with eigenvalue classification.
* Impulse-response functions with confidence bands.
* Dynamic multi-step forecasting from filtered states.
* `predict()` with one-step-ahead and filtered/smoothed modes.

## Plotting

* `plot.dsge_irf()` for multi-panel IRF plots with CI bands.
* `plot.dsge_forecast()` for forecast visualization.
