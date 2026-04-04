# dsge 1.1.0

## New features

### Bayes factor model comparison

* New `bayes_factor()` for pairwise Bayesian model comparison.  Accepts
  two or more `dsge_bayes` or `dsge_marginal_likelihood` objects, computes
  log Bayes factors, posterior model probabilities, and Kass-Raftery (1995)
  evidence labels.  Custom prior model probabilities supported via
  `prior_odds`.

### Parallel MCMC chains

* `bayes_dsge()` gains an `n_cores` argument.  When `n_cores > 1` chains
  run simultaneously using `parallel::mclapply` (POSIX) or a PSOCK cluster
  (Windows), with automatic fallback to sequential execution.

### Third-order perturbation

* `solve_dsge()` now accepts `order = 3` for nonlinear (`dsgenl_model`)
  models, computing cubic state coefficients (`g_xxx`, `h_xxx`),
  state-sigma^2 cross terms (`g_xss`, `h_xss`), and sigma^3 corrections
  (`g_sss`, `h_sss`) following Schmitt-Grohe and Uribe (2004).
* New `simulate_3rd_order()` for pruned simulation of third-order
  approximations (Andreasen, Fernandez-Villaverde and Rubio-Ramirez, 2018).

### Bootstrap particle filter

* New `particle_filter()`: sequential importance resampling with systematic
  resampling and numerically stable log-sum-exp weights.
* New `particle_filter_loglik()`: convenience wrapper accepting a
  `dsge_solution` object.
* New `bayes_particle()`: Particle Marginal Metropolis-Hastings (PMMH;
  Andrieu, Doucet and Holenstein, 2010) for fully nonlinear Bayesian
  estimation without any linearization.  Returns a `dsge_particle` object
  that inherits from `dsge_bayes`, so all existing diagnostics apply.

### Ramsey optimal policy

* New `ramsey_policy()`: computes the welfare-maximising linear feedback
  rule by solving the discrete-time algebraic Riccati equation (DARE) via
  value-function iteration.  Accepts quadratic welfare weights on states
  (`Q_xx`), controls (`Q_yy`), and cross terms (`Q_xy`).
* New `welfare_loss()`: evaluates expected welfare loss under an arbitrary
  (possibly non-optimal) feedback policy for comparison with the Ramsey
  solution.

# dsge 1.0.0

## Major release: complete DSGE estimation and analysis toolkit

### Model-implied covariance reporting (v0.10.0)

* New `model_covariance()` for unconditional covariance and correlation
  matrices of model observables and controls.
* Supports ML fits, Bayesian fits, and raw solution objects.

### Prediction and reporting tools (v0.10.0)

* New `prediction_interval()` for one-step-ahead prediction bands using
  Kalman filter innovation variance.
* New `prediction_accuracy()` for RMSE, MAE, and mean bias statistics.
* Enhanced `fitted()` method for ML-estimated models.

### Robust standard errors (v0.10.0)

* New `robust_vcov()` for sandwich (Huber-White) variance-covariance
  estimation. Provides standard errors robust to model misspecification.

### Bayesian diagnostics and model comparison (v0.10.0)

* New `posterior_predictive()` for posterior predictive checks with
  variance and autocorrelation statistics.
* New `marginal_likelihood()` via modified harmonic mean estimator
  for Bayesian model comparison.
* New `geweke_test()` for Geweke (1992) convergence diagnostic.
* New `mcmc_diagnostics()` for comprehensive MCMC health summary.

# dsge 0.9.0

## Occasionally binding constraints

* New `simulate_occbin()` for piecewise-linear simulation under
  inequality constraints (e.g., zero lower bound).
* New `obc_constraint()` helper for constraint specification.
* Shadow-shock algorithm with iterative regime detection.
* Plot method with constraint-binding period shading.

# dsge 0.8.0

## Second-order perturbation

* `solve_dsge()` now accepts `order = 2` for second-order approximation
  of nonlinear models.
* New `simulate_2nd_order()` for pruned second-order simulation.
* New `irf_2nd_order()` for generalised impulse-response functions
  with asymmetric shock effects.
* Risk/precautionary corrections to the steady state (sigma-correction).

# dsge 0.7.0

## Perfect foresight / deterministic transition paths

* New `perfect_foresight()` for deterministic simulation of transition
  paths after temporary or permanent shocks.
* Supports both linear and linearized nonlinear models.
* Plot method with deviation-from-steady-state and level views.

# dsge 0.6.0

## Identification and sensitivity analysis

* New `check_identification()` for local identification diagnostics
  via Jacobian SVD of the autocovariance mapping.
* New `parameter_sensitivity()` for sensitivity of likelihood, IRFs,
  steady state, and policy matrices to parameter perturbations.
* New `prior_posterior_update()` for Bayesian informativeness diagnostics
  comparing posterior concentration to prior width.

# dsge 0.5.0

## Kalman smoother and shock decomposition

* New `smooth_states()` for Rauch-Tung-Striebel Kalman smoothing.
* New `smooth_shocks()` for extraction of smoothed structural shocks.
* New `shock_decomposition()` for historical decomposition of observed
  variables into individual shock contributions.
* Plot methods for smoothed states and shock decomposition.

# dsge 0.4.0

## Bayesian nonlinear DSGE estimation

* `bayes_dsge()` now accepts `dsgenl_model` objects for Bayesian estimation
  of nonlinear DSGE models via first-order perturbation.
* For each MCMC draw, the model's steady state is re-solved and the system
  re-linearized at the candidate parameters, ensuring correct parameter-
  dependent dynamics.
* Failed steady-state solves, linearization failures, and Blanchard-Kahn
  violations are handled gracefully --- proposals are rejected without
  crashing the sampler. Failure counts are reported in the output.
* Data is automatically centered around the parameter-specific steady
  state for nonlinear models.
* Posterior IRFs work for nonlinear models via the same `irf()` interface.

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
* Posterior impulse-response functions via `irf()` with pointwise credible
  bands computed from posterior draws.

# dsge 0.2.0

## Nonlinear DSGE support

* New `dsgenl_model()` constructor for nonlinear DSGE models defined via
  string-based equations with `VAR(+1)` lead notation.
* `steady_state()` generic for computing the deterministic steady state
  via Newton-Raphson with damped line search.
* `linearize()` computes first-order Taylor expansion around steady state.
* `solve_dsge()` and `estimate()` now accept `dsgenl_model` objects.

# dsge 0.1.0

Initial release.

* Formula-based linear model specification.
* Klein (2000) solution via method of undetermined coefficients.
* Maximum likelihood estimation via Kalman filter.
* Policy/transition matrix extraction with delta-method standard errors.
* Impulse-response functions, forecasting, stability diagnostics.
