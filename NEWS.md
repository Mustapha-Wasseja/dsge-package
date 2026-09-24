# dsge 1.2.0

This release collects the additions made between May and September 2026;
the May-August ones were previously listed under 1.1.1.  Dates are taken
from the git history.  Entries within each section are in chronological
order.

## New features

### Nonlinear perfect foresight (2026-05-20)

* New `perfect_foresight_nonlinear()`: stacked-time Newton (LBJ) solver for
  deterministic transition paths of nonlinear DSGE models.  Solves the full
  nonlinear equilibrium conditions simultaneously over the entire horizon,
  giving paths exact to Newton tolerance even for large shocks.  Uses a
  block-bidiagonal Jacobian assembled by numerical finite differences and
  solved by O(T n³) block back-substitution; Armijo backtracking stabilises
  convergence.  Returns the same `dsge_perfect_foresight` class as the
  linearized `perfect_foresight()`, so all existing `plot`, `print`, and
  `summary` methods apply unchanged.  Reference: Juillard et al. (1998),
  *Journal of Economic Dynamics and Control*, 22, 1291–1318.

### Variance decomposition (2026-05-21)

* New `variance_decomposition()` function with methods for
  `dsge_solution`, `dsge_fit`, and `dsge_bayes`.  Computes either:
  - the **unconditional** decomposition (default), giving each shock's
    share of the long-run steady-state variance of every observable
    via per-shock discrete-Lyapunov solves; or
  - the **forecast-error variance decomposition (FEVD)** when a vector
    of integer horizons is supplied, giving each shock's share of the
    h-step-ahead forecast-error variance.
* Returns a `dsge_variance_decomposition` object with both raw
  contributions (in squared units) and percent shares; row sums equal
  100 by construction.
* New `plot.dsge_variance_decomposition()` draws horizontal stacked
  bars for the unconditional case and one stacked-bar panel per
  observable (horizons on the x-axis) for FEVD.  Uses the standard
  package theme palette and styling.

### Closing the most-cited Dynare gaps: five features (2026-05-21)

* **`osr()`** -- Optimal Simple Rules.  Finds the parameters of a
  user-specified, restricted policy rule that minimise an unconditional
  quadratic welfare loss subject to the model's rational-expectations
  equilibrium.  Complements the existing fully-flexible `ramsey_policy()`.
  Implements Dennis (2007).
* **`conditional_forecast()`** -- forecasts conditional on a pre-specified
  path for a subset of observables (e.g. holding the policy rate fixed
  for k periods).  Implements the Waggoner & Zha (1999) minimum-norm
  shock approach.  Returns a `dsge_conditional_forecast` object that
  inherits from `dsge_forecast` so existing plot methods work.
* **`irf_match()`** -- impulse-response matching estimation.  Estimates
  the model's structural parameters and shock SDs by minimising the
  weighted squared distance between model IRFs and a target IRF data
  frame (typically from a VAR).  Follows Christiano, Eichenbaum & Evans
  (1999).
* **`bayes_dsge_var()`** -- Bayesian VAR with DSGE-implied prior
  (Del Negro & Schorfheide 2004).  Combines a Bayesian VAR(p) with prior
  moments centred on the DSGE's second-moment implications, controlled
  by a single hyperparameter lambda.  Returns posterior draws of the
  VAR coefficient matrix and innovation covariance and an approximate
  log marginal likelihood for choosing lambda.
* **Multi-period / "news-shock" paths in perfect foresight**.  Both
  `perfect_foresight()` and `perfect_foresight_nonlinear()` already
  accept vector-valued shock paths; documentation now explicitly calls
  out this capability.  The nonlinear stacked-time solver correctly
  models anticipated shocks (agents adjust at t=1 in expectation of a
  future shock); the linearised version uses the recursive policy and
  treats the path as a sequence of period-by-period surprises.

### Full DSGE-VAR workflow, Dynare parity (2026-05-21)

* **`bayes_dsge_var_mh()`** -- joint Metropolis-Hastings estimation of
  the DSGE structural parameters, shock standard deviations, and the
  DSGE-prior weight lambda jointly.  The VAR coefficients are
  analytically marginalised at every iteration via the closed-form
  Normal-inverse-Wishart conjugate posterior, so only
  (theta_DSGE, sigma, lambda) are sampled.  Reuses the same adaptive
  RWMH sampler as `bayes_dsge()`.  This brings the package to
  feature-parity with Dynare's `estimation(..., dsge_var)` command.
* **`forecast.dsge_dsgevar()`** and **`forecast.dsge_dsgevar_mh()`** --
  unconditional fan-chart forecasts from a DSGE-VAR posterior.
  Iterates the VAR forward for each posterior draw, drawing innovations
  from N(0, Sigma) and aggregating across draws to produce posterior
  summary statistics.
* **`conditional_forecast.dsge_dsgevar()`** and the matching
  `_mh()` method -- DSGE-VAR conditional forecasts on a user-specified
  path for a subset of variables.  Each posterior draw is conditioned
  via a per-period minimum-norm innovation injection (Sigma-metric
  minimum-norm gap closure), parallel to the Dynare
  `Dvars_forecast.m` workflow.

### Derived parameters in dsge_model() (2026-05-22)

* `dsge_model()` gains a new `derived = function(p) list(...)` argument
  that maps primitive parameters to derived ones (analogue of Dynare's
  `# macro` substitutions).  Called by `solve_dsge()` at every solve,
  so derived parameters track the current primitives during estimation.
  This enables encoding rich structural DSGEs (e.g. Smets-Wouters 2007)
  that depend on derived steady-state ratios like `Rk`, `W`, `K_Y`,
  `beta_bar`, etc.

### Closing 9 more Dynare-feature gaps (2026-05-26)

Each is implemented as a standalone function with full documentation
and tests:

* **`calibrated_smoother()`** -- run the Kalman smoother on a
  calibrated (un-estimated) model.  Convenience wrapper around new
  `smooth_states.dsge_solution()`, `smooth_shocks.dsge_solution()`,
  `shock_decomposition.dsge_solution()` methods.
* **`predetermined()`** -- alias for `state(..., shock = FALSE)`,
  matching Dynare's `predetermined_variables` declaration for users
  porting models across implementations.
* **`perfect_foresight_expect_err()`** -- perfect-foresight simulation
  with expectation errors: agents form subjective expectations at each
  period and nature delivers the actual shock path (potentially
  different).  Analogue of Dynare's
  `perfect_foresight_with_expectation_errors_solver`.
* **`extended_path()`** -- Fair-Taylor / Adjemian-Juillard extended-path
  stochastic simulation.  Works on `dsge_solution` (linear) and
  `dsgenl_model` (nonlinear) objects.
* **`endogenous_prior()`** -- Christiano-Trabandt-Walentin endogenous
  prior on model-implied second moments.  Plugs into `bayes_dsge()`
  via a new `endogenous_prior` argument.
* **`global_sensitivity()`** -- global sensitivity analysis with
  Sobol' indices (first-order and total-effect) and Morris elementary
  effects.  Targets include model-implied moments, IRF magnitudes, or
  log-likelihood.  Companion to local sensitivity via
  `parameter_sensitivity()`.
* **`discretionary_policy()`** -- time-consistent (no-commitment)
  optimal policy via the Soederlind / Dennis fixed-point iteration.
  Companion to `ramsey_policy()` (commitment) and `osr()` (restricted
  simple rules).
* **`gmm_estimate()`** and **`smm_estimate()`** -- generalised /
  simulated method-of-moments estimation matching empirical to model
  moments (SDs, variances, covariances, autocorrelations, correlations).
  Supports one-step and two-step (optimal) weighting matrices.
* **`bayes_smc()`** -- tempered Sequential Monte Carlo sampler
  (Chopin 2002; Herbst & Schorfheide 2014).  Robust alternative to
  `bayes_dsge()` (RWMH) and `bayes_particle()` (PMMH), especially for
  multimodal posteriors.  Returns a `dsge_smc` object inheriting from
  `dsge_bayes` so all existing diagnostics and post-estimation
  methods work.

### Four further Dynare-feature gaps closed (2026-08-20)

* **`model_latex()`** -- export a model's equations as LaTeX (the
  analogue of Dynare's `write_latex_dynamic_model`).  Handles both
  linear and nonlinear models, substitutes Greek-letter parameter
  names, puts time subscripts on variables and wraps leads in a
  conditional expectation operator.  Options cover the align/gather
  environments, numbering, `\label` emission and a standalone
  compilable document.  Verified by compiling the generated output
  with `pdflatex`.
* **`kalman_filter_skewed()`** -- likelihood evaluation when the
  structural shocks are skew-normal rather than Gaussian.  The Kalman
  gain remains the optimal linear filter, and the filter additionally
  propagates the third cumulant of the state exactly (the observation
  equation carries no measurement error, so the update is a
  deterministic linear projection).  The predictive density is built
  from univariate skew-normals matched to the exact model-implied
  marginal skewness, and reduces exactly to the Gaussian filter when
  all skewness is zero.
* **`ms_filter()`** -- Markov-switching shock volatility via the Kim
  (1994) filter, with Hamilton regime probabilities, the K^2-to-K
  collapse, the Kim backward smoother, and print/plot methods.  At
  first order the policy functions are certainty-equivalent, so the
  regime rescales only the shock loading and a single solve suffices.
* **`pac_weights()`, `pac_target_loading()`, `pac_simulate()`** --
  FRB/US-style polynomial adjustment cost equations.  Factors the
  Euler equation's characteristic polynomial into stable and unstable
  roots, returns the implied lag coefficients and forward weights
  (normalised so long-run homogeneity holds exactly), and collapses
  the infinite forward sum into closed form when the target follows a
  linear state process.

### Importing Dynare .mod files (2026-09-24)

* New **`read_dynare()`** reads a Dynare `.mod` file (or model code passed
  as text) and translates it into a `dsgenl_model`, together with the
  calibration, shock standard deviations and priors it declares.  No
  Dynare, MATLAB or Octave installation is needed.  `solve_dsge()`,
  `estimate()` and `bayes_dsge()` accept the imported object directly,
  using its calibration, `shocks` block and translated priors by default.
* Supports `var`, `varexo`, `parameters`, `predetermined_variables`,
  `varobs`, parameter assignments, the `model` block (including
  `model(linear)`, equation tags and `#` model-local variables), leads
  and lags of any length, `initval`, `steady_state_model`, the `shocks`
  block (including correlated shocks, orthogonalised by Cholesky
  factorisation in declaration order as in Dynare's impulse responses)
  and `estimated_params` / `estimated_params_init`.
* Dynare timing needs no manual re-timing: lags become auxiliary state
  variables (`x_lag1`, ...), leads beyond one period become auxiliary
  controls (`x_lead1`, ...), and each shock becomes an exogenous state
  holding the current innovation, the same approach Dynare uses
  internally.
* Dynare's mean/standard-deviation priors are converted to dsge's
  parameterisation: exactly for `normal_pdf`, `beta_pdf`, `gamma_pdf`,
  `uniform_pdf` and `inv_gamma2_pdf`, and by moment matching for
  `inv_gamma_pdf`.  Anything not translated (macro directives, shifted
  priors, `weibull_pdf`, deterministic shock paths, other blocks) is
  reported, never silently dropped.
* Validated against Dynare 6.0: first-order impulse responses of 8
  models (134 responses, including Dynare's `example1`/`example2`, the
  macro-expanded `bkk` and `agtrend` examples, and models with
  two-period leads and lags, lagged shocks, correlated shocks and
  `predetermined_variables`) agree to within 2e-6, and to within 1e-7
  for all but `bkk`.  The validation script is in
  `dev/dynare-validation/`.
* New example file `inst/examples/rbc.mod`.

## Improvements

### Publication-ready plot theme (2026-05-21)

* All `plot.dsge_*` methods now use a unified visual theme: consistent
  navy/brick/olive palette, light dotted gridlines, thin solid zero
  reference lines, semi-transparent confidence bands, and sans-serif
  typography.  Centralised in a new internal `R/plot-theme.R`, so the
  look-and-feel of every figure -- IRFs, forecasts, smoothed states and
  fit, historical shock decomposition, perfect-foresight paths,
  occasionally-binding-constraint comparisons, identification
  diagnostics, sensitivity tornadoes, prior-posterior overlays, MCMC
  traces, posterior densities, running means, ACFs and posterior
  predictive checks -- is now uniform and publication-ready.  No API
  changes; all existing user code keeps working.

### Forecast fan chart + history overlay (2026-05-21)

* `forecast.dsge_fit()` now also returns the **forecast standard
  deviation** at each horizon (`sd` column in `forecasts`), computed by
  iterating the state covariance forward analytically.  The result also
  carries the in-sample observed data as `history`.
* `plot.dsge_forecast()` now shows the historical series in grey,
  followed by the point forecast in navy with three nested fan bands at
  the 50/80/95% levels and a vertical separator at the forecast origin.

### Smoothed state uncertainty bands (2026-05-21)

* The Kalman smoother now propagates and stores the smoothed state
  covariances.  `smooth_states()` returns an additional matrix
  `smoothed_states_var` (T x n_states) with the per-period state
  variances.
* `plot.dsge_smoothed(..., type = "states")` automatically draws
  semi-transparent +/- 2 sigma bands around each smoothed state when
  the variance information is available.

### Perfect-foresight comparison overlay (2026-05-21)

* `plot.dsge_perfect_foresight()` gains a new `compare = ...` argument
  that overlays a second perfect-foresight path on the same panels
  (typical use: pass the linearized `perfect_foresight()` result
  alongside a `perfect_foresight_nonlinear()` result to visualise the
  nonlinearity premium for large shocks).

## Bug fixes

### Robust Lyapunov solver (2026-05-22)

* `compute_unconditional_P()` now falls back to the doubling
  algorithm (Smith / Anderson) when the Kronecker-form system
  `(I - H ⊗ H)` is near-singular -- typical for medium-scale DSGEs
  with highly persistent shocks (eigenvalues close to 1, as in
  Smets-Wouters).  Previously the solver returned a fake fallback of
  `diag(1e6)`, which inflated `model_covariance()` and DSGE-VAR
  prior moments by many orders of magnitude.  Affected functions
  (`model_covariance`, `variance_decomposition`, `bayes_dsge_var`,
  `bayes_dsge_var_mh`) now produce sensible moments on highly
  persistent models.

### VAR stability filtering in DSGE-VAR forecasts (2026-05-22)

* `forecast.dsge_dsgevar()`, `forecast.dsge_dsgevar_mh()` and the
  matching conditional-forecast methods now drop any posterior draw
  whose VAR companion matrix has eigenvalues outside the unit disc,
  emitting a message reporting the number of skipped draws.  Prevents
  occasional explosive paths when the DSGE-VAR posterior has long tails.

## Package maintenance

* `.claude/` (local Claude Code settings) added to `.Rbuildignore` so
  it is never bundled into the source tarball.

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
