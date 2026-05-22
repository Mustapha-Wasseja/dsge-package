# dsge

Dynamic Stochastic General Equilibrium Models for R.

[![CRAN status](https://www.r-pkg.org/badges/version/dsge)](https://CRAN.R-project.org/package=dsge)

## Overview

The `dsge` package provides a comprehensive framework for specifying,
solving, and estimating DSGE models entirely in R. No external software
(Dynare, MATLAB, Octave) is required.

**Key capabilities:**

- **Linear models** via formula interface (`obs()`, `unobs()`, `state()`)
- **Nonlinear models** via string-based equations with perturbation up to
  third order (`dsgenl_model()`, `solve_dsge(order = 1, 2, 3)`)
- **Maximum likelihood estimation** via the Kalman filter
- **Bayesian estimation** via adaptive RWMH (`bayes_dsge()`) with optional
  parallel chain execution (`n_cores`)
- **Particle filter and PMMH** (`bayes_particle()`) for fully nonlinear
  Bayesian estimation without linearization
- **Bayes factor model comparison** (`bayes_factor()`) with Kass-Raftery
  evidence scales and posterior model probabilities
- **Variance decomposition** (`variance_decomposition()`) -- both
  unconditional steady-state shares and forecast-error variance
  decomposition at user-supplied horizons
- **Optimal Simple Rules** (`osr()`) -- finds the optimal coefficients
  in a user-specified, restricted policy rule (e.g. Taylor rule
  coefficients) given a quadratic welfare loss
- **Conditional forecasts** (`conditional_forecast()`) -- forecasts
  conditional on a pre-specified path for a subset of observables
  (Waggoner-Zha 1999 minimum-norm shocks)
- **IRF matching estimation** (`irf_match()`) -- estimate structural
  parameters by matching DSGE impulse responses to a user-supplied
  target (e.g. VAR-estimated IRFs), Christiano-Eichenbaum-Evans style
- **DSGE-VAR** (`bayes_dsge_var()`, `bayes_dsge_var_mh()`) -- Bayesian
  VAR with DSGE-implied prior (Del Negro & Schorfheide 2004).  The MH
  variant jointly estimates structural parameters, shock SDs, and the
  prior weight lambda; `forecast()` and `conditional_forecast()`
  methods support fan-chart projections and path-conditioned forecasts
  -- Dynare-parity workflow
- **Anticipated / news shocks** in `perfect_foresight_nonlinear()`
- **Ramsey optimal policy** via linear-quadratic regulator
  (`ramsey_policy()`, `welfare_loss()`)
- **Second- and third-order perturbation** with pruned simulation
- **Occasionally binding constraints** (OccBin) for ZLB and other bounds
- **Perfect foresight** deterministic transition paths -- both linearized
  (`perfect_foresight()`) and fully nonlinear via stacked-time Newton
  (`perfect_foresight_nonlinear()`)
- **Kalman smoothing** and historical shock decomposition
- **Local identification diagnostics** and parameter sensitivity analysis
- **Robust (sandwich) standard errors** for ML estimation
- **Posterior predictive checks** and marginal likelihood
- **Model-implied covariance matrices** and prediction tools

## Installation

```r
# Install from CRAN
install.packages("dsge")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("Mustapha-Wasseja/dsge-package")
```

## Quick Start: Maximum Likelihood

```r
library(dsge)

# Define a simple New Keynesian model
nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),       # Phillips curve
  unobs(x ~ lead(x) - (r - lead(p) - g)),      # IS curve
  obs(r   ~ psi * p + u),                       # Taylor rule
  state(u ~ rhou * u),                          # Monetary shock
  state(g ~ rhog * g),                          # Demand shock
  fixed = list(beta = 0.99),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)

# Estimate by maximum likelihood
fit <- estimate(nk, data = your_data)
summary(fit)

# Postestimation
irf(fit, periods = 20) |> plot()               # Impulse responses
forecast(fit, horizon = 12) |> plot()           # Forecasts
smooth_states(fit) |> plot()                    # Kalman-smoothed states
shock_decomposition(fit) |> plot()              # Historical decomposition
check_identification(fit)                       # Identification diagnostics
robust_vcov(fit)                                # Sandwich standard errors
model_covariance(fit)                           # Model-implied moments
```

## Bayesian Estimation

```r
# Specify priors
my_priors <- list(
  kappa = prior("beta", shape1 = 30, shape2 = 70),
  psi   = prior("gamma", shape = 184, rate = 122.7),
  rhou  = prior("beta", shape1 = 70, shape2 = 20),
  rhog  = prior("beta", shape1 = 70, shape2 = 20)
)

# Run MCMC
fit_bayes <- bayes_dsge(nk, data = your_data, priors = my_priors,
                        chains = 2, iter = 10000, warmup = 5000)

# Diagnostics and results
summary(fit_bayes)
plot(fit_bayes, type = "trace")
plot(fit_bayes, type = "prior_posterior")
plot(fit_bayes, type = "irf", periods = 20)
```

Supported priors: `normal`, `beta`, `gamma`, `uniform`, `inv_gamma`.

## Nonlinear DSGE Models

```r
rbc <- dsgenl_model(
  "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
  "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
  "Z(+1) = rho * Z",
  observed = "C",
  endo_state = "K",
  exo_state = "Z",
  fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
  start = list(rho = 0.9),
  ss_guess = c(C = 2, K = 30, Z = 0)
)

sol <- solve_dsge(rbc, params = c(alpha = 0.33, beta = 0.99,
                                   delta = 0.025, rho = 0.9),
                  shock_sd = c(Z = 0.01))
irf(sol, periods = 40) |> plot()
```

## Advanced Features

### Second- and Third-Order Perturbation

```r
sol2 <- solve_dsge(rbc, params = params, shock_sd = sd, order = 2)
simulate_2nd_order(sol2, periods = 200)
irf_2nd_order(sol2, periods = 40)

sol3 <- solve_dsge(rbc, params = params, shock_sd = sd, order = 3)
simulate_3rd_order(sol3, periods = 200)
```

### Occasionally Binding Constraints

```r
# ZLB constraint on the interest rate
obc <- simulate_occbin(sol,
  constraints = list("r >= 0"),
  shocks = list(g = -0.05),
  horizon = 40)
plot(obc)
```

### Perfect Foresight Paths

```r
# Linearized perfect foresight (fast, small shocks)
pf <- perfect_foresight(sol,
  shocks = list(Z = c(-0.05, -0.03, -0.01)),
  horizon = 60)
plot(pf)

# Fully nonlinear perfect foresight via stacked-time Newton
# (recommended for large shocks where nonlinearities matter)
pf_nl <- perfect_foresight_nonlinear(rbc,
  params   = c(rho = 0.9),
  shock_sd = c(Z = 0.01),
  shocks   = list(Z = 0.10),
  horizon  = 40)
plot(pf_nl)
```

### Particle Filter and PMMH

```r
# Bootstrap particle filter likelihood for nonlinear models
ll <- particle_filter_loglik(sol2, data = your_data, n_particles = 1000)

# Particle Marginal Metropolis-Hastings -- fully nonlinear Bayesian
fit_pmmh <- bayes_particle(rbc, data = your_data, priors = my_priors,
                           n_particles = 500, chains = 2, iter = 5000)
```

### Ramsey Optimal Policy

```r
# Quadratic welfare loss on inflation and output gap
rp <- ramsey_policy(sol,
  Q_xx = diag(c(p = 1, x = 0.5)),
  Q_yy = diag(c(r = 0.1)))
welfare_loss(sol, rp$F)   # evaluate welfare under the optimal rule
```

### Bayes Factor Model Comparison

```r
# Compare two Bayesian fits
bf <- bayes_factor(fit_bayes_A, fit_bayes_B,
                   prior_odds = c(0.5, 0.5))
print(bf)   # log Bayes factor + Kass-Raftery evidence label
```

### Variance Decomposition

```r
# Unconditional decomposition (long-run shares)
vd <- variance_decomposition(sol)
print(vd)
plot(vd)

# Forecast-error variance decomposition at multiple horizons
fevd <- variance_decomposition(sol, horizon = c(1, 4, 8, 20))
plot(fevd)
```

### Optimal Simple Rules

```r
# Optimal Taylor-rule coefficient
res <- osr(nk,
  params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
  shock_sd    = c(e.u = 1.0, e.g = 0.5),
  osr_params  = c(psi = 1.5),
  welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.1)),
  lower = 1.01, upper = 5.0)
print(res)
```

### Conditional Forecasts

```r
# Hold the policy rate at 0.5 for the next 4 periods
cf <- conditional_forecast(fit, horizon = 12,
  condition = list(r = c(0.5, 0.5, 0.5, 0.5, rep(NA, 8))))
plot(cf)
```

### IRF Matching Estimation

```r
# Match the DSGE IRF to an externally estimated target
est <- irf_match(nk,
  params_start   = c(kappa = 0.15, psi = 2.0, rhou = 0.6, rhog = 0.8),
  shock_sd_start = c(e.u = 0.8, e.g = 0.7),
  target         = target_irf_dataframe)
```

### DSGE-VAR

```r
# (a) Conditional-on-(theta,lambda) Bayesian VAR with DSGE prior
fit_dv <- bayes_dsge_var(sol, data = your_data,
                         p = 4, lambda = 1.0, n_draws = 1000)
print(fit_dv)
# Compare different lambdas by marginal likelihood
sapply(c(0.5, 1, 2, 5),
       function(l) bayes_dsge_var(sol, your_data, p=4, lambda=l,
                                  n_draws=200)$log_marg_lik)

# (b) Joint MH estimation of (theta_DSGE, sigma, lambda) -- Dynare-parity
priors <- list(kappa = prior("beta",  shape1 = 2, shape2 = 8),
               psi   = prior("normal", mean = 1.5, sd = 0.25),
               rhou  = prior("beta",  shape1 = 2, shape2 = 2),
               rhog  = prior("beta",  shape1 = 8, shape2 = 2))
fit_mh <- bayes_dsge_var_mh(nk, data = your_data, priors = priors,
                            p = 4, chains = 2, iter = 2000)
print(fit_mh)

# Unconditional + conditional forecasts from the DSGE-VAR posterior
fc_unc  <- forecast(fit_mh, horizon = 12)
fc_cond <- conditional_forecast(fit_mh, horizon = 12,
            condition = list(r = c(0.5, 0.5, 0.5, 0.5, rep(NA, 8))))
plot(fc_unc); plot(fc_cond)
```

### Anticipated / News Shocks (Perfect Foresight)

```r
# Nonlinear PF correctly anticipates a future TFP shock
pf_news <- perfect_foresight_nonlinear(rbc,
  params   = c(rho = 0.9),
  shock_sd = c(Z = 0.01),
  shocks   = list(Z = c(0, 0, 0.05)),  # announced at t=1, arrives at t=3
  horizon  = 40)
plot(pf_news)
```

### Parallel MCMC Chains

```r
# Run chains in parallel across cores (PSOCK on Windows, fork on POSIX)
fit_bayes <- bayes_dsge(nk, data = your_data, priors = my_priors,
                        chains = 4, iter = 10000, warmup = 5000,
                        n_cores = 4)
```

## Feature Comparison

| Feature | dsge (R) | DynareR (R) | Dynare (MATLAB) |
|---------|----------|-------------|-----------------|
| Native R implementation | Yes | No (wrapper) | No |
| External software required | None | Dynare + Octave | MATLAB/Octave |
| CRAN package | Yes | Yes | N/A |
| Linear DSGE | Yes | Via Dynare | Yes |
| Nonlinear DSGE | Yes | Via Dynare | Yes |
| ML estimation | Yes | Via Dynare | Yes |
| Bayesian estimation (RWMH) | Yes | Via Dynare | Yes |
| Particle filter / PMMH | Yes | Via Dynare | Yes |
| Parallel MCMC chains | Yes | Via Dynare | Yes |
| 2nd-order perturbation | Yes | Via Dynare | Yes |
| 3rd-order perturbation | Yes | Via Dynare | Yes |
| OccBin / ZLB | Yes | Via Dynare | Yes |
| Linear perfect foresight | Yes | Via Dynare | Yes |
| Nonlinear perfect foresight (LBJ) | Yes | Via Dynare | Yes |
| Ramsey optimal policy | Yes | Via Dynare | Yes |
| Bayes factor model comparison | Yes | No | Partial |
| R model interface (coef, vcov, plot) | Yes | No | No |
| Formula-based specification | Yes | No | No |

## Documentation

- **CRAN:** https://CRAN.R-project.org/package=dsge
- **Vignette:** `vignette("introduction", package = "dsge")`
- **Examples:** `system.file("examples", package = "dsge")`

## References

- Andrieu, C., Doucet, A. and Holenstein, R. (2010). "Particle Markov
  chain Monte Carlo methods." *Journal of the Royal Statistical Society:
  Series B*, 72(3), 269-342.
- Gordon, N. J., Salmond, D. J. and Smith, A. F. M. (1993). "Novel approach
  to nonlinear/non-Gaussian Bayesian state estimation." *IEE Proceedings F*,
  140(2), 107-113.
- Juillard, M., Laxton, D., McAdam, P. and Pioro, H. (1998). "An algorithm
  competition: First-order iterations versus Newton-based techniques."
  *Journal of Economic Dynamics and Control*, 22, 1291-1318.
- Kass, R. E. and Raftery, A. E. (1995). "Bayes factors." *Journal of the
  American Statistical Association*, 90(430), 773-795.
- Klein, P. (2000). "Using the generalized Schur form to solve a
  multivariate linear rational expectations model." *Journal of Economic
  Dynamics and Control*, 24(10), 1405-1423.
- Schmitt-Grohe, S. and Uribe, M. (2004). "Solving dynamic general
  equilibrium models using a second-order approximation to the policy
  function." *Journal of Economic Dynamics and Control*, 28(4), 755-775.

## License

MIT
