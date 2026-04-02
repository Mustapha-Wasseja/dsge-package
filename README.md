# dsge

Dynamic Stochastic General Equilibrium Models for R.

[![CRAN status](https://www.r-pkg.org/badges/version/dsge)](https://CRAN.R-project.org/package=dsge)

## Overview

The `dsge` package provides a comprehensive framework for specifying,
solving, and estimating DSGE models entirely in R. No external software
(Dynare, MATLAB, Octave) is required.

**Key capabilities:**

- **Linear models** via formula interface (`obs()`, `unobs()`, `state()`)
- **Nonlinear models** via string-based equations with first-order
  perturbation (`dsgenl_model()`)
- **Maximum likelihood estimation** via the Kalman filter
- **Bayesian estimation** via adaptive RWMH (`bayes_dsge()`)
- **Second-order perturbation** with pruned simulation
- **Occasionally binding constraints** (OccBin) for ZLB and other bounds
- **Perfect foresight** deterministic transition paths
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

### Second-Order Perturbation

```r
sol2 <- solve_dsge(rbc, params = params, shock_sd = sd, order = 2)
simulate_2nd_order(sol2, periods = 200)
irf_2nd_order(sol2, periods = 40)
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
pf <- perfect_foresight(sol,
  shocks = list(Z = c(-0.05, -0.03, -0.01)),
  horizon = 60)
plot(pf)
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
| Bayesian estimation | Yes | Via Dynare | Yes |
| 2nd-order perturbation | Yes | Via Dynare | Yes |
| OccBin / ZLB | Yes | Via Dynare | Yes |
| R model interface (coef, vcov, plot) | Yes | No | No |
| Formula-based specification | Yes | No | No |

## Documentation

- **CRAN:** https://CRAN.R-project.org/package=dsge
- **Vignette:** `vignette("introduction", package = "dsge")`
- **Examples:** `system.file("examples", package = "dsge")`

## References

- Klein, P. (2000). "Using the generalized Schur form to solve a
  multivariate linear rational expectations model." *Journal of Economic
  Dynamics and Control*, 24(10), 1405-1423.

## License

MIT
