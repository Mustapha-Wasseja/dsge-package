# dsge

Dynamic Stochastic General Equilibrium Models for R.

## Overview

The `dsge` package provides tools for specifying, solving, and estimating
DSGE models by maximum likelihood and Bayesian methods. It supports:

- **Linear models** via a formula-based interface (`obs()`, `unobs()`, `state()`)
- **Nonlinear models** via string-based equations with first-order
  perturbation around steady state (`dsgenl_model()`)
- **Bayesian estimation** of linear models via adaptive MCMC (`bayes_dsge()`)

Postestimation tools include policy and transition matrix extraction,
stability diagnostics, impulse-response functions with confidence/credible
bands, and forecasting with fan charts.

## Installation

```r
# Install from source
install.packages("path/to/dsge", repos = NULL, type = "source")

# Or using devtools
devtools::install("path/to/dsge")
```

## Quick Start: Maximum Likelihood

```r
library(dsge)

# Define a simple New Keynesian model
nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),
  unobs(x ~ lead(x) - (r - lead(p) - g)),
  obs(r   ~ psi * p + u),
  state(u ~ rhou * u),
  state(g ~ rhog * g),
  fixed = list(beta = 0.96),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)

# Estimate the model
fit <- estimate(nk, data = your_data)
summary(fit)

# Postestimation
policy_matrix(fit)
transition_matrix(fit)
stability(fit)

# Impulse-response functions (with plot)
irfs <- irf(fit, periods = 20)
plot(irfs)

# Forecast (with plot)
fc <- forecast(fit, horizon = 12)
plot(fc)
```

## Bayesian Estimation

Bayesian estimation is available for linear DSGE models via adaptive
Random-Walk Metropolis-Hastings (RWMH).

```r
library(dsge)

# Define the model (same syntax as ML)
nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),
  unobs(x ~ lead(x) - (r - lead(p) - g)),
  obs(r   ~ psi * p + u),
  state(u ~ rhou * u),
  state(g ~ rhog * g),
  start = list(beta = 0.95, kappa = 0.15, psi = 1.5, rhou = 0.7, rhog = 0.9)
)

# Specify priors
my_priors <- list(
  beta  = prior("beta", shape1 = 95, shape2 = 5),
  kappa = prior("beta", shape1 = 30, shape2 = 70),
  psi   = prior("gamma", shape = 184, rate = 122.7),
  rhou  = prior("beta", shape1 = 70, shape2 = 20),
  rhog  = prior("beta", shape1 = 70, shape2 = 20)
  # shock SDs default to inv_gamma(0.01, 0.01)
)

# Run MCMC
fit <- bayes_dsge(nk, data = your_data, priors = my_priors,
                  chains = 2, iter = 10000, warmup = 5000)

# Results
summary(fit)                        # posterior table with ESS, R-hat

# MCMC diagnostics
plot(fit, type = "trace")            # trace plots (all chains)
plot(fit, type = "density")          # posterior density + prior overlay
plot(fit, type = "prior_posterior")   # dedicated prior-vs-posterior
plot(fit, type = "running_mean")     # cumulative mean convergence
plot(fit, type = "acf")              # autocorrelation diagnostics
plot(fit, type = "pairs")            # pairwise scatter + correlations
plot(fit, type = "all")              # combined trace + density panel

# Parameter selection (works with all diagnostic types)
plot(fit, type = "trace", pars = c("kappa", "psi"))

# Posterior IRFs with credible bands
plot(fit, type = "irf", periods = 12)
```

Supported prior distributions: `normal`, `beta`, `gamma`, `uniform`,
`inv_gamma`.

**Plotting:** All plot types support a `pars` argument for parameter
selection (e.g., `pars = c("kappa", "psi")`). Pagination handles any
number of parameters. Forecast plotting is not yet available for
Bayesian fits.

## Nonlinear DSGE Example

```r
library(dsge)

# Define a nonlinear RBC model
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

# Compute steady state
ss <- steady_state(rbc, params = c(alpha = 0.33, beta = 0.99,
                                    delta = 0.025, rho = 0.9))

# Solve (linearize + Klein solver)
sol <- solve_dsge(rbc, params = c(alpha = 0.33, beta = 0.99,
                                   delta = 0.025, rho = 0.9))

# IRFs and postestimation work identically to linear models
irfs <- irf(sol, periods = 40, se = FALSE)
plot(irfs)
```

## Try It Yourself in RStudio

The package ships with self-contained demo scripts. After installing:

```r
library(dsge)

# AR(1) model demo (simplest case)
source(system.file("examples", "demo_ar1.R", package = "dsge"))

# New Keynesian model demo (multivariate)
source(system.file("examples", "demo_nk.R", package = "dsge"))

# Nonlinear RBC model demo (first-order perturbation)
source(system.file("examples", "demo_rbc.R", package = "dsge"))

# Bayesian NK with real FRED data (requires internet; ~10 min)
source(system.file("examples", "bayesian_nk_fred.R", package = "dsge"))

# Bayesian nonlinear NK model (~15-25 min)
source(system.file("examples", "bayesian_nk_nonlinear.R", package = "dsge"))

# Bayesian nonlinear RBC with capital — quick demo (~2.5 min)
source(system.file("examples", "bayesian_rbc_capital_demo.R", package = "dsge"))
```

## What This Package Currently Supports

**Included in v0.4.0:**

- Linear DSGE models with rational expectations
- Nonlinear DSGE models via first-order perturbation (linearization around
  deterministic steady state)
- Bayesian estimation of both linear and nonlinear DSGE models via
  adaptive RWMH
- For nonlinear models, the steady state is re-solved and the model
  re-linearized at each candidate parameter vector
- Prior specification: normal, beta, gamma, uniform, inverse-gamma
- MCMC diagnostics: ESS, R-hat, MCSE, acceptance rates
- Posterior impulse-response functions with credible bands
- Formula-based specification for linear models (`obs()`, `unobs()`, `state()`)
- String-based specification for nonlinear models (`dsgenl_model()`)
- Deterministic steady-state solver (Newton-Raphson or user-supplied)
- Method of undetermined coefficients solver (Klein 2000)
- Maximum likelihood estimation via the Kalman filter
- Delta-method standard errors for all structural parameters
- Impulse-response functions with confidence bands
- Multi-step forecasting with fan charts
- Policy matrix, transition matrix, and Blanchard-Kahn stability diagnostics
- Predictions (fitted values) and residuals

**Not yet supported (planned for future releases):**

- Higher-order perturbation (second-order, third-order)
- Nonlinear filtering (particle filter)
- Model-implied covariance matrices
- Variance decomposition
- Historical shock decomposition

## Validation

The package has been validated on:

- **Synthetic data:** parameter recovery tests for AR(1), NK, and RBC models
  (both linear and nonlinear)
- **Real data:** Bayesian NK estimation on US macroeconomic data from FRED
  (1955Q1–2015Q4), producing economically plausible posteriors and IRFs
- **Nonlinear Bayesian:** parameter recovery on nonlinear RBC and NK models
  with all true parameters within 95% credible intervals
- **Solver accuracy:** policy and transition matrices verified against
  analytical solutions and published benchmark results

## Key Features

- **Formula-based model specification** with `obs()`, `unobs()`, and `state()`
  equation wrappers
- **Forward expectations** via `lead(x)` or `E(x)` operators
- **Klein (2000) solution method** via undetermined coefficients
- **Maximum likelihood** and **Bayesian estimation** via Kalman filter
- **Prior specification** with `prior()` constructor
- **MCMC diagnostics**: ESS, R-hat, MCSE, acceptance rates
- **Postestimation tools**: policy matrix, transition matrix, stability check,
  predictions, impulse-response functions, forecasting
- **Built-in plotting**: IRF plots, forecast fan charts, trace plots,
  posterior density plots, prior-vs-posterior, running mean, ACF, pairs
- **CRAN-compatible** design with minimal dependencies

## Model Syntax

Variables are declared by their role in the model:

| Wrapper | Role | Example |
|---------|------|---------|
| `obs()` | Observed control variable | `obs(p ~ beta * lead(p) + kappa * x)` |
| `unobs()` | Unobserved control variable | `unobs(x ~ lead(x) - r)` |
| `state()` | State variable (with shock) | `state(u ~ rho * u)` |
| `state(, shock = FALSE)` | Endogenous state (no shock) | `state(k ~ delta * k, shock = FALSE)` |

## References

- Klein, P. (2000). "Using the generalized Schur form to solve a multivariate
  linear rational expectations model." *Journal of Economic Dynamics and Control*,
  24(10), 1405-1423.
- Iskrev, N. (2010). "Local identification in DSGE models." *Journal of Monetary
  Economics*, 57(2), 189-202.

## License

MIT
