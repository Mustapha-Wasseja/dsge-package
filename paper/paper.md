---
title: 'dsge: An R package for specifying, solving, and estimating Dynamic Stochastic General Equilibrium models'
tags:
  - R
  - econometrics
  - macroeconomics
  - DSGE
  - Bayesian estimation
  - Kalman filter
  - rational expectations
authors:
  - name: Mustapha Wasseja Mohammed
    orcid: 0000-0001-6693-359X
    affiliation: 1
affiliations:
  - name: World Bank Group
    index: 1
date: 8 April 2026
bibliography: paper.bib
---

# Summary

Dynamic Stochastic General Equilibrium (DSGE) models are the workhorse
framework of modern macroeconomics, combining microeconomic foundations
with stochastic shocks to describe the dynamics of output, inflation,
interest rates, and other aggregate variables. They are routinely used
by central banks and international organizations for forecasting,
policy analysis, and the evaluation of structural reforms.

The `dsge` R package provides an end-to-end, native-R workflow for
specifying, solving, and estimating DSGE models. Linear models are
defined through an equation-based formula interface, and nonlinear
models through string-based equations that are automatically
linearized around their deterministic steady state. Solution uses the
method of undetermined coefficients [@klein2000]; second-order
perturbation is also supported [@schmitt2004]. Models can be estimated
by maximum likelihood using the Kalman filter or by Bayesian methods
using an adaptive Random-Walk Metropolis-Hastings sampler with
user-specified priors. The package also provides a broad set of
post-estimation tools: impulse-response functions with delta-method or
posterior credible bands, dynamic forecasting, Kalman smoothing,
historical shock decomposition, local identification diagnostics,
parameter sensitivity analysis, perfect-foresight transition paths,
occasionally binding constraints (OccBin), model-implied covariance
matrices, robust (sandwich) standard errors, prior-versus-posterior
informativeness diagnostics, posterior predictive checks, and marginal
likelihood for model comparison. Results are returned as standard S3
objects that respond to `coef()`, `vcov()`, `logLik()`, `predict()`,
`fitted()`, `residuals()`, `summary()`, and `plot()`.

# Statement of Need

DSGE modeling in open-source software has long been dominated by
Dynare [@adjemian2011], a MATLAB/Octave-based toolbox. While Dynare is
mature and widely used, it depends on a proprietary or
semi-proprietary numerical platform and operates through its own
domain-specific language. For researchers, graduate students, and
practitioners whose primary workflow lives in R --- a common
situation at central banks, international organizations, statistical
agencies, and academic departments --- DSGE analysis has until now
required stepping outside R, managing an external installation, and
transferring data and results between systems.

The only prior R interface to DSGE modeling was DynareR [@mati2024],
which is a wrapper that executes Dynare code inside R Markdown or
Quarto documents. DynareR still depends on a working Dynare and
Octave installation and does not provide native R objects, model
specification, or estimation routines. Outside R, MacroModelling.jl
[@kockerols2023] and DSGE.jl [@frbnydsge2015] offer native
implementations in Julia but require users to leave the R ecosystem.
Within R, the `gEcon` package [@gecon2019] provides a DSGE modeling
language on R-Forge (not CRAN) but does not implement maximum
likelihood or Bayesian estimation. Table 1 summarizes the comparison.

Table 1: Comparison of `dsge` with existing DSGE software.

| Feature                 | `dsge` (R) | DynareR (R)   | Dynare (MATLAB) | MacroModelling.jl | gEcon (R)      |
|-------------------------|------------|---------------|-----------------|-------------------|----------------|
| Language                | R          | R (wrapper)   | MATLAB/Octave   | Julia             | R              |
| External software       | None       | Dynare+Octave | MATLAB/Octave   | None              | None           |
| On CRAN                 | Yes        | Yes           | N/A             | N/A               | No (R-Forge)   |
| Linear DSGE             | Yes        | Via Dynare    | Yes             | Yes               | Yes            |
| Nonlinear DSGE          | Yes        | Via Dynare    | Yes             | Yes               | Yes            |
| Maximum likelihood      | Yes        | Via Dynare    | Yes             | Yes               | No             |
| Bayesian estimation     | Yes        | Via Dynare    | Yes             | Yes               | No             |
| 2nd-order perturbation  | Yes        | Via Dynare    | Yes             | Yes               | No             |
| Occasionally binding    | Yes        | Via Dynare    | Yes             | Yes               | No             |
| R formula interface     | Yes        | No            | No              | No                | No             |

`dsge` closes this gap by providing the first comprehensive,
CRAN-available, native-R implementation of a DSGE toolkit. It targets
academic researchers in macroeconomics and applied econometrics;
graduate students learning rational-expectations and Bayesian methods;
and central-bank and policy-institution modelers who already work in R
and want a reproducible, open-source environment for DSGE analysis.

# Key Features

- **Formula-based model specification** for linear models using
  `obs()`, `unobs()`, `state()`, and a `lead()` (or `E()`) expectation
  operator that mirrors R's familiar syntax.
- **String-based specification** for nonlinear models via
  `dsgenl_model()`, with automatic Newton-Raphson solution of the
  deterministic steady state and automatic linearization.
- **Solution** by the method of undetermined coefficients [@klein2000],
  with Blanchard-Kahn stability diagnostics and support for
  second-order perturbation [@schmitt2004].
- **Maximum likelihood estimation** via the Kalman filter with
  numerical BFGS optimization, Hessian-based standard errors, and a
  sandwich (Huber-White) robust variance estimator.
- **Bayesian estimation** via an adaptive Random-Walk
  Metropolis-Hastings sampler supporting normal, beta, gamma, uniform,
  and inverse-gamma priors; multiple chains; and automatic proposal
  covariance tuning.
- **Convergence diagnostics** including effective sample size,
  Gelman-Rubin $\hat{R}$, Monte Carlo standard errors, Geweke tests,
  and prior-versus-posterior informativeness diagnostics.
- **Impulse responses** with delta-method confidence bands for ML fits
  and posterior credible bands for Bayesian fits.
- **Kalman smoothing**, **historical shock decomposition**, and
  **forecasting** with uncertainty quantification.
- **Identification diagnostics** based on the Jacobian of the
  autocovariance-to-parameter map, and **parameter sensitivity
  analysis** for the log-likelihood, IRFs, and steady state.
- **Perfect-foresight transition paths** and **occasionally binding
  constraints** (OccBin) for deterministic policy experiments such as
  zero lower bounds on interest rates.
- **Marginal likelihood** and **posterior predictive checks** for
  Bayesian model comparison.
- **Native R interface**: all results are S3 objects that respond to
  `coef()`, `vcov()`, `logLik()`, `predict()`, `fitted()`,
  `residuals()`, `summary()`, and `plot()`.

# Example Usage

Install the package from CRAN and load it.

```r
install.packages("dsge")
library(dsge)
```

Specify a simple three-equation New Keynesian model. The model has
inflation (`p`), an unobserved output gap (`x`), and a nominal
interest rate (`r`), driven by a monetary shock (`u`) and a demand
shock (`g`).

```r
nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),      # Phillips curve
  unobs(x ~ lead(x) - (r - lead(p) - g)),     # IS curve
  obs(r   ~ psi * p + u),                      # Taylor rule
  state(u ~ rhou * u),                         # Monetary shock
  state(g ~ rhog * g),                         # Demand shock
  fixed = list(beta = 0.99),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)
```

Estimate by maximum likelihood on a data frame `macro` with columns
`p` and `r`.

```r
fit <- estimate(nk, data = macro)
summary(fit)
```

Produce impulse responses, a forecast, Kalman-smoothed states, a
historical shock decomposition, and robust standard errors using R's
standard interface.

```r
plot(irf(fit, periods = 20))
plot(forecast(fit, horizon = 12))
plot(smooth_states(fit))
plot(shock_decomposition(fit))
robust_vcov(fit)
```

Estimate the same model by Bayesian methods with user-specified
priors, then inspect the posterior and plot posterior impulse
responses with credible bands.

```r
my_priors <- list(
  kappa = prior("beta", shape1 = 30, shape2 = 70),
  psi   = prior("gamma", shape = 184, rate = 122.7),
  rhou  = prior("beta", shape1 = 70, shape2 = 20),
  rhog  = prior("beta", shape1 = 70, shape2 = 20)
)

fit_bayes <- bayes_dsge(nk, data = macro,
                        priors = my_priors,
                        chains = 2, iter = 10000, warmup = 5000)
summary(fit_bayes)
plot(fit_bayes, type = "trace")
plot(fit_bayes, type = "prior_posterior")
plot(fit_bayes, type = "irf", periods = 20)
```

Nonlinear models are supported through `dsgenl_model()` with automatic
linearization around the deterministic steady state.

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
```

For more complete examples, including Bayesian estimation of a New
Keynesian model on U.S. macroeconomic data from FRED, deterministic
transition paths, and simulations under a zero lower bound on the
interest rate, see the package vignette (`vignette("introduction",
package = "dsge")`) and the example scripts under
`system.file("examples", package = "dsge")`.

# Acknowledgements

The views expressed in this paper are those of the author and do not
necessarily reflect the views of the World Bank Group, its Executive
Directors, or the countries they represent. The author thanks the R
Core Team and the CRAN maintainers for their ongoing stewardship of
the R ecosystem.

# References
