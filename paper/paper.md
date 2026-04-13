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
user-specified priors. Post-estimation tools include impulse-response
functions, dynamic forecasting, Kalman smoothing, historical shock
decomposition, local identification diagnostics, parameter sensitivity
analysis, perfect-foresight transition paths, occasionally binding
constraints, model-implied covariance matrices, robust standard
errors, posterior predictive checks, and marginal likelihood for model
comparison.

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

The `dsge` package closes this gap by providing the first
comprehensive, CRAN-available, native-R implementation of a DSGE
toolkit. It targets academic researchers in macroeconomics and applied
econometrics; graduate students learning rational-expectations and
Bayesian methods; and central-bank and policy-institution modelers
who already work in R and want a reproducible, open-source environment
for DSGE analysis.

# State of the Field

Existing DSGE software falls into three broad categories. First, full
standalone toolboxes running on MATLAB or Octave, most prominently
Dynare [@adjemian2011], which remains the reference implementation.
Second, native implementations in Julia, including MacroModelling.jl
[@kockerols2023] and DSGE.jl [@frbnydsge2015] from the Federal
Reserve Bank of New York, which take advantage of Julia's numerical
performance but require users to work outside the R ecosystem. Third,
R-based wrappers such as DynareR [@mati2024], which call Dynare from
R but do not themselves implement model specification, solution, or
estimation.

Within R, the `gEcon` package [@gecon2019] provides a DSGE modeling
language and supports calibration and impulse responses, but is
distributed only on R-Forge (not CRAN) and does not implement maximum
likelihood or Bayesian estimation.

: Comparison of `dsge` with existing DSGE software. \label{comparison}

| Feature                 | `dsge` (R) | DynareR (R)   | Dynare (MATLAB) | MacroModelling.jl | gEcon (R)      |
|-------------------------|------------|---------------|-----------------|-------------------|----------------|
| Language                | R          | R (wrapper)   | MATLAB/Octave   | Julia             | R              |
| External software       | None       | Dynare+Octave | MATLAB/Octave   | None              | None           |
| On CRAN                 | Yes        | Yes           | N/A             | N/A               | No (R-Forge)   |
| ML estimation           | Yes        | Via Dynare    | Yes             | Yes               | No             |
| Bayesian estimation     | Yes        | Via Dynare    | Yes             | Yes               | No             |
| 2nd-order perturbation  | Yes        | Via Dynare    | Yes             | Yes               | No             |
| Occasionally binding    | Yes        | Via Dynare    | Yes             | Yes               | No             |
| R formula interface     | Yes        | No            | No              | No                | No             |

The `dsge` package is the first comprehensive, CRAN-available,
native-R DSGE toolkit. Its main differentiators are: a formula-based
model specification that mirrors R's familiar syntax; no external
numerical software dependencies; deep integration with R's S3
statistical model interface; and a broad feature set spanning both
maximum-likelihood and Bayesian estimation, identification
diagnostics, perfect-foresight paths, occasionally binding
constraints, and second-order perturbation.

# Software Design

The package solves the canonical linear rational-expectations
state-space representation

$$y_t = G\,x_t, \qquad x_{t+1} = H\,x_t + M\,\varepsilon_{t+1},$$

where $y_t$ collects the control (jump) variables, $x_t$ collects the
predetermined states, and $\varepsilon_t$ is a vector of structural
shocks. The policy matrix $G$ and transition matrix $H$ are obtained
via the method of undetermined coefficients of @klein2000, which uses
the generalized Schur (QZ) decomposition to enforce the
Blanchard-Kahn saddle-path condition.

Nonlinear models are specified using string-based equations. Their
deterministic steady state is computed by a Newton-Raphson solver, and
the model is automatically linearized by numerical differentiation.
Second-order perturbation following @schmitt2004 is also available.

The log-likelihood is evaluated by the Kalman filter using the
standard prediction-error decomposition. Maximum-likelihood estimates
are obtained via BFGS optimization. Bayesian estimation uses an
adaptive RWMH sampler with automatic proposal covariance tuning,
multiple chains, and convergence monitoring via effective sample size,
the Gelman-Rubin $\hat{R}$ statistic, and Monte Carlo standard errors.

The user-facing API follows R conventions. A three-equation New
Keynesian model is specified as:

```r
library(dsge)
nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),
  unobs(x ~ lead(x) - (r - lead(p) - g)),
  obs(r   ~ psi * p + u),
  state(u ~ rhou * u),
  state(g ~ rhog * g),
  fixed = list(beta = 0.99),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)
fit <- estimate(nk, data = macro)
summary(fit)
plot(irf(fit, periods = 20))
plot(shock_decomposition(fit))
```

All results are S3 objects that respond to standard generics
(`coef()`, `vcov()`, `logLik()`, `predict()`, `residuals()`,
`summary()`, `plot()`), so that applied users can work with familiar
R idioms.

# Research Impact Statement

The package has been validated on a sequence of increasingly complex
models, culminating in a full 37-equation specification in the style
of @smets2007, featuring external habit formation, investment
adjustment costs, capacity utilization, Calvo price and wage setting
with indexation, price and wage dispersion, a Taylor rule with
interest-rate smoothing, and seven structural shocks. In
simulated-data experiments, all structural parameters were recovered
inside 95% posterior credible intervals, indicating that the package
can handle realistic, policy-scale DSGE models. The package ships with
a test suite of 570 unit and integration tests covering model
specification, solution, estimation, smoothing, identification,
sensitivity, perfect foresight, occasionally binding constraints, and
Bayesian diagnostics.

By providing a self-contained DSGE toolkit on CRAN, `dsge` lowers the
entry cost for researchers and students who wish to work with DSGE
models entirely in R. The most immediate applications include
reproducible applied macroeconomic research; graduate teaching of
rational-expectations and Bayesian econometric methods; and rapid
prototyping of policy models at central banks, ministries of finance,
and international organizations.

# AI Usage Disclosure

No generative AI tools were used in the creation of the software,
its documentation, or this paper.

# Acknowledgements

The views expressed in this paper are those of the author and do not
necessarily reflect the views of the World Bank Group, its Executive
Directors, or the countries they represent. The author thanks the R
Core Team and the CRAN maintainers for their ongoing stewardship of
the R ecosystem.

# References
