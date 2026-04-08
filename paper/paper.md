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
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: World Bank Group
    index: 1
date: 8 April 2026
bibliography: paper.bib
---

# Summary

Dynamic Stochastic General Equilibrium (DSGE) models are the standard
framework in modern macroeconomics for combining microeconomic foundations
with stochastic shocks to describe aggregate fluctuations. These models
are routinely used by central banks and international organizations for
policy analysis, forecasting, and the evaluation of structural reforms.

The `dsge` R package provides an end-to-end workflow for specifying,
solving, and estimating DSGE models natively in R. It supports both
linear models, through an equation-based formula interface, and nonlinear
models, through string-based equations with first-order and second-order
perturbation around the deterministic steady state. Models can be
estimated by maximum likelihood using the Kalman filter or by Bayesian
methods using an adaptive Random-Walk Metropolis-Hastings sampler with
user-specified priors. The package also includes a wide range of
post-estimation tools: impulse-response functions with delta-method or
posterior credible bands, dynamic forecasting, Kalman smoothing,
historical shock decomposition, local identification diagnostics,
parameter sensitivity analysis, perfect-foresight transition paths,
occasionally binding constraints via the OccBin algorithm, model-implied
covariance matrices, robust (sandwich) standard errors, prior-versus-
posterior informativeness comparisons, posterior predictive checks, and
marginal likelihood computation for model comparison.

The package is implemented entirely in base R and requires no external
numerical software. It integrates naturally with R's statistical model
interface through standard S3 methods (`coef`, `vcov`, `logLik`,
`predict`, `fitted`, `residuals`, `summary`, and `plot`).

# Statement of Need

DSGE modeling in open-source ecosystems has long been dominated by
Dynare [@adjemian2011], a MATLAB/Octave-based toolbox. While Dynare is
powerful and widely used, it requires a proprietary or semi-proprietary
numerical platform to run and operates through its own domain-specific
language. For researchers, graduate students, and practitioners whose
primary workflow lives in R --- a common situation at central banks,
international organizations, and statistical agencies --- DSGE analysis
has until now required stepping outside R, managing an external
installation, and transferring data and results between systems.

The only existing R interface to DSGE modeling before `dsge` was the
DynareR package [@mati2024], which is a wrapper that executes Dynare
code from within R Markdown or Quarto documents. DynareR still depends
on a working Dynare and Octave installation and does not provide native
R objects, model specification, or estimation routines. There was, until
the release of `dsge`, no comprehensive, self-contained DSGE package
available on the Comprehensive R Archive Network (CRAN).

The `dsge` package closes this gap. It enables users to (i) specify
models in natural R syntax without leaving the R environment, (ii)
estimate them by maximum likelihood or Bayesian methods using the
package's own Kalman filter and MCMC sampler, and (iii) perform the
full range of post-estimation analysis --- impulse responses, shock
decompositions, identification checks, and so on --- using standard R
interfaces. The target audience includes academic researchers working
in macroeconomics, monetary policy, and applied econometrics; graduate
students learning DSGE methods; central-bank and policy-institution
modelers who already work in R; and teachers who want a reproducible,
open-source environment for the classroom.

# State of the Field

Existing DSGE software falls into three broad categories. First, full
standalone toolboxes running on MATLAB or Octave, most prominently
Dynare [@adjemian2011]. Dynare remains the reference implementation
and supports a wide feature set, but depends on a proprietary or
semi-proprietary numerical platform and uses its own domain-specific
language. Second, native implementations in Julia, including
MacroModelling.jl [@kockerols2023] and DSGE.jl [@frbnydsge2015] from
the Federal Reserve Bank of New York. These take advantage of Julia's
numerical performance but require users to work outside the R
ecosystem. Third, R-based wrappers such as DynareR [@mati2024], which
call Dynare from R but do not themselves implement model
specification, solution, or estimation.

Within R, prior work is limited. The `gEcon` package [@gecon2019]
provides a DSGE modeling language and supports calibration and impulse
responses, but is distributed only on R-Forge (not CRAN), does not
offer a native formula interface, and does not implement maximum-
likelihood or Bayesian estimation. General state-space packages such
as `KFAS` implement Kalman filtering and smoothing for user-supplied
state-space models but are not designed for DSGE specification,
rational-expectations solution, or the specific diagnostic workflows
that DSGE practice requires.

The `dsge` package is the first comprehensive, CRAN-available, native-R
implementation of a DSGE toolkit. Its main differentiators are: a
formula-based model specification that mirrors R's familiar syntax; no
external numerical software dependencies (the only hard import is
`numDeriv`); deep integration with R's S3 statistical model interface;
and a broad feature set that includes both maximum-likelihood and
Bayesian estimation, perfect-foresight paths, occasionally binding
constraints, identification and sensitivity diagnostics, robust
standard errors, and second-order perturbation.

# Software Design

The package solves the canonical linear rational-expectations
state-space representation
$$
y_t = G\,x_t, \qquad x_{t+1} = H\,x_t + M\,\varepsilon_{t+1},
$$
where $y_t$ collects the control (jump) variables, $x_t$ collects the
predetermined states, and $\varepsilon_t$ is a vector of structural
shocks. The policy matrix $G$ and transition matrix $H$ are obtained
via the method of undetermined coefficients of @klein2000, which uses
the generalized Schur (QZ) decomposition to separate stable from
unstable eigenvalues and enforce the Blanchard-Kahn condition.

Nonlinear models are specified using string-based equations, their
deterministic steady state is computed by a Newton-Raphson solver, and
the model is automatically linearized by numerical differentiation.
Second-order perturbation following @schmitt2004 is also available,
providing pruned simulation and generalized impulse responses that
capture asymmetric effects.

The log-likelihood of the linearized state-space model is evaluated by
the Kalman filter using the standard prediction-error decomposition.
Maximum-likelihood estimates are obtained by numerical optimization
(BFGS), with shock standard deviations estimated in log space to
enforce positivity. Bayesian estimation uses an adaptive Random-Walk
Metropolis-Hastings sampler with automatic proposal covariance tuning,
multiple chains, and prior-posterior Jacobian corrections. Convergence
is monitored with effective sample size, the Gelman-Rubin $\hat{R}$
statistic, Geweke diagnostics, and Monte Carlo standard errors.

The user-facing API follows R conventions. A linear New Keynesian model
can be specified as follows:

```r
library(dsge)

nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),       # Phillips curve
  unobs(x ~ lead(x) - (r - lead(p) - g)),      # IS curve
  obs(r   ~ psi * p + u),                       # Taylor rule
  state(u ~ rhou * u),                          # Monetary shock
  state(g ~ rhog * g),                          # Demand shock
  fixed = list(beta = 0.99),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)

fit <- estimate(nk, data = macro)
summary(fit)
plot(irf(fit, periods = 20))
plot(shock_decomposition(fit))
```

The same model can be estimated by Bayesian methods by providing a
`priors` argument to `bayes_dsge()`. Result objects carry standard S3
classes (`dsge_fit`, `dsge_bayes`) that respond to generic methods, so
that applied users can apply familiar R idioms without learning a
package-specific interface.

# Research Impact Statement

The `dsge` package has been validated on a sequence of increasingly
complex models, culminating in a full 37-equation specification in the
style of @smets2007, featuring external habit formation, investment
adjustment costs, capacity utilization, Calvo price and wage setting
with indexation, explicit price and wage dispersion, a Taylor rule
with interest-rate smoothing, and seven structural shocks. In
simulated-data experiments, all structural parameters were recovered
inside 95% posterior credible intervals, indicating that the package
is capable of handling realistic, policy-scale DSGE models rather than
only textbook examples. The package also ships with a test suite of
570 unit and integration tests covering model specification, solution,
estimation, smoothing, identification, sensitivity, perfect foresight,
occasionally binding constraints, and Bayesian diagnostics.

By providing a self-contained DSGE toolkit on CRAN, `dsge` lowers the
entry cost for researchers and students who wish to work with DSGE
models entirely in R. The most immediate applications include
reproducible applied macroeconomic research; graduate teaching of
rational-expectations and Bayesian econometric methods; and rapid
prototyping of small- to medium-scale policy models at central banks,
ministries of finance, and international organizations.

# AI Usage Disclosure

Generative AI tools were used only for light editorial assistance in
structuring sentences and improving the readability of written text.
All mathematical content, model specifications, solution and estimation
algorithms, code, validation designs, and research decisions are the
author's own work.

# Acknowledgements

The views expressed in this paper are those of the author and do not
necessarily reflect the views of the World Bank Group, its Executive
Directors, or the countries they represent. No external funding was
received for the development of this package. The author thanks the R
Core Team and the CRAN maintainers for their ongoing stewardship of the
R ecosystem.

# References
