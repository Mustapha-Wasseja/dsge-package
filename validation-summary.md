# Validation Summary — dsge v0.4.0

## Overview

This document summarizes the validation tests performed on the `dsge`
package across its four estimation modes: ML linear, ML nonlinear,
Bayesian linear, and Bayesian nonlinear.

## Test Suite

167 unit tests passing (testthat). Tests cover:

- Model parsing and specification (linear and nonlinear)
- Klein solver correctness
- Kalman filter likelihood
- ML estimation parameter recovery (AR(1), NK)
- Bayesian estimation parameter recovery (AR(1), NK — linear)
- Bayesian nonlinear estimation parameter recovery
- Graceful handling of steady-state/linearization failures
- Nonlinear posterior IRFs
- Steady-state computation
- Linearization
- IRF computation
- Forecasting
- Policy/transition matrix extraction
- MCMC diagnostics (ESS, R-hat)

## Solver Validation — ML Estimation

### Linear NK Model (Policy Matrices)
- **Result:** Policy matrix G and transition matrix H verified against
  analytical solutions and published benchmark results.
- **Accuracy:** Differences at machine precision (~1e-7).

### Nonlinear NK Model (Policy Matrices)
- **Result:** Linearized policy matrices verified against published
  benchmark results, accounting for the level-deviation vs
  percentage-deviation convention.
- **Script:** `inst/examples/nonlinear_example.R`

## Bayesian Estimation — Linear Models

### AR(1) Model (Synthetic)
- **Result:** True rho = 0.8 recovered within posterior credible interval.
- **Diagnostics:** R-hat < 1.1, adequate ESS.

### NK Model (Real FRED Data, 2 Observables)
- **Data:** GDPDEF inflation + federal funds rate, 1955Q1–2015Q4, 244 obs
- **Result:** Economically plausible posteriors; all IRF signs correct.
- **Diagnostics:** R-hat < 1.01, ESS 573–955, acceptance 24–29%.
- **Script:** `inst/examples/bayesian_nk_fred.R`

### NK Model (Real FRED Data, 3 Observables)
- **Data:** Inflation + FFR + HP-filtered output gap
- **Result:** All results economically plausible; IRF signs correct.
- **Script:** `inst/examples/bayesian_nk_extended.R`

## Bayesian Estimation — Nonlinear Models

### Nonlinear AR(1) (Synthetic, Parameter Recovery)
- **Model:** y = z, z(+1) = z^rho
- **True:** rho = 0.85, sd = 0.50
- **Result:** rho posterior mean = 0.801, 95% CI [0.735, 0.867] — true value in CI.
- **Diagnostics:** R-hat 1.003, ESS 574–617, acceptance 24–26%.
- **Solve failures:** 0 out of ~10,000 evaluations.
- **Script:** `inst/examples/bayesian_rbc_nonlinear.R`

### Nonlinear NK (Synthetic, 5-equation, Parameter Recovery)
- **Model:** Euler, Rotemberg Phillips, Taylor rule, 2 AR(1) shocks
- **True:** phi = 47, psi = 1.9, rhom = 0.7, rhoz = 0.95
- **Result:** All 4 true parameters within 95% credible intervals.
- **Diagnostics:** R-hat < 1.02, ESS 176–207, acceptance 19–24%.
- **Solve failures:** 9 out of ~10,000 evaluations.
- **IRFs:** Monetary shock contractionary (p↓ x↓ r↑), demand shock
  expansionary (p↑ x↑ r↑) — economically correct.
- **Script:** `inst/examples/bayesian_nk_nonlinear.R`

## What Is Not Yet Validated

- Higher-order perturbation (not implemented)
- Nonlinear filtering / particle filter (not implemented)
- Variance decomposition (not implemented)
- Models with more than ~10 parameters (not tested at scale)

## R CMD Check

- 0 errors, 0 warnings, 0 notes (on code)
- Vignette build warnings are environment-dependent (missing pandoc/LaTeX),
  not code issues
