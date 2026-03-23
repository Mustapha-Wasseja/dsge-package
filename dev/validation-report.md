# Validation Report — dsge 0.1.0

Date: 2026-03-21
Platform: Windows 11, R 4.5.2

## Overview

This report documents the empirical validation of the `dsge` R package
against three test cases: one analytical benchmark and two worked examples
from the Stata 19 DSGE Reference Manual (Release 19).

---

## Example 1: AR(1) Model (Analytical Benchmark)

**Type:** Exact reproduction of analytical solution
**Source:** Textbook baseline / package vignette
**Script:** `inst/examples/validate_ar1.R`

### Model

```
y_t = z_t
z_{t+1} = rho * z_t + e_{t+1},  e ~ N(0, sigma^2)
```

### Data

Synthetic: `set.seed(42)`, n=300, true rho=0.8, true sigma=1.0

### Results

| Check                   | Result | Detail                            |
|-------------------------|--------|-----------------------------------|
| Policy matrix G         | PASS   | G=1 (exact)                       |
| Transition matrix H     | PASS   | H=0.8 (exact)                     |
| Shock matrix M          | PASS   | M=1 (exact)                       |
| Stability               | PASS   | Stable, eigenvalue = 0.8          |
| IRFs (periods 0-10)     | PASS   | 0.8^k (exact to machine epsilon)  |
| Estimation convergence  | PASS   | BFGS converged (code 0)           |
| rho recovery            | PASS   | Estimated 0.745 (true 0.8)        |
| sigma recovery          | PASS   | Estimated 0.978 (true 1.0)        |
| predict()               | PASS   | Correct dimensions (300 rows)     |
| residuals()             | PASS   | Correct dimensions (300 rows)     |
| forecast()              | PASS   | Produces 8-period forecast         |

**Classification:** Exact reproduction for solver, approximate for estimation
(sampling variability with n=300).

---

## Example 2: Linear New Keynesian Model (Stata Intro 3a)

**Type:** Approximate reproduction of documented Stata example
**Source:** Stata 19 DSGE Reference Manual, Introduction 3a, pp. 36-44
**Script:** `inst/examples/validate_nk_linear.R`
**Stata dataset:** `usmacro2.dta` (not available; synthetic data used for estimation)

### Model

```
p  = beta * E[p'] + kappa * x           (Phillips curve)
x  = E[x'] - (r - E[p'] - g)           (IS curve)
r  = psi * p + u                        (Taylor rule)
u' = rhou * u + e_u                     (monetary shock)
g' = rhog * g + e_g                     (demand shock)
```

Observed: p (inflation), r (interest rate). Unobserved: x (output gap).
Fixed: beta = 0.96.

### Stata documented parameter values

| Parameter | Stata estimate | Std. Err. |
|-----------|---------------|-----------|
| kappa     | 0.0849631     | 0.0287693 |
| psi       | 1.943004      | 0.2957869 |
| rhou      | 0.7005483     | 0.0452604 |
| rhog      | 0.9545257     | 0.0186424 |
| sd(e.u)   | 2.318204      | 0.3047434 |
| sd(e.g)   | 0.5689891     | 0.0982975 |

Log-likelihood: -753.57131

### Solver validation (at Stata parameter values)

| Matrix element | R dsge     | Stata      | Abs. diff  |
|----------------|-----------|-----------|------------|
| G[p,u]         | -0.4172519| -0.4172521| 2e-7       |
| G[p,g]         |  0.9678175|  0.9678177| 2e-7       |
| G[x,u]         | -1.6082158| -1.6082160| 2e-7       |
| G[x,g]         |  0.9529206|  0.9529203| 3e-7       |
| G[r,u]         |  0.1892778|  0.1892776| 2e-7       |
| G[r,g]         |  1.8804733|  1.8804740| 7e-7       |

**Max policy matrix difference: 7e-7**

| Matrix element | R dsge     | Stata      | Abs. diff  |
|----------------|-----------|-----------|------------|
| H[u,u]         | 0.7005483 | 0.7005483 | 0          |
| H[g,g]         | 0.9545257 | 0.9545257 | 0          |

**Max transition matrix difference: 0 (exact)**

### Results

| Check                    | Result | Detail                                 |
|--------------------------|--------|----------------------------------------|
| Policy matrix G          | PASS   | Max diff 7e-7 vs Stata                 |
| Transition matrix H      | PASS   | Exact match                            |
| Stability                | PASS   | Stable, eigenvalues 0.955, 0.701       |
| IRFs                     | PASS   | Computed, geometrically decaying       |
| Estimation (synth data)  | PASS   | Converged, rhou within 0.15 tolerance  |
| Forecast                 | PASS   | 8-period forecast produced             |

**Classification:** Exact reproduction of Stata solver output (to 7 decimal
places). Estimation not compared directly because the exact Stata dataset
is not available — synthetic data validates the estimation pipeline works.

---

## Example 3: Financial Frictions Model (Stata Intro 3c)

**Type:** Approximate reproduction of documented Stata example
**Source:** Stata 19 DSGE Reference Manual, Introduction 3c, pp. 55-62
**Script:** `inst/examples/validate_financial_frictions.R`
**Stata dataset:** `usmacro2.dta` (not available)

### Model

```
p  = beta * E[p'] + kappa * x           (Phillips curve)
x  = E[x'] - (i - E[p'] - g)           (IS curve)
i  = chi * r + e                        (Interest rate spread)
r  = psi * p + u                        (Taylor rule)
e' = rhoe * e + e_e                     (Financial shock)
u' = rhou * u + e_u                     (Monetary shock)
g' = rhoz * g + e_g                     (Demand shock)
```

Observed: p, i, r. Unobserved: x. Fixed: beta = 0.96.

### Stata documented parameter values

| Parameter | Stata estimate | Std. Err. |
|-----------|---------------|-----------|
| kappa     | 0.0503257     | 0.0184939 |
| chi       | 0.9067394     | 0.1249768 |
| psi       | 6.332377      | 2.567617  |
| rhoe      | 0.8478222     | 0.0301221 |
| rhou      | 0.815346      | 0.033254  |
| rhoz      | 0.9861866     | 0.0099767 |
| sd(e.e)   | 0.8857071     | 0.1343023 |
| sd(e.u)   | 7.160717      | 2.933797  |
| sd(e.g)   | 0.3911862     | 0.0358352 |

Log-likelihood: -882.13195

### Solver validation (at Stata parameter values)

**Max policy matrix difference: 4e-7** (vs Stata documentation)
**Max transition matrix difference: 0** (exact match)

### Results

| Check                    | Result | Detail                                 |
|--------------------------|--------|----------------------------------------|
| Model parsing (7 eqs)    | PASS   | 4 control + 3 state equations parsed   |
| Policy matrix G          | PASS   | Max diff 4e-7 vs Stata                 |
| Transition matrix H      | PASS   | Exact match                            |
| Stability                | PASS   | Stable, 3 eigenvalues inside unit circle|
| IRFs                     | PASS   | Financial shock IRFs computed           |

**Classification:** Exact reproduction of Stata solver output (to 7 decimal
places). The 3-observed, 4-control, 3-state model validates correct handling
of larger systems.

---

## Summary

| Example              | Solver   | Policy G | Trans H  | Stability | Estimation | IRFs   | Forecast |
|----------------------|----------|----------|----------|-----------|------------|--------|----------|
| AR(1)                | Exact    | Exact    | Exact    | PASS      | PASS       | Exact  | PASS     |
| NK Linear (Intro 3a) | Exact*   | Exact*   | Exact    | PASS      | PASS†      | PASS   | PASS     |
| Fin. Frictions (3c)  | Exact*   | Exact*   | Exact    | PASS      | N/A‡       | PASS   | N/A‡     |

\* Exact to 7 decimal places vs Stata documented output
† Estimation validated on synthetic data (Stata dataset not available)
‡ Not tested (would require exact Stata dataset for comparison)

## Known differences from Stata

1. **Solution method:** We use the method of undetermined coefficients
   (iterative fixed-point). Stata uses QZ (generalized Schur) decomposition.
   Both methods produce mathematically equivalent results for well-conditioned
   linear models; numerical agreement is within 1e-6.

2. **Standard errors:** We use observed information (Hessian) standard errors.
   Stata uses OPG (outer product of gradients) by default. Both are
   asymptotically equivalent but may differ in finite samples.

3. **Scope:** This release covers linear DSGE only. Stata also supports
   nonlinear DSGE (`dsgenl`), Bayesian estimation (`bayes: dsge`), and
   additional features (steady-state computation, model-implied covariances).
   These are planned for future releases.

4. **Control variable ordering:** Our package orders controls as
   (observed, unobserved); Stata orders them by equation appearance.
   This affects row ordering in G but not the numerical values.
