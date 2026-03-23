# Bayesian RBC with Capital Accumulation — Walkthrough

## Model

This example estimates a standard Real Business Cycle (RBC) model with
endogenous capital accumulation. The model is specified in nonlinear form
and estimated via Bayesian methods using first-order perturbation around
the parameter-specific deterministic steady state.

### Equations

| Equation | Role |
|----------|------|
| `1/C = beta / C(+1) * (alpha * exp(Z) * K(+1)^(alpha-1) + 1 - delta)` | Euler equation (intertemporal consumption choice) |
| `K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K` | Capital accumulation (output minus consumption plus undepreciated capital) |
| `Z(+1) = rho * Z` | Technology process (AR(1) in log-productivity) |

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| C | Observed control | Consumption |
| K | Endogenous state | Capital stock (predetermined, not directly observed) |
| Z | Exogenous state | Log-technology shock |

### Parameters

| Parameter | Value | Status |
|-----------|-------|--------|
| alpha | 0.33 | Fixed (capital share) |
| beta | 0.99 | Fixed (quarterly discount factor) |
| delta | 0.025 | Fixed (quarterly depreciation rate) |
| rho | — | **Estimated** (technology persistence) |
| sd(e.Z) | — | **Estimated** (technology shock size) |

### Steady State

The deterministic steady state has an analytical solution:

- Z_ss = 0 (log-technology at its mean)
- K_ss = (alpha * beta / (1 - beta * (1 - delta)))^(1/(1 - alpha)) ≈ 28.35
- C_ss = K_ss^alpha - delta * K_ss ≈ 2.31

The numerical steady-state solver matches these values exactly.

## Data Generation

Data is simulated from the linearized model at the true parameters:

- **True rho:** 0.9 (highly persistent technology process)
- **True sd(e.Z):** 0.007 (standard TFP shock calibration)
- **Sample size:** 300 observations (full version), 200 (demo)
- **Seed:** 42

The simulation uses the policy matrix G and transition matrix H from the
first-order perturbation solution. Data is in levels (C_ss + deviation).

## Prior Specification

| Parameter | Prior | Mean | SD | Rationale |
|-----------|-------|------|-----|-----------|
| rho | beta(10, 2) | 0.833 | 0.103 | Favors persistence; excludes unit root |
| sd(e.Z) | inv_gamma(0.01, 0.01) | (diffuse) | — | Standard diffuse prior for shock scale |

## Estimation

The Bayesian estimator uses adaptive Random-Walk Metropolis-Hastings.
At each MCMC draw, the nonlinear pipeline is:

1. Propose new (rho, sd) in unconstrained space
2. Transform to natural space with Jacobian correction
3. Evaluate log-prior
4. Solve steady state at candidate rho (Newton-Raphson)
5. Linearize (first-order Taylor expansion around steady state)
6. Solve linear RE system (Klein undetermined coefficients)
7. Evaluate Kalman filter likelihood (data centered at steady-state C)
8. Accept/reject via Metropolis ratio

### MCMC Settings (Full Version)

| Setting | Value |
|---------|-------|
| Chains | 2 |
| Iterations | 5,000 per chain |
| Warmup | 2,500 per chain |
| Post-warmup draws | 2,500 per chain (5,000 total) |
| Proposal | Adaptive RWMH (Hessian-informed) |

### Runtime

~9 minutes on a standard laptop (full version).
~2.5 minutes for the lighter demo version.

## Results (Full Version, seed = 42)

### Posterior Summary

| Parameter | True | Posterior Mean | Posterior SD | 95% CI |
|-----------|------|---------------|-------------|--------|
| rho | 0.900 | 0.890 | 0.034 | [0.816, 0.952] |
| sd(e.Z) | 0.007 | 0.0074 | 0.0016 | [0.0045, 0.0108] |

**Both true values fall within their 95% credible intervals.**

### Convergence Diagnostics

| Metric | rho | sd(e.Z) |
|--------|-----|---------|
| R-hat | 1.003 | 1.003 |
| ESS | 594 | 579 |
| Acceptance rate | 25.4% / 24.7% | — |
| Solve failures | 0 | — |

- R-hat < 1.01 for both parameters → strong evidence of convergence
- ESS > 500 → adequate effective samples
- Acceptance near 25% → well-tuned proposal
- Zero solve failures → the steady-state solver is robust across the
  posterior distribution

## Impulse-Response Functions

A one-standard-deviation positive technology shock produces:

### At Impact (period 0)

| Variable | Median Response | 95% Band | Interpretation |
|----------|----------------|----------|----------------|
| C | +0.0032 | [0.0030, 0.0036] | Consumption rises immediately |
| K | 0.0000 | [0.0000, 0.0000] | Capital is predetermined (correct) |
| Z | +0.0073 | [0.0046, 0.0108] | Technology jumps by one SD |

### At Period 5

| Variable | Median Response | 95% Band |
|----------|----------------|----------|
| C | +0.0052 | [0.0046, 0.0059] |
| K | +0.069 | [0.044, 0.093] |

### Economic Interpretation

1. **Consumption rises at impact** — higher productivity raises output,
   and households consume part of the windfall immediately.

2. **Capital is zero at impact** — capital is a predetermined (endogenous
   state) variable. It cannot jump. This is the correct RBC response.

3. **Capital accumulates gradually** — the technology shock raises output,
   and part of the extra output is invested, building the capital stock
   over time. The capital response is hump-shaped, peaking after several
   periods.

4. **Technology mean-reverts** — with rho = 0.9, the shock decays at
   10% per period, returning to zero over ~30 periods.

5. **Consumption response is persistent** — because the capital stock
   builds up, the consumption response actually increases after impact
   before eventually decaying. This is the standard RBC hump-shape.

All of these responses are consistent with standard RBC theory.

## Scripts

| Script | Purpose | Runtime |
|--------|---------|---------|
| `inst/examples/bayesian_rbc_capital_full.R` | Full validation (5,000 iter) | ~9 min |
| `inst/examples/bayesian_rbc_capital_demo.R` | Quick demo (2,000 iter) | ~2.5 min |

### Running the Examples

```r
library(dsge)

# Quick demo (interactive, ~2.5 min)
source(system.file("examples", "bayesian_rbc_capital_demo.R", package = "dsge"))

# Full validation (detailed output, ~9 min)
source(system.file("examples", "bayesian_rbc_capital_full.R", package = "dsge"))
```

## Limitations

- Uses first-order perturbation. Certainty equivalence holds, so
  risk premia and precautionary savings motives are absent.
- Only one observed variable (consumption). With more observables
  (e.g., output, investment, hours), identification of additional
  parameters (alpha, beta, delta) would be possible.
- The technology process is linear in logs (`Z(+1) = rho * Z`).
  The nonlinearity comes from the production function and Euler equation.
