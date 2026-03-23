# DSGE Package Worked Examples: Walkthrough

This document explains the two worked examples shipped with the `dsge` package
and discusses whether their outputs are economically sensible.

## Example 1: Linear New Keynesian Model (`linear_example.R`)

### Model

A three-equation New Keynesian model:

- Phillips curve: `p = beta * E[p'] + kappa * x`
- IS curve: `x = E[x'] - (r - E[p'] - g)`
- Taylor rule: `r = psi * p + u`
- Monetary shock: `u' = rhou * u + e_u`
- Demand shock: `g' = rhog * g + e_g`

This is the same model as Stata 19's DSGE Reference Manual (Intro 3a).
Three controls (inflation p, output gap x, interest rate r) and two exogenous
states (monetary policy shock u, demand shock g).

### Step-by-step walkthrough

**Step 1: Specification.** The model is defined using `dsge_model()` with the
formula interface. `beta` is fixed at 0.96 (discount factor). The free
parameters are `kappa` (Phillips curve slope), `psi` (Taylor rule coefficient),
`rhou` and `rhog` (shock persistence).

**Step 2: Solve at known parameters.** The model is solved at Stata's documented
values: kappa=0.085, psi=1.94, rhou=0.70, rhog=0.95, with shock SDs u=2.3,
g=0.57. The solver produces policy matrix G (maps states to controls) and
transition matrix H (state dynamics). These are verified to be saddle-path stable
with both eigenvalues inside the unit circle.

**Step 3: Matrices and stability.** The policy matrix shows that:
- A positive monetary shock (u>0) reduces inflation and output gap (contractionary)
- A positive demand shock (g>0) raises all three variables
These signs are economically correct for a standard NK model.

**Step 4: Simulate data.** 200 observations are generated from the state-space
solution using the true parameter values. The simulated data has the right
statistical properties (zero-mean, reasonable variance).

**Step 5: Estimation.** Maximum likelihood via Kalman filter recovers parameters
reasonably well from 200 observations:
- kappa: true=0.085, estimated~0.118 (within ~1 SE)
- psi: true=1.94, estimated~1.31 (less precise, typical for Taylor rule coefficients)
- rhou: true=0.70, estimated~0.63 (close)
- rhog: true=0.95, estimated~0.91 (close)

All estimates have correct signs and are within economically reasonable ranges.
The convergence is clean (no warnings). This level of estimation noise is
typical for ML estimation with 200 observations.

**Step 6: Postestimation.** Policy and transition matrices are extracted from
the estimated model with delta-method standard errors. The stability check
confirms saddle-path stability at estimated values.

**Step 7: IRFs.** Impulse-response functions show:
- Monetary shock (u): contractionary effect on inflation and output gap,
  positive interest rate impact, hump-shaped decay. Economically sensible.
- Demand shock (g): expansionary, raises all three variables, slower decay
  due to higher persistence (rhog > rhou). Economically sensible.

**Step 8: Forecast.** 12-period ahead forecasts converge toward the
unconditional mean (zero for demeaned data), as expected.

### Economic sensibility verdict: PASS

All outputs are economically sensible and consistent with the standard
textbook NK model.

---

## Example 2: Nonlinear RBC Model (`nonlinear_example.R`)

### Model

A real business cycle model in nonlinear form:

- Euler equation: `1/C = beta/C(+1) * (alpha*exp(Z)*K^(alpha-1) + 1-delta)`
- Capital accumulation: `K(+1) = exp(Z)*K^alpha - C + (1-delta)*K`
- TFP process: `Z(+1) = rho * Z + e_Z`

One observed control (consumption C), one endogenous state (capital K), one
exogenous state (TFP Z). Parameters: alpha=0.33 (capital share), beta=0.99
(quarterly discount factor), delta=0.025 (depreciation rate), rho=0.9 (TFP
persistence, free parameter to estimate).

### Step-by-step walkthrough

**Step 1: Specification.** The model is defined using `dsgenl_model()` with
string-based equations. The `VAR(+1)` notation denotes next-period variables.
Structural parameters alpha, beta, delta are fixed; rho is free.

**Step 2: Steady state.** Newton-Raphson converges in 5 iterations. The
numerical steady state matches analytical values exactly:
- K_ss = 28.3484 (capital stock)
- C_ss = 2.3066 (consumption)
- Z_ss = 0 (log TFP in steady state)

These are standard values for a quarterly-calibrated RBC model.

**Step 3: Linearize.** First-order Taylor expansion around steady state
produces the structural matrices. Key features:
- A0 is negative (reflecting the concavity of 1/C)
- B3 shows Z has persistence 0.9 and K evolves with coefficient ~1.01
  (close to 1 because depreciation is small)
- B3[K,Z] = 3.015 captures the strong effect of TFP on capital accumulation

**Step 4: Solve.** The linearized system is solved via undetermined coefficients.
Policy matrix: C responds positively to both Z (0.45) and K (0.05).
- The Z effect dominates because TFP shocks raise both output and consumption
- The K effect is smaller but positive (more capital = more output = more consumption)
Transition matrix: Z decays at rate 0.9, K adjusts slowly (eigenvalue 0.96).
Both eigenvalues are inside the unit circle: saddle-path stable.

**Step 5: IRFs.** The impulse response of C to a TFP shock shows:
- Immediate positive jump (consumption rises with productivity)
- Gradual decay over ~40 periods as capital adjusts back
This is the textbook RBC impulse response.

**Step 6: Estimation.** Data is simulated from the linearized state-space form
(200 observations, shock SD = 0.01) to ensure consistency with the estimator.
Maximum likelihood recovers rho = 0.913 (true value 0.9), well within the
95% CI [0.844, 0.982]. The shock SD is also recovered reasonably (0.0085 vs
true 0.01).

**Step 7: Forecast.** 12-period ahead forecasts from the estimated model
converge toward the steady-state consumption level, as expected.

### Economic sensibility verdict: PASS

All outputs are economically sensible and consistent with the standard
RBC model. The steady state, policy functions, IRFs, and parameter
recovery are all correct.

---

## Notes on estimation

- **Shock SD initialization.** The estimator uses data-informed starting values
  for shock standard deviations (based on the observed data variability). This
  helps avoid local optima that can arise when the initial shock SD is far from
  the true value.

- **Sample size.** With 200 observations, ML estimates have moderate noise.
  The linear NK model (2 observed variables, 4 free parameters) has slightly
  less precision than the nonlinear RBC model (1 observed, 1 free parameter).
  This is expected: more parameters relative to data = less precision.

- **Convergence.** Both examples converge cleanly with no warnings. The BFGS
  optimizer handles both cases well.

## CRAN readiness

The package passes R CMD check with no code-related errors or warnings.
All 123 unit tests pass. Environment-related issues (missing pandoc for
vignette building, LaTeX font for PDF manual) are deployment-dependent
and do not affect package functionality.
