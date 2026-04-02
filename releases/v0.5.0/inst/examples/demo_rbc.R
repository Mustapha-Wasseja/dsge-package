# =============================================================================
# Demo: Nonlinear RBC Model — First-Order Approximation
# =============================================================================
#
# This script demonstrates the nonlinear DSGE workflow on a simple
# Real Business Cycle (RBC) model solved via first-order perturbation.
#
# To run in RStudio:
#   library(dsge)
#   source(system.file("examples", "demo_rbc.R", package = "dsge"))
# =============================================================================

library(dsge)

# --- 1. Specify the nonlinear model ------------------------------------------
# Euler equation:       1/C = beta/C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)
# Capital accumulation: K(+1) = exp(Z) * K^alpha - C + (1-delta)*K
# TFP process:          Z(+1) = rho * Z   (shock enters linearly at first order)
#
# Observed: C (consumption)
# Endogenous state: K (capital)
# Exogenous state: Z (TFP, with shock)

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
print(rbc)

# --- 2. Compute steady state --------------------------------------------------
params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, rho = 0.9)
ss <- steady_state(rbc, params = params)
print(ss)

# --- 3. Solve (linearize + Klein solver) --------------------------------------
sol <- solve_dsge(rbc, params = params, shock_sd = c(Z = 0.01))
cat("\nSolution:\n")
print(sol)

cat("\nPolicy matrix G (maps states to controls):\n")
print(policy_matrix(sol, se = FALSE))

cat("\nTransition matrix H:\n")
print(transition_matrix(sol, se = FALSE))

cat("\nStability:\n")
print(stability(sol))

# --- 4. Impulse-response functions --------------------------------------------
irfs <- irf(sol, periods = 40, se = FALSE)
plot(irfs)

# --- 5. Simulate data and estimate --------------------------------------------
cat("\n--- Simulating data from the RBC model ---\n")
set.seed(123)
n <- 200
z_sim <- numeric(n)
k_sim <- numeric(n)
c_sim <- numeric(n)
k_sim[1] <- ss$values[["K"]]
c_sim[1] <- ss$values[["C"]]

for (t in 2:n) {
  z_sim[t] <- 0.9 * z_sim[t - 1] + 0.01 * rnorm(1)
  # Use the linear approximation for simulation
  x_dev <- c(z_sim[t], k_sim[t - 1] - ss$values[["K"]])
  c_dev <- as.numeric(sol$G %*% x_dev)
  c_sim[t] <- ss$values[["C"]] + c_dev
  k_sim[t] <- exp(z_sim[t]) * k_sim[t - 1]^0.33 - c_sim[t] +
              (1 - 0.025) * k_sim[t - 1]
}
dat <- data.frame(C = c_sim)

cat("Estimating rho from simulated data...\n")
fit <- estimate(rbc, data = dat,
                start = list(rho = 0.5),
                control = list(maxit = 300))
summary(fit)

cat("\nTrue rho: 0.9\n")
cat("Estimated rho:", round(coef(fit)["rho"], 4), "\n")

# --- 6. Forecast from estimated model ----------------------------------------
fc <- forecast(fit, horizon = 12)
plot(fc)

cat("\nDemo complete. Plots displayed:\n")
cat("  1. IRFs of consumption to TFP shock (40 periods)\n")
cat("  2. 12-period forecast with confidence bands\n")
