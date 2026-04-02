# =============================================================================
# Nonlinear DSGE Example: Real Business Cycle (RBC) Model
# =============================================================================
#
# Model (in nonlinear form):
#   Euler equation:        1/C = beta/C(+1) * (alpha*exp(Z)*K^(alpha-1) + 1-delta)
#   Capital accumulation:  K(+1) = exp(Z)*K^alpha - C + (1-delta)*K
#   TFP process:           Z(+1) = rho * Z + e_Z
#
# Observed:        C (consumption)
# Endogenous state: K (capital stock)
# Exogenous state:  Z (total factor productivity, with shock)
#
# Fixed:  alpha = 0.33, beta = 0.99, delta = 0.025
# Free:   rho (TFP persistence)
#
# The model is solved by first-order perturbation: compute the deterministic
# steady state, linearize around it, and solve the resulting linear system.
#
# To run in RStudio:
#   library(dsge)
#   source(system.file("examples", "nonlinear_example.R", package = "dsge"))
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Nonlinear DSGE Example: RBC Model\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Step 1: Define the nonlinear model
# -----------------------------------------------------------------------------
cat("--- Step 1: Model specification ---\n\n")

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

# -----------------------------------------------------------------------------
# Step 2: Compute deterministic steady state
# -----------------------------------------------------------------------------
cat("\n--- Step 2: Steady state ---\n\n")

params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, rho = 0.9)
ss <- steady_state(rbc, params = params)
print(ss)

# Verify analytically
K_ss_exact <- (0.33 / (1 / 0.99 - 1 + 0.025))^(1 / 0.67)
C_ss_exact <- K_ss_exact^0.33 - 0.025 * K_ss_exact
cat("\nAnalytical check:\n")
cat("  K_ss (analytical):", round(K_ss_exact, 4),
    "  (numerical):", round(ss$values[["K"]], 4), "\n")
cat("  C_ss (analytical):", round(C_ss_exact, 4),
    "  (numerical):", round(ss$values[["C"]], 4), "\n")
cat("  Z_ss (analytical): 0",
    "       (numerical):", round(ss$values[["Z"]], 8), "\n")

# -----------------------------------------------------------------------------
# Step 3: Linearize around steady state
# -----------------------------------------------------------------------------
cat("\n--- Step 3: Linearize ---\n\n")

lin <- linearize(rbc, ss, params = params)
cat("Structural matrices from linearization:\n")
cat("  A0 (control eq, current controls):\n")
print(round(lin$A0, 6))
cat("  A1 (control eq, lead controls):\n")
print(round(lin$A1, 6))
cat("  A3 (control eq, current states):\n")
print(round(lin$A3, 6))
cat("  B0 (state eq, lead states):\n")
print(round(lin$B0, 6))
cat("  B3 (state eq, current states):\n")
print(round(lin$B3, 6))

# -----------------------------------------------------------------------------
# Step 4: Solve the linearized system
# -----------------------------------------------------------------------------
cat("\n--- Step 4: Solve ---\n\n")

sol <- solve_dsge(rbc, params = params, shock_sd = c(Z = 0.01))
print(sol)

cat("\nPolicy matrix G (maps states to controls):\n")
print(round(policy_matrix(sol, se = FALSE), 6))

cat("\nTransition matrix H (state dynamics):\n")
print(round(transition_matrix(sol, se = FALSE), 6))

cat("\nStability:\n")
print(stability(sol))

# -----------------------------------------------------------------------------
# Step 5: Impulse-response functions
# -----------------------------------------------------------------------------
cat("\n--- Step 5: IRFs ---\n\n")

irfs <- irf(sol, periods = 40, se = FALSE)
print(irfs)
plot(irfs)

# -----------------------------------------------------------------------------
# Step 6: Simulate data and estimate
# -----------------------------------------------------------------------------
cat("\n--- Step 6: Estimation ---\n\n")
cat("Simulating 200 observations from the RBC model...\n")

set.seed(123)
n <- 200

# Simulate from the LINEAR state-space form to ensure consistency
# with the linearized estimator:
#   x_{t+1} = H * x_t + M * e_{t+1}   (state deviations)
#   y_t     = G * x_t                  (control deviations)
x_dev <- matrix(0, nrow = n, ncol = ncol(sol$H))
colnames(x_dev) <- colnames(sol$H)

for (t in 2:n) {
  e <- rnorm(1)
  x_dev[t, ] <- as.numeric(sol$H %*% x_dev[t - 1, ] + sol$M %*% e)
}

# Map state deviations to observed control deviations
c_dev <- x_dev %*% t(sol$G)
c_sim <- ss$values[["C"]] + c_dev[, 1]
dat <- data.frame(C = c_sim)

cat("  C: mean =", round(mean(dat$C), 3),
    ", sd =", round(sd(dat$C), 3), "\n")
cat("  C_ss =", round(ss$values[["C"]], 3), "\n\n")

cat("Estimating rho by maximum likelihood...\n")
fit <- estimate(rbc, data = dat,
                start = list(rho = 0.5),
                control = list(maxit = 300))
summary(fit)

cat("True rho = 0.9,  Estimated rho =",
    round(coef(fit)["rho"], 4), "\n")

# -----------------------------------------------------------------------------
# Step 7: Forecast from estimated model
# -----------------------------------------------------------------------------
cat("\n--- Step 7: Forecast ---\n\n")

fc <- forecast(fit, horizon = 12)
print(fc)
plot(fc)

cat("\n============================================================\n")
cat(" Nonlinear DSGE example complete.\n")
cat(" Two sets of plots should be displayed:\n")
cat("   1. IRF: consumption response to TFP shock (40 periods)\n")
cat("   2. Forecast: 12-period ahead consumption forecast\n")
cat("============================================================\n")
