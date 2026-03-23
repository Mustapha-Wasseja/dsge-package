# =============================================================================
# Validation: AR(1) Model (Basic Backward-Looking Example)
# =============================================================================
#
# This script validates the dsge package on the simplest possible DSGE
# model: a single observed variable driven by a single autoregressive
# state. This is a synthetic test with known analytical solution.
#
# Source: Package vignette / textbook baseline
# =============================================================================

library(dsge)

cat("=== Validation: AR(1) Model ===\n\n")

# --- Model ---
# y_t = z_t
# z_{t+1} = rho * z_t + e_{t+1}
# Analytical solution: G = 1, H = rho, M = sigma

m <- dsge_model(
  obs(y ~ z),
  state(z ~ rho * z),
  start = list(rho = 0.5)
)
print(m)

# --- Solve at known parameters ---
true_rho <- 0.8
true_sd  <- 1.0

sol <- solve_dsge(m, params = c(rho = true_rho),
                  shock_sd = c(z = true_sd))

cat("\n--- Analytical vs computed solution ---\n")
cat("G (should be 1):", as.numeric(sol$G), "\n")
cat("H (should be 0.8):", as.numeric(sol$H), "\n")
cat("M (should be 1):", as.numeric(sol$M), "\n")

stopifnot(abs(as.numeric(sol$G) - 1.0) < 1e-10)
stopifnot(abs(as.numeric(sol$H) - 0.8) < 1e-10)
stopifnot(abs(as.numeric(sol$M) - 1.0) < 1e-10)
cat("Solver: PASS (exact match)\n")

# --- Stability ---
stab <- stability(sol)
cat("\nStability:", stab$classification, "\n")
stopifnot(stab$stable)
cat("Stability: PASS\n")

# --- IRFs ---
irfs <- irf(sol, periods = 10, se = FALSE)
y_irf <- irfs$data[irfs$data$response == "y", "value"]
expected_irf <- 0.8^(0:10)
cat("\nIRF comparison (y to z shock):\n")
cat("  Computed:", round(y_irf, 6), "\n")
cat("  Expected:", round(expected_irf, 6), "\n")
stopifnot(max(abs(y_irf - expected_irf)) < 1e-10)
cat("IRFs: PASS (exact match)\n")

# --- Estimation ---
cat("\n--- Estimation on simulated data ---\n")
set.seed(42)
n <- 300
z <- numeric(n)
for (i in 2:n) z[i] <- true_rho * z[i - 1] + true_sd * rnorm(1)
dat <- data.frame(y = z)

fit <- estimate(m, data = dat, control = list(maxit = 200))
cat("Convergence:", fit$convergence, "\n")
cat("Estimated rho:", round(coef(fit)["rho"], 4),
    "(true: 0.8)\n")
cat("Estimated sd:", round(coef(fit)["sd(e.z)"], 4),
    "(true: 1.0)\n")
cat("Log-likelihood:", round(fit$loglik, 2), "\n")

stopifnot(fit$convergence == 0)
stopifnot(abs(coef(fit)["rho"] - 0.8) < 0.15)
cat("Estimation: PASS\n")

# --- Postestimation ---
cat("\n--- Postestimation ---\n")
cat("Policy matrix:\n")
print(policy_matrix(fit, se = FALSE))
cat("Transition matrix:\n")
print(transition_matrix(fit, se = FALSE))
cat("Stability:\n")
print(stability(fit))

# --- Forecast ---
fc <- forecast(fit, horizon = 8)
cat("\nForecast (8 periods):\n")
print(fc)

# --- Predict ---
pred <- predict(fit, type = "observed")
resid <- residuals(fit)
cat("\nPredictions: ", nrow(pred), "rows\n")
cat("Residuals: ", nrow(resid), "rows\n")
stopifnot(nrow(pred) == n)
stopifnot(nrow(resid) == n)
cat("Predict/residuals: PASS\n")

# --- Plots ---
cat("\n--- Generating plots ---\n")
irfs_plot <- irf(fit, periods = 20)
plot(irfs_plot)

fc_plot <- forecast(fit, horizon = 12)
plot(fc_plot)

cat("\n=== AR(1) validation complete: ALL PASS ===\n")
