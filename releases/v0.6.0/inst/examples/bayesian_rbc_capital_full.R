# =============================================================================
# Bayesian Nonlinear RBC with Capital Accumulation
# =============================================================================
#
# A standard Real Business Cycle model with endogenous capital, estimated
# via Bayesian methods. This is a more realistic test of nonlinear Bayesian
# estimation because it features:
#   - an endogenous state variable (capital stock K)
#   - nonlinear production and Euler equation
#   - interaction between consumption, investment, and capital
#
# MODEL:
#   Euler equation:
#     1/C = beta / C(+1) * (alpha * exp(Z) * K(+1)^(alpha-1) + 1 - delta)
#   Capital accumulation:
#     K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K
#   Technology process:
#     Z(+1) = rho * Z
#
#   Variables:
#     C  — consumption (observed)
#     K  — capital stock (endogenous state, unobserved)
#     Z  — log-technology (exogenous state)
#
#   Parameters:
#     alpha — capital share (fixed = 0.33)
#     beta  — discount factor (fixed = 0.99)
#     delta — depreciation rate (fixed = 0.025)
#     rho   — technology persistence (estimated)
#
# STEADY STATE:
#   Z_ss = 0  (A_ss = exp(0) = 1)
#   K_ss = (alpha * beta / (1 - beta*(1-delta)))^(1/(1-alpha))
#   C_ss = K_ss^alpha - delta * K_ss
#
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Bayesian Nonlinear RBC with Capital Accumulation\n")
cat("============================================================\n\n")

# =============================================================================
# Step 1: Define the nonlinear RBC model
# =============================================================================
cat("--- Step 1: Define RBC model ---\n\n")

rbc <- dsgenl_model(
  # Euler equation
  "1/C = beta / C(+1) * (alpha * exp(Z) * K(+1)^(alpha-1) + 1 - delta)",
  # Capital accumulation
  "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
  # Technology process (linear in logs)
  "Z(+1) = rho * Z",
  observed = "C",
  endo_state = "K",
  exo_state = "Z",
  fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
  start = list(rho = 0.9),
  ss_guess = c(C = 2, K = 30, Z = 0)
)
print(rbc)

# =============================================================================
# Step 2: Verify steady state and solve at true parameters
# =============================================================================
cat("\n--- Step 2: Steady state and solution ---\n\n")

true_rho <- 0.9
true_sd <- 0.007  # Standard TFP shock size
true_params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, rho = true_rho)

ss <- steady_state(rbc, params = true_params)
cat("Steady state:\n")
cat(sprintf("  C  = %.4f\n", ss$values["C"]))
cat(sprintf("  K  = %.4f\n", ss$values["K"]))
cat(sprintf("  Z  = %.4f\n", ss$values["Z"]))

# Analytical check
K_ss_analytical <- (0.33 * 0.99 / (1 - 0.99 * (1 - 0.025)))^(1 / (1 - 0.33))
C_ss_analytical <- K_ss_analytical^0.33 - 0.025 * K_ss_analytical
cat(sprintf("\n  Analytical K_ss = %.4f (numerical: %.4f)\n",
            K_ss_analytical, ss$values["K"]))
cat(sprintf("  Analytical C_ss = %.4f (numerical: %.4f)\n",
            C_ss_analytical, ss$values["C"]))

sol <- solve_dsge(rbc, params = true_params, shock_sd = c(Z = true_sd))
cat("\nPolicy matrix G (controls = C):\n")
print(round(sol$G, 6))
cat("\nTransition matrix H (states = Z, K):\n")
print(round(sol$H, 6))
cat("\nStability:", ifelse(sol$stable, "STABLE", "UNSTABLE"), "\n")

# =============================================================================
# Step 3: Simulate data
# =============================================================================
cat("\n--- Step 3: Simulate data ---\n\n")

set.seed(42)
n <- 300
n_states <- ncol(sol$H)
states <- matrix(0, n, n_states)
colnames(states) <- colnames(sol$H)

for (t in 2:n) {
  eps <- rnorm(1)  # single shock (Z only)
  states[t, ] <- as.numeric(sol$H %*% states[t - 1, ] + sol$M %*% eps)
}

# Observed consumption in levels
C_dev <- as.numeric(states %*% t(sol$G))
C_level <- ss$values["C"] + C_dev
dat <- data.frame(C = C_level)

cat(sprintf("  Simulated %d observations\n", n))
cat(sprintf("  C: mean = %.4f (ss = %.4f), sd = %.4f\n",
            mean(dat$C), ss$values["C"], sd(dat$C)))

# =============================================================================
# Step 4: Specify priors
# =============================================================================
cat("\n--- Step 4: Prior specification ---\n\n")

priors <- list(
  rho = prior("beta", shape1 = 10, shape2 = 2)  # mean = 0.833, favors high persistence
)

cat("  Parameter  Prior           Mean    SD\n")
cat("  ---------  --------------  ------  -----\n")
cat("  rho        beta(10, 2)     0.833   0.103\n")
cat("  sd(e.Z)    inv_gamma(0.01, 0.01)   (diffuse)\n")

# =============================================================================
# Step 5: Bayesian estimation
# =============================================================================
cat("\n--- Step 5: MCMC estimation ---\n\n")
cat("  Chains: 2, Iterations: 5000, Warmup: 2500\n")
cat("  Model has endogenous capital (K) — each draw requires\n")
cat("  steady state + linearization + Kalman filter.\n\n")

t0 <- proc.time()
fit <- bayes_dsge(rbc, data = dat,
                  priors = priors,
                  chains = 2,
                  iter = 5000,
                  warmup = 2500,
                  seed = 42)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Elapsed: %.0f seconds (%.1f minutes)\n\n", elapsed, elapsed / 60))

# =============================================================================
# Step 6: Posterior summary
# =============================================================================
cat("--- Step 6: Posterior summary ---\n\n")
summary(fit)

# =============================================================================
# Step 7: Parameter recovery
# =============================================================================
cat("\n--- Step 7: Parameter recovery ---\n\n")

est <- coef(fit)
rho_draws <- as.numeric(fit$posterior[, "rho", ])
sd_draws <- as.numeric(fit$posterior[, "sd_e.Z", ])
rho_ci <- quantile(rho_draws, c(0.025, 0.975))
sd_ci <- quantile(sd_draws, c(0.025, 0.975))

cat(sprintf("  rho:\n"))
cat(sprintf("    True:           %.3f\n", true_rho))
cat(sprintf("    Posterior mean:  %.3f\n", est["rho"]))
cat(sprintf("    Posterior SD:    %.3f\n", sd(rho_draws)))
cat(sprintf("    95%% CI:         [%.3f, %.3f]\n", rho_ci[1], rho_ci[2]))
cat(sprintf("    True in CI:     %s\n",
            ifelse(true_rho >= rho_ci[1] & true_rho <= rho_ci[2], "YES", "no")))

cat(sprintf("\n  sd(e.Z):\n"))
cat(sprintf("    True:           %.4f\n", true_sd))
cat(sprintf("    Posterior mean:  %.4f\n", est["sd_e.Z"]))
cat(sprintf("    95%% CI:         [%.4f, %.4f]\n", sd_ci[1], sd_ci[2]))
cat(sprintf("    True in CI:     %s\n",
            ifelse(true_sd >= sd_ci[1] & true_sd <= sd_ci[2], "YES", "no")))

# Diagnostics
cat(sprintf("\n  Acceptance rates:  %s\n",
            paste(round(fit$acceptance_rates, 3), collapse = ", ")))
cat(sprintf("  ESS:               rho = %.0f, sd = %.0f\n",
            fit$diagnostics$ess[fit$diagnostics$parameter == "rho"],
            fit$diagnostics$ess[fit$diagnostics$parameter == "sd_e.Z"]))
cat(sprintf("  R-hat:             rho = %.4f, sd = %.4f\n",
            fit$diagnostics$rhat[fit$diagnostics$parameter == "rho"],
            fit$diagnostics$rhat[fit$diagnostics$parameter == "sd_e.Z"]))

if (!is.null(fit$solve_failures) && fit$solve_failures > 0) {
  cat(sprintf("  Solve failures:    %d\n", fit$solve_failures))
}

# =============================================================================
# Step 8: Impulse-response functions
# =============================================================================
cat("\n--- Step 8: Posterior IRFs ---\n\n")

irfs <- irf(fit, periods = 40, n_draws = 200)
print(irfs)

# Check IRF signs at impact
cat("\nTechnology shock (Z) at impact:\n")
for (v in c("C", "K", "Z")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "Z" &
                   irfs$data$response == v, ]
  if (nrow(val) > 0) {
    cat(sprintf("  %s: median = %+.6f  [%.6f, %.6f]\n",
                v, val$value, val$lower, val$upper))
  }
}

cat("\nTechnology shock (Z) at peak (period 5):\n")
for (v in c("C", "K")) {
  val <- irfs$data[irfs$data$period == 5 & irfs$data$impulse == "Z" &
                   irfs$data$response == v, ]
  if (nrow(val) > 0) {
    cat(sprintf("  %s: median = %+.6f  [%.6f, %.6f]\n",
                v, val$value, val$lower, val$upper))
  }
}

cat("\nEconomic interpretation:\n")
cat("  A positive technology shock should:\n")
cat("  - Raise consumption (C > 0) at impact\n")
cat("  - Build capital (K > 0) gradually via higher investment\n")
cat("  - Show hump-shaped capital response (peaks after several periods)\n")
cat("  - Show persistent consumption response (due to high rho)\n")

cat("\n============================================================\n")
cat(" Bayesian RBC with capital accumulation complete.\n")
cat("============================================================\n")
