# =============================================================================
# Bayesian Estimation of a Nonlinear RBC Model
# =============================================================================
#
# This example demonstrates Bayesian estimation of a nonlinear DSGE model
# using the dsge package. The model is a simple Real Business Cycle (RBC)
# model specified in nonlinear form, estimated via first-order perturbation
# around the parameter-specific steady state.
#
# MODEL:
#   The model has one observed variable (output C), one unobserved control,
#   two state variables, and a standard RBC structure:
#
#   C = K^alpha * A - K(+1) + (1-delta)*K    [resource constraint]
#   1 = beta * (alpha * K(+1)^(alpha-1) * A(+1) + 1 - delta) * (C / C(+1))
#                                              [Euler equation, simplified]
#   K(+1) = ...                               [implied by resource constraint]
#   A(+1) = A^rho                             [technology process]
#
#   For simplicity, we use a compact 2-variable version:
#     C = A * K^alpha                          [production, all consumed]
#     A(+1) = A^rho                            [technology AR(1) in levels]
#
#   with K as an endogenous state following capital accumulation.
#
#   However, for the clearest demonstration, we use an even simpler
#   nonlinear model:
#     Y = Z                                    [output = technology]
#     Z(+1) = Z^rho                            [nonlinear AR(1)]
#
#   This is the simplest nonlinear model: a single observed variable
#   driven by a nonlinear autoregressive technology process.
#   At steady state, Z_ss = 1 (for rho < 1).
#
# ESTIMATION:
#   The Bayesian estimator re-solves the steady state and re-linearizes
#   the model at each candidate parameter draw. This makes it more
#   expensive than linear Bayesian estimation, but correctly accounts
#   for the nonlinear structure.
#
# REQUIREMENTS:
#   - dsge package (>= 0.4.0)
#
# TO RUN:
#   library(dsge)
#   source(system.file("examples", "bayesian_rbc_nonlinear.R", package = "dsge"))
#
# EXPECTED RUNTIME: ~2-5 minutes depending on hardware
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Bayesian Nonlinear DSGE — Parameter Recovery Example\n")
cat("============================================================\n\n")

# =============================================================================
# Step 1: Define the nonlinear model
# =============================================================================
cat("--- Step 1: Define nonlinear model ---\n\n")

nl_model <- dsgenl_model(
  "y = z",
  "z(+1) = z^rho",
  observed = "y",
  unobserved = character(0),
  exo_state = "z",
  start = list(rho = 0.5),
  ss_guess = c(y = 1, z = 1)
)
print(nl_model)

# =============================================================================
# Step 2: Simulate data from known parameters
# =============================================================================
cat("\n--- Step 2: Simulate data ---\n\n")

true_rho <- 0.85
true_sd <- 0.5

set.seed(123)
sol <- solve_dsge(nl_model, params = c(rho = true_rho),
                  shock_sd = c(z = true_sd))
cat(sprintf("  True parameters: rho = %.2f, sd = %.2f\n", true_rho, true_sd))
cat(sprintf("  Steady state: y = %.2f, z = %.2f\n",
            sol$steady_state["y"], sol$steady_state["z"]))
cat(sprintf("  H (transition): %.4f\n", sol$H[1,1]))
cat(sprintf("  G (policy):     %.4f\n", sol$G[1,1]))

n <- 300
z_dev <- numeric(n)
for (i in 2:n) {
  z_dev[i] <- sol$H[1,1] * z_dev[i - 1] + sol$M[1,1] * rnorm(1)
}
# Data in levels (around steady state = 1)
dat <- data.frame(y = sol$steady_state["y"] + z_dev)

cat(sprintf("  Simulated %d observations\n", n))
cat(sprintf("  y: mean = %.3f, sd = %.3f\n", mean(dat$y), sd(dat$y)))

# =============================================================================
# Step 3: Specify priors
# =============================================================================
cat("\n--- Step 3: Prior specification ---\n\n")

priors <- list(
  rho = prior("beta", shape1 = 5, shape2 = 2)  # Prior mean = 0.71, favors persistence
)

cat("  Parameter  Prior           Mean   SD\n")
cat("  ---------  --------------  -----  -----\n")
cat("  rho        beta(5, 2)      0.714  0.160\n")
cat("  sd(e.z)    inv_gamma(0.01, 0.01)  (diffuse)\n")

# =============================================================================
# Step 4: Bayesian estimation
# =============================================================================
cat("\n--- Step 4: MCMC estimation ---\n\n")
cat("  Chains: 2, Iterations: 5000, Warmup: 2500\n")
cat("  Note: Nonlinear models require re-solving the steady state\n")
cat("  and re-linearizing at each draw, so this is slower than\n")
cat("  linear Bayesian estimation.\n\n")

t0 <- proc.time()
fit <- bayes_dsge(nl_model, data = dat,
                  priors = priors,
                  chains = 2,
                  iter = 5000,
                  warmup = 2500,
                  seed = 42)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Elapsed: %.0f seconds (%.1f minutes)\n\n", elapsed, elapsed / 60))

# =============================================================================
# Step 5: Posterior summary
# =============================================================================
cat("--- Step 5: Posterior summary ---\n\n")
summary(fit)

# =============================================================================
# Step 6: Parameter recovery check
# =============================================================================
cat("\n--- Step 6: Parameter recovery ---\n\n")

est <- coef(fit)
diag_df <- fit$diagnostics

cat(sprintf("  True rho:      %.3f\n", true_rho))
cat(sprintf("  Posterior mean: %.3f\n", est["rho"]))
cat(sprintf("  Posterior SD:   %.3f\n", diag_df$mcse[diag_df$parameter == "rho"]))
cat(sprintf("  95%% CI:        [%.3f, %.3f]\n",
            quantile(fit$posterior[, "rho", ], 0.025),
            quantile(fit$posterior[, "rho", ], 0.975)))

rho_in_ci <- (true_rho >= quantile(fit$posterior[, "rho", ], 0.025) &&
              true_rho <= quantile(fit$posterior[, "rho", ], 0.975))
cat(sprintf("  True rho in 95%% CI: %s\n", ifelse(rho_in_ci, "YES", "NO")))

cat(sprintf("\n  True sd:       %.3f\n", true_sd))
cat(sprintf("  Posterior mean: %.3f\n", est["sd_e.z"]))

if (!is.null(fit$solve_failures) && fit$solve_failures > 0) {
  cat(sprintf("\n  Solve failures during MCMC: %d\n", fit$solve_failures))
}

# =============================================================================
# Step 7: Impulse-response functions
# =============================================================================
cat("\n--- Step 7: Posterior IRFs ---\n\n")

irfs <- irf(fit, periods = 12, n_draws = 200)
print(irfs)

cat("\nTechnology shock at impact:\n")
val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "z" &
                 irfs$data$response == "y", "value"]
cat(sprintf("  y = %+.4f  (positive: output rises with technology)\n", val))

cat("\n============================================================\n")
cat(" Nonlinear Bayesian estimation complete.\n")
cat("============================================================\n")
