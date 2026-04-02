# =============================================================================
# Bayesian Nonlinear RBC — Quick Demo
# =============================================================================
#
# A lighter version of the RBC-with-capital example for interactive use.
# Uses shorter chains to run in ~2-4 minutes rather than ~9 minutes.
#
# MODEL:
#   Euler:    1/C = beta/C(+1) * (alpha*exp(Z)*K(+1)^(alpha-1) + 1 - delta)
#   Capital:  K(+1) = exp(Z)*K^alpha - C + (1 - delta)*K
#   Tech:     Z(+1) = rho * Z
#
#   Observed: C (consumption)
#   States:   K (endogenous capital), Z (exogenous technology)
#   Fixed:    alpha = 0.33, beta = 0.99, delta = 0.025
#   Estimated: rho (technology persistence), sd(e.Z) (shock size)
#
# For the full validation version with longer chains and detailed
# diagnostics, see bayesian_rbc_capital_full.R.
#
# TO RUN:
#   library(dsge)
#   source(system.file("examples", "bayesian_rbc_capital_demo.R", package = "dsge"))
#
# EXPECTED RUNTIME: ~2-4 minutes
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Bayesian Nonlinear RBC — Quick Demo\n")
cat("============================================================\n\n")

# --- Define model ---
rbc <- dsgenl_model(
  "1/C = beta / C(+1) * (alpha * exp(Z) * K(+1)^(alpha-1) + 1 - delta)",
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

# --- Simulate data ---
set.seed(42)
true_params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, rho = 0.9)
sol <- solve_dsge(rbc, params = true_params, shock_sd = c(Z = 0.007))
ss <- sol$steady_state

n <- 200
states <- matrix(0, n, ncol(sol$H))
colnames(states) <- colnames(sol$H)
for (t in 2:n) {
  states[t, ] <- as.numeric(sol$H %*% states[t - 1, ] + sol$M %*% rnorm(1))
}
C_dev <- as.numeric(states %*% t(sol$G))
dat <- data.frame(C = ss["C"] + C_dev)

cat(sprintf("\nSteady state: C = %.3f, K = %.3f\n", ss["C"], ss["K"]))
cat(sprintf("Simulated %d observations (sd = %.4f)\n\n", n, sd(dat$C)))

# --- Estimate (short chains for demo) ---
cat("Running MCMC (2 chains x 2000 iter, 1000 warmup)...\n")
t0 <- proc.time()
fit <- bayes_dsge(rbc, data = dat,
                  priors = list(rho = prior("beta", shape1 = 10, shape2 = 2)),
                  chains = 2, iter = 2000, warmup = 1000,
                  seed = 42)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("Elapsed: %.0f seconds\n\n", elapsed))

# --- Results ---
summary(fit)

est <- coef(fit)
cat(sprintf("\nTrue rho = 0.900, estimated = %.3f\n", est["rho"]))
cat(sprintf("True sd  = 0.007, estimated = %.4f\n", est["sd_e.Z"]))

# --- IRFs ---
cat("\nComputing posterior IRFs...\n")
irfs <- irf(fit, periods = 20, n_draws = 100)

# Plot
plot(irfs)

cat("\nTechnology shock at impact:\n")
for (v in c("C", "K", "Z")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "Z" &
                   irfs$data$response == v, ]
  if (nrow(val) > 0) {
    cat(sprintf("  %s = %+.5f  [%+.5f, %+.5f]\n",
                v, val$value, val$lower, val$upper))
  }
}

cat("\nDone. For full validation, run bayesian_rbc_capital_full.R\n")
