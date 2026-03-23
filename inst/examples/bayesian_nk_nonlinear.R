# =============================================================================
# Bayesian Estimation of a Nonlinear New Keynesian Model
# =============================================================================
#
# This example estimates a 3-equation New Keynesian model written in
# nonlinear form. The model is linearized internally at each MCMC draw,
# producing parameter-specific steady states and policy functions.
#
# MODEL (Nonlinear New Keynesian):
#   Euler equation:    1 = beta * (x / x(+1)) * (1/z) * (r / p(+1))
#   Phillips curve:    (theta-1) + phi*(p-1)*p = theta*x + phi*beta*(p(+1)-1)*p(+1)
#   Taylor rule:       beta * r = p^psi * m
#   Monetary shock:    m(+1) = m^rhom
#   Demand shock:      z(+1) = z^rhoz
#
# VARIABLES:
#   Observed:    p (inflation gross rate), r (interest rate gross)
#   Unobserved:  x (output gap)
#   Exo states:  z (demand/preference shock), m (monetary policy shock)
#
# STEADY STATE:
#   p_ss = 1 (zero inflation), r_ss = 1/beta, x_ss = (theta-1)/theta,
#   z_ss = 1, m_ss = 1
#
# ESTIMATION:
#   Parameters estimated: phi, psi, rhom, rhoz (+ 2 shock SDs)
#   Parameters fixed: beta = 0.96, theta = 5
#
# This demonstrates the package's ability to handle genuine nonlinear
# DSGE models with Bayesian methods. The steady state and linearization
# are recomputed for each candidate parameter vector.
#
# REQUIREMENTS:
#   - dsge package (>= 0.4.0)
#
# TO RUN:
#   library(dsge)
#   source(system.file("examples", "bayesian_nk_nonlinear.R", package = "dsge"))
#
# EXPECTED RUNTIME: ~15-25 minutes depending on hardware
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Bayesian Nonlinear NK Model — Simulated Data\n")
cat("============================================================\n\n")

# =============================================================================
# Step 1: Define the nonlinear NK model
# =============================================================================
cat("--- Step 1: Define nonlinear NK model ---\n\n")

nk_nl <- dsgenl_model(
  # Euler equation
  "1 = beta * (x / x(+1)) * (1/z) * (r / p(+1))",
  # Phillips curve (Rotemberg pricing)
  "(theta - 1) + phi * (p - 1) * p = theta * x + phi * beta * (p(+1) - 1) * p(+1)",
  # Taylor rule
  "beta * r = p^psi * m",
  # State equations
  "m(+1) = m^rhom",
  "z(+1) = z^rhoz",
  observed = c("p", "r"),
  unobserved = "x",
  exo_state = c("z", "m"),
  fixed = list(beta = 0.96, theta = 5),
  start = list(phi = 47, psi = 1.9, rhom = 0.7, rhoz = 0.95),
  ss_guess = c(x = 0.8, p = 1, r = 1.04, z = 1, m = 1)
)
print(nk_nl)

# =============================================================================
# Step 2: Simulate data from known parameters
# =============================================================================
cat("\n--- Step 2: Simulate data ---\n\n")

true_params <- c(beta = 0.96, theta = 5, phi = 47, psi = 1.9,
                 rhom = 0.7, rhoz = 0.95)
true_shock_sd <- c(z = 0.01, m = 0.01)

sol <- solve_dsge(nk_nl, params = true_params, shock_sd = true_shock_sd)
cat("  Steady state:\n")
cat(sprintf("    x = %.4f, p = %.4f, r = %.4f, z = %.4f, m = %.4f\n",
            sol$steady_state["x"], sol$steady_state["p"],
            sol$steady_state["r"], sol$steady_state["z"],
            sol$steady_state["m"]))

set.seed(42)
n <- 200
n_states <- ncol(sol$H)
states <- matrix(0, n, n_states)
colnames(states) <- colnames(sol$H)
for (t in 2:n) {
  eps <- rnorm(length(true_shock_sd))
  states[t, ] <- as.numeric(sol$H %*% states[t - 1, ] + sol$M %*% eps)
}
# Observed variables in levels
obs_dev <- states %*% t(sol$G)
ss_obs <- sol$steady_state[c("p", "r")]
obs_levels <- sweep(obs_dev[, 1:2, drop = FALSE], 2, ss_obs, "+")
colnames(obs_levels) <- c("p", "r")
dat <- as.data.frame(obs_levels)

cat(sprintf("  Simulated %d observations\n", n))
cat(sprintf("  p: mean = %.4f, sd = %.4f\n", mean(dat$p), sd(dat$p)))
cat(sprintf("  r: mean = %.4f, sd = %.4f\n", mean(dat$r), sd(dat$r)))

# =============================================================================
# Step 3: Specify priors
# =============================================================================
cat("\n--- Step 3: Prior specification ---\n\n")

priors_nl <- list(
  phi   = prior("gamma", shape = 16, rate = 0.34),   # mean ~47, moderately informative
  psi   = prior("gamma", shape = 10, rate = 5.26),   # mean ~1.9
  rhom  = prior("beta", shape1 = 7, shape2 = 3),     # mean 0.7
  rhoz  = prior("beta", shape1 = 19, shape2 = 1)     # mean 0.95
)

cat("  Parameter  Prior                   Mean   SD\n")
cat("  ---------  ----------------------  -----  -----\n")
cat("  phi        gamma(16, 0.34)         47.1   11.8\n")
cat("  psi        gamma(10, 5.26)         1.90   0.60\n")
cat("  rhom       beta(7, 3)              0.70   0.14\n")
cat("  rhoz       beta(19, 1)             0.95   0.05\n")
cat("  sd(e.z)    inv_gamma(0.01, 0.01)   (diffuse)\n")
cat("  sd(e.m)    inv_gamma(0.01, 0.01)   (diffuse)\n")

# =============================================================================
# Step 4: Bayesian estimation
# =============================================================================
cat("\n--- Step 4: MCMC estimation ---\n\n")
cat("  Chains: 2, Iterations: 5000, Warmup: 2500\n")
cat("  This is a nonlinear model — each draw requires steady-state\n")
cat("  solving and re-linearization, so expect longer runtime.\n\n")

t0 <- proc.time()
fit_nl <- bayes_dsge(nk_nl, data = dat,
                     priors = priors_nl,
                     chains = 2,
                     iter = 5000,
                     warmup = 2500,
                     seed = 123)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Elapsed: %.0f seconds (%.1f minutes)\n\n", elapsed, elapsed / 60))

# =============================================================================
# Step 5: Posterior summary
# =============================================================================
cat("--- Step 5: Posterior summary ---\n\n")
summary(fit_nl)

# =============================================================================
# Step 6: Parameter recovery
# =============================================================================
cat("\n--- Step 6: Parameter recovery ---\n\n")

est <- coef(fit_nl)
true_vals <- c(phi = 47, psi = 1.9, rhom = 0.7, rhoz = 0.95)

cat("  Parameter  True    Posterior   In 95% CI?\n")
cat("  ---------  ------  ---------  ----------\n")
for (nm in names(true_vals)) {
  post_draws <- as.numeric(fit_nl$posterior[, nm, ])
  ci <- quantile(post_draws, c(0.025, 0.975))
  in_ci <- true_vals[nm] >= ci[1] && true_vals[nm] <= ci[2]
  cat(sprintf("  %-8s   %6.3f  %9.3f  %s\n",
              nm, true_vals[nm], est[nm], ifelse(in_ci, "YES", "no")))
}

if (!is.null(fit_nl$solve_failures) && fit_nl$solve_failures > 0) {
  cat(sprintf("\n  Solve failures during MCMC: %d\n", fit_nl$solve_failures))
}

# =============================================================================
# Step 7: Impulse-response functions
# =============================================================================
cat("\n--- Step 7: Posterior IRFs ---\n\n")

irfs <- irf(fit_nl, periods = 8, n_draws = 200)
print(irfs)

cat("\nMonetary shock (m) at impact:\n")
for (v in c("p", "x", "r")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "m" &
                   irfs$data$response == v, "value"]
  if (length(val) > 0) cat(sprintf("  %s = %+.4f\n", v, val))
}

cat("\nDemand shock (z) at impact:\n")
for (v in c("p", "x", "r")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "z" &
                   irfs$data$response == v, "value"]
  if (length(val) > 0) cat(sprintf("  %s = %+.4f\n", v, val))
}

cat("\n============================================================\n")
cat(" Nonlinear NK Bayesian estimation complete.\n")
cat("============================================================\n")
