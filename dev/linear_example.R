# =============================================================================
# Linear DSGE Example: Three-Equation New Keynesian Model
# =============================================================================
#
# Model:
#   Phillips curve:   p  = beta * E[p'] + kappa * x
#   IS curve:         x  = E[x'] - (r - E[p'] - g)
#   Taylor rule:      r  = psi * p + u
#   Monetary shock:   u' = rhou * u + e_u
#   Demand shock:     g' = rhog * g + e_g
#
# Observed:   p (inflation), r (interest rate)
# Unobserved: x (output gap)
# States:     u (monetary policy shock), g (demand shock)
# Fixed:      beta = 0.96
#
# This is the same model as Stata 19 DSGE Reference Manual, Intro 3a.
#
# To run in RStudio:
#   library(dsge)
#   source(system.file("examples", "linear_example.R", package = "dsge"))
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Linear DSGE Example: New Keynesian Model\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Step 1: Define the model
# -----------------------------------------------------------------------------
cat("--- Step 1: Model specification ---\n\n")

nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),
  unobs(x ~ lead(x) - (r - lead(p) - g)),
  obs(r   ~ psi * p + u),
  state(u ~ rhou * u),
  state(g ~ rhog * g),
  fixed = list(beta = 0.96),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)
print(nk)

# -----------------------------------------------------------------------------
# Step 2: Solve at known parameter values
# -----------------------------------------------------------------------------
cat("\n--- Step 2: Solve at Stata-documented parameter values ---\n\n")

true_params <- c(beta = 0.96, kappa = 0.085, psi = 1.94,
                 rhou = 0.70, rhog = 0.95)
true_sd <- c(u = 2.3, g = 0.57)

sol <- solve_dsge(nk, params = true_params, shock_sd = true_sd)
print(sol)

# -----------------------------------------------------------------------------
# Step 3: Policy matrix, transition matrix, stability
# -----------------------------------------------------------------------------
cat("\n--- Step 3: Structural matrices and stability ---\n\n")

cat("Policy matrix G (maps states to controls):\n")
G <- policy_matrix(sol, se = FALSE)
print(round(G, 6))

cat("\nTransition matrix H (state dynamics):\n")
H <- transition_matrix(sol, se = FALSE)
print(round(H, 6))

cat("\nStability check:\n")
stab <- stability(sol)
print(stab)

# -----------------------------------------------------------------------------
# Step 4: Simulate data
# -----------------------------------------------------------------------------
cat("\n--- Step 4: Simulate data from the model ---\n\n")

set.seed(42)
n <- 200
states <- matrix(0, n, 2)
colnames(states) <- c("u", "g")
for (t in 2:n) {
  states[t, "u"] <- 0.70 * states[t - 1, "u"] + 2.3 * rnorm(1)
  states[t, "g"] <- 0.95 * states[t - 1, "g"] + 0.57 * rnorm(1)
}

# Map states to observables
Z <- sol$D %*% sol$G
obs_data <- states %*% t(Z)
colnames(obs_data) <- c("p", "r")
dat <- as.data.frame(obs_data)
cat("Simulated", n, "observations of inflation (p) and interest rate (r).\n")
cat("  p: mean =", round(mean(dat$p), 3),
    ", sd =", round(sd(dat$p), 3), "\n")
cat("  r: mean =", round(mean(dat$r), 3),
    ", sd =", round(sd(dat$r), 3), "\n")

# -----------------------------------------------------------------------------
# Step 5: Estimate the model
# -----------------------------------------------------------------------------
cat("\n--- Step 5: Maximum likelihood estimation ---\n\n")

fit <- estimate(nk, data = dat,
                start = list(kappa = 0.1, psi = 1.5,
                             rhou = 0.7, rhog = 0.9),
                control = list(maxit = 500))
summary(fit)

cat("True vs estimated:\n")
cat(sprintf("  %-8s  true = %6.3f  estimated = %6.3f\n",
            "kappa", 0.085, coef(fit)["kappa"]))
cat(sprintf("  %-8s  true = %6.3f  estimated = %6.3f\n",
            "psi", 1.94, coef(fit)["psi"]))
cat(sprintf("  %-8s  true = %6.3f  estimated = %6.3f\n",
            "rhou", 0.70, coef(fit)["rhou"]))
cat(sprintf("  %-8s  true = %6.3f  estimated = %6.3f\n",
            "rhog", 0.95, coef(fit)["rhog"]))

# -----------------------------------------------------------------------------
# Step 6: Postestimation — policy and transition matrices with SEs
# -----------------------------------------------------------------------------
cat("\n--- Step 6: Postestimation matrices (with standard errors) ---\n\n")

cat("Policy matrix with delta-method SEs:\n")
print(policy_matrix(fit))

cat("\nTransition matrix with delta-method SEs:\n")
print(transition_matrix(fit))

cat("\nStability from estimated model:\n")
print(stability(fit))

# -----------------------------------------------------------------------------
# Step 7: Impulse-response functions
# -----------------------------------------------------------------------------
cat("\n--- Step 7: Impulse-response functions ---\n\n")

irfs <- irf(fit, periods = 20)
print(irfs)
plot(irfs)

# -----------------------------------------------------------------------------
# Step 8: Forecasting
# -----------------------------------------------------------------------------
cat("\n--- Step 8: Multi-step forecast ---\n\n")

fc <- forecast(fit, horizon = 12)
print(fc)
plot(fc)

cat("\n============================================================\n")
cat(" Linear DSGE example complete.\n")
cat(" Two sets of plots should be displayed:\n")
cat("   1. IRFs: p, r, x responses to u and g shocks\n")
cat("   2. Forecast: 12-period ahead for p and r\n")
cat("============================================================\n")
