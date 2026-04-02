# =============================================================================
# Demo: New Keynesian Model — Full Workflow with Plots
# =============================================================================
#
# This script demonstrates the dsge package on a three-equation New
# Keynesian model with two observed variables and one unobserved variable.
#
# To run in RStudio:
#   library(dsge)
#   source(system.file("examples", "demo_nk.R", package = "dsge"))
# =============================================================================

library(dsge)

# --- 1. Specify the model ----------------------------------------------------
# Phillips curve:  p  = beta * E[p'] + kappa * x
# IS curve:        x  = E[x'] - (r - E[p'] - g)
# Taylor rule:     r  = psi * p + u
# Monetary shock:  u' = rhou * u + e_u
# Demand shock:    g' = rhog * g + e_g
#
# Observed: p (inflation), r (interest rate)
# Unobserved: x (output gap)

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

# --- 2. Simulate data from known parameters ----------------------------------
# Solve at true parameters to get the state-space matrices
true_params <- c(beta = 0.96, kappa = 0.085, psi = 1.94,
                 rhou = 0.70, rhog = 0.95)
true_sd <- c(u = 2.3, g = 0.57)

sol <- solve_dsge(nk, params = true_params, shock_sd = true_sd)

set.seed(123)
n <- 200
states <- matrix(0, n, 2)
colnames(states) <- c("u", "g")
for (t in 2:n) {
  states[t, "u"] <- 0.70 * states[t - 1, "u"] + 2.3 * rnorm(1)
  states[t, "g"] <- 0.95 * states[t - 1, "g"] + 0.57 * rnorm(1)
}

# Map states to observables through the state-space solution
Z <- sol$D %*% sol$G
obs_data <- states %*% t(Z)
colnames(obs_data) <- c("p", "r")
dat <- as.data.frame(obs_data)

cat("Simulated", n, "observations of inflation (p) and interest rate (r)\n\n")

# --- 3. Estimate by maximum likelihood ----------------------------------------
fit <- estimate(nk, data = dat,
                start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
                control = list(maxit = 500))
summary(fit)

# --- 4. Postestimation -------------------------------------------------------
cat("\nPolicy matrix:\n")
print(policy_matrix(fit))

cat("\nTransition matrix:\n")
print(transition_matrix(fit))

cat("\nStability diagnostics:\n")
print(stability(fit))

# --- 5. Impulse-response functions (with plots) ------------------------------
irfs <- irf(fit, periods = 20)
plot(irfs)

# --- 6. Forecasting (with plot) -----------------------------------------------
fc <- forecast(fit, horizon = 12)
plot(fc)

cat("\nDemo complete. Two sets of plots should be displayed:\n")
cat("  1. IRFs: responses of p, r, x to monetary (u) and demand (g) shocks\n")
cat("  2. Forecast: 12-period ahead forecasts for p and r\n")
