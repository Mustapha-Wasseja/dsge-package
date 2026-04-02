# =============================================================================
# Demo: AR(1) Model — Full Workflow with Plots
# =============================================================================
#
# This script demonstrates the complete dsge workflow on the simplest
# possible model: a single observed variable driven by an AR(1) state.
#
# To run in RStudio:
#   library(dsge)
#   source(system.file("examples", "demo_ar1.R", package = "dsge"))
# =============================================================================

library(dsge)

# --- 1. Simulate data --------------------------------------------------------
set.seed(42)
true_rho <- 0.8
true_sd  <- 1.0
n <- 200

z <- numeric(n)
for (i in 2:n) z[i] <- true_rho * z[i - 1] + true_sd * rnorm(1)
dat <- data.frame(y = z)

cat("Simulated", n, "observations with rho =", true_rho,
    "and sigma =", true_sd, "\n\n")

# --- 2. Specify the model ----------------------------------------------------
m <- dsge_model(
  obs(y ~ z),
  state(z ~ rho * z),
  start = list(rho = 0.5)
)
print(m)

# --- 3. Estimate by maximum likelihood ----------------------------------------
fit <- estimate(m, data = dat)
summary(fit)

# --- 4. Postestimation diagnostics --------------------------------------------
cat("\nPolicy matrix (maps states to controls):\n")
print(policy_matrix(fit))

cat("\nTransition matrix (state dynamics):\n")
print(transition_matrix(fit))

cat("\nStability check:\n")
print(stability(fit))

# --- 5. Impulse-response functions (with plot) --------------------------------
irfs <- irf(fit, periods = 20)
plot(irfs)

# --- 6. Forecasting (with plot) -----------------------------------------------
fc <- forecast(fit, horizon = 12)
plot(fc)

cat("\nDemo complete. Two plots should be displayed:\n")
cat("  1. IRF: response of y to a shock to z\n")
cat("  2. Forecast: 12-period ahead forecast with confidence bands\n")
