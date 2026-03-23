# =============================================================================
# Stata Comparison: Nonlinear New Keynesian Model (Intro 3d)
# =============================================================================
#
# This replicates Stata 19's dsgenl example from DSGE Reference Manual,
# Intro 3d: "Nonlinear New Keynesian model".
#
# Model (5 equations):
#   Euler:        1 = beta * (x/x(+1)) * (1/z) * (r/p(+1))
#   Phillips:     (theta-1) + phi*(p-1)*p = theta*x + phi*beta*(p(+1)-1)*p(+1)
#   Taylor rule:  beta * r = p^psi * m
#   Monetary:     log(m(+1)) = rhom * log(m)
#   Demand:       log(z(+1)) = rhoz * log(z)
#
# Variables:
#   Observed controls:   p (inflation), r (interest rate)
#   Unobserved controls: x (output gap)
#   Exogenous states:    z (demand shock), m (monetary shock)
#
# Stata's estimated parameters (from Intro 3d, p. 46):
#   beta = 0.96 (constrained), theta = 5 (constrained)
#   phi = 47.07939, psi = 1.943008
#   rhom = 0.7005489, rhoz = 0.9545255
#   sd(e.z) = 0.5689908, sd(e.m) = 2.318208
#
# Stata's results to compare against:
#   Steady state:  x = 0.8, p = 1, r = 1/0.96, z = 1, m = 1
#   Policy matrix:
#     G[x,z] =  0.952921    G[x,m] = -1.608216
#     G[p,z] =  0.9678135   G[p,m] = -0.4172515
#     G[r,z] =  1.880469    G[r,m] =  0.189277
#   Transition matrix:
#     H[z,z] = 0.9545255    H[m,m] = 0.7005489
#
# To run:
#   library(dsge)
#   source(system.file("examples", "stata_comparison.R", package = "dsge"))
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Stata Comparison: Nonlinear NK Model (Intro 3d)\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Step 1: Define the model (matching Stata's dsgenl specification)
# -----------------------------------------------------------------------------
cat("--- Step 1: Model specification ---\n\n")

nk_nl <- dsgenl_model(
  # Euler equation
  "1 = beta * (x / x(+1)) * (1/z) * (r / p(+1))",
  # Phillips curve
  "(theta - 1) + phi * (p - 1) * p = theta * x + phi * beta * (p(+1) - 1) * p(+1)",
  # Taylor rule
  "beta * r = p^psi * m",
  # Monetary shock: log(m(+1)) = rhom*log(m)
  # Rewritten as: m(+1) = m^rhom  (equivalent in levels)
  "m(+1) = m^rhom",
  # Demand shock: log(z(+1)) = rhoz*log(z)
  # Rewritten as: z(+1) = z^rhoz  (equivalent in levels)
  "z(+1) = z^rhoz",
  observed = c("p", "r"),
  unobserved = "x",
  exo_state = c("z", "m"),
  fixed = list(beta = 0.96, theta = 5),
  start = list(phi = 47, psi = 1.9, rhom = 0.7, rhoz = 0.95),
  ss_guess = c(x = 0.8, p = 1, r = 1.04, z = 1, m = 1)
)
print(nk_nl)

# -----------------------------------------------------------------------------
# Step 2: Compute steady state
# -----------------------------------------------------------------------------
cat("\n--- Step 2: Steady state ---\n\n")

# Use Stata's estimated parameter values
stata_params <- c(
  beta = 0.96, theta = 5,
  phi = 47.07939, psi = 1.943008,
  rhom = 0.7005489, rhoz = 0.9545255
)

ss <- steady_state(nk_nl, params = stata_params)
print(ss)

# Analytical steady state
cat("\nAnalytical check:\n")
cat("  x_ss = (theta-1)/theta =", (5 - 1) / 5, "\n")
cat("  p_ss = 1\n")
cat("  r_ss = 1/beta =", 1 / 0.96, "\n")
cat("  z_ss = 1\n")
cat("  m_ss = 1\n")

# -----------------------------------------------------------------------------
# Step 3: Solve at Stata's estimated values
# -----------------------------------------------------------------------------
cat("\n--- Step 3: Solve ---\n\n")

stata_sd <- c(z = 0.5689908, m = 2.318208)
sol <- solve_dsge(nk_nl, params = stata_params, shock_sd = stata_sd)
print(sol)

# -----------------------------------------------------------------------------
# Step 4: Compare policy matrix with Stata
# -----------------------------------------------------------------------------
cat("\n--- Step 4: Compare with Stata's results ---\n\n")

G <- policy_matrix(sol, se = FALSE)

# IMPORTANT: Stata's dsgenl reports the policy matrix in PERCENTAGE deviation
# units (log-linearization), while our package uses LEVEL deviations.
#
# Conversion: G_stata[v,s] = G_ours[v,s] / ss_v * ss_s
# When state steady states are 1 (as here): G_stata = diag(1/ss_controls) * G_ours
#
# Both are mathematically correct -- just different normalization.

ss_controls <- c(p = ss$values[["p"]], r = ss$values[["r"]], x = ss$values[["x"]])
G_pct <- diag(1 / ss_controls) %*% G
rownames(G_pct) <- c("p", "r", "x")
colnames(G_pct) <- colnames(G)

cat("Our policy matrix (level deviations):\n")
print(round(G, 6))

cat("\nConverted to percentage deviations (Stata convention):\n")
print(round(G_pct, 6))

# Stata results from Intro 3d
stata_G <- matrix(c(
  0.952921, 0.9678135, 1.880469,    # z column
  -1.608216, -0.4172515, 0.189277   # m column
), nrow = 3, ncol = 2)
rownames(stata_G) <- c("x", "p", "r")
colnames(stata_G) <- c("z", "m")

cat("\nStata's policy matrix (percentage deviations):\n")
print(round(stata_G, 6))

cat("\nComparison (R pct-deviation vs Stata):\n\n")
cat(sprintf("  %-8s  %12s  %12s  %12s  %12s\n",
            "", "R (z)", "Stata (z)", "R (m)", "Stata (m)"))
cat(sprintf("  %-8s  %12s  %12s  %12s  %12s\n",
            "", "------", "------", "------", "------"))

for (v in c("x", "p", "r")) {
  cat(sprintf("  %-8s  %12.6f  %12.6f  %12.6f  %12.6f\n",
              v, G_pct[v, "z"], stata_G[v, "z"], G_pct[v, "m"], stata_G[v, "m"]))
}

cat("\nMaximum absolute difference (after unit conversion):",
    round(max(abs(G_pct[c("x", "p", "r"), ] - stata_G)), 6), "\n")

# -----------------------------------------------------------------------------
# Step 5: Compare transition matrix
# -----------------------------------------------------------------------------
cat("\n--- Step 5: Transition matrix comparison ---\n\n")

H <- transition_matrix(sol, se = FALSE)
cat("Transition matrix:\n")
print(round(H, 6))

cat("\nStata transition:\n")
cat("  z -> z:", 0.9545255, "\n")
cat("  m -> m:", 0.7005489, "\n")
cat("\nMax difference:", round(max(abs(H - diag(c(0.9545255, 0.7005489)))), 8), "\n")

# -----------------------------------------------------------------------------
# Step 6: Stability check
# -----------------------------------------------------------------------------
cat("\n--- Step 6: Stability ---\n\n")
print(stability(sol))

# -----------------------------------------------------------------------------
# Step 7: IRFs (matching Stata's irf graph)
# -----------------------------------------------------------------------------
cat("\n--- Step 7: IRFs ---\n\n")
cat("Stata describes: monetary shock (m) is contractionary,\n")
cat("demand shock (z) is expansionary.\n\n")

irfs <- irf(sol, periods = 20, se = FALSE)
print(irfs)
plot(irfs)

cat("\nIRF sign check at impact (period 0):\n")
irf_data <- irfs$data

# Check signs for m shock (contractionary monetary policy)
m_p0 <- irf_data[irf_data$period == 0 & irf_data$impulse == "m", ]
cat("  m shock:\n")
for (v in c("p", "x", "r")) {
  val <- m_p0[m_p0$response == v, "value"]
  cat(sprintf("    %s = %8.4f  (%s)\n", v, val,
              ifelse(val > 0, "positive", "negative")))
}
cat("  Stata: p < 0 (inflation falls), x < 0 (output gap falls),\n")
cat("         r > 0 (interest rate rises) -- contractionary\n\n")

# Check signs for z shock (expansionary)
z_p0 <- irf_data[irf_data$period == 0 & irf_data$impulse == "z", ]
cat("  z shock:\n")
for (v in c("p", "x", "r")) {
  val <- z_p0[z_p0$response == v, "value"]
  cat(sprintf("    %s = %8.4f  (%s)\n", v, val,
              ifelse(val > 0, "positive", "negative")))
}
cat("  Stata: p > 0, x > 0, r > 0 -- expansionary\n")

cat("\n============================================================\n")
cat(" Comparison complete.\n")
cat(" If max difference is < 0.01, the R package matches Stata.\n")
cat("============================================================\n")
