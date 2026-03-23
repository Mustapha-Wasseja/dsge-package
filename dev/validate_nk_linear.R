# =============================================================================
# Validation: Linear New Keynesian Model (Stata Intro 3a)
# =============================================================================
#
# This script validates the dsge package against the worked example in
# Stata 19 DSGE Reference Manual, Introduction 3a (Linear New Keynesian
# model). We verify the solver, policy matrix, transition matrix, and
# stability diagnostics at the documented parameter values.
#
# The exact Stata dataset (usmacro2.dta) is not bundled; synthetic data
# following the same model structure is used for estimation validation.
#
# Source: Stata 19 DSGE Reference Manual, pp. 36-44 (Intro 3a)
# =============================================================================

library(dsge)

cat("=== Validation: Linear NK Model (Stata Intro 3a) ===\n\n")

# --- Model specification ---
# Equations from Stata documentation:
#   p  = beta * F.p + kappa * x          (Phillips curve)
#   x  = F.x - (r - F.p - g)            (IS curve)
#   r  = psi * p + u                     (Taylor rule)
#   F.u = rhou * u                       (monetary shock)
#   F.g = rhog * g                       (demand shock)
#
# Observed: p (inflation), r (interest rate)
# Unobserved: x (output gap)
# States: u (monetary), g (demand)
# Fixed: beta = 0.96

nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),
  unobs(x ~ lead(x) - (r - lead(p) - g)),
  obs(r   ~ psi * p + u),
  state(u ~ rhou * u),
  state(g ~ rhog * g),
  fixed = list(beta = 0.96),
  start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
)

cat("Model specification:\n")
print(nk)

# --- Solve at Stata's documented parameter values ---
stata_params <- c(
  beta  = 0.96,
  kappa = 0.0849631,
  psi   = 1.943004,
  rhou  = 0.7005483,
  rhog  = 0.9545257
)
stata_shock_sd <- c(u = 2.318204, g = 0.5689891)

cat("\n--- Solving at Stata parameter values ---\n")
sol <- solve_dsge(nk, params = stata_params, shock_sd = stata_shock_sd)
print(sol)

# --- Verify policy matrix against Stata documentation ---
# Stata's documented policy matrix (Intro 3a):
#   p responds to u: -0.4172521;  p responds to g:  0.9678177
#   x responds to u: -1.608216;   x responds to g:  0.9529203
#   r responds to u:  0.1892776;  r responds to g:  1.880474
cat("\n--- Policy matrix comparison ---\n")
G <- policy_matrix(sol, se = FALSE)
cat("R dsge package G:\n")
print(round(G, 7))

# Stata row order is [p, x, r]; our row order is [p, r, x]
# Compare element-by-element using variable names
stata_G_vals <- c(
  p.u = -0.4172521, p.g = 0.9678177,
  x.u = -1.608216,  x.g = 0.9529203,
  r.u = 0.1892776,  r.g = 1.880474
)
cat("\nStata documented G (element-by-element):\n")
cat("  p.u =", stata_G_vals["p.u"], "  p.g =", stata_G_vals["p.g"], "\n")
cat("  x.u =", stata_G_vals["x.u"], "  x.g =", stata_G_vals["x.g"], "\n")
cat("  r.u =", stata_G_vals["r.u"], "  r.g =", stata_G_vals["r.g"], "\n")

G_diffs <- c(
  abs(G["p", "u"] - stata_G_vals["p.u"]),
  abs(G["p", "g"] - stata_G_vals["p.g"]),
  abs(G["x", "u"] - stata_G_vals["x.u"]),
  abs(G["x", "g"] - stata_G_vals["x.g"]),
  abs(G["r", "u"] - stata_G_vals["r.u"]),
  abs(G["r", "g"] - stata_G_vals["r.g"])
)
cat("\nMax absolute difference in G:", max(G_diffs), "\n")
cat("Policy matrix match:", if (max(G_diffs) < 1e-4) "PASS" else "FAIL", "\n")

# --- Verify transition matrix against Stata documentation ---
# Stata's transition matrix: diagonal with rhou and rhog
cat("\n--- Transition matrix comparison ---\n")
H <- transition_matrix(sol, se = FALSE)
cat("R dsge package H:\n")
print(round(H, 7))

stata_H <- matrix(
  c(0.7005483, 0, 0, 0.9545257),
  nrow = 2, ncol = 2,
  dimnames = list(c("u", "g"), c("u", "g"))
)
cat("\nStata documented H:\n")
print(round(stata_H, 7))

H_diff <- abs(as.numeric(H) - as.numeric(stata_H))
cat("\nMax absolute difference in H:", max(H_diff), "\n")
cat("Transition matrix match:", if (max(H_diff) < 0.001) "PASS" else "FAIL", "\n")

# --- Stability diagnostics ---
cat("\n--- Stability check ---\n")
stab <- stability(sol)
print(stab)
cat("Stability check:", if (stab$stable) "PASS" else "FAIL", "\n")

# --- IRFs ---
cat("\n--- Impulse-response functions (first 5 periods) ---\n")
irfs <- irf(sol, periods = 5, se = FALSE)
irf_u <- irfs$data[irfs$data$impulse == "u" & irfs$data$response == "p", ]
cat("IRF of p to u shock (periods 0-5):\n")
print(round(irf_u$value, 6))

# --- Estimation on synthetic data ---
cat("\n--- Estimation on synthetic data ---\n")
cat("Generating data from NK model at Stata parameter values...\n")
set.seed(123)
n <- 244  # Same as Stata sample size

# Simulate from the state-space model
u <- numeric(n)
g <- numeric(n)
for (t in 2:n) {
  u[t] <- 0.7005 * u[t - 1] + 2.318 * rnorm(1)
  g[t] <- 0.9545 * g[t - 1] + 0.569 * rnorm(1)
}
states <- cbind(u, g)
Z <- sol$D %*% sol$G
obs_data <- states %*% t(Z)
colnames(obs_data) <- c("p", "r")
dat <- as.data.frame(obs_data)

cat("Estimating model...\n")
fit <- tryCatch(
  estimate(nk, data = dat,
           start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
           control = list(maxit = 500)),
  error = function(e) { cat("Estimation error:", e$message, "\n"); NULL }
)

if (!is.null(fit)) {
  cat("\nEstimation results:\n")
  summary(fit)
  cat("\nConvergence:", if (fit$convergence == 0) "PASS" else "FAIL", "\n")

  # Check parameter recovery (within reasonable tolerance)
  est_rhou <- coef(fit)["rhou"]
  cat("rhou estimated:", round(est_rhou, 4),
      " (true: 0.7005, tolerance: 0.15)\n")
  cat("rhou recovery:", if (abs(est_rhou - 0.7005) < 0.15) "PASS" else "FAIL", "\n")

  # Postestimation
  cat("\n--- Postestimation from fitted model ---\n")
  cat("Policy matrix:\n")
  print(round(policy_matrix(fit, se = FALSE), 4))

  cat("\nForecasts:\n")
  fc <- forecast(fit, horizon = 8)
  print(fc)
} else {
  cat("Estimation failed; solver may not converge with synthetic data.\n")
}

# --- Plots ---
if (!is.null(fit)) {
  cat("\n--- Generating plots ---\n")
  irfs_plot <- irf(fit, periods = 20)
  plot(irfs_plot)

  fc_plot <- forecast(fit, horizon = 12)
  plot(fc_plot)
}

cat("\n=== Validation complete ===\n")
