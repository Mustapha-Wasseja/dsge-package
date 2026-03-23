# =============================================================================
# Validation: Financial Frictions Model (Stata Intro 3c)
# =============================================================================
#
# This script validates the dsge package against the worked example in
# Stata 19 DSGE Reference Manual, Introduction 3c (Financial Frictions
# model). We verify the solver, policy matrix, and transition matrix at
# the documented parameter values.
#
# Source: Stata 19 DSGE Reference Manual, pp. 55-62 (Intro 3c)
# =============================================================================

library(dsge)

cat("=== Validation: Financial Frictions Model (Stata Intro 3c) ===\n\n")

# --- Model specification ---
# Equations from Stata documentation:
#   p  = beta * F.p + kappa * x            (Phillips curve)
#   x  = F.x - (i - F.p - g)              (IS curve)
#   i  = chi * r + e                       (Interest rate spread)
#   r  = psi * p + u                       (Taylor rule)
#   F.e = rhoe * e                         (Financial shock)
#   F.u = rhou * u                         (Monetary shock)
#   F.g = rhoz * g                         (Demand shock)
#
# Observed: p (inflation), i (market rate), r (safe rate)
# Unobserved: x (output gap)
# States: e (financial), u (monetary), g (demand)
# Fixed: beta = 0.96

ff <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),
  unobs(x ~ lead(x) - (i - lead(p) - g)),
  obs(i   ~ chi * r + e),
  obs(r   ~ psi * p + u),
  state(e ~ rhoe * e),
  state(u ~ rhou * u),
  state(g ~ rhoz * g),
  fixed = list(beta = 0.96),
  start = list(kappa = 0.05, chi = 0.9, psi = 6,
               rhoe = 0.85, rhou = 0.82, rhoz = 0.99)
)

cat("Model specification:\n")
print(ff)

# --- Solve at Stata's documented parameter values ---
stata_params <- c(
  beta  = 0.96,
  kappa = 0.0503257,
  chi   = 0.9067394,
  psi   = 6.332377,
  rhoe  = 0.8478222,
  rhou  = 0.815346,
  rhoz  = 0.9861866
)
stata_shock_sd <- c(e = 0.8857071, u = 7.160717, g = 0.3911862)

cat("\n--- Solving at Stata parameter values ---\n")
sol <- solve_dsge(ff, params = stata_params, shock_sd = stata_shock_sd)

if (sol$stable) {
  print(sol)

  # --- Verify policy matrix ---
  # Stata's documented policy matrix (Intro 3c):
  #   p responds to e: -0.1832608; u: -0.1584194; g: 0.2096327
  #   x responds to e: -0.6776486; u: -0.683934;  g: 0.2218595
  #   i responds to e: -0.0522495; u: -0.0028755; g: 1.203672
  #   r responds to e: -1.160476;  u: -0.0031712; g: 1.327473
  cat("\n--- Policy matrix comparison ---\n")
  G <- policy_matrix(sol, se = FALSE)
  cat("R dsge package G:\n")
  print(round(G, 7))

  stata_G <- matrix(
    c(-0.1832608, -0.6776486, -0.0522495, -1.160476,
      -0.1584194, -0.683934,  -0.0028755, -0.0031712,
       0.2096327,  0.2218595,  1.203672,   1.327473),
    nrow = 4, ncol = 3,
    dimnames = list(c("p", "x", "i", "r"), c("e", "u", "g"))
  )
  cat("\nStata documented G:\n")
  print(round(stata_G, 7))

  # Compare element-by-element using row/column names
  common_rows <- intersect(rownames(G), rownames(stata_G))
  common_cols <- intersect(colnames(G), colnames(stata_G))
  G_diffs <- numeric(0)
  for (r in common_rows) {
    for (cc in common_cols) {
      G_diffs <- c(G_diffs, abs(G[r, cc] - stata_G[r, cc]))
    }
  }
  cat("\nMax absolute difference in G (by name):", max(G_diffs), "\n")
  cat("Policy matrix match:",
      if (max(G_diffs) < 1e-4) "PASS" else "FAIL", "\n")

  # --- Verify transition matrix ---
  cat("\n--- Transition matrix comparison ---\n")
  H <- transition_matrix(sol, se = FALSE)
  cat("R dsge package H:\n")
  print(round(H, 7))

  # Stata's documented H is diagonal:
  stata_H <- diag(c(0.8478222, 0.815346, 0.9861866))
  rownames(stata_H) <- colnames(stata_H) <- c("e", "u", "g")
  cat("\nStata documented H:\n")
  print(round(stata_H, 7))

  H_diff <- abs(as.numeric(H) - as.numeric(stata_H))
  cat("\nMax absolute difference in H:", max(H_diff), "\n")
  cat("Transition matrix match:",
      if (max(H_diff) < 0.001) "PASS" else "FAIL", "\n")

  # --- Stability ---
  cat("\n--- Stability check ---\n")
  stab <- stability(sol)
  print(stab)

  # --- IRFs ---
  cat("\n--- IRFs for financial shock (first 5 periods) ---\n")
  irfs <- irf(sol, periods = 5, impulse = "e", se = FALSE)
  for (resp in c("p", "x", "i", "r")) {
    vals <- irfs$data[irfs$data$response == resp, "value"]
    cat(sprintf("  %s: %s\n", resp,
                paste(round(vals, 4), collapse = ", ")))
  }

} else {
  cat("Model is NOT stable at documented parameters.\n")
  cat("This may indicate a solver discrepancy.\n")
  cat("Eigenvalues:\n")
  print(sol$eigenvalues)
}

# --- Plots ---
if (sol$stable) {
  cat("\n--- Generating IRF plots ---\n")
  irfs_plot <- irf(sol, periods = 20, se = FALSE)
  plot(irfs_plot)
}

cat("\n=== Validation complete ===\n")
