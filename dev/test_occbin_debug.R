devtools::load_all(".")

sol <- build_ar1_sol <- (function() {
  m <- dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.9))
  solve_dsge(m, params = c(rho = 0.9), shock_sd = c(z = 1.0))
})()

cat("H:", sol$H, "\nG:", sol$G, "\nM:", sol$M, "\n")

# Unconstrained path with negative shock
pf <- perfect_foresight(sol, shocks = list(z = -1.0), horizon = 10)
cat("\nUnconstrained y:\n")
print(round(pf$controls[,1], 4))

# Now try OBC
obc <- simulate_occbin(sol,
  constraints = list("y >= -0.5"),
  shocks = list(z = -1.0),
  horizon = 10)

cat("\nShadow shocks:", round(obc$shadow_shocks[,1], 4), "\n")
cat("Binding:", obc$binding[,1], "\n")
cat("Constrained y:\n")
print(round(obc$controls[,1], 4))
cat("Unconstrained y:\n")
print(round(obc$controls_unc[,1], 4))
cat("Converged:", obc$converged, "\n")
cat("Iterations:", obc$n_iter, "\n")
