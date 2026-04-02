devtools::load_all(".")

m <- dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.9))
sol <- solve_dsge(m, params = c(rho = 0.9), shock_sd = c(z = 1.0))

# Manually run the OccBin algorithm
esol <- dsge:::.extract_solution(sol)
H <- esol$H; G <- esol$G; M <- esol$M

horizon <- 10
n_s <- 1; n_c <- 1; n_shocks <- 1

shock_path <- matrix(0, horizon, 1)
shock_path[1, 1] <- -1.0
x0 <- 0

# Constraint info
bound_dev <- -0.5
var_idx <- 1
shock_idx <- 1
impact <- as.numeric(G[var_idx, ] %*% M[, shock_idx])
tol <- 1e-8

cat("Impact:", impact, "\n")

shadow <- matrix(0, horizon, 1)
binding_prev <- matrix(FALSE, horizon, 1)

for (iter in 1:5) {
  cat("\n=== Iteration", iter, "===\n")

  aug_shock <- shock_path
  aug_shock[, shock_idx] <- aug_shock[, shock_idx] + shadow[, 1]
  cat("Aug shock:", round(aug_shock[,1], 4), "\n")

  sim <- dsge:::.forward_simulate(H, G, M, x0, aug_shock, horizon)
  cat("Controls:", round(sim$controls[,1], 4), "\n")

  binding_new <- matrix(FALSE, horizon, 1)
  shadow_new <- matrix(0, horizon, 1)

  y_vals <- sim$controls[, var_idx]
  for (t in 1:horizon) {
    violated <- y_vals[t] < bound_dev - tol
    if (violated) {
      binding_new[t, 1] <- TRUE
      shadow_new[t, 1] <- (bound_dev - y_vals[t]) / impact
    }
  }

  cat("Binding:", binding_new[,1], "\n")
  cat("Shadow:", round(shadow_new[,1], 4), "\n")

  if (identical(binding_new, binding_prev)) {
    cat("CONVERGED at iter", iter, "\n")
    break
  }

  binding_prev <- binding_new
  shadow <- shadow_new
}
