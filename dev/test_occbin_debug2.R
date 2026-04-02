devtools::load_all(".")

m <- dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.9))
sol <- solve_dsge(m, params = c(rho = 0.9), shock_sd = c(z = 1.0))

H <- sol$H; G <- sol$G; M <- sol$M
n_s <- 1; n_c <- 1; n_shocks <- 1
horizon <- 10

shock_path <- matrix(0, horizon, 1)
shock_path[1, 1] <- -1.0
x0 <- 0

# Simulate unconstrained
sim <- dsge:::.forward_simulate(H, G, M, x0, shock_path, horizon)
cat("Unconstrained controls:\n")
print(round(sim$controls[,1], 4))

# Check constraint info
bound_dev <- -0.5
impact <- as.numeric(G[1, ] %*% M[, 1])
cat("Impact coef:", impact, "\n")
cat("Bound dev:", bound_dev, "\n")

# Manual check
y_vals <- sim$controls[, 1]
for (t in 1:horizon) {
  violated <- y_vals[t] < bound_dev - 1e-8
  cat("t=", t, " y=", round(y_vals[t],4), " violated=", violated, "\n")
}
