devtools::load_all(".")

mod <- dsgenl_model(
  "Y = A * H",
  "1 / C = beta * R / C(+1)",
  "Y = C",
  "R = 1 / beta",
  "A(+1) = A^rhoA",
  observed = "Y",
  unobserved = c("C", "H", "R"),
  exo_state = "A",
  fixed = list(beta = 0.99),
  start = list(rhoA = 0.90),
  ss_guess = c(Y = 1, C = 1, H = 1, R = 1.0101, A = 1)
)
sol <- solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
                  shock_sd = c(A = 0.01))

cat("G:\n"); print(sol$G)
cat("M:\n"); print(sol$M)
cat("G %*% M:\n"); print(sol$G %*% sol$M)
cat("Control names:", rownames(sol$G), "\n")
cat("State names:", rownames(sol$H), "\n")
cat("SS:", sol$steady_state, "\n")
