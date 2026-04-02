devtools::load_all(".")

m <- dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.9))
sol <- solve_dsge(m, params = c(rho = 0.9), shock_sd = c(z = 1.0))

# Check what .extract_solution returns
esol <- dsge:::.extract_solution(sol)
cat("H:", esol$H, "\n")
cat("G:", esol$G, "\n")
cat("M:", esol$M, "\n")
cat("SS:", esol$steady_state, "\n")
cat("shock_sd:", esol$shock_sd, "\n")

# Check con_info
control_names <- rownames(esol$G)
state_names <- rownames(esol$H)
shock_names <- colnames(esol$M)
cat("control_names:", control_names, "\n")
cat("state_names:", state_names, "\n")
cat("shock_names:", shock_names, "\n")

# Parse constraint
con <- dsge:::.parse_constraint_string("y >= -0.5")
cat("Constraint: ", con$variable, con$type, con$bound, "\n")

var_idx <- match(con$variable, control_names)
cat("var_idx:", var_idx, "\n")

# bound_dev: no SS for linear model
ss <- esol$steady_state
cat("SS is null:", is.null(ss), "\n")
if (!is.null(ss)) cat("SS:", ss, "\n")

# This is the issue: for linear models, ss might be non-null but y might not be in ss
bound_dev <- con$bound
if (!is.null(ss) && con$variable %in% names(ss)) {
  cat("Adjusting bound by SS value:", ss[con$variable], "\n")
  bound_dev <- con$bound - ss[con$variable]
}
cat("bound_dev:", bound_dev, "\n")
