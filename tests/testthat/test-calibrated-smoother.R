# Tests for calibrated_smoother() and dsge_solution smoother methods

make_setup <- function(TT = 100L, seed = 1L) {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    fixed = list(rho = 0.8)
  )
  sol <- solve_dsge(m, params = c(rho = 0.8), shock_sd = c(z = 1.0))
  set.seed(seed)
  e <- rnorm(TT)
  z_path <- numeric(TT)
  for (i in 2:TT) z_path[i] <- 0.8 * z_path[i - 1] + e[i]
  list(sol = sol, data = data.frame(y = z_path))
}


test_that("calibrated_smoother returns expected components", {
  s <- make_setup()
  out <- calibrated_smoother(s$sol, data = s$data)
  expect_s3_class(out, "dsge_calibrated_smoother")
  expect_named(out, c("states", "shocks", "decomposition"))
  expect_s3_class(out$states,        "dsge_smoothed")
  expect_s3_class(out$shocks,        "dsge_smoothed_shocks")
  expect_s3_class(out$decomposition, "dsge_decomposition")
})


test_that("calibrated_smoother subset via `what` works", {
  s <- make_setup()
  out_only_states <- calibrated_smoother(s$sol, data = s$data, what = "states")
  expect_named(out_only_states, "states")
})


test_that("smooth_states.dsge_solution reproduces estimate-based smoothing", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )
  set.seed(42)
  e <- rnorm(150)
  z_path <- numeric(150)
  for (i in 2:150) z_path[i] <- 0.8 * z_path[i - 1] + e[i]
  dat <- data.frame(y = z_path)

  # Calibrated smoother at the truth
  sol <- solve_dsge(m, params = c(rho = 0.8), shock_sd = c(z = 1.0))
  sm_cal <- smooth_states(sol, data = dat)

  # Estimate-based smoother (estimated at the truth too)
  fit <- estimate(m, data = dat)
  sm_est <- smooth_states(fit)

  # Smoothed states should be highly correlated (different estimation
  # paths but same model + data)
  expect_true(cor(sm_cal$smoothed_states[, 1],
                  sm_est$smoothed_states[, 1]) > 0.95)
})


test_that("calibrated_smoother errors if data lacks observables", {
  s <- make_setup()
  bad <- s$data
  colnames(bad) <- "foo"
  expect_error(calibrated_smoother(s$sol, data = bad),
               "not in model observables")
})


test_that("calibrated_smoother errors when data missing", {
  s <- make_setup()
  expect_error(smooth_states(s$sol), "required")
})


test_that("print method runs", {
  s <- make_setup()
  out <- calibrated_smoother(s$sol, data = s$data, what = "states")
  expect_output(print(out), "Calibrated DSGE smoother results")
})
