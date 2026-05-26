# Tests for extended_path()

make_linear <- function() {
  m <- dsge_model(
    obs(y ~ rho * y_lag + e),
    state(y_lag ~ y, shock = FALSE),
    state(e ~ phi * e),
    fixed = list(rho = 0.5, phi = 0.5))
  solve_dsge(m, params = c(rho = 0.5, phi = 0.5), shock_sd = c(e = 1))
}

test_that("extended_path on a linear model returns expected shape", {
  sol <- make_linear()
  sim <- extended_path(sol, periods = 50L, seed = 1L)
  expect_s3_class(sim, "dsge_extended_path")
  expect_equal(nrow(sim$states),   50L)
  expect_equal(nrow(sim$controls), 50L)
  expect_equal(nrow(sim$shocks),   50L)
  expect_false(sim$nonlinear)
})

test_that("burn-in trims initial periods", {
  sol <- make_linear()
  sim <- extended_path(sol, periods = 30L, burn = 20L, seed = 2L)
  expect_equal(nrow(sim$states), 30L)
})

test_that("results are reproducible with the same seed", {
  sol <- make_linear()
  sim1 <- extended_path(sol, periods = 20L, seed = 7L)
  sim2 <- extended_path(sol, periods = 20L, seed = 7L)
  expect_identical(sim1$shocks, sim2$shocks)
  expect_identical(sim1$states, sim2$states)
})

test_that("nonzero variance in simulated controls", {
  sol <- make_linear()
  sim <- extended_path(sol, periods = 200L, seed = 1L)
  expect_true(var(sim$controls[, "y"]) > 0)
})

test_that("print.dsge_extended_path works", {
  sol <- make_linear()
  sim <- extended_path(sol, periods = 10L, seed = 1L)
  expect_output(print(sim), "Extended-path simulation")
})

test_that("requires params + shock_sd for dsgenl_model", {
  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C", endo_state = "K", exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9))
  expect_error(extended_path(rbc, periods = 5L), "params.*shock_sd")
})
