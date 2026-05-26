# Tests for perfect_foresight_expect_err()

make_ar1 <- function() {
  m <- dsge_model(
    obs(y ~ rho * y_lag + e),
    state(y_lag ~ y, shock = FALSE),
    state(e ~ phi * e),
    fixed = list(rho = 0.5, phi = 0.5))
  solve_dsge(m, params = c(rho = 0.5, phi = 0.5), shock_sd = c(e = 1))
}


test_that("returns expected class and components", {
  sol <- make_ar1()
  pf <- perfect_foresight_expect_err(sol,
    actual_shocks = list(e = c(0, 0, 1)),
    horizon = 10)
  expect_s3_class(pf, "dsge_perfect_foresight_expecterr")
  expect_s3_class(pf, "dsge_perfect_foresight")
  expect_true("expectation_paths" %in% names(pf))
  expect_equal(length(pf$expectation_paths), 10L)
  expect_equal(dim(pf$states)[1], 10L)
})


test_that("when expected == actual, result matches perfect_foresight()", {
  sol <- make_ar1()
  shk <- list(e = c(0.5, 0.3, 0.2))
  pf_pf <- perfect_foresight(sol, shocks = shk, horizon = 20L)
  pf_ee <- perfect_foresight_expect_err(sol,
    actual_shocks   = shk,
    expected_shocks = shk,
    horizon = 20L)
  # In a backward-looking AR(1) the realised paths should match.
  expect_equal(pf_ee$states[, "e"], pf_pf$states[, "e"],
               tolerance = 1e-8)
})


test_that("unanticipated shock at period k moves state at k+1, not before", {
  sol <- make_ar1()
  pf <- perfect_foresight_expect_err(sol,
    actual_shocks = list(e = c(0, 0, 0, 1, 0)),
    horizon = 10)
  # Periods 1..3: no shock arrived yet, no propagation
  expect_lt(max(abs(pf$states[1:3, "e"])), 1e-12)
  # Period 4: shock arrives, state jumps to 1
  expect_equal(as.numeric(pf$states[4, "e"]), 1, tolerance = 1e-8)
})


test_that("expectation_paths length matches horizon and shapes are right", {
  sol <- make_ar1()
  pf <- perfect_foresight_expect_err(sol,
    actual_shocks = list(e = 1),
    horizon = 5)
  expect_equal(length(pf$expectation_paths), 5L)
  expect_equal(dim(pf$expectation_paths[[1]]$state),
               c(5L, length(pf$state_names)))
  expect_equal(dim(pf$expectation_paths[[3]]$state),
               c(3L, length(pf$state_names)))
})


test_that("errors on unknown shock name", {
  sol <- make_ar1()
  expect_error(
    perfect_foresight_expect_err(sol,
      actual_shocks = list(not_a_shock = 1), horizon = 5),
    "Unknown shock"
  )
})
