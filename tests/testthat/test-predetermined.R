# Tests for predetermined() alias

test_that("predetermined() is a synonym for state(shock = FALSE)", {
  eq1 <- predetermined(k ~ k + i)
  eq2 <- state(k ~ k + i, shock = FALSE)
  expect_equal(eq1$type,    eq2$type)
  expect_equal(eq1$shock,   eq2$shock)
  expect_equal(eq1$formula, eq2$formula)
})


test_that("model built with predetermined() solves the same as state(shock=FALSE)", {
  m_pred <- dsge_model(
    obs(y ~ rho * y_lag + e),
    predetermined(y_lag ~ y),
    state(e ~ phi * e),
    fixed = list(rho = 0.5, phi = 0.5))

  m_state <- dsge_model(
    obs(y ~ rho * y_lag + e),
    state(y_lag ~ y, shock = FALSE),
    state(e ~ phi * e),
    fixed = list(rho = 0.5, phi = 0.5))

  sol1 <- solve_dsge(m_pred,  params = c(rho = 0.5, phi = 0.5),
                     shock_sd = c(e = 1))
  sol2 <- solve_dsge(m_state, params = c(rho = 0.5, phi = 0.5),
                     shock_sd = c(e = 1))

  expect_equal(sol1$G, sol2$G, tolerance = 1e-10)
  expect_equal(sol1$H, sol2$H, tolerance = 1e-10)
})


test_that("predetermined() errors on non-formula input", {
  expect_error(predetermined("k ~ k + i"), "requires a formula")
  expect_error(predetermined(42),          "requires a formula")
})
