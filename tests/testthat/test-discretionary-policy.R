# Tests for discretionary_policy()

make_nk <- function() {
  dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
}


test_that("discretionary_policy returns expected fields", {
  m <- make_nk()
  out <- discretionary_policy(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1, e.g = 0.5),
    instruments = "r",
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5)),
    beta = 0.99)
  expect_s3_class(out, "dsge_discretionary")
  expect_true("F" %in% names(out))
  expect_true(out$converged)
  expect_true(is.finite(out$welfare_loss))
})


test_that("discretionary_policy iteration delivers a stabilising rule", {
  m <- make_nk()
  out <- discretionary_policy(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1, e.g = 0.5),
    instruments = "r",
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.05)))
  ev <- abs(eigen(out$H_ram, only.values = TRUE)$values)
  expect_true(max(ev) < 1)
})


test_that("discretion welfare loss is finite and non-negative", {
  m <- make_nk()
  out <- discretionary_policy(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1, e.g = 0.5),
    instruments = "r",
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.05)))
  expect_true(is.finite(out$welfare_loss))
  expect_true(out$welfare_loss >= 0)
})


test_that("print method runs", {
  m <- make_nk()
  out <- discretionary_policy(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1, e.g = 0.5),
    instruments = "r",
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5)))
  expect_output(print(out), "Discretionary")
})


test_that("invalid instrument errors", {
  m <- make_nk()
  expect_error(
    discretionary_policy(m,
      params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
      shock_sd    = c(e.u = 1, e.g = 0.5),
      instruments = "not_a_var",
      welfare_weights = list(Q_yy = c(p = 1))),
    "Instruments must be"
  )
})
