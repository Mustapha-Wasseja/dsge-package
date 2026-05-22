# Tests for osr() -- Optimal Simple Rules

make_nk <- function() {
  dsge_model(
    obs(p   ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r   ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )
}


test_that("osr returns dsge_osr object with expected fields", {
  m <- make_nk()
  res <- osr(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1.0, e.g = 0.5),
    osr_params  = c(psi = 1.5),
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.1)),
    lower = 1.01, upper = 5.0)

  expect_s3_class(res, "dsge_osr")
  expect_named(res$optimal, "psi")
  expect_true(is.finite(res$loss))
  expect_true(is.finite(res$loss_at_start))
  expect_true(res$converged || is.finite(res$loss))
})


test_that("osr improves on the starting value (loss decreases)", {
  m <- make_nk()
  res <- osr(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1.0, e.g = 0.5),
    osr_params  = c(psi = 1.5),
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.1)),
    lower = 1.01, upper = 5.0)

  # Should reduce loss vs starting value (or at least not worsen it)
  expect_lte(res$loss, res$loss_at_start + 1e-8)
})


test_that("osr respects bounds", {
  m <- make_nk()
  res <- osr(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1.0, e.g = 0.5),
    osr_params  = c(psi = 2.0),
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.1)),
    lower = 1.5, upper = 3.0)

  expect_true(res$optimal["psi"] >= 1.5 - 1e-8)
  expect_true(res$optimal["psi"] <= 3.0 + 1e-8)
})


test_that("osr optimises multiple parameters", {
  # Extend the model with a response to the output gap
  m <- dsge_model(
    obs(p   ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r   ~ psi_p * p + psi_x * x + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi_p = 1.5, psi_x = 0.5,
                 rhou = 0.7, rhog = 0.9)
  )
  res <- osr(m,
    params      = c(kappa = 0.1, psi_p = 1.5, psi_x = 0.5,
                    rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1.0, e.g = 0.5),
    osr_params  = c(psi_p = 1.5, psi_x = 0.5),
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.1)),
    lower = c(1.01, 0.0),
    upper = c(5.0,  3.0))

  expect_named(res$optimal, c("psi_p", "psi_x"))
  expect_lte(res$loss, res$loss_at_start + 1e-8)
})


test_that("osr errors on unknown parameter name", {
  m <- make_nk()
  expect_error(
    osr(m,
      params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
      shock_sd    = c(e.u = 1.0, e.g = 0.5),
      osr_params  = c(not_a_param = 1.0),
      welfare_weights = list(Q_yy = c(p = 1))),
    "not found in 'params'"
  )
})


test_that("print.dsge_osr produces expected output", {
  m <- make_nk()
  res <- osr(m,
    params      = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd    = c(e.u = 1.0, e.g = 0.5),
    osr_params  = c(psi = 1.5),
    welfare_weights = list(Q_yy = c(p = 1, x = 0.5, r = 0.1)),
    lower = 1.01, upper = 5.0)
  expect_output(print(res), "Optimal Simple Rule")
  expect_output(print(res), "psi")
})
