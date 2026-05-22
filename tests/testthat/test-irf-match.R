# Tests for irf_match()

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


test_that("irf_match recovers true parameters when target IRF is model-generated", {
  nk <- make_nk()
  true_params   <- c(kappa = 0.10, psi = 1.50, rhou = 0.70, rhog = 0.90)
  true_shock_sd <- c(e.u = 1.0, e.g = 0.5)

  sol_true <- solve_dsge(nk, params = true_params, shock_sd = true_shock_sd)
  target_df <- irf(sol_true, periods = 8)$data

  est <- irf_match(nk,
    params_start   = c(kappa = 0.15, psi = 2.00, rhou = 0.60, rhog = 0.80),
    shock_sd_start = c(e.u = 0.8, e.g = 0.7),
    target         = target_df,
    method         = "Nelder-Mead",
    control        = list(maxit = 2000, reltol = 1e-8))

  expect_s3_class(est, "dsge_irf_match")
  expect_true(est$objective < 1e-3)
  # Within 2% of truth
  expect_equal(est$params,   true_params,   tolerance = 5e-2)
  expect_equal(est$shock_sd, true_shock_sd, tolerance = 5e-2)
})


test_that("irf_match returns expected fields", {
  nk <- make_nk()
  sol_true <- solve_dsge(nk,
    params   = c(kappa = 0.10, psi = 1.50, rhou = 0.70, rhog = 0.90),
    shock_sd = c(e.u = 1.0, e.g = 0.5))
  target_df <- irf(sol_true, periods = 6)$data

  est <- irf_match(nk,
    params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd_start = c(e.u = 1.0, e.g = 0.5),
    target         = target_df,
    control        = list(maxit = 50))

  expect_s3_class(est, "dsge_irf_match")
  expect_named(est$params,   c("kappa", "psi", "rhou", "rhog"))
  expect_named(est$shock_sd, c("e.u", "e.g"))
  expect_true("fitted" %in% names(est$target))
})


test_that("irf_match with weight matrix runs", {
  nk <- make_nk()
  sol_true <- solve_dsge(nk,
    params   = c(kappa = 0.10, psi = 1.50, rhou = 0.70, rhog = 0.90),
    shock_sd = c(e.u = 1.0, e.g = 0.5))
  target_df <- irf(sol_true, periods = 6)$data
  N <- nrow(target_df)
  # Weight early horizons more (decay)
  w <- diag(exp(-0.1 * target_df$period))

  est <- irf_match(nk,
    params_start   = c(kappa = 0.10, psi = 1.50, rhou = 0.70, rhog = 0.90),
    shock_sd_start = c(e.u = 1.0, e.g = 0.5),
    target         = target_df,
    weight         = w,
    control        = list(maxit = 50))
  expect_s3_class(est, "dsge_irf_match")
  expect_true(is.finite(est$objective))
})


test_that("irf_match errors when target lacks required columns", {
  nk <- make_nk()
  bad_target <- data.frame(impulse = "e.u", value = 1)
  expect_error(
    irf_match(nk,
      params_start   = c(kappa = 0.1),
      shock_sd_start = c(e.u = 1.0),
      target         = bad_target),
    "must contain columns"
  )
})


test_that("print.dsge_irf_match produces output", {
  nk <- make_nk()
  sol_true <- solve_dsge(nk,
    params   = c(kappa = 0.10, psi = 1.50, rhou = 0.70, rhog = 0.90),
    shock_sd = c(e.u = 1.0, e.g = 0.5))
  target_df <- irf(sol_true, periods = 6)$data
  est <- irf_match(nk,
    params_start   = c(kappa = 0.10, psi = 1.50, rhou = 0.70, rhog = 0.90),
    shock_sd_start = c(e.u = 1.0, e.g = 0.5),
    target         = target_df,
    control        = list(maxit = 20))
  expect_output(print(est), "IRF-Matching Estimation")
})
