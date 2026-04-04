# tests/testthat/test-ramsey-policy.R

# Simple AR(1) model: one state, one control, one instrument
make_ar1_ramsey <- function() {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )
  list(model = m,
       params   = c(rho = 0.8),
       shock_sd = c(z = 0.1))
}

# Welfare weights: penalise output variance
simple_weights <- function() {
  list(Q_yy = matrix(1, 1, 1, dimnames = list("y", "y")))
}

test_that("ramsey_policy returns correct class and fields", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())

  expect_s3_class(ram, "dsge_ramsey")
  expect_true(is.matrix(ram$F))
  expect_true(is.matrix(ram$H_ram))
  expect_true(is.matrix(ram$G_ram))
  expect_true(is.matrix(ram$P))
  expect_equal(ram$instruments, "y")
})

test_that("optimal feedback rule dimensions are correct", {
  d <- make_ar1_ramsey()
  sol <- solve_dsge(d$model, d$params, d$shock_sd)
  n_s <- ncol(sol$H)
  n_u <- 1L

  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())

  expect_equal(dim(ram$F), c(n_u, n_s))
  expect_equal(dim(ram$H_ram), c(n_s, n_s))
  expect_equal(dim(ram$G_ram), c(nrow(sol$G), n_s))
})

test_that("Riccati equation converges for stable model", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())
  expect_true(ram$converged)
})

test_that("P matrix is positive semi-definite (up to numerical noise)", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())
  ev <- eigen(ram$P, symmetric = TRUE, only.values = TRUE)$values
  # Allow small negative eigenvalues due to floating-point accumulation
  expect_true(all(ev >= -1e-6 * max(abs(ev))))
})

test_that("welfare_loss under optimal policy is finite", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())
  expect_true(is.numeric(ram$welfare_loss))
  expect_true(is.finite(ram$welfare_loss))
  # Welfare loss may be slightly negative due to floating-point noise
  expect_true(ram$welfare_loss >= -1e-6)
})

test_that("welfare_loss function works with optimal rule", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())
  wl <- welfare_loss(ram)
  expect_true(is.list(wl))
  expect_true("welfare_loss" %in% names(wl))
  expect_true(is.numeric(wl$welfare_loss))
})

test_that("welfare_loss with alternative (zero) policy", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())
  n_s <- ncol(ram$F)
  F_zero <- matrix(0, 1L, n_s)
  wl_zero <- welfare_loss(ram, F_alt = F_zero)
  expect_true(is.list(wl_zero))
})

test_that("unknown instrument errors gracefully", {
  d <- make_ar1_ramsey()
  expect_error(
    ramsey_policy(d$model, d$params, d$shock_sd,
                  instruments     = "not_a_variable",
                  welfare_weights = simple_weights()),
    "not found"
  )
})

test_that("Q_yy named vector is accepted", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = list(Q_yy = c(y = 2)))
  expect_s3_class(ram, "dsge_ramsey")
  expect_true(is.finite(ram$welfare_loss))
})

test_that("print method works", {
  d <- make_ar1_ramsey()
  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = simple_weights())
  expect_output(print(ram), "Ramsey Optimal Policy")
  expect_output(print(ram), "Feedback Rule")
  expect_output(print(ram), "beta")
})

test_that("ramsey_policy works with Q_xx weight", {
  d <- make_ar1_ramsey()
  sol <- solve_dsge(d$model, d$params, d$shock_sd)
  n_s <- ncol(sol$H)
  Q_xx_mat <- diag(0.5, n_s)
  dimnames(Q_xx_mat) <- list(colnames(sol$H), colnames(sol$H))

  ram <- ramsey_policy(d$model, d$params, d$shock_sd,
                       instruments     = "y",
                       welfare_weights = list(Q_xx = Q_xx_mat,
                                              Q_yy = matrix(1, 1, 1,
                                                dimnames = list("y","y"))))
  expect_s3_class(ram, "dsge_ramsey")
  expect_true(ram$converged)
})
