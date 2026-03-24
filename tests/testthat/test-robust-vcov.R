# tests/testthat/test-robust-vcov.R

# Helper: create and estimate a simple AR(1) model
make_ar1_fit_for_robust <- function() {
  m <- dsge_model(
    obs(y ~ lead(y) + u),
    state(u ~ rho_u * u),
    fixed = list(),
    start = list(rho_u = 0.5)
  )

  set.seed(42)
  n <- 100
  u <- numeric(n)
  for (t in 2:n) u[t] <- 0.5 * u[t - 1] + rnorm(1, sd = 0.1)
  dat <- data.frame(y = u)

  estimate(m, data = dat, hessian = TRUE)
}

test_that("robust_vcov returns correct structure", {
  fit <- make_ar1_fit_for_robust()
  rv <- robust_vcov(fit)

  expect_s3_class(rv, "dsge_robust_vcov")
  expect_true(is.matrix(rv$vcov))

  # Square matrix with correct dimensions
  n_params <- length(rv$param_names)
  expect_equal(dim(rv$vcov), c(n_params, n_params))

  # Symmetric
  expect_equal(rv$vcov, t(rv$vcov), tolerance = 1e-10)

  # Positive diagonal (SEs^2 > 0)
  expect_true(all(diag(rv$vcov) >= 0))

  # SE vector matches diagonal
  expect_equal(rv$se^2, diag(rv$vcov), tolerance = 1e-10)
})

test_that("robust SEs are similar to conventional for well-specified model", {
  fit <- make_ar1_fit_for_robust()
  rv <- robust_vcov(fit)

  # For a correctly specified model, robust and conventional SEs
  # should be similar (ratio near 1)
  ratio <- rv$se / rv$se_conventional
  expect_true(all(ratio > 0.3 & ratio < 3.0))
})

test_that("vcov with type = 'robust' returns robust vcov", {
  fit <- make_ar1_fit_for_robust()
  v_conv <- vcov(fit, type = "conventional")
  v_robust <- vcov(fit, type = "robust")

  # Both should be matrices of same dimensions
  expect_equal(dim(v_conv), dim(v_robust))

  # They should differ
  expect_false(all(abs(v_conv - v_robust) < 1e-15))
})

test_that("vcov default is conventional", {
  fit <- make_ar1_fit_for_robust()
  v1 <- vcov(fit)
  v2 <- vcov(fit, type = "conventional")
  expect_equal(v1, v2)
})

test_that("robust_vcov print method works", {
  fit <- make_ar1_fit_for_robust()
  rv <- robust_vcov(fit)
  expect_output(print(rv), "Robust")
  expect_output(print(rv), "Conventional")
  expect_output(print(rv), "Ratio")
})

test_that("robust_vcov has correct parameter names", {
  fit <- make_ar1_fit_for_robust()
  rv <- robust_vcov(fit)

  expect_true("rho_u" %in% rv$param_names)
  expect_true(any(grepl("sd\\(e\\.", rv$param_names)))
})
