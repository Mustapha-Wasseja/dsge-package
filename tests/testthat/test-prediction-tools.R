# tests/testthat/test-prediction-tools.R

# Helper: create and estimate a simple AR(1) model
make_ar1_fit <- function() {
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

  estimate(m, data = dat)
}

test_that("fitted.dsge_fit returns filtered values", {
  fit <- make_ar1_fit()
  fv <- fitted(fit)

  expect_true(is.matrix(fv))
  expect_equal(ncol(fv), 1)
  expect_equal(nrow(fv), fit$nobs)
  expect_equal(colnames(fv), "y")
})

test_that("fitted values differ from one-step predictions", {
  fit <- make_ar1_fit()
  fv <- fitted(fit)
  pred <- predict(fit, type = "observed", method = "onestep")

  # Filtered and one-step should differ (filtered uses current obs)
  expect_false(all(abs(fv - pred) < 1e-12))
})

test_that("prediction_accuracy returns sensible results", {
  fit <- make_ar1_fit()
  acc <- prediction_accuracy(fit)

  expect_s3_class(acc, "dsge_prediction_accuracy")
  expect_true(all(acc$rmse > 0))
  expect_true(all(acc$mae > 0))
  expect_equal(length(acc$rmse), 1)
  expect_equal(acc$n, fit$nobs)

  # RMSE >= MAE (by Jensen's inequality)
  expect_true(all(acc$rmse >= acc$mae))
})

test_that("prediction_accuracy RMSE is correct", {
  fit <- make_ar1_fit()
  res <- residuals(fit)
  acc <- prediction_accuracy(fit)

  expected_rmse <- sqrt(mean(res[, 1]^2))
  expect_equal(acc$rmse[[1]], expected_rmse, tolerance = 1e-10)
})

test_that("prediction_accuracy print method works", {
  fit <- make_ar1_fit()
  acc <- prediction_accuracy(fit)
  expect_output(print(acc), "Prediction Accuracy")
  expect_output(print(acc), "RMSE")
})

test_that("prediction_interval returns correct structure", {
  fit <- make_ar1_fit()
  pi <- prediction_interval(fit, level = 0.95)

  expect_s3_class(pi, "dsge_prediction_interval")
  expect_equal(dim(pi$fit), c(fit$nobs, 1))
  expect_equal(dim(pi$se), c(fit$nobs, 1))
  expect_equal(dim(pi$lower), c(fit$nobs, 1))
  expect_equal(dim(pi$upper), c(fit$nobs, 1))
  expect_equal(pi$level, 0.95)

  # Lower < fit < upper
  expect_true(all(pi$lower <= pi$fit + 1e-10))
  expect_true(all(pi$upper >= pi$fit - 1e-10))
})

test_that("prediction_interval SEs are positive", {
  fit <- make_ar1_fit()
  pi <- prediction_interval(fit)

  expect_true(all(pi$se > 0, na.rm = TRUE))
})

test_that("prediction_interval width increases with level", {
  fit <- make_ar1_fit()
  pi90 <- prediction_interval(fit, level = 0.90)
  pi99 <- prediction_interval(fit, level = 0.99)

  # 99% interval should be wider than 90%
  width90 <- mean(pi90$upper - pi90$lower)
  width99 <- mean(pi99$upper - pi99$lower)
  expect_true(width99 > width90)
})

test_that("prediction_interval print method works", {
  fit <- make_ar1_fit()
  pi <- prediction_interval(fit)
  expect_output(print(pi), "Prediction Intervals")
  expect_output(print(pi), "Mean prediction SE")
})

test_that("residuals dimensions match predictions", {
  fit <- make_ar1_fit()
  res <- residuals(fit)
  pred <- predict(fit, type = "observed", method = "onestep")

  expect_equal(dim(res), dim(pred))
})
