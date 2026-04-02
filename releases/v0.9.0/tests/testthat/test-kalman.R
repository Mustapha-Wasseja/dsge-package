test_that("Kalman filter produces finite log-likelihood for AR(1)", {
  # Simple AR(1): y_t = z_t, z_{t+1} = rho * z_t + e_t
  G <- matrix(1)
  H <- matrix(0.8)
  M <- matrix(1)
  D <- matrix(1)

  set.seed(42)
  n <- 100
  z <- numeric(n)
  for (i in 2:n) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  y <- matrix(z, ncol = 1)

  kf <- kalman_filter(y, G, H, M, D)

  expect_true(is.finite(kf$loglik))
  expect_true(kf$loglik < 0)
  expect_equal(nrow(kf$filtered_states), n)
  expect_equal(nrow(kf$predicted_obs), n)
})

test_that("Kalman filter handles 2-obs model", {
  G <- diag(2)
  H <- diag(c(0.7, 0.9))
  M <- diag(c(1, 0.5))
  D <- diag(2)

  set.seed(123)
  n <- 100
  z <- matrix(0, n, 2)
  for (i in 2:n) {
    z[i, ] <- H %*% z[i - 1, ] + M %*% rnorm(2)
  }
  y <- z

  kf <- kalman_filter(y, G, H, M, D)

  expect_true(is.finite(kf$loglik))
  expect_equal(ncol(kf$filtered_states), 2)
})

test_that("unconditional P solves Lyapunov equation", {
  H <- matrix(0.8)
  Q <- matrix(1)

  P <- compute_unconditional_P(H, Q)

  # Should satisfy: P = H * P * H' + Q
  check <- H %*% P %*% t(H) + Q
  expect_equal(as.numeric(P), as.numeric(check), tolerance = 1e-10)
})

test_that("Kalman smoother returns valid states", {
  G <- matrix(1)
  H <- matrix(0.8)
  M <- matrix(1)
  D <- matrix(1)

  set.seed(42)
  n <- 50
  z <- numeric(n)
  for (i in 2:n) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  y <- matrix(z, ncol = 1)

  sm <- kalman_smoother(y, G, H, M, D)

  expect_equal(nrow(sm$smoothed_states), n)
  expect_true(all(is.finite(sm$smoothed_states)))
})
