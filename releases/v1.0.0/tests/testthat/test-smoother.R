# Tests for smooth_states(), smooth_shocks(), shock_decomposition()

test_that("smooth_states works on ML fit", {
  # Simple AR(1) model
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  e <- rnorm(100)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i - 1] + e[i]
  fit <- estimate(m, data = data.frame(y = z))

  sm <- smooth_states(fit)

  expect_s3_class(sm, "dsge_smoothed")
  expect_equal(nrow(sm$smoothed_states), 100)
  expect_equal(ncol(sm$smoothed_states), 1)  # 1 state
  expect_equal(nrow(sm$smoothed_obs), 100)
  expect_equal(ncol(sm$smoothed_obs), 1)  # 1 observable
  expect_equal(nrow(sm$residuals), 100)

  # Smoothed obs should approximate the data
  expect_true(cor(sm$smoothed_obs[, 1], z) > 0.95)
})

test_that("smooth_states dimensions match for multi-state model", {
  m <- dsge_model(
    obs(y1 ~ a),
    obs(y2 ~ b),
    state(a ~ rho_a * a),
    state(b ~ rho_b * b),
    start = list(rho_a = 0.5, rho_b = 0.3)
  )

  set.seed(42)
  e1 <- rnorm(80)
  e2 <- rnorm(80)
  a <- numeric(80); b <- numeric(80)
  for (i in 2:80) {
    a[i] <- 0.8 * a[i-1] + e1[i]
    b[i] <- 0.5 * b[i-1] + e2[i]
  }
  fit <- estimate(m, data = data.frame(y1 = a, y2 = b))

  sm <- smooth_states(fit)

  expect_equal(nrow(sm$smoothed_states), 80)
  expect_equal(ncol(sm$smoothed_states), 2)  # 2 states
  expect_equal(ncol(sm$smoothed_obs), 2)     # 2 observables
})

test_that("smooth_shocks works on ML fit", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  e <- rnorm(100)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i - 1] + e[i]
  fit <- estimate(m, data = data.frame(y = z))

  sh <- smooth_shocks(fit)

  expect_s3_class(sh, "dsge_smoothed_shocks")
  expect_equal(nrow(sh$shocks), 99)  # T-1 shocks
  expect_equal(ncol(sh$shocks), 1)   # 1 shock
  expect_equal(length(sh$shock_names), 1)

  # Shock SD should be approximately 1 (or close to estimated sd)
  expect_true(sd(sh$shocks[, 1]) > 0.5 && sd(sh$shocks[, 1]) < 2.0)
})

test_that("shock_decomposition sums to observed data", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  e <- rnorm(100)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i - 1] + e[i]
  fit <- estimate(m, data = data.frame(y = z))

  hd <- shock_decomposition(fit)

  expect_s3_class(hd, "dsge_decomposition")
  expect_equal(dim(hd$decomposition)[1], 100)  # T
  expect_equal(dim(hd$decomposition)[2], 1)    # n_obs
  expect_equal(dim(hd$decomposition)[3], 2)    # 1 shock + initial

  # Key test: sum of contributions = smoothed observables
  recon <- apply(hd$decomposition, c(1, 2), sum)
  # Should match the smoothed observables (from smooth_states)
  sm <- smooth_states(fit)
  # Decomposition sums to smoothed obs
  expect_true(max(abs(recon[, 1] - sm$smoothed_obs[, 1])) < 1e-6)
})

test_that("shock_decomposition with multiple shocks", {
  m <- dsge_model(
    obs(y1 ~ a),
    obs(y2 ~ b),
    state(a ~ rho_a * a),
    state(b ~ rho_b * b),
    start = list(rho_a = 0.5, rho_b = 0.3)
  )

  set.seed(42)
  e1 <- rnorm(80)
  e2 <- rnorm(80)
  a <- numeric(80); b <- numeric(80)
  for (i in 2:80) {
    a[i] <- 0.8 * a[i-1] + e1[i]
    b[i] <- 0.5 * b[i-1] + e2[i]
  }
  fit <- estimate(m, data = data.frame(y1 = a, y2 = b))

  hd <- shock_decomposition(fit)

  expect_equal(dim(hd$decomposition)[3], 3)  # 2 shocks + initial

  # Reconstruction test for both observables
  recon <- apply(hd$decomposition, c(1, 2), sum)
  sm <- smooth_states(fit)
  expect_true(max(abs(recon - sm$smoothed_obs)) < 1e-6)
})

test_that("print methods work", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  e <- rnorm(50)
  z <- numeric(50)
  for (i in 2:50) z[i] <- 0.8 * z[i - 1] + e[i]
  fit <- estimate(m, data = data.frame(y = z))

  expect_output(print(smooth_states(fit)), "DSGE Smoothed States")
  expect_output(print(smooth_shocks(fit)), "DSGE Smoothed Structural Shocks")
  expect_output(print(shock_decomposition(fit)), "DSGE Historical Shock Decomposition")
})

test_that("plot methods execute without error", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  e <- rnorm(50)
  z <- numeric(50)
  for (i in 2:50) z[i] <- 0.8 * z[i - 1] + e[i]
  fit <- estimate(m, data = data.frame(y = z))

  sm <- smooth_states(fit)
  hd <- shock_decomposition(fit)

  # Plot should not error
  pdf(NULL)
  expect_no_error(plot(sm, type = "states"))
  expect_no_error(plot(sm, type = "fit"))
  expect_no_error(plot(hd))
  dev.off()
})

test_that("analytical check: AR(1) smoothed states match true states", {
  # For a simple AR(1) with known rho, the Kalman smoother should
  # recover the true states well
  set.seed(123)
  rho_true <- 0.9
  sd_true <- 1.0
  n <- 200

  z_true <- numeric(n)
  for (i in 2:n) z_true[i] <- rho_true * z_true[i-1] + rnorm(1, sd = sd_true)

  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  fit <- estimate(m, data = data.frame(y = z_true))
  sm <- smooth_states(fit)

  # Smoothed states should correlate strongly with true states
  r <- cor(sm$smoothed_states[, 1], z_true)
  expect_true(r > 0.95)
})
