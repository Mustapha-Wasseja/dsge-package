test_that("policy_matrix returns correct values for AR(1)", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )

  sol <- solve_dsge(m, params = c(rho = 0.8))
  G <- policy_matrix(sol, se = FALSE)

  expect_equal(as.numeric(G), 1, tolerance = 1e-10)
})

test_that("transition_matrix returns correct values for AR(1)", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )

  sol <- solve_dsge(m, params = c(rho = 0.8))
  H <- transition_matrix(sol, se = FALSE)

  expect_equal(as.numeric(H), 0.8, tolerance = 1e-10)
})

test_that("stability returns correct classification", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )

  # Stable model
  sol <- solve_dsge(m, params = c(rho = 0.8))
  stab <- stability(sol)

  expect_true(stab$stable)
  expect_equal(stab$n_stable, 1L)
  expect_equal(stab$classification, "stable")
})

test_that("irf computes correct impulse responses for AR(1)", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )

  sol <- solve_dsge(m, params = c(rho = 0.8))
  irfs <- irf(sol, periods = 5, se = FALSE)

  expect_s3_class(irfs, "dsge_irf")

  # For AR(1) with G=1, H=0.8, M=1:
  # Response of y to shock e.z:
  #   period 0: 1 * 0.8^0 * 1 = 1
  #   period 1: 1 * 0.8^1 * 1 = 0.8
  #   period 2: 1 * 0.8^2 * 1 = 0.64
  y_resp <- irfs$data[irfs$data$response == "y", ]
  expected <- 0.8^(0:5)
  expect_equal(y_resp$value, expected, tolerance = 1e-10)
})

test_that("forecast produces reasonable output", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  dat <- data.frame(y = z)

  fit <- estimate(m, data = dat)
  fc <- forecast(fit, horizon = 10)

  expect_s3_class(fc, "dsge_forecast")
  expect_equal(nrow(fc$forecasts), 10)
  expect_true(all(is.finite(fc$forecasts$value)))

  # Forecast should converge toward the mean
  vals <- fc$forecasts$value
  expect_true(abs(vals[10]) < abs(vals[1]) + abs(mean(dat$y)) + 5)
})

test_that("predict returns correct dimensions", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  dat <- data.frame(y = z)

  fit <- estimate(m, data = dat)

  pred_obs <- predict(fit, type = "observed")
  expect_equal(nrow(pred_obs), 100)
  expect_equal(ncol(pred_obs), 1)

  pred_state <- predict(fit, type = "state")
  expect_equal(nrow(pred_state), 100)
  expect_equal(ncol(pred_state), 1)
})
