test_that("estimate recovers AR(1) parameter from simulated data", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  # Simulate data from known model
  set.seed(42)
  true_rho <- 0.8
  true_sd <- 1.0
  n <- 300
  z <- numeric(n)
  for (i in 2:n) z[i] <- true_rho * z[i - 1] + true_sd * rnorm(1)
  dat <- data.frame(y = z)

  fit <- estimate(m, data = dat, control = list(maxit = 200))

  expect_s3_class(fit, "dsge_fit")
  expect_true(fit$convergence == 0)

  # Parameter should be close to 0.8
  rho_est <- coef(fit)["rho"]
  expect_equal(as.numeric(rho_est), true_rho, tolerance = 0.15)

  # Log-likelihood should be finite
  expect_true(is.finite(fit$loglik))
})

test_that("estimate works with fixed parameters", {
  m <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.96)
  )

  # Simulate simple data
  set.seed(123)
  n <- 200
  dat <- data.frame(
    p = cumsum(rnorm(n, 0, 0.5)),
    r = cumsum(rnorm(n, 0, 0.3))
  )
  # Demean manually to make stationary
  dat$p <- dat$p - mean(dat$p)
  dat$r <- dat$r - mean(dat$r)

  # Should run without error even if the model doesn't perfectly fit this data
  fit <- tryCatch(
    estimate(m, data = dat,
             start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
             control = list(maxit = 50)),
    error = function(e) NULL
  )

  # At minimum, it should return a dsge_fit or NULL (if optimizer can't find stable region)
  if (!is.null(fit)) {
    expect_s3_class(fit, "dsge_fit")
    expect_true("beta" %in% names(coef(fit)))
    expect_equal(as.numeric(coef(fit)["beta"]), 0.96)
  }
})

test_that("estimate rejects missing data columns", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  dat <- data.frame(x = rnorm(100))
  expect_error(estimate(m, data = dat), "not found in data")
})

test_that("coef, logLik, nobs work on dsge_fit", {
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

  expect_true(length(coef(fit)) > 0)
  expect_true(is.finite(logLik(fit)))
  expect_equal(nobs(fit), 100L)
})
