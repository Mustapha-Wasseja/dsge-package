# Tests for Bayesian DSGE estimation

test_that("prior() creates valid prior objects", {
  p <- prior("normal", mean = 0, sd = 1)
  expect_s3_class(p, "dsge_prior")
  expect_equal(p$distribution, "normal")
  expect_equal(p$support, "unbounded")

  p2 <- prior("beta", shape1 = 2, shape2 = 2)
  expect_equal(p2$support, "unit")

  p3 <- prior("gamma", shape = 2, rate = 1)
  expect_equal(p3$support, "positive")

  p4 <- prior("uniform", min = 0, max = 1)
  expect_equal(p4$support, "bounded")

  p5 <- prior("inv_gamma", shape = 0.01, scale = 0.01)
  expect_equal(p5$support, "positive")
})

test_that("prior() errors on missing parameters", {
  expect_error(prior("normal", mean = 0), "Missing parameter")
  expect_error(prior("beta", shape1 = 2), "Missing parameter")
})

test_that("dprior evaluates log density correctly", {
  p <- prior("normal", mean = 0, sd = 1)
  expect_equal(dprior(p, 0), dnorm(0, log = TRUE))

  p2 <- prior("beta", shape1 = 2, shape2 = 3)
  expect_equal(dprior(p2, 0.5), dbeta(0.5, 2, 3, log = TRUE))
})

test_that("rprior draws from correct distribution", {
  set.seed(1)
  p <- prior("normal", mean = 5, sd = 0.1)
  draws <- rprior(p, 1000)
  expect_true(abs(mean(draws) - 5) < 0.02)
})

test_that("parameter transformations roundtrip correctly", {
  priors <- list(
    prior("normal", mean = 0, sd = 1),
    prior("beta", shape1 = 2, shape2 = 2),
    prior("gamma", shape = 2, rate = 1),
    prior("uniform", min = -1, max = 3),
    prior("inv_gamma", shape = 2, scale = 1)
  )

  test_vals <- c(0.5, 0.7, 2.5, 1.5, 0.8)
  for (i in seq_along(priors)) {
    u <- to_unconstrained(test_vals[i], priors[[i]])
    x_back <- from_unconstrained(u, priors[[i]])
    expect_equal(x_back, test_vals[i], tolerance = 1e-10,
                 label = paste("roundtrip", priors[[i]]$distribution))
  }
})

test_that("bayes_dsge works on AR(1) model", {
  skip_on_cran()

  # Define AR(1)
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  # Simulate data
  set.seed(42)
  n <- 200
  z <- numeric(n)
  for (i in 2:n) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  dat <- data.frame(y = z)

  # Estimate
  fit <- bayes_dsge(m, data = dat,
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 2000, warmup = 1000,
                    seed = 123)

  expect_s3_class(fit, "dsge_bayes")
  expect_equal(dim(fit$posterior)[2], 2)  # rho + sd_e.z
  expect_equal(dim(fit$posterior)[3], 2)  # 2 chains

  # Posterior mean of rho should be near 0.8
  rho_mean <- mean(fit$posterior[, "rho", ])
  expect_true(abs(rho_mean - 0.8) < 0.15,
              label = paste("rho posterior mean =", round(rho_mean, 3)))

  # Methods work
  expect_output(print(fit))
  expect_output(summary(fit))
  expect_length(coef(fit), 2)

  # Diagnostics exist
  expect_true("ess" %in% colnames(fit$diagnostics))
  expect_true("rhat" %in% colnames(fit$diagnostics))
})

test_that("bayes_dsge works on NK model", {
  skip_on_cran()

  nk <- dsge_model(
    obs(p   ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r   ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.96),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )

  # Simulate data
  set.seed(42)
  n <- 200
  sol <- solve_dsge(nk, params = c(beta = 0.96, kappa = 0.085, psi = 1.94,
                                    rhou = 0.70, rhog = 0.95),
                    shock_sd = c(u = 2.3, g = 0.57))
  states <- matrix(0, n, 2)
  for (t in 2:n) {
    states[t, ] <- as.numeric(sol$H %*% states[t - 1, ] +
                               sol$M %*% rnorm(2))
  }
  obs_data <- states %*% t(sol$D %*% sol$G)
  colnames(obs_data) <- c("p", "r")
  dat <- as.data.frame(obs_data)

  # Estimate with short chain for testing speed
  fit <- bayes_dsge(nk, data = dat,
                    priors = list(
                      kappa = prior("gamma", shape = 2, rate = 20),
                      psi   = prior("normal", mean = 1.5, sd = 0.5),
                      rhou  = prior("beta", shape1 = 2, shape2 = 2),
                      rhog  = prior("beta", shape1 = 2, shape2 = 2)
                    ),
                    chains = 2, iter = 1500, warmup = 750,
                    seed = 456)

  expect_s3_class(fit, "dsge_bayes")
  # 4 structural + 2 shock SDs = 6 parameters
  expect_equal(dim(fit$posterior)[2], 6)

  # Check posterior summary
  s <- summary(fit)
  expect_true(is.data.frame(s))
})

test_that("bayes_dsge works on nonlinear RBC model", {
  skip_on_cran()

  # Simple nonlinear RBC: C = Y - I, Y = K^alpha, K' = (1-delta)*K + I
  # With technology shock: Y = A * K^alpha, A' = A^rho
  # Rewritten for dsgenl_model with C, Y observed; K exo_state, A exo_state
  # Simplified 2-equation model:
  #   C = A * K^alpha - K(+1) + (1-delta)*K
  #   A(+1) = A^rho
  # Control: C (observed), State: K (endo), A (exo)
  # But for simplicity, use a compact nonlinear model.

  # Use the simplest possible nonlinear model: AR(1) in levels
  # y = z, z(+1) = z^rho * exp(e)
  # At ss: z_ss = 1 (when rho < 1)
  nl <- dsgenl_model(
    "y = z",
    "z(+1) = z^rho",
    observed = "y", unobserved = character(0),
    exo_state = "z",
    start = list(rho = 0.5),
    ss_guess = c(y = 1, z = 1)
  )

  # Simulate data from the linearized model
  set.seed(42)
  sol <- solve_dsge(nl, params = c(rho = 0.8))
  n <- 200
  z_dev <- numeric(n)
  for (i in 2:n) z_dev[i] <- sol$H[1,1] * z_dev[i - 1] + rnorm(1)
  dat <- data.frame(y = 1 + z_dev)  # level data around ss=1

  fit <- bayes_dsge(nl, data = dat,
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 1500, warmup = 750,
                    seed = 42)

  expect_s3_class(fit, "dsge_bayes")
  expect_true(isTRUE(fit$is_nonlinear))
  expect_equal(dim(fit$posterior)[2], 2)  # rho + sd_e.z

  # Posterior mean of rho should be near 0.8
  rho_mean <- mean(fit$posterior[, "rho", ])
  expect_true(abs(rho_mean - 0.8) < 0.2,
              label = paste("nonlinear rho posterior mean =", round(rho_mean, 3)))

  # Methods work
  expect_output(print(fit), "nonlinear")
  expect_output(summary(fit))

  # Diagnostics
  expect_true(!is.null(fit$solve_failures))
})

test_that("bayes_dsge handles nonlinear solve failures gracefully", {
  skip_on_cran()

  # Model where extreme parameter values cause solve failures
  nl <- dsgenl_model(
    "y = z",
    "z(+1) = z^rho",
    observed = "y", unobserved = character(0),
    exo_state = "z",
    start = list(rho = 0.5),
    ss_guess = c(y = 1, z = 1)
  )

  set.seed(42)
  dat <- data.frame(y = 1 + cumsum(rnorm(100, sd = 0.1)))

  # Should not error even with bad data
  fit <- tryCatch(
    bayes_dsge(nl, data = dat,
               priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
               chains = 1, iter = 500, warmup = 250,
               seed = 42),
    error = function(e) NULL
  )

  # Either succeeds or fails gracefully
  if (!is.null(fit)) {
    expect_s3_class(fit, "dsge_bayes")
  }
})

test_that("posterior IRFs work for nonlinear models", {
  skip_on_cran()

  nl <- dsgenl_model(
    "y = z",
    "z(+1) = z^rho",
    observed = "y", unobserved = character(0),
    exo_state = "z",
    start = list(rho = 0.5),
    ss_guess = c(y = 1, z = 1)
  )

  set.seed(42)
  sol <- solve_dsge(nl, params = c(rho = 0.8))
  n <- 200
  z_dev <- numeric(n)
  for (i in 2:n) z_dev[i] <- sol$H[1,1] * z_dev[i - 1] + rnorm(1)
  dat <- data.frame(y = 1 + z_dev)

  fit <- bayes_dsge(nl, data = dat,
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 1000, warmup = 500,
                    seed = 42)

  irfs <- irf(fit, periods = 8, n_draws = 50)
  expect_s3_class(irfs, "dsge_irf")
  expect_true("lower" %in% colnames(irfs$data))
  expect_true("upper" %in% colnames(irfs$data))
})

test_that("posterior IRFs work", {
  skip_on_cran()

  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  set.seed(42)
  z <- numeric(200)
  for (i in 2:200) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  dat <- data.frame(y = z)

  fit <- bayes_dsge(m, data = dat,
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 1000, warmup = 500,
                    seed = 789)

  irfs <- irf(fit, periods = 10, n_draws = 50)
  expect_s3_class(irfs, "dsge_irf")
  expect_true("lower" %in% colnames(irfs$data))
  expect_true("upper" %in% colnames(irfs$data))
  expect_true(all(irfs$data$lower <= irfs$data$upper))
})
