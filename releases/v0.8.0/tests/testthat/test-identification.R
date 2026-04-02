# Tests for check_identification(), parameter_sensitivity(), prior_posterior_update()

# ---- Helper: create and estimate a simple AR(1) model ----
make_ar1_fit <- function(seed = 42, n = 100) {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )
  set.seed(seed)
  e <- rnorm(n)
  z <- numeric(n)
  for (i in 2:n) z[i] <- 0.8 * z[i - 1] + e[i]
  estimate(m, data = data.frame(y = z))
}

make_2state_fit <- function(seed = 42, n = 100) {
  m <- dsge_model(
    obs(y1 ~ a),
    obs(y2 ~ b),
    state(a ~ rho_a * a),
    state(b ~ rho_b * b),
    start = list(rho_a = 0.5, rho_b = 0.3)
  )
  set.seed(seed)
  e1 <- rnorm(n); e2 <- rnorm(n)
  a <- numeric(n); b <- numeric(n)
  for (i in 2:n) {
    a[i] <- 0.8 * a[i-1] + e1[i]
    b[i] <- 0.5 * b[i-1] + e2[i]
  }
  estimate(m, data = data.frame(y1 = a, y2 = b))
}


# ===========================================================================
# check_identification() tests
# ===========================================================================

test_that("check_identification works on AR(1) fit", {
  fit <- make_ar1_fit()
  id <- check_identification(fit)

  expect_s3_class(id, "dsge_identification")
  expect_true(id$identified)
  expect_equal(id$n_params, 1)
  expect_equal(id$rank, 1)
  expect_true(id$condition_number > 0)
  expect_true(all(id$singular_values > 0))
})

test_that("check_identification detects well-identified 2-param model", {
  fit <- make_2state_fit()
  id <- check_identification(fit)

  expect_s3_class(id, "dsge_identification")
  expect_true(id$identified)
  expect_equal(id$n_params, 2)
  expect_equal(id$rank, 2)
  expect_equal(nrow(id$summary), 2)
  # Both parameters should be at least moderately identified
  expect_true(all(id$summary$status %in% c("strong", "moderate")))
})

test_that("check_identification summary has expected columns", {
  fit <- make_ar1_fit()
  id <- check_identification(fit)

  expect_true(all(c("parameter", "strength", "rel_strength", "status")
                  %in% colnames(id$summary)))
})

test_that("check_identification print works", {
  fit <- make_ar1_fit()
  id <- check_identification(fit)
  expect_output(print(id), "DSGE Local Identification Diagnostics")
  expect_output(print(id), "ALL parameters locally identified")
})

test_that("check_identification plot runs without error", {
  fit <- make_ar1_fit()
  id <- check_identification(fit)
  pdf(NULL)
  expect_no_error(plot(id))
  dev.off()
})

test_that("check_identification n_lags parameter works", {
  fit <- make_ar1_fit()
  id2 <- check_identification(fit, n_lags = 2)
  id8 <- check_identification(fit, n_lags = 8)

  # More lags = more moments
  expect_true(id8$n_moments > id2$n_moments)
  # Both should identify the single parameter
  expect_true(id2$identified)
  expect_true(id8$identified)
})


# ===========================================================================
# parameter_sensitivity() tests
# ===========================================================================

test_that("parameter_sensitivity works on AR(1) fit (loglik)", {
  fit <- make_ar1_fit()
  sa <- parameter_sensitivity(fit, what = "loglik")

  expect_s3_class(sa, "dsge_sensitivity")
  expect_true(!is.null(sa$loglik))
  expect_equal(nrow(sa$loglik), 1)  # 1 parameter
  expect_true(is.finite(sa$loglik$derivative[1]))
  expect_true(is.finite(sa$loglik$elasticity[1]))
})

test_that("parameter_sensitivity works on AR(1) fit (irf)", {
  fit <- make_ar1_fit()
  sa <- parameter_sensitivity(fit, what = "irf", irf_horizon = 10)

  expect_s3_class(sa, "dsge_sensitivity")
  expect_true(!is.null(sa$irf))
  expect_equal(nrow(sa$irf), 1)
  expect_true(sa$irf$max_sensitivity[1] > 0)
})

test_that("parameter_sensitivity with multiple params", {
  fit <- make_2state_fit()
  sa <- parameter_sensitivity(fit, what = c("loglik", "irf"))

  expect_equal(nrow(sa$loglik), 2)  # 2 parameters
  expect_equal(nrow(sa$irf), 2)
  expect_true(all(is.finite(sa$loglik$derivative)))
})

test_that("parameter_sensitivity print works", {
  fit <- make_ar1_fit()
  sa <- parameter_sensitivity(fit, what = c("loglik", "irf"))
  expect_output(print(sa), "DSGE Parameter Sensitivity Analysis")
})

test_that("parameter_sensitivity plot runs without error", {
  fit <- make_ar1_fit()
  sa <- parameter_sensitivity(fit, what = c("loglik", "irf"))
  pdf(NULL)
  expect_no_error(plot(sa))
  dev.off()
})

test_that("parameter_sensitivity policy works", {
  # For simple AR models, G may not depend on rho (identity mapping).
  # Just verify the function runs and returns the correct structure.
  fit <- make_2state_fit()
  sa <- parameter_sensitivity(fit, what = "policy")

  expect_true(!is.null(sa$policy))
  expect_equal(nrow(sa$policy), 2)
  expect_true(all(c("parameter", "max_G_sensitivity") %in% colnames(sa$policy)))
})


# ===========================================================================
# prior_posterior_update() tests
# ===========================================================================

test_that("prior_posterior_update works on Bayesian fit", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )
  set.seed(42)
  e <- rnorm(200)
  z <- numeric(200)
  for (i in 2:200) z[i] <- 0.8 * z[i - 1] + e[i]

  fit <- bayes_dsge(m, data = data.frame(y = z),
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 1500, seed = 1)

  pp <- prior_posterior_update(fit)

  expect_s3_class(pp, "dsge_prior_posterior")
  expect_true(nrow(pp$summary) >= 1)

  # rho should be strongly updated (prior mean = 0.5, true = 0.8)
  rho_row <- pp$summary[pp$summary$parameter == "rho", ]
  expect_true(rho_row$sd_ratio < 1)  # posterior tighter than prior
  expect_true(rho_row$mean_shift > 0)  # posterior moved from prior mean
})

test_that("prior_posterior_update flags weak updates", {
  # With a very tight prior centered on truth, update should be weak
  # (posterior can't move much from the prior)
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )
  set.seed(42)
  e <- rnorm(50)
  z <- numeric(50)
  for (i in 2:50) z[i] <- 0.8 * z[i - 1] + e[i]

  # Very tight prior near truth
  fit <- bayes_dsge(m, data = data.frame(y = z),
                    priors = list(rho = prior("beta", shape1 = 80, shape2 = 20)),
                    chains = 2, iter = 1500, seed = 1)

  pp <- prior_posterior_update(fit)
  rho_row <- pp$summary[pp$summary$parameter == "rho", ]

  # SD ratio should be close to 1 (prior dominates)
  expect_true(rho_row$sd_ratio > 0.5)
})

test_that("prior_posterior_update print works", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )
  set.seed(42)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i-1] + rnorm(1)

  fit <- bayes_dsge(m, data = data.frame(y = z),
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 1500, seed = 1)

  pp <- prior_posterior_update(fit)
  expect_output(print(pp), "Prior vs Posterior Update Diagnostics")
})

test_that("prior_posterior_update plot runs without error", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )
  set.seed(42)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i-1] + rnorm(1)

  fit <- bayes_dsge(m, data = data.frame(y = z),
                    priors = list(rho = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 1500, seed = 1)

  pp <- prior_posterior_update(fit)
  pdf(NULL)
  expect_no_error(plot(pp))
  dev.off()
})
