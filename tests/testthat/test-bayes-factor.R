# tests/testthat/test-bayes-factor.R

# Helper: two tiny Bayesian fits for comparison tests
make_two_fits <- function() {
  m1 <- dsge_model(
    obs(y ~ lead(y) + u),
    state(u ~ rho_u * u),
    start = list(rho_u = 0.5)
  )
  m2 <- dsge_model(
    obs(y ~ u),
    state(u ~ rho_u * u),
    start = list(rho_u = 0.5)
  )
  set.seed(7)
  n <- 50
  u <- numeric(n)
  for (t in 2:n) u[t] <- 0.5 * u[t - 1] + rnorm(1, sd = 0.1)
  dat <- data.frame(y = u)
  pr  <- list(rho_u = prior("beta", shape1 = 5, shape2 = 5))
  fit1 <- bayes_dsge(m1, dat, pr, chains = 1L, iter = 200L, warmup = 100L, seed = 7L)
  fit2 <- bayes_dsge(m2, dat, pr, chains = 1L, iter = 200L, warmup = 100L, seed = 7L)
  list(fit1 = fit1, fit2 = fit2)
}

test_that("bayes_factor errors on fewer than two models", {
  skip_on_cran()
  fits <- make_two_fits()
  expect_error(bayes_factor(fits$fit1), "at least two")
})

test_that("bayes_factor returns correct structure", {
  skip_on_cran()
  fits <- make_two_fits()
  bf <- bayes_factor(M1 = fits$fit1, M2 = fits$fit2)

  expect_s3_class(bf, "dsge_bayes_factor")
  expect_equal(bf$model_names, c("M1", "M2"))
  expect_length(bf$log_ml, 2L)
  expect_true(all(is.finite(bf$log_ml)))
  expect_equal(dim(bf$bf_matrix), c(2L, 2L))
  # Diagonal must be zero
  expect_equal(diag(bf$bf_matrix), c(0, 0), ignore_attr = TRUE)
  # Anti-symmetric
  expect_equal(bf$bf_matrix[1, 2], -bf$bf_matrix[2, 1])
})

test_that("bayes_factor log_ml matches standalone marginal_likelihood()", {
  skip_on_cran()
  fits <- make_two_fits()
  ml1 <- marginal_likelihood(fits$fit1)
  ml2 <- marginal_likelihood(fits$fit2)
  bf  <- bayes_factor(fits$fit1, fits$fit2)

  expect_equal(bf$log_ml[1], ml1$logml, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(bf$log_ml[2], ml2$logml, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("bayes_factor accepts dsge_marginal_likelihood objects", {
  skip_on_cran()
  fits <- make_two_fits()
  ml1 <- marginal_likelihood(fits$fit1)
  ml2 <- marginal_likelihood(fits$fit2)
  bf  <- bayes_factor(ml1, ml2)

  expect_s3_class(bf, "dsge_bayes_factor")
  expect_equal(bf$log_ml[1], ml1$logml, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("bayes_factor posterior probs sum to 1", {
  skip_on_cran()
  fits <- make_two_fits()
  bf <- bayes_factor(fits$fit1, fits$fit2)
  expect_equal(sum(bf$posterior_probs), 1, tolerance = 1e-10)
  expect_true(all(bf$posterior_probs >= 0))
})

test_that("bayes_factor respects custom prior_odds", {
  skip_on_cran()
  fits <- make_two_fits()
  bf_eq  <- bayes_factor(fits$fit1, fits$fit2)
  bf_pri <- bayes_factor(fits$fit1, fits$fit2, prior_odds = c(0.8, 0.2))

  # Posterior probs differ when priors differ
  expect_false(isTRUE(all.equal(bf_eq$posterior_probs,
                                bf_pri$posterior_probs)))
  expect_equal(sum(bf_pri$posterior_probs), 1, tolerance = 1e-10)
})

test_that("bayes_factor prior_odds length mismatch errors", {
  skip_on_cran()
  fits <- make_two_fits()
  expect_error(bayes_factor(fits$fit1, fits$fit2, prior_odds = c(1, 2, 3)),
               "same length")
})

test_that("bayes_factor handles three models", {
  skip_on_cran()
  fits <- make_two_fits()
  bf <- bayes_factor(fits$fit1, fits$fit2, fits$fit1)  # M1 twice is fine
  expect_equal(dim(bf$bf_matrix), c(3L, 3L))
  expect_equal(sum(bf$posterior_probs), 1, tolerance = 1e-10)
})

test_that("bayes_factor print method works", {
  skip_on_cran()
  fits <- make_two_fits()
  bf <- bayes_factor(A = fits$fit1, B = fits$fit2)
  expect_output(print(bf), "Bayesian Model Comparison")
  expect_output(print(bf), "Log marginal likelihoods")
  expect_output(print(bf), "Posterior model probabilities")
  expect_output(print(bf), "BF")
})

test_that(".kass_raftery returns correct labels", {
  expect_match(dsge:::.kass_raftery(-5),   "Negative")
  expect_match(dsge:::.kass_raftery(1),    "bare mention")
  expect_match(dsge:::.kass_raftery(4),    "Positive")
  expect_match(dsge:::.kass_raftery(8),    "Strong")
  expect_match(dsge:::.kass_raftery(15),   "Very strong")
})
