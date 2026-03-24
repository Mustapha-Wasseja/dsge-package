# tests/testthat/test-bayes-diagnostics.R

# Helper: create a minimal Bayesian fit for testing
# Uses a pre-built mock object to avoid long MCMC runs
make_mock_bayes <- function() {
  m <- dsge_model(
    obs(y ~ lead(y) + u),
    state(u ~ rho_u * u),
    fixed = list(),
    start = list(rho_u = 0.5)
  )

  # Simulate data
  set.seed(42)
  n <- 50
  u_vals <- numeric(n)
  for (t in 2:n) u_vals[t] <- 0.5 * u_vals[t - 1] + rnorm(1, sd = 0.1)
  dat <- data.frame(y = u_vals)

  # Run a very short MCMC (just enough for diagnostics)
  fit <- bayes_dsge(m, data = dat,
    priors = list(rho_u = prior("beta", shape1 = 5, shape2 = 5)),
    chains = 2L, iter = 200L, warmup = 100L,
    seed = 42L)
  fit
}

test_that("posterior_predictive returns correct structure", {
  skip_on_cran()
  fit <- make_mock_bayes()
  ppc <- posterior_predictive(fit, n_draws = 20, seed = 42)

  expect_s3_class(ppc, "dsge_ppc")
  expect_true("variance" %in% ppc$statistics)
  expect_true("acf1" %in% ppc$statistics)
  expect_equal(length(ppc$variables), 1)

  # P-values should be between 0 and 1
  for (stat in ppc$statistics) {
    expect_true(all(ppc$pvalues[[stat]] >= 0 & ppc$pvalues[[stat]] <= 1,
                    na.rm = TRUE))
  }
})

test_that("posterior_predictive print method works", {
  skip_on_cran()
  fit <- make_mock_bayes()
  ppc <- posterior_predictive(fit, n_draws = 20, seed = 42)
  expect_output(print(ppc), "Posterior Predictive Check")
  expect_output(print(ppc), "variance")
})

test_that("marginal_likelihood returns a number", {
  skip_on_cran()
  fit <- make_mock_bayes()
  ml <- marginal_likelihood(fit)

  expect_s3_class(ml, "dsge_marginal_likelihood")
  expect_true(is.numeric(ml$logml))
  expect_true(is.finite(ml$logml))
})

test_that("marginal_likelihood print method works", {
  skip_on_cran()
  fit <- make_mock_bayes()
  ml <- marginal_likelihood(fit)
  expect_output(print(ml), "Marginal Likelihood")
  expect_output(print(ml), "harmonic_mean")
})

test_that("geweke_test returns correct structure", {
  skip_on_cran()
  fit <- make_mock_bayes()
  gk <- geweke_test(fit)

  expect_s3_class(gk, "dsge_geweke")
  expect_true(is.matrix(gk$z_scores))
  expect_true(is.matrix(gk$p_values))

  # Correct dimensions
  expect_equal(nrow(gk$z_scores), length(fit$param_names))
  expect_equal(ncol(gk$z_scores), fit$n_chains)

  # P-values between 0 and 1
  expect_true(all(gk$p_values >= 0 & gk$p_values <= 1, na.rm = TRUE))
})

test_that("geweke_test print method works", {
  skip_on_cran()
  fit <- make_mock_bayes()
  gk <- geweke_test(fit)
  expect_output(print(gk), "Geweke")
  expect_output(print(gk), "Z-scores")
})

test_that("mcmc_diagnostics returns comprehensive summary", {
  skip_on_cran()
  fit <- make_mock_bayes()
  diag <- mcmc_diagnostics(fit)

  expect_s3_class(diag, "dsge_mcmc_summary")
  expect_true(is.data.frame(diag$table))
  expect_true("ESS" %in% colnames(diag$table))
  expect_true("Rhat" %in% colnames(diag$table))
  expect_true("Geweke_z" %in% colnames(diag$table))
  expect_true("Status" %in% colnames(diag$table))
})

test_that("mcmc_diagnostics print method works", {
  skip_on_cran()
  fit <- make_mock_bayes()
  diag <- mcmc_diagnostics(fit)
  expect_output(print(diag), "MCMC Diagnostic Summary")
  expect_output(print(diag), "Acceptance")
})

test_that("mcmc_diagnostics flags are sensible", {
  skip_on_cran()
  fit <- make_mock_bayes()
  diag <- mcmc_diagnostics(fit)

  # All status values should be character
  expect_true(is.character(diag$table$Status))

  # All should be either "OK" or contain flag keywords
  valid_flags <- c("OK", "low_ESS", "high_Rhat", "Geweke_fail")
  for (s in diag$table$Status) {
    parts <- unlist(strsplit(s, ","))
    expect_true(all(parts %in% valid_flags))
  }
})
