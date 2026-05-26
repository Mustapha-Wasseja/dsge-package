# Tests for bayes_smc()

make_ar1_dat <- function(TT = 80L, seed = 1L) {
  set.seed(seed)
  e <- rnorm(TT)
  z <- numeric(TT)
  for (i in 2:TT) z[i] <- 0.6 * z[i - 1] + e[i]
  list(model = dsge_model(
         obs(y ~ rho * y_lag + e),
         state(y_lag ~ y, shock = FALSE),
         state(e ~ phi * e),
         start = list(rho = 0.5, phi = 0.5)),
       data = data.frame(y = z))
}


test_that("bayes_smc returns expected structure", {
  d <- make_ar1_dat()
  fit <- bayes_smc(d$model, d$data,
                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2),
                                 phi = prior("beta", shape1 = 2, shape2 = 2)),
                   n_particles = 80L, n_phi = 8L,
                   seed = 1L)
  expect_s3_class(fit, "dsge_smc")
  expect_s3_class(fit, "dsge_bayes")
  expect_equal(dim(fit$posterior), c(80, 3, 1))
  expect_true(is.finite(fit$log_marg_lik))
  expect_equal(length(fit$ess_path), 8L)
})


test_that("posterior particles approach the data-generating process", {
  d <- make_ar1_dat(TT = 200L, seed = 2L)
  fit <- bayes_smc(d$model, d$data,
                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2),
                                 phi = prior("beta", shape1 = 2, shape2 = 2)),
                   n_particles = 150L, n_phi = 15L,
                   seed = 2L)
  # Posterior mean of rho should be close to 0.6 (true DGP)
  rho_mean <- mean(fit$posterior[, "rho", 1])
  expect_lt(abs(rho_mean - 0.6), 0.4)
})


test_that("ESS path is finite and bounded", {
  d <- make_ar1_dat()
  fit <- bayes_smc(d$model, d$data,
                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2),
                                 phi = prior("beta", shape1 = 2, shape2 = 2)),
                   n_particles = 50L, n_phi = 6L,
                   seed = 7L)
  expect_true(all(is.finite(fit$ess_path)))
  expect_true(all(fit$ess_path >= 1))
  expect_true(all(fit$ess_path <= 50))
})


test_that("acceptance path within [0, 1]", {
  d <- make_ar1_dat()
  fit <- bayes_smc(d$model, d$data,
                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2),
                                 phi = prior("beta", shape1 = 2, shape2 = 2)),
                   n_particles = 50L, n_phi = 5L,
                   seed = 4L)
  expect_true(all(fit$acceptance_path >= 0 & fit$acceptance_path <= 1))
})


test_that("print method works", {
  d <- make_ar1_dat()
  fit <- bayes_smc(d$model, d$data,
                   priors = list(rho = prior("beta", shape1 = 2, shape2 = 2),
                                 phi = prior("beta", shape1 = 2, shape2 = 2)),
                   n_particles = 40L, n_phi = 4L,
                   seed = 8L)
  expect_output(print(fit), "Tempered Sequential Monte Carlo")
})
