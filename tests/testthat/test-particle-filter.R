# tests/testthat/test-particle-filter.R

# Helper: simple AR(1) model and data
make_ar1_pf <- function(seed = 10L) {
  m <- dsge_model(
    obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.8)
  )
  set.seed(seed)
  n <- 80
  z <- numeric(n)
  for (t in 2:n) z[t] <- 0.8 * z[t - 1] + rnorm(1, sd = 0.2)
  dat <- data.frame(y = z - mean(z))
  sol <- solve_dsge(m, params = c(rho = 0.8), shock_sd = c(z = 0.2))
  list(model = m, data = dat, sol = sol,
       priors = list(rho = prior("beta", shape1 = 4, shape2 = 2)))
}

# --------------------------------------------------------------------------
# particle_filter()
# --------------------------------------------------------------------------

test_that("particle_filter returns finite log-likelihood", {
  d <- make_ar1_pf()
  y <- as.matrix(d$data)
  sol <- d$sol
  Z   <- sol$D %*% sol$G
  pf  <- particle_filter(y, H = sol$H, M = sol$M, Z = Z,
                          n_particles = 200L, meas_sd = 0.01, seed = 1L)

  expect_true(is.numeric(pf$loglik))
  expect_true(is.finite(pf$loglik))
})

test_that("particle_filter returns correct dimensions", {
  d   <- make_ar1_pf()
  y   <- as.matrix(d$data)
  sol <- d$sol
  Z   <- sol$D %*% sol$G
  N   <- nrow(y)
  n_s <- ncol(sol$H)

  pf <- particle_filter(y, H = sol$H, M = sol$M, Z = Z,
                         n_particles = 300L, meas_sd = 0.01, seed = 2L)

  expect_equal(dim(pf$filtered_states), c(N, n_s))
  expect_equal(length(pf$ess), N)
  expect_true(all(pf$ess > 0))
  expect_equal(pf$n_particles, 300L)
})

test_that("particle_filter log-likelihood is close to Kalman (linear model)", {
  d   <- make_ar1_pf()
  y   <- as.matrix(d$data)
  sol <- d$sol
  Z   <- sol$D %*% sol$G

  kf_ll <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D)$loglik
  pf_ll <- particle_filter(y, H = sol$H, M = sol$M, Z = Z,
                            n_particles = 2000L, meas_sd = 0.005,
                            seed = 3L)$loglik

  # With many particles and tiny measurement noise the two should be close
  expect_true(abs(pf_ll - kf_ll) < 5)  # generous tolerance for stochastic test
})

test_that("particle_filter is reproducible with seed", {
  d   <- make_ar1_pf()
  y   <- as.matrix(d$data)
  sol <- d$sol; Z <- sol$D %*% sol$G

  ll1 <- particle_filter(y, sol$H, sol$M, Z, 300L, 0.01, seed = 99L)$loglik
  ll2 <- particle_filter(y, sol$H, sol$M, Z, 300L, 0.01, seed = 99L)$loglik
  expect_equal(ll1, ll2)
})

test_that("particle_filter meas_sd length error", {
  d <- make_ar1_pf()
  y <- as.matrix(d$data)
  sol <- d$sol; Z <- sol$D %*% sol$G
  expect_error(
    particle_filter(y, sol$H, sol$M, Z, 100L, meas_sd = c(0.01, 0.02)),
    "meas_sd"
  )
})

# --------------------------------------------------------------------------
# particle_filter_loglik()
# --------------------------------------------------------------------------

test_that("particle_filter_loglik accepts dsge_solution and returns scalar", {
  d  <- make_ar1_pf()
  ll <- particle_filter_loglik(d$sol, d$data, n_particles = 300L,
                               meas_sd = 0.01, seed = 5L)
  expect_true(is.numeric(ll) && length(ll) == 1L)
  expect_true(is.finite(ll))
})

test_that("particle_filter_loglik errors on non-solution input", {
  expect_error(particle_filter_loglik(list(), data.frame(y = 1)), "dsge_solution")
})

# --------------------------------------------------------------------------
# bayes_particle()
# --------------------------------------------------------------------------

test_that("bayes_particle returns dsge_particle object", {
  skip_on_cran()
  d   <- make_ar1_pf()
  fit <- bayes_particle(d$model, d$data, d$priors,
                        chains = 1L, iter = 150L, warmup = 75L,
                        n_particles = 100L, meas_sd = 0.01,
                        seed = 11L)

  expect_s3_class(fit, "dsge_particle")
  expect_s3_class(fit, "dsge_bayes")
  expect_equal(fit$estimator, "pmmh")
  expect_equal(fit$n_particles, 100L)
})

test_that("bayes_particle posterior has correct dimensions", {
  skip_on_cran()
  d   <- make_ar1_pf()
  fit <- bayes_particle(d$model, d$data, d$priors,
                        chains = 1L, iter = 150L, warmup = 75L,
                        n_particles = 100L, meas_sd = 0.01, seed = 12L)

  # [iteration, parameter, chain]
  expect_equal(dim(fit$posterior)[2], length(fit$param_names))
  expect_equal(dim(fit$posterior)[3], 1L)
})

test_that("bayes_particle print method works", {
  skip_on_cran()
  d   <- make_ar1_pf()
  fit <- bayes_particle(d$model, d$data, d$priors,
                        chains = 1L, iter = 100L, warmup = 50L,
                        n_particles = 80L, seed = 13L)
  expect_output(print(fit), "Particle MCMC")
  expect_output(print(fit), "Particles")
})

test_that("marginal_likelihood works on dsge_particle output", {
  skip_on_cran()
  d   <- make_ar1_pf()
  fit <- bayes_particle(d$model, d$data, d$priors,
                        chains = 1L, iter = 150L, warmup = 75L,
                        n_particles = 100L, meas_sd = 0.01, seed = 14L)
  ml  <- marginal_likelihood(fit)
  expect_s3_class(ml, "dsge_marginal_likelihood")
  expect_true(is.finite(ml$logml))
})
