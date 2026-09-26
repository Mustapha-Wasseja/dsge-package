# The Kalman filter and Lyapunov solver run in C++ (src/kalman.cpp); the
# original R versions (kalman_filter_r() and friends) are kept as a
# reference. These tests check that both give the same results.

random_system <- function(seed, n_s, n_obs, n_shocks, rho = 0.95) {
  set.seed(seed)
  H <- matrix(stats::rnorm(n_s^2), n_s)
  H <- H * rho / max(abs(eigen(H, only.values = TRUE)$values))
  list(H = H,
       M = matrix(stats::rnorm(n_s * n_shocks), n_s),
       G = matrix(stats::rnorm(n_obs * n_s), n_obs),
       D = diag(n_obs),
       y = matrix(stats::rnorm(80 * n_obs), 80))
}

expect_same_filter <- function(a, b, tol = 1e-9) {
  expect_equal(a$loglik, b$loglik, tolerance = tol)
  for (k in c("filtered_states", "predicted_states", "prediction_errors",
              "predicted_obs", "filtered_P", "innovation_var")) {
    expect_equal(a[[k]], b[[k]], tolerance = tol, label = k)
  }
}

test_that("the C++ filter matches the R filter", {
  # (seed, states, observables, shocks), with as many shocks as states so
  # that the covariances stay well conditioned
  for (cfg in list(c(1, 3, 2, 3), c(2, 6, 3, 6), c(3, 12, 4, 12),
                   c(4, 1, 1, 1), c(5, 8, 1, 8))) {
    s <- random_system(cfg[1], cfg[2], cfg[3], cfg[4])
    for (ps in c(0L, 3L)) {
      a <- dsge:::kalman_filter(s$y, s$G, s$H, s$M, s$D, presample = ps)
      b <- dsge:::kalman_filter_r(s$y, s$G, s$H, s$M, s$D, presample = ps)
      expect_true(is.finite(a$loglik))
      expect_same_filter(a, b)
    }
  }
})

test_that("the filters agree to rounding error with fewer shocks than states", {
  # the state covariance is then nearly singular and rounding differences
  # are amplified, so only a looser agreement can be expected
  s <- random_system(1, 6, 2, 3)
  a <- dsge:::kalman_filter(s$y, s$G, s$H, s$M, s$D)
  b <- dsge:::kalman_filter_r(s$y, s$G, s$H, s$M, s$D)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
  expect_equal(a$filtered_states, b$filtered_states, tolerance = 1e-6)
})

test_that("loglik_only returns the same log-likelihood", {
  s <- random_system(11, 5, 2, 3)
  full <- dsge:::kalman_filter(s$y, s$G, s$H, s$M, s$D, presample = 2L)
  quick <- dsge:::kalman_filter(s$y, s$G, s$H, s$M, s$D, presample = 2L,
                                loglik_only = TRUE)
  expect_identical(names(quick), "loglik")
  expect_identical(quick$loglik, full$loglik)
})

test_that("failures give -Inf as in the R filter", {
  # singular innovation variance from the first period
  H <- diag(0.5, 2); M <- diag(2); D <- diag(2)
  G <- rbind(c(1, 0), c(0, 0))
  y <- matrix(stats::rnorm(20), 10)
  a <- dsge:::kalman_filter(y, G, H, M, D)
  expect_equal(a$loglik, -Inf)
  expect_equal(a, dsge:::kalman_filter_r(y, G, H, M, D))
  # missing data
  y[4, 1] <- NA
  expect_equal(dsge:::kalman_filter(y, diag(2), H, M, D)$loglik, -Inf)
})

test_that("the Lyapunov solver matches the Kronecker solution", {
  for (seed in 1:5) {
    set.seed(seed)
    n <- 3L * seed
    H <- matrix(stats::rnorm(n^2), n)
    H <- H * 0.999 / max(abs(eigen(H, only.values = TRUE)$values))
    M <- matrix(stats::rnorm(n * 2), n)
    Q <- M %*% t(M)
    P <- dsge:::compute_unconditional_P(H, Q)
    expect_equal(P, H %*% P %*% t(H) + Q, tolerance = 1e-9)
    expect_equal(P, dsge:::compute_unconditional_P_r(H, Q), tolerance = 1e-9)
    expect_true(isSymmetric(P))
  }
  # a unit root falls back to the same large diagonal as before
  H <- diag(c(1, 0.5))
  expect_identical(dsge:::compute_unconditional_P(H, diag(2)),
                   dsge:::compute_unconditional_P_r(H, diag(2)))
})

test_that("the lik_init = 2 filter matches the R version", {
  m <- read_dynare(text = "
    var y pi r g u;
    varexo eg eu er;
    parameters beta sigma kappa phi_pi rho_r rho_g rho_u;
    beta = 0.99; sigma = 1; kappa = 0.1; phi_pi = 1.5; rho_r = 0.8;
    rho_g = 0.9; rho_u = 0.7;
    model(linear);
      y = y(+1) - 1/sigma * (r - pi(+1)) + g;
      pi = beta * pi(+1) + kappa * y + u;
      r = rho_r * r(-1) + (1 - rho_r) * phi_pi * pi + er;
      g = rho_g * g(-1) + eg;
      u = rho_u * u(-1) + eu;
    end;
    shocks; var eg; stderr 1; var eu; stderr 0.5; var er; stderr 0.25; end;
    shocks; var y; stderr 0.2; end;
    varobs y pi r;
    estimated_params; rho_g, 0.9; end;
    estimation(datafile = 'data.csv', mode_compute = 0, lik_init = 2,
               presample = 2);
  ")
  sol <- solve_dsge(m)
  set.seed(8)
  y <- matrix(stats::rnorm(150), 50)
  a <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D, presample = 2L,
                            init = m$model$kalman_init)
  b <- dsge:::kalman_filter_dynare_state_r(y, sol$G, sol$H, sol$M, sol$D,
                                           2L, m$model$kalman_init)
  expect_true(is.finite(a$loglik))
  expect_equal(a$loglik, b$loglik, tolerance = 1e-10)
  expect_equal(a$prediction_errors, b$prediction_errors, tolerance = 1e-10)
  expect_identical(a$state_names, b$state_names)
})
