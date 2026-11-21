# The exact diffuse Kalman filter (Dynare's diffuse_filter, lik_init = 3)
# and the beta prior on an interval other than (0, 1). The reference
# log-likelihoods were computed by Dynare 7.0 (MATLAB R2023b) on the same
# models and data, with its steady-state Kalman shortcut switched off
# (options_.riccati_tol = 0), which otherwise changes them by up to 1e-5.

diffuse_loglik <- function(m, y, theta) {
  sd <- m$shock_sd
  for (k in names(theta)[startsWith(names(theta), "sd_e.")]) {
    sd[sub("sd_e.", "", k, fixed = TRUE)] <- theta[[k]]
  }
  pars <- c(theta[!startsWith(names(theta), "sd_e.")], unlist(m$model$fixed))
  sol <- solve_dsge(m$model, params = pars, shock_sd = sd[names(m$shock_sd)])
  kf <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D, init = m$model$kalman_init,
                             loglik_only = TRUE)
  kf$loglik
}

dk_data <- function() {
  set.seed(7)
  n <- 60
  x <- cumsum(rnorm(n, sd = 0.5))
  z <- as.numeric(stats::filter(rnorm(n, sd = 0.3), 0.6, "recursive"))
  slope <- cumsum(rnorm(n, sd = 0.05))
  level <- cumsum(slope + rnorm(n, sd = 0.2))
  list(a = cbind(yobs = x + z, zobs = z), b = cbind(yobs = level + rnorm(n, sd = 0.4)))
}

test_that("read_dynare() reads diffuse_filter as lik_init = 3", {
  m <- read_dynare(text = "var x y; varexo e; parameters r; r = 0.5;
    model(linear); x = x(-1) + e; y = x; end;
    shocks; var e; stderr 1; end;
    estimated_params; stderr e, inv_gamma_pdf, 1, 2; end;
    varobs y;
    estimation(diffuse_filter, mode_compute = 0);")
  expect_identical(m$estimation$lik_init, 3L)
  expect_identical(m$model$kalman_init$type, "diffuse")
  expect_false(any(grepl("lik_init", m$notes)))
})

test_that("the diffuse filter reproduces Dynare: observed random walk plus AR(1)", {
  # one unit root, seen by one of two observed variables: Finf is singular,
  # so the observations are processed one at a time (Dynare's univariate
  # diffuse filter); one diffuse period
  m <- read_dynare(text = "var x z yobs zobs; varexo ex ez; parameters rho; rho = 0.6;
    model(linear); x = x(-1) + ex; z = rho*z(-1) + ez; yobs = x + z; zobs = z; end;
    shocks; var ex; stderr 0.5; var ez; stderr 0.3; end;
    estimated_params; rho, 0.6, beta_pdf, 0.5, 0.2;
      stderr ex, 0.5, inv_gamma_pdf, 0.5, 2; stderr ez, 0.3, inv_gamma_pdf, 0.3, 2; end;
    varobs yobs zobs;
    estimation(diffuse_filter, mode_compute = 0);")
  y <- dk_data()$a
  pts <- rbind(c(0.5, 0.3, 0.6), c(0.4, 0.35, 0.8), c(0.7, 0.2, 0.3))
  dynare <- c(-49.008825848585, -54.156603215530, -77.973862921475)
  for (j in 1:3) {
    th <- stats::setNames(pts[j, ], c("sd_e.ex", "sd_e.ez", "rho"))
    expect_equal(diffuse_loglik(m, y, th), dynare[j], tolerance = 1e-9)
  }
})

test_that("the diffuse filter reproduces Dynare: local linear trend with measurement error", {
  # a double unit root (level and slope): two diffuse periods, and the slope
  # reaches the observed variable only one period later
  m <- read_dynare(text = "var mu beta yobs; varexo eta zeta; parameters c; c = 1;
    model(linear); mu = mu(-1) + c*beta(-1) + eta; beta = beta(-1) + zeta; yobs = mu; end;
    shocks; var eta; stderr 0.2; var zeta; stderr 0.05; end;
    estimated_params; stderr eta, 0.2, inv_gamma_pdf, 0.2, 2;
      stderr zeta, 0.05, inv_gamma_pdf, 0.05, 2; stderr yobs, 0.4, inv_gamma_pdf, 0.4, 2; end;
    varobs yobs;
    estimation(diffuse_filter, mode_compute = 0);")
  y <- dk_data()$b
  pts <- rbind(c(0.2, 0.05, 0.4), c(0.3, 0.02, 0.3), c(0.1, 0.1, 0.6))
  dynare <- c(-55.001445505149, -54.953286115717, -60.322729542138)
  for (j in 1:3) {
    th <- stats::setNames(pts[j, ], c("sd_e.eta", "sd_e.zeta", "sd_e.yobs_me"))
    expect_equal(diffuse_loglik(m, y, th), dynare[j], tolerance = 1e-9)
  }
  sol <- solve_dsge(m$model, shock_sd = c(eta = 0.2, zeta = 0.05, yobs_me = 0.4))
  sys <- dsge:::dynare_state_system(sol$G, sol$H, sol$M, sol$D, m$model$kalman_init)
  kf <- dsge:::kalman_filter_diffuse(y, sys, 0L)
  expect_identical(kf$diffuse_periods, 2L)
})

test_that("without unit roots the diffuse filter is the stationary filter", {
  txt <- "var x z y1 y2; varexo e1 e2; parameters a b; a = 0.8; b = 0.4;
    model(linear); x = a*x(-1) + e1; z = b*z(-1) + 0.5*x(-1) + e2; y1 = x + z; y2 = z - x(-1); end;
    shocks; var e1; stderr 0.3; var e2; stderr 0.2; end;
    varobs y1 y2;"
  m3 <- read_dynare(text = paste(txt, "estimation(diffuse_filter);"))
  m1 <- read_dynare(text = paste(txt, "estimation(mode_compute = 0);"))
  set.seed(3)
  y <- matrix(rnorm(80), 40, 2, dimnames = list(NULL, c("y1", "y2")))
  sol <- solve_dsge(m1$model)
  ll1 <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D)$loglik
  ll3 <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D, init = m3$model$kalman_init)$loglik
  expect_equal(ll3, ll1, tolerance = 1e-9)
})

test_that("a beta prior on [lower, upper] is a shifted, scaled beta", {
  p <- prior("beta", shape1 = 3, shape2 = 12, lower = 1, upper = 2)
  expect_identical(p$support, "bounded")
  expect_equal(dsge:::dprior(p, 1.3), stats::dbeta(0.3, 3, 12, log = TRUE) - log(1))
  p2 <- prior("beta", shape1 = 3, shape2 = 12, lower = -1, upper = 3)
  expect_equal(dsge:::dprior(p2, 0.2), stats::dbeta(0.3, 3, 12, log = TRUE) - log(4))
  expect_identical(dsge:::dprior(p2, 3.5), -Inf)
  u <- dsge:::to_unconstrained(0.2, p2)
  expect_equal(dsge:::from_unconstrained(u, p2), 0.2)
  set.seed(1)
  expect_true(all(dsge:::rprior(p2, 200) > -1 & dsge:::rprior(p2, 200) < 3))
  expect_error(prior("beta", shape1 = 2, shape2 = 2, lower = 2, upper = 1), "lower < upper")
})

test_that("read_dynare() imports beta_pdf with a third and fourth parameter", {
  # Jermann and Quadrini (2012): eta_bar, beta_pdf with mean 1.2, sd 0.1 on [1, 2]
  m <- read_dynare(text = "var y; varexo e; parameters eta_bar; eta_bar = 1.2;
    model; y = eta_bar*y(-1)/2 + e; end;
    shocks; var e; stderr 1; end;
    estimated_params; eta_bar, beta_pdf, 1.2, 0.1, 1, 2; end;")
  p <- m$priors$eta_bar
  expect_equal(unlist(p$params[c("lower", "upper")]), c(lower = 1, upper = 2))
  st <- dsge:::compute_prior_stats(p)
  expect_equal(st$mean, 1.2)
  expect_equal(st$sd, 0.1)
  # Dynare's log density at 1.5 (lpdfgbeta)
  a <- p$params$shape1; b <- p$params$shape2
  expect_equal(dsge:::dprior(p, 1.5),
               -lbeta(a, b) + (a - 1) * log(0.5) + (b - 1) * log(0.5) - (a + b - 1) * log(1))
})
