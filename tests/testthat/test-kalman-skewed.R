# Tests for kalman_filter_skewed()

make_sol <- function(rho = 0.8) {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    fixed = list(rho = rho))
  solve_dsge(m, params = c(rho = rho), shock_sd = c(z = 1))
}

# Draw a standardised (mean 0, variance 1) skew-normal with skewness g
rsn_std <- function(n, g) {
  p <- dsge:::.sn_params_from_skew(g)
  delta <- p$alpha / sqrt(1 + p$alpha^2)
  u0 <- abs(stats::rnorm(n))
  u1 <- stats::rnorm(n)
  p$xi + p$omega * (delta * u0 + sqrt(1 - delta^2) * u1)
}


# --------------------------------------------------------------------
# Core mathematics
# --------------------------------------------------------------------

test_that(".cum3_map matches a brute-force triple sum", {
  set.seed(1)
  A <- matrix(rnorm(6), 3, 2)
  K <- array(rnorm(8), c(2, 2, 2))
  out <- dsge:::.cum3_map(A, K)

  ref <- array(0, c(3, 3, 3))
  for (i in 1:3) for (j in 1:3) for (k in 1:3) {
    s <- 0
    for (a in 1:2) for (b in 1:2) for (cc in 1:2)
      s <- s + A[i, a] * A[j, b] * A[k, cc] * K[a, b, cc]
    ref[i, j, k] <- s
  }
  expect_equal(out, ref, tolerance = 1e-12)
})


test_that(".sn_params_from_skew inverts the skewness formula", {
  for (g in c(-0.9, -0.5, -0.1, 0.1, 0.5, 0.9)) {
    p <- dsge:::.sn_params_from_skew(g)
    delta <- p$alpha / sqrt(1 + p$alpha^2)
    b <- sqrt(2 / pi)
    mu <- p$xi + p$omega * delta * b
    v  <- p$omega^2 * (1 - 2 * delta^2 / pi)
    sk <- ((4 - pi) / 2) * (delta * b)^3 / (1 - 2 * delta^2 / pi)^(3 / 2)
    expect_equal(mu, 0, tolerance = 1e-10)
    expect_equal(v,  1, tolerance = 1e-10)
    expect_equal(sk, g, tolerance = 1e-8)
  }
})


test_that("standardised skew-normal density integrates to the right moments", {
  g <- 0.6
  f <- function(z) exp(vapply(z, dsge:::.dsn_log_std, numeric(1), g = g))
  expect_equal(stats::integrate(f, -12, 12)$value, 1, tolerance = 1e-6)
  expect_equal(stats::integrate(function(z) z * f(z), -12, 12)$value,
               0, tolerance = 1e-6)
  expect_equal(stats::integrate(function(z) z^2 * f(z), -12, 12)$value,
               1, tolerance = 1e-6)
  expect_equal(stats::integrate(function(z) z^3 * f(z), -12, 12)$value,
               g, tolerance = 1e-6)
})


test_that("zero-skew density equals the standard normal density", {
  zs <- c(-2.5, -1, 0, 0.4, 3)
  expect_equal(vapply(zs, dsge:::.dsn_log_std, numeric(1), g = 0),
               stats::dnorm(zs, log = TRUE), tolerance = 1e-12)
})


# --------------------------------------------------------------------
# Reduction to the Gaussian filter
# --------------------------------------------------------------------

test_that("zero skewness reproduces the Gaussian Kalman filter exactly", {
  sol <- make_sol()
  set.seed(42)
  dat <- matrix(rnorm(120), 120, 1, dimnames = list(NULL, "y"))

  kf  <- kalman_filter(dat, sol$G, sol$H, sol$M, sol$D)
  kfs <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                              shock_skew = 0)

  expect_equal(kfs$loglik, kf$loglik, tolerance = 1e-10)
  expect_equal(kfs$filtered_states,   kf$filtered_states,   tolerance = 1e-10)
  expect_equal(kfs$predicted_states,  kf$predicted_states,  tolerance = 1e-10)
  expect_equal(kfs$prediction_errors, kf$prediction_errors, tolerance = 1e-10)
  # Innovations carry no skewness when the shocks are Gaussian
  expect_true(all(abs(kfs$innovation_skew) < 1e-12))
})


test_that("filtered means/covariances do not depend on the skewness", {
  # The Kalman gain is the optimal linear filter regardless of shape,
  # so only the likelihood should change with `shock_skew`.
  sol <- make_sol()
  set.seed(7)
  dat <- matrix(rnorm(80), 80, 1, dimnames = list(NULL, "y"))
  a <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D, shock_skew = 0)
  b <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D, shock_skew = 0.8)
  expect_equal(a$filtered_states, b$filtered_states, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(a$loglik, b$loglik)))
})


# --------------------------------------------------------------------
# Does it actually detect skewness?
# --------------------------------------------------------------------

test_that("the likelihood prefers the correct skew sign on skewed data", {
  sol <- make_sol(rho = 0.5)
  H <- sol$H; G <- sol$G; M <- sol$M
  g_true <- -0.85

  set.seed(11)
  TT <- 400L
  x <- 0; y <- numeric(TT)
  eps <- rsn_std(TT, g_true)
  for (t in seq_len(TT)) {
    x <- as.numeric(H %*% x + M %*% eps[t])
    y[t] <- as.numeric(G %*% x)
  }
  dat <- matrix(y, TT, 1, dimnames = list(NULL, "y"))

  ll_right <- kalman_filter_skewed(dat, G, H, M, sol$D,
                                   shock_skew = g_true)$loglik
  ll_wrong <- kalman_filter_skewed(dat, G, H, M, sol$D,
                                   shock_skew = -g_true)$loglik
  ll_gauss <- kalman_filter_skewed(dat, G, H, M, sol$D,
                                   shock_skew = 0)$loglik

  expect_gt(ll_right, ll_wrong)   # correct sign beats the mirror image
  expect_gt(ll_right, ll_gauss)   # and beats assuming Gaussian shocks
})


test_that("innovation skewness has the sign of the shock skewness", {
  sol <- make_sol(rho = 0.5)
  set.seed(3)
  dat <- matrix(rnorm(50), 50, 1, dimnames = list(NULL, "y"))
  neg <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                              shock_skew = -0.7)
  pos <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                              shock_skew = 0.7)
  expect_true(all(neg$innovation_skew <= 0))
  expect_true(all(pos$innovation_skew >= 0))
  expect_equal(neg$innovation_skew, -pos$innovation_skew, tolerance = 1e-10)
})


# --------------------------------------------------------------------
# Multi-shock model and input validation
# --------------------------------------------------------------------

test_that("works on a multi-shock, multi-observable model", {
  nk <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
  sol <- solve_dsge(nk,
    params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd = c(e.u = 1, e.g = 0.5))
  set.seed(5)
  dat <- matrix(rnorm(200), 100, 2, dimnames = list(NULL, c("p", "r")))

  kf  <- kalman_filter(dat, sol$G, sol$H, sol$M, sol$D)
  z   <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                              shock_skew = c(0, 0))
  expect_equal(z$loglik, kf$loglik, tolerance = 1e-10)

  # Per-shock skewness vector is accepted
  s <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                            shock_skew = c(-0.6, 0.3))
  expect_true(is.finite(s$loglik))
  expect_equal(dim(s$innovation_skew), c(100L, 2L))
})


test_that("input validation", {
  sol <- make_sol()
  dat <- matrix(rnorm(20), 20, 1, dimnames = list(NULL, "y"))
  expect_error(kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D),
               "required")
  expect_error(kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                                    shock_skew = c(0.1, 0.2)),
               "length 1 or")
  expect_error(kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                                    shock_skew = NA_real_),
               "finite")
})


test_that("extreme skewness is clamped rather than failing", {
  sol <- make_sol()
  set.seed(9)
  dat <- matrix(rnorm(40), 40, 1, dimnames = list(NULL, "y"))
  out <- kalman_filter_skewed(dat, sol$G, sol$H, sol$M, sol$D,
                              shock_skew = 5)
  expect_equal(out$shock_skew, 0.99)
  expect_true(is.finite(out$loglik))
})
