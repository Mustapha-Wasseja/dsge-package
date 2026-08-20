# Tests for ms_filter() -- Markov-switching volatility (Kim 1994)

make_sol <- function(rho = 0.7) {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    fixed = list(rho = rho))
  solve_dsge(m, params = c(rho = rho), shock_sd = c(z = 1))
}

P2 <- function(p11 = 0.95, p22 = 0.90) {
  matrix(c(p11, 1 - p11, 1 - p22, p22), 2, 2, byrow = TRUE)
}


# --------------------------------------------------------------------
# Reduction to the single-regime Kalman filter
# --------------------------------------------------------------------

test_that("a single regime with unit scale reproduces the Kalman filter", {
  sol <- make_sol()
  set.seed(1)
  dat <- matrix(rnorm(150), 150, 1, dimnames = list(NULL, "y"))

  kf  <- kalman_filter(dat, sol$G, sol$H, sol$M, sol$D)
  ms  <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                   regime_scale = 1, P_trans = matrix(1, 1, 1))

  expect_equal(ms$loglik, kf$loglik, tolerance = 1e-10)
  expect_equal(ms$filtered_states, kf$filtered_states, tolerance = 1e-10)
  expect_true(all(abs(ms$filtered_probs - 1) < 1e-12))
})


test_that("two identical regimes also reproduce the Kalman filter", {
  # If both regimes have the same volatility the mixture is degenerate,
  # so the likelihood must be unchanged whatever the transition matrix.
  sol <- make_sol()
  set.seed(2)
  dat <- matrix(rnorm(120), 120, 1, dimnames = list(NULL, "y"))

  kf <- kalman_filter(dat, sol$G, sol$H, sol$M, sol$D)
  ms <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                  regime_scale = c(1, 1), P_trans = P2(0.8, 0.6))

  expect_equal(ms$loglik, kf$loglik, tolerance = 1e-10)
  expect_equal(ms$filtered_states, kf$filtered_states, tolerance = 1e-10)
})


# --------------------------------------------------------------------
# Probability bookkeeping
# --------------------------------------------------------------------

test_that("filtered, predicted and smoothed probabilities are valid", {
  sol <- make_sol()
  set.seed(3)
  dat <- matrix(rnorm(100), 100, 1, dimnames = list(NULL, "y"))
  out <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                   regime_scale = c(1, 3), P_trans = P2())

  for (nm in c("filtered_probs", "predicted_probs", "smoothed_probs")) {
    p <- out[[nm]]
    expect_true(all(p >= -1e-12 & p <= 1 + 1e-12), info = nm)
    expect_equal(rowSums(p), rep(1, nrow(p)), tolerance = 1e-10, info = nm)
  }
})


test_that("ergodic distribution solves pi' P = pi'", {
  P <- P2(0.95, 0.90)
  e <- dsge:::.ms_ergodic(P)
  expect_equal(sum(e), 1, tolerance = 1e-12)
  expect_equal(as.numeric(t(P) %*% e), e, tolerance = 1e-10)
  # Known closed form for a 2-state chain
  expect_equal(e[1], (1 - 0.90) / ((1 - 0.95) + (1 - 0.90)), tolerance = 1e-10)
})


# --------------------------------------------------------------------
# Does it actually detect a volatility regime change?
# --------------------------------------------------------------------

test_that("the filter recovers a known volatility break", {
  sol <- make_sol(rho = 0.5)
  H <- sol$H; G <- sol$G; M <- sol$M

  set.seed(123)
  TT <- 400L
  # First half calm (scale 1), second half turbulent (scale 4)
  true_regime <- c(rep(1L, TT / 2), rep(2L, TT / 2))
  scales <- c(1, 4)
  x <- 0; yv <- numeric(TT)
  for (t in seq_len(TT)) {
    e <- rnorm(1) * scales[true_regime[t]]
    x <- as.numeric(H %*% x + M %*% e)
    yv[t] <- as.numeric(G %*% x)
  }
  dat <- matrix(yv, TT, 1, dimnames = list(NULL, "y"))

  out <- ms_filter(dat, G, H, M, sol$D,
                   regime_scale = scales, P_trans = P2(0.98, 0.98))

  # Smoothed probability of the high-volatility regime should be low in
  # the calm half and high in the turbulent half.
  p_high <- out$smoothed_probs[, 2]
  expect_lt(mean(p_high[1:(TT / 2)]), 0.35)
  expect_gt(mean(p_high[(TT / 2 + 1):TT]), 0.65)

  # Classification accuracy against the truth
  guess <- ifelse(p_high > 0.5, 2L, 1L)
  expect_gt(mean(guess == true_regime), 0.8)
})


test_that("the correct volatility scales beat badly wrong ones", {
  sol <- make_sol(rho = 0.5)
  H <- sol$H; G <- sol$G; M <- sol$M
  set.seed(321)
  TT <- 300L
  true_regime <- c(rep(1L, TT / 2), rep(2L, TT / 2))
  x <- 0; yv <- numeric(TT)
  for (t in seq_len(TT)) {
    e <- rnorm(1) * c(1, 4)[true_regime[t]]
    x <- as.numeric(H %*% x + M %*% e)
    yv[t] <- as.numeric(G %*% x)
  }
  dat <- matrix(yv, TT, 1, dimnames = list(NULL, "y"))

  ll_right <- ms_filter(dat, G, H, M, sol$D, regime_scale = c(1, 4),
                        P_trans = P2(0.98, 0.98))$loglik
  ll_flat  <- ms_filter(dat, G, H, M, sol$D, regime_scale = c(1, 1),
                        P_trans = P2(0.98, 0.98))$loglik
  expect_gt(ll_right, ll_flat)
})


# --------------------------------------------------------------------
# Interface details
# --------------------------------------------------------------------

test_that("per-shock scale matrix is accepted on a multi-shock model", {
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
  set.seed(4)
  dat <- matrix(rnorm(160), 80, 2, dimnames = list(NULL, c("p", "r")))

  sc <- rbind(c(1, 1), c(3, 0.5))    # regime 2: shock 1 up, shock 2 down
  out <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                   regime_scale = sc, P_trans = P2())
  expect_s3_class(out, "dsge_ms_filter")
  expect_true(is.finite(out$loglik))
  expect_equal(dim(out$filtered_probs), c(80L, 2L))
})


test_that("smooth = FALSE skips the smoother", {
  sol <- make_sol()
  set.seed(5)
  dat <- matrix(rnorm(60), 60, 1, dimnames = list(NULL, "y"))
  out <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                   regime_scale = c(1, 2), P_trans = P2(),
                   smooth = FALSE)
  expect_null(out$smoothed_probs)
})


test_that("input validation", {
  sol <- make_sol()
  dat <- matrix(rnorm(30), 30, 1, dimnames = list(NULL, "y"))
  expect_error(ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                         regime_scale = c(1, 2),
                         P_trans = matrix(1, 1, 1)),
               "to match the number of regimes")
  expect_error(ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                         regime_scale = c(1, 2),
                         P_trans = matrix(c(0.5, 0.4, 0.1, 0.9), 2, 2,
                                          byrow = TRUE)),
               "must sum to 1")
  expect_error(ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                         regime_scale = c(1, -2), P_trans = P2()),
               "strictly positive")
  expect_error(ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                         regime_scale = c(1, 2), P_trans = P2(),
                         initial_probs = c(0.5, 0.9)),
               "sum to 1")
})


test_that("print and plot methods run", {
  sol <- make_sol()
  set.seed(6)
  dat <- matrix(rnorm(60), 60, 1, dimnames = list(NULL, "y"))
  out <- ms_filter(dat, sol$G, sol$H, sol$M, sol$D,
                   regime_scale = c(1, 2.5), P_trans = P2())
  expect_output(print(out), "Markov-switching volatility filter")
  expect_silent({ pdf(tempfile()); plot(out); dev.off() })
})
