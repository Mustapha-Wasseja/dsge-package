# Tests for PAC (polynomial adjustment cost) equations

# Apply the Euler operator  [1 + sum_i d_i (1-L)^i (1-beta F)^i]  to a
# deterministic path `y` at time index `t` (perfect foresight, so
# E_t y_{t+q} = y[t+q]).
euler_lhs <- function(y, t, beta, d) {
  out <- y[t]
  for (i in seq_along(d)) {
    if (d[i] == 0) next
    acc <- 0
    for (p in 0:i) for (q in 0:i) {
      acc <- acc + choose(i, p) * (-1)^p *
                   choose(i, q) * (-beta)^q * y[t - p + q]
    }
    out <- out + d[i] * acc
  }
  out
}


# --------------------------------------------------------------------
# The m = 1 closed form
# --------------------------------------------------------------------

test_that("m = 1 stable root matches the analytic quadratic root", {
  for (bet in c(0.95, 0.99, 1)) for (dd in c(0.25, 1.5, 8)) {
    p <- pac_weights(beta = bet, d = dd)
    # d beta z^2 - (1 + d + d beta) z + d = 0, smaller root
    A <- dd * bet; B <- -(1 + dd + dd * bet); C <- dd
    analytic <- (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    expect_equal(Re(p$roots[1]), analytic, tolerance = 1e-8,
                 info = paste("beta", bet, "d", dd))
  }
})


test_that("m = 1 weights match the textbook closed form", {
  bet <- 0.99; dd <- 2
  p   <- pac_weights(beta = bet, d = dd, horizon = 40)
  phi <- Re(p$roots[1])
  expect_equal(p$lag_coef[1], phi, tolerance = 1e-10)
  j   <- 0:40
  expect_equal(p$fwd_weight,
               (1 - phi) * (1 - bet * phi) * (bet * phi)^j,
               tolerance = 1e-10)
})


# --------------------------------------------------------------------
# THE decisive test: does the solution satisfy the Euler equation?
# --------------------------------------------------------------------

test_that("the PAC solution satisfies the original Euler equation", {
  for (spec in list(list(b = 0.99, d = 1.5),
                    list(b = 0.97, d = 5),
                    list(b = 0.99, d = c(1.0, 0.5)),
                    list(b = 0.98, d = c(2.0, 0.3, 0.1)))) {
    bet <- spec$b; dd <- spec$d
    p <- pac_weights(beta = bet, d = dd, horizon = 400)

    # AR(1) target decaying to zero, so terminal padding is harmless
    TT  <- 300L
    rho <- 0.85
    ystar <- rho^(seq_len(TT) - 1L)
    y <- pac_simulate(p, ystar)

    # Residual of  [1 + sum d_i (1-L)^i (1-beta F)^i] y_t - ystar_t
    m <- length(dd)
    idx <- (m + 5L):(TT - m - 60L)     # interior points only
    resid <- vapply(idx, function(t)
      euler_lhs(y, t, bet, dd) - ystar[t], numeric(1))

    expect_lt(max(abs(resid)) / max(abs(ystar[idx])), 1e-6,
              label = paste("relative Euler residual, d =",
                            paste(dd, collapse = ",")))
  }
})


# --------------------------------------------------------------------
# Long-run properties
# --------------------------------------------------------------------

test_that("long-run homogeneity holds: lag + forward weights sum to 1", {
  for (dd in list(1, 5, c(1, 0.5), c(2, 0.3, 0.1))) {
    p <- pac_weights(beta = 0.99, d = dd, horizon = 800)
    expect_equal(p$homogeneity, 1, tolerance = 1e-6,
                 info = paste(dd, collapse = ","))
  }
})


test_that("a permanent step in the target is eventually matched one-for-one", {
  p <- pac_weights(beta = 0.99, d = 3, horizon = 400)
  y <- pac_simulate(p, ystar = c(rep(0, 5), rep(1, 400)))
  expect_equal(tail(y, 1), 1, tolerance = 1e-5)

  # Adjustment is gradual and, under perfect foresight, begins *before*
  # the step lands: the agent already knows it is coming.  Five periods
  # ahead the analytic anticipatory response is (1-phi)(beta*phi)^5.
  phi <- Re(p$roots[1])
  expect_equal(y[1], (1 - phi) * (0.99 * phi)^5, tolerance = 1e-8)
  expect_lt(y[1], 0.05)                      # only a small early move
  # On impact the response is phi*y[5] + (1-phi)
  expect_equal(y[6], phi * y[5] + (1 - phi), tolerance = 1e-8)
  expect_lt(y[6], 1)                          # still incomplete
  expect_true(all(diff(y[1:300]) > -1e-12))  # monotone approach
})


test_that("costless adjustment collapses to y_t = ystar_t", {
  p <- pac_weights(beta = 0.99, d = 1e-8)
  expect_lt(abs(p$lag_coef[1]), 1e-3)
  expect_equal(p$fwd_weight[1], 1, tolerance = 1e-3)
  y <- pac_simulate(p, ystar = c(0, 1, 3, 2))
  expect_equal(y, c(0, 1, 3, 2), tolerance = 1e-3)
})


test_that("larger adjustment costs give more sluggish adjustment", {
  phis <- vapply(c(0.5, 2, 10, 100),
                 function(dd) Re(pac_weights(0.99, dd)$roots[1]),
                 numeric(1))
  expect_true(all(diff(phis) > 0))          # increasing in d
  expect_lt(phis[4], 1)
  expect_gt(phis[4], 0.8)                   # approaching 1
})


# --------------------------------------------------------------------
# Closed-form forward sum on a state process
# --------------------------------------------------------------------

test_that("pac_target_loading matches the truncated forward sum (AR(1))", {
  p   <- pac_weights(beta = 0.99, d = 2, horizon = 600)
  rho <- 0.8
  Psi <- matrix(rho, 1, 1)
  loading <- pac_target_loading(p, Psi, e = 1)
  brute   <- sum(p$fwd_weight * rho^(seq_along(p$fwd_weight) - 1L))
  expect_equal(loading, brute, tolerance = 1e-8)
})


test_that("pac_target_loading works for a multivariate state", {
  p <- pac_weights(beta = 0.99, d = c(1.5, 0.4), horizon = 800)
  Psi <- matrix(c(0.7, 0.1, 0.0, 0.5), 2, 2, byrow = TRUE)
  e   <- c(1, 0.5)
  loading <- pac_target_loading(p, Psi, e)

  # Brute force: sum_j w_j e' Psi^j
  acc <- rep(0, 2); Pk <- diag(2)
  for (j in seq_along(p$fwd_weight)) {
    acc <- acc + p$fwd_weight[j] * as.numeric(e %*% Pk)
    Pk <- Pk %*% Psi
  }
  expect_equal(loading, acc, tolerance = 1e-7)
})


test_that("a constant target gives a loading equal to the total weight", {
  # Psi = I means the target never changes, so the discounted sum of
  # weights should equal the total forward weight.
  p <- pac_weights(beta = 0.99, d = 2, horizon = 2000)
  loading <- pac_target_loading(p, Psi = matrix(1, 1, 1), e = 1)
  expect_equal(loading, sum(p$fwd_weight), tolerance = 1e-5)
  # And by homogeneity that equals 1 - sum(a)
  expect_equal(loading, 1 - sum(p$lag_coef), tolerance = 1e-5)
})


# --------------------------------------------------------------------
# Interface
# --------------------------------------------------------------------

test_that("second-order costs produce two lags and two roots", {
  p <- pac_weights(beta = 0.99, d = c(1.0, 0.5))
  expect_equal(p$m, 2L)
  expect_equal(length(p$roots), 2L)
  expect_equal(length(p$lag_coef), 2L)
  expect_true(all(Mod(p$roots) < 1))
})


test_that("pac_simulate honours initial lags", {
  p <- pac_weights(beta = 0.99, d = 2)
  y0 <- pac_simulate(p, ystar = rep(0, 20), y_init = 1)
  # With a zero target the path decays from the initial condition
  expect_equal(y0[1], p$lag_coef[1], tolerance = 1e-10)
  expect_true(all(abs(diff(abs(y0))) <= 1e-9 | diff(abs(y0)) < 0))
  expect_lt(abs(y0[20]), abs(y0[1]))
})


test_that("print method runs", {
  p <- pac_weights(beta = 0.99, d = c(1, 0.5))
  expect_output(print(p), "Polynomial Adjustment Cost")
  expect_output(print(p), "homogeneity")
})


test_that("input validation", {
  expect_error(pac_weights(beta = 1.5, d = 1),  "in \\(0, 1\\]")
  expect_error(pac_weights(beta = 0.99, d = -1), "non-negative")
  expect_error(pac_weights(beta = 0.99, d = 0),  "must be positive")
  p <- pac_weights(beta = 0.99, d = 1)
  expect_error(pac_target_loading(p, matrix(1:6, 2, 3), e = 1), "square")
  expect_error(pac_target_loading(p, matrix(0.5, 1, 1), e = c(1, 2)),
               "length nrow")
  expect_error(pac_simulate(p, ystar = numeric(0)), "non-empty")
  expect_error(pac_simulate(p, ystar = 1:5, y_init = c(1, 2)), "length m")
  expect_error(pac_target_loading("not a pac", matrix(1), 1), "dsge_pac")
})
