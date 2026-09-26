# The tensor algebra of the higher-order perturbation runs in C++
# (src/perturbation.cpp) on sparse derivative entries. These tests compare
# it with direct dense computations in R.

random_sparse_hessians <- function(n_eq, n, density = 0.3) {
  H <- lapply(seq_len(n_eq), function(k) {
    A <- matrix(stats::rnorm(n * n), n) * (stats::runif(n * n) < density)
    A + t(A)
  })
  entries <- do.call(rbind, lapply(seq_len(n_eq), function(k) {
    w <- which(H[[k]] != 0, arr.ind = TRUE)
    data.frame(eq = rep(k, nrow(w)), i = w[, 1L], j = w[, 2L],
               val = H[[k]][w])
  }))
  list(H = H, entries = entries)
}

test_that("contract2 gives vec(P' H_k Q) for every equation", {
  set.seed(1)
  h <- random_sparse_hessians(5, 7)
  P <- matrix(stats::rnorm(7 * 3), 7)
  Q <- matrix(stats::rnorm(7 * 4), 7)
  out <- dsge:::.contract2(h$entries, P, Q, 5L)
  ref <- t(vapply(h$H, function(Hk) as.vector(crossprod(P, Hk %*% Q)),
                  numeric(12)))
  expect_equal(out, ref, tolerance = 1e-12)
})

test_that("contract3 contracts the third-derivative tensor", {
  set.seed(2)
  n <- 4
  t3 <- data.frame(eq = c(1L, 1L, 2L, 3L), i = c(1L, 2L, 3L, 4L),
                   j = c(2L, 2L, 1L, 4L), l = c(3L, 1L, 4L, 4L),
                   val = stats::rnorm(4))
  P <- matrix(stats::rnorm(n * 2), n)
  Q <- matrix(stats::rnorm(n * 3), n)
  R <- matrix(stats::rnorm(n * 2), n)
  out <- dsge:::.contract3(t3, P, Q, R, 3L)
  ref <- matrix(0, 3, 2 * 3 * 2)
  for (r in seq_len(nrow(t3))) {
    ref[t3$eq[r], ] <- ref[t3$eq[r], ] + t3$val[r] *
      as.vector(P[t3$i[r], ] %o% Q[t3$j[r], ] %o% R[t3$l[r], ])
  }
  expect_equal(out, ref, tolerance = 1e-12)
})

test_that("apply_hx equals multiplication by the Kronecker power of hx", {
  set.seed(3)
  n <- 3
  hx <- matrix(stats::rnorm(n * n), n)
  for (k in 1:3) {
    X <- matrix(stats::rnorm(4 * n^k), 4)
    K <- hx
    if (k >= 2L) for (m in 2:k) K <- kronecker(hx, K)
    expect_equal(dsge:::.apply_hx(X, hx, k), X %*% K, tolerance = 1e-12)
  }
})

test_that("the C++ doubling matches the direct Kronecker solution", {
  set.seed(4)
  n_eq <- 6; n <- 3
  A <- diag(n_eq) + 0.1 * matrix(stats::rnorm(n_eq^2), n_eq)
  B <- 0.3 * matrix(stats::rnorm(n_eq^2), n_eq)
  hx <- 0.5 * matrix(stats::rnorm(n * n), n)
  hx <- hx * 0.8 / max(abs(eigen(hx, only.values = TRUE)$values))
  for (k in 1:3) {
    D <- matrix(stats::rnorm(n_eq * n^k), n_eq)
    direct <- dsge:::.solve_gen_sylvester(A, B, hx, D, k, max_direct = 1e6)
    doubling <- dsge:::.solve_gen_sylvester(A, B, hx, D, k, max_direct = 0)
    expect_equal(doubling, direct, tolerance = 1e-9)
  }
})

test_that("derivative expressions are cached as plain calls, not byte-code", {
  # byte-compiling large derivative expressions took tens of seconds
  m <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  sol <- solve_dsge(m, order = 3)
  expect_true(is.call(m$model$.cache$jac_code$code))
  expect_true(is.call(m$model$.cache$higher_code$h_code))
  expect_true(is.call(m$model$.cache$higher_code$t_code))
})
