# tests/testthat/test-third-order.R

build_nl_ar1_3 <- function(rhoA = 0.90) {
  dsgenl_model(
    "Y = A",
    "A(+1) = A^rhoA",
    observed   = "Y",
    unobserved = character(0),
    exo_state  = "A",
    fixed      = list(),
    start      = list(rhoA = rhoA),
    ss_guess   = c(Y = 1, A = 1)
  )
}

build_nl_consumption_3 <- function() {
  dsgenl_model(
    "Y = A * H",
    "1 / C = beta * R / C(+1)",
    "Y = C",
    "R = 1 / beta",
    "A(+1) = A^rhoA",
    observed   = "Y",
    unobserved = c("C", "H", "R"),
    exo_state  = "A",
    fixed      = list(beta = 0.99),
    start      = list(rhoA = 0.90),
    ss_guess   = c(Y = 1, C = 1, H = 1, R = 1.0101, A = 1)
  )
}

test_that("solve_dsge with order=3 returns third-order solution fields", {
  mod <- build_nl_ar1_3()
  sol <- solve_dsge(mod, params = c(rhoA = 0.90),
                    shock_sd = c(A = 0.01), order = 3)

  expect_s3_class(sol, "dsge_solution")
  expect_equal(sol$order, 3L)
  # Must have second-order fields
  expect_false(is.null(sol$g_xx))
  expect_false(is.null(sol$h_xx))
  expect_false(is.null(sol$g_ss))
  expect_false(is.null(sol$h_ss))
  # Must have third-order fields
  expect_false(is.null(sol$g_xxx))
  expect_false(is.null(sol$h_xxx))
  expect_false(is.null(sol$g_xss))
  expect_false(is.null(sol$h_xss))
  expect_false(is.null(sol$g_sss))
  expect_false(is.null(sol$h_sss))
})

test_that("third-order solution preserves first-order matrices", {
  mod <- build_nl_ar1_3()
  sol1 <- solve_dsge(mod, params = c(rhoA = 0.90),
                     shock_sd = c(A = 0.01), order = 1)
  sol3 <- solve_dsge(mod, params = c(rhoA = 0.90),
                     shock_sd = c(A = 0.01), order = 3)

  expect_equal(sol3$G, sol1$G, tolerance = 1e-5)
  expect_equal(sol3$H, sol1$H, tolerance = 1e-5)
  expect_equal(sol3$M, sol1$M, tolerance = 1e-5)
})

test_that("g_xxx and h_xxx have correct dimensions", {
  mod <- build_nl_ar1_3()
  sol <- solve_dsge(mod, params = c(rhoA = 0.90),
                    shock_sd = c(A = 0.01), order = 3)
  n_c <- nrow(sol$G); n_s <- ncol(sol$G)

  expect_equal(dim(sol$g_xxx), c(n_c, n_s, n_s, n_s))
  expect_equal(dim(sol$h_xxx), c(n_s, n_s, n_s, n_s))
  expect_equal(dim(sol$g_xss), c(n_c, n_s))
  expect_equal(dim(sol$h_xss), c(n_s, n_s))
  expect_equal(length(sol$g_sss), n_c)
  expect_equal(length(sol$h_sss), n_s)
})

test_that("third-order coefficients are finite", {
  mod <- build_nl_ar1_3()
  sol <- solve_dsge(mod, params = c(rhoA = 0.90),
                    shock_sd = c(A = 0.01), order = 3)

  expect_true(all(is.finite(sol$g_xxx)))
  expect_true(all(is.finite(sol$h_xxx)))
  expect_true(all(is.finite(sol$g_ss)))
  expect_true(all(is.finite(sol$h_ss)))
  expect_true(all(is.finite(sol$g_sss)))
  expect_true(all(is.finite(sol$h_sss)))
})

test_that("order=3 on linear model errors", {
  m <- dsge_model(
    obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.9)
  )
  expect_error(solve_dsge(m, params = c(rho = 0.9), order = 3),
               "dsgenl_model")
})

test_that("order=4 errors with updated message", {
  mod <- build_nl_ar1_3()
  expect_error(solve_dsge(mod, params = c(rhoA = 0.9), order = 4),
               "order must be 1, 2, or 3")
})

test_that("simulate_3rd_order produces output of correct dimensions", {
  mod <- build_nl_consumption_3()
  sol <- suppressWarnings(
    solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
               shock_sd = c(A = 0.01), order = 3)
  )
  sim <- simulate_3rd_order(sol, n = 30, n_burn = 10, seed = 7)

  expect_equal(nrow(sim$states),   30L)
  expect_equal(nrow(sim$controls), 30L)
  expect_equal(sim$order, 3L)
  expect_false(is.null(sim$state_levels))
})

test_that("simulate_3rd_order errors on non-third-order solution", {
  mod <- build_nl_ar1_3()
  sol2 <- solve_dsge(mod, params = c(rhoA = 0.90),
                     shock_sd = c(A = 0.01), order = 2)
  expect_error(simulate_3rd_order(sol2), "order = 3")
})

test_that("simulate_3rd_order is reproducible", {
  mod <- build_nl_consumption_3()
  sol <- suppressWarnings(
    solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
               shock_sd = c(A = 0.01), order = 3)
  )
  s1 <- simulate_3rd_order(sol, n = 20, n_burn = 5, seed = 42)
  s2 <- simulate_3rd_order(sol, n = 20, n_burn = 5, seed = 42)
  expect_equal(s1$states, s2$states)
  expect_equal(s1$controls, s2$controls)
})

test_that("third-order adds correction beyond second-order (non-trivial nonlinearity)", {
  mod <- build_nl_consumption_3()
  sol2 <- suppressWarnings(
    solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
               shock_sd = c(A = 0.01), order = 2)
  )
  sol3 <- suppressWarnings(
    solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
               shock_sd = c(A = 0.01), order = 3)
  )
  # g_xxx should not all be zero for a model with genuine nonlinearity
  expect_false(all(sol3$g_xxx == 0))
})
