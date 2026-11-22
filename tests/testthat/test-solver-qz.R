# The ordered QZ solve, the first of klein_solve()'s methods, against the
# fallbacks and known solutions.

qz_inputs <- function(m, params) {
  list(sm = dsge:::build_structural_matrices(m, c(params, unlist(m$fixed))),
       n_c = m$n_controls, n_s = m$n_states)
}
run_method <- function(f, x) {
  sm <- x$sm
  f(sm$A0, sm$A1, sm$A2, sm$A3, sm$B0, sm$B1, sm$B2, sm$B3, NULL,
    x$n_c, x$n_s)
}

test_that("the QZ solve matches cyclic reduction and the spectral divide", {
  nk <- dsge_model(
    obs(pi ~ beta * lead(pi) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(pi) - g)),
    obs(r ~ psi * pi + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99, kappa = 0.1, psi = 1.5),
    start = list(rhou = 0.5, rhog = 0.8)
  )
  x <- qz_inputs(nk, c(rhou = 0.5, rhog = 0.8))
  qz <- run_method(dsge:::klein_qz_solve, x)
  cr <- run_method(dsge:::klein_cyclic_reduction, x)
  sd <- run_method(dsge:::klein_spectral_divide, x)
  expect_false(is.null(qz))
  expect_false(is.null(cr))
  expect_false(is.null(sd))
  expect_equal(unname(qz$G), unname(cr$G), tolerance = 1e-10)
  expect_equal(unname(qz$H), unname(cr$H), tolerance = 1e-10)
  expect_equal(unname(qz$G), unname(sd$G), tolerance = 1e-10)
})

test_that("the QZ solve counts a unit root as stable, as Dynare does", {
  # y_t = beta E_t y_{t+1} + z_t with z a random walk: G = 1 / (1 - beta)
  m <- dsge_model(
    obs(y ~ beta * lead(y) + z),
    state(z ~ rho * z),
    fixed = list(beta = 0.9),
    start = list(rho = 1)
  )
  qz <- run_method(dsge:::klein_qz_solve, qz_inputs(m, c(rho = 1)))
  expect_false(is.null(qz))
  expect_equal(as.numeric(qz$G), 10, tolerance = 1e-10)
  expect_equal(as.numeric(qz$H), 1, tolerance = 1e-12)
  sol <- solve_dsge(m, params = c(rho = 1))
  expect_true(sol$stable)
  expect_equal(as.numeric(sol$G), 10, tolerance = 1e-10)
})

test_that("without a unique stable solution the QZ solve leaves it to the fallbacks", {
  m <- dsge_model(
    obs(y ~ beta * lead(y) + z),
    state(z ~ rho * z),
    fixed = list(beta = 0.9),
    start = list(rho = 0.5)
  )
  # an explosive state: one stable root too few
  expect_null(run_method(dsge:::klein_qz_solve, qz_inputs(m, c(rho = 1.5))))
  expect_false(solve_dsge(m, params = c(rho = 1.5))$stable)
  # beta > 1 makes the forward-looking root stable too: one too many
  m2 <- dsge_model(
    obs(y ~ beta * lead(y) + z),
    state(z ~ rho * z),
    fixed = list(beta = 1.5),
    start = list(rho = 0.5)
  )
  expect_null(run_method(dsge:::klein_qz_solve, qz_inputs(m2, c(rho = 0.5))))
})
