test_that("AR(1) model solves correctly", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )

  sol <- solve_dsge(m, params = c(rho = 0.8))

  expect_s3_class(sol, "dsge_solution")
  expect_true(sol$stable)

  # For y = z with state z_{t+1} = rho * z_t + e_{t+1}:
  # Policy matrix G should be [1] (y = 1 * z)
  # Transition matrix H should be [0.8] (z_{t+1} = 0.8 * z_t)
  expect_equal(as.numeric(sol$G), 1, tolerance = 1e-10)
  expect_equal(as.numeric(sol$H), 0.8, tolerance = 1e-10)
})

test_that("unstable parameter returns stable=FALSE", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.5)
  )

  # rho > 1 means the state is explosive
  sol <- solve_dsge(m, params = c(rho = 1.5))

  expect_false(sol$stable)
})

test_that("two-state AR model solves correctly", {
  m <- dsge_model(
    obs(y ~ u),
    obs(p ~ g),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    start = list(rhou = 0.7, rhog = 0.9)
  )

  sol <- solve_dsge(m, params = c(rhou = 0.7, rhog = 0.9))

  expect_true(sol$stable)

  # G should be identity-like (y = u, p = g)
  # H should be diagonal with rhou, rhog
  expect_equal(sol$G[1, 1], 1, tolerance = 1e-10)
  expect_equal(sol$G[2, 2], 1, tolerance = 1e-10)
  expect_equal(sol$H[1, 1], 0.7, tolerance = 1e-10)
  expect_equal(sol$H[2, 2], 0.9, tolerance = 1e-10)
})

test_that("solve_dsge rejects missing parameters", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z)
  )

  expect_error(solve_dsge(m, params = c()), "Missing parameter")
})
