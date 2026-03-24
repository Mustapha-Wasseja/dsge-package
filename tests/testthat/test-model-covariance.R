# tests/testthat/test-model-covariance.R

# Helper: create a simple AR(1) model and solve it
# y_t = rho * y_{t-1} + e_t, implemented as obs(y ~ lead(y) + u) with state(u ~ rho * u)
make_ar1_solution <- function(rho = 0.5, sigma = 0.1) {
  m <- dsge_model(
    obs(y ~ lead(y) + u),
    state(u ~ rho_u * u),
    fixed = list(),
    start = list(rho_u = rho)
  )
  solve_dsge(m, params = c(rho_u = rho), shock_sd = c(u = sigma))
}

test_that("model_covariance works on dsge_solution", {
  sol <- make_ar1_solution(0.5, 0.1)
  mc <- model_covariance(sol)

  expect_s3_class(mc, "dsge_covariance")
  expect_true(is.matrix(mc$covariance))
  expect_true(is.matrix(mc$correlation))

  # Symmetric
  expect_equal(mc$covariance, t(mc$covariance))
  expect_equal(mc$correlation, t(mc$correlation))

  # Positive diagonal
  expect_true(all(diag(mc$covariance) >= 0))

  # Correlation diagonal = 1
  expect_equal(unname(diag(mc$correlation)), rep(1, length(mc$variables)))

  # Standard deviations consistent
  expect_equal(mc$std_dev^2, diag(mc$covariance), tolerance = 1e-10)
})

test_that("model_covariance with autocovariances", {
  sol <- make_ar1_solution(0.5, 0.1)
  mc <- model_covariance(sol, n_lags = 3)

  expect_equal(mc$n_lags, 3)
  expect_length(mc$autocovariances, 3)
  expect_true("lag_1" %in% names(mc$autocovariances))
  expect_true("lag_3" %in% names(mc$autocovariances))

  # Autocovariances should decay (AR(1) model)
  acov1 <- abs(mc$autocovariances$lag_1[1, 1])
  acov3 <- abs(mc$autocovariances$lag_3[1, 1])
  expect_true(acov1 >= acov3)
})

test_that("model_covariance AR(1) analytical check", {
  # State u follows AR(1): u_{t+1} = 0.9*u_t + e_t, var(e)=0.01
  # y_t = G * u_t where G is the policy coefficient
  # Var(u) = 0.01 / (1 - 0.81), Var(y) = G^2 * Var(u)
  # Gamma_y(k) = G^2 * 0.9^k * Var(u)
  sol <- make_ar1_solution(0.9, 0.1)
  mc <- model_covariance(sol)

  G_val <- sol$G[1, 1]
  var_u <- 0.01 / (1 - 0.81)
  expected_var <- G_val^2 * var_u
  expect_equal(mc$covariance[1, 1], expected_var, tolerance = 1e-4)

  # Autocovariance at lag k: Gamma_y(k) = rho^k * Gamma_y(0)
  mc2 <- model_covariance(sol, n_lags = 2)
  expect_equal(mc2$autocovariances$lag_1[1, 1],
               0.9 * mc2$covariance[1, 1], tolerance = 1e-4)
  expect_equal(mc2$autocovariances$lag_2[1, 1],
               0.81 * mc2$covariance[1, 1], tolerance = 1e-4)
})

test_that("model_covariance print method works", {
  sol <- make_ar1_solution(0.5, 0.1)
  mc <- model_covariance(sol)

  expect_output(print(mc), "Model-Implied Covariance Matrix")
  expect_output(print(mc), "Covariance")
  expect_output(print(mc), "Correlation")
})

test_that("model_covariance works on dsge_fit", {
  m <- dsge_model(
    obs(y ~ lead(y) + u),
    state(u ~ rho_u * u),
    fixed = list(),
    start = list(rho_u = 0.5)
  )

  # Simulate data
  set.seed(42)
  n <- 100
  u <- numeric(n)
  for (t in 2:n) u[t] <- 0.5 * u[t - 1] + rnorm(1, sd = 0.1)
  dat <- data.frame(y = u)

  fit <- estimate(m, data = dat)
  mc <- model_covariance(fit)

  expect_s3_class(mc, "dsge_covariance")
  expect_true(mc$covariance[1, 1] > 0)
})

test_that("model_covariance multi-variable NK model", {
  m <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.96, kappa = 0.1),
    start = list(psi = 1.5, rhou = 0.7, rhog = 0.9)
  )

  sol <- solve_dsge(m,
    params = c(beta = 0.96, kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd = c(u = 0.01, g = 0.01))

  mc <- model_covariance(sol)

  # Should have 2 observables (p, r)
  expect_equal(length(mc$variables), 2)
  expect_true(all(c("p", "r") %in% mc$variables))

  # 2x2 covariance matrix
  expect_equal(dim(mc$covariance), c(2, 2))

  # Symmetric and positive diagonal
  expect_equal(mc$covariance, t(mc$covariance), tolerance = 1e-12)
  expect_true(all(diag(mc$covariance) > 0))

  # Correlations between -1 and 1
  expect_true(all(mc$correlation >= -1 - 1e-10))
  expect_true(all(mc$correlation <= 1 + 1e-10))
})
