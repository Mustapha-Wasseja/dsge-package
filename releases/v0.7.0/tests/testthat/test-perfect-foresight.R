# Tests for perfect_foresight()
# Feature 3: Deterministic Transition Paths (v0.7.0)

# --- Helper: build a simple AR(1) model and solve ---
build_ar1 <- function(rho = 0.9, sd = 0.01) {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = rho)
  )
  solve_dsge(m, params = c(rho = rho), shock_sd = c(z = sd))
}

test_that("perfect_foresight works with AR(1) model", {
  sol <- build_ar1(0.9, 0.01)
  pf <- perfect_foresight(sol, shocks = list(z = 0.01), horizon = 40)

  expect_s3_class(pf, "dsge_perfect_foresight")
  expect_equal(pf$horizon, 40L)
  expect_equal(nrow(pf$states), 40L)
  expect_equal(ncol(pf$states), 1L)
  expect_equal(ncol(pf$controls), 1L)
})

test_that("AR(1) impulse response matches analytical solution", {
  rho <- 0.8
  sol <- build_ar1(rho, 1.0)

  pf <- perfect_foresight(sol, shocks = list(z = 1.0), horizon = 20)

  # Analytical: z_1 = M*1 = 1, z_t = rho^(t-1) for t=1,2,...
  for (t in 1:20) {
    expect_equal(unname(pf$states[t, 1]), rho^(t - 1), tolerance = 1e-10,
                 label = paste("state at t =", t))
  }
})

test_that("initial conditions produce correct convergence", {
  rho <- 0.9
  sol <- build_ar1(rho, 0.01)

  pf <- perfect_foresight(sol, initial = c(z = 1.0), horizon = 60)

  # z_t = rho^t * z_0 (initial condition decays)
  for (t in 1:60) {
    expect_equal(unname(pf$states[t, 1]), rho^t, tolerance = 1e-10,
                 label = paste("IC at t =", t))
  }

  # Should converge to zero
  expect_lt(abs(pf$states[60, 1]), 0.01)
})

test_that("temporary multi-period shock works", {
  rho <- 0.5
  sol <- build_ar1(rho, 1.0)

  pf <- perfect_foresight(sol, shocks = list(z = c(1, 1, 1)), horizon = 20)

  # Period 1: rho*0 + 1 = 1
  expect_equal(unname(pf$states[1, 1]), 1.0, tolerance = 1e-10)
  # Period 2: rho*1 + 1 = 1.5
  expect_equal(unname(pf$states[2, 1]), 1.5, tolerance = 1e-10)
  # Period 3: rho*1.5 + 1 = 1.75
  expect_equal(unname(pf$states[3, 1]), 1.75, tolerance = 1e-10)
  # Period 4: rho*1.75 + 0 = 0.875
  expect_equal(unname(pf$states[4, 1]), 0.875, tolerance = 1e-10)
})

test_that("permanent shock produces new steady state", {
  rho <- 0.5
  sol <- build_ar1(rho, 1.0)

  perm_shock <- rep(1.0, 100)
  pf <- perfect_foresight(sol, shocks = list(z = perm_shock), horizon = 100)

  # New SS: x = rho*x + 1 => x = 1/(1-rho) = 2
  expect_equal(unname(pf$states[100, 1]), 1 / (1 - rho), tolerance = 0.001)
})

test_that("controls track states via G matrix", {
  sol <- build_ar1(0.9, 0.01)
  pf <- perfect_foresight(sol, shocks = list(z = 0.01), horizon = 20)

  for (t in 1:20) {
    expected <- as.numeric(sol$G %*% pf$states[t, ])
    expect_equal(unname(pf$controls[t, 1]), expected, tolerance = 1e-10)
  }
})

test_that("no shocks and no initial conditions gives zero paths", {
  sol <- build_ar1(0.9, 0.01)
  pf <- perfect_foresight(sol, horizon = 10)

  expect_true(all(pf$states == 0))
  expect_true(all(pf$controls == 0))
})

test_that("in_sd mode scales shocks by standard deviation", {
  sol <- build_ar1(0.9, 0.02)

  pf_sd <- perfect_foresight(sol, shocks = list(z = 1.0), in_sd = TRUE, horizon = 5)
  pf_lev <- perfect_foresight(sol, shocks = list(z = 0.02), in_sd = FALSE, horizon = 5)

  expect_equal(pf_sd$states, pf_lev$states, tolerance = 1e-12)
})

test_that("shock path matrix input works", {
  sol <- build_ar1(0.9, 0.01)

  shock_mat <- matrix(c(0.01, 0.005, 0, 0, 0), ncol = 1)
  colnames(shock_mat) <- "z"
  pf_mat <- perfect_foresight(sol, shocks = shock_mat, horizon = 5)

  pf_list <- perfect_foresight(sol, shocks = list(z = c(0.01, 0.005)), horizon = 5)

  expect_equal(pf_mat$states, pf_list$states, tolerance = 1e-12)
})

test_that("print method works", {
  sol <- build_ar1(0.9, 0.01)
  pf <- perfect_foresight(sol, shocks = list(z = 0.01), horizon = 20)
  expect_output(print(pf), "Perfect Foresight")
})

test_that("summary method works", {
  sol <- build_ar1(0.9, 0.01)
  pf <- perfect_foresight(sol, shocks = list(z = 0.01), horizon = 40)
  expect_output(summary(pf), "Impact effects")
})

test_that("plot method runs without error", {
  sol <- build_ar1(0.9, 0.01)
  pf <- perfect_foresight(sol, shocks = list(z = 0.01), horizon = 20)
  expect_silent(plot(pf))
})

test_that("error on unknown shock name", {
  sol <- build_ar1(0.9, 0.01)
  expect_error(perfect_foresight(sol, shocks = list(bad = 0.01)), "Unknown shock name")
})

test_that("error on unknown state name in initial", {
  sol <- build_ar1(0.9, 0.01)
  expect_error(perfect_foresight(sol, initial = c(bad = 1.0)), "Unknown state name")
})

test_that("horizon validation works", {
  sol <- build_ar1(0.9, 0.01)
  expect_error(perfect_foresight(sol, horizon = 0), "at least 1")
})

test_that("combined initial conditions and shocks work", {
  rho <- 0.8
  sol <- build_ar1(rho, 1.0)

  pf <- perfect_foresight(sol, initial = c(z = 0.5),
                          shocks = list(z = 1.0), horizon = 10)

  # Period 1: H * 0.5 + M * 1.0 = 0.8*0.5 + 1.0 = 1.4
  expect_equal(unname(pf$states[1, 1]), rho * 0.5 + 1.0, tolerance = 1e-10)
  # Period 2: H * 1.4 + 0 = 1.12
  expect_equal(unname(pf$states[2, 1]), rho * 1.4, tolerance = 1e-10)
})

test_that("two-state model with multiple controls works", {
  # Two AR(1) states, two observables
  m <- dsge_model(
    obs(y ~ u),
    obs(p ~ g),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    start = list(rhou = 0.7, rhog = 0.9)
  )
  sol <- solve_dsge(m, params = c(rhou = 0.7, rhog = 0.9),
                    shock_sd = c(u = 0.01, g = 0.01))

  # Shock only u
  pf <- perfect_foresight(sol, shocks = list(u = 0.01), horizon = 30)

  expect_equal(nrow(pf$states), 30L)
  expect_equal(ncol(pf$states), 2L)
  expect_equal(ncol(pf$controls), 2L)

  # u shock should not affect g state
  expect_true(all(abs(pf$states[, "g"]) < 1e-12))
  # u state should decay
  expect_lt(abs(pf$states[30, "u"]), 0.001)
})

test_that("nonlinear model works", {
  mod <- dsgenl_model(
    "Y = A * H",
    "1 / C = beta * R / C(+1)",
    "Y = C",
    "R = 1 / beta",
    "A(+1) = A^rhoA",
    observed = "Y",
    unobserved = c("C", "H", "R"),
    exo_state = "A",
    fixed = list(beta = 0.99),
    start = list(rhoA = 0.90),
    ss_guess = c(Y = 1, C = 1, H = 1, R = 1.0101, A = 1)
  )
  sol <- solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
                    shock_sd = c(A = 0.01))

  pf <- perfect_foresight(sol, shocks = list(A = 0.01), horizon = 30)

  expect_s3_class(pf, "dsge_perfect_foresight")
  expect_false(is.null(pf$state_levels))
  expect_false(is.null(pf$control_levels))
  expect_lt(max(abs(pf$states[30, ])), 0.001)
})

test_that("levels output adds steady state correctly", {
  mod <- dsgenl_model(
    "Y = A * H",
    "1 / C = beta * R / C(+1)",
    "Y = C",
    "R = 1 / beta",
    "A(+1) = A^rhoA",
    observed = "Y",
    unobserved = c("C", "H", "R"),
    exo_state = "A",
    fixed = list(beta = 0.99),
    start = list(rhoA = 0.90),
    ss_guess = c(Y = 1, C = 1, H = 1, R = 1.0101, A = 1)
  )
  sol <- solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
                    shock_sd = c(A = 0.01))

  pf <- perfect_foresight(sol, shocks = list(A = 0.01), horizon = 10)

  ss <- sol$steady_state
  for (nm in pf$state_names) {
    if (nm %in% names(ss)) {
      expect_equal(pf$state_levels[, nm],
                   pf$states[, nm] + ss[nm],
                   tolerance = 1e-12)
    }
  }
})

test_that("plot with type='level' works for nonlinear model", {
  mod <- dsgenl_model(
    "Y = A * H",
    "1 / C = beta * R / C(+1)",
    "Y = C",
    "R = 1 / beta",
    "A(+1) = A^rhoA",
    observed = "Y",
    unobserved = c("C", "H", "R"),
    exo_state = "A",
    fixed = list(beta = 0.99),
    start = list(rhoA = 0.90),
    ss_guess = c(Y = 1, C = 1, H = 1, R = 1.0101, A = 1)
  )
  sol <- solve_dsge(mod, params = c(beta = 0.99, rhoA = 0.90),
                    shock_sd = c(A = 0.01))

  pf <- perfect_foresight(sol, shocks = list(A = 0.01), horizon = 20)
  expect_silent(plot(pf, type = "level"))
})
