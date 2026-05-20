# tests/testthat/test-perfect-foresight-nonlinear.R

# ---- Shared model fixture -----------------------------------------------
# Simple nonlinear model: Euler equation + capital accumulation + AR(1) TFP
make_rbc <- function() {
  dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed   = "C",
    endo_state = "K",
    exo_state  = "Z",
    fixed  = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start  = list(rho = 0.9)
  )
}

# ---- Tests ------------------------------------------------------------------

test_that("returns dsge_perfect_foresight class", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           horizon  = 20L)

  expect_s3_class(pf, "dsge_perfect_foresight")
})

test_that("output dimensions are correct", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           horizon  = 20L)

  expect_equal(nrow(pf$states),   20L)
  expect_equal(ncol(pf$states),   2L)   # K, Z
  expect_equal(nrow(pf$controls), 20L)
  expect_equal(ncol(pf$controls), 1L)   # C
  expect_equal(pf$horizon, 20L)
})

test_that("nonlinear flag and newton_iters are present", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           horizon  = 20L)

  expect_true(isTRUE(pf$nonlinear))
  expect_true(is.integer(pf$newton_iters))
  expect_true(pf$newton_iters >= 0L)
})

test_that("converges for a small shock", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           shocks   = list(Z = 0.01),
           horizon  = 30L,
           tol      = 1e-7)

  expect_true(pf$converged)
  expect_true(pf$max_residual < 1e-7)
})

test_that("zero shock leaves all variables at steady state", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           shocks   = NULL,
           initial  = NULL,
           horizon  = 20L,
           tol      = 1e-8)

  # Deviations from SS should all be essentially zero
  expect_true(max(abs(pf$states))   < 1e-6)
  expect_true(max(abs(pf$controls)) < 1e-6)
})

test_that("large shock produces non-trivial response", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           shocks   = list(Z = 0.10),
           horizon  = 40L)

  # Some variable should respond
  expect_true(max(abs(pf$controls)) > 1e-4)
  expect_true(max(abs(pf$states))   > 1e-4)
})

test_that("nonlinear path differs from linear for a large shock", {
  mod <- make_rbc()
  p   <- c(rho = 0.9)
  sd  <- c(Z = 0.01)
  sh  <- list(Z = 0.20)    # large shock: nonlinear effects should be visible

  pf_nl <- perfect_foresight_nonlinear(mod, params = p, shock_sd = sd,
             shocks = sh, horizon = 30L, tol = 1e-7)

  sol1  <- solve_dsge(mod, params = c(rho = 0.9, alpha = 0.33,
                                      beta = 0.99, delta = 0.025),
                      shock_sd = sd, order = 1L)
  pf_lin <- perfect_foresight(sol1, shocks = sh, horizon = 30L)

  # Paths should exist for both
  expect_equal(dim(pf_nl$controls), dim(pf_lin$controls))

  # For a sizeable shock, the nonlinear path should deviate from the linear one
  max_diff <- max(abs(pf_nl$controls - pf_lin$controls))
  expect_true(max_diff > 0)   # strictly different (even numerically)
})

test_that("initial condition displaces path from SS", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           initial  = c(K = 0.1),   # capital 10% above SS
           horizon  = 30L)

  # K should start above SS and return toward it
  expect_true(pf$states[1L, "K"] > 0)
})

test_that("levels are consistent with deviations plus steady state", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           shocks   = list(Z = 0.05),
           horizon  = 20L)

  ss_C <- pf$steady_state["C"]
  ss_K <- pf$steady_state["K"]

  expect_equal(pf$control_levels[, "C"],
               pf$controls[, "C"] + ss_C,
               tolerance = 1e-10)
  expect_equal(pf$state_levels[, "K"],
               pf$states[, "K"] + ss_K,
               tolerance = 1e-10)
})

test_that("print method works (inherited)", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           horizon  = 10L)

  expect_output(print(pf), "Perfect Foresight")
})

test_that("summary method works (inherited)", {
  mod <- make_rbc()
  pf  <- perfect_foresight_nonlinear(mod,
           params   = c(rho = 0.9),
           shock_sd = c(Z = 0.01),
           shocks   = list(Z = 0.05),
           horizon  = 10L)

  expect_output(summary(pf), "Impact")
})

test_that("error on non-dsgenl_model input", {
  lin_mod <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.9)
  )
  expect_error(
    perfect_foresight_nonlinear(lin_mod, params = c(rho = 0.9),
                                 shock_sd = c(z = 0.01)),
    "dsgenl_model"
  )
})

test_that("unknown shock name in shocks list errors", {
  mod <- make_rbc()
  expect_error(
    perfect_foresight_nonlinear(mod, params = c(rho = 0.9),
                                 shock_sd = c(Z = 0.01),
                                 shocks   = list(BadName = 0.1)),
    "Unknown exogenous state"
  )
})

test_that("unknown variable in initial errors", {
  mod <- make_rbc()
  expect_error(
    perfect_foresight_nonlinear(mod, params = c(rho = 0.9),
                                 shock_sd = c(Z = 0.01),
                                 initial  = c(BadState = 0.1)),
    "Unknown state"
  )
})
