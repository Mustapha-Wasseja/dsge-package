# Tests for OccBin (Occasionally Binding Constraints) v0.9.0

# --- Helper: AR(1) model ---
build_ar1_sol <- function(rho = 0.9, sd = 0.01) {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = rho)
  )
  solve_dsge(m, params = c(rho = rho), shock_sd = c(z = sd))
}

# --- Helper: two-variable model with controls ---
build_2var_sol <- function() {
  m <- dsge_model(
    obs(y ~ u),
    obs(p ~ g),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    start = list(rhou = 0.7, rhog = 0.9)
  )
  solve_dsge(m, params = c(rhou = 0.7, rhog = 0.9),
             shock_sd = c(u = 0.01, g = 0.01))
}

# ===== Constraint specification tests =====

test_that("obc_constraint creates valid object", {
  con <- obc_constraint("r", ">=", 0)
  expect_s3_class(con, "obc_constraint")
  expect_equal(con$variable, "r")
  expect_equal(con$type, ">=")
  expect_equal(con$bound, 0)
  expect_null(con$shock)
})

test_that("obc_constraint validates inputs", {
  expect_error(obc_constraint(123, ">=", 0), "character")
  expect_error(obc_constraint("r", ">=", "abc"), "numeric")
  expect_error(obc_constraint("r", ">", 0))  # invalid type
})

test_that("string constraint parsing works", {
  con <- dsge:::.parse_constraint_string("r >= 0")
  expect_equal(con$variable, "r")
  expect_equal(con$type, ">=")
  expect_equal(con$bound, 0)

  con2 <- dsge:::.parse_constraint_string("b <= 0.6")
  expect_equal(con2$variable, "b")
  expect_equal(con2$type, "<=")
  expect_equal(con2$bound, 0.6)
})

test_that("invalid constraint strings give errors", {
  expect_error(dsge:::.parse_constraint_string("bad string"), "Cannot parse")
  expect_error(dsge:::.parse_constraint_string("r > 0"), "Cannot parse")
})

# ===== simulate_occbin basic tests =====

test_that("simulate_occbin returns correct object", {
  sol <- build_ar1_sol()
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.005"),
                         shocks = list(z = -0.01),
                         horizon = 20)

  expect_s3_class(obc, "dsge_occbin")
  expect_equal(obc$horizon, 20L)
  expect_equal(nrow(obc$states), 20L)
  expect_equal(nrow(obc$controls), 20L)
  expect_equal(nrow(obc$states_unc), 20L)
  expect_true(is.logical(obc$binding))
  expect_true(obc$converged)
})

test_that("constraint never binding gives same path as unconstrained", {
  sol <- build_ar1_sol()

  # Positive shock with lower bound — should never bind
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -1"),
                         shocks = list(z = 0.01),
                         horizon = 20)

  expect_false(any(obc$binding))
  expect_equal(obc$states, obc$states_unc, tolerance = 1e-10)
  expect_equal(obc$controls, obc$controls_unc, tolerance = 1e-10)
})

test_that("constraint binds when violated", {
  sol <- build_ar1_sol(0.9, 1.0)

  # Large negative shock pushes y below -0.5
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.5"),
                         shocks = list(z = -1.0),
                         horizon = 20)

  # Should bind in early periods
  expect_true(any(obc$binding))

  # Constrained y should never go below bound
  y_con <- obc$controls[, 1]
  expect_true(all(y_con >= -0.5 - 1e-6))

  # Unconstrained y should violate bound
  y_unc <- obc$controls_unc[, 1]
  expect_true(any(y_unc < -0.5))
})

test_that("constraint eventually stops binding", {
  sol <- build_ar1_sol(0.9, 1.0)

  # One-time shock, constraint should bind early then release
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.3"),
                         shocks = list(z = -1.0),
                         horizon = 30)

  # Should bind in early periods
  bind <- obc$binding[, 1]
  expect_true(any(bind))

  # Should stop binding later (as shock decays)
  last_bind <- max(which(bind))
  expect_lt(last_bind, 30L)

  # After constraint releases, paths should converge
  late <- 25:30
  expect_true(all(abs(obc$controls[late, 1]) < 0.3))
})

test_that("upper bound constraint works", {
  sol <- build_ar1_sol(0.9, 1.0)

  obc <- simulate_occbin(sol,
                         constraints = list("y <= 0.5"),
                         shocks = list(z = 1.0),
                         horizon = 20)

  # Should bind
  expect_true(any(obc$binding))

  # y should not exceed 0.5
  y_con <- obc$controls[, 1]
  expect_true(all(y_con <= 0.5 + 1e-6))
})

test_that("multiple constraints work", {
  sol <- build_2var_sol()

  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.005", "p >= -0.005"),
                         shocks = list(u = -0.01, g = -0.01),
                         horizon = 20)

  expect_equal(ncol(obc$binding), 2L)
  expect_true(obc$converged)
})

test_that("constraint with obc_constraint object works", {
  sol <- build_ar1_sol(0.9, 1.0)

  con <- obc_constraint("y", ">=", -0.5)
  obc <- simulate_occbin(sol,
                         constraints = list(con),
                         shocks = list(z = -1.0),
                         horizon = 20)

  expect_true(any(obc$binding))
  expect_true(all(obc$controls[, 1] >= -0.5 - 1e-6))
})

test_that("initial conditions work with OBC", {
  sol <- build_ar1_sol(0.9, 1.0)

  # Displaced initial condition
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.3"),
                         initial = c(z = -1.0),
                         horizon = 20)

  expect_true(any(obc$binding))
  expect_true(obc$converged)
})

test_that("print method works", {
  sol <- build_ar1_sol(0.9, 1.0)
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.5"),
                         shocks = list(z = -1.0),
                         horizon = 20)

  expect_output(print(obc), "OccBin")
  expect_output(print(obc), "binds")
})

test_that("summary method works", {
  sol <- build_ar1_sol(0.9, 1.0)
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.5"),
                         shocks = list(z = -1.0),
                         horizon = 20)

  expect_output(summary(obc), "Constraint")
  expect_output(summary(obc), "Binding")
})

test_that("plot method works without error", {
  sol <- build_ar1_sol(0.9, 1.0)
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.5"),
                         shocks = list(z = -1.0),
                         horizon = 20)

  expect_silent(plot(obc))
})

test_that("plot with selected variables works", {
  sol <- build_2var_sol()
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.005"),
                         shocks = list(u = -0.01),
                         horizon = 20)

  expect_silent(plot(obc, vars = "y"))
})

test_that("error on invalid constraint variable", {
  sol <- build_ar1_sol()
  expect_error(
    simulate_occbin(sol, constraints = list("bad_var >= 0")),
    "not found"
  )
})

test_that("error on zero-impact shadow shock", {
  sol <- build_2var_sol()
  # p is driven by g shock only; u shock has zero impact on p
  # Specifying u as shadow shock for p should error
  con <- obc_constraint("p", ">=", -0.005, shock = "u")
  expect_error(
    simulate_occbin(sol, constraints = list(con), shocks = list(g = -0.01)),
    "zero impact"
  )
})

test_that("nonlinear model works with OBC", {
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

  # Constraint on H (which has nonzero impact from A shock)
  # H_ss = 1, constrain H >= 0.99 in levels → bound_dev = 0.99 - 1 = -0.01
  obc <- simulate_occbin(sol,
                         constraints = list("H >= 0.99"),
                         shocks = list(A = -0.02),
                         horizon = 20)

  expect_s3_class(obc, "dsge_occbin")
  expect_true(obc$converged)
})

test_that("convergence iteration count is reasonable", {
  sol <- build_ar1_sol(0.9, 1.0)
  obc <- simulate_occbin(sol,
                         constraints = list("y >= -0.5"),
                         shocks = list(z = -1.0),
                         horizon = 20)

  # Should converge quickly (typically 2-5 iterations)
  expect_lt(obc$n_iter, 20L)
})

test_that("in_sd mode works with OBC", {
  sol <- build_ar1_sol(0.9, 0.01)

  obc_sd <- simulate_occbin(sol,
                            constraints = list("y >= -0.005"),
                            shocks = list(z = -1),
                            horizon = 20, in_sd = TRUE)
  obc_lev <- simulate_occbin(sol,
                             constraints = list("y >= -0.005"),
                             shocks = list(z = -0.01),
                             horizon = 20, in_sd = FALSE)

  expect_equal(obc_sd$states_unc, obc_lev$states_unc, tolerance = 1e-10)
})
