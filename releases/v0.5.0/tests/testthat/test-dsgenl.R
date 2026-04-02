# Tests for nonlinear DSGE functionality

# --- Parser tests ---

test_that("parse_nl_equation handles simple equation", {
  eq <- parse_nl_equation("Y = K^alpha")
  expect_equal(eq$lhs_str, "Y")
  expect_equal(eq$rhs_str, "K^alpha")
  expect_false(eq$is_state_eq)
  expect_true("K" %in% eq$current_refs)
  expect_true("alpha" %in% eq$current_refs)
  expect_true("Y" %in% eq$current_refs)
})

test_that("parse_nl_equation detects lead variables", {
  eq <- parse_nl_equation("1/C = beta/C(+1)")
  expect_true("C" %in% eq$lead_var_names)
  expect_true("C__f" %in% eq$all_names)
  expect_false(eq$is_state_eq)
})

test_that("parse_nl_equation detects state equations", {
  eq <- parse_nl_equation("K(+1) = K^alpha - C")
  expect_true(eq$is_state_eq)
  expect_equal(eq$lhs_state_var, "K")
  expect_true("K" %in% eq$lead_var_names)
})

test_that("parse_nl_equation rejects invalid equations", {
  expect_error(parse_nl_equation("K == 0"), "==")
  expect_error(parse_nl_equation("no equals here"), "must contain")
  expect_error(parse_nl_equation(42), "character string")
})

test_that("build_nl_eval_function evaluates correctly", {
  eqs <- list(
    parse_nl_equation("Y = K^alpha"),
    parse_nl_equation("K(+1) = K - C")
  )
  fn <- build_nl_eval_function(eqs)
  vals <- c(Y = 2, K = 8, alpha = 0.5, C = 1, K__f = 7)
  res <- fn(vals)
  expect_equal(res[1], 2 - 8^0.5)  # Y - K^alpha
  expect_equal(res[2], 7 - (8 - 1))  # K__f - (K - C)
})

# --- Model construction tests ---

test_that("dsgenl_model constructs correctly", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.5)
  )
  expect_s3_class(m, "dsgenl_model")
  expect_equal(m$n_controls, 1L)
  expect_equal(m$n_states, 1L)
  expect_equal(m$variables$observed, "Y")
  expect_equal(m$variables$exo_state, "Z")
  expect_true("rho" %in% m$parameters)
})

test_that("dsgenl_model validates equation counts", {
  expect_error(
    dsgenl_model(
      "Y = Z",
      observed = "Y",
      exo_state = "Z"
    ),
    "state equation"
  )
})

test_that("dsgenl_model validates lead variable declarations", {
  expect_error(
    dsgenl_model(
      "Y = W(+1)",
      "Z(+1) = rho * Z",
      observed = "Y",
      exo_state = "Z"
    ),
    "not declared"
  )
})

test_that("dsgenl_model detects parameters", {
  m <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * K^(alpha-1) + 1 - delta)",
    "K(+1) = K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C",
    endo_state = "K",
    exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9)
  )
  expect_true(all(c("alpha", "beta", "delta", "rho") %in% m$parameters))
  expect_equal(m$free_parameters, "rho")
})

# --- Steady-state tests ---

test_that("steady_state solves AR(1) model", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.8)
  )
  ss <- steady_state(m, params = c(rho = 0.8),
                     guess = c(Y = 0.1, Z = 0.1))
  expect_true(ss$converged)
  expect_equal(ss$values[["Y"]], 0, tolerance = 1e-8)
  expect_equal(ss$values[["Z"]], 0, tolerance = 1e-8)
})

test_that("steady_state solves RBC model", {
  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * K^(alpha-1) + 1 - delta)",
    "K(+1) = K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C",
    endo_state = "K",
    exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9)
  )

  params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, rho = 0.9)
  K_ss <- (0.33 / (1 / 0.99 - 1 + 0.025))^(1 / 0.67)
  C_ss <- K_ss^0.33 - 0.025 * K_ss

  ss <- steady_state(rbc, params = params,
                     guess = c(C = C_ss, K = K_ss, Z = 0))
  expect_true(ss$converged)
  expect_equal(ss$values[["Z"]], 0, tolerance = 1e-8)
  expect_equal(ss$values[["K"]], K_ss, tolerance = 1e-4)
  expect_equal(ss$values[["C"]], C_ss, tolerance = 1e-4)
})

test_that("steady_state accepts user function", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.8),
    ss_function = function(params) c(Y = 0, Z = 0)
  )
  ss <- steady_state(m, params = c(rho = 0.8))
  expect_true(ss$converged)
  expect_equal(ss$iterations, 0L)
})

# --- Linearization tests ---

test_that("linearize produces correct matrices for AR(1)", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.8)
  )
  ss <- steady_state(m, params = c(rho = 0.8),
                     guess = c(Y = 0.1, Z = 0.1))
  lin <- linearize(m, ss)

  # A0 = f_y = [1] (dY/dY = 1 from residual Y - Z)
  expect_equal(as.numeric(lin$A0), 1, tolerance = 1e-6)
  # A3 = -f_x = -(-1) = 1 (dZ coefficient, flipped)
  expect_equal(as.numeric(lin$A3), 1, tolerance = 1e-6)
  # B3 = -g_x = -(−rho) = rho = 0.8
  expect_equal(as.numeric(lin$B3), 0.8, tolerance = 1e-6)
  # B0 = g_xf = 1
  expect_equal(as.numeric(lin$B0), 1, tolerance = 1e-6)
})

# --- Solver integration tests ---

test_that("solve_dsge works with dsgenl_model (AR(1))", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.8),
    ss_guess = c(Y = 0.1, Z = 0.1)
  )
  sol <- solve_dsge(m, params = c(rho = 0.8))
  expect_true(sol$stable)
  expect_equal(as.numeric(sol$G), 1, tolerance = 1e-6)
  expect_equal(as.numeric(sol$H), 0.8, tolerance = 1e-6)
})

test_that("solve_dsge works with RBC model", {
  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * K^(alpha-1) + 1 - delta)",
    "K(+1) = K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C",
    endo_state = "K",
    exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9),
    ss_guess = c(C = 2, K = 30, Z = 0)
  )
  sol <- solve_dsge(rbc, params = c(alpha = 0.33, beta = 0.99,
                                     delta = 0.025, rho = 0.9))
  expect_true(sol$stable)
  expect_equal(nrow(sol$G), 1L)  # 1 control (C)
  expect_equal(ncol(sol$G), 2L)  # 2 states (Z, K)
  expect_true(all(Mod(sol$eigenvalues) < 1))
})

# --- IRF and forecast from nonlinear model ---

test_that("irf works after dsgenl solve", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.8),
    ss_guess = c(Y = 0.1, Z = 0.1)
  )
  sol <- solve_dsge(m, params = c(rho = 0.8))
  irfs <- irf(sol, periods = 5, se = FALSE)
  expect_s3_class(irfs, "dsge_irf")

  # IRF of Y to Z shock should be 0.8^k
  y_vals <- irfs$data[irfs$data$response == "Y", "value"]
  expected <- 0.8^(0:5)
  expect_equal(y_vals, expected, tolerance = 1e-4)
})

# --- Estimation test ---

test_that("estimate works with dsgenl_model", {
  m <- dsgenl_model(
    "Y = Z",
    "Z(+1) = rho * Z",
    observed = "Y",
    exo_state = "Z",
    start = list(rho = 0.5),
    ss_guess = c(Y = 0.1, Z = 0.1)
  )

  # Simulate data
  set.seed(42)
  n <- 150
  z <- numeric(n)
  for (i in 2:n) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  dat <- data.frame(Y = z)

  fit <- estimate(m, data = dat, control = list(maxit = 200))
  expect_s3_class(fit, "dsge_fit")
  expect_equal(fit$convergence, 0L)
  expect_true(abs(coef(fit)["rho"] - 0.8) < 0.2)
})

# --- Error handling tests ---

test_that("dsgenl_model rejects non-string equations", {
  expect_error(dsgenl_model(42, observed = "Y", exo_state = "Z"),
               "equation string")
})

test_that("dsgenl_model rejects duplicate variables", {
  expect_error(
    dsgenl_model(
      "Y = Z",
      "Z(+1) = rho * Z",
      observed = "Y",
      exo_state = "Y"
    ),
    "Duplicate"
  )
})

test_that("dsgenl_model rejects unknown fixed parameters", {
  expect_error(
    dsgenl_model(
      "Y = Z",
      "Z(+1) = rho * Z",
      observed = "Y",
      exo_state = "Z",
      fixed = list(gamma = 2)
    ),
    "not found"
  )
})
