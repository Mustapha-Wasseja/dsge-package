# linearize() uses exact symbolic first derivatives, compiled once per model
# and cached; numerical derivatives remain the fallback.

rbc_model <- function(eq1 = "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)") {
  dsgenl_model(
    eq1,
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C", endo_state = "K", exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9), ss_guess = c(C = 2, K = 30, Z = 0)
  )
}
rbc_params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, rho = 0.9)

test_that("the symbolic Jacobian matches numerical derivatives", {
  m <- rbc_model()
  ss <- steady_state(m, params = rbc_params)
  timed <- c(m$controls, m$states, paste0(m$controls, "__f"),
             paste0(m$states, "__f"))
  point <- stats::setNames(c(ss$values[m$controls], ss$values[m$states],
                             ss$values[m$controls], ss$values[m$states]),
                           timed)
  J <- dsge:::.symbolic_jacobian(m, timed, point, rbc_params)
  num <- numDeriv::jacobian(function(z) {
    m$eval_fn(c(stats::setNames(z, timed), rbc_params))
  }, point)
  expect_equal(dim(J), c(3L, length(timed)))
  expect_equal(J, num, tolerance = 1e-8)
})

test_that("the compiled Jacobian is cached and rebuilt when equations change", {
  m <- rbc_model()
  s1 <- solve_dsge(m, params = rbc_params, shock_sd = c(Z = 0.01))
  expect_false(is.null(m$.cache$jac_code))
  code1 <- m$.cache$jac_code
  s2 <- solve_dsge(m, params = rbc_params, shock_sd = c(Z = 0.01))
  expect_identical(m$.cache$jac_code, code1)
  expect_equal(s2$G, s1$G)
  # a different model sharing the cache environment gets its own Jacobian
  m2 <- rbc_model("1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - 2 * delta)")
  m2$.cache <- m$.cache
  s3 <- solve_dsge(m2, params = rbc_params, shock_sd = c(Z = 0.01))
  expect_false(isTRUE(all.equal(s3$G, s1$G)))
  expect_equal(s3$G, solve_dsge(rbc_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - 2 * delta)"),
    params = rbc_params, shock_sd = c(Z = 0.01))$G)
})

test_that("equations stats::D cannot differentiate fall back to numerical derivatives", {
  # abs() is not in stats::D's derivatives table; C > 0 so it is harmless
  m <- rbc_model("1/abs(C) = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)")
  expect_null(dsge:::.symbolic_jacobian_code(m, c("C", "K", "Z", "C__f", "K__f", "Z__f")))
  s <- solve_dsge(m, params = rbc_params, shock_sd = c(Z = 0.01))
  ref <- solve_dsge(rbc_model(), params = rbc_params, shock_sd = c(Z = 0.01))
  expect_equal(s$G, ref$G, tolerance = 1e-6)
  expect_equal(s$H, ref$H, tolerance = 1e-6)
})

test_that("models without a cache (e.g. saved by older versions) still solve", {
  m <- rbc_model()
  m$.cache <- NULL
  s <- solve_dsge(m, params = rbc_params, shock_sd = c(Z = 0.01))
  ref <- solve_dsge(rbc_model(), params = rbc_params, shock_sd = c(Z = 0.01))
  expect_equal(s$G, ref$G)
})

test_that("the steady-state solver uses the exact Jacobian", {
  m <- rbc_model()
  n_eval <- 0L
  inner <- m$eval_fn
  m$eval_fn <- function(values) {
    n_eval <<- n_eval + 1L
    inner(values)
  }
  ss <- steady_state(m, params = rbc_params)
  # analytical steady state of the RBC model
  a <- 0.33; b <- 0.99; d <- 0.025
  K <- ((1 / b - 1 + d) / a)^(1 / (a - 1))
  expect_equal(unname(ss$values["K"]), K, tolerance = 1e-10)
  expect_equal(unname(ss$values["C"]), K^a - d * K, tolerance = 1e-10)
  # Newton with the symbolic Jacobian evaluates the equations only for the
  # residual and the line search; finite differences would need many more
  expect_lte(n_eval, 3L * ss$iterations)
})
