test_that("dsge_model creates a valid model object", {
  m <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.96),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )

  expect_s3_class(m, "dsge_model")
  expect_equal(m$variables$observed, c("p", "r"))
  expect_equal(m$variables$unobserved, "x")
  expect_equal(m$variables$exo_state, c("u", "g"))
  expect_equal(m$n_controls, 3L)
  expect_equal(m$n_states, 2L)
  expect_true("beta" %in% m$parameters)
  expect_false("beta" %in% m$free_parameters)
  expect_true("kappa" %in% m$free_parameters)
})

test_that("E() is accepted as alias for lead()", {
  m <- dsge_model(
    obs(p ~ beta * E(p) + kappa * x),
    unobs(x ~ E(x) - (r - E(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g)
  )

  expect_s3_class(m, "dsge_model")
  # Check that E() terms were parsed as leads
  p_eq <- m$equations[[1]]
  lead_terms <- Filter(function(t) t$is_lead, p_eq$rhs_terms)
  expect_true(length(lead_terms) > 0)
})

test_that("obs(), unobs(), state() require formulas", {
  expect_error(obs("not a formula"), "requires a formula")
  expect_error(unobs(42), "requires a formula")
  expect_error(state("bad"), "requires a formula")
})

test_that("dsge_model validates matching shock and observed counts", {
  expect_error(
    dsge_model(
      obs(y ~ z),
      obs(p ~ z),
      state(z ~ rho * z)
    ),
    "must equal"
  )
})

test_that("dsge_model rejects duplicate LHS variables", {
  expect_error(
    dsge_model(
      obs(y ~ z),
      obs(y ~ z),
      state(z ~ rho * z),
      state(w ~ rho2 * w)
    ),
    "more than one equation"
  )
})

test_that("dsge_model detects undeclared variables", {
  expect_error(
    dsge_model(
      obs(y ~ beta * x),
      state(z ~ rho * z)
    ),
    "not declared"
  )
})

test_that("dsge_model validates fixed parameter names", {
  expect_error(
    dsge_model(
      obs(y ~ z),
      state(z ~ rho * z),
      fixed = list(nonexistent = 1)
    ),
    "not found"
  )
})

test_that("simple AR(1) model parses correctly", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    start = list(rho = 0.8)
  )

  expect_equal(length(m$equations), 2L)
  expect_equal(m$parameters, "rho")
  expect_equal(m$free_parameters, "rho")
})

test_that("state with shock=FALSE creates endogenous state", {
  m <- dsge_model(
    obs(y ~ z),
    state(z ~ rho * z),
    state(k ~ delta * k, shock = FALSE)
  )

  expect_equal(m$variables$exo_state, "z")
  expect_equal(m$variables$endo_state, "k")
  expect_equal(m$n_exo_states, 1L)
  expect_equal(m$n_endo_states, 1L)
})
