# Tests for endogenous_prior()

make_model_and_sol <- function() {
  m <- dsge_model(
    obs(y ~ rho * y_lag + e),
    state(y_lag ~ y, shock = FALSE),
    state(e ~ phi * e),
    fixed = list(rho = 0.5, phi = 0.5))
  sol <- solve_dsge(m, params = c(rho = 0.5, phi = 0.5), shock_sd = c(e = 1))
  list(model = m, sol = sol)
}

test_that("endogenous_prior returns the expected structure", {
  ep <- endogenous_prior(c(sd_y = 1.2, ac1_y = 0.5))
  expect_s3_class(ep, "dsge_endog_prior")
  expect_true(is.function(ep$moments_fn))
  expect_true(is.function(ep$log_density))
})

test_that("rejects bad moment names", {
  expect_error(endogenous_prior(c(foo_y = 1)),
               "must start with 'sd_', 'ac1_', or 'cor_'")
  expect_error(endogenous_prior(c(1.2)),
               "fully-named")
})

test_that("moment vector at the solution matches expected dimension", {
  ms <- make_model_and_sol()
  ep <- endogenous_prior(c(sd_y = 1.2, ac1_y = 0.5))
  m  <- ep$moments_fn(ms$sol)
  expect_named(m, c("sd_y","ac1_y"))
  expect_true(all(is.finite(m)))
})

test_that("log density is maximised when moments equal target", {
  ms <- make_model_and_sol()
  ep <- endogenous_prior(c(sd_y = 1.2, ac1_y = 0.5))
  m_at_sol <- ep$moments_fn(ms$sol)
  # Target equal to model-implied -> log density 0; offset target gives <0
  ep_match <- endogenous_prior(m_at_sol)
  lp_match  <- ep_match$log_density(ms$sol)
  lp_off    <- ep$log_density(ms$sol)
  expect_equal(lp_match, 0, tolerance = 1e-10)
  expect_lt(lp_off, lp_match)
})

test_that("bayes_dsge runs with endogenous_prior plumbed in", {
  m <- dsge_model(
    obs(y ~ rho * y_lag + e),
    state(y_lag ~ y, shock = FALSE),
    state(e ~ phi * e),
    start = list(rho = 0.6, phi = 0.5))
  set.seed(1)
  ee <- rnorm(80)
  z <- numeric(80); for (i in 2:80) z[i] <- 0.7 * z[i-1] + ee[i]

  emp_sd  <- stats::sd(z)
  emp_ac1 <- cor(z[-1], z[-length(z)])
  # Use a loose weight (50% SE) so the test is robust
  ep <- endogenous_prior(c(sd_y = emp_sd, ac1_y = emp_ac1),
                        weight = diag(c(1/(emp_sd*0.5)^2,
                                         1/(emp_ac1*0.5)^2)))
  fit <- bayes_dsge(m, data = data.frame(y = z),
                    priors = list(
                      rho = prior("beta",   shape1 = 2, shape2 = 2),
                      phi = prior("beta",   shape1 = 2, shape2 = 2)),
                    chains = 1, iter = 200, warmup = 100,
                    seed = 1, endogenous_prior = ep)
  expect_s3_class(fit, "dsge_bayes")
})
