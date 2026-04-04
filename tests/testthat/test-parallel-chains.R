# tests/testthat/test-parallel-chains.R

make_simple_model_and_data <- function(seed = 99L) {
  m <- dsge_model(
    obs(y ~ lead(y) + u),
    state(u ~ rho_u * u),
    start = list(rho_u = 0.5)
  )
  set.seed(seed)
  n <- 60
  u <- numeric(n)
  for (t in 2:n) u[t] <- 0.5 * u[t - 1] + rnorm(1, sd = 0.1)
  list(model = m, data = data.frame(y = u),
       priors = list(rho_u = prior("beta", shape1 = 5, shape2 = 5)))
}

test_that("n_cores = 1 (sequential) still works", {
  skip_on_cran()
  d <- make_simple_model_and_data()
  fit <- bayes_dsge(d$model, d$data, d$priors,
                    chains = 2L, iter = 200L, warmup = 100L,
                    seed = 1L, n_cores = 1L)
  expect_s3_class(fit, "dsge_bayes")
  expect_equal(fit$n_chains, 2L)
})

test_that("n_cores > chains is silently clamped", {
  skip_on_cran()
  d <- make_simple_model_and_data()
  # n_cores = 10 but chains = 2; should work without error
  fit <- bayes_dsge(d$model, d$data, d$priors,
                    chains = 2L, iter = 200L, warmup = 100L,
                    seed = 2L, n_cores = 10L)
  expect_s3_class(fit, "dsge_bayes")
})

test_that("sequential and parallel give same posterior dimensions", {
  skip_on_cran()
  # On this platform parallel may fall back to sequential (Windows without
  # installed package, or single-core CI).  Either way the dimensions must
  # match those of the sequential run.
  d <- make_simple_model_and_data()
  fit_seq <- bayes_dsge(d$model, d$data, d$priors,
                        chains = 2L, iter = 200L, warmup = 100L,
                        seed = 3L, n_cores = 1L)
  fit_par <- bayes_dsge(d$model, d$data, d$priors,
                        chains = 2L, iter = 200L, warmup = 100L,
                        seed = 3L, n_cores = 2L)

  expect_equal(dim(fit_seq$posterior), dim(fit_par$posterior))
})

test_that(".run_chain_worker returns required fields", {
  skip_on_cran()
  d <- make_simple_model_and_data()
  m <- d$model
  pr <- list(rho_u = prior("beta", shape1 = 5, shape2 = 5))
  pr_full <- dsge:::validate_priors(pr, "rho_u", m$variables$exo_state)

  set.seed(5L)
  y  <- as.matrix(d$data)
  y  <- sweep(y, 2, colMeans(y))

  # Build a minimal chain_args list
  args <- list(
    model          = m,
    y              = y,
    prior_list     = pr_full,
    free_params    = "rho_u",
    shock_names    = m$variables$exo_state,
    all_fixed      = m$fixed,
    is_nonlinear   = FALSE,
    obs_vars       = m$variables$observed,
    start_u        = c(0, 0),    # unconstrained: logit(0.5)=0 for rho, log(1)=0 for sd
    iter           = 150L,
    warmup         = 75L,
    thin           = 1L,
    proposal_scale = 0.1,
    mode_hessian   = NULL,
    chain_seed     = 42L
  )

  res <- dsge:::.run_chain_worker(args)
  expect_true(is.matrix(res$draws))
  expect_true(is.numeric(res$acceptance_rate))
  expect_true(res$acceptance_rate >= 0 && res$acceptance_rate <= 1)
  expect_true(is.integer(res$solve_failures))
})
