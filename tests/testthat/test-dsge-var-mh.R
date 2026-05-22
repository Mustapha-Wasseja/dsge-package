# Tests for bayes_dsge_var_mh() -- joint MH over (theta, sigma, lambda)

setup_mh <- function(TT = 150L, seed = 1L) {
  nk <- dsge_model(
    obs(pinfobs ~ beta * lead(pinfobs) + kappa * dy + epsP),
    obs(dy ~ lead(dy) + b_r * robs + b_pi * lead(pinfobs) + epsD),
    obs(robs ~ psi_p * pinfobs + psi_y * dy + epsR),
    state(epsP ~ rho_p * epsP),
    state(epsD ~ rho_d * epsD),
    state(epsR ~ rho_r * epsR),
    fixed = list(beta = 0.99, b_r = -1.0, b_pi = 1.0),
    start = list(kappa = 0.10, psi_p = 1.5, psi_y = 0.125,
                 rho_p = 0.50, rho_d = 0.90, rho_r = 0.50)
  )
  sol <- solve_dsge(nk,
    params   = c(kappa = 0.10, psi_p = 1.5, psi_y = 0.125,
                 rho_p = 0.50, rho_d = 0.90, rho_r = 0.50),
    shock_sd = c(e.epsP = 0.15, e.epsD = 0.50, e.epsR = 0.25))

  set.seed(seed)
  H <- sol$H; G <- sol$G; M <- sol$M
  xst <- matrix(0, TT, nrow(H))
  y   <- matrix(0, TT, nrow(G))
  for (t in 2:TT) {
    e <- rnorm(ncol(M)) * c(0.15, 0.50, 0.25)
    xst[t, ] <- as.numeric(H %*% xst[t-1, ] + M %*% e)
    y[t, ]   <- as.numeric(G %*% xst[t, ])
  }
  colnames(y) <- rownames(G)
  obs_only <- y[, nk$variables$observed, drop = FALSE]
  list(model = nk, sol = sol, data = as.data.frame(obs_only))
}


test_that("bayes_dsge_var_mh returns expected structure", {
  td <- setup_mh()
  priors <- list(
    kappa = prior("beta",  shape1 = 2, shape2 = 8),
    psi_p = prior("normal", mean = 1.5, sd = 0.25),
    psi_y = prior("normal", mean = 0.125, sd = 0.05),
    rho_p = prior("beta",  shape1 = 2, shape2 = 2),
    rho_d = prior("beta",  shape1 = 8, shape2 = 2),
    rho_r = prior("beta",  shape1 = 2, shape2 = 2)
  )
  fit <- bayes_dsge_var_mh(td$model, data = td$data,
                           priors = priors,
                           p = 2L,
                           chains = 1L, iter = 200L, warmup = 100L,
                           seed = 42L)
  expect_s3_class(fit, "dsge_dsgevar_mh")
  expect_true("lambda" %in% fit$param_names)
  expect_equal(dim(fit$posterior)[3], 1)
  expect_true(is.finite(fit$lambda_posterior_mean))
  expect_true(fit$lambda_posterior_mean > 0)
})


test_that("posterior draws of lambda respect the prior support", {
  td <- setup_mh()
  priors <- list(
    kappa = prior("beta",  shape1 = 2, shape2 = 8),
    psi_p = prior("normal", mean = 1.5, sd = 0.25),
    psi_y = prior("normal", mean = 0.125, sd = 0.05),
    rho_p = prior("beta",  shape1 = 2, shape2 = 2),
    rho_d = prior("beta",  shape1 = 8, shape2 = 2),
    rho_r = prior("beta",  shape1 = 2, shape2 = 2)
  )
  fit <- bayes_dsge_var_mh(td$model, data = td$data,
                           priors = priors,
                           lambda_prior = prior("uniform",
                                                min = 0, max = 5),
                           p = 2L,
                           chains = 1L, iter = 200L, warmup = 100L,
                           seed = 42L)
  lambda_draws <- as.numeric(fit$posterior[, "lambda", ])
  expect_true(all(lambda_draws >= 0))
  expect_true(all(lambda_draws <= 5))
})


test_that("print.dsge_dsgevar_mh works", {
  td <- setup_mh()
  priors <- list(
    kappa = prior("beta",  shape1 = 2, shape2 = 8),
    psi_p = prior("normal", mean = 1.5, sd = 0.25),
    psi_y = prior("normal", mean = 0.125, sd = 0.05),
    rho_p = prior("beta",  shape1 = 2, shape2 = 2),
    rho_d = prior("beta",  shape1 = 8, shape2 = 2),
    rho_r = prior("beta",  shape1 = 2, shape2 = 2)
  )
  fit <- bayes_dsge_var_mh(td$model, data = td$data,
                           priors = priors,
                           p = 2L,
                           chains = 1L, iter = 150L, warmup = 75L,
                           seed = 7L)
  expect_output(print(fit), "DSGE-VAR with joint MH")
  expect_output(print(fit), "lambda")
})
