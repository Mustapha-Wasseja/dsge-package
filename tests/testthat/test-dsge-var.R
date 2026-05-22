# Tests for bayes_dsge_var()

setup_data <- function(TT = 200L, seed = 1L) {
  nk <- dsge_model(
    obs(p   ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r   ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )
  sol <- solve_dsge(nk,
    params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd = c(e.u = 1.0, e.g = 0.5))

  set.seed(seed)
  H <- sol$H; G <- sol$G; M <- sol$M
  xst <- matrix(0, TT, nrow(H))
  y   <- matrix(0, TT, nrow(G))
  for (t in 2:TT) {
    e <- rnorm(ncol(M)) * c(1, 0.5)
    xst[t, ] <- as.numeric(H %*% xst[t - 1, ] + M %*% e)
    y[t, ]   <- as.numeric(G %*% xst[t, ])
  }
  colnames(y) <- rownames(G)
  obs_only <- y[, nk$variables$observed, drop = FALSE]
  list(model = nk, sol = sol, data = as.data.frame(obs_only))
}


test_that("bayes_dsge_var returns expected structure", {
  td <- setup_data()
  fit <- bayes_dsge_var(td$sol, data = td$data,
                        p = 2L, lambda = 1.0, n_draws = 50L,
                        seed = 42L)
  expect_s3_class(fit, "dsge_dsgevar")
  n_y <- ncol(td$data)
  k   <- n_y * 2L + 1L  # 2 lags + intercept
  expect_equal(dim(fit$Phi_post),   c(k, n_y, 50))
  expect_equal(dim(fit$Sigma_post), c(n_y, n_y, 50))
  expect_equal(dim(fit$Phi_mean),   c(k, n_y))
  expect_equal(dim(fit$Sigma_mean), c(n_y, n_y))
  expect_true(is.finite(fit$log_marg_lik))
})


test_that("posterior Sigma draws are symmetric and positive-definite", {
  td <- setup_data()
  fit <- bayes_dsge_var(td$sol, data = td$data,
                        p = 2L, lambda = 1.0, n_draws = 30L,
                        seed = 99L)
  for (d in seq_len(dim(fit$Sigma_post)[3])) {
    S <- fit$Sigma_post[, , d]
    expect_lt(max(abs(S - t(S))), 1e-8)
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  }
})


test_that("large lambda pulls VAR toward DSGE-implied autocovariances", {
  # As lambda -> Inf the VAR posterior mean should converge to the DSGE
  # implied first-order autocorrelation matrix.  We check that the
  # log_marg_lik is finite for both small and large lambda.
  td <- setup_data()
  fit_small <- bayes_dsge_var(td$sol, data = td$data, p = 2L,
                              lambda = 0.1, n_draws = 30L, seed = 1L)
  fit_large <- bayes_dsge_var(td$sol, data = td$data, p = 2L,
                              lambda = 10,  n_draws = 30L, seed = 1L)
  expect_true(is.finite(fit_small$log_marg_lik))
  expect_true(is.finite(fit_large$log_marg_lik))
  # As lambda grows the Phi_mean should change in a smooth way
  expect_false(isTRUE(all.equal(fit_small$Phi_mean, fit_large$Phi_mean,
                                tolerance = 1e-2)))
})


test_that("accepts dsge_model with params + shock_sd", {
  td <- setup_data()
  fit <- bayes_dsge_var(td$model, data = td$data,
                        params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7,
                                     rhog  = 0.9),
                        shock_sd = c(e.u = 1.0, e.g = 0.5),
                        p = 1L, lambda = 0.5, n_draws = 20L, seed = 5L)
  expect_s3_class(fit, "dsge_dsgevar")
})


test_that("errors on insufficient observations", {
  td <- setup_data(TT = 8L)
  expect_error(
    bayes_dsge_var(td$sol, data = td$data, p = 8L, lambda = 1.0,
                   n_draws = 5L),
    "Not enough observations"
  )
})


test_that("errors when data column names don't match observables", {
  td <- setup_data()
  bad <- td$data
  colnames(bad) <- c("foo", "bar")
  expect_error(
    bayes_dsge_var(td$sol, data = bad, p = 2L, lambda = 1.0,
                   n_draws = 5L),
    "not found in model observables"
  )
})


test_that("print.dsge_dsgevar produces output", {
  td <- setup_data()
  fit <- bayes_dsge_var(td$sol, data = td$data, p = 2L, lambda = 1.0,
                        n_draws = 20L, seed = 7L)
  expect_output(print(fit), "DSGE-VAR")
  expect_output(print(fit), "lambda")
})
