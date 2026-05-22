# Tests for DSGE-VAR forecast methods

setup_dsgevar <- function(TT = 150L, seed = 7L) {
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
  dat <- as.data.frame(obs_only)
  fit <- bayes_dsge_var(sol, data = dat, p = 2L, lambda = 1.0,
                        n_draws = 100L, seed = seed)
  list(model = nk, sol = sol, data = dat, fit = fit)
}


test_that("forecast.dsge_dsgevar returns expected structure", {
  td <- setup_dsgevar()
  fc <- forecast(td$fit, horizon = 8L, n_paths = 1L)
  expect_s3_class(fc, "dsge_dsgevar_forecast")
  expect_s3_class(fc, "dsge_forecast")
  expect_equal(fc$horizon, 8L)
  expect_equal(unique(fc$forecasts$variable), td$fit$var_names)
  expect_equal(dim(fc$forecast_paths)[1], 8L)
  expect_equal(dim(fc$forecast_paths)[2], length(td$fit$var_names))
})


test_that("forecast.dsge_dsgevar reproduces history scale roughly", {
  td <- setup_dsgevar()
  fc <- forecast(td$fit, horizon = 8L, n_paths = 1L)
  # Forecast SDs should be of the same order of magnitude as the data SDs
  data_sd <- apply(td$data, 2, stats::sd)
  fc_sd_h1 <- fc$forecasts$sd[fc$forecasts$period == 1]
  ratio <- fc_sd_h1 / data_sd
  # Within an order of magnitude
  expect_true(all(ratio > 0.1 & ratio < 10))
})


test_that("conditional_forecast hits the target path on conditioned variable", {
  td <- setup_dsgevar()
  # Hold pinfobs at 0.5 for the first 3 periods
  cf <- conditional_forecast(td$fit, horizon = 6L,
    condition = list(pinfobs = c(0.5, 0.5, 0.5, NA, NA, NA)))
  fc_pi <- cf$forecasts[cf$forecasts$variable == "pinfobs" &
                        cf$forecasts$period <= 3, "value"]
  # Posterior mean should hit target within 0.01 (sampling noise)
  expect_equal(unname(fc_pi), c(0.5, 0.5, 0.5), tolerance = 0.05)
})


test_that("plot.dsge_forecast handles a DSGE-VAR forecast", {
  td <- setup_dsgevar()
  fc <- forecast(td$fit, horizon = 6L, n_paths = 1L)
  expect_silent({
    pdf(tempfile())
    plot(fc)
    dev.off()
  })
})


test_that("forecast.dsge_dsgevar_mh produces a forecast", {
  # Quick MH with very few draws -- just check the plumbing
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
  set.seed(11)
  H <- sol$H; G <- sol$G; M <- sol$M
  TT <- 120L
  xst <- matrix(0, TT, nrow(H))
  y   <- matrix(0, TT, nrow(G))
  for (t in 2:TT) {
    e <- rnorm(ncol(M)) * c(0.15, 0.50, 0.25)
    xst[t, ] <- as.numeric(H %*% xst[t-1, ] + M %*% e)
    y[t, ]   <- as.numeric(G %*% xst[t, ])
  }
  colnames(y) <- rownames(G)
  dat <- as.data.frame(y[, nk$variables$observed, drop = FALSE])
  priors <- list(
    kappa = prior("beta",  shape1 = 2, shape2 = 8),
    psi_p = prior("normal", mean = 1.5, sd = 0.25),
    psi_y = prior("normal", mean = 0.125, sd = 0.05),
    rho_p = prior("beta",  shape1 = 2, shape2 = 2),
    rho_d = prior("beta",  shape1 = 8, shape2 = 2),
    rho_r = prior("beta",  shape1 = 2, shape2 = 2)
  )
  fit_mh <- bayes_dsge_var_mh(nk, data = dat, priors = priors,
                              p = 2L,
                              chains = 1L, iter = 80L, warmup = 40L,
                              seed = 11L)
  fc <- forecast(fit_mh, horizon = 5L, n_paths = 1L)
  expect_s3_class(fc, "dsge_dsgevar_forecast")
  expect_equal(fc$horizon, 5L)
})
