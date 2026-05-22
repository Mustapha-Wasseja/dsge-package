# Tests for conditional_forecast()

setup_fit <- function() {
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

  set.seed(42)
  H <- sol$H; G <- sol$G; M <- sol$M
  TT <- 100
  xst <- matrix(0, TT, nrow(H))
  y   <- matrix(0, TT, nrow(G))
  for (t in 2:TT) {
    e <- rnorm(ncol(M)) * c(1, 0.5)
    xst[t, ] <- as.numeric(H %*% xst[t - 1, ] + M %*% e)
    y[t, ]   <- as.numeric(G %*% xst[t, ])
  }
  colnames(y) <- rownames(G)
  dat <- as.data.frame(y)
  estimate(nk, data = dat)
}


test_that("conditional_forecast returns expected structure", {
  fit <- setup_fit()
  cf <- conditional_forecast(fit, horizon = 12,
    condition = list(r = c(0.5, 0.5, 0.5, 0.5, rep(NA, 8))))

  expect_s3_class(cf, "dsge_conditional_forecast")
  expect_s3_class(cf, "dsge_forecast")
  expect_equal(cf$horizon, 12L)
  expect_true("conditioned" %in% names(cf$forecasts))
})


test_that("conditioned variable hits the target path", {
  fit <- setup_fit()
  target_path <- c(0.5, 0.5, 0.5, 0.5)
  cf <- conditional_forecast(fit, horizon = 12,
    condition = list(r = c(target_path, rep(NA, 8))))

  # Pull conditioned-period forecast values for 'r'
  fc_r <- cf$forecasts[cf$forecasts$variable == "r" &
                       cf$forecasts$period <= 4, "value"]
  expect_equal(unname(fc_r), target_path, tolerance = 1e-6)
})


test_that("unconditioned variables move endogenously to support the path", {
  fit <- setup_fit()
  cf <- conditional_forecast(fit, horizon = 8,
    condition = list(r = c(1, 1, 1, rep(NA, 5))))

  # Implied shocks should be non-trivial when conditioning
  expect_true(max(abs(cf$shocks)) > 1e-6)
  # Inflation forecasts should respond (non-zero)
  fc_p <- cf$forecasts[cf$forecasts$variable == "p", "value"]
  expect_true(max(abs(fc_p)) > 1e-6)
})


test_that("all-NA condition reduces to ordinary forecast", {
  fit <- setup_fit()
  cf <- conditional_forecast(fit, horizon = 6,
    condition = list(r = rep(NA, 6)))
  fc <- forecast(fit, horizon = 6)
  # The function should fall back to standard forecast in this case
  expect_equal(cf$forecasts$value, fc$forecasts$value, tolerance = 1e-8)
})


test_that("unknown observable in condition errors", {
  fit <- setup_fit()
  expect_error(
    conditional_forecast(fit, horizon = 6,
                         condition = list(not_an_obs = c(1, 1, 1, NA, NA, NA))),
    "Unknown observable"
  )
})


test_that("condition vector longer than horizon errors", {
  fit <- setup_fit()
  expect_error(
    conditional_forecast(fit, horizon = 4,
                         condition = list(r = c(1, 1, 1, 1, 1))),
    "longer than horizon"
  )
})


test_that("conditioned column has correct TRUE/FALSE pattern", {
  fit <- setup_fit()
  cf <- conditional_forecast(fit, horizon = 8,
    condition = list(r = c(0.5, NA, 0.5, NA, NA, NA, NA, NA)))
  cond_r <- cf$forecasts$conditioned[cf$forecasts$variable == "r"]
  expect_equal(cond_r, c(TRUE, FALSE, TRUE, FALSE,
                         FALSE, FALSE, FALSE, FALSE))
  cond_p <- cf$forecasts$conditioned[cf$forecasts$variable == "p"]
  expect_true(all(!cond_p))
})


test_that("print.dsge_conditional_forecast works", {
  fit <- setup_fit()
  cf <- conditional_forecast(fit, horizon = 8,
    condition = list(r = c(0.5, 0.5, 0.5, rep(NA, 5))))
  expect_output(print(cf), "Conditional DSGE Forecast")
})


test_that("plot.dsge_forecast accepts a conditional forecast object", {
  fit <- setup_fit()
  cf <- conditional_forecast(fit, horizon = 8,
    condition = list(r = c(0.5, 0.5, 0.5, rep(NA, 5))))
  expect_silent({
    pdf(tempfile())
    plot(cf)
    dev.off()
  })
})
