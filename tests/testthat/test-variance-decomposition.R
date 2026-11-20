# Tests for variance_decomposition() (unconditional + FEVD)

make_test_model <- function() {
  m <- dsge_model(
    obs(p   ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r   ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )
  sol <- solve_dsge(m,
    params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd = c(e.u = 1.0, e.g = 0.5))
  list(model = m, sol = sol)
}


test_that("unconditional variance_decomposition returns correct shape", {
  tm  <- make_test_model()
  vd  <- variance_decomposition(tm$sol)

  expect_s3_class(vd, "dsge_variance_decomposition")
  expect_equal(vd$type, "unconditional")
  expect_identical(vd$horizon, Inf)
  expect_equal(dim(vd$contribution),     c(2, 2))   # 2 obs, 2 shocks
  expect_equal(dim(vd$contribution_pct), c(2, 2))
  expect_identical(vd$obs_names,   c("p", "r"))
  expect_identical(vd$shock_names, c("e.u", "e.g"))
})


test_that("unconditional shares sum to 100 per variable", {
  tm  <- make_test_model()
  vd  <- variance_decomposition(tm$sol)
  totals <- rowSums(vd$contribution_pct)
  expect_true(all(abs(totals - 100) < 1e-8))
})


test_that("unconditional contribution_pct entries are in [0, 100]", {
  tm <- make_test_model()
  vd <- variance_decomposition(tm$sol)
  expect_true(all(vd$contribution_pct >= 0))
  expect_true(all(vd$contribution_pct <= 100 + 1e-10))
})


test_that("unconditional contributions sum to model_covariance variances", {
  tm  <- make_test_model()
  vd  <- variance_decomposition(tm$sol)
  mc  <- model_covariance(tm$sol)
  total_var <- rowSums(vd$contribution)
  # Diagonal of model-implied covariance equals total variance
  expect_equal(unname(total_var), unname(diag(mc$covariance)),
               tolerance = 1e-8)
})


test_that("fevd returns correct shape", {
  tm  <- make_test_model()
  fv  <- variance_decomposition(tm$sol, horizon = c(1, 4, 8, 20))

  expect_s3_class(fv, "dsge_variance_decomposition")
  expect_equal(fv$type, "fevd")
  expect_identical(fv$horizon, c(1L, 4L, 8L, 20L))
  expect_equal(dim(fv$contribution),     c(4, 2, 2))   # 4 horizons, 2 obs, 2 shocks
  expect_equal(dim(fv$contribution_pct), c(4, 2, 2))
})


test_that("fevd shares sum to 100 for each (horizon, variable)", {
  tm  <- make_test_model()
  fv  <- variance_decomposition(tm$sol, horizon = c(1, 4, 8, 20))
  for (h in seq_len(dim(fv$contribution_pct)[1])) {
    for (i in seq_len(dim(fv$contribution_pct)[2])) {
      total <- sum(fv$contribution_pct[h, i, ])
      expect_true(abs(total - 100) < 1e-8)
    }
  }
})


test_that("fevd at large horizon approaches unconditional decomposition", {
  tm   <- make_test_model()
  vd   <- variance_decomposition(tm$sol)
  fv   <- variance_decomposition(tm$sol, horizon = 300L)
  uncon_pct <- vd$contribution_pct
  fevd_pct  <- fv$contribution_pct[1, , ]
  expect_true(max(abs(uncon_pct - fevd_pct)) < 0.5)  # within 0.5 pp
})


test_that("invalid horizon errors", {
  tm <- make_test_model()
  expect_error(variance_decomposition(tm$sol, horizon = 0L),
               "must be >= 1")
})


test_that("dsge_fit method works", {
  tm <- make_test_model()
  # Generate some data and fit
  set.seed(42)
  H <- tm$sol$H; G <- tm$sol$G; M <- tm$sol$M
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
  fit <- estimate(tm$model, data = dat)

  vd  <- variance_decomposition(fit)
  expect_s3_class(vd, "dsge_variance_decomposition")
  expect_true(all(abs(rowSums(vd$contribution_pct) - 100) < 1e-8))
})


test_that("plot.dsge_variance_decomposition does not error (unconditional)", {
  tm <- make_test_model()
  vd <- variance_decomposition(tm$sol)
  expect_silent({
    pdf(tempfile())
    plot(vd)
    dev.off()
  })
})


test_that("plot.dsge_variance_decomposition does not error (FEVD)", {
  tm <- make_test_model()
  fv <- variance_decomposition(tm$sol, horizon = c(1, 4, 12))
  expect_silent({
    pdf(tempfile())
    plot(fv)
    dev.off()
  })
})


test_that("print.dsge_variance_decomposition works for both types", {
  tm <- make_test_model()
  vd <- variance_decomposition(tm$sol)
  fv <- variance_decomposition(tm$sol, horizon = c(1, 4))
  expect_output(print(vd), "Unconditional Variance Decomposition")
  expect_output(print(fv), "Forecast Error Variance Decomposition")
})


test_that("dsge_bayes method solves at the posterior mean", {
  skip_on_cran()
  tm <- make_test_model()
  sol <- tm$sol
  set.seed(7)
  TT <- 120
  xs <- matrix(0, TT, nrow(sol$H))
  y  <- matrix(0, TT, nrow(sol$G), dimnames = list(NULL, rownames(sol$G)))
  for (t in 2:TT) {
    xs[t, ] <- sol$H %*% xs[t - 1, ] + sol$M %*% rnorm(ncol(sol$M))
    y[t, ]  <- sol$G %*% xs[t, ]
  }
  fit <- bayes_dsge(tm$model, data = as.data.frame(y[, c("p", "r")]),
                    priors = list(kappa = prior("gamma", shape = 2, rate = 20),
                                  psi   = prior("normal", mean = 1.5, sd = 0.25),
                                  rhou  = prior("beta", shape1 = 2, shape2 = 2),
                                  rhog  = prior("beta", shape1 = 2, shape2 = 2)),
                    chains = 2, iter = 300, warmup = 150, seed = 1)

  vd <- variance_decomposition(fit)
  expect_s3_class(vd, "dsge_variance_decomposition")
  expect_identical(vd$obs_names, c("p", "r"))
  expect_identical(vd$shock_names, c("u", "g"))
  expect_equal(unname(rowSums(vd$contribution_pct)) / 100, c(1, 1),
               tolerance = 1e-10)

  fv <- variance_decomposition(fit, horizon = c(1, 8))
  expect_equal(unname(apply(fv$contribution_pct, c(1, 2), sum)) / 100,
               matrix(1, 2, 2), tolerance = 1e-10)

  # Same as solving at the posterior mean, shock SDs included
  pm <- apply(fit$posterior, 2, function(col) mean(as.numeric(col)))
  sol_pm <- solve_dsge(tm$model, params = pm[fit$free_parameters],
                       shock_sd = c(u = pm[["sd_e.u"]], g = pm[["sd_e.g"]]))
  expect_equal(vd$contribution, variance_decomposition(sol_pm)$contribution)

  # `params` replaces the posterior mean of the parameters it names
  vd_u <- variance_decomposition(fit, params = c(sd_e.u = 3 * pm[["sd_e.u"]]))
  expect_gt(vd_u$contribution_pct["p", "u"], vd$contribution_pct["p", "u"])
  expect_equal(variance_decomposition(fit, params = pm)$contribution,
               vd$contribution)
  expect_error(variance_decomposition(fit, params = c(beta = 0.98)),
               "only set estimated parameters")
  expect_error(variance_decomposition(fit, params = c(sd_e.u = -1)),
               "must be positive")
  expect_error(variance_decomposition(fit, params = c(rhou = 1.5)),
               "no stable solution")

  expect_silent({
    pdf(tempfile())
    plot(vd)
    plot(fv)
    dev.off()
  })
})
