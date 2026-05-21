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
