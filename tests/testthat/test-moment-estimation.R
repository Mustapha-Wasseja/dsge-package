# Tests for gmm_estimate() and smm_estimate()

make_nk_data <- function(TT = 200L, seed = 1L) {
  nk <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
  sol_true <- solve_dsge(nk,
    params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd = c(e.u = 1, e.g = 0.5))
  set.seed(seed)
  H <- sol_true$H; G <- sol_true$G; M <- sol_true$M
  xst <- matrix(0, TT, nrow(H))
  y <- matrix(0, TT, nrow(G))
  for (t in 2:TT) {
    e <- rnorm(ncol(M)) * c(1, 0.5)
    xst[t, ] <- as.numeric(H %*% xst[t - 1, ] + M %*% e)
    y[t, ]   <- as.numeric(G %*% xst[t, ])
  }
  colnames(y) <- rownames(G)
  list(model = nk, data = as.data.frame(y[, nk$variables$observed]))
}


test_that("gmm_estimate runs and returns expected structure", {
  d <- make_nk_data()
  est <- gmm_estimate(d$model, d$data,
    moments = c("sd:p","sd:r","ac1:p","ac1:r","cor:p:r"),
    params_start   = c(kappa = 0.15, psi = 2.0, rhou = 0.6, rhog = 0.8),
    shock_sd_start = c(e.u = 0.8, e.g = 0.7),
    control = list(maxit = 200))
  expect_s3_class(est, "dsge_gmm")
  expect_named(est$params, c("kappa","psi","rhou","rhog"))
  expect_named(est$shock_sd, c("e.u","e.g"))
  expect_true(is.finite(est$objective))
  expect_true(nrow(est$moments) == 5)
})


test_that("gmm recovers true parameters within a reasonable tolerance", {
  d <- make_nk_data(TT = 500L, seed = 2L)
  est <- gmm_estimate(d$model, d$data,
    moments = c("sd:p","sd:r","var:p","var:r","ac1:p","ac1:r","cor:p:r"),
    params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd_start = c(e.u = 1.0, e.g = 0.5),
    control = list(maxit = 500))
  expect_lt(est$objective, 1.0)
  # rhou, rhog should be close (these are well-identified by ac1)
  expect_lt(abs(est$params["rhou"] - 0.7), 0.3)
  expect_lt(abs(est$params["rhog"] - 0.9), 0.3)
})


test_that("two-step GMM works", {
  d <- make_nk_data(TT = 200L)
  est <- gmm_estimate(d$model, d$data,
    moments = c("sd:p","sd:r","ac1:p","ac1:r"),
    params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd_start = c(e.u = 1, e.g = 0.5),
    weight = "optimal",
    control = list(maxit = 100))
  expect_true(est$two_step)
  expect_true(is.finite(est$objective))
})


test_that("smm_estimate runs", {
  d <- make_nk_data(TT = 200L)
  est <- smm_estimate(d$model, d$data,
    moments = c("sd:p","sd:r","ac1:p"),
    params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd_start = c(e.u = 1, e.g = 0.5),
    sim_periods = 500L, sim_replic = 3L,
    seed = 1L,
    control = list(maxit = 50))
  expect_s3_class(est, "dsge_smm")
  expect_s3_class(est, "dsge_gmm")
})


test_that("invalid moment spec errors clearly", {
  d <- make_nk_data()
  expect_error(
    gmm_estimate(d$model, d$data,
      moments = c("zzz:p"),
      params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
      shock_sd_start = c(e.u = 1, e.g = 0.5)),
    "Unrecognised moment spec"
  )
})


test_that("print method runs", {
  d <- make_nk_data()
  est <- gmm_estimate(d$model, d$data,
    moments = c("sd:p","ac1:p"),
    params_start   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
    shock_sd_start = c(e.u = 1, e.g = 0.5),
    control = list(maxit = 30))
  expect_output(print(est), "Method-of-Moments estimation")
})
