# Tests for anticipated (news) shocks via multi-period shock paths in
# perfect_foresight() and perfect_foresight_nonlinear().

make_ar_model <- function() {
  mod <- dsge_model(
    obs(p ~ x),
    state(x ~ rho * x),
    fixed = list(rho = 0.9)
  )
  sol <- solve_dsge(mod, params = c(rho = 0.9), shock_sd = c(x = 0.01))
  list(mod = mod, sol = sol)
}


test_that("vector shock path is accepted by perfect_foresight()", {
  tm <- make_ar_model()
  pf <- perfect_foresight(tm$sol,
                          shocks  = list(x = c(0, 0, 0, 0.01)),
                          horizon = 20)
  expect_s3_class(pf, "dsge_perfect_foresight")
  # Shock at period 4 -> non-zero at period 4 in shock_path
  expect_equal(pf$shock_path[1, "x"], 0, ignore_attr = TRUE)
  expect_equal(pf$shock_path[2, "x"], 0, ignore_attr = TRUE)
  expect_equal(pf$shock_path[3, "x"], 0, ignore_attr = TRUE)
  expect_equal(pf$shock_path[4, "x"], 0.01, ignore_attr = TRUE)
  expect_equal(pf$shock_path[5, "x"], 0, ignore_attr = TRUE)
})


test_that("linear PF treats shock path as period-by-period surprises", {
  # In the linearized perfect_foresight() with recursive policy, a shock at
  # t = k produces a response starting at t = k (no anticipation at t < k).
  tm <- make_ar_model()
  # Surprise at t = 1
  pf_surprise <- perfect_foresight(tm$sol,
                                   shocks  = list(x = 0.01),
                                   horizon = 20)
  # "News": shock arrives at t = 5
  pf_late <- perfect_foresight(tm$sol,
                               shocks  = list(x = c(0, 0, 0, 0, 0.01)),
                               horizon = 20)
  # Response is delayed: zero for t = 1..4 then matches surprise IRF shifted
  expect_lt(max(abs(pf_late$controls[1:4, ])), 1e-12)
  expect_gt(max(abs(pf_late$controls[5:9, "p"])), 1e-6)
  # The shifted-by-4 response should match the early-impact response
  expect_equal(pf_late$controls[5:9, "p"],
               pf_surprise$controls[1:5, "p"],
               tolerance = 1e-10)
})


test_that("multi-period shock paths work in perfect_foresight_nonlinear()", {
  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C", endo_state = "K", exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9)
  )

  # Anticipated TFP shock: announced at t=1, arrives at t=3
  pf_news <- perfect_foresight_nonlinear(rbc,
    params   = c(rho = 0.9),
    shock_sd = c(Z = 0.01),
    shocks   = list(Z = c(0, 0, 0.05)),
    horizon  = 30)
  expect_s3_class(pf_news, "dsge_perfect_foresight")
  expect_true(pf_news$converged)
  # Consumers should adjust C at t=1 in anticipation
  expect_true(abs(pf_news$controls[1, "C"]) > 1e-6)
})


test_that("zero shock vector produces zero response", {
  tm <- make_ar_model()
  pf <- perfect_foresight(tm$sol,
                          shocks  = list(x = rep(0, 10)),
                          horizon = 20)
  expect_true(max(abs(pf$states)) < 1e-12)
  expect_true(max(abs(pf$controls)) < 1e-12)
})
