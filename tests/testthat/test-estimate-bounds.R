# Parameter bounds in estimate() (Dynare's estimated_params bounds)

ar1_model <- function() {
  dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.6))
}
ar1_data <- function() {
  set.seed(42)
  z <- numeric(200)
  for (i in 2:200) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  data.frame(y = z)
}

test_that("theta_bounds() maps bounds onto the optimisation vector", {
  tb <- dsge:::theta_bounds(list(rho = c(0.2, 0.9), sd_e.z = c(0.1, 2)),
                            free_params = "rho", shock_names = "z",
                            all_params = "rho")
  expect_equal(unname(tb$lower), c(0.2, log(0.1)))
  expect_equal(unname(tb$upper), c(0.9, log(2)))
  expect_equal(tb$shock, c(FALSE, TRUE))
  # open sides, a lower bound of 0 on a standard deviation, fixed parameters
  tb <- dsge:::theta_bounds(list(rho = c(NA, 0.9), sd_e.z = c(0, NA),
                                 beta = c(0.9, 1)),
                            free_params = "rho", shock_names = "z",
                            all_params = c("rho", "beta"))
  expect_equal(unname(tb$lower), c(-Inf, -Inf))
  expect_equal(unname(tb$upper), c(0.9, Inf))
  expect_null(dsge:::theta_bounds(list(beta = c(0.9, 1)), "rho", "z",
                                  c("rho", "beta")))
  expect_null(dsge:::theta_bounds(NULL, "rho", "z", "rho"))
  expect_error(dsge:::theta_bounds(list(phi = c(0, 1)), "rho", "z", "rho"),
               "unknown")
  expect_error(dsge:::theta_bounds(list(sd_e.z = c(NA, -1)), "rho", "z",
                                   "rho"), "must be positive")
  expect_error(dsge:::theta_bounds(list(rho = c(0.9, 0.2)), "rho", "z",
                                   "rho"), "not below")
})

test_that("bounds that do not bind leave the estimate unchanged", {
  m <- ar1_model()
  dat <- ar1_data()
  f0 <- estimate(m, data = dat)
  f1 <- estimate(m, data = dat, bounds = list(rho = c(0.5, 0.95)))
  expect_equal(coef(f1), coef(f0), tolerance = 1e-5)
  expect_equal(f1$se, f0$se, tolerance = 1e-3)
  expect_length(f1$at_bound, 0L)
  expect_identical(f0$at_bound, character(0))
})

test_that("an estimate at a bound gets no standard error", {
  m <- ar1_model()
  dat <- ar1_data()
  # The data favour rho = 0.8, below the bounds; the start is outside them
  expect_warning(expect_warning(
    fit <- estimate(m, data = dat, bounds = list(rho = c(0.85, 0.95))),
    "not strictly inside"), "at a bound: rho")
  expect_equal(unname(coef(fit)["rho"]), 0.85)
  expect_identical(fit$at_bound, "rho")
  expect_true(is.na(fit$se["rho"]))
  expect_true(is.finite(fit$se["sd(e.z)"]) && fit$se["sd(e.z)"] > 0)
  expect_equal(unname(fit$vcov["rho", ]), c(0, 0))
  expect_output(print(fit), "\\(at bound\\)")

  ir <- irf(fit, periods = 4)
  expect_true(all(is.finite(ir$data$se)))
  rv <- robust_vcov(fit)
  expect_true(is.na(rv$se["rho"]))
  expect_true(is.finite(rv$se["sd(e.z)"]))

  # A binding upper bound on a shock standard deviation
  expect_warning(
    fit_sd <- estimate(m, data = dat, shock_start = c(z = 0.5),
                       bounds = list(sd_e.z = c(NA, 0.9))),
    "at a bound: sd\\(e.z\\)")
  expect_equal(unname(coef(fit_sd)["sd(e.z)"]), 0.9, tolerance = 1e-8)
  expect_true(is.na(fit_sd$se["sd(e.z)"]))
  expect_true(is.finite(fit_sd$se["rho"]))
})

test_that("bounds need method L-BFGS-B", {
  expect_error(estimate(ar1_model(), data = ar1_data(), method = "BFGS",
                        bounds = list(rho = c(0.5, 0.95))), "L-BFGS-B")
  fit <- estimate(ar1_model(), data = ar1_data(), method = "L-BFGS-B",
                  bounds = list(rho = c(0.5, 0.95)))
  expect_s3_class(fit, "dsge_fit")
})

test_that("estimate() uses the bounds of an imported Dynare model", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters rho;
    rho = 0.97;
    model(linear);
      y = rho * y(-1) + e;
    end;
    shocks; var e; stderr 1; end;
    varobs y;
    estimated_params;
      rho, 0.97, 0.95, 0.99;
      stderr e, 1, 0.1, 5;
    end;
  ")
  set.seed(3)
  dat <- data.frame(y = as.numeric(stats::arima.sim(list(ar = 0.7), n = 200)))
  expect_warning(fit <- estimate(m, data = dat), "at a bound: rho")
  expect_equal(unname(coef(fit)["rho"]), 0.95)
  expect_true(is.na(fit$se["rho"]))

  free <- estimate(m, data = dat, bounds = NULL)
  expect_lt(unname(coef(free)["rho"]), 0.9)
  expect_length(free$at_bound, 0L)
})
