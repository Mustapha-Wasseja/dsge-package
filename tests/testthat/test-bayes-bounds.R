# Parameter bounds in bayes_dsge() (Dynare's estimated_params bounds)

bounded <- function(p, b) { p$bounds <- b; p }

test_that("bounded transforms map the real line onto the allowed interval", {
  cases <- list(
    list(p = bounded(prior("beta", shape1 = 12, shape2 = 12), c(0.5, 0.95)),
         iv = c(0.5, 0.95), x = c(0.51, 0.7, 0.94)),
    list(p = bounded(prior("normal", mean = 4, sd = 1.5), c(2, NA)),
         iv = c(2, Inf), x = c(2.1, 4, 20)),
    list(p = bounded(prior("gamma", shape = 2, rate = 1), c(NA, 3)),
         iv = c(0, 3), x = c(0.01, 1, 2.9)),
    list(p = bounded(prior("normal", mean = 0, sd = 2), c(NA, 1)),
         iv = c(-Inf, 1), x = c(-5, 0, 0.9)),
    list(p = bounded(prior("inv_gamma1", s = 0.01, nu = 2), c(0.01, 3)),
         iv = c(0.01, 3), x = c(0.02, 0.5, 2.9))
  )
  for (cs in cases) {
    expect_equal(dsge:::param_interval(cs$p), cs$iv)
    u <- vapply(cs$x, dsge:::to_unconstrained, 0, prior_obj = cs$p)
    expect_equal(vapply(u, dsge:::from_unconstrained, 0, prior_obj = cs$p),
                 cs$x, tolerance = 1e-12)
    # Every real number maps inside the interval
    z <- vapply(c(-50, -5, 0, 5, 50), dsge:::from_unconstrained, 0,
                prior_obj = cs$p)
    expect_true(all(z >= cs$iv[1] & z <= cs$iv[2]))
    # Log Jacobian matches the numerical derivative
    for (uu in c(-2, 0.3, 1.7)) {
      h <- 1e-6
      d <- (dsge:::from_unconstrained(uu + h, cs$p) -
              dsge:::from_unconstrained(uu - h, cs$p)) / (2 * h)
      expect_equal(dsge:::log_jacobian(uu, cs$p), log(abs(d)),
                   tolerance = 1e-6)
    }
  }
})

test_that("priors without bounds keep their transformations", {
  p <- prior("beta", shape1 = 2, shape2 = 2)
  expect_identical(dsge:::to_unconstrained(0.3, p), stats::qlogis(0.3))
  expect_identical(dsge:::from_unconstrained(0.4, p), stats::plogis(0.4))
  g <- prior("gamma", shape = 2, rate = 1)
  expect_identical(dsge:::to_unconstrained(2, g), log(2))
  expect_identical(dsge:::log_jacobian(0.7, g), 0.7)
  n <- prior("normal", mean = 0, sd = 1)
  expect_identical(dsge:::from_unconstrained(-1.3, n), -1.3)
  expect_identical(dsge:::log_jacobian(-1.3, n), 0)
})

test_that("dprior is unchanged inside the bounds and -Inf outside", {
  p <- prior("beta", shape1 = 12, shape2 = 12)
  pb <- bounded(p, c(0.5, 0.95))
  expect_identical(dsge:::dprior(pb, 0.7), dsge:::dprior(p, 0.7))
  expect_identical(dsge:::dprior(pb, 0.45), -Inf)
  expect_identical(dsge:::dprior(pb, 0.96), -Inf)
  expect_equal(dsge:::dprior(pb, c(0.4, 0.7)), c(-Inf, dsge:::dprior(p, 0.7)))
  expect_output(print(pb), "bounded to \\[0.5, 0.95\\]")
  expect_output(print(bounded(p, c(0.5, NA))), "bounded to \\[0.5, Inf\\]")
})

test_that("apply_bounds() validates the bounds", {
  pl <- list(rho = prior("beta", shape1 = 2, shape2 = 2),
             sd_e.z = prior("inv_gamma", shape = 0.01, scale = 0.01))
  out <- dsge:::apply_bounds(pl, list(rho = c(0.2, 0.9), sd_e.z = c(NA, 2)))
  expect_equal(out$rho$bounds, c(0.2, 0.9))
  expect_equal(out$sd_e.z$bounds, c(NA, 2))
  expect_null(dsge:::apply_bounds(pl, list(rho = c(-Inf, Inf)))$rho$bounds)
  expect_identical(dsge:::apply_bounds(pl, NULL), pl)
  expect_error(dsge:::apply_bounds(pl, list(phi = c(0, 1))), "unknown")
  expect_error(dsge:::apply_bounds(pl, list(rho = c(0.9, 0.2))), "not below")
  expect_error(dsge:::apply_bounds(pl, list(rho = c(1.2, 2))), "support")
  expect_error(dsge:::apply_bounds(pl, list(rho = 0.5)), "c\\(lower, upper\\)")
  expect_error(dsge:::apply_bounds(pl, list(c(0, 1))), "named list")
})

test_that("bayes_dsge() keeps every draw within the bounds", {
  skip_on_cran()
  m <- dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.5))
  set.seed(42)
  z <- numeric(200)
  for (i in 2:200) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  dat <- data.frame(y = z)
  pri <- list(rho = prior("beta", shape1 = 2, shape2 = 2))

  # The data favour rho = 0.8, below the bounds; the start 0.5 is outside
  expect_warning(
    fit <- bayes_dsge(m, data = dat, priors = pri, chains = 1L, iter = 800L,
                      warmup = 300L, seed = 1, shock_start = c(z = 0.5),
                      bounds = list(rho = c(0.85, 0.95), sd_e.z = c(NA, 0.9))),
    "rho = 0.5 lies outside its bounds")
  rho <- fit$posterior[, "rho", 1]
  sdz <- fit$posterior[, "sd_e.z", 1]
  expect_true(all(rho > 0.85 & rho < 0.95))
  expect_true(all(sdz < 0.9))
  expect_gt(mean(rho < 0.87), 0.3)        # mass piles up at the lower bound
  expect_equal(fit$priors$rho$bounds, c(0.85, 0.95))

  fit0 <- bayes_dsge(m, data = dat, priors = pri, chains = 1L, iter = 800L,
                     warmup = 300L, seed = 1)
  expect_gt(mean(fit0$posterior[, "rho", 1] < 0.85), 0.5)
  expect_null(fit0$priors$rho$bounds)
})

test_that("read_dynare() imports the bounds and bayes_dsge() uses them", {
  m <- read_dynare(text = "
    var y x;
    varexo e;
    parameters rho;
    rho = 0.9;
    model(linear);
      y = rho * y(-1) + e;
      x = 2 * y;
    end;
    shocks;
      var e; stderr 1;
      var x; stderr 0.3;
    end;
    varobs y x;
    estimated_params;
      rho, 0.9, 0.5, 0.95, beta_pdf, 0.8, 0.1;
      stderr e, 1, 0.1, 5, inv_gamma_pdf, 1, inf;
      stderr x, 0.3, 0.01, 1, inv_gamma_pdf, 0.2, inf;
    end;
    estimated_params_bounds;
      stderr e, 0.2, 3;
    end;
  ")
  expect_equal(m$bounds, list(rho = c(0.5, 0.95), sd_e.e = c(0.2, 3),
                              sd_e.x_me = c(0.01, 1)))
  expect_setequal(names(m$bounds), names(m$priors))
  expect_false(any(grepl("estimated_params_bounds", m$notes)))
  expect_output(print(m), "Bounds:")

  m0 <- read_dynare(text = "
    var y; varexo e; parameters rho; rho = 0.5;
    model(linear); y = rho * y(-1) + e; end;
    shocks; var e; stderr 1; end;
    varobs y;
    estimated_params; rho, beta_pdf, 0.5, 0.2; end;
  ")
  expect_null(m0$bounds)

  skip_on_cran()
  set.seed(2)
  y <- as.numeric(stats::arima.sim(list(ar = 0.9), n = 120))
  dat <- data.frame(y = y + stats::rnorm(120, sd = 0.2),
                    x = 2 * y + stats::rnorm(120, sd = 0.3))
  fit <- bayes_dsge(m, data = dat, chains = 1L, iter = 400L, warmup = 200L,
                    seed = 1)
  expect_equal(fit$priors$rho$bounds, c(0.5, 0.95))
  expect_equal(fit$priors$sd_e.e$bounds, c(0.2, 3))
  rho <- fit$posterior[, "rho", 1]
  expect_true(all(rho > 0.5 & rho < 0.95))
  expect_true(all(fit$posterior[, "sd_e.e", 1] > 0.2))

  fit_free <- bayes_dsge(m, data = dat, chains = 1L, iter = 300L,
                         warmup = 150L, seed = 1, bounds = NULL)
  expect_null(fit_free$priors$rho$bounds)
})
