# Tests for read_dynare()

irf_path <- function(sol, response, impulse, periods = 10L) {
  d <- irf(sol, periods = periods, se = FALSE)$data
  d$value[d$response == response & d$impulse == impulse]
}

ar1_text <- "
  var y;
  varexo e;
  parameters rho;
  rho = 0.9;
  model;
    y = rho * y(-1) + e;
  end;
  shocks;
    var e; stderr 0.01;
  end;
"


test_that("a simple AR(1) is imported with Dynare timing", {
  m <- read_dynare(text = ar1_text)
  expect_s3_class(m, "dsge_dynare")
  expect_s3_class(m$model, "dsgenl_model")
  expect_equal(m$variables, "y")
  expect_equal(m$shocks, "e")
  expect_equal(m$shock_sd, c(e = 0.01))
  expect_equal(m$aux$name, "y_lag1")
  expect_equal(m$observed, character(0))
  expect_true(any(grepl("No varobs", m$notes)))

  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "y", "e", 5L), 0.01 * 0.9^(0:5),
               tolerance = 1e-8)
})


test_that("the bundled RBC example matches a hand-written dsge model", {
  m <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  expect_setequal(m$aux$name, c("k_lag1", "a_lag1"))
  expect_length(m$commands, 3L)

  sol <- solve_dsge(m)
  k_ss <- ((1 / 0.99 - 1 + 0.025) / 0.33)^(1 / (0.33 - 1))
  expect_equal(unname(sol$steady_state["k"]), k_ss, tolerance = 1e-8)

  hand <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z(+1)) * K(+1)^(alpha - 1) + 1 - delta)",
    "Y = exp(Z) * K^alpha",
    "I = Y - C",
    "K(+1) = I + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "Y", unobserved = c("C", "I"),
    exo_state = "Z", endo_state = "K",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025, rho = 0.95),
    ss_guess = c(Y = 3, C = 2.3, I = 0.7, K = 28, Z = 0)
  )
  hs <- solve_dsge(hand, shock_sd = c(Z = 0.01))
  for (v in c("y", "c", "i")) {
    expect_equal(irf_path(sol, v, "e"), irf_path(hs, toupper(v), "Z"),
                 tolerance = 1e-6)
  }
})


test_that("predetermined_variables gives the same dynamics as lag timing", {
  common <- "
    var y c k z;
    varexo e;
    parameters alpha beta delta rho;
    alpha = 0.33; beta = 0.99; delta = 0.025; rho = 0.9;
    initval; k = 28; c = 2.3; y = 3; z = 0; end;
    shocks; var e; stderr 0.01; end;
  "
  lagged <- read_dynare(text = paste(common, "
    model;
      1/c = beta / c(+1) * (alpha * exp(z(+1)) * k^(alpha - 1) + 1 - delta);
      y = exp(z) * k(-1)^alpha;
      k = y - c + (1 - delta) * k(-1);
      z = rho * z(-1) + e;
    end;"))
  predet <- read_dynare(text = paste(common, "
    predetermined_variables k;
    model;
      1/c = beta / c(+1) * (alpha * exp(z(+1)) * k(+1)^(alpha - 1) + 1 - delta);
      y = exp(z) * k^alpha;
      k(+1) = y - c + (1 - delta) * k;
      z = rho * z(-1) + e;
    end;"))
  s1 <- solve_dsge(lagged)
  s2 <- solve_dsge(predet)
  expect_equal(irf_path(s1, "c", "e"), irf_path(s2, "c", "e"),
               tolerance = 1e-6)
  expect_equal(irf_path(s1, "y", "e"), irf_path(s2, "y", "e"),
               tolerance = 1e-6)
})


test_that("long leads, long lags and lagged shocks get auxiliary variables", {
  m <- read_dynare(text = "
    var x y;
    varexo e;
    parameters a1 a2 b1 b2;
    a1 = 0.5; a2 = 0.2; b1 = 0.3; b2 = 0.1;
    model(linear);
      x = a1 * x(-1) + a2 * x(-2) + e + 0.5 * e(-1);
      y = b1 * y(+1) + b2 * y(+2) + x;
    end;
    shocks; var e; stderr 1; end;
  ")
  expect_setequal(m$aux$name, c("y_lead1", "x_lag1", "x_lag2", "e_lag1"))
  expect_true("y_lead1" %in% m$model$controls)
  expect_true(all(c("x_lag1", "x_lag2", "e_lag1") %in% m$model$states))

  sol <- solve_dsge(m)
  # x_t = a1 x_{t-1} + a2 x_{t-2} + e_t + 0.5 e_{t-1}
  expected <- numeric(6)
  expected[1] <- 1
  expected[2] <- 0.5 * expected[1] + 0.5
  for (t in 3:6) expected[t] <- 0.5 * expected[t - 1] + 0.2 * expected[t - 2]
  expect_equal(irf_path(sol, "x", "e", 5L), expected, tolerance = 1e-8)
})


test_that("comments, equation tags and model-local variables are handled", {
  m <- read_dynare(text = "
    /* block
       comment */
    var y pi;       // line comment
    varexo u v;     % MATLAB comment
    parameters beta kappa rho;
    beta = 0.99; kappa = 0.1; rho = 0.5;
    model(linear);
      # slope = kappa / (1 - beta * rho);
      [name = 'Phillips curve']
      pi = beta * pi(+1) + kappa * y + u;
      [name = 'demand', mcp = 'x > 0']
      y = rho * y(-1) + v;
    end;
    shocks; var u; stderr 0.1; var v; stderr 0.2; end;
  ")
  expect_equal(m$variables, c("y", "pi"))
  expect_true(any(grepl("pi(+1)", m$model$eq_strings, fixed = TRUE)))
  expect_equal(unname(m$model$ss_guess[c("y", "pi")]), c(0, 0))
  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "pi", "u", 0L), 0.1, tolerance = 1e-8)
})


test_that("model-local variables are inlined into equations", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters a b;
    a = 0.5; b = 2;
    model;
      # ab = a / b;
      y = ab * y(-1) + e;
    end;
    shocks; var e; stderr 1; end;
  ")
  expect_match(m$model$eq_strings[1], "(a / b)", fixed = TRUE)
  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "y", "e", 2L), 0.25^(0:2), tolerance = 1e-8)
})


test_that("steady_state_model becomes an ss_function that uses parameters", {
  m <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  ss_fn <- m$model$ss_function
  expect_true(is.function(ss_fn))
  p <- m$params[m$model$parameters]
  ss1 <- ss_fn(p)
  expect_equal(unname(ss1["k_lag1"]), unname(ss1["k"]))
  expect_equal(unname(ss1["e"]), 0)
  p["beta"] <- 0.98
  expect_lt(ss_fn(p)["k"], ss1["k"])
})


test_that("correlated shocks are orthogonalised as in Dynare", {
  m <- read_dynare(text = "
    var y1 y2;
    varexo e1 e2;
    parameters rho;
    rho = 0;
    model(linear);
      y1 = rho * y1(-1) + e1;
      y2 = rho * y2(-1) + e2;
    end;
    shocks;
      var e1; stderr 0.2;
      var e2; stderr 0.5;
      corr e1, e2 = 0.6;
    end;
  ")
  expect_equal(unname(m$shock_sd), c(0.2, 0.5 * sqrt(1 - 0.6^2)))
  expect_true(any(grepl("orthogonalised", m$notes)))
  sol <- solve_dsge(m)
  # Cholesky: a one-s.d. e1 shock moves y2 by corr * sd2 on impact
  expect_equal(irf_path(sol, "y2", "e1", 0L), 0.6 * 0.5, tolerance = 1e-8)
  expect_equal(irf_path(sol, "y2", "e2", 0L), 0.5 * sqrt(1 - 0.36),
               tolerance = 1e-8)
  expect_equal(irf_path(sol, "y1", "e2", 0L), 0, tolerance = 1e-10)

  # The same covariance given as a variance statement
  m2 <- read_dynare(text = sub("corr e1, e2 = 0.6", "var e1, e2 = 0.06",
                               "
    var y1 y2; varexo e1 e2; parameters rho; rho = 0;
    model(linear); y1 = rho * y1(-1) + e1; y2 = rho * y2(-1) + e2; end;
    shocks; var e1; stderr 0.2; var e2; stderr 0.5; corr e1, e2 = 0.6; end;
  ", fixed = TRUE))
  expect_equal(m2$shock_sd, m$shock_sd)
})


test_that("estimated_params are translated into dsge priors", {
  m <- read_dynare(text = "
    var y x;
    varexo e u;
    parameters rho phi kappa sig tau;
    rho = 0.9; phi = 1.5; kappa = 0.2; sig = 1; tau = 0.5;
    model(linear);
      y = rho * y(-1) + e;
      x = phi * y + kappa * x(-1) + sig * tau * u;
    end;
    shocks; var e; stderr 0.01; end;
    varobs y x;
    estimated_params;
      rho, beta_pdf, 0.7, 0.1;
      phi, 1.4, gamma_pdf, 1.5, 0.25;
      kappa, 0.2, 0, 1, normal_pdf, 0.2, 0.05;
      stderr e, inv_gamma_pdf, 0.01, inf;
      stderr u, inv_gamma2_pdf, 0.02, 0.01;
      sig, uniform_pdf, , , 0, 2;
      tau, weibull_pdf, 1, 0.5;
    end;
    estimated_params_init;
      rho, 0.8;
    end;
  ")
  pr <- m$priors
  expect_setequal(names(pr), c("rho", "phi", "kappa", "sd_e.e", "sd_e.u",
                               "sig"))
  k <- 0.7 * 0.3 / 0.01 - 1
  expect_equal(pr$rho$params$shape1, 0.7 * k)
  expect_equal(pr$rho$params$shape2, 0.3 * k)
  expect_equal(pr$phi$params$shape, 1.5^2 / 0.25^2)
  expect_equal(pr$phi$params$rate, 1.5 / 0.25^2)
  expect_equal(pr$kappa$params, list(mean = 0.2, sd = 0.05))
  expect_equal(pr[["sd_e.e"]]$params, list(shape = 2, scale = 0.01))
  a <- 2 + 0.02^2 / 0.01^2
  expect_equal(pr[["sd_e.u"]]$params, list(shape = a, scale = 0.02 * (a - 1)))
  expect_equal(pr$sig$params, list(min = 0, max = 2))

  # Free parameters are the estimated ones, with initial values as starts
  expect_setequal(m$model$free_parameters,
                  c("rho", "phi", "kappa", "sig", "tau"))
  expect_equal(m$model$start$rho, 0.8)
  expect_equal(m$model$start$phi, 1.4)
  expect_equal(m$model$start$tau, 0.5)
  expect_true(any(grepl("tau", m$notes)))
  expect_true(any(grepl("approximated", m$notes)))

  # Shock with no declared variance but an estimated stderr gets its mean
  expect_equal(unname(m$shock_sd["u"]), 0.02)
})


test_that("Dynare math functions are translated", {
  m <- read_dynare(text = "
    var y z;
    varexo e;
    parameters a;
    a = 0.5;
    model;
      z = a * z(-1) + e;
      ln(y) = normcdf(z) - 0.5 + 0 * erf(z) + 0 * exp(1);
    end;
    initval; y = 1; end;
    shocks; var e; stderr 0.1; end;
  ")
  expect_true(any(grepl("log(y)", m$model$eq_strings, fixed = TRUE)))
  expect_true(any(grepl("stats::pnorm", m$model$eq_strings, fixed = TRUE)))
  sol <- solve_dsge(m)
  expect_equal(unname(sol$steady_state["y"]), 1, tolerance = 1e-8)
  expect_equal(irf_path(sol, "y", "e", 0L), stats::dnorm(0) * 0.1,
               tolerance = 1e-6)
})


test_that("observed variables follow varobs or an override", {
  base <- "
    var a b c;
    varexo e1 e2;
    parameters r;
    r = 0.5;
    model(linear);
      a = r * a(-1) + e1;
      b = r * b(-1) + e2;
      c = a + b;
    end;
  "
  m <- read_dynare(text = base)
  expect_equal(m$observed, character(0))
  expect_equal(m$model$variables$observed, character(0))
  expect_true(any(grepl("No varobs", m$notes)))
  expect_s3_class(solve_dsge(m, shock_sd = c(e1 = 1, e2 = 1)),
                  "dsge_solution")

  m <- read_dynare(text = paste(base, "varobs c a;"))
  expect_equal(m$observed, c("c", "a"))

  # Fewer observed variables than shocks is allowed
  m <- read_dynare(text = paste(base, "varobs c;"))
  expect_equal(m$observed, "c")
  expect_equal(m$model$variables$observed, "c")

  # More observed variables than shocks is truncated with a note
  m <- read_dynare(text = paste(base, "varobs a b c;"))
  expect_equal(m$observed, c("a", "b"))
  expect_true(any(grepl("singular", m$notes)))

  m <- read_dynare(text = base, observed = c("b", "c"))
  expect_equal(m$observed, c("b", "c"))
  expect_true(any(grepl("standard deviations are 0", m$notes)))
})


test_that("constants and undeclared assignments are usable in later blocks", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters rho;
    rho = 0.9;
    scale = 0.02;
    model;
      y = rho * y(-1) + e;
    end;
    shocks; var e = scale^2; end;
  ")
  expect_equal(unname(m$shock_sd), 0.02)
  expect_false("scale" %in% names(m$params))
})


test_that("unsupported input fails with informative errors", {
  expect_error(read_dynare(text = "@#if 1\nvar y;"), "not closed")
  expect_error(read_dynare(text = "@#foo x\nvar y;"), "Unsupported macro")
  expect_error(read_dynare(text = "var y; varexo e; model; y = e"),
               "missing ';'")
  expect_error(read_dynare(text = "var y; varexo e; model; y = e;"),
               "not closed")
  expect_error(read_dynare(text = "
    var y; varexo e; parameters a; a = 1;
    model; y = a * y(-1) + b + e; end;"), "Undeclared name")
  expect_error(read_dynare(text = "
    var y; varexo e; parameters a;
    model; y = a * y(-1) + e; end;"), "No value assigned")
  expect_error(read_dynare(text = "
    var y; varexo e; model; y = EXPECTATION(-1)(y) + e; end;"), "EXPECTATION")
  expect_error(read_dynare(text = "trend_var(growth_factor = g) A; var y;"),
               "trend_var")
  expect_error(read_dynare("no/such/file.mod"), "File not found")
})


test_that("print method summarises the import", {
  m <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  out <- capture.output(print(m))
  expect_true(any(grepl("Endogenous", out)))
  expect_true(any(grepl("k_lag1", out)))
  expect_true(any(grepl("stoch_simul", out)))
})


test_that("imported models can be solved at second order", {
  m <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  sol2 <- solve_dsge(m, order = 2)
  expect_s3_class(sol2, "dsge_solution")
})


test_that("bayes_dsge uses the imported priors by default", {
  skip_on_cran()
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters rho;
    rho = 0.5;
    model(linear);
      y = rho * y(-1) + e;
    end;
    shocks; var e; stderr 1; end;
    varobs y;
    estimated_params;
      rho, beta_pdf, 0.5, 0.2;
      stderr e, inv_gamma_pdf, 1, inf;
    end;
  ")
  set.seed(1)
  y <- as.numeric(stats::arima.sim(list(ar = 0.7), n = 150))
  fit <- bayes_dsge(m, data = data.frame(y = y), chains = 1L, iter = 400L,
                    warmup = 200L, seed = 1)
  expect_s3_class(fit, "dsge_bayes")
  expect_true("rho" %in% fit$free_parameters)
})


test_that("edge cases: element-wise operators, reserved names, unused priors", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters a b c0;
    a = 0.5; b = 2; c0 = 0;
    model;
      y = a .* y(-1) + e ./ b .^ 1;
    end;
    steady_state_model;
      y = e / (1 - a);
    end;
    shocks; var e; stderr 1; end;
    estimated_params;
      a, normal_pdf, 0.5, 0.1;
      b, normal_pdf, 2, 0.1;
      c0, normal_pdf, 0, 1;
    end;
  ")
  expect_false(any(grepl(".*", m$model$eq_strings, fixed = TRUE)))
  expect_equal(unname(m$model$ss_function(c(a = 0.5, b = 2))["y"]), 0)
  expect_setequal(names(m$priors), c("a", "b"))
  expect_true(any(grepl("c0", m$notes)))
  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "y", "e", 1L), c(0.5, 0.25), tolerance = 1e-8)

  expect_error(read_dynare(text = "var in; varexo e; model; in = e; end;"),
               "reserved in R")
})


test_that("macro directives are expanded", {
  m <- read_dynare(text = "
    @#define N = 3
    @#define names = [\"a\", \"b\", \"c\"]
    @#define linear = true
    var @{names[1]}
    @#for i in 2:N
      @{names[i]}
    @#endfor
    ;
    varexo e;
    parameters rho;
    rho = 0.5;
    @#if linear
    model(linear);
    @#else
    model;
    @#endif
      a = rho * a(-1) + e;
    @#for (v, w) in [(\"b\", \"a\"), (\"c\", \"b\")]
      @{v} = 0.5 * @{w}(-1);
    @#endfor
    end;
    @#ifdef extra
    varobs a;
    @#endif
    shocks; var e; stderr 1; end;
  ")
  expect_equal(m$variables, c("a", "b", "c"))
  expect_equal(m$observed, character(0))
  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "c", "e", 3L), c(0, 0, 0.25, 0.125),
               tolerance = 1e-8)

  m2 <- read_dynare(text = "
    var y; varexo e; parameters r;
    @#ifndef r_value
    @#define r_value = 0.9
    @#endif
    r = @{r_value};
    model; y = r * y(-1) + e; end;
  ", defines = list(r_value = 0.4))
  expect_equal(unname(m2$params["r"]), 0.4)
})


test_that("measurement errors become observation shocks", {
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
      rho, beta_pdf, 0.8, 0.1;
      stderr e, inv_gamma_pdf, 1, inf;
      stderr y, inv_gamma_pdf, 0.2, inf;
    end;
  ")
  expect_equal(m$measurement_errors, c("y", "x"))
  expect_equal(m$observed, c("y", "x"))
  expect_equal(m$model$variables$observed, c("y_obs", "x_obs"))
  expect_setequal(m$model$variables$exo_state, c("e", "y_me", "x_me"))
  expect_equal(unname(m$shock_sd[c("e", "y_me", "x_me")]), c(1, 0.2, 0.3))
  expect_setequal(names(m$priors), c("rho", "sd_e.e", "sd_e.y_me"))
  expect_equal(m$data_map, c(y = "y_obs", x = "x_obs"))

  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "x_obs", "x_me", 0L), 0.3, tolerance = 1e-10)
  expect_equal(irf_path(sol, "x", "x_me", 0L), 0, tolerance = 1e-10)
  expect_equal(irf_path(sol, "x_obs", "e", 2L), irf_path(sol, "x", "e", 2L),
               tolerance = 1e-10)

  skip_on_cran()
  set.seed(2)
  y <- as.numeric(stats::arima.sim(list(ar = 0.9), n = 120))
  dat <- data.frame(y = y + stats::rnorm(120, sd = 0.2),
                    x = 2 * y + stats::rnorm(120, sd = 0.3))
  fit <- bayes_dsge(m, data = dat, chains = 1L, iter = 300L,
                    warmup = 150L, seed = 1)
  expect_s3_class(fit, "dsge_bayes")
})


test_that("varexo_det and leads of shocks are supported", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    varexo_det d;
    parameters rho;
    rho = 0.5;
    model(linear);
      y = rho * y(-1) + e + e(+1) + d;
    end;
    shocks;
      var e; stderr 1;
      var d; periods 1:2 4; values 0.5 0.25;
    end;
  ")
  expect_equal(m$shocks, "e")
  expect_equal(m$shocks_det, "d")
  expect_equal(unname(m$shock_sd["d"]), 0)
  expect_equal(unname(m$shock_paths[, "d"]), c(0.5, 0.5, 0, 0.25))
  sol <- solve_dsge(m)
  # E_t e(+1) = 0, so the IRF is that of the AR(1)
  expect_equal(irf_path(sol, "y", "e", 3L), 0.5^(0:3), tolerance = 1e-8)
})


test_that("STEADY_STATE() refers to the steady state at current parameters", {
  txt <- "
    var y z;
    varexo e;
    parameters a rho;
    a = 2; rho = 0.5;
    model;
      z = rho * z(-1) + e;
      y = STEADY_STATE(y) * exp(z) + 0 * (y - STEADY_STATE(y));
    end;
    initval; y = 1; z = 0; end;
    steady_state_model; z = 0; y = a; end;
    shocks; var e; stderr 0.1; end;
  "
  expect_error(read_dynare(text = txt), NA)
  m <- read_dynare(text = sub("y = STEADY_STATE(y) * exp(z) + 0 * (y - STEADY_STATE(y));",
                              "log(y) = log(a) + z * STEADY_STATE(y);", txt,
                              fixed = TRUE))
  expect_equal(unname(m$model$ss_links), "y")
  sol <- solve_dsge(m)
  expect_equal(unname(sol$steady_state["y"]), 2, tolerance = 1e-8)
  # d log(y) / d z = STEADY_STATE(y) = a, so dy = y_ss * a * dz
  expect_equal(irf_path(sol, "y", "e", 0L), 2 * 2 * 0.1, tolerance = 1e-6)
  sol3 <- solve_dsge(m, params = c(a = 3, rho = 0.5))
  expect_equal(irf_path(sol3, "y", "e", 0L), 3 * 3 * 0.1, tolerance = 1e-6)
})
