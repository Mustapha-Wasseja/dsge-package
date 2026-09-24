# Optimal policy and OccBin in imported Dynare models.
# Reference values were computed with Dynare 6.0 (see dev/dynare-validation).

irf_path <- function(sol, response, impulse, periods = 10L) {
  d <- irf(sol, periods = periods, se = FALSE)$data
  d$value[d$response == response & d$impulse == impulse]
}

nk_block <- "
  var y pi r u g;
  varexo eu eg;
  parameters beta sigma kappa lambda rho_u rho_g;
  beta = 0.99; sigma = 1; kappa = 0.1; lambda = 0.25;
  rho_u = 0.5; rho_g = 0.8;
  model(linear);
    y = y(+1) - 1/sigma * (r - pi(+1)) + g;
    pi = beta * pi(+1) + kappa * y + u;
    u = rho_u * u(-1) + eu;
    g = rho_g * g(-1) + eg;
  end;
  shocks;
    var eu; stderr 0.01;
    var eg; stderr 0.01;
  end;
  planner_objective pi^2 + lambda * y^2;
"


test_that("ramsey_model adds multipliers and matches Dynare", {
  m <- read_dynare(text = paste(nk_block,
    "ramsey_model(planner_discount = beta, instruments = (r));"))
  expect_equal(m$policy$type, "ramsey")
  expect_equal(m$policy$multipliers, c("MULT_1", "MULT_2", "MULT_3",
                                       "MULT_4"))
  expect_true(all(m$policy$multipliers %in% m$model$controls))
  sol <- solve_dsge(m)
  expect_true(sol$stable)
  expect_equal(irf_path(sol, "pi", "eu", 2L),
               c(0.0138780618572, 0.00447796397123, 0.000214348489257),
               tolerance = 1e-8)
  expect_equal(irf_path(sol, "y", "eu", 2L),
               c(-0.00555122474289, -0.00734241033138, -0.00742814972708),
               tolerance = 1e-8)
})


test_that("discretionary_policy closes the model with the Dennis rule", {
  m <- read_dynare(text = paste(nk_block,
    "discretionary_policy(instruments = (r), planner_discount = beta);"))
  expect_equal(m$policy$type, "discretion")
  expect_equal(rownames(m$policy$rule$F1), "r")
  sol <- solve_dsge(m)
  expect_equal(irf_path(sol, "pi", "eu", 2L),
               c(0.0183486223735, 0.00917431118674, 0.00458715559337),
               tolerance = 1e-6)
  expect_equal(irf_path(sol, "y", "eu", 2L),
               c(-0.00733944894939, -0.00366972447469, -0.00183486223735),
               tolerance = 1e-6)
  # Targeting rule under discretion: kappa * pi + lambda * y = 0
  expect_equal(0.1 * irf_path(sol, "pi", "eu", 3L) +
                 0.25 * irf_path(sol, "y", "eu", 3L), rep(0, 4),
               tolerance = 1e-10)
})


test_that("a nonlinear Ramsey problem solves for its steady state", {
  m <- read_dynare(text = "
    var c k a;
    varexo e;
    parameters alpha delta rho beta;
    alpha = 0.33; delta = 0.025; rho = 0.9; beta = 0.99;
    model;
      k = exp(a) * k(-1)^alpha + (1 - delta) * k(-1) - c;
      a = rho * a(-1) + e;
    end;
    initval; a = 0; k = 28; c = 2.3; end;
    shocks; var e; stderr 0.01; end;
    planner_objective log(c);
    ramsey_model(planner_discount = beta, instruments = (c));
  ")
  sol <- solve_dsge(m)
  k_ss <- ((1 / 0.99 - 1 + 0.025) / 0.33)^(1 / (0.33 - 1))
  expect_equal(unname(sol$steady_state["k"]), k_ss, tolerance = 1e-8)
  # Multiplier on the resource constraint equals marginal utility 1/c
  expect_equal(abs(unname(sol$steady_state["MULT_1"])),
               1 / unname(sol$steady_state["c"]), tolerance = 1e-8)
})


test_that("osr() runs the OSR problem declared in the file", {
  m <- read_dynare(text = "
    var y pi r u g;
    varexo eu eg;
    parameters beta sigma kappa rho_u rho_g phi_pi phi_y;
    beta = 0.99; sigma = 1; kappa = 0.1; rho_u = 0.5; rho_g = 0.8;
    phi_pi = 1.5; phi_y = 0.5;
    model(linear);
      y = y(+1) - 1/sigma * (r - pi(+1)) + g;
      pi = beta * pi(+1) + kappa * y + u;
      r = phi_pi * pi + phi_y * y;
      u = rho_u * u(-1) + eu;
      g = rho_g * g(-1) + eg;
    end;
    shocks; var eu; stderr 0.01; var eg; stderr 0.01; end;
    osr_params phi_pi phi_y;
    osr_params_bounds; phi_pi, 1.01, 5; phi_y, 0, 2; end;
    optim_weights; pi 1; y 0.25; r 0.1; end;
    osr;
  ")
  expect_equal(m$policy$type, "osr")
  expect_equal(unname(m$policy$weights["y", "y"]), 0.25)
  o <- osr(m)
  expect_equal(o$loss_at_start, 0.000654717, tolerance = 1e-6)
  expect_equal(o$loss, 0.000529050096519, tolerance = 1e-8)
  expect_equal(o$optimal[["phi_pi"]], 2.3362157159, tolerance = 1e-4)
  expect_equal(o$optimal[["phi_y"]], 2, tolerance = 1e-6)

  # Supplied values for fixed parameters are honoured by solve_dsge()
  s1 <- solve_dsge(m)
  s2 <- solve_dsge(m, params = c(phi_pi = 3))
  expect_false(isTRUE(all.equal(s1$G, s2$G)))
})


test_that("occbin_constraints are simulated piecewise-linearly as in Dynare", {
  m <- read_dynare(text = "
    var y pi i inot g;
    varexo eg;
    parameters beta sigma kappa phi_pi phi_y rho_g ilb;
    beta = 0.99; sigma = 1; kappa = 0.1; phi_pi = 1.5; phi_y = 0.125;
    rho_g = 0.8; ilb = -0.01;
    model(linear);
      y = y(+1) - 1/sigma * (i - pi(+1)) + g;
      pi = beta * pi(+1) + kappa * y;
      inot = phi_pi * pi + phi_y * y;
      [name = 'policy', relax = 'zlb']
      i = inot;
      [name = 'policy', bind = 'zlb']
      i = ilb;
      g = rho_g * g(-1) + eg;
    end;
    occbin_constraints;
      name 'zlb'; bind inot <= ilb; relax inot > ilb;
    end;
    shocks(surprise);
      var eg; periods 1; values -0.06;
    end;
  ")
  expect_false(is.null(m$occbin))
  expect_equal(names(m$occbin$constraints), "zlb")
  o <- simulate_occbin(m, horizon = 30)
  expect_s3_class(o, "dsge_occbin")
  expect_true(o$converged)
  expect_equal(which(o$binding[, 1]), 1:10)
  expect_equal(unname(o$controls[1:3, "y"]),
               c(-0.496464036478405, -0.342512400952604, -0.234107961123776),
               tolerance = 1e-8)
  expect_true(all(o$controls[o$binding[, 1], "i"] + 0.01 < 1e-10))
  expect_true(any(o$controls_unc[, "i"] < -0.01))
  expect_output(print(o), "inot >= ilb")
})
