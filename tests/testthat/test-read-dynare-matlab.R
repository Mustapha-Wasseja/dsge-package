# MATLAB code in Dynare files, lik_init = 2, higher-order solutions and
# the more robust first-order solvers. Reference values were computed with
# Dynare 6.0 (Octave 8.4).

irf_path <- function(sol, response, impulse, periods = 10L) {
  d <- irf(sol, periods = periods, se = FALSE)$data
  d$value[d$response == response & d$impulse == impulse]
}

test_that("the MATLAB interpreter evaluates common constructs", {
  ctx <- dsge:::mat_new_ctx()
  dsge:::mat_run_script("
    x = [1 2 3; 4 5 6];
    y = x(2, end) + x(end);
    z = [x(1, :) -1];
    s.a = 3; s.b.c = [1, 2]; s.b.c(4) = 7;
    c = {'alpha', 'beta'};
    str = ['k = ' num2str(3) ';'];
    eval(str);
    q = zeros(1, 3);
    for i = 1:3, if i == 2, continue; end, q(i) = i^2; end
    f = @(v) v.^2 - 2;
    r = fzero(f, 1);
    [r2, fv, ef] = fsolve(@(v) [v(1)^2 - 4; v(2) - v(1)], [1; 1]);
    w = -2^2;
    msg = sprintf('%d items: %4.2f', 3, pi);
  ", ctx)
  gv <- function(nm) base::get(nm, envir = ctx$vars)
  expect_equal(as.numeric(gv("y")), 12)
  expect_equal(as.numeric(gv("z")), c(1, 2, 3, -1))
  expect_equal(as.numeric(gv("s")$b$c), c(1, 2, 0, 7))
  expect_equal(gv("c")[[2]], "beta")
  expect_equal(as.numeric(gv("k")), 3)
  expect_equal(as.numeric(gv("q")), c(1, 0, 9))
  expect_equal(as.numeric(gv("r")), sqrt(2), tolerance = 1e-10)
  expect_equal(as.numeric(gv("r2")), c(2, 2), tolerance = 1e-8)
  expect_equal(as.numeric(gv("ef")), 1)
  expect_equal(as.numeric(gv("w")), -4)
  expect_equal(gv("msg"), "3 items: 3.14")
})

test_that("MATLAB statements in a .mod file are run or skipped", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters rho sig;
    Q = [0.5 0
         0.1 0.2];
    V = Q*Q';
    rho = sqrt(V(1,1));
    fprintf('Calibration: rho = %4.2f (%s)\\n', rho, 'ok');
    verbatim;
      tmp = zeros(2, 1);
      for i = 1:2
        tmp(i) = i;
      end
    end;
    set_param_value('sig', tmp(2) / 10);
    model;
      y = rho * y(-1) + sig * e;
    end;
    steady_state_model; y = 0; end;
    shocks; var e; stderr 1; end;
    stoch_simul(order = 1);
    rho = 0.99;
    disp(rho)
  ")
  expect_equal(unname(m$params["rho"]), 0.5)
  expect_equal(unname(m$params["sig"]), 0.2)
  expect_true(any(grepl("after the first computing command", m$notes)))
})

test_that("parameters set in steady_state_model follow the other parameters", {
  txt <- "
    var y k;
    varexo e;
    parameters alpha beta delta kbar;
    alpha = 0.3; beta = 0.99; delta = 0.025; kbar = 0;
    model;
      k = (1 - delta) * k(-1) + delta * kbar * exp(e);
      y = k(-1)^alpha;
    end;
    steady_state_model;
      kbar = ((1/beta - 1 + delta) / alpha)^(1 / (alpha - 1));
      k = kbar;
      y = kbar^alpha;
    end;
    shocks; var e; stderr 0.01; end;
  "
  m <- read_dynare(text = txt)
  kb <- ((1 / 0.99 - 1 + 0.025) / 0.3)^(1 / (0.3 - 1))
  expect_equal(unname(m$params["kbar"]), kb)
  sol <- solve_dsge(m, params = c(alpha = 0.35))
  kb2 <- ((1 / 0.99 - 1 + 0.025) / 0.35)^(1 / (0.35 - 1))
  expect_equal(unname(sol$steady_state["k"]), kb2, tolerance = 1e-8)
})

test_that("a _steadystate.m file is run, including recalibrated parameters", {
  dir <- tempfile("ssm")
  dir.create(dir)
  writeLines("
    var y c k a;
    varexo e;
    parameters alpha beta delta rho sigma_e psi hbar;
    alpha = 0.33; beta = 0.99; delta = 0.025; rho = 0.9; sigma_e = 0.01;
    psi = 0; hbar = 0.3;
    model;
      1/c = beta / c(+1) * (alpha * exp(a(+1)) * k^(alpha - 1) + 1 - delta);
      y = exp(a) * k(-1)^alpha;
      k = y - c + (1 - delta) * k(-1) + psi * 0;
      a = rho * a(-1) + e;
    end;
    shocks; var e; stderr sigma_e; end;
  ", file.path(dir, "rbcm.mod"))
  writeLines("
function [ys, params, check] = rbcm_steadystate(ys, exo, M_, options_)
  NumberOfParameters = M_.param_nbr;
  for ii = 1:NumberOfParameters
    paramname = M_.param_names{ii};
    eval([paramname ' = M_.params(' int2str(ii) ');']);
  end
  check = 0;
  % solve the capital Euler equation numerically
  [k, fval, exitflag] = fsolve(@(k) alpha * k^(alpha - 1) + 1 - delta - 1/beta, 20, optimset('Display', 'off'));
  if exitflag < 1
    check = 1;
    return
  end
  psi = get_psi(hbar);   % a recalibrated parameter
  y = k^alpha;
  c = y - delta * k;
  a = 0;
  params = NaN(NumberOfParameters, 1);
  for iter = 1:length(M_.params)
    eval(['params(' num2str(iter) ') = ' M_.param_names{iter} ';'])
  end
  for ii = 1:M_.orig_endo_nbr
    varname = M_.endo_names{ii};
    eval(['ys(' int2str(ii) ') = ' varname ';']);
  end
end

function psi = get_psi(h)
  psi = h / (1 - h);
end
", file.path(dir, "rbcm_steadystate.m"))
  m <- read_dynare(file.path(dir, "rbcm.mod"))
  expect_true(any(grepl("rbcm_steadystate.m", m$notes)))
  expect_equal(unname(m$params["psi"]), 0.3 / 0.7)
  sol <- solve_dsge(m)
  kss <- ((1 / 0.99 - 1 + 0.025) / 0.33)^(1 / (0.33 - 1))
  expect_equal(unname(sol$steady_state["k"]), kss, tolerance = 1e-8)
  sol2 <- solve_dsge(m, params = c(beta = 0.98))
  kss2 <- ((1 / 0.98 - 1 + 0.025) / 0.33)^(1 / (0.33 - 1))
  expect_equal(unname(sol2$steady_state["k"]), kss2, tolerance = 1e-8)
})

test_that("lik_init = 2 reproduces Dynare's likelihood", {
  mod <- "
    var y pi r g u v;
    varexo eg eu ev;
    parameters beta sigma kappa phi_pi phi_y rho_r rho_g rho_u rho_v;
    beta = 0.99; sigma = 1; kappa = 0.1;
    phi_pi = 1.5; phi_y = 0.125; rho_r = 0.8;
    rho_g = 0.9; rho_u = 0.7; rho_v = 0.5;
    model(linear);
      # rr = r - pi(+1);
      y = y(+1) - 1/sigma * rr + g;
      pi = beta * pi(+1) + kappa * y + u;
      r = rho_r * r(-1) + (1 - rho_r) * (phi_pi * pi + phi_y * y) + v;
      g = rho_g * g(-1) + eg;
      u = rho_u * u(-1) + eu;
      v = rho_v * v(-1) + ev;
    end;
    shocks;
      var eg; stderr 0.01;
      var eu = 0.005^2;
      var ev; stderr 0.0025;
    end;
    shocks;
      var r; stderr 0.002;
    end;
    varobs y pi r;
    estimated_params;
      rho_g, 0.9;
    end;
    estimation(datafile = 'data.csv', mode_compute = 0, lik_init = 2);
  "
  m <- read_dynare(text = mod)
  expect_equal(m$model$kalman_init$type, "lik_init_2")
  set.seed(5)
  dat <- data.frame(y = cumsum(stats::rnorm(60)) * 0.01,
                    pi = stats::rnorm(60) * 0.005,
                    r = stats::rnorm(60) * 0.004)
  dat$y <- dat$y - mean(dat$y)
  sol <- solve_dsge(m)
  d <- dsge:::dyn_map_data(m, dat)
  y <- as.matrix(d[, m$model$variables$observed])
  ll <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D,
                             init = m$model$kalman_init)$loglik
  # Dynare 6.0: 508.4024094091 (it switches to the steady-state Kalman gain
  # once the covariance has converged, hence the small difference)
  expect_equal(ll, 508.4024094091, tolerance = 1e-8)
  ll1 <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D)$loglik
  expect_equal(ll1, 539.544231, tolerance = 1e-8)
})

test_that("second- and third-order solutions match Dynare", {
  m <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  pt <- function(sol, order) {
    x <- stats::setNames(numeric(ncol(sol$H)), colnames(sol$H))
    x[c("k_lag1", "a_lag1", "e")] <- c(0.3, 0.01, 0.004)
    y <- sol$G %*% x + 0.5 * sol$g_ss
    for (i in seq_len(nrow(sol$G))) {
      y[i] <- y[i] + 0.5 * sum(sol$g_xx[i, , ] * (x %o% x))
      if (order == 3L) {
        y[i] <- y[i] + sum(sol$g_xxx[i, , , ] * (x %o% x %o% x)) / 6 +
          0.5 * sum(sol$g_xss[i, ] * x) + sol$g_sss[i] / 6
      }
    }
    stats::setNames(as.numeric(y) + sol$steady_state[rownames(sol$G)],
                    rownames(sol$G))
  }
  s2 <- solve_dsge(m, order = 2)
  s3 <- solve_dsge(m, order = 3)
  # Dynare: 0.5 * ghs2 and the decision rules at (k, a, e) deviations
  # (0.3, 0.01, 0.004)
  expect_equal(unname(0.5 * s2$g_ss["c"]) / -1.92611860662934e-05, 1,
               tolerance = 1e-6)
  expect_equal(unname(pt(s2, 2L)[c("y", "c", "k", "i")]),
               c(3.06694453470985, 2.33111640959941, 28.6680367096327,
                 0.73582812511044), tolerance = 1e-10)
  expect_equal(unname(pt(s3, 3L)[c("y", "c", "k", "i")]),
               c(3.06694644669871, 2.33111689572592, 28.668038135495,
                 0.735829550972783), tolerance = 1e-10)
  expect_equal(unname(s3$g_sss), numeric(length(s3$g_sss)))
})

test_that("leads of two periods in nonlinear terms become auxiliary variables", {
  r <- dsge:::dyn_substitute_leads(
    c("0 = e/(p(+1)*c(+1)) - beta*e(+1)*(1+rf)/(p(+2)*c(+2))",
      "y = x(+2) + 2*x(+3)"),
    c("e", "p", "c", "rf", "x", "y"), "u")
  expect_equal(r$equations[1], "0 = e/(p(+1) * c(+1)) - aux_lead_1(+1)")
  expect_equal(r$equations[2], "y = x(+2) + 2 * x(+3)")
  expect_equal(r$defs$aux_lead_1, "beta * e * (1 + rf(-1))/(p(+1) * c(+1))")
})

test_that("first-order solver handles unit roots, redundant states and pegs", {
  # A unit root (price level) and a forward-looking NK model closed by the
  # discretionary targeting rule (an instrument peg would be indeterminate)
  m <- read_dynare(text = "
    var pi x i p u;
    varexo eu;
    parameters beta kappa sigma rho_u lambda;
    beta = 0.99; kappa = 0.1; sigma = 1; rho_u = 0.5; lambda = 0.25;
    model(linear);
      pi = beta * pi(+1) + kappa * x + u;
      x = x(+1) - 1/sigma * (i - pi(+1));
      pi = p - p(-1);
      u = rho_u * u(-1) + eu;
    end;
    shocks; var eu; stderr 1; end;
    planner_objective pi^2 + lambda * x^2;
    discretionary_policy(instruments = (i), planner_discount = beta);
  ")
  sol <- solve_dsge(m)
  expect_true(sol$stable)
  # Gali (2015, ch. 5): x = -kappa / (kappa^2 + lambda (1 - beta rho_u)) u
  coef <- -0.1 / (0.1^2 + 0.25 * (1 - 0.99 * 0.5))
  expect_equal(irf_path(sol, "x", "eu", 0L), coef, tolerance = 1e-8)
})
