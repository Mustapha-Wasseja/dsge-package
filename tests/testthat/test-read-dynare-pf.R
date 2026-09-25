# Perfect-foresight simulation of imported models and MATLAB data/optimisation
# functions. Reference paths were computed with Dynare 6.0 (Octave 8.4),
# perfect_foresight_solver with tolf = tolx = 1e-12.

test_that("a permanent shock reproduces Dynare's perfect-foresight path", {
  m <- read_dynare(text = "
    var c k;
    varexo a;
    parameters alpha beta delta;
    alpha = 0.33; beta = 0.99; delta = 0.025;
    model;
      1/c = beta/c(+1) * (alpha * exp(a(+1)) * k^(alpha - 1) + 1 - delta);
      k = exp(a) * k(-1)^alpha + (1 - delta) * k(-1) - c;
    end;
    steady_state_model;
      k = ((1/beta - 1 + delta)/(alpha*exp(a)))^(1/(alpha - 1));
      c = exp(a)*k^alpha - delta*k;
    end;
    initval; a = 0; k = 30; c = 2.3; end;
    steady;
    endval; a = 0.1; end;
    steady;
    perfect_foresight_setup(periods = 100);
    perfect_foresight_solver;
  ")
  pf <- simulate_perfect_foresight(m)
  expect_true(pf$converged)
  expect_equal(pf$periods, 100L)
  rows <- as.character(c(0, 1, 2, 10, 50, 101))
  expect_equal(unname(pf$path[rows, "c"]),
               c(2.306617231988, 2.453436260101, 2.462025442352,
                 2.519845828378, 2.644784899034, 2.677907700536),
               tolerance = 1e-9)
  expect_equal(unname(pf$path[rows, "k"]),
               c(28.348419061048, 28.518724816338, 28.682777099393,
                 29.795089833737, 32.232999526409, 32.911593934550),
               tolerance = 1e-9)
  # a different horizon and parameters can be given
  pf2 <- simulate_perfect_foresight(m, periods = 50, params = c(delta = 0.03))
  expect_equal(nrow(pf2$path), 52L)
})

test_that("mcp tags give a zero lower bound (Dynare's lmmcp)", {
  txt <- "
    var pi x i;
    varexo r_nat;
    parameters beta kappa sigma phi;
    beta = 0.99; kappa = 0.1; sigma = 1; phi = 1.5;
    model;
      pi = beta*pi(+1) + kappa*x;
      x = x(+1) - 1/sigma*(i - pi(+1) - r_nat);
      [name = 'Taylor rule', mcp = 'i > 0']
      i = 1 + phi*pi + 0.5*x;
    end;
    initval; r_nat = 1; i = 1; end;
    steady;
    shocks;
      var r_nat; periods 1:6; values -1;
    end;
    perfect_foresight_setup(periods = 40);
    perfect_foresight_solver(lmmcp);
  "
  pf <- simulate_perfect_foresight(read_dynare(text = txt))
  expect_true(pf$converged)
  expect_equal(unname(pf$path[1:7, "i"]),
               c(0, 0, 0, 0, 0, 0.212121212121, 1), tolerance = 1e-9)
  expect_equal(unname(pf$path[1:6, "pi"]),
               c(-3.051110435533, -2.004970246667, -1.251795333333,
                 -0.718466666667, -0.353333333333, -0.121212121212),
               tolerance = 1e-9)
  # without lmmcp the rule holds and the rate goes negative
  pf0 <- simulate_perfect_foresight(read_dynare(text = txt), lmmcp = FALSE)
  expect_lt(min(pf0$path[, "i"]), 0)
})

test_that("histval and oo_.endo_simul set the initial condition", {
  m <- read_dynare(text = "
    var y;
    varexo e;
    parameters rho;
    rho = 0.5;
    model;
      y = rho * y(-1) + e;
    end;
    initval; y = 0; e = 0; end;
    histval; y(0) = 1; end;
    perfect_foresight_setup(periods = 5);
    perfect_foresight_solver;
  ")
  pf <- simulate_perfect_foresight(m)
  # no leads, so no terminal period
  expect_equal(unname(pf$path[, "y"]), c(1, 0.5^(1:5)), tolerance = 1e-12)
  m2 <- read_dynare(text = "
    var y;
    varexo e;
    parameters rho;
    rho = 0.5;
    model;
      y = rho * y(-1) + e;
    end;
    initval; y = 0; e = 0; end;
    perfect_foresight_setup(periods = 5);
    oo_.endo_simul(:, 1) = 2;
    perfect_foresight_solver;
  ")
  expect_equal(unname(simulate_perfect_foresight(m2)$path[2, "y"]), 1)
})

test_that("MATLAB code can read data files and optimise", {
  dir <- tempfile("matdata")
  dir.create(dir)
  writeLines(c("# Created by Octave 8.4.0", "# name: V", "# type: matrix",
               "# rows: 2", "# columns: 2", " 4 1", " 1 9", "", "",
               "# name: s", "# type: scalar", "0.5", "", ""),
             file.path(dir, "cal.mat"))
  writeLines(c("1 2", "3 4"), file.path(dir, "nums.txt"))
  ctx <- dsge:::mat_new_ctx()
  ctx$path <- dir
  dsge:::mat_run_script("
    load cal
    S = load('cal.mat');
    load('nums.txt');
    [x, fv, ef] = fmincon(@(z) (z(1)-3)^2 + (z(2)-1)^2, [0; 0], [1 1], 2, ...
                          [], [], [0; 0], [5; 5]);
    [tr, cy] = hpfilter((1:20)'.^1.5, 1600);
    yy = interp1([1 2 3], [10 20 30], 2.5);
    p = polyfit([1 2 3], [2 4 6], 1);
  ", ctx)
  gv <- function(nm) base::get(nm, envir = ctx$vars)
  expect_equal(gv("V"), matrix(c(4, 1, 1, 9), 2))
  expect_equal(as.numeric(gv("S")$s), 0.5)
  expect_equal(gv("nums"), matrix(c(1, 3, 2, 4), 2))
  expect_equal(as.numeric(gv("x")), c(2, 0), tolerance = 1e-6)
  expect_equal(as.numeric(gv("ef")), 1)
  expect_equal(as.numeric(gv("tr") + gv("cy")), (1:20)^1.5)
  expect_equal(as.numeric(gv("yy")), 25)
  expect_equal(as.numeric(gv("p")), c(2, 0), tolerance = 1e-10)

  skip_if_not_installed("R.matlab")
  R.matlab::writeMat(file.path(dir, "bin.mat"),
                     cov_matrix = matrix(c(2, 0.5, 0.5, 3), 2))
  dsge:::mat_run_script("B = load('bin.mat'); sx = sqrt(B.cov_matrix(2, 2));",
                        ctx)
  expect_equal(as.numeric(gv("sx")), sqrt(3))
})

test_that("xlsread reads spreadsheets (readxl)", {
  skip_if_not_installed("readxl")
  ctx <- dsge:::mat_new_ctx()
  ctx$path <- dirname(readxl::readxl_example("datasets.xlsx"))
  dsge:::mat_run_script("d = xlsread('datasets.xlsx', 1);", ctx)
  d <- base::get("d", envir = ctx$vars)
  expect_equal(dim(d), c(32L, 11L))
  expect_equal(d[1, 1], 21)
})
