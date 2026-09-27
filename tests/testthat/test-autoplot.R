# ggplot2 versions of the plots (ggplot2 is optional)

skip_if_not_installed("ggplot2")

build_ok <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
}

nk_solution <- function() {
  m <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )
  list(model = m,
       sol = solve_dsge(m, params = c(kappa = 0.1, psi = 1.5, rhou = 0.7,
                                      rhog = 0.9),
                        shock_sd = c(e.u = 1, e.g = 0.5)))
}

test_that("autoplot methods are registered with ggplot2", {
  for (cl in c("dsge_irf", "dsge_forecast", "dsge_variance_decomposition",
               "dsge_decomposition", "dsge_smoothed", "dsge_bayes")) {
    expect_false(is.null(utils::getS3method("autoplot", cl,
                                            optional = TRUE,
                                            envir = asNamespace("ggplot2"))),
                 info = cl)
  }
})

test_that("autoplot draws IRFs and variance decompositions", {
  nk <- nk_solution()
  ir <- irf(nk$sol, periods = 8)
  build_ok(ggplot2::autoplot(ir))
  build_ok(ggplot2::autoplot(ir, impulse = "e.u"))
  build_ok(ggplot2::autoplot(variance_decomposition(nk$sol,
                                                    horizon = c(1, 4))))
  build_ok(ggplot2::autoplot(variance_decomposition(nk$sol)))
  expect_error(ggplot2::autoplot(ir, response = "nope"), "Nothing to plot")
})

test_that("autoplot draws fit-based results", {
  nk <- nk_solution()
  sol <- nk$sol
  set.seed(1)
  TT <- 60
  xs <- matrix(0, TT, 2)
  y <- matrix(0, TT, nrow(sol$G))
  for (t in 2:TT) {
    xs[t, ] <- sol$H %*% xs[t - 1, ] + sol$M %*% (rnorm(2) * c(1, 0.5))
    y[t, ] <- sol$G %*% xs[t, ]
  }
  colnames(y) <- rownames(sol$G)
  fit <- estimate(nk$model, data = as.data.frame(y[, c("p", "r")]))
  build_ok(ggplot2::autoplot(irf(fit, periods = 8)))
  build_ok(ggplot2::autoplot(forecast(fit, horizon = 6)))
  build_ok(ggplot2::autoplot(smooth_states(fit)))
  build_ok(ggplot2::autoplot(shock_decomposition(fit)))
})

test_that("autoplot draws Dynare-model, OccBin and Bayesian results", {
  rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  s1 <- solve_dsge(rbc)
  p <- ggplot2::autoplot(irf(s1, periods = 10))
  build_ok(p)
  expect_setequal(unique(as.character(p$data$response)),
                  c("y", "c", "k", "i", "a"))
  build_ok(ggplot2::autoplot(irf_2nd_order(solve_dsge(rbc, order = 2),
                                           "e", 0.05, periods = 10)))
  build_ok(ggplot2::autoplot(perfect_foresight(s1, shocks = list(e = 0.01),
                                               horizon = 10)))

  nk <- dsge_model(
    obs(pi ~ beta * lead(pi) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(pi) - g)),
    obs(r ~ psi * pi + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99, kappa = 0.1, psi = 1.5),
    start = list(rhou = 0.5, rhog = 0.5)
  )
  s2 <- solve_dsge(nk, params = list(rhou = 0.5, rhog = 0.5),
                   shock_sd = c(u = 0.5, g = 0.5))
  build_ok(ggplot2::autoplot(simulate_occbin(s2, constraints = list("r >= 0"),
                                             shocks = list(g = -0.05),
                                             horizon = 20)))

  skip_on_cran()
  m <- dsge_model(obs(y ~ z), state(z ~ rho * z), start = list(rho = 0.5))
  set.seed(42)
  z <- numeric(100)
  for (i in 2:100) z[i] <- 0.8 * z[i - 1] + rnorm(1)
  bfit <- bayes_dsge(m, data = data.frame(y = z),
                     priors = list(rho = prior("beta", shape1 = 2,
                                               shape2 = 2)),
                     chains = 2, iter = 400, seed = 1)
  build_ok(ggplot2::autoplot(bfit))
  build_ok(ggplot2::autoplot(bfit, type = "density"))
  expect_error(ggplot2::autoplot(bfit, pars = "nope"), "Unknown parameter")
})

test_that("theme and scales are ggplot2 objects", {
  expect_s3_class(theme_dsge(), "theme")
  expect_true(inherits(scale_colour_dsge(), "Scale"))
  expect_true(inherits(scale_fill_dsge(), "Scale"))
})
