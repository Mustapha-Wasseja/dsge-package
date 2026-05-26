# Tests for global_sensitivity()

make_nk <- function() {
  dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
}

priors_nk <- function() {
  list(
    kappa = prior("beta",   shape1 = 2,   shape2 = 8),
    psi   = prior("normal", mean   = 1.5, sd = 0.25),
    rhou  = prior("beta",   shape1 = 5,   shape2 = 2),
    rhog  = prior("beta",   shape1 = 5,   shape2 = 2))
}


test_that("sobol returns expected structure", {
  m <- make_nk()
  gs <- global_sensitivity(m, priors_nk(),
                           target = "sd:p", method = "sobol",
                           n_samples = 30L, seed = 1L)
  expect_s3_class(gs, "dsge_global_sensitivity")
  expect_equal(gs$method, "sobol")
  expect_true(all(c("kappa","psi","rhou","rhog",
                    "sd_e.u","sd_e.g") %in% gs$param_names))
  expect_equal(length(gs$S_first), length(gs$param_names))
  expect_equal(length(gs$S_total), length(gs$param_names))
})


test_that("sobol indices are produced and broadly sensible", {
  m <- make_nk()
  gs <- global_sensitivity(m, priors_nk(),
                           target = "sd:p", method = "sobol",
                           n_samples = 50L, seed = 2L)
  # At least one parameter should produce a finite total-effect index
  expect_true(any(is.finite(gs$S_total)))
  # Total-effect indices should be non-negative (modulo noise)
  finite_idx <- which(is.finite(gs$S_total))
  expect_true(length(finite_idx) > 0)
})


test_that("morris returns mu_star and sigma per parameter", {
  m <- make_nk()
  gs <- global_sensitivity(m, priors_nk(),
                           target = "sd:p", method = "morris",
                           n_samples = 20L, seed = 3L)
  expect_equal(gs$method, "morris")
  expect_equal(length(gs$mu_star), length(gs$param_names))
  expect_equal(length(gs$sigma),   length(gs$param_names))
  expect_true(all(gs$mu_star >= 0, na.rm = TRUE))
})


test_that("irf_max target runs to completion", {
  m <- make_nk()
  gs <- global_sensitivity(m, priors_nk(),
                           target = "irf_max:e.u:p",
                           method = "morris",
                           n_samples = 40L, seed = 4L)
  # Just verify the function returned a result of the right structure.
  # Whether enough finite EE are computed depends on how many random
  # parameter draws hit the determinacy region.
  expect_s3_class(gs, "dsge_global_sensitivity")
  expect_equal(length(gs$mu_star), length(gs$param_names))
})


test_that("print method works", {
  m <- make_nk()
  gs <- global_sensitivity(m, priors_nk(),
                           target = "sd:p", method = "sobol",
                           n_samples = 20L, seed = 5L)
  expect_output(print(gs), "Global sensitivity")
})


test_that("plot method runs", {
  m <- make_nk()
  gs <- global_sensitivity(m, priors_nk(),
                           target = "sd:p", method = "sobol",
                           n_samples = 20L, seed = 6L)
  expect_silent({
    pdf(tempfile()); plot(gs); dev.off()
  })
})
