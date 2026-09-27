# Models without stochastic shocks (no varexo), e.g. deterministic growth
# models written for perfect-foresight simulation.

solow_no_shocks <- function() {
  f <- tempfile(fileext = ".mod")
  writeLines(c(
    "var k c;",
    "parameters alpha delta s;",
    "alpha = 0.33; delta = 0.1; s = 0.2;",
    "model;",
    "k = (1 - delta) * k(-1) + s * k(-1)^alpha;",
    "c = (1 - s) * k(-1)^alpha;",
    "end;",
    "initval;",
    "k = 1; c = 1;",
    "end;"
  ), f)
  suppressWarnings(read_dynare(f))
}

test_that("solve_dsge solves a model without shocks", {
  m <- solow_no_shocks()
  s <- solve_dsge(m)
  expect_true(s$stable)
  expect_equal(ncol(s$M), 0L)
  expect_equal(nrow(s$M), nrow(s$H))
})

test_that("irf and variance_decomposition explain that there are no shocks", {
  s <- solve_dsge(solow_no_shocks())
  expect_error(irf(s), "no stochastic shocks|has none")
  expect_error(variance_decomposition(s), "simulate_perfect_foresight")
})
