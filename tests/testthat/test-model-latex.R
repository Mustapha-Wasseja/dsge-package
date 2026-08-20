# Tests for model_latex()

make_nk <- function() {
  dsge_model(
    obs(pi ~ beta * lead(pi) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(pi) - g)),
    obs(r ~ psi * pi + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
}


test_that("model_latex returns one body line per equation", {
  m <- make_nk()
  tex <- model_latex(m)
  # begin + 5 equations + end
  expect_equal(length(tex), 7L)
  expect_match(tex[1], "^\\\\begin\\{align\\}$")
  expect_match(tex[length(tex)], "^\\\\end\\{align\\}$")
})


test_that("greek parameter names are substituted", {
  m <- make_nk()
  tex <- paste(model_latex(m), collapse = "\n")
  expect_match(tex, "\\\\beta", fixed = FALSE)
  expect_match(tex, "\\\\kappa", fixed = FALSE)
  expect_match(tex, "\\\\psi", fixed = FALSE)
})


test_that("greek = FALSE leaves names as plain text", {
  m <- make_nk()
  tex <- paste(model_latex(m, greek = FALSE), collapse = "\n")
  expect_false(grepl("\\\\beta", tex))
  expect_true(grepl("beta", tex, fixed = TRUE))
})


test_that("variables get time subscripts and leads get expectations", {
  m <- make_nk()
  tex <- paste(model_latex(m), collapse = "\n")
  expect_true(grepl("_{t}", tex, fixed = TRUE))
  expect_true(grepl("\\mathbb{E}_{t}", tex, fixed = TRUE))
  expect_true(grepl("_{t+1}", tex, fixed = TRUE))
})


test_that("state equations show t+1 on the LHS and an innovation on the RHS", {
  m <- make_nk()
  tex <- model_latex(m)
  # The two state rows are for u and g
  state_rows <- grep("varepsilon", tex, value = TRUE)
  expect_equal(length(state_rows), 2L)
  expect_true(all(grepl("_\\{t\\+1\\}", state_rows)))
})


test_that("env and numbered arguments control the wrapper", {
  m <- make_nk()
  t_star <- model_latex(m, numbered = FALSE)
  expect_match(t_star[1], "align\\*")

  t_gather <- model_latex(m, env = "gather")
  expect_match(t_gather[1], "^\\\\begin\\{gather\\}$")

  t_none <- model_latex(m, env = "none")
  expect_false(any(grepl("begin\\{", t_none)))
  expect_equal(length(t_none), 5L)   # bodies only
})


test_that("standalone wraps a compilable document", {
  m <- make_nk()
  tex <- model_latex(m, standalone = TRUE)
  expect_match(tex[1], "documentclass")
  expect_true(any(grepl("amsmath", tex)))
  expect_match(tex[length(tex)], "end\\{document\\}")
})


test_that("labels are emitted when requested", {
  m <- make_nk()
  tex <- paste(model_latex(m, labels = TRUE), collapse = "\n")
  expect_true(grepl("\\label{eq:", tex, fixed = TRUE))
})


test_that("file argument writes to disk", {
  m <- make_nk()
  f <- tempfile(fileext = ".tex")
  expect_message(model_latex(m, file = f), "LaTeX written")
  expect_true(file.exists(f))
  content <- readLines(f)
  expect_true(any(grepl("begin\\{align\\}", content)))
})


test_that("nonlinear models render from their equation strings", {
  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C", endo_state = "K", exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9))
  tex <- paste(model_latex(rbc), collapse = "\n")
  # Fractions, exponentials and powers should all be rendered
  expect_true(grepl("\\frac", tex, fixed = TRUE))
  expect_true(grepl("\\exp", tex, fixed = TRUE))
  expect_true(grepl("^{", tex, fixed = TRUE))
  # `beta`, `alpha`, `delta`, `rho` are Greek
  expect_true(grepl("\\beta",  tex, fixed = TRUE))
  expect_true(grepl("\\alpha", tex, fixed = TRUE))
  expect_true(grepl("\\delta", tex, fixed = TRUE))
})


test_that("bar / ss decorations are handled", {
  expect_equal(dsge:::.latex_sym("pi_bar", greek = TRUE),  "\\bar{\\pi}")
  expect_equal(dsge:::.latex_sym("Rk_ss",  greek = TRUE),  "Rk^{ss}")
  expect_match(dsge:::.latex_sym("rho_a",  greek = TRUE),  "\\\\rho_\\{a\\}")
})


test_that("errors on a non-model input", {
  expect_error(model_latex(42), "must be a dsge_model")
})


test_that("no double subscripts are emitted (would break LaTeX compilation)", {
  m <- make_nk()
  tex <- paste(model_latex(m), collapse = "\n")
  # A subscript immediately followed by another subscript,
  # e.g. varepsilon_{u}_{t+1}, is a LaTeX double-subscript
  # error.  Note `}_{` alone is fine (mathbb{E}_{t}).
  expect_false(grepl("_[{][^{}]*[}]_[{]", tex))

  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C", endo_state = "K", exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9))
  expect_false(grepl("_[{][^{}]*[}]_[{]",
                     paste(model_latex(rbc), collapse = "\n")))
})


test_that("state innovations render as a single combined subscript", {
  m <- make_nk()
  tex <- paste(model_latex(m), collapse = "\n")
  expect_true(grepl("\\varepsilon_{u,t+1}", tex, fixed = TRUE))
  expect_true(grepl("\\varepsilon_{g,t+1}", tex, fixed = TRUE))
})


test_that("nonlinear model variables carry time subscripts", {
  rbc <- dsgenl_model(
    "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
    "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
    "Z(+1) = rho * Z",
    observed = "C", endo_state = "K", exo_state = "Z",
    fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
    start = list(rho = 0.9))
  tex <- paste(model_latex(rbc), collapse = "\n")
  # Variables get _{t}
  expect_true(grepl("C_{t}", tex, fixed = TRUE))
  expect_true(grepl("Z_{t}", tex, fixed = TRUE))
  expect_true(grepl("K_{t}", tex, fixed = TRUE))
  # Parameters do NOT get a time subscript
  expect_false(grepl("\\alpha_{t}", tex, fixed = TRUE))
  expect_false(grepl("\\beta_{t}",  tex, fixed = TRUE))
  expect_false(grepl("\\delta_{t}", tex, fixed = TRUE))
})
