test_that("IRFs stay accurate when the states are linearly dependent", {
  # Kiyotaki and Moore (1997) credit cycle model: land market clearing
  # (k + m * kp = K_bar) makes the states linearly dependent, which gives a
  # transition matrix with very large entries. Forming H^k then lost
  # accuracy to cancellation (0.5% error after four periods).
  f <- tempfile(fileext = ".mod")
  writeLines(c(
    "var x xp b k kp q mu phi;",
    "varexo ed;",
    "parameters alpha m K_bar beta betap a c z;",
    "alpha = 1/3; m = 0.5; K_bar = 1; betap = 0.99; beta = 0.98;",
    "a = 0.7; c = 0.3; z = 0.01;",
    "model;",
    "1 + phi = (beta*(1+phi(+1)) + mu)/betap;",
    "q*(1+phi) + beta*c*phi(+1) = beta*(1+phi(+1))*((1+ed(+1))*(a+c) + q(+1)) + mu*q(+1);",
    "q*(k - k(-1)) + b(-1)/betap + x = (1+ed)*(a+c)*k(-1) + b;",
    "b = betap*q(+1)*k;",
    "q = betap*((1+ed(+1))*alpha*(z + kp)^(alpha-1) + q(+1));",
    "x + m*xp = (1+ed)*(a+c)*k(-1) + m*(1+ed)*(z + kp(-1))^alpha;",
    "k + m*kp = K_bar;",
    "x = c*k(-1);",
    "end;",
    "steady_state_model;",
    "q = a/(1-betap);",
    "kp = (betap*alpha/a)^(1/(1-alpha)) - z;",
    "k = K_bar - m*kp;",
    "b = betap*q*k;",
    "xp = (1/m)*(a*k + m*(z + kp)^alpha);",
    "phi = (a*(beta-1) + beta*c)/(a*(1-beta));",
    "mu = (betap-beta)*beta*c/(a*(1-beta));",
    "x = c*k;",
    "end;",
    "shocks; var ed = 0.0011^2; end;"
  ), f)
  sol <- solve_dsge(suppressWarnings(read_dynare(f)))
  expect_gt(max(abs(sol$H)), 1e5)
  ir <- irf(sol, periods = 8, response = "k", se = FALSE)$data
  # the capital response decays at the rate of the one non-zero eigenvalue
  lambda <- max(Mod(eigen(sol$H, only.values = TRUE)$values))
  ratio <- ir$value[-1] / ir$value[-nrow(ir)]
  expect_equal(ratio[-1], rep(lambda, length(ratio) - 1L), tolerance = 1e-5)
})

test_that("IRF plots hide zero, auxiliary and shock-spike responses", {
  rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
  d <- irf(solve_dsge(rbc), periods = 10)$data
  shown <- dsge:::.dsge_irf_default_responses(d)
  expect_setequal(shown, c("y", "c", "k", "i", "a"))

  pdf(NULL)
  on.exit(grDevices::dev.off())
  m <- dsge_model(
    obs(p ~ beta * lead(p) + kappa * x),
    unobs(x ~ lead(x) - (r - lead(p) - g)),
    obs(r ~ psi * p + u),
    state(u ~ rhou * u),
    state(g ~ rhog * g),
    fixed = list(beta = 0.99),
    start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
  )
  sol <- solve_dsge(m, params = c(kappa = 0.1, psi = 1.5, rhou = 0.7,
                                  rhog = 0.9))
  ir <- irf(sol, periods = 8)
  expect_silent(plot(ir))
  expect_silent(plot(ir, drop_zero = FALSE))
  expect_silent(plot(ir, ncol = 2))
  expect_silent(plot(ir, impulse = "u", col = "black", lwd = 1))
  expect_error(plot(ir, response = "nope"), "Nothing to plot")
})
