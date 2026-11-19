# Paccagnini, *Macroeconometrics*, Chapter 6 -- empirical example
# Bayesian estimation of Smets and Wouters (2007, AER) with dsge.
#
# The book runs this example in Dynare (Chapter6/Codes_Empirical_Examples):
#   sw_table1.mod  (@#include sw_defs.mod, sw_model_def.mod,
#                   sw_estimation_def.mod)
#   usmodel_data.mat   (US data, 1947Q3-2004Q4)
#   table1.m           (Table 1A / 1B: posterior mode, mean, 90% HPD)
# Here the same .mod files are imported unchanged with read_dynare() and
# estimated in R, without Dynare or MATLAB. The estimation options in the
# file are honoured: first_obs = 71, presample = 4, lik_init = 2,
# prefilter = 0.
#
# Steps
#   1. Posterior mode within the estimated_params bounds (the analogue of
#      Dynare's mode_compute = 1).
#   2. Random-walk Metropolis-Hastings with bayes_dsge(), started at the
#      mode. Dynare's file asks for 2 blocks x 2000 draws and itself warns
#      that this is too short. Random-walk MH mixes slowly on this
#      posterior (autocorrelation times of several hundred draws), so the
#      default here is 4 chains x 150000 iterations (50000 warm-up, every
#      5th draw kept), settable via SW_CHAINS, SW_ITER, SW_WARMUP, SW_THIN,
#      SW_CORES. Per-chain means and log-posterior levels are printed so
#      that a chain stuck in a minor mode is visible.
#   3. Table 1A and Table 1B (as table1.m) and Figure 6.2 (prior vs
#      posterior) from the posterior draws.
#
# Run from the package root:
#   Rscript dev/macroeconometrics-book/ch6_smets_wouters_estimation.R
# Requires the R.matlab package to read usmodel_data.mat.

suppressPackageStartupMessages(library(dsge))
local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
  source(file.path(here, "book_path.R"))
})

env_int <- function(name, default) as.integer(Sys.getenv(name, default))
n_chains <- env_int("SW_CHAINS", 4L)
n_iter   <- env_int("SW_ITER", 150000L)
n_warmup <- env_int("SW_WARMUP", 50000L)
n_thin   <- env_int("SW_THIN", 5L)
n_cores  <- env_int("SW_CORES", min(n_chains, parallel::detectCores()))
seed     <- env_int("SW_SEED", 2007L)
out      <- out_dir()

# ---------------------------------------------------------------------------
# 0. Import the book's Dynare files and data
# ---------------------------------------------------------------------------

ex_dir <- file.path(book_dir(), "Chapter6", "Codes_Empirical_Examples")
m <- read_dynare(file.path(ex_dir, "sw_table1.mod"))
print(m)
cat("\nEstimation options read from the file:\n")
str(m$estimation[c("first_obs", "presample", "lik_init", "datafile")])

raw <- R.matlab::readMat(file.path(ex_dir, "usmodel_data.mat"))
obs <- m$model$variables$observed
dat <- as.data.frame(lapply(stats::setNames(obs, obs),
                            function(v) as.numeric(raw[[v]])))
first <- m$estimation$first_obs
quarter <- function(i) {           # usmodel_data.mat starts in 1947Q3
  q <- 1947 * 4 + 2 + (i - 1)
  sprintf("%dQ%d", q %/% 4, q %% 4 + 1)
}
cat(sprintf("\nData: %d quarters (%s-%s); estimation sample from %s, ",
            nrow(dat), quarter(1), quarter(nrow(dat)), quarter(first)),
    sprintf("likelihood from %s (presample %d)\n",
            quarter(first + m$estimation$presample), m$estimation$presample),
    sep = "")

# ---------------------------------------------------------------------------
# 1. Posterior kernel and mode
#    log p(Y | theta) + log p(theta) in the natural parameterisation, as
#    Dynare maximises it, with theta kept within the estimated_params
#    bounds. The Kalman filter and prior densities are dsge's own.
# ---------------------------------------------------------------------------

est   <- m$estimated_params
key   <- ifelse(est$type == "stderr", paste0("sd_e.", est$name), est$name)
stopifnot(setequal(key, names(m$priors)))
lower <- stats::setNames(est$lower, key)
upper <- stats::setNames(est$upper, key)
init  <- stats::setNames(est$init, key)
is_sd <- startsWith(key, "sd_e.")
fixed <- unlist(m$model$fixed)
Y     <- as.matrix(dat[first:nrow(dat), obs])

split_theta <- function(theta) {
  names(theta) <- key
  list(params = c(theta[!is_sd], fixed),
       sd = stats::setNames(theta[is_sd], sub("sd_e.", "", key[is_sd], fixed = TRUE)))
}

log_lik <- function(theta) {
  p <- split_theta(theta)
  sol <- tryCatch(solve_dsge(m$model, params = p$params, shock_sd = p$sd),
                  error = function(e) NULL)
  if (is.null(sol) || !isTRUE(sol$stable)) return(-Inf)
  dsge:::kalman_filter(sweep(Y, 2, sol$steady_state[obs]), sol$G, sol$H,
                       sol$M, sol$D, presample = m$estimation$presample,
                       init = m$model$kalman_init, loglik_only = TRUE)$loglik
}
log_prior <- function(theta) {
  names(theta) <- key
  sum(vapply(key, function(nm) dsge:::dprior(m$priors[[nm]], theta[[nm]]), 0))
}
log_post <- function(theta) {
  lp <- log_prior(theta)
  if (!is.finite(lp)) return(-Inf)
  ll <- log_lik(theta)
  if (!is.finite(ll)) -Inf else ll + lp
}

cat(sprintf("\nAt the file's initial values: log-likelihood %.4f, log-prior %.4f\n",
            log_lik(init), log_prior(init)))

mode_file <- file.path(out, "sw_posterior_mode.rds")
if (file.exists(mode_file) && !nzchar(Sys.getenv("SW_REFIT_MODE"))) {
  mode <- readRDS(mode_file)
  cat("Posterior mode read from", mode_file, "\n")
} else {
  cat("\nSearching for the posterior mode (L-BFGS-B within the bounds)...\n")
  t0 <- Sys.time()
  obj <- function(x) { v <- log_post(x); if (is.finite(v)) -v else 1e10 }
  scale <- pmin(stats::setNames(est$prior_sd, key), upper - lower)
  # Central-difference gradient, coordinates evaluated in parallel
  grad <- function(x) {
    h <- 1e-5 * scale
    unlist(parallel::mclapply(seq_along(x), function(j) {
      up <- min(x[j] + h[j], upper[j]); dn <- max(x[j] - h[j], lower[j])
      xu <- x; xu[j] <- up; xd <- x; xd[j] <- dn
      (obj(xu) - obj(xd)) / (up - dn)
    }, mc.cores = n_cores))
  }
  x <- init
  for (round in 1:5) {                       # restart until no improvement
    o <- stats::optim(x, obj, grad, method = "L-BFGS-B", lower = lower,
                      upper = upper, control = list(maxit = 1000,
                      parscale = scale, factr = 1e7))
    improved <- obj(x) - o$value
    x <- o$par
    cat(sprintf("  round %d: log posterior %.4f (improved by %.4f)\n",
                round, -o$value, improved))
    if (improved < 1e-3) break
  }
  mode <- stats::setNames(x, key)
  saveRDS(mode, mode_file)
  cat(sprintf("  done in %.1f minutes\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
mode_ll <- log_lik(mode); mode_lp <- log_prior(mode)
cat(sprintf("At the mode: log-likelihood %.4f, log-prior %.4f, log posterior %.4f\n",
            mode_ll, mode_lp, mode_ll + mode_lp))
at_bound <- key[abs(mode - lower) < 1e-6 | abs(mode - upper) < 1e-6]
if (length(at_bound))
  cat("Parameters at a bound at the mode:", paste(at_bound, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 2. MCMC with bayes_dsge(), started at the mode
# ---------------------------------------------------------------------------

fit_file <- file.path(out, "sw_fit.rds")
if (file.exists(fit_file) && !nzchar(Sys.getenv("SW_REFIT_MCMC"))) {
  fit <- readRDS(fit_file)
  cat("\nPosterior draws read from", fit_file, "\n")
} else {
  m_start <- m
  m_start$model$start[key[!is_sd]] <- as.list(mode[!is_sd])
  sd_start <- split_theta(mode)$sd
  cat(sprintf("\nbayes_dsge(): %d chains x %d iterations (%d warm-up, thin %d) on %d cores...\n",
              n_chains, n_iter, n_warmup, n_thin, n_cores))
  t0 <- Sys.time()
  fit <- bayes_dsge(m_start, data = dat, chains = n_chains, iter = n_iter,
                    warmup = n_warmup, thin = n_thin, n_cores = n_cores,
                    seed = seed, shock_start = sd_start)
  cat(sprintf("  done in %.1f minutes\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit, fit_file)
}
cat("\nAcceptance rates:", paste(sprintf("%.3f", fit$acceptance_rates), collapse = ", "), "\n")

draws <- do.call(rbind, lapply(seq_len(dim(fit$posterior)[3]),
                               function(ch) fit$posterior[, , ch]))
draws <- draws[, key]
cat(sprintf("Posterior draws kept: %d (%d chains)\n", nrow(draws), fit$n_chains))

dg <- fit$diagnostics
dg <- dg[match(key, dg$parameter), ]
cat(sprintf("R-hat: max %.3f (%s); ESS: min %.0f (%s), median %.0f\n",
            max(dg$rhat), dg$parameter[which.max(dg$rhat)],
            min(dg$ess), dg$parameter[which.min(dg$ess)], stats::median(dg$ess)))

# Per chain: log-likelihood and log posterior kernel at 150 evenly spaced
# draws, and means of the parameters that separate the minor mode
n_keep <- dim(fit$posterior)[1]
pick <- round(seq(1, n_keep, length.out = min(150L, n_keep)))
per_chain <- t(sapply(seq_len(fit$n_chains), function(ch) {
  th <- fit$posterior[pick, key, ch, drop = FALSE][, , 1]
  ll <- apply(th, 1, log_lik)
  lp <- ll + apply(th, 1, log_prior)
  c(mean_loglik = mean(ll), mean_logpost = mean(lp), max_logpost = max(lp),
    colMeans(fit$posterior[, c("rho_a", "rho_b", "gamma_bar", "sigma_c"), ch]))
}))
rownames(per_chain) <- paste0("chain", seq_len(fit$n_chains))
cat(sprintf("\nPer chain (log posterior at the mode: %.2f):\n", mode_ll + mode_lp))
print(round(per_chain, 3))

outside <- colMeans(sweep(draws, 2, lower, "<") | sweep(draws, 2, upper, ">"))
if (any(outside > 0)) {
  cat("Share of draws outside the Dynare bounds:\n")
  print(round(outside[outside > 0], 4))
} else {
  cat("All posterior draws lie within the estimated_params bounds.\n")
}

# ---------------------------------------------------------------------------
# 3. Tables 1A and 1B (as in the book's table1.m)
# ---------------------------------------------------------------------------

hpd <- function(x, prob = 0.90) {
  x <- sort(x); n <- length(x); k <- floor(prob * n)
  w <- x[(k + 1):n] - x[1:(n - k)]
  i <- which.min(w)
  c(x[i], x[i + k])
}

label <- c(
  phi = "phi (investment adj. cost)", sigma_c = "sigma_c (risk aversion)",
  lambda = "h (habit)", xi_w = "xi_w (wage Calvo)",
  sigma_l = "sigma_l (labour supply elast.)", xi_p = "xi_p (price Calvo)",
  iota_w = "iota_w (wage indexation)", iota_p = "iota_p (price indexation)",
  psi = "psi (capital utilisation)", phi_p = "Phi (fixed cost)",
  r_pi = "r_pi (Taylor: inflation)", rho = "rho (interest smoothing)",
  r_y = "r_y (Taylor: output gap)", r_dy = "r_dy (Taylor: d output gap)",
  pi_bar = "pi_bar (SS inflation)", beta_const = "100(1/beta - 1)",
  l_bar = "l_bar (SS hours)", gamma_bar = "gamma_bar (trend growth)",
  alpha = "alpha (capital share)",
  sd_e.eta_a = "sigma_a (productivity)", sd_e.eta_b = "sigma_b (risk premium)",
  sd_e.eta_g = "sigma_g (exogenous spending)", sd_e.eta_i = "sigma_i (investment)",
  sd_e.eta_r = "sigma_r (monetary policy)", sd_e.eta_p = "sigma_p (price mark-up)",
  sd_e.eta_w = "sigma_w (wage mark-up)",
  rho_a = "rho_a", rho_b = "rho_b", rho_g = "rho_g", rho_i = "rho_i",
  rho_r = "rho_r", rho_p = "rho_p", rho_w = "rho_w", Mu_p = "mu_p",
  Mu_w = "mu_w", rho_ga = "rho_ga")

prior_txt <- function(nm) {
  e <- est[match(nm, key), ]
  sprintf("%s(%.3g, %.3g)", c(normal = "N", beta = "B", gamma = "G",
          inv_gamma = "IG", inv_gamma1 = "IG", uniform = "U")[[e$shape]],
          e$prior_mean, e$prior_sd)
}

make_table <- function(names) {
  data.frame(
    parameter = names,
    description = unname(label[names]),
    prior = vapply(names, prior_txt, ""),
    mode = unname(mode[names]),
    mean = unname(colMeans(draws[, names, drop = FALSE])),
    hpd_5 = vapply(names, function(nm) hpd(draws[, nm])[1], 0),
    hpd_95 = vapply(names, function(nm) hpd(draws[, nm])[2], 0),
    rhat = dg$rhat[match(names, dg$parameter)],
    ess = round(dg$ess[match(names, dg$parameter)]),
    row.names = NULL)
}

structural <- c("phi", "sigma_c", "lambda", "xi_w", "sigma_l", "xi_p",
                "iota_w", "iota_p", "psi", "phi_p", "r_pi", "rho", "r_y",
                "r_dy", "pi_bar", "beta_const", "l_bar", "gamma_bar", "alpha")
shocks <- c(paste0("sd_e.eta_", c("a", "b", "g", "i", "r", "p", "w")),
            "rho_a", "rho_b", "rho_g", "rho_i", "rho_r", "rho_p", "rho_w",
            "Mu_p", "Mu_w", "rho_ga")
t1a <- make_table(structural)
t1b <- make_table(shocks)

show <- function(tab, title) {
  cat("\n", title, "\n", strrep("=", 98), "\n", sep = "")
  cat(sprintf("%-31s %-14s %7s %7s %7s %7s %6s %6s\n", "Parameter", "Prior",
              "Mode", "Mean", "5%", "95%", "R-hat", "ESS"))
  cat(strrep("-", 98), "\n")
  for (i in seq_len(nrow(tab)))
    cat(sprintf("%-31s %-14s %7.3f %7.3f %7.3f %7.3f %6.3f %6d\n",
                tab$description[i], tab$prior[i], tab$mode[i], tab$mean[i],
                tab$hpd_5[i], tab$hpd_95[i], tab$rhat[i], as.integer(tab$ess[i])))
  cat(strrep("-", 98), "\n")
}
show(t1a, "Table 1A -- Prior and posterior distribution of structural parameters (dsge)")
show(t1b, "Table 1B -- Prior and posterior distribution of shock processes (dsge)")
cat("5% and 95%: bounds of the 90% highest posterior density interval, as in Dynare.\n")

utils::write.csv(t1a, file.path(out, "table_1a_structural.csv"), row.names = FALSE)
utils::write.csv(t1b, file.path(out, "table_1b_shocks.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Figure 6.2 -- prior vs posterior for the six parameters in the book
#    The book's figure script draws normal approximations around the
#    published SW (2007) mode and 90% interval; here the posterior is the
#    kernel density of the MCMC draws and the prior is the file's prior.
#    The published values are read from the book's script and marked.
# ---------------------------------------------------------------------------

book_fig <- file.path(book_dir(), "Chapter6", "Codes_Figures",
                      "chapter6_dsge_companion.R")
sw_pub <- local({
  e <- new.env()
  for (x in parse(book_fig))
    if (is.call(x) && identical(x[[1]], as.name("<-")) &&
        identical(x[[2]], as.name("params"))) eval(x, e)
  e$params
})
fig_map <- c(xi_p = "xi_p", xi_w = "xi_w", h = "lambda", phi = "phi",
             r_pi = "r_pi", rho = "rho")

cmp <- data.frame(
  parameter = names(fig_map),
  sw2007_mode = vapply(names(fig_map), function(k) sw_pub[[k]]$mode, 0),
  sw2007_5 = vapply(names(fig_map), function(k) sw_pub[[k]]$p5, 0),
  sw2007_95 = vapply(names(fig_map), function(k) sw_pub[[k]]$p95, 0),
  dsge_mode = unname(mode[fig_map]),
  dsge_mean = unname(colMeans(draws[, fig_map])),
  dsge_5 = vapply(fig_map, function(nm) hpd(draws[, nm])[1], 0),
  dsge_95 = vapply(fig_map, function(nm) hpd(draws[, nm])[2], 0),
  row.names = NULL)
cat("\nFigure 6.2 parameters: published SW (2007) values (from the book's script) vs dsge\n")
print(format(cmp, digits = 3), row.names = FALSE)
utils::write.csv(cmp, file.path(out, "figure_6_2_comparison.csv"), row.names = FALSE)

prior_density <- function(p, x) exp(vapply(x, function(v) dsge:::dprior(p, v), 0))

col_post <- "#2a78d6"; col_pub <- "#eb6834"; col_prior <- "#8a8984"
col_grid <- "#e4e3df"
titles <- c(xi_p = "xi[p]~~'(price stickiness)'", xi_w = "xi[w]~~'(wage stickiness)'",
            h = "h~~'(habit formation)'", phi = "varphi~~'(investment adj. cost)'",
            r_pi = "r[pi]~~'(Taylor: inflation)'", rho = "rho~~'(interest smoothing)'")

png(file.path(out, "figure_6_2_prior_posterior_dsge.png"),
    width = 1500, height = 950, res = 150)
par(mfrow = c(2, 3), mar = c(3, 1.5, 2.5, 1), oma = c(5, 0, 2.5, 0),
    las = 1, bty = "l", col.axis = "#52514e")
for (k in names(fig_map)) {
  nm <- fig_map[[k]]
  x_draw <- draws[, nm]
  dens <- stats::density(x_draw, adjust = 1.2)
  pr <- m$priors[[nm]]
  rng <- range(stats::quantile(x_draw, c(0.0005, 0.9995)),
               if (pr$distribution == "beta") stats::qbeta(c(0.005, 0.995),
                 pr$params$shape1, pr$params$shape2)
               else pr$params$mean + c(-2.6, 2.6) * pr$params$sd,
               dens$x)
  xs <- seq(rng[1], rng[2], length.out = 400)
  pd <- prior_density(pr, xs)
  ylim <- c(0, max(dens$y, pd) * 1.08)
  plot(xs, pd, type = "n", ylim = ylim, yaxt = "n", xlab = "", ylab = "",
       main = parse(text = titles[[k]])[[1]], font.main = 1, cex.main = 1.15)
  abline(v = pretty(rng), col = col_grid, lwd = 0.6)
  polygon(c(xs, rev(xs)), c(pd, rep(0, length(xs))),
          col = grDevices::adjustcolor(col_prior, 0.25), border = NA)
  lines(xs, pd, col = col_prior, lwd = 1.2)
  polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$x))),
          col = grDevices::adjustcolor(col_post, 0.30), border = NA)
  lines(dens$x, dens$y, col = col_post, lwd = 2)
  abline(v = mode[[nm]], col = col_post, lty = 2, lwd = 1.3)
  abline(v = sw_pub[[k]]$mode, col = col_pub, lty = 1, lwd = 1.6)
}
mtext("Figure 6.2  Bayesian estimation of Smets-Wouters (2007): prior vs posterior (dsge)",
      outer = TRUE, cex = 1.05, font = 2)
# Legend and source line in the outer bottom margin
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
legend("bottom", inset = c(0, 0.045), horiz = TRUE, bty = "n", cex = 0.85,
       c("Prior", "Posterior (dsge MCMC)", "Posterior mode (dsge)",
         "SW (2007) published mode"),
       col = c(col_prior, col_post, col_post, col_pub),
       lwd = c(7, 7, 1.3, 1.6), lty = c(1, 1, 2, 1), seg.len = 1.6,
       text.col = "#0b0b0b")
text(0.5, 0.018, cex = 0.72, font = 3, col = "#52514e",
     sprintf("Book's sw_table1.mod and usmodel_data.mat estimated with dsge::bayes_dsge(): %d chains, %d kept draws each, %s-%s.",
             fit$n_chains, dim(fit$posterior)[1],
             quarter(first + m$estimation$presample), quarter(nrow(dat))))
invisible(dev.off())
cat("\nSaved", file.path(out, "figure_6_2_prior_posterior_dsge.png"), "\n")

# All 36 estimated parameters, with the package's own plot method
pdf(file.path(out, "sw_prior_posterior_all.pdf"), width = 11, height = 8.5)
plot(fit, type = "prior_posterior")
invisible(dev.off())
cat("Saved", file.path(out, "sw_prior_posterior_all.pdf"), "\n")
