# Paccagnini, *Macroeconometrics*, Chapter 6 -- Smets-Wouters (2007)
# Compare dsge's estimation of the book's sw_table1.mod with Dynare's own
# estimation of the same file.
#
# Dynare side (run separately; see README "Dynare comparison"):
#   - the book's run: posterior mode, then 2 blocks x 2000 MH draws
#     (exports <book>_summary.csv, _scalars.csv, _draws.csv);
#   - long MH chains started at Dynare's mode (exports sw_long*_draws.csv).
#   export_dynare.m writes those files, and evaluates Dynare's own
#   log-likelihood and log-prior at the parameter vectors this script
#   writes to points.csv (-> <book>_points_eval.csv).
#
# This script
#   1. writes points.csv: the file's initial values, both modes and both
#      posterior means;
#   2. evaluates dsge's log-likelihood and log-prior at the same points and
#      compares them with Dynare's values (same model, same data?);
#   3. compares the two posterior modes, and the Laplace approximations of
#      the marginal data density;
#   4. compares the posterior distributions: dsge's bayes_dsge() chains
#      against Dynare's long chains (and the book's 2 x 2000 run), with
#      Monte Carlo standard errors so that differences can be judged
#      against sampling noise;
#   5. draws the posterior densities from both programs.
#
# Requires the outputs of ch6_smets_wouters_estimation.R (output/
# sw_posterior_mode.rds; output/sw_fit.rds for the posterior sections).
# SW_DYNARE_DIR is the folder holding the MATLAB exports (default
# output/dynare/matlab), SW_DYNARE_OCTAVE_DIR the Octave exports (default
# output/dynare/octave), SW_DYNARE_BOOK the book run's file prefix (default
# sw_table1). The *_draws.csv files are not committed (about 90 MB).
#
# Run from the package root:
#   Rscript dev/macroeconometrics-book/ch6_compare_dynare.R

suppressPackageStartupMessages(library(dsge))
local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
  source(file.path(here, "book_path.R"))
})
out     <- out_dir()
dyn_dir <- Sys.getenv("SW_DYNARE_DIR", file.path(out, "dynare", "matlab"))
book    <- Sys.getenv("SW_DYNARE_BOOK", "sw_table1")
dyn_label <- Sys.getenv("SW_DYNARE_LABEL", "MATLAB, fmincon")
dir.create(dyn_dir, showWarnings = FALSE, recursive = TRUE)
dyn_file <- function(...) file.path(dyn_dir, paste0(...))

# ---------------------------------------------------------------------------
# 0. The book's model, data and posterior kernel in dsge
#    (the same construction as ch6_smets_wouters_estimation.R)
# ---------------------------------------------------------------------------

ex_dir <- file.path(book_dir(), "Chapter6", "Codes_Empirical_Examples")
m      <- read_dynare(file.path(ex_dir, "sw_table1.mod"))
est    <- m$estimated_params
keys   <- ifelse(est$type == "stderr", paste0("sd_e.", est$name), est$name)
is_sd  <- startsWith(keys, "sd_e.")
fixed  <- unlist(m$model$fixed)
obs    <- m$model$variables$observed
raw    <- R.matlab::readMat(file.path(ex_dir, "usmodel_data.mat"))
Y      <- sapply(obs, function(v) as.numeric(raw[[v]]))
Y      <- Y[m$estimation$first_obs:nrow(Y), , drop = FALSE]

log_lik <- function(theta) {
  theta <- theta[keys]
  sol <- tryCatch(solve_dsge(m$model, params = c(theta[!is_sd], fixed),
                             shock_sd = stats::setNames(theta[is_sd],
                               sub("sd_e.", "", keys[is_sd], fixed = TRUE))),
                  error = function(e) NULL)
  if (is.null(sol) || !isTRUE(sol$stable)) return(-Inf)
  dsge:::kalman_filter(sweep(Y, 2, sol$steady_state[obs]), sol$G, sol$H,
                       sol$M, sol$D, presample = m$estimation$presample,
                       init = m$model$kalman_init, loglik_only = TRUE)$loglik
}
log_prior <- function(theta)
  sum(vapply(keys, function(k) dsge:::dprior(m$priors[[k]], theta[[k]]), 0))
log_post <- function(theta) {
  names(theta) <- keys
  lp <- log_prior(theta)
  if (!is.finite(lp)) return(-Inf)
  lp + log_lik(theta)
}

hpd <- function(x, prob = 0.90) {
  x <- sort(unname(x)); n <- length(x); k <- floor(prob * n)
  w <- x[(k + 1):n] - x[1:(n - k)]
  i <- which.min(w)
  c(x[i], x[i + k])
}

# ---------------------------------------------------------------------------
# 1. Inputs: dsge's mode and draws, Dynare's exports
# ---------------------------------------------------------------------------

mode_path <- file.path(out, "sw_posterior_mode.rds")
fit_path  <- file.path(out, "sw_fit.rds")
if (!file.exists(mode_path))
  stop("Run ch6_smets_wouters_estimation.R first (needs ", mode_path, ").")
dsge_mode <- readRDS(mode_path)[keys]
# The draws are optional: without them only sections 2, 3 and the modes in
# section 6 are produced
fit <- dsge_chains <- dsge_draws <- NULL
if (file.exists(fit_path)) {
  fit <- readRDS(fit_path)
  dsge_chains <- lapply(seq_len(dim(fit$posterior)[3]),
                        function(ch) fit$posterior[, keys, ch])
  dsge_draws <- do.call(rbind, dsge_chains)
} else {
  cat("(", fit_path, " not found: posterior sections skipped)\n\n", sep = "")
}

read_csv <- function(f) utils::read.csv(f, check.names = FALSE,
                                        stringsAsFactors = FALSE)
dyn_sum <- read_csv(dyn_file(book, "_summary.csv"))
stopifnot(setequal(dyn_sum$parameter, keys))
dyn_sum <- dyn_sum[match(keys, dyn_sum$parameter), ]
dyn_sc  <- read_csv(dyn_file(book, "_scalars.csv"))
dyn_scalar <- function(q) {
  v <- dyn_sc$value[dyn_sc$quantity == q]
  if (length(v)) as.numeric(v) else NA_real_
}
dyn_mode <- stats::setNames(dyn_sum$mode, keys)

# The bounds Dynare used are the bounds dsge enforces
stopifnot(isTRUE(all.equal(dyn_sum$lower, est$lower)),
          isTRUE(all.equal(dyn_sum$upper, est$upper)))

# Chains from a Dynare draws file: one per block, the first mh_drop share
# of each block dropped (Dynare's mh_drop)
dynare_chains <- function(files, drop) {
  ch <- list()
  for (f in files) {
    d <- read_csv(f)
    for (b in unique(d$block)) {
      db <- d[d$block == b, ]
      n_total <- max(db$iter)
      keep <- db$iter > floor(drop * n_total)
      ch[[length(ch) + 1]] <- as.matrix(db[keep, keys])
    }
  }
  # equal lengths for R-hat
  n <- min(vapply(ch, nrow, 0L))
  lapply(ch, function(x) x[seq_len(n), , drop = FALSE])
}

book_chains <- NULL
if (file.exists(dyn_file(book, "_draws.csv")))
  book_chains <- dynare_chains(dyn_file(book, "_draws.csv"),
                               drop = dyn_scalar("mh_drop"))

# Optional second Dynare platform (e.g. Octave), exported to its own folder
oct_dir <- Sys.getenv("SW_DYNARE_OCTAVE_DIR", file.path(out, "dynare", "octave"))
oct_label <- Sys.getenv("SW_DYNARE_OCTAVE_LABEL", "Dynare (Octave)")
oct <- NULL
if (nzchar(oct_dir) && file.exists(file.path(oct_dir, paste0(book, "_summary.csv")))) {
  of <- function(...) file.path(oct_dir, paste0(book, ...))
  s <- read_csv(of("_summary.csv")); s <- s[match(keys, s$parameter), ]
  sc <- read_csv(of("_scalars.csv"))
  oct <- list(summary = s,
              scalar = function(q) { v <- sc$value[sc$quantity == q]
                                     if (length(v)) as.numeric(v) else NA_real_ },
              eval = if (file.exists(of("_points_eval.csv"))) read_csv(of("_points_eval.csv")),
              chains = if (file.exists(of("_draws.csv")))
                dynare_chains(of("_draws.csv"), drop = as.numeric(sc$value[sc$quantity == "mh_drop"])))
  oct$mode <- stats::setNames(s$mode, keys)
}
long_files <- list.files(dyn_dir, pattern = "^sw_long.*_draws\\.csv$",
                         full.names = TRUE)
long_drop <- as.numeric(Sys.getenv("SW_DYNARE_LONG_DROP", "0.3333333333"))
long_chains <- if (length(long_files)) dynare_chains(long_files, long_drop)
long_draws  <- if (length(long_files)) do.call(rbind, long_chains)

# ---------------------------------------------------------------------------
# 2. The posterior kernel at common points, in both programs
# ---------------------------------------------------------------------------

points <- list(initial_values = stats::setNames(est$init, keys),
               dynare_mode = dyn_mode,
               dsge_mode = dsge_mode)
if (!is.null(dsge_draws)) points$dsge_posterior_mean <- colMeans(dsge_draws)
if (!is.null(oct)) points$octave_mode <- oct$mode
if (!is.null(long_draws)) points$dynare_posterior_mean <- colMeans(long_draws)
pts <- data.frame(point = names(points), do.call(rbind, points),
                  check.names = FALSE, row.names = NULL)
utils::write.csv(pts, dyn_file("points.csv"), row.names = FALSE, quote = FALSE)
if (nzchar(oct_dir))
  utils::write.csv(pts, file.path(oct_dir, "points.csv"), row.names = FALSE,
                   quote = FALSE)

dsge_eval <- t(vapply(points, function(th)
  c(loglik = log_lik(th), logprior = log_prior(th)), c(loglik = 0, logprior = 0)))

cat("Posterior kernel at common points (book's sample, 1966Q1-2004Q4)\n")
cat(strrep("=", 100), "\n")
eval_file <- dyn_file(book, "_points_eval.csv")
kernel_tab <- data.frame(point = names(points),
                         dsge_loglik = dsge_eval[, "loglik"],
                         dsge_logprior = dsge_eval[, "logprior"],
                         row.names = NULL)
if (file.exists(eval_file)) {
  de <- read_csv(eval_file)
  de <- de[match(kernel_tab$point, de$point), ]
  kernel_tab$dynare_loglik <- de$loglik
  kernel_tab$dynare_logprior <- de$logprior
  kernel_tab$diff_loglik <- kernel_tab$dsge_loglik - kernel_tab$dynare_loglik
  kernel_tab$diff_logprior <- kernel_tab$dsge_logprior - kernel_tab$dynare_logprior
  cat(sprintf("%-22s %15s %15s %10s %13s %13s %10s\n", "Point", "loglik dsge",
              "loglik Dynare", "diff", "prior dsge", "prior Dynare", "diff"))
  for (i in seq_len(nrow(kernel_tab))) with(kernel_tab[i, ], cat(sprintf(
    "%-22s %15.6f %15.6f %10.1e %13.6f %13.6f %10.1e\n", point, dsge_loglik,
    dynare_loglik, diff_loglik, dsge_logprior, dynare_logprior, diff_logprior)))
} else {
  cat("(", eval_file, " not found: run export_dynare.m with points.csv in place\n",
      " to get Dynare's values at these points)\n", sep = "")
  print(kernel_tab, digits = 10, row.names = FALSE)
}
if (!is.null(oct$eval)) {
  oe <- oct$eval[match(kernel_tab$point, oct$eval$point), ]
  kernel_tab$octave_loglik <- oe$loglik
  kernel_tab$octave_logprior <- oe$logprior
  cat(sprintf("\n%s at the same points:\n", oct_label))
  cat(sprintf("%-22s %15s %10s %13s %10s\n", "Point", "loglik", "vs dsge",
              "logprior", "vs dsge"))
  for (i in seq_len(nrow(kernel_tab))) with(kernel_tab[i, ], if (!is.na(octave_loglik))
    cat(sprintf("%-22s %15.6f %10.1e %13.6f %10.1e\n", point, octave_loglik,
                dsge_loglik - octave_loglik, octave_logprior,
                dsge_logprior - octave_logprior)))
}
cat(sprintf("Dynare's reported log posterior at its mode: %.4f\n",
            dyn_scalar("logpost_at_mode")))
utils::write.csv(kernel_tab, file.path(out, "dynare_kernel_at_points.csv"),
                 row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Posterior modes and Laplace approximations
# ---------------------------------------------------------------------------

lp_dyn_mode  <- sum(dsge_eval["dynare_mode", ])
lp_dsge_mode <- sum(dsge_eval["dsge_mode", ])

# Laplace approximation at dsge's mode with a numerical Hessian (dsge side)
laplace_file <- file.path(out, "sw_laplace.rds")
lap <- if (file.exists(laplace_file)) readRDS(laplace_file)
if (is.null(lap) || !isTRUE(all.equal(lap$at, dsge_mode))) {   # cache is per mode
  cat("\nNumerical Hessian of the log posterior at dsge's mode...\n")
  H <- numDeriv::hessian(log_post, dsge_mode,
                         method.args = list(eps = 1e-4, d = 1e-3, r = 4))
  V <- solve(-H)
  ld <- as.numeric(determinant(-H, logarithm = TRUE)$modulus)
  lap <- list(at = dsge_mode, H = H, sd = sqrt(diag(V)),
              log_mdd = lp_dsge_mode + length(keys) / 2 * log(2 * pi) - 0.5 * ld)
  saveRDS(lap, laplace_file)
}

mode_tab <- data.frame(
  parameter = keys,
  dynare_mode = dyn_sum$mode,
  dynare_sd_at_mode = dyn_sum$sd_at_mode,
  dsge_mode = unname(dsge_mode),
  dsge_sd_at_mode = unname(lap$sd),
  difference = unname(dsge_mode) - dyn_sum$mode,
  row.names = NULL)
mode_tab$diff_in_sd <- mode_tab$difference / mode_tab$dynare_sd_at_mode

cat("\nPosterior mode: Dynare (", dyn_label, ") vs dsge (L-BFGS-B)\n", sep = "")
cat(strrep("=", 100), "\n")
cat(sprintf("%-12s %12s %10s %12s %10s %11s %9s\n", "Parameter", "Dynare mode",
            "Dynare sd", "dsge mode", "dsge sd", "difference", "diff/sd"))
for (i in seq_len(nrow(mode_tab))) with(mode_tab[i, ], cat(sprintf(
  "%-12s %12.4f %10.4f %12.4f %10.4f %+11.4f %+9.2f\n", parameter, dynare_mode,
  dynare_sd_at_mode, dsge_mode, dsge_sd_at_mode, difference, diff_in_sd)))
cat(sprintf("|difference| / sd at mode: median %.3f, max %.3f (%s)\n",
            stats::median(abs(mode_tab$diff_in_sd)), max(abs(mode_tab$diff_in_sd)),
            keys[which.max(abs(mode_tab$diff_in_sd))]))
cat(sprintf("Log posterior (dsge's evaluation): at Dynare's mode %.4f, at dsge's mode %.4f\n",
            lp_dyn_mode, lp_dsge_mode))
if (!is.null(oct)) {
  d_oct <- (oct$mode - dyn_mode) / dyn_sum$sd_at_mode
  cat(sprintf(paste0("%s mode vs %s mode: |difference| / sd median %.3f, max %.3f (%s); ",
                     "log posterior at it (dsge's evaluation) %.4f\n"),
              oct_label, dyn_label, stats::median(abs(d_oct)), max(abs(d_oct)),
              keys[which.max(abs(d_oct))], sum(dsge_eval["octave_mode", ])))
}
cat(sprintf("Laplace log marginal data density: Dynare %.4f, dsge %.4f\n",
            dyn_scalar("log_mdd_laplace"), lap$log_mdd))
utils::write.csv(mode_tab, file.path(out, "dynare_mode_comparison.csv"),
                 row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Posterior distributions
# ---------------------------------------------------------------------------

summarise <- function(chains) {
  dr <- do.call(rbind, chains)
  t(vapply(keys, function(k) {
    ch <- lapply(chains, function(x) x[, k])
    ess <- dsge:::compute_ess(ch)
    h <- hpd(dr[, k])
    c(mean = mean(dr[, k]), sd = stats::sd(dr[, k]), hpd_5 = h[1], hpd_95 = h[2],
      ess = ess, mcse = stats::sd(dr[, k]) / sqrt(ess),
      rhat = dsge:::compute_rhat(ch))
  }, numeric(7)))
}
if (!is.null(dsge_chains)) {        # (section 4 needs dsge's draws)
s_dsge <- summarise(dsge_chains)

post_tab <- data.frame(parameter = keys, row.names = NULL)
add <- function(tab, s, prefix) {
  for (nm in colnames(s)) tab[[paste0(prefix, "_", nm)]] <- unname(s[, nm])
  tab
}
post_tab <- add(post_tab, s_dsge, "dsge")
ref <- NULL
if (!is.null(long_chains)) {
  s_long <- summarise(long_chains)
  post_tab <- add(post_tab, s_long, "dynare")
  post_tab$diff_mean_in_sd <- (post_tab$dsge_mean - post_tab$dynare_mean) /
    post_tab$dynare_sd
  post_tab$diff_mean_in_mcse <- (post_tab$dsge_mean - post_tab$dynare_mean) /
    sqrt(post_tab$dsge_mcse^2 + post_tab$dynare_mcse^2)
  post_tab$sd_ratio <- post_tab$dsge_sd / post_tab$dynare_sd
  ref <- "dynare"
}
if (!is.null(book_chains)) {
  s_book <- summarise(book_chains)
  post_tab <- add(post_tab, s_book, "book")
}

cat("\nPosterior distributions\n")
cat(strrep("=", 100), "\n")
cat(sprintf("dsge: %d chains x %d kept draws (bayes_dsge(), thin %d)\n",
            length(dsge_chains), nrow(dsge_chains[[1]]), fit$thin %||% NA))
if (!is.null(long_chains))
  cat(sprintf("Dynare long run: %d chains x %d kept draws (after dropping %.0f%%)\n",
              length(long_chains), nrow(long_chains[[1]]), 100 * long_drop))
if (!is.null(book_chains))
  cat(sprintf("Dynare book run: %d blocks x %d kept draws (mh_drop %.2f)\n",
              length(book_chains), nrow(book_chains[[1]]), dyn_scalar("mh_drop")))

if (!is.null(ref)) {
  cat(sprintf("\n%-12s %8s %8s %8s %8s | %8s %8s %8s %8s | %7s %7s %6s %6s\n",
              "Parameter", "dsge", "5%", "95%", "sd", "Dynare", "5%", "95%", "sd",
              "d/sd", "d/mcse", "Rh-ds", "Rh-dy"))
  for (i in seq_len(nrow(post_tab))) with(post_tab[i, ], cat(sprintf(
    "%-12s %8.4f %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %8.4f | %+7.3f %+7.2f %6.3f %6.3f\n",
    parameter, dsge_mean, dsge_hpd_5, dsge_hpd_95, dsge_sd, dynare_mean,
    dynare_hpd_5, dynare_hpd_95, dynare_sd, diff_mean_in_sd, diff_mean_in_mcse,
    dsge_rhat, dynare_rhat)))
  cat(sprintf(paste0("\nMean difference / posterior sd: median %.3f, max %.3f (%s)\n",
                     "Mean difference / MCSE: median %.2f, max %.2f (%s); ",
                     "%d of %d within 2 MCSE\n",
                     "sd ratio dsge/Dynare: %.3f to %.3f\n",
                     "HPD bounds, |difference| / sd: median %.3f, max %.3f\n"),
              stats::median(abs(post_tab$diff_mean_in_sd)),
              max(abs(post_tab$diff_mean_in_sd)),
              keys[which.max(abs(post_tab$diff_mean_in_sd))],
              stats::median(abs(post_tab$diff_mean_in_mcse)),
              max(abs(post_tab$diff_mean_in_mcse)),
              keys[which.max(abs(post_tab$diff_mean_in_mcse))],
              sum(abs(post_tab$diff_mean_in_mcse) <= 2), nrow(post_tab),
              min(post_tab$sd_ratio), max(post_tab$sd_ratio),
              stats::median(abs(c(post_tab$dsge_hpd_5 - post_tab$dynare_hpd_5,
                                  post_tab$dsge_hpd_95 - post_tab$dynare_hpd_95) /
                                rep(post_tab$dynare_sd, 2))),
              max(abs(c(post_tab$dsge_hpd_5 - post_tab$dynare_hpd_5,
                        post_tab$dsge_hpd_95 - post_tab$dynare_hpd_95) /
                      rep(post_tab$dynare_sd, 2)))))
}
if (!is.null(book_chains)) {
  base <- if (!is.null(ref)) post_tab$dynare_mean else post_tab$dsge_mean
  base_sd <- if (!is.null(ref)) post_tab$dynare_sd else post_tab$dsge_sd
  d <- (post_tab$book_mean - base) / base_sd
  cat(sprintf(paste0("Book's short Dynare run vs the %s long run: mean difference / sd ",
                     "median %.2f, max %.2f (%s); its R-hat up to %.2f, ESS %d to %d\n"),
              if (!is.null(ref)) "Dynare" else "dsge", stats::median(abs(d)),
              max(abs(d)), keys[which.max(abs(d))], max(post_tab$book_rhat),
              as.integer(min(post_tab$book_ess)), as.integer(max(post_tab$book_ess))))
}
for (nm in c("log_mdd_laplace", "log_mdd_mhm"))
  if (!is.na(dyn_scalar(nm))) cat(sprintf("Dynare %s (book run): %.4f\n", nm, dyn_scalar(nm)))
utils::write.csv(post_tab, file.path(out, "dynare_posterior_comparison.csv"),
                 row.names = FALSE)
}                                   # end of section 4

# ---------------------------------------------------------------------------
# 5. Posterior densities, dsge vs Dynare
# ---------------------------------------------------------------------------

if (!is.null(long_draws) && !is.null(dsge_draws)) {
  col_dsge <- "#2a78d6"; col_dyn <- "#eb6834"; col_prior <- "#8a8984"
  col_grid <- "#e4e3df"; col_ink2 <- "#52514e"
  prior_density <- function(p, x) exp(vapply(x, function(v) dsge:::dprior(p, v), 0))

  panel <- function(k, title = k, show_prior = TRUE) {
    a <- dsge_draws[, k]; b <- long_draws[, k]
    da <- stats::density(a, adjust = 1.2); db <- stats::density(b, adjust = 1.2)
    rng <- range(da$x, db$x)
    xs <- seq(rng[1], rng[2], length.out = 300)
    pd <- if (show_prior) prior_density(m$priors[[k]], xs) else 0 * xs
    ylim <- c(0, max(da$y, db$y) * 1.1)
    plot(NA, xlim = rng, ylim = ylim, yaxt = "n", xlab = "", ylab = "",
         main = title, font.main = 1, cex.main = 1.0)
    abline(v = pretty(rng), col = col_grid, lwd = 0.6)
    if (show_prior) {
      polygon(c(xs, rev(xs)), c(pmin(pd, ylim[2]), rep(0, length(xs))),
              col = grDevices::adjustcolor(col_prior, 0.18), border = NA)
      lines(xs, pmin(pd, ylim[2]), col = col_prior, lwd = 1)
    }
    lines(db$x, db$y, col = col_dyn, lwd = 2)
    lines(da$x, da$y, col = col_dsge, lwd = 2, lty = 2)
    abline(v = dyn_mode[[k]], col = col_dyn, lwd = 1, lty = 3)
    abline(v = dsge_mode[[k]], col = col_dsge, lwd = 1, lty = 3)
  }
  legend_row <- function(cex = 0.85) {
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    legend("bottom", inset = c(0, 0.01), horiz = TRUE, bty = "n", cex = cex,
           c("Prior", "Posterior, Dynare", "Posterior, dsge",
             "Mode (dotted, same colours)"),
           col = c(col_prior, col_dyn, col_dsge, col_ink2),
           lwd = c(6, 2, 2, 1), lty = c(1, 1, 2, 3), seg.len = 2,
           text.col = "#0b0b0b")
  }

  six <- c(xi_p = "xi_p (price stickiness)", xi_w = "xi_w (wage stickiness)",
           lambda = "h (habit)", phi = "phi (investment adj. cost)",
           r_pi = "r_pi (Taylor: inflation)", rho = "rho (interest smoothing)")
  png(file.path(out, "dynare_vs_dsge_posteriors.png"), width = 1500,
      height = 950, res = 150)
  par(mfrow = c(2, 3), mar = c(3, 1.5, 2.5, 1), oma = c(3.5, 0, 2.5, 0),
      las = 1, bty = "l", col.axis = col_ink2)
  for (k in names(six)) panel(k, six[[k]])
  mtext("Smets-Wouters (2007), book's sw_table1.mod: posterior from Dynare and from dsge",
        outer = TRUE, cex = 1.0, font = 2)
  legend_row()
  invisible(dev.off())

  pdf(file.path(out, "dynare_vs_dsge_posteriors_all.pdf"), width = 11, height = 8.5)
  for (page in split(keys, ceiling(seq_along(keys) / 12))) {
    par(mfrow = c(3, 4), mar = c(3, 1.5, 2.5, 1), oma = c(3, 0, 2, 0),
        las = 1, bty = "l", col.axis = col_ink2)
    for (k in page) panel(k)
    mtext("Posterior densities: Dynare (solid) vs dsge (dashed), prior in grey",
          outer = TRUE, cex = 0.95, font = 2)
    legend_row(0.8)
  }
  invisible(dev.off())
  cat("\nSaved", file.path(out, "dynare_vs_dsge_posteriors.png"), "and",
      file.path(out, "dynare_vs_dsge_posteriors_all.pdf"), "\n")
}

# ---------------------------------------------------------------------------
# 6. One table: every estimated parameter on every platform
#    (Table 1A/1B order of the book's table1.m)
# ---------------------------------------------------------------------------

order_keys <- c("phi", "sigma_c", "lambda", "xi_w", "sigma_l", "xi_p", "iota_w",
                "iota_p", "psi", "phi_p", "r_pi", "rho", "r_y", "r_dy", "pi_bar",
                "beta_const", "l_bar", "gamma_bar", "alpha",
                paste0("sd_e.eta_", c("a", "b", "g", "i", "r", "p", "w")),
                "rho_a", "rho_b", "rho_g", "rho_i", "rho_r", "rho_p", "rho_w",
                "Mu_p", "Mu_w", "rho_ga")
stopifnot(setequal(order_keys, keys))
shape_abbr <- c(normal = "N", beta = "B", gamma = "G", inv_gamma = "IG",
                inv_gamma1 = "IG", uniform = "U")
prior_txt <- function(k) {
  e <- est[match(k, keys), ]
  sprintf("%s(%.3g, %.3g)", shape_abbr[[e$shape]], e$prior_mean, e$prior_sd)
}
stat_of <- function(chains, k, what) {
  if (is.null(chains)) return(NA_real_)
  x <- unlist(lapply(chains, function(ch) ch[, k]))
  switch(what, mean = mean(x), lo = hpd(x)[1], hi = hpd(x)[2])
}
plat <- data.frame(parameter = order_keys,
                   prior = vapply(order_keys, prior_txt, ""),
                   mode_dynare_matlab = unname(dyn_mode[order_keys]),
                   mode_dynare_octave = if (!is.null(oct)) unname(oct$mode[order_keys]) else NA,
                   mode_dsge = unname(dsge_mode[order_keys]),
                   row.names = NULL, stringsAsFactors = FALSE)
for (src in list(list("dynare_long", long_chains), list("dsge", dsge_chains),
                 list("dynare_matlab_book", book_chains),
                 list("dynare_octave_book", oct$chains))) {
  for (w in c("mean", "lo", "hi"))
    plat[[paste0(src[[1]], "_", c(mean = "mean", lo = "hpd5", hi = "hpd95")[[w]])]] <-
      vapply(order_keys, function(k) stat_of(src[[2]], k, w), 0)
}
utils::write.csv(plat, file.path(out, "platform_comparison.csv"), row.names = FALSE)

fmt <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
cell <- function(m, lo, hi) ifelse(is.na(m), "", sprintf("%s [%s, %s]", fmt(m), fmt(lo), fmt(hi)))
md <- c(
  "| Parameter | Prior | Mode: Dynare (MATLAB, fmincon) | Mode: Dynare (Octave, csminwel) | Mode: dsge (R, L-BFGS-B) | Posterior mean [90% HPD]: Dynare (MATLAB), 4 x 150k | Posterior mean [90% HPD]: dsge, 4 x 150k | Book settings, 2 x 2000: Dynare (MATLAB) | Book settings, 2 x 2000: Dynare (Octave) |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(plat))) with(plat[i, ], {
  md <<- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    sub("sd_e.", "sd ", parameter, fixed = TRUE), prior, fmt(mode_dynare_matlab),
    fmt(mode_dynare_octave), fmt(mode_dsge),
    cell(dynare_long_mean, dynare_long_hpd5, dynare_long_hpd95),
    cell(dsge_mean, dsge_hpd5, dsge_hpd95),
    cell(dynare_matlab_book_mean, dynare_matlab_book_hpd5, dynare_matlab_book_hpd95),
    cell(dynare_octave_book_mean, dynare_octave_book_hpd5, dynare_octave_book_hpd95)))
})
lp_row <- function(lab, a, b, c) sprintf("| %s | | %s | %s | %s | | | | |", lab, a, b, c)
md <- c(md,
  lp_row("Log posterior at the mode (dsge's evaluation)", fmt(lp_dyn_mode, 4),
         if (!is.null(oct)) fmt(sum(dsge_eval["octave_mode", ]), 4) else "",
         fmt(lp_dsge_mode, 4)),
  lp_row("Log posterior at the mode (as reported)", fmt(dyn_scalar("logpost_at_mode"), 4),
         if (!is.null(oct)) fmt(oct$scalar("logpost_at_mode"), 4) else "",
         fmt(lp_dsge_mode, 4)),
  lp_row("Laplace log marginal data density",
         fmt(dyn_scalar("log_mdd_laplace"), 4),
         if (!is.null(oct)) fmt(oct$scalar("log_mdd_laplace"), 4) else "",
         fmt(lap$log_mdd, 4)))
writeLines(md, file.path(out, "platform_comparison.md"))
cat("\nSaved", file.path(out, "platform_comparison.csv"), "and",
    file.path(out, "platform_comparison.md"), "\n")
