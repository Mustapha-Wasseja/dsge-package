# Paccagnini, *Macroeconometrics*, Chapter 6 -- Smets-Wouters (2007)
# Compare dsge's posterior mode (from the book's sw_table1.mod) with the
# Smets-Wouters mode file distributed with Johannes Pfeifer's replication
# (DSGE_mod/Smets_Wouters_2007/usmodel_mode.mat, GPL, not copied here).
#
# 1. Both .mod files are imported with read_dynare(); the parameters are
#    matched by their estimated_params entries (initial value, bounds,
#    prior shape, mean and sd), which the two files share under different
#    names (the book's file renames Smets and Wouters' parameters).
# 2. The 36 values of the mode file (ordered as Pfeifer's estimated_params)
#    are compared with dsge's mode, and the differences are scaled by the
#    posterior standard deviations from the dsge MCMC run.
# 3. The book's posterior kernel (log-likelihood + log-prior, its sample
#    1966Q1-2004Q4) is evaluated at both points.
# 4. dsge's mode search is restarted from the Smets-Wouters mode.
#
# Requires the outputs of ch6_smets_wouters_estimation.R (output/
# sw_posterior_mode.rds and, for posterior standard deviations,
# output/sw_fit.rds). Set PFEIFER_DSGE_MOD_DIR to a clone of
# https://github.com/JohannesPfeifer/DSGE_mod, or leave it unset to clone
# the Smets_Wouters_2007 folder into a temporary directory. Set
# SW_SKIP_RESTART=1 to skip step 4 (about 6 minutes on 4 cores).
#
# Run from the package root:
#   Rscript dev/macroeconometrics-book/ch6_compare_sw_mode.R

suppressPackageStartupMessages(library(dsge))
local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
  source(file.path(here, "book_path.R"))
})
out <- out_dir()
n_cores <- as.integer(Sys.getenv("SW_CORES", min(4L, parallel::detectCores())))

pfeifer_dir <- function() {
  d <- Sys.getenv("PFEIFER_DSGE_MOD_DIR", "")
  if (!nzchar(d)) {
    d <- file.path(tempdir(), "DSGE_mod")
    if (!dir.exists(d)) {
      ok <- system2("git", c("clone", "--quiet", "--depth", "1",
                             "--filter=blob:none", "--sparse",
                             "https://github.com/JohannesPfeifer/DSGE_mod.git", d)) == 0L &&
        system2("git", c("-C", d, "sparse-checkout", "set",
                         "Smets_Wouters_2007")) == 0L
      if (!ok) stop("Could not clone DSGE_mod.", call. = FALSE)
    }
  }
  file.path(normalizePath(d), "Smets_Wouters_2007")
}

# ---------------------------------------------------------------------------
# 1. Match the parameters of the two files
# ---------------------------------------------------------------------------

book_mod <- file.path(book_dir(), "Chapter6", "Codes_Empirical_Examples",
                      "sw_table1.mod")
pf <- pfeifer_dir()
m  <- read_dynare(book_mod)
mp <- suppressWarnings(read_dynare(file.path(pf, "Smets_Wouters_2007.mod")))

key_of <- function(tab) ifelse(tab$type == "stderr", paste0("sd_e.", tab$name),
                               tab$name)
sig <- function(tab) paste(tab$type, signif(tab$init, 10), tab$lower, tab$upper,
                           tab$shape, tab$prior_mean, tab$prior_sd)
eb <- m$estimated_params
ep <- mp$estimated_params
stopifnot(!anyDuplicated(sig(eb)), !anyDuplicated(sig(ep)))
map <- data.frame(book = key_of(eb)[match(sig(ep), sig(eb))],
                  sw = key_of(ep), stringsAsFactors = FALSE)
if (anyNA(map$book) || nrow(map) != 36L)
  stop("The estimated_params entries of the two files do not match one to one.")
cat("All 36 estimated_params entries match one to one (same initial value,",
    "bounds and prior).\n\n")

# Dynare orders the estimated vector: shock standard deviations, then the
# other parameters, each in the order of estimated_params
ord <- c(which(ep$type == "stderr"), which(ep$type == "param"))
mode_file <- R.matlab::readMat(file.path(pf, "usmodel_mode.mat"))
stopifnot(length(mode_file$xparam1) == 36L)
sw_mode <- stats::setNames(as.numeric(mode_file$xparam1), map$book[ord])

# ---------------------------------------------------------------------------
# 2. dsge's mode and posterior standard deviations
# ---------------------------------------------------------------------------

mode_path <- file.path(out, "sw_posterior_mode.rds")
if (!file.exists(mode_path))
  stop("Run ch6_smets_wouters_estimation.R first (needs ", mode_path, ").")
dsge_mode <- readRDS(mode_path)
keys <- names(m$priors)
post_sd <- NULL
fit_path <- file.path(out, "sw_fit.rds")
if (file.exists(fit_path)) {
  P <- readRDS(fit_path)$posterior
  draws <- do.call(rbind, lapply(seq_len(dim(P)[3]), function(ch) P[, , ch]))
  post_sd <- apply(draws[, keys], 2, stats::sd)
}

# ---------------------------------------------------------------------------
# 3. The book's posterior kernel at both points
# ---------------------------------------------------------------------------

raw <- R.matlab::readMat(file.path(dirname(book_mod), "usmodel_data.mat"))
obs <- m$model$variables$observed
Y <- sapply(obs, function(v) as.numeric(raw[[v]]))[m$estimation$first_obs:230, ]
fixed <- unlist(m$model$fixed)
is_sd <- startsWith(keys, "sd_e.")
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
kernel <- function(theta) c(loglik = log_lik(theta), logprior = log_prior(theta))

k_sw <- kernel(sw_mode)
k_dsge <- kernel(dsge_mode)

# ---------------------------------------------------------------------------
# 4. Restart the mode search from the Smets-Wouters mode
# ---------------------------------------------------------------------------

restart <- NULL
restart_path <- file.path(out, "sw_mode_from_published.rds")
if (!nzchar(Sys.getenv("SW_SKIP_RESTART"))) {
  if (file.exists(restart_path)) {
    restart <- readRDS(restart_path)
  } else {
    lower <- stats::setNames(eb$lower, key_of(eb))[keys]
    upper <- stats::setNames(eb$upper, key_of(eb))[keys]
    scale <- pmin(stats::setNames(eb$prior_sd, key_of(eb))[keys], upper - lower)
    obj <- function(x) {
      names(x) <- keys
      v <- sum(kernel(x))
      if (is.finite(v)) -v else 1e10
    }
    grad <- function(x) {
      h <- 1e-5 * scale
      unlist(parallel::mclapply(seq_along(x), function(j) {
        up <- min(x[j] + h[j], upper[j]); dn <- max(x[j] - h[j], lower[j])
        xu <- x; xu[j] <- up; xd <- x; xd[j] <- dn
        (obj(xu) - obj(xd)) / (up - dn)
      }, mc.cores = n_cores))
    }
    cat("Restarting the mode search from the Smets-Wouters mode...\n")
    x <- pmin(pmax(sw_mode[keys], lower), upper)
    for (round in 1:5) {
      o <- stats::optim(x, obj, grad, method = "L-BFGS-B", lower = lower,
                        upper = upper, control = list(maxit = 1000,
                        parscale = scale, factr = 1e7))
      improved <- obj(x) - o$value
      x <- o$par
      if (improved < 1e-3) break
    }
    restart <- stats::setNames(x, keys)
    saveRDS(restart, restart_path)
  }
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

tab <- data.frame(
  parameter = keys,
  sw_name = map$sw[match(keys, map$book)],
  sw2007_mode = unname(sw_mode[keys]),
  dsge_mode = unname(dsge_mode[keys]),
  difference = unname(dsge_mode[keys] - sw_mode[keys]),
  stringsAsFactors = FALSE)
if (!is.null(post_sd)) {
  tab$posterior_sd <- unname(post_sd[keys])
  tab$diff_in_sd <- tab$difference / tab$posterior_sd
}
if (!is.null(restart)) tab$dsge_mode_from_sw <- unname(restart[keys])

cat("Posterior mode: Smets-Wouters mode file vs dsge (book's sw_table1.mod)\n")
print(format(tab, digits = 3), row.names = FALSE)
utils::write.csv(tab, file.path(out, "sw_mode_vs_published.csv"),
                 row.names = FALSE)

cat("\nBook's posterior kernel (sample 1966Q1-2004Q4):\n")
kt <- rbind(`Smets-Wouters mode file` = k_sw, `dsge mode` = k_dsge)
if (!is.null(restart)) kt <- rbind(kt, `dsge, restarted from SW mode` = kernel(restart))
kt <- cbind(kt, logpost = rowSums(kt))
print(round(kt, 4))
if (!is.null(post_sd)) {
  cat(sprintf("\n|difference| / posterior sd: median %.2f, max %.2f (%s)\n",
              stats::median(abs(tab$diff_in_sd)), max(abs(tab$diff_in_sd)),
              tab$parameter[which.max(abs(tab$diff_in_sd))]))
}
if (!is.null(restart)) {
  cat(sprintf("Restarted search vs dsge mode: max |difference| %.2g (%s)\n",
              max(abs(restart[keys] - dsge_mode[keys])),
              keys[which.max(abs(restart[keys] - dsge_mode[keys]))]))
}
