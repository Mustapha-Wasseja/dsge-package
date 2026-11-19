# Paccagnini, *Macroeconometrics*, Chapter 6 -- Smets-Wouters (2007)
# Run the book's sw_table1.mod in Dynare and compare with dsge.
#
# Dynare 6.0 runs in Octave (`apt-get install octave dynare` on Ubuntu
# 24.04). The book's files are copied to a work directory and one change is
# made: the file's mode search, mode_compute = 1 with
# optim = ('Algorithm', 'active-set'), uses MATLAB's fmincon, which Octave
# does not have, so it becomes mode_compute = 4 (csminwel, Dynare's
# default). Everything else is the book's: priors, bounds, sample
# (first_obs = 71, presample = 4), lik_init = 2, and the Metropolis-Hastings
# settings (2 blocks of 2000 draws, mh_jscale = 0.2, mh_drop = 0.2).
#
# Compared with dsge (from ch6_smets_wouters_estimation.R's cached output):
#   1. the log posterior kernel, evaluated by each program at both modes;
#   2. the posterior modes of all 36 parameters;
#   3. posterior means and 90% HPD intervals (Dynare's 2 x 2000 draws
#      against dsge's 4 x 20000 kept draws).
#
# Set DYNARE_WORK_DIR to reuse a finished run; DYNARE_MATLAB if Dynare's
# matlab directory is not /usr/lib/dynare/matlab.
#
# Run from the package root:
#   Rscript dev/macroeconometrics-book/ch6_compare_dynare.R

suppressPackageStartupMessages(library(dsge))
local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
  source(file.path(here, "book_path.R"))
})
out <- out_dir()
ex_dir <- file.path(book_dir(), "Chapter6", "Codes_Empirical_Examples")
dynare_matlab <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")

mode_path <- file.path(out, "sw_posterior_mode.rds")
fit_path <- file.path(out, "sw_fit.rds")
if (!file.exists(mode_path) || !file.exists(fit_path))
  stop("Run ch6_smets_wouters_estimation.R first (needs ", mode_path,
       " and ", fit_path, ").")
dsge_mode <- readRDS(mode_path)

# ---------------------------------------------------------------------------
# Run Dynare (unless a finished run is given)
# ---------------------------------------------------------------------------

octave_driver <- c(
  sprintf("addpath('%s');", dynare_matlab),
  "t0 = tic;",
  "dynare sw_table1 noclearall",
  "fprintf('Dynare estimation took %.1f minutes\\n', toc(t0) / 60);",
  "names = bayestopt_.name; n = numel(names); xmode = zeros(n, 1);",
  "fid = fopen('dynare_posterior.csv', 'w');",
  "fprintf(fid, 'name,mode,mean,hpd_inf,hpd_sup\\n');",
  "for i = 1:n",
  "  nm = names{i};",
  "  if isfield(oo_.posterior_mode.parameters, nm), f = 'parameters'; else, f = 'shocks_std'; end",
  "  xmode(i) = oo_.posterior_mode.(f).(nm);",
  "  fprintf(fid, '%s,%.12g,%.12g,%.12g,%.12g\\n', nm, xmode(i), oo_.posterior_mean.(f).(nm), oo_.posterior_hpdinf.(f).(nm), oo_.posterior_hpdsup.(f).(nm));",
  "end",
  "fclose(fid);",
  "if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end",
  "bnds = prior_bounds(bayestopt_, options_.prior_trunc);",
  "fid = fopen('dsge_mode.csv'); c = textscan(fid, '%s %f', 'Delimiter', ',', 'HeaderLines', 1); fclose(fid);",
  "xdsge = zeros(n, 1);",
  "for i = 1:n, j = find(strcmp(c{1}, names{i})); if numel(j) ~= 1, error('No dsge value for %s', names{i}); end, xdsge(i) = c{2}(j); end",
  "fid = fopen('dynare_kernel.csv', 'w'); fprintf(fid, 'point,loglik,logprior,logpost\\n');",
  "pts = {'dynare_mode', xmode; 'dsge_mode', xdsge};",
  "for k = 1:2",
  "  x = pts{k, 2};",
  "  fval = dsge_likelihood(x, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bnds, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);",
  "  lp = priordens(x, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, bayestopt_.p3, bayestopt_.p4);",
  "  fprintf(fid, '%s,%.10f,%.10f,%.10f\\n', pts{k, 1}, -fval - lp, lp, -fval);",
  "end",
  "fclose(fid);")

work <- Sys.getenv("DYNARE_WORK_DIR", "")
if (!nzchar(work) || !file.exists(file.path(work, "dynare_posterior.csv"))) {
  if (!nzchar(work)) work <- file.path(tempdir(), "dynare_sw")
  dir.create(work, showWarnings = FALSE, recursive = TRUE)
  files <- c("sw_defs.mod", "sw_model_def.mod", "sw_estimation_def.mod",
             "sw_table1.mod", "usmodel_data.mat")
  stopifnot(all(file.copy(file.path(ex_dir, files), work, overwrite = TRUE)))
  mod <- readLines(file.path(work, "sw_table1.mod"))
  drop <- grepl("optim=('Algorithm','active-set')", mod, fixed = TRUE)
  stopifnot(sum(drop) == 1L, sum(grepl("mode_compute=1,", mod, fixed = TRUE)) == 1L)
  mod <- sub("mode_compute=1,", "mode_compute=4,", mod[!drop], fixed = TRUE)
  writeLines(mod, file.path(work, "sw_table1.mod"))
  utils::write.csv(data.frame(name = sub("^sd_e\\.", "", names(dsge_mode)),
                              value = sprintf("%.17g", dsge_mode)),
                   file.path(work, "dsge_mode.csv"), row.names = FALSE,
                   quote = FALSE)
  writeLines(octave_driver, file.path(work, "run_book_sw.m"))
  cat("Running Dynare in", work, "...\n")
  old <- setwd(work)
  status <- system2("octave", c("--no-gui", "--quiet", "run_book_sw.m"),
                    stdout = "octave_run.log", stderr = "octave_run.log")
  setwd(old)
  if (status != 0L || !file.exists(file.path(work, "dynare_posterior.csv")))
    stop("Dynare failed; see ", file.path(work, "octave_run.log"))
}
dyn <- utils::read.csv(file.path(work, "dynare_posterior.csv"))
dyn_kernel <- utils::read.csv(file.path(work, "dynare_kernel.csv"))
for (f in c("octave_run.log", "sw_table1.log")) {
  if (file.exists(file.path(work, f)))
    file.copy(file.path(work, f), file.path(out, paste0("dynare_", f)),
              overwrite = TRUE)
}

# ---------------------------------------------------------------------------
# dsge's kernel at both modes, posterior summaries
# ---------------------------------------------------------------------------

m <- read_dynare(file.path(ex_dir, "sw_table1.mod"))
keys <- names(m$priors)
dyn$key <- ifelse(paste0("sd_e.", dyn$name) %in% keys, paste0("sd_e.", dyn$name),
                  dyn$name)
stopifnot(setequal(dyn$key, keys))
dyn <- dyn[match(keys, dyn$key), ]
dyn_mode <- stats::setNames(dyn$mode, keys)

raw <- R.matlab::readMat(file.path(ex_dir, "usmodel_data.mat"))
obs <- m$model$variables$observed
Y <- sapply(obs, function(v) as.numeric(raw[[v]]))[m$estimation$first_obs:230, ]
fixed <- unlist(m$model$fixed)
is_sd <- startsWith(keys, "sd_e.")
dsge_kernel <- function(theta) {
  theta <- theta[keys]
  sol <- solve_dsge(m$model, params = c(theta[!is_sd], fixed),
                    shock_sd = stats::setNames(theta[is_sd],
                      sub("sd_e.", "", keys[is_sd], fixed = TRUE)))
  ll <- dsge:::kalman_filter(sweep(Y, 2, sol$steady_state[obs]), sol$G, sol$H,
                             sol$M, sol$D, presample = m$estimation$presample,
                             init = m$model$kalman_init, loglik_only = TRUE)$loglik
  lp <- sum(vapply(keys, function(k) dsge:::dprior(m$priors[[k]], theta[[k]]), 0))
  c(loglik = ll, logprior = lp, logpost = ll + lp)
}

hpd <- function(x, prob = 0.90) {
  x <- sort(x); n <- length(x); k <- floor(prob * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)])
  c(x[i], x[i + k])
}
P <- readRDS(fit_path)$posterior
draws <- do.call(rbind, lapply(seq_len(dim(P)[3]), function(ch) P[, , ch]))[, keys]
post_sd <- apply(draws, 2, stats::sd)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

kt <- rbind(
  `Dynare mode, Dynare's kernel` = unlist(dyn_kernel[dyn_kernel$point == "dynare_mode", 2:4]),
  `Dynare mode, dsge's kernel` = dsge_kernel(dyn_mode),
  `dsge mode, Dynare's kernel` = unlist(dyn_kernel[dyn_kernel$point == "dsge_mode", 2:4]),
  `dsge mode, dsge's kernel` = dsge_kernel(dsge_mode))
cat("Posterior kernel of the book's model (1966Q1-2004Q4):\n")
print(round(kt, 6))
utils::write.csv(data.frame(point = rownames(kt), kt, row.names = NULL),
                 file.path(out, "dynare_kernel_comparison.csv"), row.names = FALSE)

tab <- data.frame(
  parameter = keys,
  dynare_mode = dyn$mode,
  dsge_mode = unname(dsge_mode[keys]),
  mode_diff_sd = (unname(dsge_mode[keys]) - dyn$mode) / post_sd,
  dynare_mean = dyn$mean,
  dsge_mean = unname(colMeans(draws)),
  dynare_hpd_inf = dyn$hpd_inf,
  dynare_hpd_sup = dyn$hpd_sup,
  dsge_hpd_inf = apply(draws, 2, function(x) hpd(x)[1]),
  dsge_hpd_sup = apply(draws, 2, function(x) hpd(x)[2]),
  posterior_sd = unname(post_sd),
  row.names = NULL)
tab$mean_diff_sd <- (tab$dsge_mean - tab$dynare_mean) / tab$posterior_sd
cat("\nDynare (book's file, mode_compute = 4, 2 x 2000 MH draws) vs dsge:\n")
print(format(tab[, c("parameter", "dynare_mode", "dsge_mode", "mode_diff_sd",
                     "dynare_mean", "dsge_mean", "mean_diff_sd")], digits = 3),
      row.names = FALSE)
utils::write.csv(tab, file.path(out, "dynare_vs_dsge.csv"), row.names = FALSE)
cat(sprintf("\n|mode difference| / posterior sd: median %.3f, max %.3f (%s)\n",
            stats::median(abs(tab$mode_diff_sd)), max(abs(tab$mode_diff_sd)),
            keys[which.max(abs(tab$mode_diff_sd))]))
cat(sprintf("|mean difference| / posterior sd: median %.3f, max %.3f (%s)\n",
            stats::median(abs(tab$mean_diff_sd)), max(abs(tab$mean_diff_sd)),
            keys[which.max(abs(tab$mean_diff_sd))]))
