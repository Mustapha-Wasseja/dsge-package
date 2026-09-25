# Full Bayesian estimation: dsge versus Dynare on the same data.
#
# Uses Johannes Pfeifer's RBC_baseline_first_diff_bayesian.mod (GPL, not
# included; from https://github.com/JohannesPfeifer/DSGE_mod): a nonlinear
# RBC model estimated on simulated output and consumption growth, with two
# autocorrelations (beta priors) and two shock standard deviations
# (inverse gamma priors). Dynare simulates the data (stoch_simul with
# periods = 200, fixed seed), finds the posterior mode (csminwel) and runs
# two Metropolis-Hastings chains. dsge imports the unchanged file, uses the
# same data and priors, maximises the same posterior and runs its own
# sampler (bayes_dsge()). Compared: the log posterior kernel at Dynare's
# mode, the posterior mode, and posterior means and standard deviations.
#
# Usage: Rscript dev/dynare-validation/extra/validate_estimation.R \
#          path/to/DSGE_mod/RBC_baseline [MH_DRAWS]

suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
src_dir <- args[1]
mh <- if (length(args) > 1L) as.integer(args[2]) else 20000L
dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")
mod <- file.path(src_dir, "RBC_baseline_first_diff_bayesian.mod")

work <- Sys.getenv("EST_WORK", tempfile("est"))
dir.create(work, showWarnings = FALSE)
lines <- readLines(mod, warn = FALSE)
cut <- grep("^estimation\\(", lines)
tail_end <- cut + grep(";", lines[cut:length(lines)])[1] - 1L
est_cmd <- sprintf(paste0(
  "estimation(datafile=first_diff_data, mode_compute=4, mh_replic=%d, ",
  "mh_nblocks=2, mh_jscale=1.5, nograph, nodisplay, nodiagnostic);"), mh)
writeLines(c(lines[seq_len(cut - 1L)], est_cmd), file.path(work, "rbc_est.mod"))
writeLines(c(
  sprintf("addpath('%s');", dynare_path),
  "dynare rbc_est noclearall nolog",
  "fid = fopen('data.csv', 'w'); fprintf(fid, 'g_obs,c_obs\\n');",
  "d = dataset_.data; for t = 1:size(d, 1), fprintf(fid, '%.17g,%.17g\\n', d(t, 1), d(t, 2)); end; fclose(fid);",
  "names = bayestopt_.name;",
  "fid = fopen('post.csv', 'w'); fprintf(fid, 'name,mode,mean,sd,hpd_lo,hpd_hi\\n');",
  "for i = 1:numel(names)",
  "  nm = names{i}; key = strrep(nm, 'SE_', '');",
  "  if any(strcmp(fieldnames(oo_.posterior_mean.parameters), key))",
  "    mo = oo_.posterior_mode.parameters.(key); me = oo_.posterior_mean.parameters.(key);",
  "    sd = oo_.posterior_std.parameters.(key); lo = oo_.posterior_hpdinf.parameters.(key); hi = oo_.posterior_hpdsup.parameters.(key);",
  "  else",
  "    mo = oo_.posterior_mode.shocks_std.(key); me = oo_.posterior_mean.shocks_std.(key);",
  "    sd = oo_.posterior_std.shocks_std.(key); lo = oo_.posterior_hpdinf.shocks_std.(key); hi = oo_.posterior_hpdsup.shocks_std.(key);",
  "  end",
  "  fprintf(fid, '%s,%.12g,%.12g,%.12g,%.12g,%.12g\\n', key, mo, me, sd, lo, hi);",
  "end; fclose(fid);",
  "xparam1 = get_posterior_parameters('mode', M_, estim_params_, oo_, options_);",
  "if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end",
  "fval = dsge_likelihood(xparam1, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, prior_bounds(bayestopt_, options_.prior_trunc), oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);",
  "lnprior = priordens(xparam1, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, bayestopt_.p3, bayestopt_.p4);",
  "fid = fopen('kernel.csv', 'w'); fprintf(fid, '%.12f,%.12f\\n', -fval - lnprior, lnprior); fclose(fid);",
  "fid = fopen('mdd.csv', 'w'); fprintf(fid, '%.8f,%.8f\\n', oo_.MarginalDensity.LaplaceApproximation, oo_.MarginalDensity.ModifiedHarmonicMean); fclose(fid);"
), file.path(work, "run.m"))
# DYNARE_WORK=<dir> reuses the output of an earlier Dynare run
reuse <- Sys.getenv("DYNARE_WORK")
t_dynare <- NA_real_
if (nzchar(reuse)) {
  work <- reuse
} else {
  old <- setwd(work)
  t0 <- Sys.time()
  system2("octave", c("--no-gui", "--quiet", "run.m"), stdout = "octave.log",
          stderr = "octave.log")
  t_dynare <- as.numeric(Sys.time() - t0, units = "mins")
  setwd(old)
}
if (!file.exists(file.path(work, "post.csv"))) {
  cat(utils::tail(readLines(file.path(work, "octave.log")), 30), sep = "\n")
  stop("Dynare failed")
}
post <- utils::read.csv(file.path(work, "post.csv"))
kern <- scan(file.path(work, "kernel.csv"), sep = ",", quiet = TRUE)
mdd <- scan(file.path(work, "mdd.csv"), sep = ",", quiet = TRUE)
dat <- utils::read.csv(file.path(work, "data.csv"))

m <- read_dynare(mod)
params <- setdiff(post$name, m$shocks)
shocks <- intersect(post$name, m$shocks)
theta_names <- c(params, paste0("sd_e.", shocks))

# log posterior kernel (natural parameterisation, as in Dynare)
kernel <- function(theta) {
  names(theta) <- theta_names
  lp <- sum(vapply(theta_names, function(nm) dprior(m$priors[[nm]], theta[[nm]]), 0))
  if (!is.finite(lp)) return(-Inf)
  sd <- m$shock_sd
  sd[shocks] <- theta[paste0("sd_e.", shocks)]
  pv <- c(theta[params], unlist(m$model$fixed))
  sol <- tryCatch(solve_dsge(m, params = theta[params], shock_sd = sd),
                  error = function(e) NULL)
  if (is.null(sol) || !isTRUE(sol$stable)) return(-Inf)
  obs <- m$model$variables$observed
  y <- sweep(as.matrix(dat[, obs]), 2, sol$steady_state[obs])
  ll <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D,
                      presample = m$estimation$presample,
                      init = m$model$kalman_init)$loglik
  ll + lp
}
dyn_key <- ifelse(post$name %in% shocks, paste0("sd_e.", post$name), post$name)
dyn_mode <- stats::setNames(post$mode, dyn_key)[theta_names]
k_dsge <- kernel(dyn_mode)

# posterior mode, started from the prior means (not Dynare's mode); rho's
# on the logit scale and standard deviations on the log scale
is_rho <- theta_names %in% params
to_nat <- function(z) ifelse(is_rho, stats::plogis(z), exp(z))
start <- vapply(theta_names, function(nm) {
  p <- m$priors[[nm]]
  if (!is.null(p$mean)) p$mean else mean(stats::na.omit(dyn_mode[nm]))
}, 0)
z0 <- ifelse(is_rho, stats::qlogis(start), log(start))
negk <- function(z) {
  v <- -kernel(to_nat(z))
  if (is.finite(v)) v else 1e10
}
opt <- stats::optim(z0, negk, method = "Nelder-Mead",
                    control = list(maxit = 4000, reltol = 1e-12))
opt <- stats::optim(opt$par, negk, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-14))
opt$par <- to_nat(opt$par)
mode_dsge <- stats::setNames(opt$par, theta_names)

t0 <- Sys.time()
fit_file <- file.path(work, "dsge_fit.rds")
fit <- if (nzchar(reuse) && file.exists(fit_file)) readRDS(fit_file) else bayes_dsge(m, data = dat, chains = 2L, iter = mh + mh %/% 4L,
                  warmup = mh %/% 4L, seed = 1, n_cores = 2L)
t_dsge <- as.numeric(Sys.time() - t0, units = "mins")
saveRDS(fit, file.path(work, "dsge_fit.rds"))
# posterior: draws x parameters x chains
draws <- do.call(rbind, lapply(seq_len(dim(fit$posterior)[3]), function(ch) {
  fit$posterior[, , ch]
}))
colnames(draws) <- dimnames(fit$posterior)[[2]]

tab <- data.frame(
  parameter = theta_names,
  mode_dynare = signif(dyn_mode, 6), mode_dsge = signif(mode_dsge, 6),
  mean_dynare = signif(post$mean[match(sub("^sd_e\\.", "", theta_names), post$name)], 4),
  mean_dsge = signif(colMeans(draws[, theta_names, drop = FALSE]), 4),
  sd_dynare = signif(post$sd[match(sub("^sd_e\\.", "", theta_names), post$name)], 3),
  sd_dsge = signif(apply(draws[, theta_names, drop = FALSE], 2, stats::sd), 3),
  row.names = NULL)
print(tab, row.names = FALSE)
cat(sprintf("\nLog posterior kernel at Dynare's mode: Dynare %.9f, dsge %.9f\n",
            kern[1] + kern[2], k_dsge))
cat(sprintf("Log posterior kernel at dsge's mode: %.9f\n", -opt$value))
cat(sprintf("Dynare log marginal density: Laplace %.4f, modified harmonic mean %.4f\n",
            mdd[1], mdd[2]))
cat(sprintf("Draws per chain: %d (Dynare), %d after warmup (dsge); time: Dynare %.1f min, dsge %.1f min\n",
            mh, mh, t_dynare, t_dsge))
