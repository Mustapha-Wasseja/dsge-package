# End-to-end check of read_dynare() on Smets & Wouters (2007) against Dynare.
#
# Uses Johannes Pfeifer's replication file (GPL, not included here):
#   git clone https://github.com/JohannesPfeifer/DSGE_mod
# Usage (from the package root):
#   Rscript dev/dynare-validation/extra/validate_smets_wouters.R \
#           path/to/DSGE_mod/Smets_Wouters_2007
#
# Dynare runs the file's estimation command at the published posterior mode
# (mode_file = usmodel_mode, mode_compute = 0, mh_replic = 0), with
# lik_init = 1 because dsge initialises the Kalman filter at the stationary
# distribution. It exports the parameters, shock standard deviations, the
# data, the log-likelihood and log-prior at the mode and all first-order
# IRFs; the same quantities are then computed with dsge.

suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
sw_dir <- if (length(args) > 0L) args[1] else stop("Give the SW directory.")
dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")
mod <- file.path(sw_dir, "Smets_Wouters_2007.mod")

work <- tempfile("sw")
dir.create(work)
invisible(file.copy(file.path(sw_dir, c("usmodel_data.mat", "usmodel_mode.mat")),
                     work))
src <- readLines(mod)
src <- src[!grepl("^shock_decomposition", src)]
src <- sub("lik_init=2", "lik_init=1", src, fixed = TRUE)
src <- sub(", tex);", ");", src, fixed = TRUE)
writeLines(c(src, "stoch_simul(order = 1, irf = 20, nograph, noprint);"),
           file.path(work, "sw_est.mod"))
writeLines(c(
  sprintf("addpath('%s');", dynare_path),
  "dynare sw_est noclearall nolog",
  "fid = fopen('params.csv', 'w');",
  "for i = 1:M_.param_nbr, fprintf(fid, '%s,%.17g\\n', M_.param_names{i}, M_.params(i)); end",
  "fclose(fid);",
  "fid = fopen('sd.csv', 'w');",
  "for i = 1:M_.exo_nbr, fprintf(fid, '%s,%.17g\\n', M_.exo_names{i}, sqrt(M_.Sigma_e(i, i))); end",
  "fclose(fid);",
  "d = load('usmodel_data.mat');",
  "obs = {'dy','dc','dinve','labobs','pinfobs','dw','robs'};",
  "fid = fopen('data.csv', 'w'); fprintf(fid, '%s,', obs{1:end-1}); fprintf(fid, '%s\\n', obs{end});",
  "for t = 1:numel(d.dy)",
  "  for j = 1:numel(obs), v = d.(obs{j}); fprintf(fid, '%.17g', v(t)); if j < numel(obs), fprintf(fid, ','); end, end",
  "  fprintf(fid, '\\n');",
  "end",
  "fclose(fid);",
  "if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end",
  "xparam1 = get_all_parameters(estim_params_, M_);",
  "fval = dsge_likelihood(xparam1, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, prior_bounds(bayestopt_, options_.prior_trunc), oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);",
  "lnprior = priordens(xparam1, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, bayestopt_.p3, bayestopt_.p4);",
  "fid = fopen('lik.csv', 'w'); fprintf(fid, '%.12f,%.12f\\n', -fval - lnprior, lnprior); fclose(fid);",
  "fid = fopen('irfs.csv', 'w'); f = fieldnames(oo_.irfs);",
  "for j = 1:numel(f), fprintf(fid, '%s', f{j}); fprintf(fid, ',%.15g', oo_.irfs.(f{j})); fprintf(fid, '\\n'); end",
  "fclose(fid);"
), file.path(work, "run.m"))
old <- setwd(work)
system2("octave", c("--no-gui", "--quiet", "run.m"), stdout = FALSE,
        stderr = FALSE)
setwd(old)

m <- read_dynare(mod)
pr <- utils::read.csv(file.path(work, "params.csv"), header = FALSE)
mode <- stats::setNames(pr$V2, pr$V1)
sdv <- utils::read.csv(file.path(work, "sd.csv"), header = FALSE)
sd <- stats::setNames(sdv$V2, sdv$V1)
lik <- scan(file.path(work, "lik.csv"), sep = ",", quiet = TRUE)

sol <- solve_dsge(m, params = mode[m$model$parameters], shock_sd = sd)

lines <- readLines(file.path(work, "irfs.csv"))
dyn <- lapply(strsplit(lines, ","), function(z) as.numeric(z[-1]))
names(dyn) <- vapply(strsplit(lines, ","), `[`, "", 1)
ir <- irf(sol, periods = 19, se = FALSE)$data
worst <- 0
for (v in m$variables) for (e in m$shocks) {
  ours <- ir$value[ir$response == v & ir$impulse == e]
  th <- dyn[[paste0(v, "_", e)]]
  if (is.null(th)) th <- 0 * ours
  worst <- max(worst, max(abs(ours - th)))
}

dat <- utils::read.csv(file.path(work, "data.csv"))
obs <- m$model$variables$observed
y <- sweep(as.matrix(dat[, obs]), 2, sol$steady_state[obs])
ll <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D,
                           presample = m$estimation$presample)$loglik
lp <- sum(vapply(names(m$priors), function(nm) {
  x <- if (startsWith(nm, "sd_e.")) sd[[sub("sd_e.", "", nm, fixed = TRUE)]]
       else mode[[nm]]
  dsge:::dprior(m$priors[[nm]], x)
}, 0))

res <- data.frame(
  check = c("IRFs (40 variables x 7 shocks), max abs diff",
            "log-likelihood at posterior mode (presample 4)",
            "log-prior at posterior mode",
            "log-posterior kernel at posterior mode"),
  dynare = c("", sprintf("%.9f", lik[1]), sprintf("%.9f", lik[2]),
             sprintf("%.9f", lik[1] + lik[2])),
  dsge = c(sprintf("%.3g", worst), sprintf("%.9f", ll), sprintf("%.9f", lp),
           sprintf("%.9f", ll + lp)))
print(res, row.names = FALSE)
