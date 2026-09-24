# Further checks of read_dynare() against Dynare (via Octave):
#   * Kalman-filter log-likelihood with a measurement error and with fewer
#     observables than shocks (Dynare: estimation with mode_compute = 0)
#   * optimal simple rule (Dynare: osr)
#   * OccBin piecewise-linear simulation (Dynare: occbin_solver), with one
#     and with two surprise shocks
#
# Usage (from the package root):
#   Rscript dev/dynare-validation/extra/validate_extra.R

suppressMessages(devtools::load_all(quiet = TRUE))
dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")
here <- "dev/dynare-validation/extra"

run_octave <- function(mod, extra_lines, files = character(0)) {
  work <- tempfile("dyn")
  dir.create(work)
  file.copy(mod, work)
  for (f in files) file.copy(f, work)
  base <- tools::file_path_sans_ext(basename(mod))
  writeLines(c(sprintf("addpath('%s');", dynare_path),
               sprintf("dynare %s noclearall nolog", base), extra_lines),
             file.path(work, "run.m"))
  old <- setwd(work)
  on.exit(setwd(old))
  out <- system2("octave", c("--no-gui", "--quiet", "run.m"),
                 stdout = TRUE, stderr = TRUE)
  list(dir = work, out = out)
}

results <- list()

# --- log-likelihood -------------------------------------------------------
set.seed(5)
dat <- data.frame(y = cumsum(rnorm(60)) * 0.01, pi = rnorm(60) * 0.005,
                  r = rnorm(60) * 0.004)
dat$y <- dat$y - mean(dat$y)
csv <- file.path(tempdir(), "data.csv")
utils::write.csv(dat, csv, row.names = FALSE)
for (f in c("nk_me.mod", "nk_fewobs.mod")) {
  mod <- file.path(here, f)
  r <- run_octave(mod, character(0), files = csv)
  line <- grep("Initial value of the log posterior", r$out, value = TRUE)
  dyn_ll <- as.numeric(sub(".*:\\s*", "", line))
  m <- read_dynare(mod)
  sol <- solve_dsge(m)
  d <- dsge:::dyn_map_data(m, dat)
  y <- as.matrix(d[, m$model$variables$observed])
  ours <- dsge:::kalman_filter(y, sol$G, sol$H, sol$M, sol$D)$loglik
  results[[length(results) + 1L]] <- data.frame(
    check = paste("log-likelihood,", f), dynare = sprintf("%.4f", dyn_ll),
    dsge = sprintf("%.6f", ours), abs_diff = signif(abs(ours - dyn_ll), 3))
}

# --- OSR ------------------------------------------------------------------
mod <- file.path(here, "nk_osr.mod")
r <- run_octave(mod, c(
  "fprintf('OSR %.12f %.12f %.15g\\n', oo_.osr.optim_params.phi_pi, ",
  "        oo_.osr.optim_params.phi_y, oo_.osr.objective_function);"))
vals <- as.numeric(strsplit(sub("^OSR ", "", grep("^OSR ", r$out,
                                                  value = TRUE)), " ")[[1]])
o <- osr(read_dynare(mod))
results[[length(results) + 1L]] <- data.frame(
  check = "OSR loss at optimum", dynare = sprintf("%.10g", vals[3]),
  dsge = sprintf("%.10g", o$loss), abs_diff = signif(abs(o$loss - vals[3]), 3))
results[[length(results) + 1L]] <- data.frame(
  check = "OSR phi_pi (phi_y at bound 2)", dynare = sprintf("%.6f", vals[1]),
  dsge = sprintf("%.6f", o$optimal[["phi_pi"]]),
  abs_diff = signif(abs(o$optimal[["phi_pi"]] - vals[1]), 3))

# --- OccBin ---------------------------------------------------------------
occ_lines <- c(
  "pw = oo_.occbin.simul.piecewise; names = cellstr(M_.endo_names);",
  "fid = fopen('occ.csv', 'w');",
  "for j = 1:numel(names)",
  "  fprintf(fid, '%s', names{j}); fprintf(fid, ',%.15g', pw(:, j));",
  "  fprintf(fid, '\\n');",
  "end",
  "fclose(fid);")
base_occ <- file.path(here, "nk_zlb_occbin.mod")
two <- file.path(tempdir(), "nk_zlb_two_shocks.mod")
writeLines(sub("var eg; periods 1; values -0.06;",
               "var eg; periods 1 6; values -0.04 -0.05;",
               readLines(base_occ), fixed = TRUE), two)
for (mod in c(base_occ, two)) {
  r <- run_octave(mod, occ_lines)
  lines <- readLines(file.path(r$dir, "occ.csv"))
  dyn <- lapply(strsplit(lines, ","), function(z) as.numeric(z[-1]))
  names(dyn) <- vapply(strsplit(lines, ","), `[`, "", 1)
  o <- simulate_occbin(read_dynare(mod), horizon = length(dyn[[1]]))
  diff <- max(vapply(names(dyn), function(v) max(abs(o$controls[, v] -
                                                      dyn[[v]])), 0))
  results[[length(results) + 1L]] <- data.frame(
    check = paste("OccBin path,", basename(mod)),
    dynare = paste(sum(o$binding[, 1]), "periods binding"),
    dsge = paste(sum(o$binding[, 1]), "periods binding"),
    abs_diff = signif(diff, 3))
}

print(do.call(rbind, results), row.names = FALSE)
