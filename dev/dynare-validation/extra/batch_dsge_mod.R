# Batch test of read_dynare() on a collection of real .mod files.
#
# For every .mod file under ROOT (e.g. a clone of
# https://github.com/JohannesPfeifer/DSGE_mod), the file is truncated at its
# first computational command (stoch_simul, estimation, ...), so Dynare and
# dsge see exactly the same model and calibration. Dynare computes first-
# order IRFs; read_dynare() + solve_dsge() must reproduce them.
#
# Usage: Rscript dev/dynare-validation/extra/batch_dsge_mod.R ROOT [OUT.csv] [LOGDIR]

suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
out_csv <- if (length(args) > 1L) args[2] else "batch_results.csv"
dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")
horizon <- 20L
log_dir <- if (length(args) > 2L) args[3] else ""
if (nzchar(log_dir)) dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "mod_tools.R"))

run_one <- function(mod) {
  res <- data.frame(model = sub(paste0("^", root, "/?"), "", mod),
                    status = NA_character_, irfs = NA_integer_,
                    max_abs_diff = NA_real_, scale = NA_real_,
                    message = "", stringsAsFactors = FALSE)
  if (is_perfect_foresight(mod)) {
    res$status <- "skipped"
    res$message <- "perfect foresight model"
    return(res)
  }
  prep <- prepare_mod(mod, "bt", order = 1L, irf = horizon)
  work <- prep$dir
  base <- sub("_bt$", "", prep$name)
  test_mod <- prep$mod
  writeLines(c(
    sprintf("addpath('%s');", dynare_path),
    sprintf("dynare %s_bt noclearall nolog", base),
    "% shocks without a variance get sd 0.01 (in dsge too), so IRFs exist",
    "z = find(diag(M_.Sigma_e) == 0);",
    "fid = fopen('sd_override.csv', 'w'); fprintf(fid, '%s\\n', M_.exo_names{z}); fclose(fid);",
    "if ~isempty(z)",
    "  for k = z', M_.Sigma_e(k, k) = 1e-4; end",
    "  [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, {});",
    "end",
    "fid = fopen('irfs.csv', 'w');",
    "if isfield(oo_, 'irfs'), f = fieldnames(oo_.irfs); else, f = {}; end",
    "for j = 1:numel(f)",
    "  fprintf(fid, '%s', f{j}); fprintf(fid, ',%.15g', oo_.irfs.(f{j}));",
    "  fprintf(fid, '\\n');",
    "end",
    "fclose(fid);"), file.path(work, "run.m"))
  old <- setwd(work)
  on.exit(setwd(old))
  system2("timeout", c("600", "octave", "--no-gui", "--quiet", "run.m"),
          stdout = "octave.log", stderr = "octave.log")
  if (nzchar(log_dir)) file.copy("octave.log", file.path(log_dir, paste0(base, ".log")),
                                 overwrite = TRUE)
  if (!file.exists("irfs.csv") || file.size("irfs.csv") == 0) {
    log <- readLines("octave.log", warn = FALSE)
    err <- grep("error|ERROR", log, value = TRUE)
    if (length(err) == 0L) err <- utils::tail(log[nzchar(log)], 2)
    res$status <- "dynare_failed"
    res$message <- substr(paste(utils::head(err, 2), collapse = " | "), 1, 200)
    return(res)
  }
  lines <- readLines("irfs.csv")
  dyn <- lapply(strsplit(lines, ","), function(z) as.numeric(z[-1]))
  names(dyn) <- vapply(strsplit(lines, ","), `[`, "", 1)

  imp <- tryCatch(read_dynare(test_mod), error = function(e) e)
  if (inherits(imp, "error")) {
    res$status <- "import_failed"
    res$message <- substr(conditionMessage(imp), 1, 200)
    return(res)
  }
  sol <- tryCatch({
    setTimeLimit(elapsed = 600)
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
    sd <- imp$shock_sd
    ov <- if (file.exists("sd_override.csv")) readLines("sd_override.csv") else character(0)
    sd[intersect(ov[nzchar(ov)], names(sd))] <- 0.01
    solve_dsge(imp, shock_sd = sd)
  }, error = function(e) e)
  if (inherits(sol, "error") || !isTRUE(sol$stable)) {
    res$status <- "solve_failed"
    res$message <- if (inherits(sol, "error")) {
      substr(conditionMessage(sol), 1, 200)
    } else "no stable solution"
    return(res)
  }
  ir <- irf(sol, periods = horizon - 1L, se = FALSE)$data
  worst <- 0
  scale <- 0
  n <- 0L
  for (v in imp$variables) for (e in imp$shocks) {
    ours <- ir$value[ir$response == v & ir$impulse == e]
    th <- dyn[[paste0(v, "_", e)]]
    if (is.null(th)) th <- 0 * ours
    if (nzchar(Sys.getenv("DEBUG_BATCH")) && max(abs(ours - th)) > 1e-6) {
      cat(sprintf("  %s <- %s: max diff %.3g (dynare max %.3g)\n", v, e,
                  max(abs(ours - th)), max(abs(th))))
    }
    worst <- max(worst, max(abs(ours - th)))
    scale <- max(scale, max(abs(th)))
    n <- n + 1L
  }
  res$irfs <- n
  res$max_abs_diff <- signif(worst, 3)
  res$scale <- signif(scale, 3)
  res$status <- if (worst <= 1e-6 * max(1, scale)) "match" else "mismatch"
  res
}

mods <- sort(list.files(root, "\\.mod$", recursive = TRUE, full.names = TRUE))
results <- list()
for (m in mods) {
  r <- tryCatch(run_one(m), error = function(e) {
    data.frame(model = m, status = "harness_error", irfs = NA, max_abs_diff = NA,
               scale = NA, message = substr(conditionMessage(e), 1, 200))
  })
  results[[length(results) + 1L]] <- r
  cat(sprintf("%-60s %-14s %s %s\n", r$model, r$status,
              format(r$max_abs_diff), r$message))
  utils::write.csv(do.call(rbind, results), out_csv, row.names = FALSE)
}
tab <- do.call(rbind, results)
print(table(tab$status))
