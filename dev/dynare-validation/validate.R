# Validate read_dynare() against Dynare itself.
#
# For each .mod file: run Dynare (via Octave) for first-order IRFs of every
# endogenous variable to every shock, then import the same file with
# read_dynare(), solve with dsge, and compare the IRFs.
#
# Requires Octave and Dynare (e.g. `apt-get install octave dynare`).
# Usage (from the package root):
#   Rscript dev/dynare-validation/validate.R [file.mod ...]

suppressMessages(devtools::load_all(quiet = TRUE))

dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")
horizon <- 20L

args <- commandArgs(trailingOnly = TRUE)
mods <- if (length(args) > 0L) args else c(
  "inst/examples/rbc.mod",
  list.files("dev/dynare-validation", "\\.mod$", full.names = TRUE)
)

run_dynare_irfs <- function(mod) {
  work <- tempfile("dyn")
  dir.create(work)
  base <- gsub("[^A-Za-z0-9_]", "_", tools::file_path_sans_ext(basename(mod)))
  src <- readLines(mod)
  # Replace any stoch_simul with our own, covering all variables
  src <- src[!grepl("^\\s*stoch_simul", src)]
  if (!any(grepl("^\\s*discretionary_policy", src))) {
    src <- c(src, sprintf("stoch_simul(order = 1, irf = %d, nograph, noprint);",
                          horizon))
  }
  writeLines(src, file.path(work, paste0(base, ".mod")))
  writeLines(c(
    sprintf("addpath('%s');", dynare_path),
    sprintf("dynare %s noclearall nolog", base),
    "fid = fopen('irfs.csv', 'w');",
    "f = fieldnames(oo_.irfs);",
    "for j = 1:numel(f)",
    "  fprintf(fid, '%s', f{j}); fprintf(fid, ',%.15g', oo_.irfs.(f{j}));",
    "  fprintf(fid, '\\n');",
    "end",
    "fclose(fid);"
  ), file.path(work, "run.m"))
  old <- setwd(work)
  on.exit(setwd(old))
  status <- system2("octave", c("--no-gui", "--quiet", "run.m"),
                    stdout = "octave.log", stderr = "octave.log")
  if (!file.exists("irfs.csv")) {
    stop("Dynare failed for ", mod, "; see ", file.path(work, "octave.log"))
  }
  lines <- readLines("irfs.csv")
  vals <- lapply(strsplit(lines, ","), function(x) as.numeric(x[-1]))
  names(vals) <- vapply(strsplit(lines, ","), `[`, "", 1)
  vals
}

compare <- function(mod) {
  dyn <- run_dynare_irfs(mod)
  imp <- read_dynare(mod)
  sol <- solve_dsge(imp)
  ir <- irf(sol, periods = horizon - 1L, se = FALSE)$data
  worst <- 0
  n <- 0L
  for (v in imp$variables) {
    for (e in imp$shocks) {
      key <- paste0(v, "_", e)
      ours <- ir$value[ir$response == v & ir$impulse == e]
      theirs <- dyn[[key]]
      if (is.null(theirs)) {
        # Dynare drops IRFs that are numerically zero
        theirs <- numeric(length(ours))
      }
      worst <- max(worst, max(abs(ours - theirs)))
      n <- n + 1L
    }
  }
  data.frame(model = basename(mod), irfs_compared = n,
             max_abs_diff = signif(worst, 3))
}

res <- do.call(rbind, lapply(mods, compare))
print(res, row.names = FALSE)
