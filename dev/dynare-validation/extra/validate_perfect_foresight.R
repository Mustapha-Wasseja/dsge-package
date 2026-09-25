# Perfect-foresight simulations: simulate_perfect_foresight() versus Dynare.
#
# Each file is cut after its first perfect_foresight_solver (or simul)
# command. Dynare runs it and exports oo_.endo_simul; read_dynare() imports
# the same cut file and simulate_perfect_foresight() computes the path. The
# largest absolute difference over all variables and periods is reported,
# relative to the largest deviation of the path from its terminal value.
#
# Usage: Rscript dev/dynare-validation/extra/validate_perfect_foresight.R file.mod ...

suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")

cut_after_solver <- function(lines) {
  hit <- grep("^\\s*(perfect_foresight_solver|simul)\\b", lines)
  if (length(hit) == 0L) return(NULL)
  j <- hit[1]
  while (!grepl(";", lines[j]) && j < length(lines)) j <- j + 1L
  # a second solver call right after the first (e.g. retried with lmmcp)
  nxt <- j + 1L
  while (nxt <= length(lines) && !nzchar(trimws(sub("%.*$", "", lines[nxt])))) nxt <- nxt + 1L
  if (nxt <= length(lines) && grepl("^\\s*perfect_foresight_solver\\b", lines[nxt])) {
    j <- nxt
    while (!grepl(";", lines[j]) && j < length(lines)) j <- j + 1L
  }
  kept <- lines[seq_len(j)]
  # ask Dynare for a tight solution (its default tolerance is 1e-5)
  last <- gsub("\\b(tolf|tolx)\\s*=\\s*[^,)]+,?\\s*", "", kept[j])
  last <- sub(",\\s*\\)", ")", sub("\\(\\s*\\)", "", last))
  if (grepl("perfect_foresight_solver\\s*\\(", last)) {
    kept[j] <- sub("perfect_foresight_solver\\s*\\(",
                   "perfect_foresight_solver(tolf = 1e-12, tolx = 1e-12, ", last)
  } else {
    kept[j] <- sub("(perfect_foresight_solver|simul)\\s*;",
                   "\\1(tolf = 1e-12, tolx = 1e-12);", last)
  }
  stack <- character(0)
  for (l in kept) {
    if (grepl("^\\s*@#\\s*(if|ifdef|ifndef)\\b", l)) stack <- c(stack, "@#endif")
    if (grepl("^\\s*@#\\s*for\\b", l)) stack <- c(stack, "@#endfor")
    if (grepl("^\\s*@#\\s*(endif|endfor)\\b", l)) stack <- stack[-length(stack)]
  }
  c(kept, rev(stack))
}

run_one <- function(mod) {
  work <- tempfile("pf")
  dir.create(work)
  file.copy(list.files(dirname(mod), full.names = TRUE), work, recursive = TRUE)
  base <- paste0(gsub("[^A-Za-z0-9_]", "_",
                      tools::file_path_sans_ext(basename(mod))), "_pf")
  # expand macros first, so the cut is made in the branch actually used
  raw <- paste(readLines(mod, warn = FALSE), collapse = "\n")
  if (grepl("@#", raw)) {
    raw <- dyn_macro_expand(dyn_strip_comments(raw), dir = dirname(mod),
                            defines = list())
  }
  src <- cut_after_solver(strsplit(raw, "\n", fixed = TRUE)[[1]])
  if (is.null(src)) return(list(status = "no perfect-foresight command"))
  writeLines(src, file.path(work, paste0(base, ".mod")))
  ssm <- sub("\\.mod$", "_steadystate.m", mod)
  if (file.exists(ssm)) file.copy(ssm, file.path(work, paste0(base, "_steadystate.m")))
  writeLines(c(
    sprintf("addpath('%s');", dynare_path),
    sprintf("dynare %s noclearall nolog", base),
    "dlmwrite('endo_simul.csv', oo_.endo_simul, 'precision', '%.15g');",
    "fid = fopen('names.csv', 'w');",
    "for i = 1:M_.orig_endo_nbr, fprintf(fid, '%s\\n', M_.endo_names{i}); end; fclose(fid);"
  ), file.path(work, "run.m"))
  old <- setwd(work)
  on.exit(setwd(old))
  system2("timeout", c("900", "octave", "--no-gui", "--quiet", "run.m"),
          stdout = "octave.log", stderr = "octave.log")
  if (!file.exists("endo_simul.csv")) {
    err <- grep("error", readLines("octave.log", warn = FALSE), value = TRUE)
    return(list(status = "dynare_failed", msg = paste(utils::head(err, 1))))
  }
  dyn <- as.matrix(utils::read.csv("endo_simul.csv", header = FALSE))
  nm <- readLines("names.csv")
  dyn <- dyn[seq_along(nm), , drop = FALSE]
  rownames(dyn) <- nm
  imp <- tryCatch(read_dynare(file.path(work, paste0(base, ".mod"))),
                  error = function(e) e)
  if (inherits(imp, "error")) return(list(status = "import_failed",
                                          msg = conditionMessage(imp)))
  pf <- tryCatch(simulate_perfect_foresight(imp), error = function(e) e)
  if (inherits(pf, "error")) return(list(status = "solve_failed",
                                         msg = conditionMessage(pf)))
  vars <- intersect(colnames(pf$path), nm)
  ours <- t(pf$path[, vars, drop = FALSE])
  theirs <- dyn[vars, , drop = FALSE]
  if (ncol(ours) != ncol(theirs)) {
    return(list(status = "mismatch", msg = sprintf(
      "periods: dsge %d, Dynare %d", ncol(ours), ncol(theirs))))
  }
  if (nzchar(Sys.getenv("DEBUG_PF"))) {
    d <- apply(abs(ours - theirs), 1, max)
    print(signif(d[d > 1e-6], 3))
    w <- names(which.max(d))
    print(rbind(dsge = round(ours[w, 1:min(12, ncol(ours))], 5),
                dynare = round(theirs[w, 1:min(12, ncol(ours))], 5)))
  }
  diff <- max(abs(ours - theirs))
  scale <- max(abs(theirs - theirs[, ncol(theirs)]))
  list(status = if (diff <= 1e-6 * max(1, scale)) "match" else "mismatch",
       diff = diff, scale = scale, vars = length(vars),
       periods = pf$periods, converged = pf$converged)
}

out <- list()
for (mod in args) {
  r <- tryCatch(run_one(mod), error = function(e) list(status = "harness_error",
                                                     msg = conditionMessage(e)))
  row <- data.frame(model = basename(mod), status = r$status,
                    vars = if (is.null(r$vars)) NA else r$vars,
                    periods = if (is.null(r$periods)) NA else r$periods,
                    max_abs_diff = if (is.null(r$diff)) NA else signif(r$diff, 3),
                    max_deviation = if (is.null(r$scale)) NA else signif(r$scale, 3),
                    message = if (is.null(r$msg)) "" else substr(r$msg, 1, 90))
  print(row, row.names = FALSE)
  out[[length(out) + 1L]] <- row
}
print(do.call(rbind, out), row.names = FALSE)
