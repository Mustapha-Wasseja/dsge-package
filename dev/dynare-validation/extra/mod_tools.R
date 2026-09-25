# Helpers shared by the batch and higher-order comparison scripts.

stop_cmds <- c("stoch_simul", "estimation", "shock_decomposition",
               "realtime_shock_decomposition", "identification", "osr",
               "perfect_foresight_setup", "perfect_foresight_solver",
               "simul", "forecast", "conditional_forecast", "calib_smoother",
               "dynare_sensitivity", "extended_path", "occbin_setup",
               "occbin_solver", "ramsey_policy", "discretionary_policy",
               "write_latex_dynamic_model", "write_latex_static_model",
               "write_latex_original_model", "write_latex_parameter_table",
               "write_latex_prior_table", "collect_latex_files",
               "evaluate_planner_objective", "model_diagnostics",
               "method_of_moments", "mode_check", "smoother2histval",
               "histval_file", "initval_file", "model_info", "save_params_and_steady_state",
               "load_params_and_steady_state", "planner_objective_value")

policy_cmds <- c("ramsey_policy", "discretionary_policy")

# Cut the file at its first computing command, close any macro @#if/@#for
# left open by the cut, and append a first-order stoch_simul (or the policy
# command itself, with IRF options, for ramsey_policy/discretionary_policy).
truncate_mod <- function(lines, order = 1L, irf = 20L) {
  pat <- paste0("^\\s*(", paste(stop_cmds, collapse = "|"), ")\\b")
  hit <- grep(pat, lines)
  in_block <- FALSE
  first <- NA
  for (i in seq_along(lines)) {
    l <- lines[i]
    if (!in_block && i %in% hit) { first <- i; break }
    opens <- lengths(regmatches(l, gregexpr("/\\*", l)))
    closes <- lengths(regmatches(l, gregexpr("\\*/", l)))
    if (opens > closes) in_block <- TRUE
    if (closes > opens) in_block <- FALSE
  }
  irf_opts <- sprintf("irf = %d, nograph, noprint, nomoments, nocorr", irf)
  tail_cmd <- sprintf("stoch_simul(order = %d, %s);", order, irf_opts)
  if (is.na(first)) {
    kept <- lines
  } else {
    kept <- lines[seq_len(first - 1L)]
    cmd <- sub("^\\s*([a-z_]+).*$", "\\1", lines[first])
    if (cmd %in% policy_cmds) {
      j <- first
      stmt <- lines[j]
      while (!grepl(";", stmt) && j < length(lines)) {
        j <- j + 1L
        stmt <- paste(stmt, lines[j])
      }
      stmt <- sub(";.*$", "", stmt)
      opts <- if (grepl("\\(", stmt)) sub("^[^(]*\\((.*)\\).*$", "\\1", stmt) else ""
      opts <- trimws(strsplit(opts, ",")[[1]])
      opts <- opts[nzchar(opts) & !grepl(
        "^(irf|order|nograph|noprint|nomoments|nocorr|graph_format|periods|hp_filter|tex|nodisplay|irf_shocks|conditional_variance_decomposition)\\b",
        opts)]
      tail_cmd <- sprintf("%s(%s);", cmd,
                          paste(c(opts, sprintf("order = %d", order), irf_opts),
                                collapse = ", "))
    }
  }
  # close macro blocks left open by the cut
  stack <- character(0)
  for (l in kept) {
    if (grepl("^\\s*@#\\s*(if|ifdef|ifndef)\\b", l)) stack <- c(stack, "@#endif")
    if (grepl("^\\s*@#\\s*for\\b", l)) stack <- c(stack, "@#endfor")
    if (grepl("^\\s*@#\\s*(endif|endfor)\\b", l)) stack <- stack[-length(stack)]
  }
  c(kept, rev(stack), tail_cmd)
}


# Copy a model's folder to a fresh directory with the truncated file as
# <base>_<suffix>.mod (and its steady-state file renamed to match).
prepare_mod <- function(mod, suffix, order = 1L, irf = 20L) {
  work <- tempfile("bm")
  dir.create(work)
  file.copy(list.files(dirname(mod), full.names = TRUE), work, recursive = TRUE)
  base <- gsub("[^A-Za-z0-9_]", "_", tools::file_path_sans_ext(basename(mod)))
  name <- paste0(base, "_", suffix)
  writeLines(truncate_mod(readLines(mod, warn = FALSE), order, irf),
             file.path(work, paste0(name, ".mod")))
  ssm <- file.path(dirname(mod), paste0(
    tools::file_path_sans_ext(basename(mod)), "_steadystate.m"))
  if (file.exists(ssm)) {
    file.copy(ssm, file.path(work, paste0(name, "_steadystate.m")))
  }
  list(dir = work, name = name, mod = file.path(work, paste0(name, ".mod")))
}

is_perfect_foresight <- function(mod) {
  any(grepl("^\\s*(simul|perfect_foresight_setup|perfect_foresight_solver)\\b",
            readLines(mod, warn = FALSE)))
}
