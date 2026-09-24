# Import Dynare .mod files
#
# Translates the declarations, model block, calibration, steady state,
# shocks and estimation blocks of a Dynare .mod file into a dsgenl_model
# plus the calibrated parameter values, shock standard deviations and
# priors needed to solve or estimate it.
#
# Dynare timing is mapped onto the dsge state/control split mechanically,
# in the same way Dynare itself handles long leads and lags:
#   * every Dynare endogenous variable becomes a control;
#   * a lag x(-k) becomes an auxiliary state x_lagk, with
#     x_lag1(+1) = x and x_lagj(+1) = x_lag(j-1);
#   * a lead x(+k), k >= 2, becomes x_lead(k-1)(+1), with auxiliary
#     controls x_lead1 = x(+1) and x_leadj = x_lead(j-1)(+1);
#   * every varexo e becomes an exogenous state with e(+1) = 0, so that
#     e_t is the period-t innovation and carries the shock.

#' Import a Dynare .mod File
#'
#' Reads a Dynare model file and translates it into a nonlinear dsge
#' model, together with the calibration, shock standard deviations and
#' priors it declares, so that the model can be solved, simulated or
#' estimated in R without Dynare, MATLAB or Octave.
#'
#' @param file Path to a `.mod` file.
#' @param text Alternatively, the model code as a character vector (one
#'   element per line, or a single string). Used instead of `file` when
#'   supplied.
#' @param observed Optional character vector of observed variables. Overrides
#'   the file's `varobs` declaration. dsge requires as many observed
#'   variables as shocks; see Details.
#'
#' @return An object of class `"dsge_dynare"`, a list with components:
#' \describe{
#'   \item{model}{The translated `dsgenl_model`. Parameters listed in
#'     `estimated_params` are free (with starting values); all other
#'     parameters are fixed at their calibrated values.}
#'   \item{params}{Named numeric vector of calibrated parameter values.}
#'   \item{shock_sd}{Named numeric vector of shock standard deviations from
#'     the `shocks` block (0 for shocks with no declared variance, as in
#'     Dynare).}
#'   \item{priors}{Named list of [prior()] objects translated from
#'     `estimated_params`, ready for [bayes_dsge()], or `NULL`.}
#'   \item{estimated_params}{Data frame describing each `estimated_params`
#'     entry and how it was translated.}
#'   \item{observed}{Observed variables used for the model.}
#'   \item{variables, shocks, parameters}{Names declared in the file.}
#'   \item{aux}{Data frame of auxiliary variables created for leads and
#'     lags.}
#'   \item{commands}{List of Dynare commands found in the file (such as
#'     `stoch_simul` or `estimation`), recorded but not executed.}
#'   \item{notes}{Character vector of translation notes, including
#'     anything that was ignored or approximated.}
#' }
#'
#' @details
#' **Supported:** `var`, `varexo`, `parameters`,
#' `predetermined_variables`, `varobs`, parameter assignments, the
#' `model` block (including `model(linear)`, equation tags and `#`
#' model-local variables), leads and lags of any length, `initval`,
#' `steady_state_model`, the variance part of the `shocks` block, and
#' `estimated_params` / `estimated_params_init`.
#'
#' **Timing.** No manual re-timing is needed. Every Dynare variable
#' becomes a control; lags become auxiliary state variables named
#' `x_lag1`, `x_lag2`, ...; leads beyond one period become auxiliary
#' controls `x_lead1`, ...; and each shock becomes an exogenous state that
#' holds the current innovation. Impulse responses therefore have the
#' same timing as in Dynare. The auxiliary variables also appear in
#' solution and IRF output.
#'
#' **Observed variables.** dsge requires the number of observed variables
#' to equal the number of shocks. If `varobs` (or `observed`) has that
#' many variables it is used as is; otherwise, and when no `varobs` is
#' given, the first declared endogenous variables are used and a note is
#' recorded. The choice only matters for estimation.
#'
#' **Priors.** Dynare's mean/standard-deviation prior parameterisation is
#' converted to dsge's: `normal_pdf`, `beta_pdf`, `gamma_pdf`,
#' `uniform_pdf` and `inv_gamma2_pdf` are translated exactly.
#' `inv_gamma_pdf` / `inv_gamma1_pdf` (a prior on a standard deviation) has
#' no exact dsge counterpart and is approximated by an inverse gamma with
#' the same mean and standard deviation (an infinite standard deviation
#' maps to shape 2). Shifted or generalised priors and `weibull_pdf` are
#' not translated and are listed in `notes`.
#'
#' **Not supported:** macro-processor directives (`@#`; expand them first
#' with `dynare model.mod savemacro onlymacro` and import the generated
#' `model-macroexp.mod`), `varexo_det`, `external_function`,
#' `STEADY_STATE()`, `EXPECTATION()`, `diff()`, `adl()`, PAC and VAR
#' expectation operators, and correlated shocks (correlations are listed
#' in `notes` and ignored). Other blocks and commands are recorded but not
#' run.
#'
#' @seealso [dsgenl_model()], [solve_dsge()], [bayes_dsge()]
#'
#' @examples
#' rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
#' rbc
#' sol <- solve_dsge(rbc)
#' irf(sol, periods = 20, se = FALSE)
#'
#' # Model code can also be passed as text
#' ar1 <- read_dynare(text = "
#'   var y;
#'   varexo e;
#'   parameters rho;
#'   rho = 0.9;
#'   model;
#'     y = rho * y(-1) + e;
#'   end;
#'   shocks;
#'     var e; stderr 0.01;
#'   end;
#' ")
#' solve_dsge(ar1)
#'
#' @export
read_dynare <- function(file, text = NULL, observed = NULL) {
  if (is.null(text)) {
    if (missing(file) || !is.character(file) || length(file) != 1L) {
      stop("Supply either `file` (a path) or `text`.", call. = FALSE)
    }
    if (!file.exists(file)) {
      stop("File not found: ", file, call. = FALSE)
    }
    src <- readLines(file, warn = FALSE)
    source_name <- file
  } else {
    src <- as.character(text)
    source_name <- "<text>"
  }

  src <- paste(src, collapse = "\n")
  dyn_check_macros(src)
  src <- dyn_strip_comments(src)
  statements <- dyn_split_statements(src)

  parsed <- dyn_parse_statements(statements)
  build <- dyn_build(parsed, observed = observed)
  build$file <- source_name
  structure(build, class = "dsge_dynare")
}

#' @export
print.dsge_dynare <- function(x, ...) {
  cat("Dynare model imported from", x$file, "\n")
  cat("  Endogenous:  ", length(x$variables), "\n")
  cat("  Shocks:      ", length(x$shocks), "\n")
  cat("  Parameters:  ", length(x$parameters),
      "(", length(x$model$free_parameters), "estimated )\n")
  if (nrow(x$aux) > 0L) {
    cat("  Auxiliary:   ", nrow(x$aux), "(",
        paste(x$aux$name, collapse = ", "), ")\n")
  }
  cat("  Observed:    ", paste(x$observed, collapse = ", "), "\n")
  if (!is.null(x$priors)) {
    cat("  Priors:      ", length(x$priors), "\n")
  }
  if (length(x$commands) > 0L) {
    cat("  Commands:    ",
        paste(vapply(x$commands, `[[`, "", "name"), collapse = ", "),
        "(recorded, not run)\n")
  }
  if (length(x$notes) > 0L) {
    cat("\nNotes:\n")
    for (n in x$notes) cat("  -", n, "\n")
  }
  cat("\nUse solve_dsge() on this object to solve at the calibration.\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Lexing
# ---------------------------------------------------------------------------

#' @noRd
dyn_check_macros <- function(src) {
  lines <- strsplit(src, "\n", fixed = TRUE)[[1]]
  macro <- grep("^\\s*@#", lines, value = TRUE)
  if (length(macro) > 0L || grepl("@\\{", src)) {
    stop("Dynare macro-processor directives (@#...) are not supported. ",
         "Expand them first by running `dynare model.mod savemacro onlymacro` ",
         "and import the generated model-macroexp.mod file.", call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
dyn_strip_comments <- function(src) {
  src <- gsub("(?s)/\\*.*?\\*/", " ", src, perl = TRUE)
  src <- gsub("//[^\n]*", "", src, perl = TRUE)
  src <- gsub("%[^\n]*", "", src, perl = TRUE)
  src
}

#' Split source text on ';' outside quotes
#' @noRd
dyn_split_statements <- function(src) {
  if (!grepl("['\"]", src)) {
    parts <- strsplit(src, ";", fixed = TRUE)[[1]]
    if (!grepl(";\\s*$", src) && nzchar(trimws(parts[length(parts)]))) {
      stop("Unterminated statement at end of file (missing ';'): ",
           substr(trimws(parts[length(parts)]), 1L, 60L), call. = FALSE)
    }
    parts <- trimws(gsub("\\s+", " ", parts))
    return(parts[nzchar(parts)])
  }
  chars <- strsplit(src, "", fixed = TRUE)[[1]]
  out <- character(0)
  buf <- character(0)
  quote <- ""
  for (ch in chars) {
    if (nzchar(quote)) {
      if (ch == quote) quote <- ""
      buf <- c(buf, ch)
    } else if (ch == "'" || ch == "\"") {
      quote <- ch
      buf <- c(buf, ch)
    } else if (ch == ";") {
      out <- c(out, paste(buf, collapse = ""))
      buf <- character(0)
    } else {
      buf <- c(buf, ch)
    }
  }
  rest <- trimws(paste(buf, collapse = ""))
  if (nzchar(rest)) {
    stop("Unterminated statement at end of file (missing ';'): ",
         substr(rest, 1L, 60L), call. = FALSE)
  }
  out <- trimws(gsub("\\s+", " ", out))
  out[nzchar(out)]
}

#' Split on commas that are not inside parentheses
#' @noRd
dyn_split_commas <- function(s) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  out <- character(0)
  buf <- character(0)
  depth <- 0L
  for (ch in chars) {
    if (ch == "(") depth <- depth + 1L
    if (ch == ")") depth <- depth - 1L
    if (ch == "," && depth == 0L) {
      out <- c(out, paste(buf, collapse = ""))
      buf <- character(0)
    } else {
      buf <- c(buf, ch)
    }
  }
  trimws(c(out, paste(buf, collapse = "")))
}

# ---------------------------------------------------------------------------
# Statement parsing
# ---------------------------------------------------------------------------

dyn_block_names <- c(
  "model", "initval", "endval", "histval", "steady_state_model", "shocks",
  "mshocks", "estimated_params", "estimated_params_init",
  "estimated_params_bounds", "estimated_params_remove", "observation_trends",
  "optim_weights", "osr_params_bounds", "occbin_constraints", "verbatim",
  "conditional_forecast_paths", "moment_calibration", "irf_calibration",
  "filter_initial_state", "homotopy_setup", "ramsey_constraints",
  "shock_groups", "init2shocks", "heteroskedastic_shocks",
  "matched_moments", "model_replace", "model_options", "perfect_foresight_controlled_paths",
  "svar_identification", "markov_switching", "var_model", "pac_model",
  "generate_irfs", "epilogue", "declare_optimal_policy_discretionary"
)

#' @noRd
dyn_parse_statements <- function(statements) {
  p <- list(
    endo = character(0), exo = character(0), params = character(0),
    predetermined = character(0), varobs = NULL,
    param_exprs = list(), model = character(0), model_linear = FALSE,
    blocks = list(), commands = list(), notes = character(0)
  )

  i <- 1L
  n <- length(statements)
  while (i <= n) {
    st <- statements[i]
    kw <- tolower(sub("^([A-Za-z_][A-Za-z0-9_]*).*$", "\\1", st))
    head_only <- grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*(\\(.*\\))?$", st)

    if (kw %in% dyn_block_names && head_only) {
      opts <- if (grepl("\\(", st)) sub("^[^(]*\\((.*)\\)$", "\\1", st) else ""
      j <- i + 1L
      body <- character(0)
      while (j <= n && tolower(statements[j]) != "end") {
        body <- c(body, statements[j])
        j <- j + 1L
      }
      if (j > n) {
        stop("Block '", kw, "' is not closed with 'end;'.", call. = FALSE)
      }
      if (kw == "model") {
        if (length(p$model) > 0L) {
          stop("Multiple model blocks are not supported.", call. = FALSE)
        }
        p$model <- body
        p$model_linear <- grepl("\\blinear\\b", opts)
      } else {
        p$blocks[[kw]] <- c(p$blocks[[kw]], body)
        if (kw %in% c("steady_state_model", "shocks", "initval")) {
          p$blocks[[paste0(kw, "_opts")]] <- opts
        }
      }
      i <- j + 1L
      next
    }

    if (kw %in% c("var", "varexo", "parameters", "predetermined_variables",
                  "varobs", "model_local_variable") &&
        !grepl("^[A-Za-z_]+\\s*=", st)) {
      names_found <- dyn_declared_names(st, kw)
      if (kw == "var") {
        if (grepl("^var\\s*\\([^)]*\\blog\\b", st)) {
          stop("var(log) declarations are not supported.", call. = FALSE)
        }
        p$endo <- c(p$endo, names_found)
      } else if (kw == "varexo") {
        p$exo <- c(p$exo, names_found)
      } else if (kw == "parameters") {
        p$params <- c(p$params, names_found)
      } else if (kw == "predetermined_variables") {
        p$predetermined <- c(p$predetermined, names_found)
      } else if (kw == "varobs") {
        p$varobs <- c(p$varobs, names_found)
      }
      i <- i + 1L
      next
    }

    if (kw %in% c("varexo_det", "trend_var", "log_trend_var",
                  "external_function", "change_type")) {
      stop("'", kw, "' is not supported by read_dynare().", call. = FALSE)
    }

    # Parameter assignment: name = expr
    if (grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*=[^=]", st)) {
      lhs <- trimws(sub("=.*$", "", st))
      rhs <- trimws(sub("^[^=]*=", "", st))
      # Assignments to undeclared names are plain MATLAB variables in
      # Dynare; keep them as constants usable in later expressions.
      p$param_exprs[[length(p$param_exprs) + 1L]] <- list(name = lhs,
                                                          expr = rhs)
      i <- i + 1L
      next
    }

    # MATLAB-style statements (options_.x = ..., set_param_value, ...)
    if (grepl("^[A-Za-z_][A-Za-z0-9_]*\\.", st) ||
        grepl("^set_param_value", st)) {
      p$notes <- c(p$notes, paste0("Ignored MATLAB statement: ",
                                   substr(st, 1L, 60L)))
      i <- i + 1L
      next
    }

    # Anything else is a command, e.g. stoch_simul(order=1) y c;
    opts <- ""
    if (grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*\\(", st)) {
      opts <- dyn_paren_content(st, regexpr("\\(", st))
    }
    p$commands[[length(p$commands) + 1L]] <- list(name = kw, options = opts,
                                                  statement = st)
    i <- i + 1L
  }

  if (length(p$endo) == 0L) {
    stop("No endogenous variables declared ('var').", call. = FALSE)
  }
  if (length(p$model) == 0L) {
    stop("No model block found.", call. = FALSE)
  }
  all_names <- c(p$endo, p$exo, p$params)
  reserved <- c("if", "else", "repeat", "while", "function", "for", "next",
                "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "in",
                "T", "F")
  bad <- intersect(all_names, reserved)
  if (length(bad) > 0L) {
    stop("Name(s) reserved in R cannot be imported; rename them in the ",
         ".mod file: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(all_names)) {
    stop("Name declared more than once: ",
         paste(unique(all_names[duplicated(all_names)]), collapse = ", "),
         call. = FALSE)
  }
  p
}

#' Names from a declaration statement, dropping LaTeX and options
#' @noRd
dyn_declared_names <- function(st, kw) {
  body <- sub(paste0("^", kw, "\\s*(\\([^)]*\\))?"), "", st, ignore.case = TRUE)
  body <- gsub("\\$[^$]*\\$", " ", body)
  body <- gsub("\\((?:[^()'\"]|'[^']*'|\"[^\"]*\")*\\)", " ", body, perl = TRUE)
  toks <- strsplit(trimws(body), "[[:space:],]+")[[1]]
  toks <- toks[nzchar(toks)]
  bad <- toks[!grepl("^[A-Za-z_][A-Za-z0-9_]*$", toks)]
  if (length(bad) > 0L) {
    stop("Cannot parse '", kw, "' declaration near: ",
         paste(bad, collapse = " "), call. = FALSE)
  }
  toks
}

#' Content of the parenthesised group starting at position `start`
#' @noRd
dyn_paren_content <- function(s, start) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  depth <- 0L
  for (k in seq(start, length(chars))) {
    if (chars[k] == "(") depth <- depth + 1L
    if (chars[k] == ")") {
      depth <- depth - 1L
      if (depth == 0L) {
        return(substr(s, start + 1L, k - 1L))
      }
    }
  }
  stop("Unbalanced parentheses in: ", s, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Expression rewriting
# ---------------------------------------------------------------------------

#' Rewrite identifiers (optionally followed by an integer time index)
#'
#' `fn(name, lag, has_index)` returns the replacement text, or NULL to
#' leave the token unchanged.
#' @noRd
dyn_rewrite_ids <- function(s, fn) {
  pat <- paste0("(?<![A-Za-z0-9_.])([A-Za-z_][A-Za-z0-9_]*)",
                "(\\s*\\(\\s*[+-]?\\s*[0-9]+\\s*\\))?")
  m <- gregexpr(pat, s, perl = TRUE)
  toks <- regmatches(s, m)[[1]]
  if (length(toks) == 0L) return(s)
  new <- vapply(toks, function(tok) {
    name <- sub("^([A-Za-z_][A-Za-z0-9_]*).*$", "\\1", tok)
    has_index <- grepl("(", tok, fixed = TRUE)
    lag <- if (has_index) {
      as.integer(gsub("[^0-9+-]", "", sub("^[^(]*\\(", "", tok)))
    } else {
      0L
    }
    r <- fn(name, lag, has_index)
    if (is.null(r)) tok else r
  }, character(1), USE.NAMES = FALSE)
  regmatches(s, m) <- list(new)
  s
}

#' Replace calls f(ARG) by a template in which ARG is substituted
#' @noRd
dyn_rewrite_calls <- function(s, fname, template) {
  pat <- paste0("(?<![A-Za-z0-9_.])", fname, "\\s*\\(")
  repeat {
    pos <- regexpr(pat, s, perl = TRUE)
    if (pos == -1L) return(s)
    open <- pos + attr(pos, "match.length") - 1L
    arg <- dyn_paren_content(s, open)
    end <- open + nchar(arg) + 1L
    repl <- gsub("ARG", arg, template, fixed = TRUE)
    s <- paste0(substr(s, 1L, pos - 1L), repl, substr(s, end + 1L, nchar(s)))
  }
}

#' Translate Dynare functions into R
#' @noRd
dyn_translate_math <- function(s) {
  unsupported <- c("STEADY_STATE", "EXPECTATION", "diff", "adl",
                   "var_expectation", "pac_expectation",
                   "pac_target_nonstationary")
  for (f in unsupported) {
    if (grepl(paste0("(?<![A-Za-z0-9_.])", f, "\\s*\\("), s, perl = TRUE)) {
      stop("The Dynare operator ", f, "() is not supported by read_dynare().",
           call. = FALSE)
    }
  }
  s <- gsub(".^", "^", s, fixed = TRUE)
  s <- gsub(".*", "*", s, fixed = TRUE)
  s <- gsub("./", "/", s, fixed = TRUE)
  s <- gsub("(?<![A-Za-z0-9_.])ln\\s*\\(", "log(", s, perl = TRUE)
  s <- gsub("(?<![A-Za-z0-9_.])normcdf\\s*\\(", "stats::pnorm(", s, perl = TRUE)
  s <- gsub("(?<![A-Za-z0-9_.])normpdf\\s*\\(", "stats::dnorm(", s, perl = TRUE)
  s <- dyn_rewrite_calls(s, "erfc", "(2 * stats::pnorm(-sqrt(2) * (ARG)))")
  s <- dyn_rewrite_calls(s, "erf", "(2 * stats::pnorm(sqrt(2) * (ARG)) - 1)")
  s <- gsub("(?<![A-Za-z0-9_.])(inf|Inf)(?![A-Za-z0-9_.])", "Inf", s,
            perl = TRUE)
  s
}

#' Evaluate a Dynare expression in an environment
#' @noRd
dyn_eval <- function(expr, env, what) {
  if (!nzchar(trimws(expr))) return(NA_real_)
  val <- tryCatch(
    eval(parse(text = dyn_translate_math(expr)), envir = env),
    error = function(e) {
      stop("Cannot evaluate ", what, ": ", expr, "\n  ", conditionMessage(e),
           call. = FALSE)
    }
  )
  as.numeric(val)
}

# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

#' @noRd
dyn_build <- function(p, observed = NULL) {
  notes <- p$notes
  endo <- p$endo
  exo <- p$exo

  # --- calibration --------------------------------------------------------
  cal_env <- new.env(parent = baseenv())
  for (pe in p$param_exprs) {
    assign(pe$name, dyn_eval(pe$expr, cal_env,
                             paste0("parameter '", pe$name, "'")),
           envir = cal_env)
  }
  params <- unlist(mget(intersect(p$params, ls(cal_env)), envir = cal_env))
  if (is.null(params)) params <- numeric(0)
  params <- params[intersect(p$params, names(params))]

  # --- model block --------------------------------------------------------
  eqs <- dyn_model_equations(p$model, p$predetermined, endo, exo)

  # --- shocks (correlated shocks are orthogonalised in the equations) -----
  shk <- dyn_parse_shocks(p$blocks$shocks, exo, cal_env)
  notes <- c(notes, shk$notes)
  if (length(shk$cross) > 0L) {
    orth <- dyn_orthogonalise(shk$sd, shk$cross, exo)
    eqs <- vapply(eqs, function(eq) {
      dyn_rewrite_ids(eq, function(name, lag, has_index) {
        if (!name %in% names(orth$subst)) return(NULL)
        parts <- vapply(orth$subst[[name]], function(t) {
          dyn_rewrite_ids(t, function(n2, l2, h2) {
            if (n2 %in% exo) dyn_timed(n2, lag) else NULL
          })
        }, character(1))
        paste0("(", paste(parts, collapse = " + "), ")")
      })
    }, character(1), USE.NAMES = FALSE)
    shk$sd <- orth$sd
    notes <- c(notes, paste0(
      "Correlated shocks orthogonalised by Cholesky factorisation in the ",
      "order ", paste(exo, collapse = ", "), " (as in Dynare's IRFs); shock ",
      "standard deviations refer to the orthogonal shocks."))
  }

  # Leads and lags present in the model
  max_lead <- stats::setNames(integer(length(endo)), endo)
  max_lag <- stats::setNames(integer(length(c(endo, exo))), c(endo, exo))
  for (eq in eqs) {
    dyn_rewrite_ids(eq, function(name, lag, has_index) {
      if (name %in% endo) {
        if (lag > max_lead[name]) max_lead[name] <<- lag
        if (-lag > max_lag[name]) max_lag[name] <<- -lag
      } else if (name %in% exo) {
        if (lag > 0L) {
          stop("Lead of shock '", name, "' is not supported.", call. = FALSE)
        }
        if (-lag > max_lag[name]) max_lag[name] <<- -lag
      }
      NULL
    })
  }

  lag_name <- function(v, j) paste0(v, "_lag", j)
  lead_name <- function(v, j) paste0(v, "_lead", j)

  rewrite_timing <- function(name, lag, has_index) {
    if (name %in% endo) {
      if (lag == 0L) return(name)
      if (lag == 1L) return(paste0(name, "(+1)"))
      if (lag > 1L) return(paste0(lead_name(name, lag - 1L), "(+1)"))
      return(lag_name(name, -lag))
    }
    if (name %in% exo) {
      if (lag == 0L) return(name)
      return(lag_name(name, -lag))
    }
    if (has_index && !exists(name, envir = baseenv(), mode = "function")) {
      stop("Time index on '", name, "', which is not a declared variable.",
           call. = FALSE)
    }
    NULL
  }
  model_eqs <- vapply(eqs, function(eq) {
    dyn_rewrite_ids(eq, rewrite_timing)
  }, character(1), USE.NAMES = FALSE)

  # Auxiliary lead controls and lag states
  aux <- data.frame(name = character(0), base = character(0),
                    type = character(0), shift = integer(0),
                    stringsAsFactors = FALSE)
  aux_ctrl_eqs <- character(0)
  lead_ctrls <- character(0)
  for (v in endo) {
    if (max_lead[v] >= 2L) {
      for (j in seq_len(max_lead[v] - 1L)) {
        nm <- lead_name(v, j)
        prev <- if (j == 1L) v else lead_name(v, j - 1L)
        aux_ctrl_eqs <- c(aux_ctrl_eqs, paste0(nm, " = ", prev, "(+1)"))
        lead_ctrls <- c(lead_ctrls, nm)
        aux <- rbind(aux, data.frame(name = nm, base = v, type = "lead",
                                     shift = j, stringsAsFactors = FALSE))
      }
    }
  }
  aux_state_eqs <- character(0)
  lag_states <- character(0)
  for (v in c(endo, exo)) {
    if (max_lag[v] >= 1L) {
      for (j in seq_len(max_lag[v])) {
        nm <- lag_name(v, j)
        prev <- if (j == 1L) v else lag_name(v, j - 1L)
        aux_state_eqs <- c(aux_state_eqs, paste0(nm, "(+1) = ", prev))
        lag_states <- c(lag_states, nm)
        aux <- rbind(aux, data.frame(name = nm, base = v, type = "lag",
                                     shift = -j, stringsAsFactors = FALSE))
      }
    }
  }
  clash <- intersect(aux$name, c(endo, exo, p$params))
  if (length(clash) > 0L) {
    stop("Auxiliary variable name(s) clash with declared names: ",
         paste(clash, collapse = ", "), call. = FALSE)
  }

  # --- steady state ------------------------------------------------------
  exo_ss <- stats::setNames(numeric(length(exo)), exo)
  init_vals <- list()
  if (!is.null(p$blocks$initval)) {
    init_env <- new.env(parent = cal_env)
    for (st in p$blocks$initval) {
      if (!grepl("=", st)) next
      nm <- trimws(sub("=.*$", "", st))
      val <- dyn_eval(sub("^[^=]*=", "", st), init_env,
                      paste0("initval for '", nm, "'"))
      assign(nm, val, envir = init_env)
      init_vals[[nm]] <- val
    }
    for (e in intersect(names(init_vals), exo)) exo_ss[e] <- init_vals[[e]]
  }

  shock_state_eqs <- paste0(exo, "(+1) = ", format_num(exo_ss))

  default_guess <- if (p$model_linear) 0 else 1
  base_guess <- stats::setNames(rep(default_guess, length(endo)), endo)
  for (v in intersect(names(init_vals), endo)) base_guess[v] <- init_vals[[v]]
  fill_aux <- function(base_vals) {
    out <- c(base_vals, exo_ss)
    for (k in seq_len(nrow(aux))) out[aux$name[k]] <- out[[aux$base[k]]]
    out
  }
  ss_guess <- fill_aux(base_guess)

  ss_function <- NULL
  if (!is.null(p$blocks$steady_state_model)) {
    ss_function <- dyn_ss_function(p$blocks$steady_state_model, cal_env,
                                   endo, exo_ss, fill_aux)
  }

  # --- observed variables ------------------------------------------------
  obs_source <- if (!is.null(observed)) "observed" else "varobs"
  obs <- if (!is.null(observed)) observed else p$varobs
  unknown_obs <- setdiff(obs, endo)
  if (length(unknown_obs) > 0L) {
    stop("Observed variable(s) not declared with 'var': ",
         paste(unknown_obs, collapse = ", "), call. = FALSE)
  }
  n_exo <- length(exo)
  if (n_exo > length(endo)) {
    stop("dsge needs at least as many endogenous variables as shocks.",
         call. = FALSE)
  }
  if (length(obs) != n_exo) {
    chosen <- c(obs, setdiff(endo, obs))[seq_len(n_exo)]
    if (is.null(obs)) {
      notes <- c(notes, paste0(
        "No varobs: '", paste(chosen, collapse = "', '"),
        "' marked as observed (dsge needs one observed variable per shock; ",
        "this only matters for estimation)."))
    } else {
      notes <- c(notes, paste0(
        obs_source, " lists ", length(obs), " variable(s) but the model has ",
        n_exo, " shock(s); dsge needs equal numbers, so observed was set to '",
        paste(chosen, collapse = "', '"), "'."))
    }
    obs <- chosen
  }

  # --- estimated parameters and priors -----------------------------------
  est <- dyn_parse_estimated(p$blocks$estimated_params,
                             p$blocks$estimated_params_init,
                             p$params, exo, cal_env)
  notes <- c(notes, est$notes)

  # --- assemble dsgenl_model ---------------------------------------------
  controls <- c(endo, lead_ctrls)
  unobs <- setdiff(controls, obs)
  all_eqs <- c(model_eqs, aux_ctrl_eqs, shock_state_eqs, aux_state_eqs)
  make_model <- function(fixed, start) {
    do.call(dsgenl_model, c(
      as.list(all_eqs),
      list(observed = obs, unobserved = unobs, exo_state = exo,
           endo_state = lag_states, fixed = fixed, start = start,
           ss_guess = ss_guess, ss_function = ss_function)
    ))
  }
  probe <- make_model(list(), list())
  model_params <- probe$parameters

  undeclared <- setdiff(model_params, p$params)
  if (length(undeclared) > 0L) {
    stop("Undeclared name(s) in the model block: ",
         paste(undeclared, collapse = ", "), call. = FALSE)
  }

  est_all <- est$table$name[est$table$type == "param"]
  est_names <- intersect(est_all, model_params)
  unused <- setdiff(est_all, model_params)
  if (length(unused) > 0L) {
    notes <- c(notes, paste0(
      "Estimated parameter(s) not used in the model equations were dropped: ",
      paste(unused, collapse = ", "), "."))
    est$priors <- est$priors[setdiff(names(est$priors), unused)]
    if (length(est$priors) == 0L) est$priors <- NULL
  }
  start <- list()
  for (nm in est_names) {
    init <- est$table$init[est$table$name == nm & est$table$type == "param"]
    if (is.na(init)) init <- params[nm]
    if (is.na(init)) init <- est$table$prior_mean[est$table$name == nm &
                                                    est$table$type == "param"]
    if (is.na(init)) {
      stop("No calibrated or initial value for estimated parameter '", nm,
           "'.", call. = FALSE)
    }
    start[[nm]] <- unname(init)
  }
  fixed_names <- setdiff(model_params, est_names)
  missing_cal <- fixed_names[!fixed_names %in% names(params)]
  if (length(missing_cal) > 0L) {
    stop("No value assigned to parameter(s): ",
         paste(missing_cal, collapse = ", "), call. = FALSE)
  }
  fixed <- as.list(params[fixed_names])

  model <- make_model(fixed, start)

  # Calibrated values of estimated parameters complete the calibration
  for (nm in est_names) {
    if (!nm %in% names(params)) params[nm] <- start[[nm]]
  }

  shock_sd <- shk$sd
  for (k in seq_len(nrow(est$table))) {
    if (est$table$type[k] == "stderr" && is.na(shk$declared[est$table$name[k]])) {
      v <- est$table$init[k]
      if (is.na(v)) v <- est$table$prior_mean[k]
      if (!is.na(v) && is.finite(v)) shock_sd[est$table$name[k]] <- v
    }
  }
  zero_sd <- names(shock_sd)[shock_sd == 0]
  if (length(zero_sd) > 0L && length(zero_sd) < length(shock_sd)) {
    notes <- c(notes, paste0("No variance declared for shock(s) ",
                             paste(zero_sd, collapse = ", "),
                             "; standard deviation set to 0 as in Dynare."))
  } else if (length(zero_sd) > 0L) {
    notes <- c(notes, paste0("No shock variances declared; all standard ",
                             "deviations are 0 as in Dynare. Pass shock_sd ",
                             "to solve_dsge() to simulate."))
  }

  for (cmd in p$commands) {
    if (cmd$name %in% c("ramsey_model", "ramsey_policy", "planner_objective",
                        "discretionary_policy", "osr", "occbin_setup")) {
      notes <- c(notes, paste0("Command '", cmd$name, "' recorded but not ",
                               "translated; see the matching dsge function."))
    }
  }
  ignored_blocks <- setdiff(names(p$blocks),
                            c("initval", "initval_opts", "steady_state_model",
                              "steady_state_model_opts", "shocks",
                              "shocks_opts", "estimated_params",
                              "estimated_params_init"))
  if (length(ignored_blocks) > 0L) {
    notes <- c(notes, paste0("Block(s) not translated: ",
                             paste(ignored_blocks, collapse = ", "), "."))
  }

  list(
    model = model,
    params = params,
    shock_sd = shock_sd,
    priors = est$priors,
    estimated_params = est$table,
    observed = obs,
    variables = endo,
    shocks = exo,
    parameters = p$params,
    aux = aux,
    commands = p$commands,
    notes = unique(notes)
  )
}

#' Model block statements -> equation strings in Dynare timing
#' @noRd
dyn_model_equations <- function(statements, predetermined, endo, exo) {
  locals <- list()
  eqs <- character(0)
  inline_locals <- function(s) {
    if (length(locals) == 0L) return(s)
    dyn_rewrite_ids(s, function(name, lag, has_index) {
      if (!name %in% names(locals)) return(NULL)
      if (has_index) {
        stop("Model-local variable '", name, "' used with a time index.",
             call. = FALSE)
      }
      paste0("(", locals[[name]], ")")
    })
  }
  for (st in statements) {
    st <- sub("^\\[[^]]*\\]\\s*", "", st)
    if (!nzchar(st)) next
    if (startsWith(st, "#")) {
      def <- sub("^#\\s*", "", st)
      nm <- trimws(sub("=.*$", "", def))
      if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", nm) || !grepl("=", def)) {
        stop("Cannot parse model-local variable: ", st, call. = FALSE)
      }
      locals[[nm]] <- inline_locals(trimws(sub("^[^=]*=", "", def)))
      next
    }
    eqs <- c(eqs, inline_locals(st))
  }

  if (length(predetermined) > 0L) {
    unknown <- setdiff(predetermined, endo)
    if (length(unknown) > 0L) {
      stop("predetermined_variables not declared with 'var': ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    # Dynare convention: x is beginning-of-period, x(+1) chosen today.
    eqs <- vapply(eqs, function(eq) {
      dyn_rewrite_ids(eq, function(name, lag, has_index) {
        if (!name %in% predetermined) return(NULL)
        dyn_timed(name, lag - 1L)
      })
    }, character(1), USE.NAMES = FALSE)
  }

  vapply(eqs, function(eq) {
    eq <- dyn_translate_math(eq)
    n_eq <- lengths(regmatches(eq, gregexpr("=", eq, fixed = TRUE)))
    if (n_eq == 0L) {
      eq <- paste0(eq, " = 0")
    } else if (n_eq > 1L) {
      stop("Equation has more than one '=': ", eq, call. = FALSE)
    }
    eq
  }, character(1), USE.NAMES = FALSE)
}

#' Dynare-timed reference text
#' @noRd
dyn_timed <- function(name, lag) {
  if (lag == 0L) name else sprintf("%s(%+d)", name, lag)
}

#' @noRd
format_num <- function(x) {
  sprintf("%.17g", x)
}

#' Build a steady-state function from a steady_state_model block
#' @noRd
dyn_ss_function <- function(statements, cal_env, endo, exo_ss, fill_aux) {
  assigns <- lapply(statements, function(st) {
    if (!grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*=[^=]", st)) {
      stop("Unsupported statement in steady_state_model: ", st,
           call. = FALSE)
    }
    list(name = trimws(sub("=.*$", "", st)),
         expr = parse(text = dyn_translate_math(sub("^[^=]*=", "", st))))
  })
  defined <- vapply(assigns, `[[`, "", "name")
  missing_v <- setdiff(endo, defined)
  if (length(missing_v) > 0L) {
    stop("steady_state_model does not define: ",
         paste(missing_v, collapse = ", "), call. = FALSE)
  }
  force(cal_env)
  force(exo_ss)
  force(fill_aux)
  function(params) {
    env <- new.env(parent = baseenv())
    for (nm in ls(cal_env)) assign(nm, get(nm, envir = cal_env), envir = env)
    for (nm in names(exo_ss)) assign(nm, exo_ss[[nm]], envir = env)
    for (nm in names(params)) assign(nm, params[[nm]], envir = env)
    for (a in assigns) {
      assign(a$name, eval(a$expr, envir = env), envir = env)
    }
    fill_aux(unlist(mget(endo, envir = env)))
  }
}

#' Parse the shocks block (variances only)
#' @noRd
dyn_parse_shocks <- function(statements, exo, cal_env) {
  sd <- stats::setNames(numeric(length(exo)), exo)
  declared <- stats::setNames(rep(NA_real_, length(exo)), exo)
  notes <- character(0)
  cross <- list()
  current <- NULL
  for (st in statements) {
    if (grepl("^var\\s", st)) {
      body <- sub("^var\\s+", "", st)
      if (grepl("=", body)) {
        lhs <- strsplit(trimws(sub("=.*$", "", body)), "[[:space:],]+")[[1]]
        if (length(lhs) == 1L) {
          dyn_check_shock(lhs, exo)
          val <- dyn_eval(sub("^[^=]*=", "", body), cal_env,
                          "shock variance")
          sd[lhs] <- sqrt(val)
          declared[lhs] <- sd[lhs]
        } else if (length(lhs) == 2L) {
          for (nm in lhs) dyn_check_shock(nm, exo)
          val <- dyn_eval(sub("^[^=]*=", "", body), cal_env,
                          "shock covariance")
          cross[[length(cross) + 1L]] <- list(a = lhs[1], b = lhs[2],
                                              type = "cov", value = val)
        } else {
          stop("Cannot parse shocks statement: ", st, call. = FALSE)
        }
        current <- NULL
      } else {
        current <- trimws(body)
        dyn_check_shock(current, exo)
      }
    } else if (grepl("^stderr\\s", st)) {
      if (is.null(current)) {
        stop("'stderr' without a preceding 'var' in the shocks block.",
             call. = FALSE)
      }
      sd[current] <- dyn_eval(sub("^stderr\\s+", "", st), cal_env,
                              "shock standard deviation")
      declared[current] <- sd[current]
    } else if (grepl("^corr\\s", st)) {
      lhs <- strsplit(trimws(sub("=.*$", "", sub("^corr\\s+", "", st))),
                      "[[:space:],]+")[[1]]
      if (length(lhs) != 2L || !grepl("=", st)) {
        stop("Cannot parse shocks statement: ", st, call. = FALSE)
      }
      for (nm in lhs) dyn_check_shock(nm, exo)
      cross[[length(cross) + 1L]] <- list(
        a = lhs[1], b = lhs[2], type = "corr",
        value = dyn_eval(sub("^[^=]*=", "", st), cal_env, "shock correlation"))
    } else if (grepl("^(periods|values)\\s", st)) {
      if (!is.null(current)) {
        notes <- c(notes, paste0("Deterministic shock path for '", current,
                                 "' ignored; use perfect_foresight()."))
      }
    } else {
      notes <- c(notes, paste0("Unrecognised shocks statement ignored: ",
                               st, "."))
    }
  }
  list(sd = sd, declared = declared, cross = cross, notes = unique(notes))
}

#' Cholesky rewrite for correlated shocks
#'
#' With Sigma = L L', shock i is replaced in the equations by
#' e_i + sum_{j<i} (L_ij / L_jj) e_j, where the e_j are orthogonal with
#' standard deviations L_jj. This reproduces the covariance and Dynare's
#' Cholesky-ordered impulse responses (ordering = varexo declaration order).
#' @noRd
dyn_orthogonalise <- function(sd, cross, exo) {
  Sigma <- diag(sd^2, nrow = length(exo))
  dimnames(Sigma) <- list(exo, exo)
  for (cr in cross) {
    v <- if (cr$type == "cov") cr$value else cr$value * sd[cr$a] * sd[cr$b]
    Sigma[cr$a, cr$b] <- Sigma[cr$b, cr$a] <- v
  }
  L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
  if (is.null(L)) {
    stop("Shock covariance matrix is not positive definite.", call. = FALSE)
  }
  subst <- list()
  for (i in seq_along(exo)) {
    terms <- exo[i]
    for (j in seq_len(i - 1L)) {
      if (abs(L[i, j]) > 0) {
        terms <- c(terms, paste0(sprintf("%.17g", L[i, j] / L[j, j]), " * ",
                                 exo[j]))
      }
    }
    if (length(terms) > 1L) subst[[exo[i]]] <- terms
  }
  list(sd = stats::setNames(diag(L), exo), subst = subst)
}

#' @noRd
dyn_check_shock <- function(name, exo) {
  if (!name %in% exo) {
    stop("'", name, "' in the shocks block is not declared with 'varexo'.",
         call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# estimated_params
# ---------------------------------------------------------------------------

dyn_prior_shapes <- c(
  beta_pdf = "beta", gamma_pdf = "gamma", normal_pdf = "normal",
  uniform_pdf = "uniform", inv_gamma_pdf = "inv_gamma1",
  inv_gamma1_pdf = "inv_gamma1", inv_gamma2_pdf = "inv_gamma2",
  weibull_pdf = "weibull"
)

#' @noRd
dyn_parse_estimated <- function(statements, init_statements, params, exo,
                                cal_env) {
  table <- data.frame(
    name = character(0), type = character(0), init = numeric(0),
    lower = numeric(0), upper = numeric(0), shape = character(0),
    prior_mean = numeric(0), prior_sd = numeric(0), translation = character(0),
    stringsAsFactors = FALSE
  )
  priors <- list()
  notes <- character(0)
  num <- function(x, what) {
    if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_real_)
    dyn_eval(x, cal_env, what)
  }

  target <- function(fields, st) {
    f1 <- fields[1]
    if (grepl("^stderr\\s", f1)) {
      nm <- trimws(sub("^stderr\\s+", "", f1))
      if (!nm %in% exo) {
        stop("stderr of undeclared shock '", nm, "' in: ", st, call. = FALSE)
      }
      return(list(name = nm, type = "stderr", rest = fields[-1]))
    }
    if (grepl("^corr\\s", f1)) {
      return(list(name = paste(trimws(sub("^corr\\s+", "", f1)), fields[2]),
                  type = "corr", rest = fields[-(1:2)]))
    }
    if (!f1 %in% params) {
      stop("Estimated parameter '", f1, "' is not declared with 'parameters'.",
           call. = FALSE)
    }
    list(name = f1, type = "param", rest = fields[-1])
  }

  for (st in statements) {
    fields <- dyn_split_commas(st)
    tg <- target(fields, st)
    rest <- tg$rest
    shape_pos <- which(rest %in% names(dyn_prior_shapes))
    init <- lower <- upper <- mean <- sd <- NA_real_
    shape <- NA_character_
    p3 <- p4 <- NA_real_
    what <- paste0("estimated_params entry '", tg$name, "'")
    if (length(shape_pos) == 0L) {
      init <- num(rest[1], what)
      lower <- num(rest[2], what)
      upper <- num(rest[3], what)
    } else {
      sp <- shape_pos[1]
      pre <- rest[seq_len(sp - 1L)]
      if (length(pre) >= 1L) init <- num(pre[1], what)
      if (length(pre) >= 3L) {
        lower <- num(pre[2], what)
        upper <- num(pre[3], what)
      }
      shape <- dyn_prior_shapes[[rest[sp]]]
      post <- rest[-seq_len(sp)]
      mean <- num(post[1], what)
      sd <- num(post[2], what)
      p3 <- num(post[3], what)
      p4 <- num(post[4], what)
    }

    translation <- if (is.na(shape)) "no prior (ML entry)" else ""
    if (tg$type == "corr") {
      translation <- "ignored (shock correlations not supported)"
      notes <- c(notes, paste0("Estimated correlation '", tg$name,
                               "' ignored."))
    } else if (!is.na(shape)) {
      pr <- dyn_translate_prior(shape, mean, sd, p3, p4)
      translation <- pr$translation
      if (!is.null(pr$prior)) {
        key <- if (tg$type == "stderr") paste0("sd_e.", tg$name) else tg$name
        priors[[key]] <- pr$prior
      } else {
        notes <- c(notes, paste0("Prior for '", tg$name, "' not translated: ",
                                 pr$translation, "."))
      }
      if (!is.null(pr$note)) {
        notes <- c(notes, paste0("Prior for '", tg$name, "': ", pr$note, "."))
      }
    }

    table <- rbind(table, data.frame(
      name = tg$name, type = tg$type, init = init, lower = lower,
      upper = upper, shape = if (is.na(shape)) NA_character_ else shape,
      prior_mean = mean, prior_sd = sd, translation = translation,
      stringsAsFactors = FALSE
    ))
  }

  # estimated_params_init overrides initial values
  for (st in init_statements) {
    fields <- dyn_split_commas(st)
    tg <- target(fields, st)
    idx <- which(table$name == tg$name & table$type == tg$type)
    if (length(idx) == 1L) {
      table$init[idx] <- num(tg$rest[1], "estimated_params_init entry")
    }
  }

  n_param <- sum(table$type == "param")
  n_with_prior <- sum(table$type == "param" & !is.na(table$shape))
  if (n_with_prior > 0L && n_with_prior < n_param) {
    notes <- c(notes, paste0("Some estimated parameters have no prior; ",
                             "add them before calling bayes_dsge()."))
  }
  if (length(priors) == 0L) priors <- NULL

  list(table = table, priors = priors, notes = notes)
}

#' Convert a Dynare prior (mean, sd, p3, p4) into a dsge prior
#' @noRd
dyn_translate_prior <- function(shape, mean, sd, p3, p4) {
  out <- list(prior = NULL, translation = "", note = NULL)
  shifted <- !is.na(p3) && p3 != 0
  if (shape == "weibull") {
    out$translation <- "weibull_pdf has no dsge counterpart"
    return(out)
  }
  if (shape == "uniform") {
    lo <- if (!is.na(p3)) p3 else mean - sqrt(3) * sd
    hi <- if (!is.na(p4)) p4 else mean + sqrt(3) * sd
    if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
      out$translation <- "invalid uniform bounds"
      return(out)
    }
    out$prior <- prior("uniform", min = lo, max = hi)
    out$translation <- sprintf("uniform(min = %g, max = %g)", lo, hi)
    return(out)
  }
  if (is.na(mean) || is.na(sd)) {
    out$translation <- "missing prior mean or standard deviation"
    return(out)
  }
  if (shape == "normal") {
    out$prior <- prior("normal", mean = mean, sd = sd)
    out$translation <- sprintf("normal(mean = %g, sd = %g)", mean, sd)
    return(out)
  }
  if (shape == "beta") {
    if ((!is.na(p3) && p3 != 0) || (!is.na(p4) && p4 != 1)) {
      out$translation <- "generalised beta on a non-unit interval"
      return(out)
    }
    k <- mean * (1 - mean) / sd^2 - 1
    if (mean <= 0 || mean >= 1 || k <= 0) {
      out$translation <- "beta mean/sd are not feasible"
      return(out)
    }
    out$prior <- prior("beta", shape1 = mean * k, shape2 = (1 - mean) * k)
    out$translation <- sprintf("beta(shape1 = %g, shape2 = %g)", mean * k,
                               (1 - mean) * k)
    return(out)
  }
  if (shifted) {
    out$translation <- "shifted prior (non-zero lower bound)"
    return(out)
  }
  if (mean <= 0) {
    out$translation <- "prior mean must be positive"
    return(out)
  }
  if (shape == "gamma") {
    out$prior <- prior("gamma", shape = mean^2 / sd^2, rate = mean / sd^2)
    out$translation <- sprintf("gamma(shape = %g, rate = %g)", mean^2 / sd^2,
                               mean / sd^2)
    return(out)
  }
  # Inverse gamma: dsge density x^-(a+1) exp(-b/x); moment matching.
  a <- if (is.finite(sd)) 2 + mean^2 / sd^2 else 2
  b <- mean * (a - 1)
  out$prior <- prior("inv_gamma", shape = a, scale = b)
  out$translation <- sprintf("inv_gamma(shape = %g, scale = %g)", a, b)
  if (shape == "inv_gamma1") {
    out$note <- paste0("inv_gamma_pdf approximated by an inverse gamma with ",
                       "the same mean and sd")
  }
  out
}

#' Unwrap a dsge_dynare object passed to a model-taking function
#' @noRd
unwrap_dynare <- function(x) {
  if (inherits(x, "dsge_dynare")) x$model else x
}
