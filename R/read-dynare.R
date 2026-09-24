# Import Dynare .mod files
#
# Translates a Dynare .mod file into a dsgenl_model plus the calibrated
# parameter values, shock standard deviations and priors needed to solve
# or estimate it. Macro directives are expanded first (dynare-macro.R);
# optimal-policy and OccBin blocks are handled in dynare-policy.R and
# dynare-occbin.R.
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
#' model, together with the calibration, shock standard deviations,
#' measurement errors, priors, optimal-policy problem and occasionally
#' binding constraints it declares, so that the model can be solved,
#' simulated or estimated in R without Dynare, MATLAB or Octave.
#'
#' @param file Path to a `.mod` file.
#' @param text Alternatively, the model code as a character vector (one
#'   element per line, or a single string). Used instead of `file` when
#'   supplied.
#' @param observed Optional character vector of observed variables. Overrides
#'   the file's `varobs` declaration.
#' @param defines Optional named list of macro variables, the equivalent of
#'   Dynare's `-D` command-line option (e.g. `list(N = 3)`).
#'
#' @return An object of class `"dsge_dynare"`, a list with components:
#' \describe{
#'   \item{model}{The translated `dsgenl_model`. Parameters listed in
#'     `estimated_params` are free (with starting values); all other
#'     parameters are fixed at their calibrated values.}
#'   \item{params}{Named numeric vector of calibrated parameter values.}
#'   \item{shock_sd}{Named numeric vector of standard deviations of the
#'     model's shocks, including measurement errors (0 for shocks with no
#'     declared variance, as in Dynare).}
#'   \item{priors}{Named list of [prior()] objects translated from
#'     `estimated_params`, ready for [bayes_dsge()], or `NULL`.}
#'   \item{estimated_params}{Data frame describing each `estimated_params`
#'     entry and how it was translated.}
#'   \item{observed}{Observed variables (Dynare names).}
#'   \item{measurement_errors, data_map}{Observed variables with a
#'     measurement error, and the model variable (`y_obs`) each is mapped
#'     to.}
#'   \item{variables, shocks, shocks_det, parameters}{Names declared in the
#'     file.}
#'   \item{shock_paths}{Period-by-shock matrix of deterministic shock values
#'     from `shocks` blocks with `periods`/`values`, or `NULL`.}
#'   \item{aux}{Data frame of auxiliary variables created for leads, lags
#'     and measurement errors.}
#'   \item{policy}{Optimal-policy problem (`ramsey`, `discretion` or `osr`),
#'     or `NULL`.}
#'   \item{occbin}{Occasionally binding constraints and regime models, or
#'     `NULL`.}
#'   \item{estimation}{Options of the file's `estimation` command used by
#'     [estimate()] and [bayes_dsge()] (`presample`, `first_obs`, `nobs`).}
#'   \item{commands}{List of Dynare commands found in the file (such as
#'     `stoch_simul` or `estimation`), recorded but not executed.}
#'   \item{notes}{Character vector of translation notes, including
#'     anything that was ignored or approximated.}
#' }
#'
#' @details
#' **Declarations and model.** `var`, `varexo`, `varexo_det`, `parameters`,
#' `predetermined_variables`, `varobs`, parameter assignments (and
#' constants assigned to undeclared names), the `model` block (including
#' `model(linear)`, equation tags and `#` model-local variables), leads and
#' lags of any length on variables and shocks, and `STEADY_STATE(x)`,
#' which always equals the steady state of `x` at the current parameter
#' values.
#'
#' **Macro processor.** `@#define`, `@#if`/`@#elseif`/`@#else`/`@#endif`,
#' `@#ifdef`, `@#ifndef`, `@#for` (over arrays, ranges and tuples, with
#' optional `when` filters), `@#include`, `@#includepath`, `@#echo`,
#' `@#error`, simple macro functions and `@{...}` substitution are
#' expanded before the model is translated.
#'
#' **Timing.** No manual re-timing is needed. Every Dynare variable
#' becomes a control; lags become auxiliary state variables named
#' `x_lag1`, `x_lag2`, ...; leads beyond one period become auxiliary
#' controls `x_lead1`, ...; and each shock becomes an exogenous state that
#' holds the current innovation. Impulse responses therefore have the
#' same timing as in Dynare. The auxiliary variables also appear in
#' solution and IRF output.
#'
#' **Steady state and shocks.** `steady_state_model` becomes the model's
#' steady-state function (variables it does not set keep their `initval`
#' value, 0 by default, as in Dynare) and `initval` supplies starting
#' values for the numerical solver. Models declared `model(linear)` are
#' linearised with an exact Jacobian. In the `shocks` block, standard deviations,
#' variances, covariances and correlations are supported; correlated shocks
#' are orthogonalised by Cholesky factorisation in `varexo` order, which
#' reproduces Dynare's impulse responses. Deterministic paths
#' (`periods`/`values`) are returned in `shock_paths`.
#'
#' **Observed variables and measurement errors.** A model may have fewer
#' observed variables than shocks. A `stderr` on an observed endogenous
#' variable (in `shocks` or `estimated_params`) is a measurement error: the
#' variable `y` is observed as `y_obs = y + y_me`, where `y_me` is an
#' i.i.d. shock. [estimate()] and [bayes_dsge()] rename a data column `y`
#' to `y_obs` automatically. With more observed variables than shocks and
#' measurement errors, the likelihood would be singular, so the extra
#' variables are dropped with a note.
#'
#' **Priors.** Dynare's mean/standard-deviation prior parameterisation is
#' converted to dsge's, exactly: `normal_pdf`, `beta_pdf`, `gamma_pdf`,
#' `uniform_pdf`, `inv_gamma2_pdf`, and `inv_gamma_pdf` / `inv_gamma1_pdf`
#' (a prior on a standard deviation, translated to the `"inv_gamma1"`
#' distribution of [prior()] with Dynare's own parameterisation). Shape
#' names are case-insensitive. Shifted or generalised priors and
#' `weibull_pdf` are not translated and are listed in `notes`.
#'
#' **Estimation options.** `presample`, `first_obs` and `nobs` from the
#' file's `estimation` command are stored in `estimation` and used by
#' [estimate()] and [bayes_dsge()]: the data are restricted to the
#' estimation sample and the first `presample` observations only
#' initialise the Kalman filter. dsge always initialises the filter at the
#' stationary distribution (Dynare's `lik_init = 1`); other `lik_init`
#' values are reported in `notes`.
#'
#' **Optimal policy.** With `planner_objective` and `ramsey_model` or
#' `ramsey_policy`, the planner's first-order conditions are derived
#' symbolically and added to the model together with Lagrange multipliers
#' `MULT_1`, ..., as in Dynare, so [solve_dsge()] returns the Ramsey
#' equilibrium; the steady state is found with the multipliers concentrated
#' out. With `discretionary_policy` (linear models), the time-consistent
#' rule is computed with the Dennis (2007) algorithm at the calibrated
#' parameters and the model is closed with it. With `osr_params`,
#' `osr_params_bounds` and `optim_weights`, [osr()] can be called directly
#' on the imported model.
#'
#' **OccBin.** Equations tagged `bind = 'c'` / `relax = 'c'` and an
#' `occbin_constraints` block define occasionally binding constraints;
#' [simulate_occbin()] on the imported model solves them with the
#' piecewise-linear algorithm of Guerrieri and Iacoviello (2015), as
#' Dynare's `occbin_solver` does, using the file's `shocks(surprise)`
#' block by default.
#'
#' **Not supported:** `external_function`, `trend_var`, `EXPECTATION()`,
#' `diff()`, `adl()`, PAC and VAR expectation operators; these raise an
#' error. Other blocks and commands are recorded but not run.
#'
#' @references
#' Dennis, R. (2007). Optimal policy in rational expectations models: new
#' solution algorithms. \emph{Macroeconomic Dynamics}, 11(1), 31-55.
#'
#' Guerrieri, L. and Iacoviello, M. (2015). OccBin: A toolkit for solving
#' dynamic models with occasionally binding constraints easily.
#' \emph{Journal of Monetary Economics}, 70, 22-38.
#'
#' @seealso [dsgenl_model()], [solve_dsge()], [bayes_dsge()], [osr()],
#'   [simulate_occbin()]
#'
#' @examples
#' rbc <- read_dynare(system.file("examples", "rbc.mod", package = "dsge"))
#' rbc
#' sol <- solve_dsge(rbc)
#' irf(sol, periods = 20, se = FALSE)
#'
#' # Model code can also be passed as text, including macro directives
#' ar <- read_dynare(text = "
#'   @#define lags = 2
#'   var y;
#'   varexo e;
#'   parameters rho;
#'   rho = 0.5;
#'   model;
#'     y = e
#'   @#for k in 1:lags
#'       + rho^@{k} * y(-@{k})
#'   @#endfor
#'     ;
#'   end;
#'   shocks;
#'     var e; stderr 0.01;
#'   end;
#' ")
#' solve_dsge(ar)
#'
#' @export
read_dynare <- function(file, text = NULL, observed = NULL, defines = NULL) {
  if (is.null(text)) {
    if (missing(file) || !is.character(file) || length(file) != 1L) {
      stop("Supply either `file` (a path) or `text`.", call. = FALSE)
    }
    if (!file.exists(file)) {
      stop("File not found: ", file, call. = FALSE)
    }
    src <- readLines(file, warn = FALSE)
    source_name <- file
    src_dir <- dirname(file)
  } else {
    src <- as.character(text)
    source_name <- "<text>"
    src_dir <- getwd()
  }

  src <- dyn_strip_comments(paste(src, collapse = "\n"))
  if (grepl("(^|\n)\\s*@#", src) || grepl("@\\{", src) ||
      length(defines) > 0L) {
    src <- dyn_macro_expand(src, dir = src_dir,
                            defines = as.list(defines))
  }
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
    endo = character(0), exo = character(0), exo_det = character(0),
    params = character(0),
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

    if (kw %in% c("var", "varexo", "varexo_det", "parameters",
                  "predetermined_variables", "varobs",
                  "model_local_variable") &&
        !grepl("^[A-Za-z_]+\\s*=", st)) {
      names_found <- dyn_declared_names(st, kw)
      if (kw == "var") {
        if (grepl("^var\\s*\\([^)]*\\blog\\b", st)) {
          stop("var(log) declarations are not supported.", call. = FALSE)
        }
        p$endo <- c(p$endo, names_found)
      } else if (kw == "varexo") {
        p$exo <- c(p$exo, names_found)
      } else if (kw == "varexo_det") {
        p$exo_det <- c(p$exo_det, names_found)
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

    if (kw %in% c("trend_var", "log_trend_var",
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
  all_names <- c(p$endo, p$exo, p$exo_det, p$params)
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
  unsupported <- c("EXPECTATION", "diff", "adl",
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
  stoch <- p$exo
  exo <- c(p$exo, p$exo_det)

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
  meq <- dyn_model_equations(p$model, p$predetermined, endo, exo)
  eqs <- meq$equations
  ss_links <- meq$ss_links

  # --- shocks (correlated shocks are orthogonalised in the equations) -----
  shk <- dyn_parse_shocks(p$blocks$shocks, exo, cal_env, endo = endo)
  notes <- c(notes, shk$notes)
  if (length(shk$cross) > 0L) {
    orth <- dyn_orthogonalise(shk$sd[stoch], shk$cross, stoch)
    eqs <- vapply(eqs, function(eq) {
      dyn_rewrite_ids(eq, function(name, lag, has_index) {
        if (!name %in% names(orth$subst)) return(NULL)
        parts <- vapply(orth$subst[[name]], function(t) {
          dyn_rewrite_ids(t, function(n2, l2, h2) {
            if (n2 %in% stoch) dyn_timed(n2, lag) else NULL
          })
        }, character(1))
        paste0("(", paste(parts, collapse = " + "), ")")
      })
    }, character(1), USE.NAMES = FALSE)
    shk$sd[stoch] <- orth$sd
    notes <- c(notes, paste0(
      "Correlated shocks orthogonalised by Cholesky factorisation in the ",
      "order ", paste(stoch, collapse = ", "), " (as in Dynare's IRFs); ",
      "shock standard deviations refer to the orthogonal shocks."))
  }
  if (length(p$exo_det) > 0L) {
    notes <- c(notes, paste0(
      "varexo_det ", paste(p$exo_det, collapse = ", "), " imported as ",
      "exogenous variables with zero variance; give them paths with ",
      "perfect_foresight()."))
  }

  # --- OccBin: keep the relaxed equations, set aside the binding ones -----
  occ_split <- dyn_occbin_split(eqs, meq$tagged$tags)
  if (!is.null(occ_split)) {
    for (k in seq_along(occ_split$alt)) {
      occ_split$alt[[k]]$equation <- eqs[setdiff(seq_along(eqs),
                                                 occ_split$base_idx)[k]]
    }
    eqs <- eqs[occ_split$base_idx]
  }

  # --- optimal policy (Ramsey FOCs or discretionary rule) ----------------
  spec <- dyn_policy_spec(p)
  if (!is.null(occ_split) && identical(spec$type, "ramsey")) {
    stop("OccBin constraints combined with Ramsey policy are not supported.",
         call. = FALSE)
  }
  policy_extra <- list()
  mults <- character(0)
  if (identical(spec$type, "ramsey")) {
    ram <- dyn_ramsey_equations(eqs, spec$objective, endo, exo, spec$discount)
    eqs <- ram$equations
    mults <- ram$multipliers
    policy_extra$multipliers <- mults
    notes <- c(notes, paste0(
      "Ramsey problem: ", length(mults), " Lagrange multipliers (",
      paste(mults, collapse = ", "), ") and the planner's first-order ",
      "conditions were added to the model."))
  } else if (identical(spec$type, "discretion")) {
    rule <- dyn_discretion_rule(eqs, spec$objective, endo, exo,
                                spec$instruments, spec$discount, params,
                                cal_env)
    eqs <- c(eqs, dyn_rule_equations(rule, endo, exo))
    policy_extra$rule <- rule
    notes <- c(notes, paste0(
      "Discretionary policy: the model is closed with the time-consistent ",
      "rule for ", paste(spec$instruments, collapse = ", "), ", computed at ",
      "the calibrated parameters (Dennis 2007)."))
  }
  endo_model <- c(endo, mults)

  # --- estimated parameters and priors -----------------------------------
  est <- dyn_parse_estimated(p$blocks$estimated_params,
                             p$blocks$estimated_params_init,
                             p$params, exo, cal_env, endo = endo)
  notes <- c(notes, est$notes)

  # --- leads and lags -----------------------------------------------------
  timed <- c(endo_model, exo)
  max_lead <- stats::setNames(integer(length(timed)), timed)
  max_lag <- stats::setNames(integer(length(timed)), timed)
  occ_eqs <- if (is.null(occ_split)) character(0) else
    vapply(occ_split$alt, `[[`, "", "equation")
  for (eq in c(eqs, occ_eqs)) {
    dyn_rewrite_ids(eq, function(name, lag, has_index) {
      if (name %in% timed) {
        if (lag > max_lead[name]) max_lead[name] <<- lag
        if (-lag > max_lag[name]) max_lag[name] <<- -lag
      }
      NULL
    })
  }

  lag_name <- function(v, j) paste0(v, "_lag", j)
  lead_name <- function(v, j) paste0(v, "_lead", j)

  rewrite_timing <- function(name, lag, has_index) {
    if (name %in% timed) {
      if (lag == 0L) return(name)
      if (lag == 1L) return(paste0(name, "(+1)"))
      if (lag > 1L) return(paste0(lead_name(name, lag - 1L), "(+1)"))
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
  add_aux <- function(nm, base, type, shift) {
    aux <<- rbind(aux, data.frame(name = nm, base = base, type = type,
                                  shift = shift, stringsAsFactors = FALSE))
  }
  aux_ctrl_eqs <- character(0)
  lead_ctrls <- character(0)
  for (v in timed) {
    if (max_lead[v] >= 2L) {
      for (j in seq_len(max_lead[v] - 1L)) {
        nm <- lead_name(v, j)
        prev <- if (j == 1L) v else lead_name(v, j - 1L)
        aux_ctrl_eqs <- c(aux_ctrl_eqs, paste0(nm, " = ", prev, "(+1)"))
        lead_ctrls <- c(lead_ctrls, nm)
        add_aux(nm, v, "lead", j)
      }
    }
  }
  aux_state_eqs <- character(0)
  lag_states <- character(0)
  for (v in timed) {
    if (max_lag[v] >= 1L) {
      for (j in seq_len(max_lag[v])) {
        nm <- lag_name(v, j)
        prev <- if (j == 1L) v else lag_name(v, j - 1L)
        aux_state_eqs <- c(aux_state_eqs, paste0(nm, "(+1) = ", prev))
        lag_states <- c(lag_states, nm)
        add_aux(nm, v, "lag", -j)
      }
    }
  }

  # --- observed variables and measurement errors -------------------------
  obs <- if (!is.null(observed)) observed else p$varobs
  if (is.null(obs)) obs <- character(0)
  unknown_obs <- setdiff(obs, endo)
  if (length(unknown_obs) > 0L) {
    stop("Observed variable(s) not declared with 'var': ",
         paste(unknown_obs, collapse = ", "), call. = FALSE)
  }
  me_est <- est$table$name[est$table$type == "me"]
  me_vars <- unique(c(names(shk$me_sd), me_est))
  dropped_me <- setdiff(me_vars, obs)
  if (length(dropped_me) > 0L) {
    notes <- c(notes, paste0(
      "Measurement error on non-observed variable(s) ignored: ",
      paste(dropped_me, collapse = ", "), "."))
  }
  me_vars <- intersect(obs, me_vars)
  me_states <- dyn_suffix(me_vars, "_me")
  me_obs <- dyn_suffix(me_vars, "_obs")
  me_eqs <- if (length(me_vars) > 0L) {
    c(paste0(me_obs, " = ", me_vars, " + ", me_states),
      paste0(me_states, "(+1) = 0"))
  } else {
    character(0)
  }
  for (k in seq_along(me_vars)) add_aux(me_obs[k], me_vars[k], "obs", 0L)
  data_map <- stats::setNames(me_obs, me_vars)
  obs_model <- obs
  obs_model[obs %in% me_vars] <- data_map[obs[obs %in% me_vars]]

  clash <- intersect(c(aux$name, me_states, names(ss_links)),
                     c(endo, exo, p$params))
  if (length(clash) > 0L) {
    stop("Auxiliary variable name(s) clash with declared names: ",
         paste(clash, collapse = ", "), call. = FALSE)
  }

  exo_all <- c(exo, me_states)
  n_shocks <- length(exo_all)
  if (length(obs_model) > n_shocks) {
    kept <- obs_model[seq_len(n_shocks)]
    notes <- c(notes, paste0(
      "More observed variables (", length(obs_model), ") than shocks (",
      n_shocks, ") would make the likelihood singular; only '",
      paste(kept, collapse = "', '"), "' kept as observed. Add ",
      "measurement errors to observe more variables."))
    obs_model <- kept
  }
  if (length(obs_model) == 0L && length(exo_all) > 0L) {
    notes <- c(notes, paste0(
      "No varobs: the model has no observed variables. Declare varobs or ",
      "pass `observed` to estimate it."))
  }

  # --- steady state ------------------------------------------------------
  exo_ss <- stats::setNames(numeric(length(exo_all)), exo_all)
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
  shock_state_eqs <- if (length(exo) > 0L) {
    paste0(exo, "(+1) = ", format_num(exo_ss[exo]))
  } else {
    character(0)
  }

  default_guess <- if (p$model_linear) 0 else 1
  base_guess <- stats::setNames(rep(default_guess, length(endo_model)),
                                endo_model)
  base_guess[mults] <- 0
  for (v in intersect(names(init_vals), endo)) base_guess[v] <- init_vals[[v]]
  fill_aux <- function(base_vals) {
    out <- c(base_vals, exo_ss)
    for (k in seq_len(nrow(aux))) out[aux$name[k]] <- out[[aux$base[k]]]
    out
  }
  ss_guess <- fill_aux(base_guess)

  ss_function <- NULL
  ss_defaults <- stats::setNames(numeric(length(endo)), endo)
  for (v in intersect(names(init_vals), endo)) ss_defaults[v] <- init_vals[[v]]
  if (!is.null(p$blocks$steady_state_model)) {
    if (length(mults) > 0L) {
      # Ramsey: the file's steady state only covers the original variables;
      # use it as the starting guess for the augmented system.
      ssf <- dyn_ss_function(p$blocks$steady_state_model, cal_env, endo,
                             exo_ss, function(v) v, defaults = ss_defaults)
      guess <- tryCatch(ssf(params), error = function(e) NULL)
      if (!is.null(guess)) {
        base_guess[endo] <- guess[endo]
        ss_guess <- fill_aux(base_guess)
      }
    } else {
      ss_function <- dyn_ss_function(p$blocks$steady_state_model, cal_env,
                                     endo, exo_ss, fill_aux,
                                     defaults = ss_defaults)
    }
  }

  # Ramsey: steady state with the multipliers concentrated out
  ramsey_holder <- NULL
  if (length(mults) > 0L) {
    ramsey_holder <- new.env(parent = emptyenv())
    ss_function <- function(params) {
      dyn_ramsey_ss(ramsey_holder$model, params, ss_guess, endo, mults,
                    fill_aux)
    }
  }

  # --- assemble dsgenl_model ---------------------------------------------
  controls <- c(endo_model, lead_ctrls, me_obs)
  unobs <- setdiff(controls, obs_model)
  all_eqs <- c(model_eqs, aux_ctrl_eqs, me_eqs[seq_along(me_obs)],
               shock_state_eqs, me_eqs[length(me_obs) + seq_along(me_states)],
               aux_state_eqs)
  make_model <- function(eq_set, fixed, start) {
    do.call(dsgenl_model, c(
      as.list(eq_set),
      list(observed = obs_model, unobserved = unobs, exo_state = exo_all,
           endo_state = lag_states, fixed = fixed, start = start,
           ss_guess = ss_guess, ss_function = ss_function)
    ))
  }
  probe <- make_model(all_eqs, list(), list())
  model_params <- setdiff(probe$parameters, names(ss_links))

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
  }
  start <- list()
  for (nm in est_names) {
    row <- est$table$name == nm & est$table$type == "param"
    init <- est$table$init[row]
    if (is.na(init)) init <- params[nm]
    if (is.na(init)) init <- est$table$prior_mean[row]
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

  model <- make_model(all_eqs, fixed, start)
  model$linear <- p$model_linear
  if (length(ss_links) > 0L) {
    ss_eqs <- vapply(all_eqs, function(eq) {
      dyn_rewrite_ids(eq, function(name, lag, has_index) {
        if (name %in% names(ss_links)) ss_links[[name]] else NULL
      })
    }, character(1), USE.NAMES = FALSE)
    model <- dyn_link_steady_state(model, make_model(ss_eqs, fixed, start),
                                   ss_links)
  }

  if (!is.null(ramsey_holder)) ramsey_holder$model <- model

  occbin <- dyn_occbin_build(
    occ_split, p$blocks$occbin_constraints,
    function(eq) dyn_rewrite_ids(eq, rewrite_timing),
    function(eq_set) {
      prm <- setdiff(make_model(eq_set, list(), list())$parameters,
                     names(ss_links))
      miss <- setdiff(prm, c(names(params), est_names))
      if (length(miss) > 0L) {
        stop("No value assigned to parameter(s): ",
             paste(miss, collapse = ", "), call. = FALSE)
      }
      make_model(eq_set, as.list(params[setdiff(prm, est_names)]),
                 start[intersect(names(start), prm)])
    }, all_eqs)

  # Calibrated values of estimated parameters complete the calibration
  for (nm in est_names) {
    if (!nm %in% names(params)) params[nm] <- start[[nm]]
  }

  # --- shock standard deviations -----------------------------------------
  shock_sd <- c(shk$sd, stats::setNames(numeric(length(me_states)), me_states))
  for (k in seq_along(me_vars)) {
    v <- shk$me_sd[me_vars[k]]
    if (!is.na(v)) shock_sd[me_states[k]] <- v
  }
  for (k in seq_len(nrow(est$table))) {
    ty <- est$table$type[k]
    nm <- est$table$name[k]
    key <- if (ty == "stderr") nm else if (ty == "me") paste0(nm, "_me") else NA
    if (is.na(key) || !key %in% names(shock_sd)) next
    declared <- if (ty == "stderr") shk$declared[nm] else shk$me_sd[nm]
    if (is.na(declared)) {
      v <- est$table$init[k]
      if (is.na(v)) v <- est$table$prior_mean[k]
      if (!is.na(v) && is.finite(v)) shock_sd[key] <- v
    }
  }
  if (!is.null(est$priors)) {
    keep <- names(est$priors)
    me_keys <- dyn_suffix("sd_e.", me_est, prefix = TRUE)
    for (k in which(keep %in% me_keys)) {
      v <- sub("^sd_e\\.", "", keep[k])
      keep[k] <- if (v %in% me_vars) paste0("sd_e.", v, "_me") else NA
    }
    names(est$priors) <- keep
    est$priors <- est$priors[!is.na(keep)]
    if (length(est$priors) == 0L) est$priors <- NULL
  }

  zero_sd <- intersect(names(shock_sd)[shock_sd == 0], stoch)
  if (length(zero_sd) > 0L && length(zero_sd) < length(stoch)) {
    notes <- c(notes, paste0("No variance declared for shock(s) ",
                             paste(zero_sd, collapse = ", "),
                             "; standard deviation set to 0 as in Dynare."))
  } else if (length(zero_sd) > 0L) {
    notes <- c(notes, paste0("No shock variances declared; all standard ",
                             "deviations are 0 as in Dynare. Pass shock_sd ",
                             "to solve_dsge() to simulate."))
  }

  estimation <- dyn_estimation_options(p$commands, cal_env)
  notes <- c(notes, estimation$notes)
  estimation$notes <- NULL

  policy <- dyn_policy_info(p, cal_env, endo, params, spec = spec,
                            extra = policy_extra)
  notes <- c(notes, policy$notes)

  ignored_blocks <- setdiff(names(p$blocks), c(
    "initval", "initval_opts", "steady_state_model",
    "steady_state_model_opts", "shocks", "shocks_opts", "estimated_params",
    "estimated_params_init", "osr_params_bounds", "optim_weights",
    "occbin_constraints"))
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
    observed = obs[seq_along(obs_model)],
    data_map = data_map,
    measurement_errors = me_vars,
    variables = endo,
    shocks = stoch,
    shocks_det = p$exo_det,
    shock_paths = shk$det,
    parameters = p$params,
    aux = aux,
    policy = policy$policy,
    occbin = occbin,
    estimation = estimation,
    commands = p$commands,
    notes = unique(notes)
  )
}

#' Make STEADY_STATE(x) references follow the steady state
#'
#' `model` uses internal parameters (named in `links`) for STEADY_STATE(x);
#' `ss_model` is the same model with those references replaced by x itself,
#' which has the same steady state. The evaluation function fills the
#' internal parameters from ss_model's steady state at the current
#' parameter values (cached), so every solver can use the model unchanged.
#' @noRd
dyn_link_steady_state <- function(model, ss_model, links) {
  inner <- model$eval_fn
  param_names <- ss_model$parameters
  cache <- new.env(parent = emptyenv())
  model$eval_fn <- function(values) {
    pv <- values[param_names]
    if (is.null(cache$key) || !identical(cache$key, pv)) {
      ss <- steady_state(ss_model, params = pv)
      cache$key <- pv
      cache$vals <- ss$values[unname(links)]
    }
    values[names(links)] <- cache$vals
    inner(values)
  }
  model$parameters <- setdiff(model$parameters, names(links))
  model$free_parameters <- setdiff(model$free_parameters, names(links))
  model$ss_links <- links
  model
}

#' Model block statements -> equation strings in Dynare timing
#' @noRd
dyn_model_equations <- function(statements, predetermined, endo, exo) {
  locals <- list()
  eqs <- character(0)
  tags <- list()
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
    tag <- list()
    if (grepl("^\\[", st)) {
      tag <- dyn_parse_tags(sub("^\\[([^]]*)\\].*$", "\\1", st))
      st <- trimws(sub("^\\[[^]]*\\]\\s*", "", st))
    }
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
    tags[[length(eqs)]] <- tag
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

  # STEADY_STATE(x) -> internal parameter linked to x's steady state
  ss_links <- character(0)
  eqs <- vapply(eqs, function(eq) {
    repeat {
      pos <- regexpr("(?<![A-Za-z0-9_.])STEADY_STATE\\s*\\(", eq, perl = TRUE)
      if (pos == -1L) break
      open <- pos + attr(pos, "match.length") - 1L
      arg <- trimws(dyn_paren_content(eq, open))
      if (!arg %in% c(endo, exo)) {
        stop("STEADY_STATE() must be applied to a declared variable: ",
             arg, call. = FALSE)
      }
      link <- paste0(arg, "__ss")
      ss_links[link] <<- arg
      eq <- paste0(substr(eq, 1L, pos - 1L), link,
                   substr(eq, open + nchar(arg) + 2L, nchar(eq)))
    }
    eq
  }, character(1), USE.NAMES = FALSE)

  equations <- vapply(eqs, function(eq) {
    eq <- dyn_translate_math(eq)
    n_eq <- lengths(regmatches(eq, gregexpr("=", eq, fixed = TRUE)))
    if (n_eq == 0L) {
      eq <- paste0(eq, " = 0")
    } else if (n_eq > 1L) {
      stop("Equation has more than one '=': ", eq, call. = FALSE)
    }
    eq
  }, character(1), USE.NAMES = FALSE)
  list(equations = equations, ss_links = ss_links,
       tagged = list(equations = equations, tags = tags))
}

#' Parse an equation tag list such as name='x', bind='y < 0'
#' @noRd
dyn_parse_tags <- function(txt) {
  m <- gregexpr("([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*('[^']*'|\"[^\"]*\")", txt,
                perl = TRUE)
  items <- regmatches(txt, m)[[1]]
  out <- list()
  for (it in items) {
    key <- trimws(sub("=.*$", "", it))
    val <- trimws(sub("^[^=]*=", "", it))
    out[[key]] <- substr(val, 2L, nchar(val) - 1L)
  }
  out
}

#' paste0 that keeps zero-length input zero-length
#' @noRd
dyn_suffix <- function(x, add, prefix = FALSE) {
  if (prefix) {
    if (length(add) == 0L) character(0) else paste0(x, add)
  } else {
    if (length(x) == 0L) character(0) else paste0(x, add)
  }
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
dyn_ss_function <- function(statements, cal_env, endo, exo_ss, fill_aux,
                            defaults = NULL) {
  assigns <- lapply(statements, function(st) {
    if (!grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*=[^=]", st)) {
      stop("Unsupported statement in steady_state_model: ", st,
           call. = FALSE)
    }
    list(name = trimws(sub("=.*$", "", st)),
         expr = parse(text = dyn_translate_math(sub("^[^=]*=", "", st))))
  })
  # As in Dynare, variables not set in steady_state_model keep their
  # initval value (0 by default).
  if (is.null(defaults)) defaults <- stats::setNames(numeric(length(endo)), endo)
  force(defaults)
  force(cal_env)
  force(exo_ss)
  force(fill_aux)
  function(params) {
    env <- new.env(parent = baseenv())
    for (nm in ls(cal_env)) assign(nm, get(nm, envir = cal_env), envir = env)
    for (nm in names(exo_ss)) assign(nm, exo_ss[[nm]], envir = env)
    for (nm in names(params)) assign(nm, params[[nm]], envir = env)
    for (nm in endo) assign(nm, defaults[[nm]], envir = env)
    for (a in assigns) {
      assign(a$name, eval(a$expr, envir = env), envir = env)
    }
    fill_aux(unlist(mget(endo, envir = env)))
  }
}

#' Parse the shocks block (variances only)
#' @noRd
dyn_parse_shocks <- function(statements, exo, cal_env, endo = character(0)) {
  sd <- stats::setNames(numeric(length(exo)), exo)
  declared <- stats::setNames(rep(NA_real_, length(exo)), exo)
  me_sd <- numeric(0)
  notes <- character(0)
  cross <- list()
  det <- list()
  current <- NULL
  check_name <- function(nm) {
    if (!nm %in% c(exo, endo)) {
      stop("'", nm, "' in the shocks block is not declared with 'varexo' ",
           "or 'var'.", call. = FALSE)
    }
  }
  set_sd <- function(nm, value) {
    if (nm %in% exo) {
      sd[nm] <<- value
      declared[nm] <<- value
    } else {
      me_sd[nm] <<- value
    }
  }
  pair_names <- function(txt) {
    strsplit(trimws(txt), "[[:space:],]+")[[1]]
  }
  for (st in statements) {
    if (grepl("^var\\s", st)) {
      body <- sub("^var\\s+", "", st)
      if (grepl("=", body)) {
        lhs <- pair_names(sub("=.*$", "", body))
        for (nm in lhs) check_name(nm)
        val <- dyn_eval(sub("^[^=]*=", "", body), cal_env, "shock variance")
        if (length(lhs) == 1L) {
          set_sd(lhs, sqrt(val))
        } else if (length(lhs) == 2L && all(lhs %in% exo)) {
          cross[[length(cross) + 1L]] <- list(a = lhs[1], b = lhs[2],
                                              type = "cov", value = val)
        } else if (length(lhs) == 2L) {
          notes <- c(notes, paste0("Measurement-error covariance ignored: ",
                                   st, "."))
        } else {
          stop("Cannot parse shocks statement: ", st, call. = FALSE)
        }
        current <- NULL
      } else {
        current <- trimws(body)
        check_name(current)
      }
    } else if (grepl("^stderr\\s", st)) {
      if (is.null(current)) {
        stop("'stderr' without a preceding 'var' in the shocks block.",
             call. = FALSE)
      }
      set_sd(current, dyn_eval(sub("^stderr\\s+", "", st), cal_env,
                               "shock standard deviation"))
    } else if (grepl("^corr\\s", st)) {
      lhs <- pair_names(sub("=.*$", "", sub("^corr\\s+", "", st)))
      if (length(lhs) != 2L || !grepl("=", st)) {
        stop("Cannot parse shocks statement: ", st, call. = FALSE)
      }
      for (nm in lhs) check_name(nm)
      if (all(lhs %in% exo)) {
        cross[[length(cross) + 1L]] <- list(
          a = lhs[1], b = lhs[2], type = "corr",
          value = dyn_eval(sub("^[^=]*=", "", st), cal_env,
                           "shock correlation"))
      } else {
        notes <- c(notes, paste0("Measurement-error correlation ignored: ",
                                 st, "."))
      }
    } else if (grepl("^periods\\s", st)) {
      if (is.null(current)) {
        stop("'periods' without a preceding 'var' in the shocks block.",
             call. = FALSE)
      }
      det[[length(det) + 1L]] <- list(
        shock = current,
        periods = dyn_parse_periods(sub("^periods\\s+", "", st)))
    } else if (grepl("^values\\s", st)) {
      k <- length(det)
      if (k == 0L || det[[k]]$shock != current || !is.null(det[[k]]$values)) {
        stop("'values' without a preceding 'periods' in the shocks block.",
             call. = FALSE)
      }
      vals <- vapply(dyn_split_values(sub("^values\\s+", "", st)),
                     dyn_eval, numeric(1), env = cal_env,
                     what = "deterministic shock value")
      det[[k]]$values <- vals
    } else {
      notes <- c(notes, paste0("Unrecognised shocks statement ignored: ",
                               st, "."))
    }
  }
  list(sd = sd, declared = declared, me_sd = me_sd, cross = cross,
       det = dyn_det_paths(det), notes = unique(notes))
}

#' Parse a Dynare periods list such as "1 2:4 7"
#' @noRd
dyn_parse_periods <- function(txt) {
  parts <- strsplit(trimws(txt), "[[:space:],]+")[[1]]
  lapply(parts, function(p) {
    if (grepl(":", p)) {
      ab <- as.integer(strsplit(p, ":", fixed = TRUE)[[1]])
      seq(ab[1], ab[2])
    } else {
      as.integer(p)
    }
  })
}

#' Split a values list: numbers separated by spaces/commas, or
#' parenthesised expressions
#' @noRd
dyn_split_values <- function(txt) {
  txt <- trimws(txt)
  toks <- regmatches(txt, gregexpr("\\((?:[^()]|\\([^()]*\\))*\\)|[^[:space:],()]+",
                                   txt, perl = TRUE))[[1]]
  toks[nzchar(toks)]
}

#' Collect deterministic shock paths into a period-by-shock matrix
#' @noRd
dyn_det_paths <- function(det) {
  if (length(det) == 0L) return(NULL)
  shocks <- unique(vapply(det, `[[`, "", "shock"))
  horizon <- max(unlist(lapply(det, `[[`, "periods")))
  out <- matrix(0, horizon, length(shocks),
                dimnames = list(NULL, shocks))
  for (d in det) {
    if (is.null(d$values) || length(d$values) != length(d$periods)) {
      if (!is.null(d$values) && length(d$values) == 1L) {
        d$values <- rep(d$values, length(d$periods))
      } else {
        stop("Deterministic shock '", d$shock, "': 'periods' and 'values' ",
             "have different lengths.", call. = FALSE)
      }
    }
    for (k in seq_along(d$periods)) out[d$periods[[k]], d$shock] <- d$values[k]
  }
  out
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
                                cal_env, endo = character(0)) {
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
      if (nm %in% endo) {
        return(list(name = nm, type = "me", rest = fields[-1]))
      }
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
    shape_pos <- which(tolower(rest) %in% names(dyn_prior_shapes))
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
      shape <- dyn_prior_shapes[[tolower(rest[sp])]]
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
        key <- if (tg$type %in% c("stderr", "me")) paste0("sd_e.", tg$name)
               else tg$name
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
  if (shape == "inv_gamma1") {
    ig <- inv_gamma1_from_moments(mean, sd)
    out$prior <- prior("inv_gamma1", s = ig$s, nu = ig$nu)
    out$translation <- sprintf("inv_gamma1(s = %g, nu = %g)", ig$s, ig$nu)
    return(out)
  }
  # inv_gamma2_pdf: inverse gamma on the parameter itself, dsge density
  # x^-(a+1) exp(-b/x), with the same moment matching as Dynare.
  a <- if (is.finite(sd)) 2 + mean^2 / sd^2 else 2
  b <- mean * (a - 1)
  out$prior <- prior("inv_gamma", shape = a, scale = b)
  out$translation <- sprintf("inv_gamma(shape = %g, scale = %g)", a, b)
  out
}

#' Unwrap a dsge_dynare object passed to a model-taking function
#' @noRd
unwrap_dynare <- function(x) {
  if (inherits(x, "dsge_dynare")) x$model else x
}

#' Free fixed parameters of a dsgenl_model so that supplied values are used
#'
#' dsgenl models evaluate fixed parameters after supplied ones, so a value
#' passed for a fixed parameter would be ignored. For imported Dynare models
#' every non-estimated parameter is fixed; when the caller supplies values
#' for some of them, they are moved to `start`.
#' @noRd
dyn_unfix <- function(model, nms) {
  nms <- intersect(nms, names(model$fixed))
  if (length(nms) == 0L) return(model)
  for (nm in nms) model$start[[nm]] <- model$fixed[[nm]]
  model$fixed[nms] <- NULL
  model$free_parameters <- union(model$free_parameters, nms)
  model
}

#' Options of the file's estimation command used by dsge's estimators
#' @noRd
dyn_estimation_options <- function(commands, cal_env) {
  est <- Filter(function(cm) cm$name == "estimation", commands)
  out <- list(presample = 0L, first_obs = 1L, nobs = NA_integer_,
              datafile = NULL, notes = character(0))
  if (length(est) == 0L) return(out)
  opts <- dyn_split_commas(est[[length(est)]]$options)
  kv <- list()
  for (o in opts) {
    if (!grepl("=", o)) next
    kv[[tolower(trimws(sub("=.*$", "", o)))]] <- trimws(sub("^[^=]*=", "", o))
  }
  int_opt <- function(key, default) {
    if (is.null(kv[[key]]) || grepl("^\\[", kv[[key]])) return(default)
    as.integer(dyn_eval(kv[[key]], cal_env, paste0("estimation option ", key)))
  }
  out$presample <- int_opt("presample", 0L)
  out$first_obs <- int_opt("first_obs", 1L)
  out$nobs <- int_opt("nobs", NA_integer_)
  out$datafile <- kv$datafile
  lik_init <- int_opt("lik_init", 1L)
  if (lik_init != 1L) {
    out$notes <- c(out$notes, paste0(
      "estimation uses lik_init = ", lik_init, "; dsge initialises the ",
      "Kalman filter at the stationary distribution (Dynare's lik_init = 1)."))
  }
  if (int_opt("prefilter", 0L) != 0L) {
    out$notes <- c(out$notes, paste0(
      "estimation uses prefilter = 1 (demeaned data); pass demeaned data ",
      "or remove the constants to reproduce it."))
  }
  out
}

#' @noRd
dyn_estimation_option <- function(x, key) {
  v <- x$estimation[[key]]
  if (is.null(v)) 0L else v
}

#' Restrict data to the file's first_obs / nobs sample
#' @noRd
dyn_estimation_sample <- function(x, data) {
  est <- x$estimation
  if (is.null(est)) return(data)
  first <- if (is.null(est$first_obs)) 1L else est$first_obs
  n <- nrow(data)
  last <- if (is.null(est$nobs) || is.na(est$nobs)) n else first + est$nobs - 1L
  if (first == 1L && last == n) return(data)
  if (last > n) {
    stop("The data have ", n, " rows but the estimation command asks for ",
         "observations ", first, " to ", last, ".", call. = FALSE)
  }
  data[first:last, , drop = FALSE]
}

#' Rename data columns for observed variables with measurement error
#'
#' A variable y observed with measurement error is represented by the
#' observable y_obs = y + y_me; data supplied under Dynare's name y is
#' renamed to y_obs.
#' @noRd
dyn_map_data <- function(x, data) {
  if (!inherits(x, "dsge_dynare") || length(x$data_map) == 0L) return(data)
  if (is.null(colnames(data))) return(data)
  cn <- colnames(data)
  for (v in names(x$data_map)) {
    if (v %in% cn && !x$data_map[[v]] %in% cn) {
      cn[cn == v] <- x$data_map[[v]]
    }
  }
  colnames(data) <- cn
  data
}


