# Nonlinear DSGE equation parser
#
# Parses string-based nonlinear equations of the form "LHS = RHS".
# Handles lead variables via the VAR(+1) notation, replacing them with
# VAR__f for safe R expression parsing.

#' Parse a nonlinear DSGE equation string
#'
#' @param eq_str Character string of the form "LHS = RHS".
#' @return List with parsed components including the R expression,
#'   variable references, and whether it is a state equation.
#' @noRd
parse_nl_equation <- function(eq_str) {
  if (!is.character(eq_str) || length(eq_str) != 1L) {
    stop("Each equation must be a single character string.", call. = FALSE)
  }

  eq_str <- trimws(eq_str)

  if (grepl("==|!=|<=|>=", eq_str)) {
    stop("Use '=' for equations, not '==', '!=', '<=', or '>=': ",
         eq_str, call. = FALSE)
  }

  # Split on '='
  parts <- strsplit(eq_str, "=", fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    stop("Equation must contain exactly one '=' sign: ", eq_str, call. = FALSE)
  }

  lhs_str <- trimws(parts[1])
  rhs_str <- trimws(parts[2])

  if (nchar(lhs_str) == 0L || nchar(rhs_str) == 0L) {
    stop("Both sides of '=' must be non-empty: ", eq_str, call. = FALSE)
  }

  # Form residual: LHS - RHS
  residual_str <- paste0("(", lhs_str, ") - (", rhs_str, ")")

  # Detect and replace lead variables: VAR(+1) -> VAR__f
  lead_pattern <- "([A-Za-z][A-Za-z0-9_.]*)\\(\\+1\\)"
  lead_matches <- regmatches(residual_str,
                             gregexpr(lead_pattern, residual_str))[[1]]
  lead_vars <- unique(sub("\\(\\+1\\)", "", lead_matches))

  safe_str <- residual_str
  for (lv in lead_vars) {
    safe_str <- gsub(paste0(lv, "\\(\\+1\\)"),
                     paste0(lv, "__f"), safe_str)
  }

  # Parse as R expression
  expr <- tryCatch(
    parse(text = safe_str),
    error = function(e) {
      stop("Cannot parse equation as R expression: ", eq_str,
           "\n  ", e$message, call. = FALSE)
    }
  )

  # Extract all symbol names
  all_names <- all.vars(expr)

  lead_refs <- all_names[grepl("__f$", all_names)]
  current_refs <- all_names[!grepl("__f$", all_names)]
  lead_var_names <- sub("__f$", "", lead_refs)

  # Detect state equation: LHS is exactly "VAR(+1)"
  lhs_lead <- regmatches(
    lhs_str,
    regexec("^([A-Za-z][A-Za-z0-9_.]*)\\(\\+1\\)$", lhs_str)
  )[[1]]
  is_state_eq <- length(lhs_lead) == 2L
  lhs_state_var <- if (is_state_eq) lhs_lead[2L] else NULL

  list(
    original = eq_str,
    lhs_str = lhs_str,
    rhs_str = rhs_str,
    safe_str = safe_str,
    expression = expr,
    all_names = all_names,
    current_refs = current_refs,
    lead_var_names = lead_var_names,
    is_state_eq = is_state_eq,
    lhs_state_var = lhs_state_var
  )
}

#' Build evaluation function for a system of nonlinear equations
#'
#' Returns a function that takes a named numeric vector (containing variable
#' values at all timings plus parameter values) and returns the vector of
#' equation residuals.
#'
#' @param parsed_eqs List of parsed equation objects.
#' @return Function(values) -> numeric vector of residuals.
#' @noRd
build_nl_eval_function <- function(parsed_eqs) {
  exprs <- lapply(parsed_eqs, function(eq) eq$expression)
  n <- length(exprs)

  function(values) {
    env <- list2env(as.list(values), parent = baseenv())
    res <- numeric(n)
    for (i in seq_len(n)) {
      res[i] <- eval(exprs[[i]], envir = env)
    }
    res
  }
}
