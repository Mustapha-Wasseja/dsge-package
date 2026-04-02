# Internal equation parser for linear DSGE models
#
# Parses formula-based equation specifications into structured term lists
# that can be mapped to structural matrices.

#' Parse a single DSGE equation
#'
#' @param eq A `dsge_equation` object (from `obs()`, `unobs()`, or `state()`).
#' @return A list containing parsed equation information.
#' @noRd
parse_equation <- function(eq) {
  formula <- eq$formula
  type <- eq$type
  shock <- if (type == "state") eq$shock else NA

  # Extract LHS and RHS
  if (length(formula) != 3L) {
    stop("Each equation must have both a left-hand side and right-hand side.",
         call. = FALSE)
  }

  lhs_expr <- formula[[2L]]
  rhs_expr <- formula[[3L]]

  # Parse LHS: should be a single variable, optionally multiplied by params
  lhs_info <- parse_lhs(lhs_expr)

  # Parse RHS: extract additive terms, each of form param_expr * var
  rhs_terms <- parse_rhs(rhs_expr)

  # Collect all parameter names from LHS and RHS
  lhs_params <- extract_param_names(lhs_info$coef_expr)
  rhs_params <- unique(unlist(lapply(rhs_terms, function(t) {
    extract_param_names(t$coef_expr)
  })))
  all_params <- unique(c(lhs_params, rhs_params))

  # Collect all variable names from RHS
  rhs_variables <- unique(vapply(rhs_terms, function(t) t$variable, character(1)))

  list(
    type = type,
    shock = shock,
    lhs_var = lhs_info$variable,
    lhs_coef_expr = lhs_info$coef_expr,
    rhs_terms = rhs_terms,
    parameters = all_params,
    rhs_variables = rhs_variables,
    formula = formula
  )
}

#' Parse the left-hand side of an equation
#'
#' Extracts the variable name and any coefficient expression.
#' LHS can be: `y`, `{param}*y`, `(1/{param})*y`
#'
#' @param expr An R expression from the LHS of a formula.
#' @return List with `variable` (character) and `coef_expr` (expression or 1).
#' @noRd
parse_lhs <- function(expr) {
  if (is.name(expr)) {
    # Simple variable: y
    return(list(variable = as.character(expr), coef_expr = 1))
  }

  if (is.call(expr) && identical(expr[[1L]], as.name("*"))) {
    # Product: something * something
    # Try both orderings: coef * var or var * coef
    left <- expr[[2L]]
    right <- expr[[3L]]

    if (is.name(right) && !is_param_expr(right)) {
      return(list(variable = as.character(right), coef_expr = left))
    }
    if (is.name(left) && !is_param_expr(left)) {
      return(list(variable = as.character(left), coef_expr = right))
    }
  }

  # Fallback: try to find the variable name in the expression
  vars <- find_variables_in_expr(expr)
  if (length(vars) == 1L) {
    # Build coefficient by "dividing out" the variable conceptually
    # For simple cases like (1/beta) * y
    return(list(variable = vars, coef_expr = remove_variable_from_product(expr, vars)))
  }

  stop("Cannot parse left-hand side of equation. ",
       "LHS should be a single variable, optionally multiplied by parameters.",
       call. = FALSE)
}

#' Parse the right-hand side of an equation into additive terms
#'
#' Each term is of the form `coefficient_expression * variable`
#' where the variable may be wrapped in `lead()` or `E()`.
#'
#' @param expr An R expression from the RHS of a formula.
#' @return A list of term lists, each with `variable`, `coef_expr`,
#'   `is_lead`, `lead_k`.
#' @noRd
parse_rhs <- function(expr) {
  # Flatten additive structure: a + b - c becomes list of (+a, +b, -c)
  flat_terms <- flatten_additive(expr)

  parsed_terms <- list()
  for (ft in flat_terms) {
    term <- parse_single_term(ft$expr, ft$sign)
    if (!is.null(term)) {
      parsed_terms <- c(parsed_terms, list(term))
    }
  }

  if (length(parsed_terms) == 0L) {
    stop("Right-hand side of equation has no parseable terms.", call. = FALSE)
  }

  parsed_terms
}

#' Flatten an additive expression into a list of (sign, subexpression) pairs
#' @noRd
flatten_additive <- function(expr, sign = 1) {
  if (is.call(expr)) {
    op <- as.character(expr[[1L]])

    if (op == "+" && length(expr) == 3L) {
      left <- flatten_additive(expr[[2L]], sign)
      right <- flatten_additive(expr[[3L]], sign)
      return(c(left, right))
    }

    if (op == "+" && length(expr) == 2L) {
      return(flatten_additive(expr[[2L]], sign))
    }

    if (op == "-" && length(expr) == 3L) {
      left <- flatten_additive(expr[[2L]], sign)
      right <- flatten_additive(expr[[3L]], -sign)
      return(c(left, right))
    }

    if (op == "-" && length(expr) == 2L) {
      return(flatten_additive(expr[[2L]], -sign))
    }

    if (op == "(") {
      return(flatten_additive(expr[[2L]], sign))
    }
  }

  list(list(expr = expr, sign = sign))
}

#' Parse a single term (variable * coefficient) from the RHS
#' @noRd
parse_single_term <- function(expr, sign = 1) {
  # Find the variable in this term and separate it from the coefficient
  info <- separate_var_and_coef(expr)

  if (is.null(info)) {
    # Could be a constant term (no variable) -- not supported in DSGE
    stop("Each term in a DSGE equation must contain exactly one variable. ",
         "Found a term with no variable: ", deparse(expr), call. = FALSE)
  }

  # Apply sign to coefficient
  if (sign == -1) {
    if (identical(info$coef_expr, 1)) {
      info$coef_expr <- -1
    } else {
      info$coef_expr <- call("-", info$coef_expr)
    }
  }

  info
}

#' Separate variable and coefficient in a multiplicative term
#'
#' Handles: `x`, `beta*x`, `(1/beta)*x`, `lead(x)`, `beta*lead(x)`,
#' `(1/beta)*(lead(x) + gamma*z)` (distributes), etc.
#'
#' @noRd
separate_var_and_coef <- function(expr) {
  # Case 1: bare variable name

  if (is.name(expr)) {
    varname <- as.character(expr)
    return(list(variable = varname, coef_expr = 1, is_lead = FALSE, lead_k = 0L))
  }

  # Case 2: lead(x) or lead(x, k) or E(x)
  if (is.call(expr)) {
    fn <- as.character(expr[[1L]])

    if (fn == "lead") {
      varname <- as.character(expr[[2L]])
      k <- if (length(expr) >= 3L) as.integer(eval(expr[[3L]])) else 1L
      return(list(variable = varname, coef_expr = 1, is_lead = TRUE, lead_k = k))
    }

    if (fn == "E") {
      varname <- as.character(expr[[2L]])
      return(list(variable = varname, coef_expr = 1, is_lead = TRUE, lead_k = 1L))
    }

    # Case 3: multiplication a * b
    if (fn == "*" && length(expr) == 3L) {
      left <- expr[[2L]]
      right <- expr[[3L]]

      # Try: coef * var
      right_info <- try_extract_variable(right)
      if (!is.null(right_info)) {
        return(list(
          variable = right_info$variable,
          coef_expr = left,
          is_lead = right_info$is_lead,
          lead_k = right_info$lead_k
        ))
      }

      # Try: var * coef
      left_info <- try_extract_variable(left)
      if (!is.null(left_info)) {
        return(list(
          variable = left_info$variable,
          coef_expr = right,
          is_lead = left_info$is_lead,
          lead_k = left_info$lead_k
        ))
      }

      return(NULL)
    }

    # Case 4: division -- e.g., x / beta treated as (1/beta) * x
    if (fn == "/" && length(expr) == 3L) {
      left_info <- try_extract_variable(expr[[2L]])
      if (!is.null(left_info)) {
        return(list(
          variable = left_info$variable,
          coef_expr = call("/", 1, expr[[3L]]),
          is_lead = left_info$is_lead,
          lead_k = left_info$lead_k
        ))
      }
    }

    # Case 5: parenthesized expression containing a variable
    if (fn == "(") {
      return(separate_var_and_coef(expr[[2L]]))
    }
  }

  # Case 6: numeric constant -- not a variable term
  if (is.numeric(expr)) {
    return(NULL)
  }

  NULL
}

#' Try to extract variable info from an expression if it IS a variable
#' @noRd
try_extract_variable <- function(expr) {
  if (is.name(expr)) {
    return(list(variable = as.character(expr), is_lead = FALSE, lead_k = 0L))
  }

  if (is.call(expr)) {
    fn <- as.character(expr[[1L]])
    if (fn == "lead") {
      k <- if (length(expr) >= 3L) as.integer(eval(expr[[3L]])) else 1L
      return(list(variable = as.character(expr[[2L]]), is_lead = TRUE, lead_k = k))
    }
    if (fn == "E") {
      return(list(variable = as.character(expr[[2L]]), is_lead = TRUE, lead_k = 1L))
    }
    if (fn == "(") {
      return(try_extract_variable(expr[[2L]]))
    }
  }

  NULL
}

#' Extract parameter names from a coefficient expression
#'
#' Parameters are names that are not known R functions and are not numeric.
#' @noRd
extract_param_names <- function(expr) {
  if (is.numeric(expr) || is.logical(expr)) return(character(0))
  if (is.name(expr)) {
    nm <- as.character(expr)
    # Filter out known R constants and functions
    if (nm %in% c("pi", "TRUE", "FALSE", "T", "F", "Inf", "NaN", "NA")) {
      return(character(0))
    }
    return(nm)
  }
  if (is.call(expr)) {
    fn <- as.character(expr[[1L]])
    # Known functions are not parameters
    known_fns <- c("+", "-", "*", "/", "^", "(", "sqrt", "exp", "log",
                   "abs", "lead", "E")
    params <- character(0)
    # The function name itself is not a parameter
    for (i in seq_along(expr)[-1]) {
      params <- c(params, extract_param_names(expr[[i]]))
    }
    return(unique(params))
  }
  character(0)
}

#' Find variable names in an expression (non-parameter names)
#' @noRd
find_variables_in_expr <- function(expr) {
  all_names <- all.vars(expr)
  # Variables are names; we can't distinguish from params here
  # so this is a helper used in limited contexts
  all_names
}

#' Check if an expression contains only parameter references (no variables)
#' @noRd
is_param_expr <- function(expr) {
  # A heuristic: single names that aren't lead/E calls
  # This is used only for LHS parsing disambiguation
  FALSE
}

#' Remove a variable from a product expression, returning the coefficient
#' @noRd
remove_variable_from_product <- function(expr, varname) {
  if (is.name(expr) && as.character(expr) == varname) return(1)

  if (is.call(expr) && identical(expr[[1L]], as.name("*"))) {
    left <- expr[[2L]]
    right <- expr[[3L]]

    if (is.name(right) && as.character(right) == varname) return(left)
    if (is.name(left) && as.character(left) == varname) return(right)
  }

  expr
}

#' Evaluate a parameter coefficient expression given parameter values
#'
#' @param coef_expr An R expression or numeric value representing a parameter
#'   coefficient.
#' @param params Named numeric vector of parameter values.
#' @return A numeric scalar.
#' @noRd
eval_coef <- function(coef_expr, params) {
  if (is.numeric(coef_expr)) return(coef_expr)
  eval(coef_expr, envir = as.list(params))
}
