#' Forward Lead Operator for DSGE Equations
#'
#' Marks a variable as a one-period-ahead model-consistent expectation
#' in a DSGE equation formula. This is the core primitive for forward-looking
#' variables.
#'
#' @param x A variable name (unquoted) within a DSGE equation formula.
#' @param k Integer lead horizon. Currently only `k = 1` is supported.
#'   Default is 1.
#'
#' @details
#' `lead(x)` in a DSGE equation represents the expectation of variable `x`
#' one period ahead, conditional on the model. In the literature, this
#' corresponds to \eqn{E_t[x_{t+1}]}.
#'
#' This function is not meant to be called directly. It is recognized by
#' the equation parser inside [dsge_model()].
#'
#' @return This function is not meant to be called directly; it always
#'   throws an error. It is recognized as a syntactic marker by the
#'   equation parser inside [dsge_model()].
#' @seealso [E()] for a user-friendly alias, [dsge_model()]
#' @export
lead <- function(x, k = 1L) {
  stop("`lead()` should only be used inside dsge_model() equation formulas.",
       call. = FALSE)
}

#' Expectation Operator (Alias for lead)
#'
#' A user-friendly alias for `lead(x, 1)`. Represents the one-period-ahead
#' model-consistent expectation of variable `x`.
#'
#' @param x A variable name (unquoted) within a DSGE equation formula.
#'
#' @details
#' `E(x)` is equivalent to `lead(x, 1)`. The parser translates `E(x)` to
#' `lead(x, 1)` internally.
#'
#' @return This function is not meant to be called directly; it always
#'   throws an error. It is recognized as a syntactic marker by the
#'   equation parser inside [dsge_model()].
#' @seealso [lead()], [dsge_model()]
#' @export
E <- function(x) {
  stop("`E()` should only be used inside dsge_model() equation formulas.",
       call. = FALSE)
}

#' Define an Observed Control Variable Equation
#'
#' Wraps a formula to mark it as an equation for an observed control variable
#' in a linear DSGE model.
#'
#' @param formula A formula of the form `variable ~ expression`.
#'
#' @return A list with class `"dsge_equation"` containing the parsed equation
#'   and its type.
#'
#' @seealso [unobs()], [state()], [dsge_model()]
#' @export
obs <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`obs()` requires a formula (e.g., obs(y ~ beta * lead(x))).",
         call. = FALSE)
  }
  structure(list(formula = formula, type = "observed"),
            class = "dsge_equation")
}

#' Define an Unobserved Control Variable Equation
#'
#' Wraps a formula to mark it as an equation for an unobserved control variable
#' in a linear DSGE model.
#'
#' @param formula A formula of the form `variable ~ expression`.
#'
#' @return A list with class `"dsge_equation"` containing the parsed equation
#'   and its type.
#'
#' @seealso [obs()], [state()], [dsge_model()]
#' @export
unobs <- function(formula) {

  if (!inherits(formula, "formula")) {
    stop("`unobs()` requires a formula (e.g., unobs(x ~ lead(x) - r)).",
         call. = FALSE)
  }
  structure(list(formula = formula, type = "unobserved"),
            class = "dsge_equation")
}

#' Define a State Variable Equation
#'
#' Wraps a formula to mark it as an equation for a state variable
#' in a linear DSGE model. State equations describe the evolution of
#' state variables one period ahead.
#'
#' @param formula A formula of the form `variable ~ expression`. The
#'   left-hand side variable represents the one-period lead of the state.
#' @param shock Logical. If `TRUE` (default), an unobserved shock is
#'   attached to this state equation (exogenous state). If `FALSE`, no
#'   shock is attached (endogenous state / deterministic state).
#'
#' @return A list with class `"dsge_equation"` containing the parsed equation
#'   and its type.
#'
#' @seealso [obs()], [unobs()], [dsge_model()]
#' @export
state <- function(formula, shock = TRUE) {
  if (!inherits(formula, "formula")) {
    stop("`state()` requires a formula (e.g., state(u ~ rhou * u)).",
         call. = FALSE)
  }
  structure(list(formula = formula, type = "state", shock = shock),
            class = "dsge_equation")
}

#' Define a Linear DSGE Model
#'
#' Constructs a linear DSGE model object from a set of equations.
#' Each equation is wrapped in [obs()], [unobs()], or [state()] to
#' indicate the role of the left-hand-side variable.
#'
#' @param ... Equation specifications created by [obs()], [unobs()],
#'   and [state()].
#' @param fixed Named list of parameter values to hold fixed
#'   (constrained) during estimation.
#' @param start Named list of starting values for free parameters.
#'
#' @return An object of class `"dsge_model"` containing:
#'   \describe{
#'     \item{equations}{List of parsed equation objects.}
#'     \item{variables}{List with elements `observed`, `unobserved`,
#'       `exo_state`, `endo_state`.}
#'     \item{parameters}{Character vector of all parameter names.}
#'     \item{free_parameters}{Character vector of free (estimable) parameter names.}
#'     \item{fixed}{Named list of fixed parameter values.}
#'     \item{start}{Named list of starting values.}
#'   }
#'
#' @details
#' A linear DSGE model is specified as a system of equations in which
#' variables enter linearly but parameters may enter nonlinearly.
#'
#' Use `lead(x)` or `E(x)` within formulas to denote the one-period-ahead
#' model-consistent expectation of variable `x`.
#'
#' The number of exogenous state variables (those with shocks) must equal
#' the number of observed control variables.
#'
#' @examples
#' # Simple New Keynesian model
#' nk <- dsge_model(
#'   obs(p   ~ beta * lead(p) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(p) - g)),
#'   obs(r   ~ psi * p + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.96),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
#' )
#'
#' @export
dsge_model <- function(..., fixed = list(), start = list()) {
  eqs <- list(...)

  # Validate that all arguments are dsge_equation objects

  for (i in seq_along(eqs)) {
    if (!inherits(eqs[[i]], "dsge_equation")) {
      stop("Argument ", i, " is not a DSGE equation. ",
           "Use obs(), unobs(), or state() to wrap each equation.",
           call. = FALSE)
    }
  }

  if (length(eqs) == 0L) {
    stop("At least one equation must be provided.", call. = FALSE)
  }

  # Parse each equation
  parsed <- lapply(eqs, parse_equation)

  # Extract variable roles
  observed <- character(0)
  unobserved <- character(0)
  exo_state <- character(0)
  endo_state <- character(0)

  for (p in parsed) {
    if (p$type == "observed") {
      observed <- c(observed, p$lhs_var)
    } else if (p$type == "unobserved") {
      unobserved <- c(unobserved, p$lhs_var)
    } else if (p$type == "state") {
      if (p$shock) {
        exo_state <- c(exo_state, p$lhs_var)
      } else {
        endo_state <- c(endo_state, p$lhs_var)
      }
    }
  }

  all_vars <- c(observed, unobserved, exo_state, endo_state)

  # Check for duplicate LHS variables
  if (anyDuplicated(all_vars)) {
    dups <- all_vars[duplicated(all_vars)]
    stop("Variable(s) appear on the LHS of more than one equation: ",
         paste(dups, collapse = ", "), call. = FALSE)
  }

  # Check that number of observed controls equals number of shocks
  n_shocks <- length(exo_state)
  n_obs <- length(observed)
  if (n_shocks != n_obs) {
    stop("Number of exogenous state variables (", n_shocks,
         ") must equal number of observed control variables (", n_obs, ").",
         call. = FALSE)
  }

  # Collect all parameter names from all equations
  all_params <- unique(unlist(lapply(parsed, function(p) p$parameters)))

  # Validate fixed parameters exist in the model

  if (length(fixed) > 0) {
    unknown_fixed <- setdiff(names(fixed), all_params)
    if (length(unknown_fixed) > 0) {
      stop("Fixed parameter(s) not found in model equations: ",
           paste(unknown_fixed, collapse = ", "), call. = FALSE)
    }
  }

  free_params <- setdiff(all_params, names(fixed))

  # Validate start parameters
  if (length(start) > 0) {
    unknown_start <- setdiff(names(start), free_params)
    if (length(unknown_start) > 0) {
      # Could be fixed params specified in start -- warn but don't error
      unknown_start <- setdiff(unknown_start, names(fixed))
      if (length(unknown_start) > 0) {
        stop("Start value parameter(s) not found in model: ",
             paste(unknown_start, collapse = ", "), call. = FALSE)
      }
    }
  }

  # Validate that all RHS variables are declared somewhere
  rhs_vars <- unique(unlist(lapply(parsed, function(p) p$rhs_variables)))
  undeclared <- setdiff(rhs_vars, all_vars)
  if (length(undeclared) > 0) {
    stop("Variable(s) used in equations but not declared as control or state: ",
         paste(undeclared, collapse = ", "), call. = FALSE)
  }

  structure(
    list(
      equations = parsed,
      variables = list(
        observed = observed,
        unobserved = unobserved,
        exo_state = exo_state,
        endo_state = endo_state
      ),
      all_variables = all_vars,
      parameters = all_params,
      free_parameters = free_params,
      fixed = fixed,
      start = start,
      n_obs_controls = n_obs,
      n_unobs_controls = length(unobserved),
      n_exo_states = length(exo_state),
      n_endo_states = length(endo_state),
      n_controls = n_obs + length(unobserved),
      n_states = length(exo_state) + length(endo_state),
      call = match.call()
    ),
    class = "dsge_model"
  )
}
