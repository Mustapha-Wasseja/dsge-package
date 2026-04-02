# Nonlinear DSGE model specification
#
# Constructor and print method for nonlinear DSGE models defined via
# string-based equations.

#' Define a Nonlinear DSGE Model
#'
#' Constructs a nonlinear DSGE model from string-based equations.
#' Each equation is a character string of the form `"LHS = RHS"`.
#' State equations are detected automatically when the left-hand side
#' is of the form `VAR(+1)`.
#'
#' @param ... Character strings, each defining one model equation.
#'   Use `VAR(+1)` to denote the one-period-ahead value of variable `VAR`.
#' @param observed Character vector of observed control variable names.
#' @param unobserved Character vector of unobserved control variable names.
#'   Default is `character(0)`.
#' @param exo_state Character vector of exogenous state variable names.
#'   These have shocks attached. The number of exogenous states must equal
#'   the number of observed controls.
#' @param endo_state Character vector of endogenous (predetermined) state
#'   variable names. These have no shocks. Default is `character(0)`.
#' @param fixed Named list of parameter values to hold fixed during estimation.
#' @param start Named list of starting values for free parameters.
#' @param ss_guess Named numeric vector of initial guesses for steady-state
#'   solving. If `NULL`, defaults to 1 for all variables.
#' @param ss_function Optional function that computes the steady state
#'   analytically. Must accept a named parameter vector and return a named
#'   numeric vector of steady-state variable values.
#'
#' @return An object of class `"dsgenl_model"`.
#'
#' @details
#' Equations are written in standard mathematical notation with `=` as
#' the equality sign and `VAR(+1)` for leads. For example:
#'
#' ```
#' "1/C = beta / C(+1) * (alpha * K^(alpha-1) + 1 - delta)"
#' ```
#'
#' State equations must have exactly one lead variable on the left-hand side
#' (e.g., `"K(+1) = K^alpha - C + (1 - delta) * K"`). Control equations
#' are all remaining equations.
#'
#' Control equations must be provided in the order matching the controls
#' vector (observed first, then unobserved).
#'
#' @examples
#' # Simple RBC model
#' rbc <- dsgenl_model(
#'   "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
#'   "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
#'   "Z(+1) = rho * Z",
#'   observed = "C",
#'   endo_state = "K",
#'   exo_state = "Z",
#'   fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
#'   start = list(rho = 0.9)
#' )
#'
#' @export
dsgenl_model <- function(..., observed = character(0),
                         unobserved = character(0),
                         exo_state, endo_state = character(0),
                         fixed = list(), start = list(),
                         ss_guess = NULL, ss_function = NULL) {

  eq_list <- list(...)
  for (i in seq_along(eq_list)) {
    if (!is.character(eq_list[[i]]) || length(eq_list[[i]]) != 1L) {
      stop("Argument ", i, " must be a single equation string.", call. = FALSE)
    }
  }
  eq_strings <- unlist(eq_list)

  if (length(eq_strings) == 0L) {
    stop("At least one equation must be provided.", call. = FALSE)
  }

  # Parse all equations
  parsed <- lapply(eq_strings, parse_nl_equation)

  # Classify equations
  is_state <- vapply(parsed, function(p) p$is_state_eq, logical(1))
  state_eq_idx <- which(is_state)
  control_eq_idx <- which(!is_state)

  # Validate variable declarations
  controls <- c(observed, unobserved)
  states <- c(exo_state, endo_state)
  all_vars <- c(controls, states)

  if (length(all_vars) == 0L) {
    stop("At least one variable must be declared.", call. = FALSE)
  }

  if (anyDuplicated(all_vars)) {
    stop("Duplicate variable declarations: ",
         paste(all_vars[duplicated(all_vars)], collapse = ", "), call. = FALSE)
  }

  n_c <- length(controls)
  n_s <- length(states)
  n_obs <- length(observed)
  n_exo <- length(exo_state)

  # Check equation counts
  if (length(control_eq_idx) != n_c) {
    stop("Expected ", n_c, " control equation(s) but found ",
         length(control_eq_idx), ".", call. = FALSE)
  }
  if (length(state_eq_idx) != n_s) {
    stop("Expected ", n_s, " state equation(s) but found ",
         length(state_eq_idx), ".", call. = FALSE)
  }

  # Check observability condition
  if (n_obs != n_exo) {
    stop("Number of observed controls (", n_obs,
         ") must equal number of exogenous states (", n_exo, ").",
         call. = FALSE)
  }

  # Order state equations to match the states vector
  state_parsed <- parsed[state_eq_idx]
  state_lhs_vars <- vapply(state_parsed,
                           function(p) p$lhs_state_var, character(1))

  undeclared_state <- setdiff(state_lhs_vars, states)
  if (length(undeclared_state) > 0L) {
    stop("State equation LHS variable(s) not declared as state: ",
         paste(undeclared_state, collapse = ", "), call. = FALSE)
  }

  state_order <- match(states, state_lhs_vars)
  if (any(is.na(state_order))) {
    missing <- states[is.na(state_order)]
    stop("No state equation found for: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  state_parsed <- state_parsed[state_order]

  # Control equations in user-provided order
  control_parsed <- parsed[control_eq_idx]

  # Ordered: control first, then state
  ordered_parsed <- c(control_parsed, state_parsed)

  # Identify parameters: all names in equations that are not variables
  all_current_refs <- unique(unlist(lapply(ordered_parsed,
                                          function(p) p$current_refs)))
  all_params <- setdiff(all_current_refs, all_vars)

  # Validate lead variables are declared
  all_lead_vars <- unique(unlist(lapply(ordered_parsed,
                                       function(p) p$lead_var_names)))
  undeclared_leads <- setdiff(all_lead_vars, all_vars)
  if (length(undeclared_leads) > 0L) {
    stop("Lead variable(s) not declared: ",
         paste(undeclared_leads, collapse = ", "), call. = FALSE)
  }

  # Validate fixed parameters
  if (length(fixed) > 0L) {
    unknown_fixed <- setdiff(names(fixed), all_params)
    if (length(unknown_fixed) > 0L) {
      stop("Fixed parameter(s) not found in equations: ",
           paste(unknown_fixed, collapse = ", "), call. = FALSE)
    }
  }

  free_params <- setdiff(all_params, names(fixed))

  # Validate start parameters
  if (length(start) > 0L) {
    unknown_start <- setdiff(names(start), free_params)
    unknown_start <- setdiff(unknown_start, names(fixed))
    if (length(unknown_start) > 0L) {
      stop("Start value parameter(s) not found in model: ",
           paste(unknown_start, collapse = ", "), call. = FALSE)
    }
  }

  # Build evaluation function
  eval_fn <- build_nl_eval_function(ordered_parsed)

  structure(
    list(
      equations = ordered_parsed,
      eq_strings = vapply(ordered_parsed, function(p) p$original, character(1)),
      variables = list(
        observed = observed,
        unobserved = unobserved,
        exo_state = exo_state,
        endo_state = endo_state
      ),
      all_variables = all_vars,
      controls = controls,
      states = states,
      parameters = all_params,
      free_parameters = free_params,
      fixed = fixed,
      start = start,
      n_obs_controls = n_obs,
      n_unobs_controls = length(unobserved),
      n_exo_states = n_exo,
      n_endo_states = length(endo_state),
      n_controls = n_c,
      n_states = n_s,
      eval_fn = eval_fn,
      ss_guess = ss_guess,
      ss_function = ss_function,
      call = match.call()
    ),
    class = "dsgenl_model"
  )
}

#' @export
print.dsgenl_model <- function(x, ...) {
  cat("Nonlinear DSGE Model\n")
  cat("  Equations:    ", length(x$equations), "\n")
  cat("  Controls:     ", x$n_controls,
      "(", x$n_obs_controls, "observed,",
      x$n_unobs_controls, "unobserved )\n")
  cat("  States:       ", x$n_states,
      "(", x$n_exo_states, "exogenous,",
      x$n_endo_states, "endogenous )\n")
  cat("  Parameters:   ", length(x$parameters),
      "(", length(x$free_parameters), "free,",
      length(x$fixed), "fixed )\n")
  cat("\n")
  for (i in seq_along(x$eq_strings)) {
    cat("  [", i, "] ", x$eq_strings[i], "\n", sep = "")
  }
  cat("\n")
  if (length(x$fixed) > 0L) {
    cat("  Fixed: ", paste(names(x$fixed), "=",
        unlist(x$fixed), collapse = ", "), "\n")
  }
  invisible(x)
}
