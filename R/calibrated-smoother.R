# ==========================================================================
# Calibrated smoother -- Kalman smoothing without estimation
# ==========================================================================
#
# Provides Dynare's `calibrated_smoother` functionality: run the Kalman
# smoother on a model that has been solved (calibrated) but not estimated.
# All existing smoother infrastructure works on dsge_fit / dsge_bayes
# objects; this file adds methods for dsge_solution.
# ==========================================================================


#' Run the Kalman Smoother on a Calibrated Model
#'
#' Convenience wrapper that runs the Kalman smoother, smoothed-shock
#' recovery, and historical shock decomposition on a calibrated
#' \code{dsge_solution} object together with observed data.  Analogous to
#' Dynare's \code{calibrated_smoother} command.
#'
#' @param sol A \code{dsge_solution} object (from \code{\link{solve_dsge}}).
#' @param data Matrix or data frame of observed variables.  Column names
#'   must match \code{rownames(sol$D)} (the observable names).
#' @param what Character vector of outputs to compute.  Any subset of
#'   \code{c("states", "shocks", "decomposition")}.  Default is all three.
#'
#' @return A named list containing the requested smoothed-output objects.
#'   Elements are themselves objects of the same classes returned by
#'   \code{\link{smooth_states}}, \code{\link{smooth_shocks}}, and
#'   \code{\link{shock_decomposition}}, so all existing print and plot
#'   methods apply.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ z),
#'   state(z ~ rho * z),
#'   fixed = list(rho = 0.8))
#' sol <- solve_dsge(m, params = c(rho = 0.8), shock_sd = c(z = 1))
#' set.seed(1)
#' e <- rnorm(100)
#' z_path <- numeric(100)
#' for (i in 2:100) z_path[i] <- 0.8 * z_path[i - 1] + e[i]
#' out <- calibrated_smoother(sol, data = data.frame(y = z_path))
#' plot(out$states)
#' }
#'
#' @seealso \code{\link{smooth_states.dsge_solution}},
#'   \code{\link{smooth_shocks.dsge_solution}},
#'   \code{\link{shock_decomposition.dsge_solution}}.
#'
#' @export
calibrated_smoother <- function(sol, data,
                                what = c("states", "shocks", "decomposition")) {
  if (!inherits(sol, "dsge_solution"))
    stop("`sol` must be a dsge_solution object (from solve_dsge()).",
         call. = FALSE)
  what <- match.arg(what, several.ok = TRUE)
  out <- list()
  if ("states" %in% what)
    out$states <- smooth_states(sol, data = data)
  if ("shocks" %in% what)
    out$shocks <- smooth_shocks(sol, data = data)
  if ("decomposition" %in% what)
    out$decomposition <- shock_decomposition(sol, data = data)
  class(out) <- "dsge_calibrated_smoother"
  out
}

#' @export
print.dsge_calibrated_smoother <- function(x, ...) {
  cat("Calibrated DSGE smoother results\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Available components:", paste(names(x), collapse = ", "), "\n")
  for (nm in names(x)) {
    cat("\n[", nm, "]\n", sep = "")
    print(x[[nm]])
  }
  invisible(x)
}


# --------------------------------------------------------------------------
# S3 methods on dsge_solution
# --------------------------------------------------------------------------

#' @rdname smooth_states
#' @param data Matrix or data frame of observed data with columns matching
#'   the model's observable names.  Required when \code{x} is a
#'   \code{dsge_solution}; ignored for \code{dsge_fit} / \code{dsge_bayes}
#'   objects (which carry their own data).
#' @export
smooth_states.dsge_solution <- function(x, data, ...) {
  .smoother_inputs_from_solution(x, data, fn = "smoother")
}

#' @rdname smooth_shocks
#' @param data Matrix or data frame of observed data with columns matching
#'   the model's observable names.  Required when \code{x} is a
#'   \code{dsge_solution}.
#' @export
smooth_shocks.dsge_solution <- function(x, data, ...) {
  .smoother_inputs_from_solution(x, data, fn = "shocks")
}

#' @rdname shock_decomposition
#' @param data Matrix or data frame of observed data with columns matching
#'   the model's observable names.  Required when \code{x} is a
#'   \code{dsge_solution}.
#' @export
shock_decomposition.dsge_solution <- function(x, data, ...) {
  .smoother_inputs_from_solution(x, data, fn = "decomposition")
}


#' Wrap a (sol, data) pair in a temporary dsge_fit-like object so the
#' existing smoother machinery can be reused unchanged.
#' @noRd
.smoother_inputs_from_solution <- function(sol, data, fn) {
  # Validate data
  if (missing(data) || is.null(data))
    stop("`data` is required when running the smoother on a dsge_solution.",
         call. = FALSE)
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data))
    stop("`data` must be a matrix or data frame.", call. = FALSE)

  obs_names <- rownames(sol$D)
  if (is.null(obs_names) && !is.null(sol$model))
    obs_names <- sol$model$variables$observed

  if (!is.null(obs_names)) {
    if (is.null(colnames(data))) {
      if (ncol(data) != length(obs_names))
        stop("`data` has ", ncol(data),
             " columns but model has ", length(obs_names),
             " observables.", call. = FALSE)
      colnames(data) <- obs_names
    } else {
      bad <- setdiff(colnames(data), obs_names)
      if (length(bad) > 0)
        stop("Data column(s) not in model observables: ",
             paste(bad, collapse = ", "), call. = FALSE)
      data <- data[, obs_names, drop = FALSE]
    }
  }

  # Demean data (the smoother operates around steady state / zero)
  data_means <- colMeans(data, na.rm = TRUE)
  data_dem   <- sweep(data, 2, data_means, FUN = "-")

  # Build a pseudo-dsge_fit object carrying just what the impl needs
  pseudo_fit <- structure(
    list(model    = sol$model,
         solution = sol,
         data     = data_dem,
         data_means = data_means,
         kalman   = NULL),
    class = "dsge_fit")

  if (fn == "smoother")      return(smooth_states_impl(pseudo_fit))
  if (fn == "shocks")        return(smooth_shocks_impl(pseudo_fit))
  if (fn == "decomposition") return(shock_decomposition_impl(pseudo_fit))
  stop("Unknown fn: ", fn)
}
