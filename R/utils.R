# Internal utility functions

#' Check if a dsge_fit or dsge_solution is stable
#' @noRd
check_stable <- function(x) {
  sol <- if (inherits(x, "dsge_fit")) x$solution else x
  if (!sol$stable) {
    stop("Model is not saddle-path stable.", call. = FALSE)
  }
  invisible(TRUE)
}
