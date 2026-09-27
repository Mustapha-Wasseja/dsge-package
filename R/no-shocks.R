# Error for tools that need stochastic shocks, called on a solution of a
# deterministic model (one without varexo, e.g. a perfect-foresight model).
.require_shocks <- function(sol, what) {
  if (is.null(sol$M) || ncol(sol$M) == 0L) {
    stop(what, " requires a model with stochastic shocks, but this model ",
         "has none. For deterministic transitions use ",
         "simulate_perfect_foresight().", call. = FALSE)
  }
  invisible(TRUE)
}
