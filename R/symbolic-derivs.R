# Exact second and third derivatives of the model equations
#
# The equations of a dsgenl_model are R expressions in the timed variables
# (x for current values, x__f for leads). Their derivatives at the steady
# state are computed symbolically with stats::D, only for the variables that
# appear in each equation. This is exact (unlike finite differences, whose
# error enters the risk corrections of higher-order solutions) and much
# faster for large models, since each equation involves few variables. If an
# equation uses a function stats::D cannot differentiate, NULL is returned
# and the caller falls back to finite differences.

#' Symbolic Hessians (and third derivatives) of the model equations
#'
#' @param model A dsgenl_model.
#' @param timed_names Names of the timed variables, in the order used for
#'   the derivative arrays.
#' @param point Values of the timed variables (named).
#' @param params Parameter values (named).
#' @param order 2 for Hessians only, 3 to add third derivatives.
#' @return list(hess = list of n x n matrices, third = list of matrices
#'   with columns i, j, l, value holding the nonzero third derivatives with
#'   all index permutations, or NULL), or NULL when symbolic
#'   differentiation fails.
#' @noRd
.symbolic_equation_derivs <- function(model, timed_names, point, params,
                                      order = 2L) {
  vals <- c(point, params)
  vals <- vals[!duplicated(names(vals))]
  if (is.function(model$resolve_values)) vals <- model$resolve_values(vals)
  env <- list2env(as.list(vals), parent = baseenv())
  n <- length(timed_names)
  n_eq <- length(model$equations)
  hess <- vector("list", n_eq)
  third <- if (order >= 3L) vector("list", n_eq) else NULL
  ok <- tryCatch({
    for (k in seq_len(n_eq)) {
      e <- model$equations[[k]]$expression[[1L]]
      v <- intersect(all.vars(e), timed_names)
      idx <- match(v, timed_names)
      H <- matrix(0, n, n)
      Tk <- list()
      nv <- length(v)
      if (nv > 0L) {
        d1 <- lapply(v, function(a) stats::D(e, a))
        for (i in seq_len(nv)) {
          for (j in i:nv) {
            d2 <- stats::D(d1[[i]], v[j])
            h <- as.numeric(eval(d2, env))
            H[idx[i], idx[j]] <- h
            H[idx[j], idx[i]] <- h
            if (order >= 3L) {
              for (l in j:nv) {
                t3 <- as.numeric(eval(stats::D(d2, v[l]), env))
                if (t3 != 0) {
                  p <- idx[c(i, j, l)]
                  perms <- unique(list(p[c(1, 2, 3)], p[c(1, 3, 2)],
                                       p[c(2, 1, 3)], p[c(2, 3, 1)],
                                       p[c(3, 1, 2)], p[c(3, 2, 1)]))
                  for (q in perms) Tk[[length(Tk) + 1L]] <- c(q, t3)
                }
              }
            }
          }
        }
      }
      if (!all(is.finite(H))) stop("non-finite derivative")
      hess[[k]] <- H
      if (order >= 3L) {
        third[[k]] <- if (length(Tk)) do.call(rbind, Tk) else
          matrix(numeric(0), 0L, 4L)
      }
    }
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  list(hess = hess, third = third)
}
