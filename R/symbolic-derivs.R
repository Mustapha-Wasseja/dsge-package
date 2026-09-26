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

#' Compiled symbolic Jacobian of the model equations
#'
#' Differentiates each equation with stats::D with respect to the timed
#' variables it contains, and collects the nonzero derivatives into one
#' byte-compiled expression `c(d_1, ..., d_m)` with their (row, column)
#' positions. The result is cached in `model$.cache` (when the model has
#' one) and rebuilt if the equations or variable names change.
#'
#' @return list(rows, cols, code), or NULL if an equation uses a function
#'   stats::D cannot differentiate.
#' @noRd
.symbolic_jacobian_code <- function(model, timed_names) {
  exprs <- lapply(model$equations, function(eq) eq$expression[[1L]])
  key <- list(exprs, timed_names)
  cache <- model$.cache
  if (is.environment(cache) && identical(cache$jac_key, key)) {
    return(cache$jac_code)
  }
  code <- tryCatch({
    rows <- integer(0)
    cols <- integer(0)
    terms <- list()
    for (k in seq_along(exprs)) {
      v <- intersect(all.vars(exprs[[k]]), timed_names)
      for (a in v) {
        d <- stats::D(exprs[[k]], a)
        if (is.numeric(d) && length(d) == 1L && d == 0) next
        rows <- c(rows, k)
        cols <- c(cols, match(a, timed_names))
        terms[[length(terms) + 1L]] <- d
      }
    }
    call <- as.call(c(list(as.name("c")), terms))
    list(rows = rows, cols = cols,
         code = compiler::compile(call, options = list(suppressAll = TRUE)))
  }, error = function(e) NULL)
  if (is.environment(cache)) {
    cache$jac_key <- key
    cache$jac_code <- code
  }
  code
}

#' Exact Jacobian of the model equations at a point
#'
#' @param model A dsgenl_model.
#' @param timed_names Names of the timed variables (columns of the result).
#' @param point Values of the timed variables (named, in that order).
#' @param params Parameter values (named).
#' @return An n_eq x length(timed_names) matrix, or NULL when symbolic
#'   differentiation is not possible (the caller then differentiates
#'   numerically).
#' @noRd
.symbolic_jacobian <- function(model, timed_names, point, params) {
  jc <- .symbolic_jacobian_code(model, timed_names)
  if (is.null(jc)) return(NULL)
  vals <- c(point, params)
  vals <- vals[!duplicated(names(vals))]
  if (is.function(model$resolve_values)) vals <- model$resolve_values(vals)
  env <- list2env(as.list(vals), parent = baseenv())
  d <- tryCatch(as.numeric(eval(jc$code, env)), error = function(e) NULL)
  if (is.null(d) || length(d) != length(jc$rows) || !all(is.finite(d))) {
    return(NULL)
  }
  J <- matrix(0, length(model$equations), length(timed_names))
  J[cbind(jc$rows, jc$cols)] <- d
  J
}
