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
#' @return list(hess = data frame (eq, i, j, val) of the nonzero second
#'   derivatives, both (i, j) and (j, i); third = data frame
#'   (eq, i, j, l, val) of the nonzero third derivatives with all index
#'   permutations, or NULL), or NULL when symbolic differentiation fails.
#' @noRd
.symbolic_equation_derivs <- function(model, timed_names, point, params,
                                      order = 2L) {
  code <- .symbolic_higher_code(model, timed_names, order)
  if (is.null(code)) return(NULL)
  vals <- c(point, params)
  vals <- vals[!duplicated(names(vals))]
  if (is.function(model$resolve_values)) vals <- model$resolve_values(vals)
  env <- list2env(as.list(vals), parent = baseenv())
  h <- tryCatch(as.numeric(eval(code$h_code, env)), error = function(e) NULL)
  if (is.null(h) || length(h) != nrow(code$h_idx) || !all(is.finite(h))) {
    return(NULL)
  }
  third <- NULL
  if (order >= 3L) {
    t3 <- tryCatch(as.numeric(eval(code$t_code, env)),
                   error = function(e) NULL)
    if (is.null(t3) || length(t3) != nrow(code$t_idx)) return(NULL)
    third <- .expand_third(code$t_idx, t3)
  }
  list(hess = .expand_hess(code$h_idx, h), third = third)
}

#' Sparse symmetric Hessian entries: (eq, i, j, val) with both (i, j) and
#' (j, i) for i != j, from entries with i <= j
#' @noRd
.expand_hess <- function(idx, val) {
  keep <- val != 0
  idx <- idx[keep, , drop = FALSE]
  val <- val[keep]
  off <- idx[, 2L] != idx[, 3L]
  data.frame(eq = c(idx[, 1L], idx[off, 1L]),
             i = c(idx[, 2L], idx[off, 3L]),
             j = c(idx[, 3L], idx[off, 2L]),
             val = c(val, val[off]))
}

#' Sparse third-derivative entries (eq, i, j, l, val) with every distinct
#' permutation of (i, j, l), from entries with i <= j <= l
#' @noRd
.expand_third <- function(idx, val) {
  keep <- val != 0
  idx <- idx[keep, , drop = FALSE]
  val <- val[keep]
  if (!length(val)) {
    return(data.frame(eq = integer(0), i = integer(0), j = integer(0),
                      l = integer(0), val = numeric(0)))
  }
  perms <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2),
                c(3, 2, 1))
  p <- idx[, 2:4, drop = FALSE]
  all <- do.call(rbind, lapply(perms, function(q) {
    cbind(r = seq_len(nrow(p)), p[, q, drop = FALSE])
  }))
  all <- unique(all)
  data.frame(eq = idx[all[, 1L], 1L], i = all[, 2L], j = all[, 3L],
             l = all[, 4L], val = val[all[, 1L]])
}

#' Symbolic Jacobian of the model equations
#'
#' Differentiates each equation with stats::D with respect to the timed
#' variables it contains, and collects the nonzero derivatives into one
#' expression `c(d_1, ..., d_m)` with their (row, column) positions. The result is cached in `model$.cache` (when the model has
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
    # not byte-compiled: compiling large derivative expressions takes far
    # longer (tens of seconds for big models) than evaluating them
    list(rows = rows, cols = cols,
         code = as.call(c(list(as.name("c")), terms)))
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

#' Symbolic second (and third) derivatives of the model equations
#'
#' Like `.symbolic_jacobian_code()`: every nonzero second derivative (i <= j)
#' and, for order 3, third derivative (i <= j <= l) of every equation is
#' collected into one expression, cached in `model$.cache`.
#'
#' @return list(h_idx, h_code, t_idx, t_code) with index matrices
#'   (equation, i, j[, l]) in timed-variable positions, or NULL if an
#'   equation cannot be differentiated symbolically.
#' @noRd
.symbolic_higher_code <- function(model, timed_names, order) {
  exprs <- lapply(model$equations, function(eq) eq$expression[[1L]])
  key <- list(exprs, timed_names, order)
  cache <- model$.cache
  if (is.environment(cache) && identical(cache$higher_key, key)) {
    return(cache$higher_code)
  }
  is_zero <- function(d) is.numeric(d) && length(d) == 1L && d == 0
  # one c(...) call per order (not byte-compiled; see .symbolic_jacobian_code)
  make_c <- function(terms) as.call(c(list(as.name("c")), terms))
  code <- tryCatch({
    h_idx <- list(); h_terms <- list()
    t_idx <- list(); t_terms <- list()
    for (k in seq_along(exprs)) {
      e <- exprs[[k]]
      v <- intersect(all.vars(e), timed_names)
      idx <- match(v, timed_names)
      nv <- length(v)
      if (nv == 0L) next
      d1 <- lapply(v, function(a) stats::D(e, a))
      for (i in seq_len(nv)) {
        if (is_zero(d1[[i]])) next
        for (j in i:nv) {
          d2 <- stats::D(d1[[i]], v[j])
          if (is_zero(d2)) next
          h_idx[[length(h_idx) + 1L]] <- c(k, idx[i], idx[j])
          h_terms[[length(h_terms) + 1L]] <- d2
          if (order >= 3L) {
            for (l in j:nv) {
              d3 <- stats::D(d2, v[l])
              if (is_zero(d3)) next
              t_idx[[length(t_idx) + 1L]] <- c(k, idx[i], idx[j], idx[l])
              t_terms[[length(t_terms) + 1L]] <- d3
            }
          }
        }
      }
    }
    list(h_idx = if (length(h_idx)) do.call(rbind, h_idx) else
           matrix(integer(0), 0L, 3L),
         h_code = make_c(h_terms),
         t_idx = if (length(t_idx)) do.call(rbind, t_idx) else
           matrix(integer(0), 0L, 4L),
         t_code = if (order >= 3L) make_c(t_terms) else NULL)
  }, error = function(e) NULL)
  if (is.environment(cache)) {
    cache$higher_key <- key
    cache$higher_code <- code
  }
  code
}
