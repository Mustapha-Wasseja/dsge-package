# ==========================================================================
# Variance decomposition (unconditional and FEVD)
# ==========================================================================

#' Variance Decomposition
#'
#' Computes the share of each observable variable's variance attributable
#' to each structural shock.
#'
#' Two flavours are supported, controlled by the \code{horizon} argument:
#' \itemize{
#'   \item \strong{Unconditional} (default, \code{horizon = NULL}): the
#'     decomposition of the long-run / steady-state variance.  For each
#'     shock \eqn{j}, the state covariance contribution \eqn{\Sigma_x^{(j)}}
#'     solves the discrete Lyapunov equation
#'     \eqn{\Sigma_x^{(j)} = H \Sigma_x^{(j)} H' + M_j \sigma_j^2 M_j'} and
#'     the observable variance share is
#'     \eqn{\mathrm{diag}(G \Sigma_x^{(j)} G')}.
#'   \item \strong{Forecast-error variance decomposition (FEVD)}
#'     (\code{horizon = 1:H}): the share of each shock in the
#'     \eqn{h}-step-ahead forecast-error variance for a vector of
#'     horizons.  Forecast-error variance at horizon \eqn{h} is
#'     \eqn{\sum_{k=0}^{h-1} H^k M \Sigma_\varepsilon M' (H^k)'}.
#' }
#'
#' @param x A \code{dsge_solution}, \code{dsge_fit}, or \code{dsge_bayes}
#'   object.
#' @param horizon \code{NULL} (default) for the unconditional
#'   decomposition, or an integer vector of horizons for FEVD.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"dsge_variance_decomposition"}
#'   containing:
#'   \describe{
#'     \item{contribution}{Either a \eqn{n_o \times n_e} matrix (unconditional)
#'       or a \eqn{n_h \times n_o \times n_e} array (FEVD) of variance
#'       contributions in level units (squared standard deviations).}
#'     \item{contribution_pct}{Same shape as \code{contribution}, but
#'       normalised so that each variable's shares across shocks sum to
#'       100 percent.}
#'     \item{obs_names}{Character vector of observable variable names.}
#'     \item{shock_names}{Character vector of structural shock names.}
#'     \item{horizon}{The horizon argument (or \code{Inf} for the
#'       unconditional case).}
#'     \item{type}{\code{"unconditional"} or \code{"fevd"}.}
#'   }
#'
#' @examples
#' \donttest{
#' nk <- dsge_model(
#'   obs(p   ~ beta * lead(p) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(p) - g)),
#'   obs(r   ~ psi * p + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9)
#' )
#' sol <- solve_dsge(nk,
#'   params   = c(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9),
#'   shock_sd = c(e.u = 1.0, e.g = 0.5))
#' vd  <- variance_decomposition(sol)
#' print(vd)
#' plot(vd)
#'
#' fevd <- variance_decomposition(sol, horizon = c(1, 4, 8, 20))
#' plot(fevd)
#' }
#'
#' @export
variance_decomposition <- function(x, horizon = NULL, ...) {
  UseMethod("variance_decomposition")
}

#' @rdname variance_decomposition
#' @export
variance_decomposition.dsge_solution <- function(x, horizon = NULL, ...) {
  .variance_decomposition_impl(x, horizon = horizon)
}

#' @rdname variance_decomposition
#' @export
variance_decomposition.dsge_fit <- function(x, horizon = NULL, ...) {
  .variance_decomposition_impl(x$solution, horizon = horizon)
}

#' @rdname variance_decomposition
#' @export
variance_decomposition.dsge_bayes <- function(x, horizon = NULL, ...) {
  # Use the underlying solution stored on first model solve at posterior mode/mean.
  # Fall back to first stored solution if available.
  sol <- x$solution
  if (is.null(sol)) {
    stop("variance_decomposition.dsge_bayes requires the dsge_bayes object ",
         "to carry a $solution component.", call. = FALSE)
  }
  .variance_decomposition_impl(sol, horizon = horizon)
}


#' @noRd
.variance_decomposition_impl <- function(sol, horizon = NULL) {
  G <- sol$G
  H <- sol$H
  M <- sol$M
  D <- sol$D

  shock_sd <- sol$shock_sd
  if (is.null(shock_sd)) {
    stop("Solution does not carry shock_sd; cannot decompose variance.",
         call. = FALSE)
  }
  shock_names <- names(shock_sd)
  if (is.null(shock_names)) {
    shock_names <- paste0("e", seq_along(shock_sd))
  }

  Z <- D %*% G
  obs_names <- rownames(D)
  if (is.null(obs_names)) obs_names <- paste0("y", seq_len(nrow(D)))

  n_e <- ncol(M)
  n_o <- nrow(Z)

  if (length(shock_sd) != n_e) {
    stop("shock_sd length (", length(shock_sd),
         ") does not match number of shock columns in M (", n_e, ").",
         call. = FALSE)
  }

  if (is.null(horizon)) {
    # Unconditional decomposition via per-shock Lyapunov solves
    contribution <- matrix(0, n_o, n_e,
                           dimnames = list(obs_names, shock_names))
    # NOTE: in solve_dsge, M is built as B0_B1G^{-1} %*% C %*% diag(shock_sd),
    # so the shock standard deviations are already absorbed into M.  No
    # further sigma scaling is needed when forming the per-shock Q.
    for (j in seq_len(n_e)) {
      M_j <- M[, j, drop = FALSE]
      Q_j <- M_j %*% t(M_j)
      Sigma_x_j <- compute_unconditional_P(H, Q_j)
      Sigma_y_j <- Z %*% Sigma_x_j %*% t(Z)
      contribution[, j] <- pmax(diag(Sigma_y_j), 0)
    }
    contribution_pct <- .row_normalize_pct(contribution)
    return(structure(
      list(
        contribution     = contribution,
        contribution_pct = contribution_pct,
        obs_names        = obs_names,
        shock_names      = shock_names,
        horizon          = Inf,
        type             = "unconditional"
      ),
      class = "dsge_variance_decomposition"
    ))
  }

  # FEVD path
  horizons <- as.integer(horizon)
  if (any(horizons < 1L)) {
    stop("All horizons must be >= 1.", call. = FALSE)
  }
  n_h <- length(horizons)
  max_h <- max(horizons)

  fevd <- array(0, dim = c(n_h, n_o, n_e),
                dimnames = list(as.character(horizons), obs_names, shock_names))

  # Pre-compute H^(k-1) M for k = 1..max_h.  M already absorbs shock_sd,
  # so no further sigma scaling is required here.
  HkM <- array(0, dim = c(max_h, nrow(H), n_e))
  HkM[1, , ] <- M
  if (max_h > 1L) {
    for (k in 2:max_h) {
      HkM[k, , ] <- H %*% HkM[k - 1, , , drop = TRUE]
    }
  }

  for (h_idx in seq_along(horizons)) {
    h <- horizons[h_idx]
    for (j in seq_len(n_e)) {
      Sigma_x_j <- matrix(0, nrow(H), nrow(H))
      for (k in seq_len(h)) {
        v <- HkM[k, , j, drop = TRUE]
        Sigma_x_j <- Sigma_x_j + tcrossprod(v)
      }
      Sigma_y_j <- Z %*% Sigma_x_j %*% t(Z)
      fevd[h_idx, , j] <- pmax(diag(Sigma_y_j), 0)
    }
  }

  fevd_pct <- fevd
  for (h_idx in seq_len(n_h)) {
    fevd_pct[h_idx, , ] <- .row_normalize_pct(fevd[h_idx, , ])
  }

  structure(
    list(
      contribution     = fevd,
      contribution_pct = fevd_pct,
      obs_names        = obs_names,
      shock_names      = shock_names,
      horizon          = horizons,
      type             = "fevd"
    ),
    class = "dsge_variance_decomposition"
  )
}


#' Normalise rows of a non-negative matrix to sum to 100 (percent)
#' @noRd
.row_normalize_pct <- function(M) {
  totals <- rowSums(M)
  out <- M
  for (i in seq_len(nrow(M))) {
    if (totals[i] > 0) out[i, ] <- M[i, ] / totals[i] * 100
  }
  out
}


# --------------------------------------------------------------------------
# Print / plot methods
# --------------------------------------------------------------------------

#' @export
print.dsge_variance_decomposition <- function(x, digits = 1, ...) {
  if (x$type == "unconditional") {
    cat("Unconditional Variance Decomposition (% of total variance)\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")
    print(round(x$contribution_pct, digits))
  } else {
    cat("Forecast Error Variance Decomposition (% of FEV)\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")
    for (i in seq_along(x$obs_names)) {
      cat("\nVariable: ", x$obs_names[i], "\n", sep = "")
      mat_i <- x$contribution_pct[, i, , drop = FALSE]
      dim(mat_i) <- dim(mat_i)[c(1, 3)]
      dimnames(mat_i) <- list(paste0("h=", x$horizon), x$shock_names)
      print(round(mat_i, digits))
    }
  }
  invisible(x)
}


#' Plot a Variance Decomposition
#'
#' For unconditional decompositions, draws a horizontal stacked bar chart
#' showing shock shares (in percent) for each observable.  For FEVD, draws
#' one stacked-bar panel per observable with horizons on the x-axis.
#'
#' @param x A \code{dsge_variance_decomposition} object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}; called for the side effect of
#'   producing a plot.
#'
#' @export
plot.dsge_variance_decomposition <- function(x, ...) {
  if (x$type == "unconditional") {
    .plot_vd_unconditional(x)
  } else {
    .plot_vd_fevd(x)
  }
  invisible(x)
}


#' @noRd
.plot_vd_unconditional <- function(x) {
  pct <- x$contribution_pct       # n_o x n_e
  n_e <- ncol(pct)
  cols <- .dsge_palette(n_e)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_single()
  graphics::par(mar = c(4.0, 5.0, 2.4, 8.5), xpd = TRUE)

  # barplot expects rows = stack categories, columns = bars
  graphics::barplot(t(pct), beside = FALSE, horiz = TRUE,
                    col = cols, border = "white",
                    xlim = c(0, 100), las = 1,
                    main = "Unconditional Variance Decomposition",
                    xlab = "Share (%)")
  .dsge_grid(horizontal = FALSE, vertical = TRUE)
  .dsge_legend("topright", inset = c(-0.22, 0),
               legend = x$shock_names, fill = cols, title = "Shocks")
}


#' @noRd
.plot_vd_fevd <- function(x) {
  n_h <- length(x$horizon)
  n_o <- length(x$obs_names)
  n_e <- length(x$shock_names)
  cols <- .dsge_palette(n_e)

  nc <- min(2L, n_o)
  nr <- ceiling(n_o / nc)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(nr, nc)
  graphics::par(mar = c(3.6, 3.8, 2.0, 7.5), xpd = TRUE)

  for (i in seq_len(n_o)) {
    mat_i <- x$contribution_pct[, i, , drop = FALSE]
    dim(mat_i) <- dim(mat_i)[c(1, 3)]   # n_h x n_e
    graphics::barplot(t(mat_i), beside = FALSE,
                      col = cols, border = "white",
                      names.arg = x$horizon,
                      las = 1, ylim = c(0, 100),
                      main = x$obs_names[i],
                      xlab = "Horizon", ylab = "Share (%)")
    if (i == 1L) {
      .dsge_legend("topright", inset = c(-0.32, 0),
                   legend = x$shock_names, fill = cols, title = "Shocks")
    }
  }
}
