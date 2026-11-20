# ==========================================================================
# Model-implied covariance reporting
# ==========================================================================

#' Model-implied covariance and correlation matrices
#'
#' Computes the unconditional (model-implied) covariance and correlation
#' matrices of observable variables from a solved or estimated DSGE model.
#' These are the theoretical second moments implied by the model at the
#' given parameter values.
#'
#' @param x A fitted model (\code{"dsge_fit"}, \code{"dsge_bayes"}) or
#'   a solved model (\code{"dsge_solution"}).
#' @param variables Character vector of variable names to include.
#'   Default \code{NULL} returns all observable variables.
#' @param n_lags Integer. If positive, also compute autocovariances
#'   at lags 1, ..., n_lags. Default 0 (contemporaneous only).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"dsge_covariance"} containing:
#' \describe{
#'   \item{covariance}{Covariance matrix of selected variables.}
#'   \item{correlation}{Correlation matrix of selected variables.}
#'   \item{std_dev}{Standard deviations (square root of diagonal).}
#'   \item{autocovariances}{List of lagged autocovariance matrices
#'     (empty if \code{n_lags = 0}).}
#'   \item{variables}{Variable names.}
#'   \item{n_lags}{Number of autocovariance lags computed.}
#' }
#'
#' @examples
#' mod <- dsge_model(
#'   obs(pi ~ beta * lead(pi) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(pi) - g)),
#'   obs(r ~ psi * pi + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.5, rhog = 0.5)
#' )
#' p <- list(kappa = 0.1, psi = 1.5, rhou = 0.5, rhog = 0.5)
#' s <- c(u = 0.5, g = 0.5)
#' sol <- solve_dsge(mod, params = p, shock_sd = s)
#' model_covariance(sol)
#'
#' @export
model_covariance <- function(x, variables = NULL, n_lags = 0L, ...) {
  UseMethod("model_covariance")
}

#' @export
model_covariance.dsge_solution <- function(x, variables = NULL,
                                           n_lags = 0L, ...) {
  sol <- x
  .compute_model_covariance(sol, variables = variables, n_lags = n_lags)
}

#' @export
model_covariance.dsge_fit <- function(x, variables = NULL,
                                      n_lags = 0L, ...) {
  sol <- x$solution
  .compute_model_covariance(sol, variables = variables, n_lags = n_lags)
}

#' @export
model_covariance.dsge_bayes <- function(x, variables = NULL,
                                        n_lags = 0L, ...) {
  # Use posterior mean parameters to solve the model
  post_mean <- colMeans(do.call(rbind, lapply(seq_len(dim(x$posterior)[3]),
    function(ch) x$posterior[, , ch])))

  # Extract structural and shock SD parameters
  pnames <- x$param_names
  free_pars <- x$free_parameters
  shock_names <- x$shock_names

  # Build parameter vector
  model <- x$model
  all_params <- c(unlist(model$start), unlist(model$fixed))
  for (p in free_pars) {
    all_params[p] <- post_mean[p]
  }

  # Build shock_sd
  sd_names <- paste0("sd_e.", shock_names)
  shock_sd <- post_mean[sd_names]
  names(shock_sd) <- shock_names

  sol <- solve_dsge(model, params = all_params, shock_sd = shock_sd)
  .compute_model_covariance(sol, variables = variables, n_lags = n_lags)
}

#' Internal covariance computation
#' @noRd
.compute_model_covariance <- function(sol, variables = NULL, n_lags = 0L) {
  G <- sol$G
  H <- sol$H
  M <- sol$M
  D <- sol$D

  # Observable selection matrix: Z = D %*% G maps states to observables
  Z <- D %*% G
  # (paste0("y", seq_len(0)) is "y", not character(0): a model without
  # observed variables, e.g. any imported .mod file without varobs, got a
  # phantom observable and the dimnames assignment below failed)
  obs_names <- rownames(D)
  if (is.null(obs_names)) {
    obs_names <- if (nrow(D) > 0L) paste0("y", seq_len(nrow(D))) else character(0)
  }

  # All control names for optional subsetting
  ctrl_names <- rownames(G)
  if (is.null(ctrl_names)) {
    ctrl_names <- if (nrow(G) > 0L) paste0("c", seq_len(nrow(G))) else character(0)
  }

  # Rows to report: the observables by default, else the chosen controls
  if (is.null(variables)) {
    Zsel <- Z
    var_names <- obs_names
  } else {
    all_avail <- c(obs_names, ctrl_names)
    missing <- setdiff(variables, all_avail)
    if (length(missing) > 0) {
      stop("Variables not found: ", paste(missing, collapse = ", "))
    }
    idx <- match(variables, ctrl_names)
    if (any(is.na(idx))) {
      stop("Could not match all variables to control names.")
    }
    Zsel <- G[idx, , drop = FALSE]
    var_names <- variables
  }

  # With a unit root in H (a price level, a cumulated variable) the
  # unconditional state covariance does not exist, and the Lyapunov solution
  # returned a large meaningless number that also contaminated stationary
  # variables (Gali 2008, ch. 5: money growth sd 408217 instead of 33.3).
  # Covariances are then sums of products of impulse responses, which
  # converge for stationary variables; variables whose responses do not die
  # out get an infinite variance, as in Dynare.
  roots <- if (nrow(H) > 0L) abs(eigen(H, only.values = TRUE)$values) else numeric(0)
  if (length(roots) > 0L && max(roots) > 1 - 1e-6) {
    gam <- .covariance_irf_sum(Zsel, H, M, n_lags)
  } else {
    P <- compute_unconditional_P(H, M %*% t(M))
    gam <- list(Zsel %*% P %*% t(Zsel))
    HkP <- P  # H^k P, updated as H (H^{k-1} P) to avoid forming H^k
    for (k in seq_len(n_lags)) {
      HkP <- H %*% HkP
      gam[[k + 1L]] <- Zsel %*% HkP %*% t(Zsel)
    }
  }
  for (k in seq_along(gam)) dimnames(gam[[k]]) <- list(var_names, var_names)

  cov_mat <- gam[[1L]]
  # Ensure symmetry (Inf and NA entries are symmetric already)
  fin <- is.finite(cov_mat) & is.finite(t(cov_mat))
  cov_mat[fin] <- ((cov_mat + t(cov_mat)) / 2)[fin]

  # Standard deviations (Inf for a non-stationary variable)
  sds <- sqrt(pmax(diag(cov_mat), 0))

  # Correlation matrix (NA where a variance is infinite)
  cor_mat <- cov_mat
  for (i in seq_along(sds)) {
    for (j in seq_along(sds)) {
      cor_mat[i, j] <- if (!is.finite(sds[i]) || !is.finite(sds[j])) NA_real_
        else if (sds[i] > 0 && sds[j] > 0) cov_mat[i, j] / (sds[i] * sds[j])
        else 0
    }
  }
  diag(cor_mat) <- ifelse(is.finite(sds), 1, NA_real_)

  # Autocovariances at lags 1..K
  autocov <- list()
  for (k in seq_len(n_lags)) autocov[[paste0("lag_", k)]] <- gam[[k + 1L]]

  structure(
    list(
      covariance = cov_mat,
      correlation = cor_mat,
      std_dev = sds,
      autocovariances = autocov,
      variables = var_names,
      n_lags = n_lags
    ),
    class = "dsge_covariance"
  )
}

#' Covariances and autocovariances as sums of impulse-response products
#'
#' Gamma(j) = Cov(y_{t+j}, y_t) = sum_k R_{k+j} R_k' with R_k = Zsel H^k M,
#' for a transition matrix with a unit root. Rows whose responses do not
#' die out (non-stationary) get variance Inf and NA covariances.
#' @noRd
.covariance_irf_sum <- function(Zsel, H, M, n_lags = 0L, max_iter = 200000L,
                                tol = 1e-20) {
  n <- nrow(Zsel)
  # Responses after 2^24 periods: a stable root of at most 1 - 1e-6
  # (Dynare's qz_criterium) has decayed below 1e-6 by then, while a unit
  # root keeps the response alive. (2^30 periods was too far: rounding
  # moves a double unit root, a Jordan block at 1 as in an HP-filter trend,
  # about 1e-8 inside the unit circle, and over 1e9 periods it died out.)
  Hp <- H
  for (i in seq_len(24L)) Hp <- Hp %*% Hp
  r_inf <- sqrt(rowSums((Zsel %*% Hp %*% M)^2))
  r_0 <- sqrt(rowSums((Zsel %*% M)^2))
  nonstat <- !is.finite(r_inf) | r_inf > 1e-6 * max(1, max(r_0))
  stat <- !nonstat

  gam <- replicate(n_lags + 1L, matrix(0, n, n), simplify = FALSE)
  buf <- vector("list", n_lags + 1L)   # buf[[j + 1]] = R_{k - j}
  X <- M
  converged <- FALSE
  for (k in 0:max_iter) {
    R <- Zsel %*% X
    buf <- c(list(R), buf[seq_len(n_lags)])
    for (j in 0:n_lags) {
      if (!is.null(buf[[j + 1L]])) gam[[j + 1L]] <- gam[[j + 1L]] + R %*% t(buf[[j + 1L]])
    }
    X <- H %*% X
    if (k >= n_lags + 10L) {
      inc <- rowSums(R^2)[stat]
      if (all(inc <= tol * pmax(diag(gam[[1L]])[stat], 1e-300))) {
        converged <- TRUE
        break
      }
    }
  }
  if (!converged) {
    warning("Model covariance with a unit root: the impulse-response sums ",
            "did not converge in ", max_iter, " periods.", call. = FALSE)
  }
  for (j in seq_along(gam)) {
    gam[[j]][nonstat, ] <- NA_real_
    gam[[j]][, nonstat] <- NA_real_
  }
  diag(gam[[1L]])[nonstat] <- Inf
  gam
}

#' @export
print.dsge_covariance <- function(x, digits = 4, ...) {
  cat("Model-Implied Covariance Matrix\n")
  cat("  Variables:", paste(x$variables, collapse = ", "), "\n\n")

  cat("Covariance:\n")
  print(round(x$covariance, digits))

  cat("\nCorrelation:\n")
  print(round(x$correlation, digits))

  cat("\nStandard deviations:\n")
  sd_vec <- x$std_dev
  names(sd_vec) <- x$variables
  print(round(sd_vec, digits))

  if (x$n_lags > 0) {
    cat("\nAutocovariances computed at lags 1 to", x$n_lags, "\n")
  }
  invisible(x)
}

#' @export
summary.dsge_covariance <- function(object, ...) {
  print(object, ...)
}
