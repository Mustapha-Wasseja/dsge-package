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
#' \dontrun{
#'   sol <- solve_dsge(model, params = p, shock_sd = s)
#'   mc <- model_covariance(sol)
#'   print(mc)
#'
#'   # With autocovariances
#'   mc2 <- model_covariance(sol, n_lags = 4)
#' }
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
  obs_names <- rownames(D)
  if (is.null(obs_names)) obs_names <- paste0("y", seq_len(nrow(D)))

  # All control names for optional subsetting
  ctrl_names <- rownames(G)
  if (is.null(ctrl_names)) ctrl_names <- paste0("c", seq_len(nrow(G)))

  # Shock covariance
  Q <- M %*% t(M)

  # Unconditional state covariance: P = H P H' + Q

  P <- compute_unconditional_P(H, Q)

  # Observable covariance: Gamma(0) = Z P Z'
  Gamma0_obs <- Z %*% P %*% t(Z)
  rownames(Gamma0_obs) <- colnames(Gamma0_obs) <- obs_names

  # Full control covariance: G P G'
  Gamma0_all <- G %*% P %*% t(G)
  rownames(Gamma0_all) <- colnames(Gamma0_all) <- ctrl_names

  # Determine which matrix to use
  if (is.null(variables)) {
    # Default: observables only
    cov_mat <- Gamma0_obs
    var_names <- obs_names
  } else {
    # User-selected variables
    # Check which are available
    all_avail <- c(obs_names, ctrl_names)
    missing <- setdiff(variables, all_avail)
    if (length(missing) > 0) {
      stop("Variables not found: ", paste(missing, collapse = ", "))
    }

    # Map to indices in the full control covariance
    idx <- match(variables, ctrl_names)
    if (any(is.na(idx))) {
      stop("Could not match all variables to control names.")
    }
    cov_mat <- Gamma0_all[idx, idx, drop = FALSE]
    var_names <- variables
  }

  # Ensure symmetry
  cov_mat <- (cov_mat + t(cov_mat)) / 2

  # Standard deviations
  sds <- sqrt(pmax(diag(cov_mat), 0))

  # Correlation matrix
  cor_mat <- cov_mat
  for (i in seq_along(sds)) {
    for (j in seq_along(sds)) {
      if (sds[i] > 0 && sds[j] > 0) {
        cor_mat[i, j] <- cov_mat[i, j] / (sds[i] * sds[j])
      } else {
        cor_mat[i, j] <- 0
      }
    }
  }
  diag(cor_mat) <- 1

  # Autocovariances at lags 1..K
  autocov <- list()
  if (n_lags > 0) {
    n_s <- ncol(H)
    Hk <- diag(n_s)

    # Determine the selection matrix for requested variables
    if (is.null(variables)) {
      Zsel <- Z
    } else {
      idx <- match(variables, ctrl_names)
      Zsel <- G[idx, , drop = FALSE]
    }

    for (k in seq_len(n_lags)) {
      Hk <- Hk %*% H
      Gk <- Zsel %*% Hk %*% P %*% t(Zsel)
      rownames(Gk) <- colnames(Gk) <- var_names
      autocov[[paste0("lag_", k)]] <- Gk
    }
  }

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
