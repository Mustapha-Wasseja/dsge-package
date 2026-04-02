# ==========================================================================
# Prediction and reporting tools
# ==========================================================================

#' Fitted values from a DSGE model
#'
#' Returns filtered fitted values of observed variables (using all
#' information up to and including time t).
#'
#' @param object A \code{"dsge_fit"} object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix of fitted values with columns named by observed variables.
#'
#' @export
fitted.dsge_fit <- function(object, ...) {
  predict(object, type = "observed", method = "filter")
}

#' Prediction accuracy measures for a fitted DSGE model
#'
#' Computes root mean squared error (RMSE), mean absolute error (MAE),
#' and mean error (bias) of one-step-ahead predictions.
#'
#' @param object A \code{"dsge_fit"} object.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"dsge_prediction_accuracy"} containing:
#' \describe{
#'   \item{rmse}{Named numeric vector of RMSE by variable.}
#'   \item{mae}{Named numeric vector of MAE by variable.}
#'   \item{bias}{Named numeric vector of mean error by variable.}
#'   \item{n}{Number of observations.}
#'   \item{variables}{Variable names.}
#' }
#'
#' @examples
#' \dontrun{
#'   fit <- estimate(model, data = dat)
#'   acc <- prediction_accuracy(fit)
#'   print(acc)
#' }
#'
#' @export
prediction_accuracy <- function(object, ...) {
  UseMethod("prediction_accuracy")
}

#' @export
prediction_accuracy.dsge_fit <- function(object, ...) {
  res <- residuals(object)
  obs_names <- colnames(res)
  n <- nrow(res)

  rmse <- sqrt(colMeans(res^2))
  mae <- colMeans(abs(res))
  bias <- colMeans(res)

  names(rmse) <- names(mae) <- names(bias) <- obs_names

  structure(
    list(
      rmse = rmse,
      mae = mae,
      bias = bias,
      n = n,
      variables = obs_names
    ),
    class = "dsge_prediction_accuracy"
  )
}

#' @export
print.dsge_prediction_accuracy <- function(x, digits = 4, ...) {
  cat("One-Step-Ahead Prediction Accuracy\n")
  cat("  Observations:", x$n, "\n\n")

  tbl <- data.frame(
    RMSE = x$rmse,
    MAE = x$mae,
    Bias = x$bias,
    row.names = x$variables
  )
  print(round(tbl, digits))
  invisible(x)
}

#' Prediction standard errors from Kalman filter
#'
#' Extracts the standard errors of one-step-ahead predictions from the
#' innovation variance (F_t) stored by the Kalman filter.
#'
#' @param object A \code{"dsge_fit"} object.
#'
#' @return A matrix of prediction standard errors (T x n_obs).
#'
#' @noRd
prediction_se <- function(object) {
  kf <- object$kalman
  if (is.null(kf$innovation_var) || length(kf$innovation_var) == 0) {
    stop("Innovation variances not available. Re-estimate the model.")
  }

  n_T <- length(kf$innovation_var)
  n_obs <- nrow(kf$innovation_var[[1]])

  se_mat <- matrix(NA_real_, n_T, n_obs)
  for (t in seq_len(n_T)) {
    F_t <- kf$innovation_var[[t]]
    if (!is.null(F_t)) {
      se_mat[t, ] <- sqrt(pmax(diag(F_t), 0))
    }
  }
  colnames(se_mat) <- object$model$variables$observed
  se_mat
}

#' Prediction intervals for DSGE models
#'
#' Computes point predictions and prediction intervals using the
#' one-step-ahead innovation variance from the Kalman filter.
#'
#' @param object A \code{"dsge_fit"} object.
#' @param level Confidence level for prediction intervals (default 0.95).
#' @param ... Additional arguments passed to \code{predict()}.
#'
#' @return An object of class \code{"dsge_prediction_interval"} with
#'   components \code{fit}, \code{lower}, \code{upper}, \code{se},
#'   \code{level}, and \code{variables}.
#'
#' @examples
#' \dontrun{
#'   fit <- estimate(model, data = dat)
#'   pi <- prediction_interval(fit, level = 0.95)
#'   print(pi)
#' }
#'
#' @export
prediction_interval <- function(object, level = 0.95, ...) {
  UseMethod("prediction_interval")
}

#' @export
prediction_interval.dsge_fit <- function(object, level = 0.95, ...) {
  pred <- predict(object, type = "observed", method = "onestep")
  se <- prediction_se(object)

  z <- qnorm((1 + level) / 2)
  lower <- pred - z * se
  upper <- pred + z * se

  structure(
    list(
      fit = pred,
      se = se,
      lower = lower,
      upper = upper,
      level = level,
      variables = colnames(pred)
    ),
    class = "dsge_prediction_interval"
  )
}

#' @export
print.dsge_prediction_interval <- function(x, digits = 4, ...) {
  cat(sprintf("One-Step-Ahead Prediction Intervals (%.0f%%)\n",
              x$level * 100))
  cat("  Variables:", paste(x$variables, collapse = ", "), "\n")
  cat("  Observations:", nrow(x$fit), "\n\n")

  # Show summary: mean SE by variable
  cat("Mean prediction SE by variable:\n")
  mean_se <- colMeans(x$se, na.rm = TRUE)
  names(mean_se) <- x$variables
  print(round(mean_se, digits))
  invisible(x)
}
