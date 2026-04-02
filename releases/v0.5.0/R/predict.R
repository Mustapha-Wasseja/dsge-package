# Prediction methods for fitted DSGE models

#' Predict Method for DSGE Models
#'
#' Computes one-step-ahead predictions or filtered state estimates
#' from a fitted DSGE model.
#'
#' @param object A `dsge_fit` object.
#' @param type Character. `"observed"` (default) for predicted values of
#'   observed control variables, or `"state"` for filtered latent state
#'   estimates.
#' @param method Character. `"onestep"` (default) for one-step-ahead
#'   predictions using only past data, or `"filter"` for filtered
#'   estimates using past and contemporaneous data.
#' @param newdata Optional new data for prediction. If `NULL`, uses
#'   the estimation data.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix of predictions or state estimates.
#'
#' @export
predict.dsge_fit <- function(object, type = c("observed", "state"),
                             method = c("onestep", "filter"),
                             newdata = NULL, ...) {
  type <- match.arg(type)
  method <- match.arg(method)

  if (!is.null(newdata)) {
    y <- prepare_data(newdata, object$model$variables$observed, demean = TRUE)
    sol <- object$solution
    kf <- kalman_filter(y, sol$G, sol$H, sol$M, sol$D)
  } else {
    kf <- object$kalman
  }

  if (type == "observed") {
    if (method == "onestep") {
      pred <- kf$predicted_obs
    } else {
      # Filtered: use filtered states to compute fitted values
      Z <- object$solution$D %*% object$solution$G
      pred <- kf$filtered_states %*% t(Z)
    }
    colnames(pred) <- object$model$variables$observed
    # Add back means
    pred <- sweep(pred, 2, object$data_means)
  } else {
    if (method == "filter") {
      # Use smoother for smoothed state estimates
      sol <- object$solution
      y_data <- if (!is.null(newdata)) {
        prepare_data(newdata, object$model$variables$observed, demean = TRUE)
      } else {
        object$data
      }
      sm <- kalman_smoother(y_data, sol$G, sol$H, sol$M, sol$D)
      pred <- sm$smoothed_states
    } else {
      pred <- kf$filtered_states
    }
    states <- c(object$model$variables$exo_state,
                object$model$variables$endo_state)
    colnames(pred) <- states
  }

  pred
}

#' Residuals from a fitted DSGE model
#'
#' Returns one-step-ahead prediction errors.
#'
#' @param object A `dsge_fit` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix of prediction errors.
#'
#' @export
residuals.dsge_fit <- function(object, ...) {
  res <- object$kalman$prediction_errors
  colnames(res) <- object$model$variables$observed
  res
}
