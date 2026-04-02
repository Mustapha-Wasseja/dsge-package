# Posterior impulse-response functions for Bayesian DSGE models

#' @describeIn irf Compute posterior IRFs from a Bayesian DSGE fit.
#'   Returns pointwise posterior median and credible bands.
#' @param n_draws Integer. Number of posterior draws to use for IRF
#'   computation. Default is 200. Set to `NULL` to use all draws.
#' @export
irf.dsge_bayes <- function(x, periods = 20L, impulse = NULL,
                           response = NULL, se = TRUE, level = 0.95,
                           n_draws = 200L, ...) {
  model <- x$model
  all_fixed <- model$fixed
  free_params <- x$free_parameters
  shock_names <- x$shock_names
  n_shocks <- length(shock_names)

  controls <- c(model$variables$observed, model$variables$unobserved)
  states <- c(model$variables$exo_state, model$variables$endo_state)
  shocks <- model$variables$exo_state
  all_vars <- c(controls, states)

  if (is.null(impulse)) impulse <- shocks
  if (is.null(response)) response <- all_vars

  # Collect all posterior draws
  post <- x$posterior
  total_draws <- dim(post)[1] * dim(post)[3]
  all_draws <- matrix(NA_real_, nrow = total_draws, ncol = dim(post)[2])
  colnames(all_draws) <- dimnames(post)[[2]]
  idx <- 0L
  for (ch in seq_len(dim(post)[3])) {
    for (i in seq_len(dim(post)[1])) {
      idx <- idx + 1L
      all_draws[idx, ] <- post[i, , ch]
    }
  }

  # Subsample if requested
  if (!is.null(n_draws) && n_draws < total_draws) {
    draw_idx <- sample.int(total_draws, n_draws)
    all_draws <- all_draws[draw_idx, , drop = FALSE]
  }
  n_use <- nrow(all_draws)

  # Template for IRF structure
  n_imp <- length(impulse)
  n_resp <- length(response)
  n_periods <- periods + 1  # 0:periods
  n_irf <- n_imp * n_resp * n_periods

  # Storage for IRF values across draws
  irf_values <- matrix(NA_real_, nrow = n_use, ncol = n_irf)

  for (d in seq_len(n_use)) {
    draw <- all_draws[d, ]
    n_free <- length(free_params)
    struct_vals <- draw[seq_len(n_free)]
    names(struct_vals) <- free_params
    all_params <- c(struct_vals, unlist(all_fixed))

    shock_sd <- draw[(n_free + 1):(n_free + n_shocks)]
    names(shock_sd) <- shock_names

    sol <- tryCatch(
      solve_dsge(model, params = all_params, shock_sd = shock_sd),
      error = function(e) NULL
    )

    if (is.null(sol) || !sol$stable) next

    irf_df <- compute_irf_values(sol$G, sol$H, sol$M, controls, states, shocks,
                                  impulse, response, periods)
    irf_values[d, ] <- irf_df$value
  }

  # Compute pointwise quantiles
  valid <- apply(irf_values, 1, function(r) all(is.finite(r)))
  if (sum(valid) < 10) {
    stop("Fewer than 10 valid posterior IRF draws. Check model specification.",
         call. = FALSE)
  }
  irf_valid <- irf_values[valid, , drop = FALSE]

  alpha <- (1 - level) / 2

  # Build output data frame from template
  # Re-create template
  template <- compute_irf_values(
    sol$G, sol$H, sol$M, controls, states, shocks,
    impulse, response, periods
  )

  template$value <- apply(irf_valid, 2, stats::median)
  template$lower <- apply(irf_valid, 2, stats::quantile, probs = alpha)
  template$upper <- apply(irf_valid, 2, stats::quantile, probs = 1 - alpha)
  template$se <- apply(irf_valid, 2, stats::sd)

  structure(
    list(data = template, periods = periods, level = level,
         n_draws = sum(valid)),
    class = "dsge_irf"
  )
}
