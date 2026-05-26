# ==========================================================================
# Global sensitivity analysis: Sobol indices and Morris elementary effects
# ==========================================================================
#
# Provides Dynare's `sensitivity` block functionality: sample parameters
# from the prior, evaluate a target output (model-implied moment, log
# likelihood, or IRF magnitude), and quantify how each parameter
# contributes to output variability:
#
#   * Sobol' (Saltelli 2010): first-order index S_i and total-effect S_T_i
#       S_i   = Var_{X_i}[E_{X_{-i}}(Y | X_i)] / Var(Y)
#       S_T_i = E_{X_{-i}}[Var_{X_i}(Y | X_{-i})] / Var(Y)
#   * Morris (1991): elementary effects mu* and sigma over r trajectories.
#
# Output: a `dsge_global_sensitivity` object with per-parameter indices
# plus a plot method (tornado of S_T_i for Sobol; mu* vs sigma scatter
# for Morris).
# ==========================================================================


#' Global Sensitivity Analysis of a DSGE Model
#'
#' Quantifies how much each free parameter contributes to variation in
#' a chosen scalar model output, sampling from the prior space.
#' Implements two complementary methods:
#' \itemize{
#'   \item \strong{Sobol' (Saltelli 2010)}: estimates first-order and
#'     total-effect variance-based sensitivity indices.  More expensive
#'     but quantitatively interpretable as "share of output variance
#'     attributable to parameter i".
#'   \item \strong{Morris (1991) elementary effects}: cheaper
#'     "screening" method that produces \eqn{\mu^*} (mean absolute
#'     elementary effect, ranking importance) and \eqn{\sigma} (spread,
#'     measuring nonlinearity / interactions).
#' }
#'
#' @param model A \code{dsge_model} object.
#' @param priors Named list of \code{dsge_prior} objects, one per free
#'   structural parameter (use the same priors you would for Bayesian
#'   estimation).  Shock SDs default to \code{inv_gamma(0.1, 2)} if not
#'   supplied; pass them under names \code{"sd_e.<shock>"} to override.
#' @param target A character string naming the scalar output to analyse:
#'   \itemize{
#'     \item \code{"sd:<var>"} -- model-implied unconditional standard
#'       deviation of observable \code{<var>}.
#'     \item \code{"cor:<v1>:<v2>"} -- model-implied unconditional
#'       correlation between \code{<v1>} and \code{<v2>}.
#'     \item \code{"irf_max:<shock>:<response>"} -- maximum absolute
#'       impulse response of \code{<response>} to \code{<shock>} over
#'       \code{horizon} periods.
#'     \item \code{"loglik"} -- Kalman-filter log likelihood (requires
#'       \code{data}).
#'   }
#' @param method \code{"sobol"} (default) or \code{"morris"}.
#' @param n_samples Integer.  Number of base samples (Sobol uses
#'   \code{n_samples * (k + 2)} model solves total; Morris uses
#'   \code{n_samples * (k + 1)}).  Default 200.
#' @param horizon Integer.  IRF horizon when \code{target} starts with
#'   \code{"irf_max:"}.  Default 20.
#' @param data Required when \code{target = "loglik"}.
#' @param seed Optional integer seed.
#'
#' @return An object of class \code{"dsge_global_sensitivity"} with the
#'   computed indices and metadata.  For Sobol the key fields are
#'   \code{S_first} (first-order) and \code{S_total} (total-effect); for
#'   Morris, \code{mu_star} and \code{sigma}.  All indices are named by
#'   parameter.
#'
#' @examples
#' \donttest{
#' nk <- dsge_model(
#'   obs(p ~ beta * lead(p) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(p) - g)),
#'   obs(r ~ psi * p + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
#' priors <- list(
#'   kappa = prior("beta",   shape1 = 2,   shape2 = 8),
#'   psi   = prior("normal", mean   = 1.5, sd = 0.25),
#'   rhou  = prior("beta",   shape1 = 5,   shape2 = 2),
#'   rhog  = prior("beta",   shape1 = 5,   shape2 = 2))
#' gs <- global_sensitivity(nk, priors, target = "sd:p",
#'                          method = "sobol", n_samples = 100, seed = 1)
#' print(gs)
#' }
#'
#' @references
#' Saltelli, A. et al. (2010). Variance based sensitivity analysis of
#'   model output. \emph{Computer Physics Communications}, 181: 259-270.
#'
#' Morris, M.D. (1991). Factorial sampling plans for preliminary
#'   computational experiments.  \emph{Technometrics}, 33(2): 161-174.
#'
#' @export
global_sensitivity <- function(model, priors, target,
                               method = c("sobol", "morris"),
                               n_samples = 200L,
                               horizon = 20L,
                               data = NULL,
                               seed = NULL) {
  if (!inherits(model, "dsge_model"))
    stop("`model` must be a dsge_model.", call. = FALSE)
  method <- match.arg(method)
  n_samples <- as.integer(n_samples)
  if (n_samples < 8L)
    stop("`n_samples` must be >= 8.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  # Parameter setup
  free_params <- model$free_parameters
  shock_names <- model$variables$exo_state
  sd_names    <- paste0("sd_e.", shock_names)
  all_names   <- c(free_params, sd_names)
  prior_list  <- validate_priors(priors, free_params, shock_names)
  k <- length(all_names)

  # Build target evaluator
  eval_target <- .gs_target_evaluator(model, target, horizon, data)

  # Sample function: takes a vector u in [0,1]^k, returns natural params
  sample_from_unit <- function(u_vec) {
    out <- numeric(k)
    for (j in seq_len(k)) {
      out[j] <- .gs_quantile_from_unit(u_vec[j], prior_list[[j]])
    }
    names(out) <- all_names
    out
  }

  eval_at <- function(u_vec) {
    p_nat <- sample_from_unit(u_vec)
    struct  <- p_nat[free_params]
    shocksd <- p_nat[sd_names]
    names(shocksd) <- shock_names
    eval_target(struct, shocksd)
  }

  if (method == "sobol") {
    res <- .gs_sobol(eval_at, k, n_samples)
    names(res$S_first) <- all_names
    names(res$S_total) <- all_names
    structure(list(method = "sobol",
                   target = target,
                   S_first = res$S_first,
                   S_total = res$S_total,
                   n_samples = n_samples,
                   n_solves  = n_samples * (k + 2),
                   param_names = all_names),
              class = "dsge_global_sensitivity")
  } else {
    res <- .gs_morris(eval_at, k, n_samples)
    names(res$mu_star) <- all_names
    names(res$sigma)   <- all_names
    structure(list(method = "morris",
                   target = target,
                   mu_star = res$mu_star,
                   sigma   = res$sigma,
                   n_samples = n_samples,
                   n_solves  = n_samples * (k + 1),
                   param_names = all_names),
              class = "dsge_global_sensitivity")
  }
}


# --------------------------------------------------------------------------
# Internal: target evaluator factory
# --------------------------------------------------------------------------

#' @noRd
.gs_target_evaluator <- function(model, target, horizon, data) {
  parts <- strsplit(target, ":", fixed = TRUE)[[1]]
  kind  <- parts[1]
  if (kind == "loglik") {
    if (is.null(data))
      stop("target = 'loglik' requires `data`.", call. = FALSE)
    obs_vars <- model$variables$observed
    y <- prepare_data(data, obs_vars, demean = TRUE)
    return(function(struct, shocksd) {
      all_params <- c(struct, unlist(model$fixed))
      sol <- tryCatch(solve_dsge(model, params = all_params,
                                 shock_sd = shocksd),
                      error = function(e) NULL)
      if (is.null(sol) || !sol$stable) return(NA_real_)
      kf <- tryCatch(kalman_filter(y, sol$G, sol$H, sol$M, sol$D),
                     error = function(e) NULL)
      if (is.null(kf)) return(NA_real_)
      kf$loglik
    })
  }

  if (kind == "sd" || kind == "cor") {
    v1 <- if (length(parts) >= 2) parts[2] else NA
    v2 <- if (length(parts) >= 3) parts[3] else NA
    return(function(struct, shocksd) {
      all_params <- c(struct, unlist(model$fixed))
      sol <- tryCatch(solve_dsge(model, params = all_params,
                                 shock_sd = shocksd),
                      error = function(e) NULL)
      if (is.null(sol) || !sol$stable) return(NA_real_)
      mc <- tryCatch(model_covariance(sol),
                     error = function(e) NULL)
      if (is.null(mc)) return(NA_real_)
      if (kind == "sd")  return(sqrt(pmax(mc$covariance[v1, v1], 0)))
      if (kind == "cor") return(mc$correlation[v1, v2])
    })
  }

  if (kind == "irf_max") {
    impulse  <- if (length(parts) >= 2) parts[2] else NA
    response <- if (length(parts) >= 3) parts[3] else NA
    return(function(struct, shocksd) {
      all_params <- c(struct, unlist(model$fixed))
      sol <- tryCatch(solve_dsge(model, params = all_params,
                                 shock_sd = shocksd),
                      error = function(e) NULL)
      if (is.null(sol) || !sol$stable) return(NA_real_)
      ir <- tryCatch(irf(sol, periods = horizon),
                     error = function(e) NULL)
      if (is.null(ir)) return(NA_real_)
      sub <- ir$data[ir$data$impulse == impulse &
                     ir$data$response == response, ]
      if (nrow(sub) == 0) return(NA_real_)
      max(abs(sub$value), na.rm = TRUE)
    })
  }
  stop("Unrecognised target: '", target, "'.", call. = FALSE)
}


#' Map a [0,1] uniform value to the natural parameter scale via the
#' prior's CDF inverse (for known parametric families) or a rough
#' substitute based on a normal approximation.
#' @noRd
.gs_quantile_from_unit <- function(u, p) {
  u <- min(max(u, 1e-6), 1 - 1e-6)
  dist <- p$distribution
  pars <- p$params
  switch(dist,
    normal    = stats::qnorm(u, pars$mean,   pars$sd),
    beta      = stats::qbeta(u, pars$shape1, pars$shape2),
    gamma     = stats::qgamma(u, shape = pars$shape, rate = pars$rate),
    uniform   = stats::qunif(u, pars$min,   pars$max),
    inv_gamma = 1 / stats::qgamma(1 - u, shape = pars$shape,
                                       scale = pars$scale),
    # Fallback: rough normal approximation around 0.5
    stats::qnorm(u))
}


# --------------------------------------------------------------------------
# Sobol indices via Saltelli's matrices
# --------------------------------------------------------------------------

#' @noRd
.gs_sobol <- function(eval_at, k, N) {
  # Two independent N x k matrices on the unit cube
  A <- matrix(stats::runif(N * k), N, k)
  B <- matrix(stats::runif(N * k), N, k)

  yA <- vapply(seq_len(N), function(n) eval_at(A[n, ]), numeric(1))
  yB <- vapply(seq_len(N), function(n) eval_at(B[n, ]), numeric(1))

  # Build C^(i): copy of B with column i from A
  S_first <- numeric(k)
  S_total <- numeric(k)
  varY    <- stats::var(c(yA, yB), na.rm = TRUE)
  if (!is.finite(varY) || varY <= 0) {
    return(list(S_first = rep(NA_real_, k),
                S_total = rep(NA_real_, k)))
  }

  for (i in seq_len(k)) {
    Ci <- B; Ci[, i] <- A[, i]
    yCi <- vapply(seq_len(N), function(n) eval_at(Ci[n, ]), numeric(1))
    keep <- is.finite(yA) & is.finite(yB) & is.finite(yCi)
    if (sum(keep) < 4) {
      S_first[i] <- NA_real_; S_total[i] <- NA_real_; next
    }
    # Saltelli 2010 estimators
    S_first[i] <- mean(yB[keep] * (yCi[keep] - yA[keep])) / varY
    S_total[i] <- 0.5 * mean((yA[keep] - yCi[keep])^2) / varY
  }
  list(S_first = S_first, S_total = S_total)
}


# --------------------------------------------------------------------------
# Morris elementary effects
# --------------------------------------------------------------------------

#' @noRd
.gs_morris <- function(eval_at, k, r) {
  delta <- 0.5
  # r trajectories, each of length k + 1
  EE <- matrix(NA_real_, r, k)
  for (t in seq_len(r)) {
    base <- stats::runif(k)
    y_base <- eval_at(base)
    if (!is.finite(y_base)) next
    perm <- sample.int(k)
    cur <- base
    y_prev <- y_base
    for (j in perm) {
      cur_new <- cur
      cur_new[j] <- min(1, max(0, cur[j] + delta))
      y_new <- eval_at(cur_new)
      if (!is.finite(y_new)) { cur <- cur_new; next }
      EE[t, j] <- (y_new - y_prev) / (cur_new[j] - cur[j])
      cur <- cur_new
      y_prev <- y_new
    }
  }
  mu_star <- apply(abs(EE), 2, mean, na.rm = TRUE)
  sigma   <- apply(EE,      2, stats::sd, na.rm = TRUE)
  list(mu_star = mu_star, sigma = sigma)
}


# --------------------------------------------------------------------------
# Print / plot
# --------------------------------------------------------------------------

#' @export
print.dsge_global_sensitivity <- function(x, digits = 4, ...) {
  cat(sprintf("Global sensitivity analysis (%s) on target '%s'\n",
              x$method, x$target))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat(sprintf("Samples: %d (total model solves: %d)\n",
              x$n_samples, x$n_solves))
  if (x$method == "sobol") {
    tab <- data.frame(parameter = x$param_names,
                      S_first   = round(x$S_first, digits),
                      S_total   = round(x$S_total, digits))
  } else {
    tab <- data.frame(parameter = x$param_names,
                      mu_star   = round(x$mu_star, digits),
                      sigma     = round(x$sigma,   digits))
  }
  print(tab, row.names = FALSE)
  invisible(x)
}


#' @export
plot.dsge_global_sensitivity <- function(x, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .dsge_par_grid(1L, 1L)
  graphics::par(mar = c(4.4, 7.8, 2.4, 0.8))
  if (x$method == "sobol") {
    vals <- x$S_total
    vals[!is.finite(vals)] <- 0
    ord <- order(vals, decreasing = FALSE)
    graphics::barplot(vals[ord], horiz = TRUE, las = 1,
                      names.arg = x$param_names[ord],
                      col = .DSGE_INK_PRIMARY, border = "white",
                      main = paste("Sobol total-effect indices --", x$target),
                      xlab = "S_total")
  } else {
    plot(x$mu_star, x$sigma,
         xlab = expression(mu^"*"), ylab = expression(sigma),
         main = paste("Morris elementary effects --", x$target),
         pch = 19, col = .DSGE_INK_PRIMARY)
    graphics::text(x$mu_star, x$sigma, labels = x$param_names,
                   pos = 4, cex = 0.7)
  }
  invisible(x)
}
