#' Perfect Foresight / Deterministic Transition Paths
#'
#' Simulate deterministic transition paths for DSGE models under perfect
#' foresight. Supports temporary shocks, permanent shocks, and initial
#' condition experiments using the linearized solution.
#'
#' @param x A solved DSGE model object. Can be a \code{dsge_solution} (from
#'   \code{solve_dsge}), \code{dsge_fit} (from \code{estimate}), or
#'   \code{dsge_bayes} (from \code{bayes_dsge}).
#' @param shocks Deterministic shock specification. Can be:
#'   \itemize{
#'     \item A named list where each element is a numeric vector giving the
#'       shock path for that shock variable (e.g., \code{list(u = c(1, 0, 0))}).
#'       Unnamed periods after the vector ends are treated as zero.
#'     \item A matrix of dimension \code{horizon x n_shocks} with column names
#'       matching shock names.
#'   }
#'   Shock values are in units of the shock variable (not standard deviations).
#'   If \code{NULL} (default), no shocks are applied (useful with
#'   \code{initial} to study convergence from displaced initial conditions).
#' @param initial Named numeric vector of initial state deviations from steady
#'   state. Names must match state variable names. Unspecified states default
#'   to zero. Default is \code{NULL} (all states start at steady state).
#' @param horizon Integer. Number of periods to simulate. Default 40.
#' @param params Named numeric vector of parameters. Required only when
#'   \code{x} is a model object that has not been solved. For solved objects,
#'   parameters are extracted automatically.
#' @param shock_sd Named numeric vector of shock standard deviations. Used
#'   only when \code{x} is a model object. Default \code{NULL}.
#' @param in_sd Logical. If \code{TRUE}, shock values in \code{shocks} are
#'   interpreted as multiples of the shock standard deviation. Default
#'   \code{FALSE} (shocks are in level units).
#'
#' @return An object of class \code{"dsge_perfect_foresight"} containing:
#'   \describe{
#'     \item{states}{Matrix (horizon x n_states) of state deviations from SS}
#'     \item{controls}{Matrix (horizon x n_controls) of control deviations}
#'     \item{state_levels}{Matrix of state levels (SS + deviation), if SS available}
#'     \item{control_levels}{Matrix of control levels, if SS available}
#'     \item{steady_state}{Named numeric vector of steady-state values}
#'     \item{shock_path}{Matrix (horizon x n_shocks) of applied shocks}
#'     \item{initial}{Named vector of initial state deviations}
#'     \item{horizon}{Integer horizon}
#'     \item{state_names}{Character vector of state names}
#'     \item{control_names}{Character vector of control names}
#'     \item{shock_names}{Character vector of shock names}
#'     \item{H}{State transition matrix used}
#'     \item{G}{Policy matrix used}
#'     \item{M}{Shock impact matrix used}
#'   }
#'
#' @details
#' The deterministic transition path is computed using the linearized
#' state-space representation:
#' \deqn{x_{t+1} = H x_t + M \varepsilon_{t+1}}
#' \deqn{y_t = G x_t}
#' where \eqn{x_t} are state deviations from steady state, \eqn{y_t} are
#' control deviations, and \eqn{\varepsilon_t} are deterministic shocks.
#'
#' This uses the first-order linearized solution, so results are approximate
#' for large shocks. For small to moderate shocks, the linearized paths are
#' accurate.
#'
#' @examples
#' # Simple AR(1) model
#' mod <- dsge_model(
#'   obs(p ~ x),
#'   state(x ~ rho * x),
#'   start = list(rho = 0.9)
#' )
#' sol <- solve_dsge(mod, params = list(rho = 0.9), shock_sd = c(x = 0.01))
#'
#' # One-time shock at period 1
#' pf <- perfect_foresight(sol, shocks = list(x = 0.01), horizon = 40)
#' plot(pf)
#'
#' # Displaced initial condition
#' pf2 <- perfect_foresight(sol, initial = c(x = 0.05), horizon = 40)
#' plot(pf2)
#'
#' @export
perfect_foresight <- function(x, shocks = NULL, initial = NULL,
                              horizon = 40L, params = NULL,
                              shock_sd = NULL, in_sd = FALSE) {
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("horizon must be at least 1")

  # --- Extract solution matrices ---
  sol <- .extract_solution(x, params, shock_sd)
  H <- sol$H
  G <- sol$G
  M <- sol$M
  ss <- sol$steady_state

  n_states <- nrow(H)
  n_controls <- nrow(G)
  n_shocks <- ncol(M)

  state_names <- rownames(H)
  control_names <- rownames(G)
  shock_names <- colnames(M)

  if (is.null(state_names)) state_names <- paste0("state_", seq_len(n_states))
  if (is.null(control_names)) control_names <- paste0("ctrl_", seq_len(n_controls))
  if (is.null(shock_names)) shock_names <- paste0("shock_", seq_len(n_shocks))

  # --- Build shock path matrix ---
  shock_path <- matrix(0, nrow = horizon, ncol = n_shocks)
  colnames(shock_path) <- shock_names

  if (!is.null(shocks)) {
    if (is.matrix(shocks)) {
      if (ncol(shocks) != n_shocks) {
        stop("shocks matrix must have ", n_shocks, " columns")
      }
      nr <- min(nrow(shocks), horizon)
      shock_path[seq_len(nr), ] <- shocks[seq_len(nr), , drop = FALSE]
      if (!is.null(colnames(shocks))) colnames(shock_path) <- colnames(shocks)
    } else if (is.list(shocks)) {
      for (nm in names(shocks)) {
        idx <- match(nm, shock_names)
        if (is.na(idx)) stop("Unknown shock name: '", nm, "'. Available: ",
                             paste(shock_names, collapse = ", "))
        vals <- as.numeric(shocks[[nm]])
        nr <- min(length(vals), horizon)
        shock_path[seq_len(nr), idx] <- vals[seq_len(nr)]
      }
    } else {
      stop("shocks must be a named list or a matrix")
    }

    # Scale by shock_sd if in_sd = TRUE
    if (in_sd) {
      sds <- sol$shock_sd
      if (is.null(sds) || length(sds) != n_shocks) {
        stop("in_sd = TRUE requires shock standard deviations")
      }
      shock_path <- sweep(shock_path, 2, sds, "*")
    }
  }

  # --- Build initial state vector ---
  x0 <- rep(0, n_states)
  names(x0) <- state_names

  if (!is.null(initial)) {
    for (nm in names(initial)) {
      idx <- match(nm, state_names)
      if (is.na(idx)) stop("Unknown state name: '", nm, "'. Available: ",
                           paste(state_names, collapse = ", "))
      x0[idx] <- initial[nm]
    }
  }

  # --- Forward simulation ---
  states_mat <- matrix(0, nrow = horizon, ncol = n_states)
  colnames(states_mat) <- state_names

  x_curr <- x0
  for (t in seq_len(horizon)) {
    # State evolution: x_t = H * x_{t-1} + M * eps_t
    # (for t=1, x_curr = x0 which is the initial condition)
    if (t == 1L) {
      # Period 1: apply initial condition + period-1 shock
      x_curr <- as.numeric(H %*% x0 + M %*% shock_path[1L, ])
    } else {
      x_curr <- as.numeric(H %*% x_curr + M %*% shock_path[t, ])
    }
    states_mat[t, ] <- x_curr
  }

  # Controls: y_t = G * x_t
  controls_mat <- states_mat %*% t(G)
  colnames(controls_mat) <- control_names

  # --- Levels (if steady state available) ---
  state_levels <- NULL
  control_levels <- NULL

  if (!is.null(ss)) {
    # Match SS names to state/control names
    ss_states <- ss[state_names]
    ss_controls <- ss[control_names]

    if (!any(is.na(ss_states))) {
      state_levels <- sweep(states_mat, 2, ss_states, "+")
    }
    if (!any(is.na(ss_controls))) {
      control_levels <- sweep(controls_mat, 2, ss_controls, "+")
    }
  }

  # --- Return object ---
  result <- list(
    states = states_mat,
    controls = controls_mat,
    state_levels = state_levels,
    control_levels = control_levels,
    steady_state = ss,
    shock_path = shock_path,
    initial = x0,
    horizon = horizon,
    state_names = state_names,
    control_names = control_names,
    shock_names = shock_names,
    H = H,
    G = G,
    M = M
  )
  class(result) <- "dsge_perfect_foresight"
  result
}


# --- Helper: extract solution from various input types ---
.extract_solution <- function(x, params = NULL, shock_sd = NULL) {
  if (inherits(x, "dsge_solution")) {
    return(list(
      H = x$H, G = x$G, M = x$M,
      steady_state = x$steady_state,
      shock_sd = x$shock_sd
    ))
  }

  if (inherits(x, "dsge_fit")) {
    sol <- x$solution
    return(list(
      H = sol$H, G = sol$G, M = sol$M,
      steady_state = sol$steady_state,
      shock_sd = sol$shock_sd
    ))
  }

  if (inherits(x, "dsge_bayes")) {
    # Use posterior mean parameters to solve
    post <- x$posterior
    pnames <- x$param_names
    post_means <- numeric(length(pnames))
    names(post_means) <- pnames
    for (p in pnames) {
      post_means[p] <- mean(as.numeric(post[, p, ]))
    }

    # Build full parameter vector
    model <- x$model
    all_params <- .build_params_from_bayes(x, post_means)
    shock_sds <- .extract_shock_sds(post_means, x$shock_names)

    sol <- solve_dsge(model, params = all_params, shock_sd = shock_sds)
    return(list(
      H = sol$H, G = sol$G, M = sol$M,
      steady_state = sol$steady_state,
      shock_sd = sol$shock_sd
    ))
  }

  if (inherits(x, "dsge_model") || inherits(x, "dsgenl_model")) {
    if (is.null(params)) stop("params required when passing a model object")
    if (is.null(shock_sd)) shock_sd <- rep(1, length(x$variables$exo_state))
    sol <- solve_dsge(x, params = params, shock_sd = shock_sd)
    return(list(
      H = sol$H, G = sol$G, M = sol$M,
      steady_state = sol$steady_state,
      shock_sd = sol$shock_sd
    ))
  }

  stop("x must be a dsge_solution, dsge_fit, dsge_bayes, dsge_model, or dsgenl_model")
}


# Build full params from bayes fit + posterior means
.build_params_from_bayes <- function(fit, post_means) {
  model <- fit$model
  free_names <- fit$free_parameters
  sd_names <- paste0("sd_e.", fit$shock_names)

  # Start with fixed params from model
  if (inherits(model, "dsgenl_model")) {
    all_p <- unlist(model$fixed)
    for (p in free_names) {
      all_p[p] <- post_means[p]
    }
  } else {
    all_p <- model$params
    for (p in free_names) {
      all_p[p] <- post_means[p]
    }
  }
  all_p
}


# Extract shock SDs from posterior means
.extract_shock_sds <- function(post_means, shock_names) {
  sd_names <- paste0("sd_e.", shock_names)
  sds <- post_means[sd_names]
  names(sds) <- shock_names
  sds
}


# ======================================================================
#  Nonlinear perfect foresight  (stacked-time Newton / LBJ algorithm)
# ======================================================================

#' Perfect Foresight for Nonlinear DSGE Models
#'
#' Computes deterministic perfect foresight transition paths for nonlinear
#' DSGE models using a stacked-time Newton solver (Juillard et al., 1998).
#' Unlike the linearized \code{perfect_foresight()}, this function solves the
#' full nonlinear equilibrium conditions simultaneously over the entire
#' horizon, giving exact (up to Newton tolerance) paths even for large shocks.
#'
#' @param model A \code{dsgenl_model} object.
#' @param params Named numeric vector or list of free parameter values.
#' @param shock_sd Named numeric vector of shock standard deviations with
#'   names matching \code{model$variables$exo_state}. Used only for the
#'   linearized warm-start step; the nonlinear solver uses \code{shocks}
#'   directly.
#' @param shocks Shock path specification (same interface as
#'   \code{\link{perfect_foresight}}): a named list where each element is
#'   a numeric vector giving the shock values for that exogenous state over
#'   time, or a matrix of dimension \code{horizon x n_exo} with column names
#'   matching exogenous state names.  Unspecified periods are zero.
#' @param initial Named numeric vector of initial DEVIATIONS from steady
#'   state for state variables.  Defaults to zero (all states start at SS).
#' @param horizon Integer. Number of periods to simulate.  Default 40.
#' @param tol Convergence tolerance (maximum absolute equation residual).
#'   Default \code{1e-8}.
#' @param max_iter Maximum Newton iterations.  Default 50.
#' @param verbose Logical.  If \code{TRUE}, prints iteration progress.
#'   Default \code{FALSE}.
#'
#' @return A \code{"dsge_perfect_foresight"} object with the same fields as
#'   \code{\link{perfect_foresight}}, so all \code{plot}, \code{print}, and
#'   \code{summary} methods apply unchanged.  Three extra fields are appended:
#'   \describe{
#'     \item{nonlinear}{\code{TRUE}}
#'     \item{newton_iters}{Number of Newton iterations taken}
#'     \item{converged}{Logical; \code{TRUE} if \code{max|R| < tol}}
#'   }
#'
#' @details
#' The algorithm (Juillard et al., 1998) stacks the \eqn{n \times T} nonlinear
#' equilibrium conditions into a single system and solves it by Newton's method.
#' The Jacobian is block-bidiagonal (each period's equations depend only on
#' the current and next-period variables) and is assembled numerically using
#' forward finite differences.  The block structure is exploited via block
#' back-substitution, giving \eqn{O(T n^3)} cost per Newton step.  Armijo
#' backtracking stabilises convergence for large initial steps.
#'
#' The initial state \eqn{x_1} is pinned to \code{steady_state + initial}
#' via equality constraints that replace the state equations at \eqn{t=1}.
#' For \eqn{t = T+1} the terminal condition \eqn{v_{T+1} = \bar{v}} (steady
#' state) is imposed.  The solver is warm-started from the linearized
#' \code{perfect_foresight()} path, which is exact for small shocks and
#' provides a good initial guess for large ones.
#'
#' @references
#' Juillard, M., Laxton, D., McAdam, P. and Pioro, H. (1998).
#' An algorithm competition: First-order iterations versus Newton-based
#' techniques. \emph{Journal of Economic Dynamics and Control}, 22, 1291--1318.
#'
#' @seealso \code{\link{perfect_foresight}} for linearized paths,
#'   \code{\link{dsgenl_model}} for nonlinear model specification.
#'
#' @examples
#' rbc <- dsgenl_model(
#'   "1/C = beta / C(+1) * (alpha * exp(Z) * K^(alpha-1) + 1 - delta)",
#'   "K(+1) = exp(Z) * K^alpha - C + (1 - delta) * K",
#'   "Z(+1) = rho * Z",
#'   observed = "C", endo_state = "K", exo_state = "Z",
#'   fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025),
#'   start = list(rho = 0.9)
#' )
#' pf <- perfect_foresight_nonlinear(
#'   rbc, params = c(rho = 0.9), shock_sd = c(Z = 0.01),
#'   shocks = list(Z = 0.1), horizon = 40
#' )
#' plot(pf)
#'
#' @importFrom stats setNames
#' @export
perfect_foresight_nonlinear <- function(model, params, shock_sd,
                                         shocks   = NULL,
                                         initial  = NULL,
                                         horizon  = 40L,
                                         tol      = 1e-8,
                                         max_iter = 50L,
                                         verbose  = FALSE) {

  if (!inherits(model, "dsgenl_model"))
    stop("model must be a dsgenl_model object", call. = FALSE)

  horizon  <- as.integer(horizon)
  max_iter <- as.integer(max_iter)
  if (horizon < 1L) stop("horizon must be at least 1", call. = FALSE)

  # ------------------------------------------------------------------
  # 1.  Merge fixed and free parameters
  # ------------------------------------------------------------------
  if (is.list(params)) params <- unlist(params)

  all_params <- unlist(model$fixed)
  for (nm in model$free_parameters) {
    if (!nm %in% names(params))
      stop("Free parameter '", nm, "' not supplied in params", call. = FALSE)
    all_params[nm] <- params[nm]
  }

  # ------------------------------------------------------------------
  # 2.  First-order solution (steady state + warm start matrices)
  # ------------------------------------------------------------------
  sol1 <- tryCatch(
    solve_dsge(model, params = all_params, shock_sd = shock_sd, order = 1L),
    error = function(e)
      stop("Linearization failed: ", e$message, call. = FALSE)
  )

  ss_vec <- sol1$steady_state   # named vector, all variables at SS

  ctrl_names  <- model$controls
  state_names <- model$states
  exo_names   <- model$variables$exo_state

  n_c    <- length(ctrl_names)
  n_s    <- length(state_names)
  n_exo  <- length(exo_names)
  n_vars <- n_c + n_s            # unknowns per period
  n_eq   <- n_vars               # equations per period (square system)

  ss_ctrl  <- ss_vec[ctrl_names]
  ss_state <- ss_vec[state_names]
  ss_v     <- c(ss_ctrl, ss_state)

  # Indices of the state components within a single period's sub-vector
  state_idx_in_v <- n_c + seq_len(n_s)   # positions n_c+1 .. n_vars

  # ------------------------------------------------------------------
  # 3.  Shock path  (horizon x n_exo, in units of shock variable)
  # ------------------------------------------------------------------
  shock_path <- matrix(0, nrow = horizon, ncol = n_exo,
                       dimnames = list(NULL, exo_names))

  if (!is.null(shocks)) {
    if (is.matrix(shocks)) {
      nr <- min(nrow(shocks), horizon)
      if (!is.null(colnames(shocks))) {
        for (nm in colnames(shocks)) {
          j <- match(nm, exo_names)
          if (is.na(j))
            stop("Unknown exogenous state '", nm, "' in shocks", call. = FALSE)
          shock_path[seq_len(nr), j] <- shocks[seq_len(nr), nm]
        }
      } else {
        nc_s <- min(ncol(shocks), n_exo)
        shock_path[seq_len(nr), seq_len(nc_s)] <-
          shocks[seq_len(nr), seq_len(nc_s), drop = FALSE]
      }
    } else if (is.list(shocks)) {
      for (nm in names(shocks)) {
        j <- match(nm, exo_names)
        if (is.na(j))
          stop("Unknown exogenous state '", nm, "' in shocks", call. = FALSE)
        vals <- as.numeric(shocks[[nm]])
        nr   <- min(length(vals), horizon)
        shock_path[seq_len(nr), j] <- vals[seq_len(nr)]
      }
    } else {
      stop("shocks must be a named list or matrix", call. = FALSE)
    }
  }

  # ------------------------------------------------------------------
  # 4.  Initial state in LEVELS
  #
  #  Convention (matches perfect_foresight()):
  #    shock_path[t, j]  = eps_t,j  which  *produces*  x_t  from x_{t-1}.
  #
  #  For exo states: x_1[j] = SS_j + shock_path[1, j]   (direct impact)
  #  For endo states: x_1[j] = SS_j + initial[j]          (user-supplied)
  #  The shock at period t >= 2 enters the state equation at eval-period t-1:
  #    x_t = f_s(y_{t-1}, x_{t-1}) + shock_t
  #  so in .eval_period called for period t-1 we pass shock_path[t, ].
  # ------------------------------------------------------------------
  x_init <- ss_state

  # User-specified initial deviations (typically for endo states)
  if (!is.null(initial)) {
    for (nm in names(initial)) {
      idx <- match(nm, state_names)
      if (is.na(idx))
        stop("Unknown state variable '", nm, "' in initial. Available: ",
             paste(state_names, collapse = ", "), call. = FALSE)
      x_init[idx] <- ss_state[idx] + initial[nm]
    }
  }

  # Period-1 shock directly displaces the exogenous states at t = 1
  for (j in seq_len(n_exo)) {
    ej         <- match(exo_names[j], state_names)
    x_init[ej] <- x_init[ej] + shock_path[1L, j]
  }

  # ------------------------------------------------------------------
  # 5.  Initialise z in LEVELS; warm-start from linear perfect_foresight()
  #     z = [v_1, ..., v_T],  v_t = c(controls_t, states_t)
  # ------------------------------------------------------------------
  z <- rep(ss_v, horizon)

  pf_lin <- tryCatch({
    shock_list_lin <- lapply(seq_len(n_exo), function(j) shock_path[, j])
    names(shock_list_lin) <- exo_names

    init_vec <- NULL
    if (!is.null(initial)) {
      init_vec <- setNames(rep(0, n_s), state_names)
      for (nm in names(initial)) init_vec[nm] <- initial[nm]
    }
    perfect_foresight(sol1,
                      shocks  = shock_list_lin,
                      initial = init_vec,
                      horizon = horizon)
  }, error = function(e) NULL)

  if (!is.null(pf_lin)) {
    for (t in seq_len(horizon)) {
      y_lev <- ss_ctrl  + pf_lin$controls[t, ctrl_names]
      x_lev <- ss_state + pf_lin$states[t,   state_names]
      z[.pf_idx(t, n_vars)] <- c(y_lev, x_lev)
    }
  }

  # Force z[state part of period 1] = x_init so the warm start is consistent
  z[state_idx_in_v] <- x_init

  # ------------------------------------------------------------------
  # 6.  Core building blocks
  # ------------------------------------------------------------------

  # .eval_period: evaluate nonlinear residuals for one period.
  #   v_t   : level vector [controls, states]        length n_vars
  #   v_t1  : level vector for NEXT period (or SS)   length n_vars
  #   sp    : additive shock to exo-state equations  length n_exo
  #           (this is shock_path[t+1, ] in the main loop, producing x_{t+1})
  .eval_period <- function(v_t, v_t1, sp) {
    vals <- c(
      setNames(v_t[seq_len(n_c)],       ctrl_names),
      setNames(v_t[state_idx_in_v],      state_names),
      setNames(v_t1[seq_len(n_c)],       paste0(ctrl_names,  "__f")),
      setNames(v_t1[state_idx_in_v],     paste0(state_names, "__f")),
      all_params
    )
    R <- model$eval_fn(vals)
    # Exo-state equations occupy rows n_c+1 .. n_c+n_exo
    if (n_exo > 0L) R[n_c + seq_len(n_exo)] <- R[n_c + seq_len(n_exo)] - sp
    R
  }

  # Shock argument to .eval_period for eval-period t:
  # shock_path[t+1, ] produces x_{t+1} (or zeros if beyond horizon)
  .sp_for_eval <- function(t) {
    if (t < horizon) shock_path[t + 1L, ] else rep(0, n_exo)
  }

  # Apply pinning override for t=1: replace state-equation rows with
  # x_1 - x_init = 0.  This is a helper used in both .residual_full and
  # .block_jacobian to ensure consistent treatment.
  .pin1 <- function(R, v_t) {
    R[state_idx_in_v] <- v_t[state_idx_in_v] - x_init
    R
  }

  # Full stacked residual vector of length n_eq * horizon.
  # State equations at t=1 are pinned to x_init (initial condition).
  .residual_full <- function(z) {
    R_full <- numeric(n_eq * horizon)
    for (t in seq_len(horizon)) {
      v_t  <- z[.pf_idx(t, n_vars)]
      v_t1 <- if (t < horizon) z[.pf_idx(t + 1L, n_vars)] else ss_v
      Rt   <- .eval_period(v_t, v_t1, .sp_for_eval(t))
      if (t == 1L) Rt <- .pin1(Rt, v_t)
      R_full[.pf_idx(t, n_eq)] <- Rt
    }
    R_full
  }

  # Block Jacobian via forward finite differences.
  # At t=1, numerical differentiation inherits the pinning override, giving:
  #   A_1[state rows, ctrl cols] = 0
  #   A_1[state rows, state cols] = I      (from pinning: d(x_1-x_init)/dx_1 = I)
  #   B_1[state rows, *]         = 0      (pinning doesn't depend on v_2)
  .block_jacobian <- function(z) {
    h <- 1e-5
    A <- vector("list", horizon)
    B <- vector("list", horizon - 1L)

    for (t in seq_len(horizon)) {
      v_t  <- z[.pf_idx(t, n_vars)]
      v_t1 <- if (t < horizon) z[.pf_idx(t + 1L, n_vars)] else ss_v
      sp   <- .sp_for_eval(t)

      R0 <- .eval_period(v_t, v_t1, sp)
      if (t == 1L) R0 <- .pin1(R0, v_t)

      # Diagonal block A_t = dF_t / dv_t
      At <- matrix(0, n_eq, n_vars)
      for (k in seq_len(n_vars)) {
        vp      <- v_t; vp[k] <- vp[k] + h
        Rp      <- .eval_period(vp, v_t1, sp)
        if (t == 1L) Rp <- .pin1(Rp, vp)   # vp used → dpin/dv_t computed correctly
        At[, k] <- (Rp - R0) / h
      }
      A[[t]] <- At

      # Super-diagonal block B_t = dF_t / dv_{t+1}  (only for t < T)
      if (t < horizon) {
        Bt <- matrix(0, n_eq, n_vars)
        for (k in seq_len(n_vars)) {
          vp      <- v_t1; vp[k] <- vp[k] + h
          Rp      <- .eval_period(v_t, vp, sp)
          if (t == 1L) Rp <- .pin1(Rp, v_t)   # v_t unchanged → B_1[state rows] = 0
          Bt[, k] <- (Rp - R0) / h
        }
        B[[t]] <- Bt
      }
    }
    list(A = A, B = B)
  }

  # Block back-substitution for the upper-bidiagonal block system:
  #   A_t * delta_t + B_t * delta_{t+1} = rhs_t  (t = 1..T-1)
  #   A_T * delta_T                      = rhs_T
  # With pinning, A_1[state, state] = I, so the solve naturally yields
  # delta_1[state] = rhs_1[state] = -R_state_1 = -(x_1 - x_init) → 0 if x_1 = x_init.
  .block_backsolve <- function(A, B, rhs) {
    delta <- numeric(n_vars * horizon)

    delta_t1 <- .safe_solve(A[[horizon]], rhs[.pf_idx(horizon, n_eq)])
    delta[.pf_idx(horizon, n_vars)] <- delta_t1

    for (t in seq(horizon - 1L, 1L)) {
      rhs_t   <- as.numeric(rhs[.pf_idx(t, n_eq)] - B[[t]] %*% delta_t1)
      delta_t <- .safe_solve(A[[t]], rhs_t)
      delta[.pf_idx(t, n_vars)] <- delta_t
      delta_t1 <- delta_t
    }
    delta
  }

  # ------------------------------------------------------------------
  # 7.  Newton iterations with Armijo backtracking
  #
  #  x_1 is pinned to x_init via the .pin1() override in the residuals
  #  and Jacobian.  The solve naturally gives delta_1[state] = 0 when
  #  z[state at t=1] = x_init, so x_1 remains fixed throughout without
  #  any manual zeroing.
  # ------------------------------------------------------------------
  R_cur    <- .residual_full(z)
  res_norm <- max(abs(R_cur))
  converged <- FALSE
  iters     <- 0L

  if (verbose) {
    cat("Nonlinear perfect foresight -- Newton iterations:\n")
    cat(sprintf("  start   : max|R| = %.3e\n", res_norm))
  }

  for (iter in seq_len(max_iter)) {
    if (res_norm < tol) { converged <- TRUE; break }

    jac   <- .block_jacobian(z)
    delta <- .block_backsolve(jac$A, jac$B, -R_cur)

    # Armijo backtracking: accept first step with ||R|| reduction
    alpha <- 1.0
    for (ls in seq_len(30L)) {
      z_try    <- z + alpha * delta
      R_try    <- .residual_full(z_try)
      norm_try <- max(abs(R_try))
      if (norm_try < res_norm || alpha < 1e-14) break
      alpha <- alpha * 0.5
    }

    z        <- z_try
    R_cur    <- R_try
    res_norm <- norm_try
    iters    <- iter

    if (verbose)
      cat(sprintf("  iter %3d: max|R| = %.3e  alpha = %.4f\n",
                  iter, res_norm, alpha))
  }

  if (res_norm < tol) converged <- TRUE

  if (!converged)
    warning("perfect_foresight_nonlinear did not converge in ", max_iter,
            " iterations.  max|R| = ", format(res_norm, digits = 4),
            "\nConsider increasing max_iter or tol, or checking your model.",
            call. = FALSE)

  # ------------------------------------------------------------------
  # 8.  Extract solution matrices (levels → deviations)
  # ------------------------------------------------------------------
  ctrl_lev  <- matrix(0, horizon, n_c,  dimnames = list(NULL, ctrl_names))
  state_lev <- matrix(0, horizon, n_s,  dimnames = list(NULL, state_names))

  for (t in seq_len(horizon)) {
    vt           <- z[.pf_idx(t, n_vars)]
    ctrl_lev[t,  ] <- vt[seq_len(n_c)]
    state_lev[t, ] <- vt[state_idx_in_v]
  }

  ctrl_dev  <- sweep(ctrl_lev,  2L, ss_ctrl,  "-")
  state_dev <- sweep(state_lev, 2L, ss_state, "-")

  # ------------------------------------------------------------------
  # 9.  Return  (same class as perfect_foresight())
  # ------------------------------------------------------------------
  result <- list(
    states         = state_dev,
    controls       = ctrl_dev,
    state_levels   = state_lev,
    control_levels = ctrl_lev,
    steady_state   = ss_vec,
    shock_path     = shock_path,
    initial        = x_init - ss_state,
    horizon        = horizon,
    state_names    = state_names,
    control_names  = ctrl_names,
    shock_names    = exo_names,
    H              = sol1$H,
    G              = sol1$G,
    M              = sol1$M,
    nonlinear      = TRUE,
    newton_iters   = iters,
    converged      = converged,
    max_residual   = res_norm
  )
  class(result) <- "dsge_perfect_foresight"
  result
}


# ---- Small utilities shared by the nonlinear PF solver ---------------

# Index into a stacked vector:  block t of length n (1-based)
.pf_idx <- function(t, n) ((t - 1L) * n + 1L):(t * n)

# Solve A x = b; fall back to SVD pseudo-inverse if A is near-singular
.safe_solve <- function(A, b) {
  tryCatch(solve(A, b), error = function(e) {
    sv  <- svd(A)
    tol <- max(dim(A)) * max(sv$d) * .Machine$double.eps
    sv$v %*% (ifelse(sv$d > tol, 1 / sv$d, 0) * (t(sv$u) %*% b))
  })
}


#' @export
print.dsge_perfect_foresight <- function(x, ...) {
  cat("Perfect Foresight Transition Path\n")
  cat("  Horizon:  ", x$horizon, " periods\n", sep = "")
  cat("  States:   ", length(x$state_names), " (", paste(x$state_names, collapse = ", "), ")\n", sep = "")
  cat("  Controls: ", length(x$control_names), " (", paste(x$control_names, collapse = ", "), ")\n", sep = "")
  cat("  Shocks:   ", length(x$shock_names), " (", paste(x$shock_names, collapse = ", "), ")\n", sep = "")

  # Active shocks
  active <- colSums(abs(x$shock_path)) > 0
  if (any(active)) {
    cat("  Active shocks: ", paste(x$shock_names[active], collapse = ", "), "\n", sep = "")
  }

  # Non-zero initial conditions
  nonzero_ic <- x$initial[abs(x$initial) > 1e-12]
  if (length(nonzero_ic) > 0) {
    ic_str <- paste(names(nonzero_ic), "=", round(nonzero_ic, 4), collapse = ", ")
    cat("  Initial conditions: ", ic_str, "\n", sep = "")
  }

  has_levels <- !is.null(x$state_levels)
  cat("  Levels available: ", has_levels, "\n", sep = "")

  invisible(x)
}


#' Plot Perfect Foresight Transition Paths
#'
#' Plot the deterministic transition paths from a \code{perfect_foresight}
#' result, optionally overlaying a second path for comparison (e.g.
#' linearized vs. nonlinear).
#'
#' @param x A \code{dsge_perfect_foresight} object.
#' @param vars Character vector of variable names to plot. If \code{NULL}
#'   (default), plots all variables with non-trivial paths.
#' @param type Character. One of \code{"deviation"} (default) or
#'   \code{"level"}. Controls whether to plot deviations from SS or levels.
#' @param compare Optional second \code{dsge_perfect_foresight} object to
#'   overlay on the same panels (e.g. pass the linearized
#'   \code{perfect_foresight()} result alongside a
#'   \code{perfect_foresight_nonlinear()} result to visualise the
#'   nonlinearity premium).  Only variables common to both objects are
#'   overlaid.
#' @param max_panels Integer. Maximum number of panels per plot page.
#'   Default 9.
#' @param ... Additional arguments (currently unused).
#'
#' @return No return value, called for the side effect of producing
#'   transition path plots on the active graphics device.
#'
#' @export
plot.dsge_perfect_foresight <- function(x, vars = NULL, type = "deviation",
                                        compare = NULL,
                                        max_panels = 9L, ...) {
  type <- match.arg(type, c("deviation", "level"))

  # Combine states and controls
  if (type == "level" && !is.null(x$state_levels)) {
    all_data <- cbind(x$control_levels, x$state_levels)
  } else {
    all_data <- cbind(x$controls, x$states)
    type <- "deviation"  # force if levels not available
  }

  all_names <- colnames(all_data)

  # Optional comparison object
  has_cmp <- !is.null(compare)
  if (has_cmp) {
    if (!inherits(compare, "dsge_perfect_foresight"))
      stop("'compare' must be a dsge_perfect_foresight object.",
           call. = FALSE)
    if (type == "level" && !is.null(compare$state_levels)) {
      cmp_data <- cbind(compare$control_levels, compare$state_levels)
    } else {
      cmp_data <- cbind(compare$controls, compare$states)
    }
  }

  # Select variables to plot
  if (is.null(vars)) {
    # Plot variables with non-trivial paths (max absolute deviation > 1e-10)
    max_abs <- apply(abs(all_data), 2, max)
    active <- max_abs > 1e-10
    vars <- all_names[active]
    if (length(vars) == 0) {
      message("No variables with non-trivial paths to plot.")
      return(invisible(x))
    }
  } else {
    bad <- setdiff(vars, all_names)
    if (length(bad) > 0) {
      stop("Unknown variable(s): ", paste(bad, collapse = ", "),
           "\nAvailable: ", paste(all_names, collapse = ", "))
    }
  }

  n_vars <- length(vars)
  n_pages <- ceiling(n_vars / max_panels)
  periods <- seq_len(x$horizon)

  # Label for the primary series in the legend
  primary_label <- if (isTRUE(x$nonlinear)) "Nonlinear" else "Linear"
  compare_label <- if (has_cmp) {
    if (isTRUE(compare$nonlinear)) "Nonlinear" else "Linear"
  } else NA_character_

  for (page in seq_len(n_pages)) {
    idx_start <- (page - 1L) * max_panels + 1L
    idx_end <- min(page * max_panels, n_vars)
    page_vars <- vars[idx_start:idx_end]
    n_page <- length(page_vars)

    # Layout
    nc <- min(3L, n_page)
    nr <- ceiling(n_page / nc)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    .dsge_par_grid(nr, nc)

    for (v_idx in seq_along(page_vars)) {
      v <- page_vars[v_idx]
      y <- all_data[, v]
      y_cmp <- if (has_cmp && v %in% colnames(cmp_data)) cmp_data[, v]
               else NULL

      # ylim accommodating both series
      ylim <- range(y, na.rm = TRUE)
      if (!is.null(y_cmp)) {
        ylim <- range(c(ylim, y_cmp), na.rm = TRUE)
      }
      if (type == "level") {
        ss_val <- x$steady_state[v]
        if (!is.na(ss_val)) ylim <- range(c(ylim, ss_val))
      }

      ylab <- if (type == "deviation") "Deviation from SS" else "Level"
      graphics::plot(periods, y, type = "n",
                     xlab = "Period", ylab = ylab,
                     main = v, ylim = ylim)
      .dsge_grid()
      if (type == "deviation") {
        .dsge_zero_line()
      } else {
        ss_val <- x$steady_state[v]
        if (!is.na(ss_val)) {
          graphics::abline(h = ss_val, lty = "dotted",
                           col = .DSGE_INK_SECONDARY, lwd = 1.0)
        }
      }

      # Comparison line first (so primary sits on top)
      if (!is.null(y_cmp)) {
        graphics::lines(periods, y_cmp,
                        col = .DSGE_INK_NEUTRAL, lwd = 1.2,
                        lty = "dashed")
      }

      # Primary path
      graphics::lines(periods, y,
                      col = .DSGE_INK_PRIMARY, lwd = 1.8)

      # Legend on the first panel only
      if (has_cmp && v_idx == 1L) {
        .dsge_legend("topright",
                     legend = c(primary_label, compare_label),
                     col    = c(.DSGE_INK_PRIMARY, .DSGE_INK_NEUTRAL),
                     lty    = c("solid", "dashed"),
                     lwd    = c(1.8, 1.2))
      }
    }
  }

  invisible(x)
}


#' Summary of Perfect Foresight Transition
#'
#' @param object A \code{dsge_perfect_foresight} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the `dsge_perfect_foresight` object. Called
#'   for the side effect of printing impact effects, peak deviations,
#'   and convergence diagnostics to the console.
#'
#' @export
summary.dsge_perfect_foresight <- function(object, ...) {
  cat("Perfect Foresight Transition Summary\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")

  # Impact effects (period 1)
  cat("\nImpact effects (period 1 deviations from SS):\n")
  all_t1 <- c(object$controls[1, ], object$states[1, ])
  nonzero <- all_t1[abs(all_t1) > 1e-12]
  if (length(nonzero) > 0) {
    for (nm in names(nonzero)) {
      cat(sprintf("  %-15s %10.6f\n", nm, nonzero[nm]))
    }
  } else {
    cat("  (no impact effects)\n")
  }

  # Peak effects
  cat("\nPeak absolute deviations:\n")
  all_data <- cbind(object$controls, object$states)
  for (j in seq_len(ncol(all_data))) {
    peak <- max(abs(all_data[, j]))
    if (peak > 1e-10) {
      peak_t <- which.max(abs(all_data[, j]))
      cat(sprintf("  %-15s %10.6f  (period %d)\n",
                  colnames(all_data)[j], peak, peak_t))
    }
  }

  # Convergence check
  last_dev <- all_data[nrow(all_data), ]
  max_last <- max(abs(last_dev))
  cat(sprintf("\nMax deviation at period %d: %.2e\n", object$horizon, max_last))
  if (max_last < 1e-4) {
    cat("  Converged to steady state.\n")
  } else {
    cat("  Has NOT converged to steady state (increase horizon).\n")
  }

  invisible(object)
}
