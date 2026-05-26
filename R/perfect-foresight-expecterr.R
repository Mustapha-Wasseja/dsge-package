# ==========================================================================
# Perfect foresight with expectation errors (Dynare:
#   `perfect_foresight_with_expectation_errors_solver`)
# ==========================================================================
#
# Standard perfect_foresight() solves once: agents at t = 1 see the entire
# future shock path and optimise accordingly.
#
# In `perfect_foresight_expect_err()` agents instead believe at every
# period t = k that no further shocks will arrive (or that only the
# announced future shocks will arrive), then nature actually delivers the
# shock path -- the realised path is *not* equal to what agents expected
# at t = 1.  This captures "surprise" shocks at later dates that agents
# did not initially anticipate.
#
# Algorithm:
#   For each k = 1 .. horizon:
#     - Construct the agent's expectation of the remaining shock path
#       (a `subjective` argument; defaults to "no further shocks
#       beyond what has already realised").
#     - Solve a perfect foresight problem with the residual horizon
#       and that expected shock path, starting from the actual state at
#       k.
#     - Apply the actual one-period shock at k (potentially different
#       from what agents expected).  This gives the realised state at
#       k + 1.
#   Final realised paths are returned.
# ==========================================================================


#' Perfect Foresight Simulation with Expectation Errors
#'
#' Generalises \code{\link{perfect_foresight}} by allowing the realised
#' shock path to differ from the path agents anticipate at each point in
#' time.  In standard perfect foresight, agents at \eqn{t = 1} see every
#' future shock; in the with-expectation-errors variant, agents form
#' subjective expectations at each period \eqn{k}, solve the residual
#' perfect-foresight problem, then nature delivers the actual one-period
#' shock (which may be a surprise).
#'
#' @param x A \code{dsge_solution} object.
#' @param actual_shocks Named list (or matrix) describing the
#'   \emph{realised} shock path -- same format as the \code{shocks}
#'   argument of \code{\link{perfect_foresight}}.
#' @param expected_shocks Optional named list (or matrix) describing
#'   what agents \emph{expect} at \eqn{t = 1}.  Default \code{NULL}
#'   meaning agents expect no further shocks beyond those already
#'   realised.  When agents expect the same path that actually
#'   materialises, this function reproduces \code{perfect_foresight()}.
#' @param initial Optional named numeric vector of initial state
#'   deviations.
#' @param horizon Integer.  Number of periods.  Default 40.
#'
#' @return An object of class \code{c("dsge_perfect_foresight_expecterr",
#'   "dsge_perfect_foresight")} containing the same fields as
#'   \code{\link{perfect_foresight}} plus an extra element
#'   \code{expectation_paths} -- a list of per-period subjective
#'   forecast paths (one matrix per starting period) so users can
#'   inspect how agent expectations evolved.
#'
#' @details
#' This is the analogue of Dynare's
#' \code{perfect_foresight_with_expectation_errors_solver} command.  A
#' typical use case: study how the economy reacts to a sequence of
#' news/MIT shocks that arrive unexpectedly even though each shock,
#' once it lands, is treated as fully credible going forward.
#'
#' @examples
#' \donttest{
#' m <- dsge_model(
#'   obs(y ~ beta * lead(y) + 0.1 * x),
#'   state(x ~ 0.9 * x),
#'   fixed = list(beta = 0.99))
#' sol <- solve_dsge(m, params = c(), shock_sd = c(x = 1))
#' # Agents expect no shocks; nature delivers a one-time shock at t = 5
#' pf <- perfect_foresight_expect_err(sol,
#'   actual_shocks = list(x = c(0, 0, 0, 0, 1)),
#'   horizon = 30)
#' plot(pf)
#' }
#'
#' @seealso \code{\link{perfect_foresight}} (no expectation errors),
#'   \code{\link{perfect_foresight_nonlinear}}.
#'
#' @export
perfect_foresight_expect_err <- function(x,
                                         actual_shocks,
                                         expected_shocks = NULL,
                                         initial = NULL,
                                         horizon = 40L) {

  if (!inherits(x, "dsge_solution"))
    stop("`x` must be a dsge_solution object.", call. = FALSE)
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("horizon must be >= 1.", call. = FALSE)

  H <- x$H;  G <- x$G;  M <- x$M
  shock_names <- colnames(M)
  state_names <- colnames(H)
  control_names <- rownames(G)
  n_s <- nrow(H);  n_c <- nrow(G);  n_e <- ncol(M)

  # ---- 1. Build actual and expected shock matrices ----
  actual_path <- .pf_shock_to_matrix(actual_shocks, shock_names, horizon)
  if (is.null(expected_shocks)) {
    expected_path <- matrix(0, horizon, n_e)
    colnames(expected_path) <- shock_names
  } else {
    expected_path <- .pf_shock_to_matrix(expected_shocks, shock_names, horizon)
  }

  # ---- 2. Initial state ----
  x_state <- numeric(n_s)
  names(x_state) <- state_names
  if (!is.null(initial)) {
    bad <- setdiff(names(initial), state_names)
    if (length(bad) > 0)
      stop("Unknown state(s) in initial: ", paste(bad, collapse = ", "),
           call. = FALSE)
    x_state[names(initial)] <- initial
  }

  # ---- 3. Forward loop with revised expectations ----
  states_realised   <- matrix(0, horizon, n_s,
                              dimnames = list(NULL, state_names))
  controls_realised <- matrix(0, horizon, n_c,
                              dimnames = list(NULL, control_names))
  expectation_paths <- vector("list", horizon)

  for (k in seq_len(horizon)) {
    # Subjective expectation from period k onward: combination of the
    # initially-expected path and what has actually arrived since.
    if (k == 1L) {
      subj_path <- expected_path
    } else {
      subj_path <- expected_path
      # The agent observes shocks 1..k-1 as they realised
      if (k > 1L) subj_path[seq_len(k - 1L), ] <- actual_path[seq_len(k - 1L), ]
    }
    # Forecasted state path from the agent's perspective starting at k
    # with state x_state and remaining shocks `subj_path[k:horizon, ]`.
    # We only need k-th step for the realisation, so simulate one step.
    # But we also store the full subjective path for inspection.
    h_remaining <- horizon - k + 1L
    subj_state <- matrix(0, h_remaining, n_s)
    subj_ctrl  <- matrix(0, h_remaining, n_c)
    x_cur <- x_state
    for (j in seq_len(h_remaining)) {
      shk <- subj_path[k + j - 1L, ]
      x_cur <- as.numeric(H %*% x_cur + M %*% shk)
      subj_state[j, ] <- x_cur
      subj_ctrl[j, ]  <- as.numeric(G %*% x_cur)
    }
    expectation_paths[[k]] <- list(state = subj_state, control = subj_ctrl)

    # Realised one-step update using the ACTUAL shock at period k
    shk_actual <- actual_path[k, ]
    x_state <- as.numeric(H %*% x_state + M %*% shk_actual)
    names(x_state) <- state_names
    states_realised[k, ]   <- x_state
    controls_realised[k, ] <- as.numeric(G %*% x_state)
  }

  # Levels (deviations + steady state if available)
  ss <- if (!is.null(x$steady_state)) x$steady_state else NULL

  structure(
    list(
      states          = states_realised,
      controls        = controls_realised,
      shock_path      = actual_path,
      expected_path   = expected_path,
      expectation_paths = expectation_paths,
      steady_state    = ss,
      initial         = initial,
      horizon         = horizon,
      state_names     = state_names,
      control_names   = control_names,
      shock_names     = shock_names,
      H = H, G = G, M = M
    ),
    class = c("dsge_perfect_foresight_expecterr",
              "dsge_perfect_foresight"))
}


#' Convert a shock specification (list or matrix) to a (horizon x n_shocks)
#' matrix, with unspecified periods set to zero.
#' @noRd
.pf_shock_to_matrix <- function(shocks, shock_names, horizon) {
  n_e <- length(shock_names)
  out <- matrix(0, horizon, n_e, dimnames = list(NULL, shock_names))
  if (is.null(shocks)) return(out)
  if (is.matrix(shocks)) {
    if (ncol(shocks) != n_e)
      stop("shocks matrix must have ", n_e, " columns.", call. = FALSE)
    nr <- min(nrow(shocks), horizon)
    out[seq_len(nr), ] <- shocks[seq_len(nr), , drop = FALSE]
    return(out)
  }
  if (is.list(shocks)) {
    for (nm in names(shocks)) {
      idx <- match(nm, shock_names)
      if (is.na(idx))
        stop("Unknown shock '", nm, "'. Available: ",
             paste(shock_names, collapse = ", "), call. = FALSE)
      v <- as.numeric(shocks[[nm]])
      nr <- min(length(v), horizon)
      out[seq_len(nr), idx] <- v[seq_len(nr)]
    }
    return(out)
  }
  stop("`shocks` must be a named list or a matrix.", call. = FALSE)
}
