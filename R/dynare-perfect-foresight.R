# Deterministic (perfect foresight) simulation of imported Dynare models
#
# Reproduces Dynare's perfect_foresight_setup / perfect_foresight_solver
# (and the older simul) on the model's equations in Dynare timing: the
# endogenous variables are stacked over periods 1..T, with initial values
# (initval, histval, or the steady state computed by steady) before period
# 1 and terminal values (endval or the steady state) after period T; the
# exogenous path comes from initval/endval and the deterministic shocks
# block. The stacked system is solved by Newton's method with a sparse
# Jacobian built from symbolic derivatives. Equations tagged mcp (Dynare's
# lmmcp option) are mixed complementarity conditions, solved as
# min(v - lower, F) = 0 (semismooth Newton).

#' Perfect-foresight simulation of an imported Dynare model
#'
#' Solves the model's deterministic path exactly as Dynare's
#' `perfect_foresight_setup` + `perfect_foresight_solver` (or `simul`) do,
#' using the file's `initval`, `endval`, `histval` and `steady` statements
#' and its deterministic `shocks` (with `periods`/`values`), and
#' complementarity conditions from `mcp` equation tags (Dynare's `lmmcp`
#' option, e.g. a zero lower bound).
#'
#' @param x A model imported with [read_dynare()].
#' @param periods Number of simulation periods. Defaults to the file's
#'   `perfect_foresight_setup(periods = )`, `simul(periods = )` or
#'   `periods` statement.
#' @param shocks Optional matrix (or data frame) of exogenous values by
#'   period (rows 1, 2, ...) with columns named after shocks, replacing the
#'   file's deterministic shocks for those shocks.
#' @param params Optional named parameter values overriding the calibration.
#' @param lmmcp Logical: treat equations tagged `mcp` as complementarity
#'   conditions (Dynare's `lmmcp` option). Defaults to `TRUE` when a
#'   `perfect_foresight_solver` command of the file uses `lmmcp`; without
#'   it Dynare ignores the tags.
#' @param tol Convergence tolerance on the largest equation residual.
#' @param max_iter Maximum number of Newton iterations.
#'
#' @return An object of class `"dsge_dynare_pf"`: a list with `path` (a
#'   matrix of the endogenous variables, one row per period from the
#'   initial condition(s) to the terminal one(s), as Dynare's
#'   `oo_.endo_simul` transposed), `exo` (the exogenous path, as
#'   `oo_.exo_simul`), `periods`, `initial` and `terminal` (steady states
#'   or initial/terminal values), `converged`, `iterations` and
#'   `max_residual`.
#'
#' @examples
#' m <- read_dynare(text = "
#'   var c k;
#'   varexo a;
#'   parameters alpha beta delta;
#'   alpha = 0.33; beta = 0.99; delta = 0.025;
#'   model;
#'     1/c = beta/c(+1) * (alpha * exp(a(+1)) * k^(alpha - 1) + 1 - delta);
#'     k = exp(a) * k(-1)^alpha + (1 - delta) * k(-1) - c;
#'   end;
#'   initval; a = 0; k = 30; c = 2.3; end;
#'   steady;
#'   endval; a = 0.1; end;
#'   steady;
#'   perfect_foresight_setup(periods = 100);
#'   perfect_foresight_solver;
#' ")
#' pf <- simulate_perfect_foresight(m)
#' head(pf$path)
#' @export
simulate_perfect_foresight <- function(x, periods = NULL, shocks = NULL,
                                       params = NULL, lmmcp = NULL,
                                       tol = 1e-10, max_iter = 100L) {
  if (!inherits(x, "dsge_dynare")) {
    stop("`x` must be a model imported with read_dynare(); for other ",
         "models use perfect_foresight() or perfect_foresight_nonlinear().",
         call. = FALSE)
  }
  spec <- x$perfect_foresight
  if (is.null(spec)) stop("The model has no perfect-foresight information.",
                          call. = FALSE)
  if (is.null(periods)) periods <- spec$periods
  if (is.null(periods) || is.na(periods)) {
    stop("Give `periods` (the file sets none).", call. = FALSE)
  }
  periods <- as.integer(periods)
  prm <- x$params
  if (!is.null(params)) prm[names(params)] <- params
  if (is.null(lmmcp)) lmmcp <- isTRUE(spec$lmmcp)
  sys <- dyn_pf_system(spec)
  if (!lmmcp) sys$mcp <- vector("list", length(sys$exprs))
  endo <- spec$endo
  exo <- spec$exo
  n <- length(endo)
  env0 <- new.env(parent = spec$cal_env)
  for (nm in names(prm)) assign(nm, prm[[nm]], envir = env0)

  # --- initial and terminal conditions (Dynare's make_y_ / make_ex_) -------
  eval_block <- function(statements, base_endo, base_exo) {
    vals <- c(base_endo, base_exo)
    env <- new.env(parent = env0)
    for (nm in names(vals)) assign(nm, vals[[nm]], envir = env)
    for (st in statements) {
      if (!grepl("=", st)) next
      nm <- trimws(sub("=.*$", "", st))
      v <- dyn_eval(sub("^[^=]*=", "", st), env, paste0("value of '", nm, "'"))
      assign(nm, v, envir = env)
      vals[nm] <- v
    }
    list(endo = vals[endo], exo = vals[exo])
  }
  zero_endo <- stats::setNames(numeric(n), endo)
  zero_exo <- stats::setNames(numeric(length(exo)), exo)
  ini <- eval_block(spec$initval, zero_endo, zero_exo)
  if (spec$steady_after_init) {
    ini$endo <- dyn_pf_steady(sys, ini$endo, ini$exo, env0, spec)
  }
  if (length(spec$endval) > 0L) {
    ter <- eval_block(spec$endval, ini$endo, ini$exo)
    if (spec$steady_after_end) {
      ter$endo <- dyn_pf_steady(sys, ter$endo, ter$exo, env0, spec)
    }
    ex_before <- ini$exo
  } else {
    ter <- ini
    ex_before <- ini$exo
  }
  lagmax <- sys$maxlag
  leadmax <- sys$maxlead
  nper <- lagmax + periods + leadmax
  Y <- matrix(rep(ter$endo, nper), n, nper, dimnames = list(endo, NULL))
  if (lagmax > 0L) Y[, seq_len(lagmax)] <- ini$endo
  # histval: x(0) is the last initial period, x(-1) the one before, ...
  for (st in spec$histval) {
    lhs <- trimws(sub("=.*$", "", st))
    nm <- sub("\\s*\\(.*$", "", lhs)
    lag <- if (grepl("\\(", lhs)) as.integer(sub("^.*\\(\\s*([-+]?[0-9]+)\\s*\\).*$", "\\1", lhs)) else 0L
    col <- lagmax + lag
    if (nm %in% endo && col >= 1L) {
      Y[nm, col] <- dyn_eval(sub("^[^=]*=", "", st), env0, "histval")
    }
  }
  # oo_.endo_simul(rows, 1) = value (rows ':' or an index vector, in
  # declaration order with Ramsey multipliers last, as in Dynare)
  for (es in spec$endo_simul_init) {
    rows <- if (trimws(es$rows) == ":") seq_len(n) else
      as.integer(dyn_eval_matlab(es$rows, env0))
    val <- dyn_eval_matlab(es$value, env0)
    if (lagmax >= 1L) Y[rows, 1L] <- val
  }
  X <- matrix(0, nper, length(exo), dimnames = list(NULL, exo))
  if (length(exo)) {
    X[] <- rep(ter$exo, each = nper)
    if (lagmax > 0L) X[seq_len(lagmax), ] <- rep(ex_before, each = lagmax)
    det <- spec$det
    if (!is.null(det)) {
      for (nm in intersect(colnames(det), exo)) {
        rows <- seq_len(min(nrow(det), periods))
        v <- det[rows, nm]
        set <- spec$det_set[[nm]]
        idx <- if (is.null(set)) rows else intersect(set, rows)
        X[lagmax + idx, nm] <- v[idx]
      }
    }
    if (!is.null(shocks)) {
      for (nm in intersect(colnames(shocks), exo)) {
        r <- seq_len(min(nrow(shocks), periods))
        X[lagmax + r, nm] <- shocks[r, nm]
      }
    }
  }

  # STEADY_STATE(x): Dynare's oo_.steady_state (the terminal one)
  for (nm in names(spec$ss_links)) {
    assign(nm, ter$endo[[spec$ss_links[[nm]]]], envir = env0)
  }

  # Dynare starts from the terminal values; if Newton fails from there
  # (e.g. purely backward-looking dynamics far from the end point), try the
  # initial values and a linear path between the two
  res <- dyn_pf_newton(sys, Y, X, periods, env0, tol, max_iter)
  if (!res$converged) {
    cols <- lagmax + seq_len(periods)
    alt <- list(ini$endo, NULL)
    for (a in seq_along(alt)) {
      Y2 <- Y
      Y2[, cols] <- if (!is.null(alt[[a]])) alt[[a]] else {
        w <- matrix(seq_len(periods) / (periods + 1), n, periods, byrow = TRUE)
        (1 - w) * ini$endo + w * ter$endo
      }
      r2 <- dyn_pf_newton(sys, Y2, X, periods, env0, tol, max_iter)
      if (r2$converged || r2$max_residual < res$max_residual) res <- r2
      if (res$converged) break
    }
  }
  path <- t(res$Y)
  rownames(path) <- as.character(seq_len(nper) - lagmax)
  exo_path <- X
  rownames(exo_path) <- rownames(path)
  structure(list(
    path = path[, spec$variables, drop = FALSE],
    path_all = path,
    exo = exo_path,
    periods = periods,
    initial = ini$endo[spec$variables],
    terminal = ter$endo[spec$variables],
    converged = res$converged,
    iterations = res$iterations,
    max_residual = res$max_residual
  ), class = "dsge_dynare_pf")
}

#' @export
print.dsge_dynare_pf <- function(x, ...) {
  cat("Perfect-foresight simulation (", x$periods, " periods)\n", sep = "")
  cat(if (x$converged) "  Converged" else "  NOT converged", " after ",
      x$iterations, " Newton iterations; max residual ",
      format(x$max_residual, digits = 3), "\n", sep = "")
  cat("  Variables:", paste(colnames(x$path), collapse = ", "), "\n")
  invisible(x)
}

#' @export
plot.dsge_dynare_pf <- function(x, vars = NULL, ...) {
  if (is.null(vars)) vars <- utils::head(colnames(x$path), 9L)
  n <- length(vars)
  nc <- min(3L, n)
  op <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(3, 3, 2, 1))
  on.exit(graphics::par(op))
  t <- as.numeric(rownames(x$path))
  for (v in vars) {
    graphics::plot(t, x$path[, v], type = "l", main = v, xlab = "",
                   ylab = "", col = "#1f3a5f", lwd = 1.5)
    graphics::abline(h = x$terminal[[v]], lty = 3, col = "grey50")
  }
  invisible(x)
}

#' Compile the Dynare-timed equations for stacked evaluation
#' @noRd
dyn_pf_system <- function(spec) {
  endo <- spec$endo
  exo <- spec$exo
  timed <- list()
  exprs <- vector("list", length(spec$equations))
  refs <- vector("list", length(spec$equations))
  for (k in seq_along(spec$equations)) {
    eq <- spec$equations[[k]]
    if (grepl("=", eq)) {
      lhs <- sub("=.*$", "", eq)
      rhs <- sub("^[^=]*=", "", eq)
      eq <- paste0("(", lhs, ") - (", rhs, ")")
    }
    used <- list()
    s <- dyn_rewrite_ids(eq, function(name, lag, has_index) {
      if (!name %in% c(endo, exo)) return(NULL)
      tn <- if (lag == 0L) paste0(name, "__t0") else
        paste0(name, if (lag < 0L) "__m" else "__p", abs(lag))
      used[[tn]] <<- list(name = name, lag = lag, endo = name %in% endo)
      tn
    })
    s <- dyn_translate_math(s)
    # evaluated on vectors over periods: element-wise max/min
    s <- gsub("(?<![A-Za-z0-9_.])max\\s*\\(", "pmax(", s, perl = TRUE)
    s <- gsub("(?<![A-Za-z0-9_.])min\\s*\\(", "pmin(", s, perl = TRUE)
    exprs[[k]] <- str2lang(s)
    refs[[k]] <- used
  }
  lags <- unlist(lapply(refs, function(r) {
    vapply(r, function(u) if (u$endo) u$lag else 0L, 0L)
  }))
  mcp <- lapply(spec$tags, function(tg) {
    if (is.null(tg) || is.null(tg$mcp)) return(NULL)
    m <- regmatches(tg$mcp, regexec(
      "^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*([<>]=?)\\s*(.+?)\\s*$", tg$mcp))[[1]]
    if (length(m) != 4L) return(NULL)
    list(var = m[2], lower = startsWith(m[3], ">"), bound = m[4])
  })
  length(mcp) <- length(exprs)
  list(exprs = exprs, refs = refs, endo = endo, exo = exo,
       maxlag = max(0L, -min(c(0L, lags))), maxlead = max(0L, lags), mcp = mcp,
       derivs = new.env(parent = emptyenv()))
}

#' Residuals (and Jacobian) of the stacked system
#' @noRd
dyn_pf_eval <- function(sys, Y, X, periods, env0, jacobian = TRUE) {
  n <- length(sys$endo)
  lagmax <- sys$maxlag
  cols <- lagmax + seq_len(periods)
  env <- new.env(parent = env0)
  for (k in seq_along(sys$exprs)) {
    for (tn in names(sys$refs[[k]])) {
      if (exists(tn, envir = env, inherits = FALSE)) next
      u <- sys$refs[[k]][[tn]]
      v <- if (u$endo) Y[u$name, cols + u$lag] else X[cols + u$lag, u$name]
      assign(tn, v, envir = env)
    }
  }
  n_eq <- length(sys$exprs)
  R <- matrix(0, periods, n_eq)
  ii <- integer(0)
  jj <- integer(0)
  xx <- numeric(0)
  for (k in seq_len(n_eq)) {
    Fk <- rep_len(as.numeric(eval(sys$exprs[[k]], env)), periods)
    rows <- (seq_len(periods) - 1L) * n_eq + k
    mc <- sys$mcp[[k]]
    use_F <- rep(TRUE, periods)
    if (!is.null(mc)) {
      b <- dyn_eval(mc$bound, env0, "mcp bound")
      vv <- Y[mc$var, cols]
      d <- if (mc$lower) vv - b else b - vv
      Fs <- if (mc$lower) Fk else -Fk
      use_F <- Fs <= d
      Fk <- ifelse(use_F, Fs, d)
      if (!mc$lower) Fk <- -Fk
    }
    R[, k] <- Fk
    if (!jacobian) next
    if (!is.null(mc) && any(!use_F)) {
      tt <- which(!use_F)
      ii <- c(ii, rows[tt])
      jj <- c(jj, (tt - 1L) * n + match(mc$var, sys$endo))
      xx <- c(xx, rep(1, length(tt)))
    }
    for (tn in names(sys$refs[[k]])) {
      u <- sys$refs[[k]][[tn]]
      if (!u$endo) next
      g <- dyn_pf_deriv(sys, k, tn, env, periods)
      t_idx <- seq_len(periods)
      ok <- use_F & (t_idx + u$lag >= 1L) & (t_idx + u$lag <= periods) & g != 0
      if (!any(ok)) next
      ii <- c(ii, rows[ok])
      jj <- c(jj, (t_idx[ok] + u$lag - 1L) * n + match(u$name, sys$endo))
      xx <- c(xx, g[ok])
    }
  }
  list(R = as.vector(t(R)), i = ii, j = jj, x = xx)
}

#' Derivative of equation k with respect to a timed variable (vector over t)
#' @noRd
dyn_pf_deriv <- function(sys, k, tn, env, periods) {
  key <- paste(k, tn)
  d <- sys$derivs[[key]]
  if (is.null(d)) {
    d <- tryCatch(stats::D(sys$exprs[[k]], tn), error = function(e) NA)
    assign(key, d, envir = sys$derivs)
  }
  if (!identical(d, NA)) {
    return(rep_len(as.numeric(eval(d, env)), periods))
  }
  # functions stats::D cannot differentiate (max, min, abs, ...)
  x0 <- get(tn, envir = env)
  h <- 1e-6 * pmax(1, abs(x0))
  assign(tn, x0 + h, envir = env)
  fp <- rep_len(as.numeric(eval(sys$exprs[[k]], env)), periods)
  assign(tn, x0 - h, envir = env)
  fm <- rep_len(as.numeric(eval(sys$exprs[[k]], env)), periods)
  assign(tn, x0, envir = env)
  (fp - fm) / (2 * h)
}

#' Newton's method on the stacked system
#' @noRd
dyn_pf_newton <- function(sys, Y, X, periods, env0, tol, max_iter) {
  # plain Newton steps first (as Dynare's solver), damped ones if that fails
  r <- dyn_pf_newton1(sys, Y, X, periods, env0, tol, max_iter, damped = FALSE)
  if (r$converged) return(r)
  r2 <- dyn_pf_newton1(sys, Y, X, periods, env0, tol, max_iter, damped = TRUE)
  if (r2$converged || !is.finite(r$max_residual) ||
      r2$max_residual < r$max_residual) r2 else r
}

#' @noRd
dyn_pf_newton1 <- function(sys, Y, X, periods, env0, tol, max_iter, damped) {
  n <- length(sys$endo)
  cols <- sys$maxlag + seq_len(periods)
  N <- n * periods
  have_matrix <- requireNamespace("Matrix", quietly = TRUE)
  res <- dyn_pf_eval(sys, Y, X, periods, env0)
  err <- max(abs(res$R))
  it <- 0L
  nonmono <- 0L
  while (it < max_iter && (!is.finite(err) || err > tol)) {
    it <- it + 1L
    step <- dyn_pf_step(res, N, have_matrix, sys, Y, X, periods, env0, cols)
    if (is.null(step) || !all(is.finite(step))) break
    # line search on the sum of squares; with complementarity conditions
    # and kinks (max, min) a full step may not decrease it, so after a
    # failed search the full step is taken anyway (non-monotone), a limited
    # number of times
    if (!damped) {
      Yn <- Y
      Yn[, cols] <- Y[, cols] + matrix(step, n, periods)
      rn <- tryCatch(dyn_pf_eval(sys, Yn, X, periods, env0),
                     error = function(e) NULL)
      if (is.null(rn) || !all(is.finite(rn$R))) break
      Y <- Yn
      res <- rn
      err <- max(abs(res$R))
      next
    }
    ssr <- sum(res$R^2)
    lam <- 1
    improved <- FALSE
    for (ls in 1:30) {
      Yn <- Y
      Yn[, cols] <- Y[, cols] + lam * matrix(step, n, periods)
      rn <- tryCatch(dyn_pf_eval(sys, Yn, X, periods, env0, jacobian = FALSE),
                     error = function(e) NULL)
      if (!is.null(rn) && all(is.finite(rn$R)) && sum(rn$R^2) < ssr) {
        improved <- TRUE
        break
      }
      lam <- lam / 2
    }
    if (!improved) {
      nonmono <- nonmono + 1L
      if (nonmono > 20L) break
      Yn <- Y
      Yn[, cols] <- Y[, cols] + matrix(step, n, periods)
      rn <- tryCatch(dyn_pf_eval(sys, Yn, X, periods, env0, jacobian = FALSE),
                     error = function(e) NULL)
      if (is.null(rn) || !all(is.finite(rn$R))) break
    }
    Y <- Yn
    res <- dyn_pf_eval(sys, Y, X, periods, env0)
    err <- max(abs(res$R))
  }
  list(Y = Y, converged = is.finite(err) && err <= tol, iterations = it,
       max_residual = err)
}

#' Steady state of the Dynare equations at given exogenous values
#' @noRd
dyn_pf_steady <- function(sys, guess, exo_vals, env0, spec) {
  n <- length(sys$endo)
  X <- matrix(exo_vals, 1L, length(exo_vals), dimnames = list(NULL, sys$exo))
  # a one-period system with every lead and lag at the same value
  f <- function(y) {
    env <- new.env(parent = env0)
    for (k in seq_along(sys$exprs)) {
      for (tn in names(sys$refs[[k]])) {
        u <- sys$refs[[k]][[tn]]
        assign(tn, if (u$endo) y[[u$name]] else X[1L, u$name], envir = env)
      }
    }
    for (nm in names(spec$ss_links)) assign(nm, y[[spec$ss_links[[nm]]]],
                                           envir = env)
    vapply(sys$exprs, function(e) as.numeric(eval(e, env))[1L], 0)
  }
  # As Dynare's `steady`, start from the steady_state_model block when the
  # file has one (evaluated at the current exogenous values); Newton then
  # only polishes it.
  y <- dyn_pf_ss_model(spec, guess, stats::setNames(exo_vals, sys$exo), env0)
  for (it in 1:100) {
    r <- f(y)
    if (!all(is.finite(r))) {
      stop("The steady state for the perfect-foresight simulation cannot be ",
           "computed from the starting values (the model equations are not ",
           "finite there). Give starting values in `initval` or a ",
           "`steady_state_model` block.", call. = FALSE)
    }
    if (max(abs(r)) < 1e-12) break
    J <- numDeriv::jacobian(function(z) f(stats::setNames(z, names(y))), y)
    step <- dyn_lstsq(J, -r)
    lam <- 1
    for (ls in 1:30) {
      yn <- y + lam * step
      rn <- tryCatch(f(yn), error = function(e) rep(NaN, length(r)))
      if (all(is.finite(rn)) && sum(rn^2) < sum(r^2)) break
      lam <- lam / 2
    }
    y <- yn
  }
  if (max(abs(f(y))) > 1e-8) {
    warning("Steady state for the perfect-foresight simulation did not ",
            "converge (max residual ", format(max(abs(f(y))), digits = 3),
            ").", call. = FALSE)
  }
  y
}

#' Steady state from the file's steady_state_model block, used as the
#' starting point of dyn_pf_steady(); variables the block does not set keep
#' their value in `guess`
#' @noRd
dyn_pf_ss_model <- function(spec, guess, exo_vals, env0) {
  stmts <- spec$steady_state_model
  if (length(stmts) == 0L) return(guess)
  env <- new.env(parent = env0)
  for (nm in names(exo_vals)) assign(nm, exo_vals[[nm]], envir = env)
  for (nm in names(guess)) assign(nm, guess[[nm]], envir = env)
  ok <- tryCatch({
    for (st in stmts) {
      if (!grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*=[^=]", st)) next
      nm <- trimws(sub("=.*$", "", st))
      assign(nm, eval(parse(text = dyn_translate_math(sub("^[^=]*=", "", st))),
                      envir = env), envir = env)
    }
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(guess)
  for (nm in names(guess)) {
    v <- get(nm, envir = env)
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) guess[[nm]] <- v
  }
  guess
}

#' Newton step, or an adaptive Levenberg-Marquardt step when the Jacobian is
#' (near) singular and the Newton step does not solve the linear system
#' @noRd
dyn_pf_step <- function(res, N, have_matrix, sys, Y, X, periods, env0, cols) {
  n <- length(sys$endo)
  if (have_matrix) {
    J <- Matrix::sparseMatrix(i = res$i, j = res$j, x = res$x, dims = c(N, N))
  } else {
    J <- matrix(0, N, N)
    idx <- (res$j - 1L) * N + res$i
    sums <- tapply(res$x, idx, sum)
    J[as.numeric(names(sums))] <- sums
  }
  slv <- if (have_matrix) Matrix::solve else solve
  st <- tryCatch(as.numeric(slv(J, -res$R)), error = function(e) NULL)
  if (!is.null(st) && all(is.finite(st))) {
    lin <- as.numeric(J %*% st) + res$R
    if (max(abs(lin)) <= 1e-6 * max(1, max(abs(res$R)))) return(st)
  }
  # Levenberg-Marquardt: increase the damping until the residual falls
  JtJ <- if (have_matrix) Matrix::crossprod(J) else crossprod(J)
  g <- -as.numeric(if (have_matrix) Matrix::crossprod(J, res$R) else
    crossprod(J, res$R))
  d <- if (have_matrix) Matrix::diag(JtJ) else diag(JtJ)
  scale <- max(1, max(abs(d)))
  ssr <- sum(res$R^2)
  I_N <- if (have_matrix) Matrix::Diagonal(N) else diag(N)
  for (mu in scale * 10^seq(-10, 4)) {
    st <- tryCatch(as.numeric(slv(JtJ + mu * I_N, g)),
                   error = function(e) NULL)
    if (is.null(st) || !all(is.finite(st))) next
    Yn <- Y
    Yn[, cols] <- Y[, cols] + matrix(st, n, periods)
    rn <- tryCatch(dyn_pf_eval(sys, Yn, X, periods, env0, jacobian = FALSE),
                   error = function(e) NULL)
    if (!is.null(rn) && all(is.finite(rn$R)) && sum(rn$R^2) < ssr) return(st)
  }
  NULL
}
