# OccBin for imported Dynare models
#
# Dynare's occbin_constraints block together with equations tagged
# [name = '...', relax = 'C'] / [name = '...', bind = 'C'] defines a model
# with regime-dependent equations. simulate_occbin() on an imported model
# solves it with the piecewise-linear algorithm of Guerrieri and Iacoviello
# (2015), as Dynare's occbin_solver does: each regime is linearised around
# the steady state of the reference (relaxed) regime, a regime sequence is
# guessed, time-varying decision rules are computed by backward recursion
# (agents anticipate how long each constraint binds), and the guess is
# updated from the bind/relax conditions until it is consistent. Shocks are
# surprises in the period they occur.

#' Split tagged model equations into base (relaxed) and binding variants
#' @noRd
dyn_occbin_split <- function(eqs, tags) {
  is_bind <- vapply(seq_along(eqs), function(i) {
    !is.null(tags[[i]]$bind) && is.null(tags[[i]]$relax)
  }, logical(1))
  if (!any(is_bind)) return(NULL)
  base_idx <- which(!is_bind)
  alt <- lapply(which(is_bind), function(i) {
    nm <- tags[[i]]$name
    partner <- base_idx[vapply(base_idx, function(j) {
      identical(tags[[j]]$name, nm) && !is.null(tags[[j]]$relax)
    }, logical(1))]
    if (is.null(nm) || length(partner) != 1L) {
      stop("OccBin: the equation tagged bind = '", tags[[i]]$bind,
           "' needs a matching equation with the same name tag and ",
           "relax = '", tags[[i]]$bind, "'.", call. = FALSE)
    }
    list(constraint = tags[[i]]$bind, equation = eqs[i],
         position = match(partner, base_idx))
  })
  list(base_idx = base_idx, alt = alt)
}

#' Parse the occbin_constraints block
#' @noRd
dyn_occbin_constraints <- function(statements) {
  cons <- list()
  current <- NULL
  for (st in statements) {
    kw <- sub("^([A-Za-z_]+).*$", "\\1", st)
    body <- trimws(sub("^[A-Za-z_]+\\s*", "", st))
    if (kw == "name") {
      current <- gsub("^['\"]|['\"]$", "", body)
      cons[[current]] <- list(name = current, bind = NULL, relax = NULL)
    } else if (kw %in% c("bind", "relax")) {
      if (is.null(current)) {
        stop("occbin_constraints: '", kw, "' before 'name'.", call. = FALSE)
      }
      cons[[current]][[kw]] <- body
    } else if (kw %in% c("error_bind", "error_relax")) {
      next
    } else {
      stop("Cannot parse occbin_constraints entry: ", st, call. = FALSE)
    }
  }
  for (c in cons) {
    if (is.null(c$bind)) {
      stop("occbin_constraints: constraint '", c$name, "' has no bind ",
           "condition.", call. = FALSE)
    }
  }
  cons
}

#' Condition string -> obc_constraint-like description for print/plot
#' @noRd
dyn_occbin_describe <- function(cond, controls) {
  m <- regmatches(cond, regexec(
    "^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*(<=|>=|<|>)\\s*(.+)$", cond))[[1]]
  if (length(m) == 4L && m[2] %in% controls) {
    type <- if (m[3] %in% c("<", "<=")) ">=" else "<="
    return(list(variable = m[2], type = type, bound = trimws(m[4])))
  }
  list(variable = controls[1], type = "", bound = cond)
}

#' Build OccBin information at import time
#'
#' @param split Result of dyn_occbin_split().
#' @param block occbin_constraints statements.
#' @param rewrite Function turning a Dynare-timed equation into dsge form.
#' @param make_model Function(eq_set) -> dsgenl_model.
#' @param all_eqs Equations of the relaxed (reference) model.
#' @noRd
dyn_occbin_build <- function(split, block, rewrite, make_model, all_eqs) {
  if (is.null(split)) {
    if (length(block) > 0L) {
      stop("occbin_constraints given but no equation is tagged with ",
           "bind = '...'.", call. = FALSE)
    }
    return(NULL)
  }
  cons <- dyn_occbin_constraints(block)
  alt_cons <- vapply(split$alt, `[[`, "", "constraint")
  unknown <- setdiff(alt_cons, names(cons))
  if (length(unknown) > 0L) {
    stop("Equation tagged bind = '", unknown[1], "' but no such constraint ",
         "in occbin_constraints.", call. = FALSE)
  }
  n_con <- length(cons)
  if (n_con > 4L) {
    stop("OccBin: at most 4 constraints are supported.", call. = FALSE)
  }
  regimes <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n_con)))
  colnames(regimes) <- names(cons)
  models <- lapply(seq_len(nrow(regimes)), function(r) {
    eq_set <- all_eqs
    for (a in split$alt) {
      if (regimes[r, a$constraint]) {
        eq_set[a$position] <- rewrite(a$equation)
      }
    }
    if (r == 1L) NULL else make_model(eq_set)
  })
  list(constraints = cons, regimes = regimes, models = models)
}

#' Piecewise-linear OccBin simulation for an imported Dynare model
#' @noRd
dyn_occbin_simulate <- function(x, shocks = NULL, horizon = 40L,
                                max_iter = 100L, lookahead = NULL,
                                params = NULL) {
  ob <- x$occbin
  if (is.null(ob)) {
    stop("The imported model has no occbin_constraints.", call. = FALSE)
  }
  horizon <- as.integer(horizon)
  if (is.null(lookahead)) lookahead <- max(horizon, 100L)
  sol <- solve_dsge(x, params = params)
  if (!isTRUE(sol$stable)) {
    stop("The reference (relaxed) regime has no stable solution.",
         call. = FALSE)
  }
  model <- sol$model
  G <- sol$G
  H <- sol$H
  ss <- sol$steady_state
  pvec <- sol$params
  controls <- rownames(G)
  states <- colnames(H)
  n_c <- length(controls)
  n_s <- length(states)
  exo <- model$variables$exo_state

  # Linearise every regime around the reference steady state
  timed_names <- c(controls, states, paste0(controls, "__f"),
                   paste0(states, "__f"))
  timed_ss <- c(ss[controls], ss[states], ss[controls], ss[states])
  names(timed_ss) <- timed_names
  pv <- assemble_params_nl(model, pvec)
  lin <- lapply(seq_len(nrow(ob$regimes)), function(r) {
    m <- if (r == 1L) model else ob$models[[r]]
    pv_r <- assemble_params_nl(m, pvec[intersect(names(pvec), m$parameters)])
    f <- function(v) {
      names(v) <- timed_names
      m$eval_fn(c(v, pv_r))
    }
    J <- numDeriv::jacobian(f, timed_ss)
    list(Jy = J[, seq_len(n_c), drop = FALSE],
         Jx = J[, n_c + seq_len(n_s), drop = FALSE],
         Jyf = J[, n_c + n_s + seq_len(n_c), drop = FALSE],
         Jxf = J[, 2L * n_c + n_s + seq_len(n_s), drop = FALSE],
         c = f(timed_ss))
  })

  # Shock paths (raw innovations, surprises in the period they occur)
  shock_path <- matrix(0, horizon, length(exo), dimnames = list(NULL, exo))
  if (is.null(shocks)) shocks <- x$shock_paths
  if (!is.null(shocks)) {
    if (is.matrix(shocks)) {
      for (nm in colnames(shocks)) {
        if (!nm %in% exo) stop("Unknown shock: ", nm, call. = FALSE)
        k <- min(nrow(shocks), horizon)
        shock_path[seq_len(k), nm] <- shocks[seq_len(k), nm]
      }
    } else {
      for (nm in names(shocks)) {
        if (!nm %in% exo) stop("Unknown shock: ", nm, call. = FALSE)
        v <- as.numeric(shocks[[nm]])
        k <- min(length(v), horizon)
        shock_path[seq_len(k), nm] <- v[seq_len(k)]
      }
    }
  }

  cons <- ob$constraints
  n_con <- length(cons)
  env <- new.env(parent = baseenv())
  for (nm in names(x$params)) assign(nm, x$params[[nm]], envir = env)
  for (nm in names(pv)) assign(nm, pv[[nm]], envir = env)
  cond_exprs <- lapply(cons, function(cn) {
    list(bind = parse(text = dyn_translate_math(cn$bind))[[1]],
         relax = if (is.null(cn$relax)) NULL
                 else parse(text = dyn_translate_math(cn$relax))[[1]])
  })
  eval_cond <- function(expr, levels) {
    e <- new.env(parent = env)
    for (nm in names(levels)) assign(nm, levels[[nm]], envir = e)
    isTRUE(as.logical(eval(expr, envir = e)))
  }
  regime_row <- function(bind_vec) {
    which(apply(ob$regimes, 1L, function(r) all(r == bind_vec)))
  }

  # Perfect-foresight path from state x0 under a regime sequence
  pf_path <- function(x0, reg) {
    T_ <- nrow(reg)
    last <- if (any(reg)) max(which(apply(reg, 1L, any))) else 0L
    P <- vector("list", T_ + 1L)
    q <- vector("list", T_ + 1L)
    Hs <- vector("list", T_)
    s <- vector("list", T_)
    P_next <- G
    q_next <- numeric(n_c)
    for (t in rev(seq_len(T_))) {
      if (t > last) {
        P[[t]] <- G
        q[[t]] <- numeric(n_c)
        Hs[[t]] <- H
        s[[t]] <- numeric(n_s)
      } else {
        L <- lin[[regime_row(reg[t, ])]]
        A <- cbind(L$Jy, L$Jyf %*% P_next + L$Jxf)
        rhs_x <- -solve(A, L$Jx)
        rhs_c <- -solve(A, L$Jyf %*% q_next + L$c)
        P[[t]] <- rhs_x[seq_len(n_c), , drop = FALSE]
        q[[t]] <- rhs_c[seq_len(n_c)]
        Hs[[t]] <- rhs_x[n_c + seq_len(n_s), , drop = FALSE]
        s[[t]] <- rhs_c[n_c + seq_len(n_s)]
      }
      P_next <- P[[t]]
      q_next <- q[[t]]
    }
    ys <- matrix(0, T_, n_c)
    xs <- matrix(0, T_ + 1L, n_s)
    xs[1, ] <- x0
    for (t in seq_len(T_)) {
      ys[t, ] <- P[[t]] %*% xs[t, ] + q[[t]]
      xs[t + 1L, ] <- Hs[[t]] %*% xs[t, ] + s[[t]]
    }
    list(y = ys, x = xs)
  }

  update_regime <- function(path, reg) {
    new <- reg
    for (t in seq_len(nrow(reg))) {
      lv <- c(ss[controls] + path$y[t, ], ss[states] + path$x[t, ])
      names(lv) <- c(controls, states)
      for (k in seq_len(n_con)) {
        if (!reg[t, k]) {
          new[t, k] <- eval_cond(cond_exprs[[k]]$bind, lv)
        } else if (is.null(cond_exprs[[k]]$relax)) {
          new[t, k] <- eval_cond(cond_exprs[[k]]$bind, lv)
        } else {
          new[t, k] <- !eval_cond(cond_exprs[[k]]$relax, lv)
        }
      }
    }
    new
  }

  simulate <- function(constrained) {
    y_out <- matrix(0, horizon, n_c, dimnames = list(NULL, controls))
    x_out <- matrix(0, horizon, n_s, dimnames = list(NULL, states))
    binding <- matrix(FALSE, horizon, n_con,
                      dimnames = list(NULL, names(cons)))
    xt <- numeric(n_s)
    reg <- matrix(FALSE, lookahead, n_con)
    iters <- 0L
    ok <- TRUE
    for (t in seq_len(horizon)) {
      xt[match(exo, states)] <- xt[match(exo, states)] + shock_path[t, ]
      if (constrained) {
        conv <- FALSE
        for (it in seq_len(max_iter)) {
          path <- pf_path(xt, reg)
          new <- update_regime(path, reg)
          if (identical(new, reg)) {
            conv <- TRUE
            break
          }
          reg <- new
        }
        iters <- max(iters, it)
        if (!conv) ok <- FALSE
        if (any(reg[lookahead, ])) {
          stop("OccBin: a constraint still binds at the end of the ",
               "look-ahead window; increase `lookahead`.", call. = FALSE)
        }
      } else {
        path <- pf_path(xt, reg)
      }
      y_out[t, ] <- path$y[1, ]
      x_out[t, ] <- path$x[1, ]
      binding[t, ] <- reg[1, ]
      xt <- path$x[2, ]
      reg <- rbind(reg[-1, , drop = FALSE], rep(FALSE, n_con))
    }
    list(y = y_out, x = x_out, binding = binding, iters = iters, ok = ok)
  }

  unc <- simulate(FALSE)
  occ <- simulate(TRUE)

  structure(
    list(
      states = occ$x, controls = occ$y,
      states_unc = unc$x, controls_unc = unc$y,
      binding = occ$binding,
      shadow_shocks = matrix(0, horizon, n_con,
                             dimnames = list(NULL, names(cons))),
      n_iter = occ$iters, converged = occ$ok,
      constraints = lapply(cons, function(cn) {
        dyn_occbin_describe(cn$bind, controls)
      }),
      steady_state = ss, horizon = horizon,
      state_names = states, control_names = controls,
      shock_names = exo, H = H, G = G, M = sol$M,
      method = "piecewise-linear (Guerrieri-Iacoviello)"
    ),
    class = "dsge_occbin"
  )
}
