# Optimal policy blocks of imported Dynare models
#
# * ramsey_model / ramsey_policy: the planner's first-order conditions are
#   derived symbolically (stats::D) from the Lagrangian
#     E_0 sum_t beta^t [ U_t + sum_i MULT_i,t f_i,t ],
#   exactly as Dynare does, and appended to the model; the augmented model
#   is an ordinary rational-expectations model solved by solve_dsge().
# * discretionary_policy: for linear models with a quadratic objective, the
#   time-consistent policy is computed with the Dennis (2007) / Soederlind
#   (1999) fixed-point algorithm and the model is closed with the resulting
#   rule for the instrument(s).
# * osr: the rule parameters, bounds and loss weights are collected so that
#   osr() can be called directly on the imported model.

#' Collect optimal-policy declarations from parsed statements
#' @noRd
dyn_policy_spec <- function(p) {
  cmds <- p$commands
  names_c <- vapply(cmds, `[[`, "", "name")
  spec <- list(type = NULL, objective = NULL, discount = NULL,
               instruments = NULL, osr_params = NULL)

  obj <- cmds[names_c == "planner_objective"]
  if (length(obj) > 0L) {
    st <- obj[[length(obj)]]$statement
    spec$objective <- trimws(sub("^planner_objective\\s*", "", st))
    if (grepl("^\\(.*\\)$", spec$objective)) {
      spec$objective <- dyn_paren_content(spec$objective, 1L)
    }
  }
  opt_value <- function(opts, key) {
    m <- regmatches(opts, regexec(paste0(key, "\\s*=\\s*(\\([^)]*\\)|[^,]+)"),
                                  opts))[[1]]
    if (length(m) == 2L) trimws(gsub("^\\(|\\)$", "", m[2])) else NULL
  }
  for (kw in c("ramsey_model", "ramsey_policy", "discretionary_policy")) {
    hit <- cmds[names_c == kw]
    if (length(hit) == 0L) next
    opts <- hit[[1]]$options
    spec$type <- if (kw == "discretionary_policy") "discretion" else "ramsey"
    spec$discount <- opt_value(opts, "planner_discount")
    instr <- opt_value(opts, "instruments")
    if (!is.null(instr)) {
      spec$instruments <- strsplit(instr, "[[:space:],]+")[[1]]
    }
  }
  osr_decl <- cmds[names_c == "osr_params"]
  if (length(osr_decl) > 0L) {
    spec$osr_params <- unlist(lapply(osr_decl, function(cm) {
      strsplit(trimws(sub("^osr_params\\s*", "", cm$statement)),
               "[[:space:],]+")[[1]]
    }))
    if (is.null(spec$type)) spec$type <- "osr"
  }
  spec
}

#' Policy information returned with the imported model
#' @noRd
dyn_policy_info <- function(p, cal_env, endo, params, spec = NULL,
                            extra = list()) {
  if (is.null(spec)) spec <- dyn_policy_spec(p)
  notes <- character(0)
  if (is.null(spec$type)) return(list(policy = NULL, notes = notes))
  if (spec$type == "osr") {
    w <- dyn_optim_weights(p$blocks$optim_weights, endo, cal_env)
    bounds <- dyn_osr_bounds(p$blocks$osr_params_bounds, spec$osr_params,
                             cal_env)
    missing_p <- setdiff(spec$osr_params, names(params))
    if (length(missing_p) > 0L) {
      stop("osr_params without a calibrated starting value: ",
           paste(missing_p, collapse = ", "), call. = FALSE)
    }
    notes <- c(notes, "OSR problem imported; run osr() on this object.")
    return(list(policy = list(type = "osr",
                              osr_params = params[spec$osr_params],
                              weights = w, lower = bounds$lower,
                              upper = bounds$upper),
                notes = notes))
  }
  list(policy = c(list(type = spec$type, objective = spec$objective,
                       discount = spec$discount,
                       instruments = spec$instruments), extra),
       notes = notes)
}

#' optim_weights block -> symmetric weight matrix on endogenous variables
#' @noRd
dyn_optim_weights <- function(statements, endo, cal_env) {
  W <- matrix(0, length(endo), length(endo), dimnames = list(endo, endo))
  for (st in statements) {
    m <- regmatches(st, regexec(
      "^([A-Za-z_][A-Za-z0-9_]*)\\s*(?:,\\s*([A-Za-z_][A-Za-z0-9_]*))?\\s+(.+)$",
      st, perl = TRUE))[[1]]
    if (length(m) != 4L) {
      stop("Cannot parse optim_weights entry: ", st, call. = FALSE)
    }
    a <- m[2]
    b <- if (nzchar(m[3])) m[3] else a
    if (!all(c(a, b) %in% endo)) {
      stop("optim_weights refers to an undeclared variable: ", st,
           call. = FALSE)
    }
    val <- dyn_eval(m[4], cal_env, "optim_weights entry")
    W[a, b] <- W[b, a] <- val
  }
  W
}

#' osr_params_bounds block
#' @noRd
dyn_osr_bounds <- function(statements, osr_params, cal_env) {
  lower <- stats::setNames(rep(-Inf, length(osr_params)), osr_params)
  upper <- stats::setNames(rep(Inf, length(osr_params)), osr_params)
  for (st in statements) {
    f <- dyn_split_commas(st)
    if (length(f) != 3L || !f[1] %in% osr_params) {
      stop("Cannot parse osr_params_bounds entry: ", st, call. = FALSE)
    }
    lower[f[1]] <- dyn_eval(f[2], cal_env, "OSR bound")
    upper[f[1]] <- dyn_eval(f[3], cal_env, "OSR bound")
  }
  list(lower = lower, upper = upper)
}

# ---------------------------------------------------------------------------
# Symbolic helpers on Dynare-timed expressions
# ---------------------------------------------------------------------------

#' Replace timed references v(k) by symbols v__L<k>/v__F<k> (and v by v)
#' @noRd
dyn_to_symbols <- function(expr, vars) {
  dyn_rewrite_ids(expr, function(name, lag, has_index) {
    if (!name %in% vars) return(NULL)
    dyn_sym(name, lag)
  })
}

#' @noRd
dyn_sym <- function(name, lag) {
  if (lag == 0L) name
  else if (lag < 0L) paste0(name, "__L", -lag)
  else paste0(name, "__F", lag)
}

#' Convert symbols back to Dynare timing, shifting every lag by `shift`
#' @noRd
dyn_from_symbols <- function(expr, vars, shift = 0L) {
  dyn_rewrite_ids(expr, function(name, lag, has_index) {
    m <- regmatches(name, regexec("^(.*)__(L|F)([0-9]+)$", name))[[1]]
    if (length(m) == 4L && m[2] %in% vars) {
      k <- as.integer(m[4]) * if (m[3] == "L") -1L else 1L
      return(dyn_timed(m[2], k + shift))
    }
    if (name %in% vars) return(dyn_timed(name, shift))
    NULL
  })
}

#' Symbolic derivative of an R expression string with respect to a symbol
#' @noRd
dyn_deriv <- function(expr, sym) {
  e <- gsub("stats::(pnorm|dnorm)", "\\1", expr)
  d <- stats::D(parse(text = e)[[1]], sym)
  out <- paste(deparse(d, width.cutoff = 500L), collapse = " ")
  out <- gsub("(?<![A-Za-z0-9_.:])(pnorm|dnorm)\\(", "stats::\\1(", out,
              perl = TRUE)
  out
}

#' Timed symbols present in an expression
#' @noRd
dyn_timed_symbols <- function(expr, vars) {
  syms <- all.vars(parse(text = expr))
  base <- sub("__(L|F)[0-9]+$", "", syms)
  syms[base %in% vars]
}

#' @noRd
dyn_sym_lag <- function(sym) {
  m <- regmatches(sym, regexec("__(L|F)([0-9]+)$", sym))[[1]]
  if (length(m) != 3L) return(0L)
  as.integer(m[3]) * if (m[2] == "L") -1L else 1L
}

#' @noRd
dyn_sym_base <- function(sym) sub("__(L|F)[0-9]+$", "", sym)

#' Residual string "(lhs) - (rhs)" of an equation
#' @noRd
dyn_residual <- function(eq) {
  parts <- strsplit(eq, "=", fixed = TRUE)[[1]]
  paste0("(", parts[1], ") - (", parts[2], ")")
}

# ---------------------------------------------------------------------------
# Ramsey (commitment)
# ---------------------------------------------------------------------------

#' Augment the model block with the Ramsey planner's first-order conditions
#'
#' @param eqs Dynare-timed equation strings (translated math).
#' @param objective Planner objective (Dynare syntax).
#' @param discount Planner discount factor as an expression string.
#' @return list(equations, multipliers)
#' @noRd
dyn_ramsey_equations <- function(eqs, objective, endo, exo, discount) {
  if (is.null(objective)) {
    stop("ramsey_model/ramsey_policy requires a planner_objective.",
         call. = FALSE)
  }
  n_mult <- length(eqs)
  mults <- paste0("MULT_", seq_len(n_mult))
  clash <- intersect(mults, c(endo, exo))
  if (length(clash) > 0L) {
    stop("Names reserved for Ramsey multipliers are already used: ",
         paste(clash, collapse = ", "), call. = FALSE)
  }
  vars <- c(endo, exo, mults)
  beta <- paste0("(", if (is.null(discount)) "1" else discount, ")")

  terms <- c(list(dyn_to_symbols(dyn_translate_math(objective), vars)),
             lapply(seq_len(n_mult), function(i) {
               paste0(mults[i], " * (", dyn_to_symbols(dyn_residual(eqs[i]),
                                                       vars), ")")
             }))

  focs <- vapply(endo, function(v) {
    pieces <- character(0)
    for (tm in terms) {
      syms <- dyn_timed_symbols(tm, endo)
      syms <- syms[dyn_sym_base(syms) == v]
      for (sy in syms) {
        tau <- dyn_sym_lag(sy)
        d <- dyn_deriv(tm, sy)
        if (identical(d, "0")) next
        d <- dyn_from_symbols(d, vars, shift = -tau)
        w <- if (tau == 0L) "" else paste0(beta, "^(", -tau, ") * ")
        pieces <- c(pieces, paste0(w, "(", d, ")"))
      }
    }
    if (length(pieces) == 0L) {
      stop("Variable '", v, "' does not appear in the Ramsey problem.",
           call. = FALSE)
    }
    paste0(paste(pieces, collapse = " + "), " = 0")
  }, character(1), USE.NAMES = FALSE)

  list(equations = c(eqs, focs), multipliers = mults)
}

# ---------------------------------------------------------------------------
# Discretion (linear-quadratic)
# ---------------------------------------------------------------------------

#' Time-consistent optimal rule for a linear model with quadratic loss
#'
#' Model: A_lag y_{t-1} + A_0 y_t + A_lead E_t y_{t+1} + B e_t = 0 (m
#' equations in n variables, n - m instruments). Loss: y_t' W y_t
#' discounted at beta. Returns the instrument rule
#' u_t = F1 y_{t-1} + F2 e_t (Dennis 2007).
#' @noRd
dyn_discretion_rule <- function(eqs, objective, endo, exo, instruments,
                                discount, params, cal_env,
                                tol = 1e-12, max_iter = 10000L) {
  if (is.null(objective)) {
    stop("discretionary_policy requires a planner_objective.", call. = FALSE)
  }
  if (is.null(instruments)) {
    stop("discretionary_policy requires instruments = (...).", call. = FALSE)
  }
  unknown <- setdiff(instruments, endo)
  if (length(unknown) > 0L) {
    stop("Unknown instrument(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  n <- length(endo)
  m <- length(eqs)
  if (n - m != length(instruments)) {
    stop("discretionary_policy: the model has ", m, " equations for ", n,
         " variables but ", length(instruments), " instrument(s).",
         call. = FALSE)
  }
  env <- new.env(parent = baseenv())
  for (nm in ls(cal_env)) assign(nm, get(nm, envir = cal_env), envir = env)
  for (nm in names(params)) assign(nm, params[[nm]], envir = env)
  beta <- if (is.null(discount)) 1 else dyn_eval(discount, env, "discount")

  vars <- c(endo, exo)
  timed <- c(paste0(endo, "__L1"), endo, paste0(endo, "__F1"), exo)
  res <- vapply(eqs, function(eq) dyn_to_symbols(dyn_residual(eq), vars),
                character(1))
  bad <- setdiff(unlist(lapply(res, dyn_timed_symbols, vars)), timed)
  if (length(bad) > 0L) {
    stop("discretionary_policy supports leads and lags of one period and ",
         "current shocks only; found: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  exprs <- lapply(res, function(r) parse(text = r)[[1]])
  f <- function(z) {
    e2 <- new.env(parent = env)
    for (k in seq_along(timed)) assign(timed[k], z[k], envir = e2)
    vapply(exprs, eval, numeric(1), envir = e2)
  }
  z0 <- numeric(length(timed))
  f0 <- f(z0)
  J <- vapply(seq_along(timed), function(k) {
    z <- z0
    z[k] <- 1
    f(z) - f0
  }, numeric(m))
  J <- matrix(J, m, length(timed))
  if (max(abs(f0)) > 1e-10) {
    stop("discretionary_policy requires a linear model in deviations ",
         "(model(linear)).", call. = FALSE)
  }
  A_lag <- J[, seq_len(n), drop = FALSE]
  A_0 <- J[, n + seq_len(n), drop = FALSE]
  A_lead <- J[, 2L * n + seq_len(n), drop = FALSE]
  B <- J[, 3L * n + seq_along(exo), drop = FALSE]

  # Quadratic loss matrix W from the objective (current variables only)
  obj <- dyn_to_symbols(dyn_translate_math(objective), vars)
  obj_syms <- dyn_timed_symbols(obj, vars)
  if (any(!obj_syms %in% endo)) {
    stop("discretionary_policy: the planner objective may only contain ",
         "current endogenous variables.", call. = FALSE)
  }
  W <- matrix(0, n, n, dimnames = list(endo, endo))
  for (a in intersect(endo, obj_syms)) {
    da <- dyn_deriv(obj, a)
    for (b in intersect(endo, obj_syms)) {
      dab <- dyn_deriv(da, b)
      W[a, b] <- eval(parse(text = dab), envir = env) / 2
    }
  }

  u_idx <- match(instruments, endo)
  z_idx <- setdiff(seq_len(n), u_idx)
  H1 <- matrix(0, n, n)
  P <- matrix(0, n, n)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    D <- A_0 + A_lead %*% H1
    Dz <- D[, z_idx, drop = FALSE]
    Du <- D[, u_idx, drop = FALSE]
    Dz_inv <- solve(Dz)
    Gam <- matrix(0, n, length(u_idx))
    Gam[z_idx, ] <- -Dz_inv %*% Du
    Gam[u_idx, ] <- diag(length(u_idx))
    Lam <- matrix(0, n, m)
    Lam[z_idx, ] <- -Dz_inv
    S <- W + beta * P
    K <- (diag(n) - Gam %*% solve(t(Gam) %*% S %*% Gam, t(Gam) %*% S)) %*% Lam
    H1_new <- K %*% A_lag
    P_new <- t(H1_new) %*% S %*% H1_new
    delta <- max(abs(H1_new - H1)) + max(abs(P_new - P))
    H1 <- H1_new
    P <- P_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    stop("discretionary_policy: the Dennis (2007) iteration did not ",
         "converge.", call. = FALSE)
  }
  H2 <- K %*% B
  dimnames(H1) <- list(endo, endo)
  dimnames(H2) <- list(endo, exo)
  list(F1 = H1[instruments, , drop = FALSE], F2 = H2[instruments, , drop = FALSE],
       H1 = H1, H2 = H2, W = W, P = P, beta = beta, iterations = it)
}

#' Equation strings for the discretionary rule u = F1 y(-1) + F2 e
#' @noRd
dyn_rule_equations <- function(rule, endo, exo) {
  vapply(rownames(rule$F1), function(u) {
    terms <- character(0)
    for (v in endo) {
      cf <- rule$F1[u, v]
      if (abs(cf) > 1e-14) terms <- c(terms, sprintf("(%.17g) * %s(-1)", cf, v))
    }
    for (e in exo) {
      cf <- rule$F2[u, e]
      if (abs(cf) > 1e-14) terms <- c(terms, sprintf("(%.17g) * %s", cf, e))
    }
    if (length(terms) == 0L) terms <- "0"
    paste0(u, " = ", paste(terms, collapse = " + "))
  }, character(1), USE.NAMES = FALSE)
}


#' Steady state of a Ramsey-augmented model with concentrated multipliers
#'
#' At the steady state the planner's first-order conditions are linear in
#' the Lagrange multipliers, so for given values of the original variables
#' the multipliers solve a least-squares problem. The remaining system in
#' the original variables is solved by Gauss-Newton (as in Dynare's
#' Ramsey steady-state computation).
#' @noRd
dyn_ramsey_ss <- function(model, params, guess, endo, mults, fill,
                          tol = 1e-10, max_iter = 200L) {
  resid_at <- function(v, lam) {
    vals <- c(v, lam)
    names(vals) <- c(endo, mults)
    eval_ss_residual(model, fill(vals)[model$all_variables], params)
  }
  concentrate <- function(v) {
    a <- resid_at(v, numeric(length(mults)))
    B <- vapply(seq_along(mults), function(i) {
      e <- numeric(length(mults))
      e[i] <- 1
      resid_at(v, e) - a
    }, numeric(length(a)))
    B <- matrix(B, length(a), length(mults))
    lam <- -qr.solve(B, a)
    list(lam = lam, r = a + B %*% lam)
  }
  v <- guess[endo]
  for (it in seq_len(max_iter)) {
    cr <- concentrate(v)
    r <- as.numeric(cr$r)
    if (max(abs(r)) < tol) {
      vals <- c(v, cr$lam)
      names(vals) <- c(endo, mults)
      return(fill(vals)[model$all_variables])
    }
    J <- numDeriv::jacobian(function(z) as.numeric(concentrate(z)$r), v)
    step <- qr.solve(J, -r)
    lambda <- 1
    for (ls in seq_len(20L)) {
      v_new <- v + lambda * step
      r_new <- tryCatch(as.numeric(concentrate(v_new)$r),
                        error = function(e) rep(Inf, length(r)))
      if (all(is.finite(r_new)) && sum(r_new^2) < sum(r^2)) break
      lambda <- lambda / 2
    }
    v <- v + lambda * step
  }
  stop("Ramsey steady state did not converge; supply better initval ",
       "values.", call. = FALSE)
}
