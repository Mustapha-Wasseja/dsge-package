# Compare second- and third-order solutions of imported models with Dynare.
#
# For each .mod file, Dynare solves the model at order 2 and 3 (after the
# file is cut at its first computing command, as in batch_dsge_mod.R) and
# exports its decision rules (ghx, ghu, ghxx, ghxu, ghuu, ghs2 and, at
# order 3, ghxxx, ghxxu, ghxuu, ghuuu, ghxss, ghuss). dsge solves the
# imported model at the same order. Both decision rules are evaluated at
# the same random points (states drawn from the first-order ergodic
# distribution, shocks from their distribution), and the largest absolute
# difference over all variables and points is reported.
#
# Usage: Rscript dev/dynare-validation/extra/validate_higher_order.R \
#          [ORDER] file1.mod file2.mod ...

suppressMessages(devtools::load_all(quiet = TRUE))
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)))
source(file.path(script_dir, "mod_tools.R"))
args <- commandArgs(trailingOnly = TRUE)
orders <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2:3
if (length(args) && grepl("^[0-9]+$", args[1])) args <- args[-1]
dynare_path <- Sys.getenv("DYNARE_MATLAB", "/usr/lib/dynare/matlab")
n_points <- 20L

export_m <- c(
  "dr = oo_.dr;",
  "wr = @(nm, x) dlmwrite([nm '.csv'], full(x), 'precision', '%.17g');",
  "wr('ys', dr.ys); wr('order_var', dr.order_var); wr('state_var', dr.order_var(M_.nstatic+(1:M_.nspred)));",
  "wr('ghx', dr.ghx); wr('ghu', dr.ghu); wr('Sigma_e', M_.Sigma_e);",
  "if isfield(dr, 'ghxx'), wr('ghxx', dr.ghxx); wr('ghxu', dr.ghxu); wr('ghuu', dr.ghuu); wr('ghs2', dr.ghs2); end",
  "if options_.order >= 3",
  "  wr('ghxxx', dr.ghxxx); wr('ghxxu', dr.ghxxu); wr('ghxuu', dr.ghxuu);",
  "  wr('ghuuu', dr.ghuuu); wr('ghxss', dr.ghxss); wr('ghuss', dr.ghuss);",
  "end",
  "fid = fopen('names.csv', 'w');",
  "for i = 1:M_.endo_nbr, fprintf(fid, '%s\\n', M_.endo_names{i}); end; fclose(fid);",
  "fid = fopen('exo.csv', 'w');",
  "for i = 1:M_.exo_nbr, fprintf(fid, '%s\\n', M_.exo_names{i}); end; fclose(fid);",
  "fid = fopen('aux.csv', 'w');",
  "for i = 1:numel(M_.aux_vars)",
  "  a = M_.aux_vars(i); ol = 0; oi = 0;",
  "  if isfield(a, 'orig_lead_lag') && ~isempty(a.orig_lead_lag), ol = a.orig_lead_lag; end",
  "  if isfield(a, 'orig_index') && ~isempty(a.orig_index), oi = a.orig_index; end",
  "  fprintf(fid, '%d,%d,%d,%d\\n', a.endo_index, a.type, oi, ol);",
  "end; fclose(fid);")

rd <- function(dir, nm) {
  f <- file.path(dir, paste0(nm, ".csv"))
  if (!file.exists(f) || file.size(f) == 0) return(matrix(numeric(0), 0, 0))
  as.matrix(utils::read.csv(f, header = FALSE))
}

kron_pow <- function(a, b = a, c = NULL) {
  out <- kronecker(a, b)
  if (!is.null(c)) out <- kronecker(out, c)
  out
}

compare_one <- function(mod, order) {
  prep <- prepare_mod(mod, paste0("ho", order), order = order, irf = 0L)
  writeLines(c(sprintf("addpath('%s');", dynare_path),
               sprintf("dynare %s noclearall nolog", prep$name), export_m),
             file.path(prep$dir, "run.m"))
  old <- setwd(prep$dir)
  on.exit(setwd(old))
  system2("timeout", c("900", "octave", "--no-gui", "--quiet", "run.m"),
          stdout = "octave.log", stderr = "octave.log")
  setwd(old)
  if (file.exists(file.path(prep$dir, "ghx.csv")) &&
      !file.exists(file.path(prep$dir, "ghxx.csv"))) {
    return(list(status = "skipped", msg = "linear model (first order only)"))
  }
  if (!file.exists(file.path(prep$dir, "ghx.csv"))) {
    err <- grep("error", readLines(file.path(prep$dir, "octave.log"), warn = FALSE),
                value = TRUE, ignore.case = TRUE)
    return(list(status = "dynare_failed", msg = paste(utils::head(err, 1))))
  }
  d <- list()
  for (nm in c("ys", "order_var", "state_var", "ghx", "ghu", "ghxx", "ghxu",
               "ghuu", "ghs2", "Sigma_e", "ghxxx", "ghxxu", "ghxuu", "ghuuu",
               "ghxss", "ghuss")) d[[nm]] <- rd(prep$dir, nm)
  if (!file.exists(file.path(prep$dir, "names.csv"))) {
    log <- readLines(file.path(prep$dir, "octave.log"), warn = FALSE)
    return(list(status = "dynare_failed",
                msg = paste(utils::tail(grep("error", log, value = TRUE), 2),
                            collapse = " | ")))
  }
  endo_names <- readLines(file.path(prep$dir, "names.csv"))
  exo_names <- readLines(file.path(prep$dir, "exo.csv"))
  aux <- if (file.size(file.path(prep$dir, "aux.csv")) > 0) {
    utils::read.csv(file.path(prep$dir, "aux.csv"), header = FALSE)
  } else NULL
  # correlated shocks: dsge's shock states are the orthogonal components,
  # u = T e with T = L diag(1/diag(L)) (unit lower-triangular Cholesky)
  Lc <- t(chol(d$Sigma_e + diag(1e-300, nrow(d$Sigma_e))))
  Tc <- Lc %*% diag(1 / pmax(diag(Lc), 1e-300), nrow(Lc))

  imp <- tryCatch(read_dynare(prep$mod), error = function(e) e)
  if (inherits(imp, "error")) return(list(status = "import_failed",
                                          msg = conditionMessage(imp)))
  sol <- tryCatch({
    setTimeLimit(elapsed = 1800)
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
    solve_dsge(imp, order = order)
  }, error = function(e) e)
  if (inherits(sol, "error")) return(list(status = "solve_failed",
                                          msg = conditionMessage(sol)))

  # Dynare state variable -> dsge state name
  sv <- as.integer(d$state_var)
  dyn_state_name <- vapply(sv, function(k) {
    if (!is.null(aux) && k %in% aux$V1) {
      a <- aux[aux$V1 == k, ]
      base <- if (a$V2 %in% c(1, 3) && a$V3 > 0) {
        if (a$V2 == 1) endo_names[a$V3] else exo_names[a$V3]
      } else NA
      if (is.na(base) || a$V4 > 0 || (a$V2 == 1 && a$V4 == 0)) {
        return(paste0("<aux type ", a$V2, " orig ", a$V3, " lag ", a$V4, ">"))
      }
      return(paste0(base, "_lag", -a$V4 + 1L))
    }
    paste0(endo_names[k], "_lag1")
  }, "")
  states <- colnames(sol$H)
  if (!all(dyn_state_name %in% states)) {
    bad <- setdiff(dyn_state_name, states)
    return(list(status = "skipped", msg = paste("state mapping:",
                                                 paste(bad, collapse = " "))))
  }
  shock_states <- intersect(exo_names, states)
  sd_u <- sqrt(diag(d$Sigma_e))
  P <- compute_unconditional_P(sol$H, sol$M %*% t(sol$M))
  dimnames(P) <- dimnames(sol$H)
  lagst <- dyn_state_name
  set.seed(1)
  ov <- as.integer(d$order_var)
  ys <- as.numeric(d$ys)
  vars <- intersect(imp$variables, rownames(sol$G))
  worst <- 0
  scale <- 0
  for (k in seq_len(n_points)) {
    Pl <- P[lagst, lagst, drop = FALSE]
    ev <- eigen((Pl + t(Pl)) / 2, symmetric = TRUE)
    xs <- as.numeric(ev$vectors %*% (sqrt(pmax(ev$values, 0)) * stats::rnorm(length(lagst))))
    u <- stats::rnorm(length(exo_names)) * sd_u
    # Dynare
    ydr <- ys[ov] + 0.5 * d$ghs2 + d$ghx %*% xs + d$ghu %*% u +
      0.5 * d$ghxx %*% kron_pow(xs) + d$ghxu %*% kronecker(xs, u) +
      0.5 * d$ghuu %*% kron_pow(u)
    if (order >= 3) {
      ydr <- ydr + (d$ghxxx %*% kron_pow(xs, xs, xs) + d$ghuuu %*% kron_pow(u, u, u) +
                    3 * d$ghxxu %*% kron_pow(xs, xs, u) +
                    3 * d$ghxuu %*% kron_pow(xs, u, u) +
                    3 * d$ghxss %*% xs + 3 * d$ghuss %*% u) / 6
    }
    y_dyn <- stats::setNames(numeric(length(ys)), endo_names)
    y_dyn[ov] <- ydr
    # dsge
    x <- stats::setNames(numeric(length(states)), states)
    x[lagst] <- xs
    x[shock_states] <- solve(Tc, u)[match(shock_states, exo_names)]
    yd <- sol$G %*% x + 0.5 * sol$g_ss
    for (i in seq_len(nrow(sol$G))) {
      yd[i] <- yd[i] + 0.5 * sum(sol$g_xx[i, , ] * (x %o% x))
      if (order >= 3) {
        yd[i] <- yd[i] + sum(sol$g_xxx[i, , , ] * (x %o% x %o% x)) / 6 +
          0.5 * sum(sol$g_xss[i, ] * x) + sol$g_sss[i] / 6
      }
    }
    y_dsge <- stats::setNames(as.numeric(yd) + sol$steady_state[rownames(sol$G)],
                              rownames(sol$G))
    if (nzchar(Sys.getenv("DEBUG_HO")) && k == 1L) {
      decl <- function(v) {
        out <- stats::setNames(numeric(length(ys)), endo_names)
        out[ov] <- v
        out[vars]
      }
      q_dyn <- decl(0.5 * d$ghxx %*% kron_pow(xs) + d$ghxu %*% kronecker(xs, u) +
                      0.5 * d$ghuu %*% kron_pow(u))
      iv <- match(vars, rownames(sol$G))
      q_dsge <- vapply(iv, function(v) 0.5 * sum(sol$g_xx[v, , ] * (x %o% x)), 0)
      print(rbind(const_dyn = decl(0.5 * d$ghs2), const_dsge = 0.5 * sol$g_ss[iv],
                  lin_dyn = decl(d$ghx %*% xs + d$ghu %*% u),
                  lin_dsge = as.numeric(sol$G[vars, ] %*% x),
                  quad_dyn = q_dyn, quad_dsge = q_dsge,
                  ss_dyn = ys[match(vars, endo_names)],
                  ss_dsge = sol$steady_state[vars]))
      if (order >= 3) {
        print(rbind(
          xss_dyn = decl(0.5 * (d$ghxss %*% xs + d$ghuss %*% u)),
          xss_dsge = vapply(iv, function(v) 0.5 * sum(sol$g_xss[v, ] * x), 0),
          cub_dyn = decl((d$ghxxx %*% kron_pow(xs, xs, xs) + d$ghuuu %*% kron_pow(u, u, u) +
                          3 * d$ghxxu %*% kron_pow(xs, xs, u) +
                          3 * d$ghxuu %*% kron_pow(xs, u, u)) / 6),
          cub_dsge = vapply(iv, function(v) sum(sol$g_xxx[v, , , ] * (x %o% x %o% x)) / 6, 0),
          sss_dsge = sol$g_sss[iv] / 6))
      }
    }
    diff <- abs(y_dsge[vars] - y_dyn[vars])
    worst <- max(worst, diff)
    scale <- max(scale, abs(y_dyn[vars] - ys[match(vars, endo_names)]))
  }
  list(status = "ok", max_abs_diff = worst, scale = scale, n_vars = length(vars))
}

out <- list()
for (mod in args) {
  if (is_perfect_foresight(mod)) next
  for (ord in orders) {
    t0 <- Sys.time()
    r <- tryCatch(compare_one(mod, ord),
                  error = function(e) list(status = "harness_error",
                                           msg = conditionMessage(e)))
    row <- data.frame(model = basename(mod), order = ord, status = r$status,
                      vars = if (is.null(r$n_vars)) NA else r$n_vars,
                      max_abs_diff = if (is.null(r$max_abs_diff)) NA else signif(r$max_abs_diff, 3),
                      max_deviation = if (is.null(r$scale)) NA else signif(r$scale, 3),
                      relative = if (is.null(r$scale)) NA else
                        signif(r$max_abs_diff / max(1, r$scale), 3),
                      seconds = round(as.numeric(Sys.time() - t0, units = "secs")),
                      message = if (is.null(r$msg)) "" else substr(r$msg, 1, 80))
    print(row, row.names = FALSE)
    out[[length(out) + 1L]] <- row
  }
}
res <- do.call(rbind, out)
print(res, row.names = FALSE)
