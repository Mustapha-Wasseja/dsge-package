# ==========================================================================
# LaTeX export of model equations (Dynare: `write_latex_dynamic_model`)
# ==========================================================================
#
# Renders the equations of a dsge_model or dsgenl_model as LaTeX, with
# Greek-letter substitution for parameter names, time subscripts on
# variables, and expectation operators on leads.
#
# Linear models are rendered from the parsed term structure produced by
# parse_equation(); nonlinear models are rendered by re-parsing the
# stored equation strings into R expressions and walking the AST.
# ==========================================================================


#' Export DSGE Model Equations to LaTeX
#'
#' Renders the equations of a \code{dsge_model} or \code{dsgenl_model} as
#' LaTeX source, suitable for pasting into a paper or writing straight to
#' a \code{.tex} file.  Parameter names that match Greek letters are
#' converted automatically (\code{beta} becomes \code{\\beta}), variables
#' receive time subscripts, and leads are wrapped in a conditional
#' expectation operator.
#'
#' @param model A \code{dsge_model} or \code{dsgenl_model} object.
#' @param file Optional path to write the LaTeX source to.  If
#'   \code{NULL} (default) the source is returned invisibly and printed
#'   to the console.
#' @param env Character.  LaTeX environment to wrap the system in.  One
#'   of \code{"align"} (default), \code{"gather"}, or \code{"none"} (bare
#'   equation lines with no wrapper).
#' @param numbered Logical.  Use the numbered environment
#'   (\code{align}) rather than the starred one (\code{align*}).
#'   Default \code{TRUE}.
#' @param greek Logical.  Substitute Greek-letter parameter names with
#'   their LaTeX commands.  Default \code{TRUE}.
#' @param standalone Logical.  Wrap the output in a minimal complete
#'   LaTeX document (with \code{\\documentclass} and \code{amsmath}) so
#'   the file compiles on its own.  Default \code{FALSE}.
#' @param labels Logical.  Emit a \code{\\label{eq:<lhs>}} on each
#'   equation so they can be cross-referenced.  Only used when
#'   \code{numbered = TRUE}.  Default \code{FALSE}.
#'
#' @return A character vector of LaTeX source lines, returned invisibly.
#'   When \code{file} is supplied the same lines are written there.
#'
#' @details
#' \subsection{Rendering conventions}{
#' \itemize{
#'   \item Variables carry a time subscript: \code{y} renders as
#'     \eqn{y_{t}}.
#'   \item Leads render inside a conditional expectation:
#'     \code{lead(y)} becomes \eqn{\mathbb{E}_{t} y_{t+1}}.
#'   \item Parameter names matching a Greek letter are converted
#'     (\code{sigma} to \code{\\sigma}); a trailing \code{_bar} becomes
#'     an overbar and a trailing \code{_ss} a steady-state superscript.
#'   \item Divisions render as \code{\\frac}, powers as superscripts,
#'     and multiplication as juxtaposition.
#' }
#' }
#'
#' @examples
#' nk <- dsge_model(
#'   obs(pi ~ beta * lead(pi) + kappa * x),
#'   unobs(x ~ lead(x) - (r - lead(pi) - g)),
#'   obs(r ~ psi * pi + u),
#'   state(u ~ rhou * u),
#'   state(g ~ rhog * g),
#'   fixed = list(beta = 0.99),
#'   start = list(kappa = 0.1, psi = 1.5, rhou = 0.7, rhog = 0.9))
#' tex <- model_latex(nk)
#' cat(head(tex, 4), sep = "\n")
#'
#' @export
model_latex <- function(model, file = NULL,
                        env = c("align", "gather", "none"),
                        numbered = TRUE,
                        greek = TRUE,
                        standalone = FALSE,
                        labels = FALSE) {
  env <- match.arg(env)
  is_nl <- inherits(model, "dsgenl_model")
  if (!is_nl && !inherits(model, "dsge_model"))
    stop("`model` must be a dsge_model or dsgenl_model object.",
         call. = FALSE)

  bodies <- if (is_nl) .latex_bodies_nl(model, greek)
            else       .latex_bodies_linear(model, greek)

  # Join the equation bodies into the requested environment
  sep <- if (env == "align") " \\\\" else if (env == "gather") " \\\\" else ""
  lines <- character(0)

  if (standalone) {
    lines <- c(lines,
               "\\documentclass{article}",
               "\\usepackage{amsmath}",
               "\\usepackage{amssymb}",
               "\\begin{document}")
  }

  env_name <- if (env == "none") NA_character_
              else if (numbered) env else paste0(env, "*")

  if (!is.na(env_name)) lines <- c(lines, paste0("\\begin{", env_name, "}"))

  n <- length(bodies)
  for (i in seq_len(n)) {
    body <- bodies[[i]]$tex
    if (labels && numbered && env != "none")
      body <- paste0(body, " \\label{eq:", bodies[[i]]$lhs, "}")
    # No line-break marker on the final equation
    suffix <- if (i < n) sep else ""
    lines <- c(lines, paste0("  ", body, suffix))
  }

  if (!is.na(env_name)) lines <- c(lines, paste0("\\end{", env_name, "}"))
  if (standalone) lines <- c(lines, "\\end{document}")

  if (!is.null(file)) {
    writeLines(lines, con = file)
    message("LaTeX written to: ", file)
  } else {
    cat(lines, sep = "\n")
    cat("\n")
  }
  invisible(lines)
}


# --------------------------------------------------------------------------
# Body builders
# --------------------------------------------------------------------------

#' Build LaTeX bodies for a linear dsge_model from its parsed terms.
#' @noRd
.latex_bodies_linear <- function(model, greek) {
  lapply(model$equations, function(p) {
    # Left-hand side
    lhs_var <- .latex_var(p$lhs_var, is_lead = FALSE, greek = greek)
    lhs <- if (identical(p$lhs_coef_expr, 1)) {
      lhs_var
    } else {
      paste0(.latex_coef(p$lhs_coef_expr, greek), " ", lhs_var)
    }
    # State equations describe x_{t+1}; show that explicitly
    if (identical(p$type, "state")) {
      lhs <- .latex_var(p$lhs_var, is_lead = FALSE, greek = greek,
                        subscript = "t+1")
    }

    # Right-hand side: accumulate signed terms
    rhs <- ""
    for (j in seq_along(p$rhs_terms)) {
      tm <- p$rhs_terms[[j]]
      piece <- .latex_term(tm, greek)
      if (j == 1L) {
        rhs <- piece$tex
        if (piece$negative) rhs <- paste0("-", rhs)
      } else {
        rhs <- paste0(rhs, if (piece$negative) " - " else " + ", piece$tex)
      }
    }

    # State equations carry an additive innovation
    if (identical(p$type, "state") && isTRUE(p$shock)) {
      eps <- if (greek) "\\varepsilon" else "varepsilon"
      rhs <- paste0(rhs, " + ", eps, "_{",
                    .latex_escape(p$lhs_var), ",t+1}")
    }

    list(lhs = p$lhs_var,
         tex = paste0(lhs, " &= ", rhs))
  })
}


#' Build LaTeX bodies for a nonlinear dsgenl_model from its equation
#' strings, by translating `X(+1)` to `lead(X)` and parsing.
#' @noRd
.latex_bodies_nl <- function(model, greek) {
  vars <- unique(unlist(model$variables, use.names = FALSE))
  lapply(model$eq_strings, function(s) {
    parts <- strsplit(s, "=", fixed = TRUE)[[1]]
    if (length(parts) < 2L) {
      return(list(lhs = "", tex = paste0("\\text{", s, "}")))
    }
    lhs_s <- trimws(parts[1])
    rhs_s <- trimws(paste(parts[-1], collapse = "="))
    lhs_tex <- .latex_nl_side(lhs_s, greek, vars)
    rhs_tex <- .latex_nl_side(rhs_s, greek, vars)
    list(lhs = gsub("[^A-Za-z0-9]", "", lhs_s),
         tex = paste0(lhs_tex, " &= ", rhs_tex))
  })
}


#' Convert one side of a nonlinear equation string to LaTeX.
#' @noRd
.latex_nl_side <- function(s, greek, vars = character(0)) {
  # Translate Dynare-ish lead syntax `X(+1)` into `lead(X)` before parsing
  s2 <- gsub("([A-Za-z][A-Za-z0-9_]*)\\s*\\(\\s*\\+\\s*1\\s*\\)",
             "lead(\\1)", s)
  e <- tryCatch(str2lang(s2), error = function(err) NULL)
  if (is.null(e)) return(paste0("\\text{", s, "}"))
  .expr_to_latex(e, greek, vars)
}


#' Render one parsed RHS term (coefficient x variable).
#' @noRd
.latex_term <- function(tm, greek) {
  coef_expr <- tm$coef_expr
  negative <- FALSE

  # Detect a leading unary minus so we can render "- x" rather than "+ -x"
  if (is.numeric(coef_expr) && length(coef_expr) == 1L && coef_expr < 0) {
    negative <- TRUE
    coef_expr <- if (coef_expr == -1) 1 else abs(coef_expr)
  } else if (is.call(coef_expr) &&
             identical(as.character(coef_expr[[1L]]), "-") &&
             length(coef_expr) == 2L) {
    negative <- TRUE
    coef_expr <- coef_expr[[2L]]
  }

  var_tex <- .latex_var(tm$variable, is_lead = isTRUE(tm$is_lead),
                        greek = greek, lead_k = tm$lead_k)

  tex <- if (identical(coef_expr, 1)) var_tex
         else paste0(.latex_coef(coef_expr, greek), " ", var_tex)

  list(tex = tex, negative = negative)
}


#' Render a coefficient expression, adding parentheses when the top-level
#' operator is additive (so `(1 - alpha) x` reads correctly).
#' @noRd
.latex_coef <- function(e, greek) {
  tex <- .expr_to_latex(e, greek)
  if (is.call(e)) {
    op <- as.character(e[[1L]])
    if (op %in% c("+", "-") && length(e) == 3L)
      return(paste0("\\left(", tex, "\\right)"))
  }
  tex
}


#' Render a model variable with a time subscript, wrapping leads in a
#' conditional expectation operator.
#' @noRd
.latex_var <- function(name, is_lead = FALSE, greek = TRUE,
                       lead_k = 1L, subscript = NULL) {
  base <- .latex_sym(name, greek)
  if (!is.null(subscript))
    return(paste0(base, "_{", subscript, "}"))
  if (isTRUE(is_lead)) {
    k <- if (is.null(lead_k) || is.na(lead_k)) 1L else as.integer(lead_k)
    sub <- if (k == 1L) "t+1" else paste0("t+", k)
    return(paste0("\\mathbb{E}_{t} ", base, "_{", sub, "}"))
  }
  paste0(base, "_{t}")
}


# --------------------------------------------------------------------------
# Expression -> LaTeX
# --------------------------------------------------------------------------

#' Recursively convert an R expression to LaTeX source.
#' @noRd
.expr_to_latex <- function(e, greek = TRUE,
                           vars = character(0)) {
  if (is.numeric(e) || is.logical(e)) return(.latex_num(e))
  if (is.character(e)) return(.latex_sym(e, greek))
  if (is.name(e)) {
    nm <- as.character(e)
    if (nm %in% vars)
      return(.latex_var(nm, is_lead = FALSE, greek = greek))
    return(.latex_sym(nm, greek))
  }

  if (is.call(e)) {
    op <- as.character(e[[1L]])

    if (op == "(" && length(e) == 2L)
      return(paste0("\\left(", .expr_to_latex(e[[2L]], greek, vars), "\\right)"))

    if (op %in% c("+", "-") && length(e) == 3L)
      return(paste(.expr_to_latex(e[[2L]], greek, vars), op,
                   .expr_to_latex(e[[3L]], greek, vars)))

    if (op == "-" && length(e) == 2L)
      return(paste0("-", .expr_to_latex(e[[2L]], greek, vars)))

    if (op == "+" && length(e) == 2L)
      return(.expr_to_latex(e[[2L]], greek, vars))

    if (op == "/" && length(e) == 3L)
      return(paste0("\\frac{", .expr_to_latex(e[[2L]], greek, vars), "}{",
                    .expr_to_latex(e[[3L]], greek, vars), "}"))

    if (op == "*" && length(e) == 3L) {
      l <- .latex_wrap_if_additive(e[[2L]], greek, vars)
      r <- .latex_wrap_if_additive(e[[3L]], greek, vars)
      return(paste0(l, " ", r))
    }

    if (op == "^" && length(e) == 3L) {
      base <- .latex_wrap_if_additive(e[[2L]], greek, vars)
      return(paste0(base, "^{", .expr_to_latex(e[[3L]], greek, vars), "}"))
    }

    if (op %in% c("lead", "E") && length(e) >= 2L) {
      k <- if (op == "lead" && length(e) >= 3L)
             as.integer(eval(e[[3L]])) else 1L
      return(.latex_var(as.character(e[[2L]]), is_lead = TRUE,
                        greek = greek, lead_k = k))
    }

    if (op == "exp")  return(paste0("\\exp\\left(",
                                    .expr_to_latex(e[[2L]], greek, vars),
                                    "\\right)"))
    if (op == "log")  return(paste0("\\log\\left(",
                                    .expr_to_latex(e[[2L]], greek, vars),
                                    "\\right)"))
    if (op == "sqrt") return(paste0("\\sqrt{",
                                    .expr_to_latex(e[[2L]], greek, vars), "}"))

    # Generic function call: f(a, b)
    args <- vapply(as.list(e)[-1L],
                   function(a) .expr_to_latex(a, greek, vars), character(1))
    return(paste0("\\mathrm{", .latex_escape(op), "}\\left(",
                  paste(args, collapse = ", "), "\\right)"))
  }

  .latex_escape(paste(deparse(e), collapse = ""))
}


#' Wrap a sub-expression in parentheses when it is an additive term.
#' @noRd
.latex_wrap_if_additive <- function(e, greek, vars = character(0)) {
  tex <- .expr_to_latex(e, greek, vars)
  if (is.call(e)) {
    op <- as.character(e[[1L]])
    if (op %in% c("+", "-") && length(e) == 3L)
      return(paste0("\\left(", tex, "\\right)"))
  }
  tex
}


#' Format a numeric literal without scientific notation clutter.
#' @noRd
.latex_num <- function(x) {
  if (is.logical(x)) return(as.character(as.integer(x)))
  out <- format(x, trim = TRUE, digits = 6, scientific = FALSE)
  sub("\\.0+$", "", out)
}


#' Convert a symbol name to LaTeX, applying Greek substitution and
#' subscript / decoration conventions.
#' @noRd
.latex_sym <- function(name, greek = TRUE) {
  name <- as.character(name)

  # Split trailing decorations first: foo_bar -> \bar{foo}, foo_ss -> foo^{ss}
  deco <- NULL
  if (grepl("_bar$", name)) {
    name <- sub("_bar$", "", name); deco <- "bar"
  } else if (grepl("_ss$", name)) {
    name <- sub("_ss$", "", name);  deco <- "ss"
  } else if (grepl("_star$", name)) {
    name <- sub("_star$", "", name); deco <- "star"
  }

  # Split remaining name into stem + subscript on the first separator
  parts <- strsplit(name, "[_.]", perl = TRUE)[[1]]
  stem  <- parts[1]
  subs  <- if (length(parts) > 1L) paste(parts[-1], collapse = ",") else NULL

  stem_tex <- if (greek) .latex_greek(stem) else .latex_escape(stem)
  out <- stem_tex
  if (!is.null(subs)) out <- paste0(out, "_{", .latex_escape(subs), "}")

  if (identical(deco, "bar"))  out <- paste0("\\bar{", stem_tex, "}",
                                             if (!is.null(subs))
                                               paste0("_{",
                                                      .latex_escape(subs),
                                                      "}") else "")
  if (identical(deco, "ss"))   out <- paste0(out, "^{ss}")
  if (identical(deco, "star")) out <- paste0(out, "^{*}")
  out
}


#' Map a bare name to its Greek LaTeX command when it matches one.
#' @noRd
.latex_greek <- function(x) {
  greek_lower <- c("alpha","beta","gamma","delta","epsilon","varepsilon",
                   "zeta","eta","theta","vartheta","iota","kappa","lambda",
                   "mu","nu","xi","rho","varrho","sigma","varsigma","tau",
                   "upsilon","phi","varphi","chi","psi","omega","pi","varpi")
  greek_upper <- c("Gamma","Delta","Theta","Lambda","Xi","Pi","Sigma",
                   "Upsilon","Phi","Psi","Omega")
  if (x %in% greek_lower || x %in% greek_upper) return(paste0("\\", x))
  # Capitalised spelling of a lower-case Greek name (e.g. "Beta") is not a
  # LaTeX command; fall through to plain text in that case.
  .latex_escape(x)
}


#' Escape the LaTeX special characters that can appear in variable names.
#' @noRd
.latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}
