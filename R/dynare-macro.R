# Dynare macro processor
#
# Expands the directives of Dynare's macro language (@#define, @#if,
# @#ifdef, @#ifndef, @#elseif, @#else, @#endif, @#for ... @#endfor,
# @#include, @#includepath, @#echo, @#error) and @{expr} substitutions,
# so that read_dynare() can import macro-based .mod files directly.
#
# Macro values are numbers, strings, booleans and arrays (R lists with
# class "dyn_array"). Expressions are translated into R calls on small
# helper functions and evaluated in a dedicated environment.

#' Expand Dynare macro directives
#'
#' @param src Source text (a single string, comments already removed).
#' @param dir Directory used to resolve @#include paths.
#' @param defines Named list of macro variables defined up front (like
#'   Dynare's -D command-line option).
#' @return Expanded source text.
#' @noRd
dyn_macro_expand <- function(src, dir = ".", defines = list()) {
  state <- new.env(parent = emptyenv())
  state$vars <- dyn_macro_env()
  state$include_path <- dir
  for (nm in names(defines)) {
    assign(nm, dyn_macro_value(defines[[nm]]), envir = state$vars)
  }
  lines <- dyn_macro_lines(src)
  out <- dyn_macro_run(lines, state, dir)
  paste(out, collapse = "\n")
}

#' Split into lines, joining directive continuation lines ending in '\'
#' @noRd
dyn_macro_lines <- function(src) {
  lines <- strsplit(src, "\n", fixed = TRUE)[[1]]
  out <- character(0)
  buf <- NULL
  for (ln in lines) {
    if (!is.null(buf)) {
      ln <- paste(buf, ln)
      buf <- NULL
    }
    if (grepl("^\\s*@#", ln) && grepl("\\\\\\s*$", ln)) {
      buf <- sub("\\\\\\s*$", "", ln)
      next
    }
    out <- c(out, ln)
  }
  if (!is.null(buf)) out <- c(out, buf)
  out
}

#' Parse a directive line into keyword and argument text
#' @noRd
dyn_macro_directive <- function(line) {
  body <- sub("^\\s*@#\\s*", "", line)
  kw <- sub("^([A-Za-z_]+).*$", "\\1", body)
  arg <- trimws(substr(body, nchar(kw) + 1L, nchar(body)))
  list(kw = kw, arg = arg)
}

#' Find the index of the directive closing a block opened at `start`
#' @noRd
dyn_macro_find_end <- function(lines, start, open_kw, close_kw,
                               mid_kw = character(0)) {
  depth <- 0L
  mids <- integer(0)
  for (k in seq(start, length(lines))) {
    if (!grepl("^\\s*@#", lines[k])) next
    kw <- dyn_macro_directive(lines[k])$kw
    if (kw %in% open_kw) depth <- depth + 1L
    if (kw == close_kw) {
      depth <- depth - 1L
      if (depth == 0L) return(list(end = k, mids = mids))
    }
    if (depth == 1L && kw %in% mid_kw && k != start) mids <- c(mids, k)
  }
  stop("Macro block starting with '", trimws(lines[start]),
       "' is not closed with @#", close_kw, ".", call. = FALSE)
}

#' Process a vector of lines
#' @noRd
dyn_macro_run <- function(lines, state, dir) {
  out <- character(0)
  i <- 1L
  n <- length(lines)
  if_kw <- c("if", "ifdef", "ifndef")
  while (i <= n) {
    ln <- lines[i]
    if (!grepl("^\\s*@#", ln)) {
      out <- c(out, dyn_macro_substitute(ln, state$vars))
      i <- i + 1L
      next
    }
    d <- dyn_macro_directive(ln)
    kw <- d$kw
    if (kw == "define") {
      dyn_macro_define(d$arg, state$vars)
      i <- i + 1L
    } else if (kw %in% if_kw) {
      blk <- dyn_macro_find_end(lines, i, if_kw, "endif", c("elseif", "else"))
      starts <- c(i, blk$mids)
      ends <- c(blk$mids, blk$end)
      for (b in seq_along(starts)) {
        dd <- dyn_macro_directive(lines[starts[b]])
        take <- switch(dd$kw,
          "if" = , "elseif" = dyn_macro_truthy(dyn_macro_eval(dd$arg,
                                                              state$vars)),
          "ifdef" = exists(dd$arg, envir = state$vars, inherits = FALSE),
          "ifndef" = !exists(dd$arg, envir = state$vars, inherits = FALSE),
          "else" = TRUE
        )
        if (take) {
          body <- lines[seq_len(ends[b] - starts[b] - 1L) + starts[b]]
          out <- c(out, dyn_macro_run(body, state, dir))
          break
        }
      }
      i <- blk$end + 1L
    } else if (kw == "for") {
      blk <- dyn_macro_find_end(lines, i, "for", "endfor")
      body <- lines[seq_len(blk$end - i - 1L) + i]
      out <- c(out, dyn_macro_for(d$arg, body, state, dir))
      i <- blk$end + 1L
    } else if (kw == "include") {
      out <- c(out, dyn_macro_include(d$arg, state, dir))
      i <- i + 1L
    } else if (kw == "includepath") {
      p <- dyn_macro_eval(d$arg, state$vars)
      state$include_path <- c(state$include_path, as.character(p))
      i <- i + 1L
    } else if (kw == "echo") {
      message("Dynare macro: ",
              dyn_macro_format(dyn_macro_eval(d$arg, state$vars)))
      i <- i + 1L
    } else if (kw == "error") {
      stop("Dynare macro @#error: ",
           dyn_macro_format(dyn_macro_eval(d$arg, state$vars)),
           call. = FALSE)
    } else if (kw %in% c("echomacrovars", "line")) {
      i <- i + 1L
    } else if (kw %in% c("endif", "endfor", "else", "elseif")) {
      stop("Unexpected @#", kw, " without a matching opening directive.",
           call. = FALSE)
    } else {
      stop("Unsupported macro directive: @#", kw, call. = FALSE)
    }
  }
  out
}

#' @#define NAME = EXPR  or  @#define f(x, y) = EXPR
#' @noRd
dyn_macro_define <- function(arg, env) {
  if (!grepl("=", arg)) {
    stop("Cannot parse @#define ", arg, call. = FALSE)
  }
  lhs <- trimws(sub("=.*$", "", arg))
  rhs <- trimws(sub("^[^=]*=", "", arg))
  if (grepl("^[A-Za-z_][A-Za-z0-9_]*$", lhs)) {
    assign(lhs, dyn_macro_eval(rhs, env), envir = env)
    return(invisible(NULL))
  }
  fm <- regmatches(lhs, regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*\\((.*)\\)$",
                                lhs))[[1]]
  if (length(fm) != 3L) {
    stop("Cannot parse @#define ", arg, call. = FALSE)
  }
  args <- trimws(strsplit(fm[3], ",", fixed = TRUE)[[1]])
  code <- dyn_macro_translate(rhs)
  f <- function(...) {
    vals <- list(...)
    local_env <- new.env(parent = env)
    for (k in seq_along(args)) assign(args[k], vals[[k]], envir = local_env)
    eval(parse(text = code), envir = local_env)
  }
  assign(fm[2], f, envir = env)
  invisible(NULL)
}

#' @#for VAR in EXPR [when COND]  /  @#for (A, B) in EXPR
#' @noRd
dyn_macro_for <- function(arg, body, state, dir) {
  m <- regmatches(arg, regexec(
    "^\\(?\\s*([A-Za-z_][A-Za-z0-9_, ]*?)\\s*\\)?\\s+in\\s+(.+)$", arg))[[1]]
  if (length(m) != 3L) {
    stop("Cannot parse @#for ", arg, call. = FALSE)
  }
  loop_vars <- trimws(strsplit(m[2], ",", fixed = TRUE)[[1]])
  rest <- m[3]
  cond <- NULL
  wm <- regexpr("\\s+when\\s+", rest)
  if (wm > 0) {
    cond <- substr(rest, wm + attr(wm, "match.length"), nchar(rest))
    rest <- substr(rest, 1L, wm - 1L)
  }
  seq_vals <- dyn_macro_as_array(dyn_macro_eval(rest, state$vars))
  out <- character(0)
  for (v in seq_vals) {
    if (length(loop_vars) == 1L) {
      assign(loop_vars, v, envir = state$vars)
    } else {
      tup <- dyn_macro_as_array(v)
      for (k in seq_along(loop_vars)) {
        assign(loop_vars[k], tup[[k]], envir = state$vars)
      }
    }
    if (!is.null(cond) &&
        !dyn_macro_truthy(dyn_macro_eval(cond, state$vars))) next
    out <- c(out, dyn_macro_run(body, state, dir))
  }
  out
}

#' @#include "file"
#' @noRd
dyn_macro_include <- function(arg, state, dir) {
  f <- as.character(dyn_macro_eval(arg, state$vars))
  cands <- unique(c(if (grepl("^(/|[A-Za-z]:)", f)) f,
                    file.path(c(dir, state$include_path), f)))
  path <- cands[file.exists(cands)][1]
  if (is.na(path)) {
    stop("@#include file not found: ", f, call. = FALSE)
  }
  src <- dyn_strip_comments(paste(readLines(path, warn = FALSE),
                                  collapse = "\n"))
  dyn_macro_run(dyn_macro_lines(src), state, dirname(path))
}

#' Replace @{expr} in a text line
#' @noRd
dyn_macro_substitute <- function(line, env) {
  repeat {
    pos <- regexpr("@\\{", line)
    if (pos < 0) return(line)
    close <- regexpr("\\}", substr(line, pos + 2L, nchar(line)))
    if (close < 0) {
      stop("Unterminated @{ in: ", line, call. = FALSE)
    }
    expr <- substr(line, pos + 2L, pos + close)
    val <- dyn_macro_format(dyn_macro_eval(expr, env))
    line <- paste0(substr(line, 1L, pos - 1L), val,
                   substr(line, pos + close + 2L, nchar(line)))
  }
}

# ---------------------------------------------------------------------------
# Values and expression evaluation
# ---------------------------------------------------------------------------

#' @noRd
dyn_array <- function(x) structure(as.list(x), class = "dyn_array")

#' @noRd
dyn_macro_value <- function(x) {
  if (is.list(x) || length(x) > 1L) dyn_array(x) else x
}

#' @noRd
dyn_macro_as_array <- function(x) {
  if (inherits(x, "dyn_array")) unclass(x) else as.list(x)
}

#' @noRd
dyn_macro_truthy <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  if (is.numeric(x)) return(length(x) == 1L && x != 0)
  stop("Macro condition is not a boolean or number.", call. = FALSE)
}

#' Format a macro value for substitution into the text
#' @noRd
dyn_macro_format <- function(x, quote = FALSE) {
  if (inherits(x, "dyn_array")) {
    return(paste0("[", paste(vapply(unclass(x), dyn_macro_format, "",
                                    quote = TRUE), collapse = ", "), "]"))
  }
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.numeric(x)) {
    if (x == round(x) && abs(x) < 1e15) return(format(x, scientific = FALSE))
    return(format(x, digits = 15))
  }
  if (quote) paste0("\"", x, "\"") else as.character(x)
}

#' Environment with the macro-language helpers
#' @noRd
dyn_macro_env <- function() {
  helpers <- new.env(parent = baseenv())
  helpers$true <- TRUE
  helpers$false <- FALSE
  helpers$.arr <- function(...) {
    parts <- list(...)
    if (length(parts) == 1L && inherits(parts[[1]], "dyn_array")) {
      return(parts[[1]])
    }
    dyn_array(parts)
  }
  helpers$.range <- function(a, b, by = 1) {
    if (by > 0 && a > b || by < 0 && a < b) return(dyn_array(list()))
    dyn_array(seq(a, b, by = by))
  }
  helpers$.idx <- function(x, i) {
    xs <- dyn_macro_as_array(x)
    if (is.character(x) && !inherits(x, "dyn_array")) {
      chars <- strsplit(x, "")[[1]]
      return(paste(chars[unlist(dyn_macro_as_array(i))], collapse = ""))
    }
    ii <- unlist(dyn_macro_as_array(i))
    if (length(ii) == 1L && !inherits(i, "dyn_array")) return(xs[[ii]])
    dyn_array(xs[ii])
  }
  helpers$.in <- function(x, arr) {
    any(vapply(dyn_macro_as_array(arr), function(el) identical(el, x),
               logical(1)))
  }
  helpers$.add <- function(e1, e2) {
    if (missing(e2)) return(e1)
    if (inherits(e1, "dyn_array") || inherits(e2, "dyn_array")) {
      return(dyn_array(c(dyn_macro_as_array(e1), dyn_macro_as_array(e2))))
    }
    if (is.character(e1) || is.character(e2)) return(paste0(e1, e2))
    e1 + e2
  }
  helpers$.eq <- function(a, b) identical(a, b) ||
    (is.numeric(a) && is.numeric(b) && isTRUE(all.equal(a, b)))
  helpers$.neq <- function(a, b) !helpers$.eq(a, b)
  helpers$length <- function(x) length(dyn_macro_as_array(x))
  helpers$isempty <- function(x) length(dyn_macro_as_array(x)) == 0L
  helpers$sum <- function(x) sum(unlist(dyn_macro_as_array(x)))
  helpers$min <- function(...) min(unlist(lapply(list(...), dyn_macro_as_array)))
  helpers$max <- function(...) max(unlist(lapply(list(...), dyn_macro_as_array)))
  helpers$floor <- floor
  helpers$ceil <- ceiling
  helpers$round <- round
  helpers$trunc <- trunc
  helpers$mod <- function(a, b) a %% b
  helpers$abs <- abs
  helpers$sqrt <- sqrt
  helpers$exp <- exp
  helpers$log <- log
  helpers$ln <- log
  helpers$str <- function(x) dyn_macro_format(x)
  helpers$defined <- function(name) {
    exists(as.character(substitute(name)), envir = parent.frame())
  }
  new.env(parent = helpers)
}

#' Evaluate a macro expression
#' @noRd
dyn_macro_eval <- function(expr, env) {
  code <- dyn_macro_translate(expr)
  tryCatch(eval(parse(text = code), envir = env), error = function(e) {
    stop("Cannot evaluate macro expression '", expr, "': ",
         conditionMessage(e), call. = FALSE)
  })
}

#' Tokenise a macro expression
#' @noRd
dyn_macro_tokens <- function(expr) {
  pat <- paste0(
    "\"[^\"]*\"|",                         # string
    "[0-9]+\\.?[0-9]*(?:[eE][+-]?[0-9]+)?|\\.[0-9]+|",  # number
    "[A-Za-z_][A-Za-z0-9_]*|",             # identifier
    "&&|\\|\\||==|!=|<=|>=|[-+*/^<>!:,()\\[\\]]"
  )
  m <- gregexpr(pat, expr, perl = TRUE)
  toks <- regmatches(expr, m)[[1]]
  leftover <- gsub("\\s", "", gsub(pat, "", expr, perl = TRUE))
  if (nzchar(leftover)) {
    stop("Unsupported characters in macro expression: ", expr, call. = FALSE)
  }
  toks
}

#' Translate a macro expression into R code
#'
#' Recursive-descent parser producing R source that calls the helpers in
#' dyn_macro_env(). Precedence (low to high): ||, &&, in, comparisons,
#' range (:), additive, multiplicative, unary, power, postfix indexing.
#' @noRd
dyn_macro_translate <- function(expr) {
  toks <- dyn_macro_tokens(expr)
  pos <- 1L
  peek <- function() if (pos <= length(toks)) toks[pos] else ""
  take <- function(expected = NULL) {
    t <- peek()
    if (!is.null(expected) && t != expected) {
      stop("Expected '", expected, "' in macro expression: ", expr,
           call. = FALSE)
    }
    pos <<- pos + 1L
    t
  }
  p_or <- function() {
    x <- p_and()
    while (peek() == "||") { take(); x <- paste0("(", x, " || ", p_and(), ")") }
    x
  }
  p_and <- function() {
    x <- p_in()
    while (peek() == "&&") { take(); x <- paste0("(", x, " && ", p_in(), ")") }
    x
  }
  p_in <- function() {
    x <- p_cmp()
    while (peek() == "in") { take(); x <- paste0(".in(", x, ", ", p_cmp(), ")") }
    x
  }
  p_cmp <- function() {
    x <- p_range()
    while (peek() %in% c("==", "!=", "<", ">", "<=", ">=")) {
      op <- take()
      y <- p_range()
      x <- switch(op,
        "==" = paste0(".eq(", x, ", ", y, ")"),
        "!=" = paste0(".neq(", x, ", ", y, ")"),
        paste0("(", x, " ", op, " ", y, ")"))
    }
    x
  }
  p_range <- function() {
    x <- p_add()
    if (peek() == ":") {
      take()
      y <- p_add()
      if (peek() == ":") {
        take()
        z <- p_add()
        return(paste0(".range(", x, ", ", z, ", ", y, ")"))
      }
      x <- paste0(".range(", x, ", ", y, ")")
    }
    x
  }
  p_add <- function() {
    x <- p_mul()
    while (peek() %in% c("+", "-")) {
      op <- take()
      y <- p_mul()
      x <- if (op == "+") paste0(".add(", x, ", ", y, ")")
           else paste0("(", x, " - ", y, ")")
    }
    x
  }
  p_mul <- function() {
    x <- p_unary()
    while (peek() %in% c("*", "/")) {
      op <- take()
      x <- paste0("(", x, " ", op, " ", p_unary(), ")")
    }
    x
  }
  p_unary <- function() {
    if (peek() == "-") { take(); return(paste0("(-", p_unary(), ")")) }
    if (peek() == "+") { take(); return(p_unary()) }
    if (peek() == "!") { take(); return(paste0("(!", p_unary(), ")")) }
    p_pow()
  }
  p_pow <- function() {
    x <- p_postfix()
    if (peek() == "^") { take(); x <- paste0("(", x, " ^ ", p_unary(), ")") }
    x
  }
  p_postfix <- function() {
    x <- p_primary()
    while (peek() == "[") {
      take()
      i <- p_or()
      take("]")
      x <- paste0(".idx(", x, ", ", i, ")")
    }
    x
  }
  p_list <- function(close) {
    items <- character(0)
    if (peek() != close) {
      items <- p_or()
      while (peek() == ",") { take(); items <- c(items, p_or()) }
    }
    take(close)
    items
  }
  p_primary <- function() {
    t <- take()
    if (t == "(") {
      items <- p_list(")")
      if (length(items) == 1L) return(paste0("(", items, ")"))
      return(paste0(".arr(", paste(items, collapse = ", "), ")"))
    }
    if (t == "[") {
      return(paste0(".arr(", paste(p_list("]"), collapse = ", "), ")"))
    }
    if (grepl("^\"", t)) return(t)
    if (grepl("^[0-9.]", t)) return(t)
    if (grepl("^[A-Za-z_]", t)) {
      if (peek() == "(") {
        take()
        return(paste0(t, "(", paste(p_list(")"), collapse = ", "), ")"))
      }
      return(t)
    }
    stop("Unexpected '", t, "' in macro expression: ", expr, call. = FALSE)
  }
  out <- p_or()
  if (pos <= length(toks)) {
    stop("Unexpected '", toks[pos], "' in macro expression: ", expr,
         call. = FALSE)
  }
  out
}
