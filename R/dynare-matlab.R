# A small MATLAB interpreter for Dynare files
#
# Dynare model files often rely on MATLAB code: a `<model>_steadystate.m`
# file that computes the steady state (and sometimes calibrates parameters
# along the way), or MATLAB statements in the .mod file itself that set up
# the calibration (vectors, structures, loops). read_dynare() runs this code
# with the interpreter below. It covers the MATLAB language used in such
# files: scalars, matrices, strings, cell arrays and structures; indexing
# and assignment; control flow; subfunctions and anonymous functions;
# eval(); and the numerical solvers fsolve, csolve and fzero. Graphics and
# file input/output are ignored.

# ---------------------------------------------------------------------------
# Lexer
# ---------------------------------------------------------------------------

mat_keywords <- c("if", "elseif", "else", "end", "for", "parfor", "while",
                  "break", "continue", "return", "function", "switch",
                  "case", "otherwise", "try", "catch", "global",
                  "persistent")

mat_ops <- c("==", "~=", "!=", "<=", ">=", "&&", "||", ".*", "./", ".\\",
             ".^", ".'")

#' Tokenise MATLAB source
#' @return A list with character vectors `type` and `val` and a logical
#'   vector `sp` (whitespace before the token).
#' @noRd
mat_lex <- function(src) {
  ch <- strsplit(src, "", fixed = TRUE)[[1]]
  n <- length(ch)
  type <- character(0)
  val <- character(0)
  sp <- logical(0)
  stack <- character(0)
  space <- FALSE
  line_start <- TRUE
  i <- 1L
  top <- function() if (length(stack)) stack[length(stack)] else ""
  last_type <- function() if (length(type)) type[length(type)] else ""
  last_val <- function() if (length(val)) val[length(val)] else ""
  value_end <- function() {
    lt <- last_type()
    lt %in% c("num", "str", "id", "endidx") ||
      (lt == "op" && last_val() %in% c(")", "]", "}", "'", ".'")) ||
      (lt == "kw" && last_val() == "end")
  }
  emit <- function(t, v, starts_value = FALSE) {
    if (top() %in% c("[", "{") && space && starts_value && value_end()) {
      type[length(type) + 1L] <<- "op"
      val[length(val) + 1L] <<- ","
      sp[length(sp) + 1L] <<- TRUE
    }
    type[length(type) + 1L] <<- t
    val[length(val) + 1L] <<- v
    sp[length(sp) + 1L] <<- space
    space <<- FALSE
    line_start <<- FALSE
  }
  rest_of_line <- function(j) {
    k <- j
    while (k <= n && ch[k] != "\n") k <- k + 1L
    k
  }
  while (i <= n) {
    c1 <- ch[i]
    c2 <- if (i < n) ch[i + 1L] else ""
    if (c1 == "%" || (c1 == "#" && line_start)) {
      k <- rest_of_line(i)
      body <- trimws(paste(ch[i:(k - 1L)], collapse = ""))
      if (line_start && body %in% c("%{", "#{")) {
        # block comment up to a line holding only %}
        repeat {
          if (k > n) break
          j <- k + 1L
          k <- rest_of_line(j)
          if (j <= n && trimws(paste(ch[j:max(j, k - 1L)], collapse = "")) %in%
              c("%}", "#}")) break
        }
      }
      i <- k
      next
    }
    if (c1 == "." && c2 == "." && i + 2L <= n && ch[i + 2L] == ".") {
      i <- rest_of_line(i) + 1L
      space <- TRUE
      next
    }
    if (c1 == " " || c1 == "\t" || c1 == "\r") {
      space <- TRUE
      i <- i + 1L
      next
    }
    if (c1 == "\n") {
      if (top() %in% c("[", "{")) {
        emit("op", ";")
      } else if (top() != "(") {
        emit("sep", "\n")
      }
      line_start <- TRUE
      i <- i + 1L
      next
    }
    if (c1 == "'" && !(value_end() && (!space || !top() %in% c("[", "{")))) {
      j <- i + 1L
      buf <- character(0)
      repeat {
        if (j > n || ch[j] == "\n") stop("Unterminated string", call. = FALSE)
        if (ch[j] == "'") {
          if (j < n && ch[j + 1L] == "'") {
            buf <- c(buf, "'")
            j <- j + 2L
            next
          }
          break
        }
        buf <- c(buf, ch[j])
        j <- j + 1L
      }
      emit("str", paste(buf, collapse = ""), TRUE)
      i <- j + 1L
      next
    }
    if (c1 == "\"") {
      j <- i + 1L
      buf <- character(0)
      repeat {
        if (j > n || ch[j] == "\n") stop("Unterminated string", call. = FALSE)
        if (ch[j] == "\"") {
          if (j < n && ch[j + 1L] == "\"") {
            buf <- c(buf, "\"")
            j <- j + 2L
            next
          }
          break
        }
        buf <- c(buf, ch[j])
        j <- j + 1L
      }
      emit("str", paste(buf, collapse = ""), TRUE)
      i <- j + 1L
      next
    }
    if (grepl("[0-9]", c1) || (c1 == "." && grepl("[0-9]", c2))) {
      s <- paste(ch[i:min(n, i + 60L)], collapse = "")
      m <- regmatches(s, regexpr(
        "^([0-9]+\\.?[0-9]*|\\.[0-9]+)([eEdD][+-]?[0-9]+)?", s))
      emit("num", sub("[dD]", "e", m), TRUE)
      i <- i + nchar(m)
      # imaginary units are not supported; a trailing i/j is an error below
      next
    }
    if (grepl("[A-Za-z_]", c1)) {
      j <- i
      while (j < n && grepl("[A-Za-z0-9_]", ch[j + 1L])) j <- j + 1L
      w <- paste(ch[i:j], collapse = "")
      if (w %in% mat_keywords) {
        if (w == "end" && top() %in% c("(", "{")) {
          emit("endidx", w, TRUE)
        } else {
          emit("kw", w)
        }
      } else {
        emit("id", w, TRUE)
      }
      i <- j + 1L
      next
    }
    two <- paste0(c1, c2)
    if (two %in% mat_ops) {
      emit("op", two)
      i <- i + 2L
      next
    }
    if (c1 == "'" ) {
      emit("op", "'")
      i <- i + 1L
      next
    }
    if (c1 %in% c("(", "[", "{")) {
      emit("op", c1, TRUE)
      stack <- c(stack, c1)
      i <- i + 1L
      next
    }
    if (c1 %in% c(")", "]", "}")) {
      if (length(stack)) stack <- stack[-length(stack)]
      emit("op", c1)
      i <- i + 1L
      next
    }
    if (c1 %in% c("+", "-", "~", "!", "@")) {
      unary_start <- !grepl("[[:space:]]", c2) && c2 != "="
      emit("op", c1, unary_start)
      i <- i + 1L
      next
    }
    if (c1 %in% c("*", "/", "\\", "^", "<", ">", "&", "|", "=", ",", ";",
                  ":", ".")) {
      emit("op", c1)
      i <- i + 1L
      next
    }
    stop("Unexpected character '", c1, "' in MATLAB code", call. = FALSE)
  }
  list(type = type, val = val, sp = sp)
}

# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

#' Parse MATLAB source into a program (script statements and functions)
#' @noRd
mat_parse <- function(src) {
  tk <- mat_lex(src)
  n <- length(tk$type)
  pos <- 1L
  peek_t <- function(k = 0L) if (pos + k <= n) tk$type[pos + k] else "eof"
  peek_v <- function(k = 0L) if (pos + k <= n) tk$val[pos + k] else ""
  is_op <- function(v, k = 0L) peek_t(k) == "op" && peek_v(k) %in% v
  is_kw <- function(v, k = 0L) peek_t(k) == "kw" && peek_v(k) %in% v
  fail <- function(msg) {
    ctx <- paste(tk$val[max(1L, pos - 5L):min(n, pos + 5L)], collapse = " ")
    stop(msg, " near: ", ctx, call. = FALSE)
  }
  expect_op <- function(v) {
    if (!is_op(v)) fail(paste0("Expected '", v, "'"))
    pos <<- pos + 1L
  }
  skip_seps <- function() {
    while (peek_t() == "sep" || is_op(c(";", ","))) pos <<- pos + 1L
  }
  at_stmt_end <- function() {
    peek_t() %in% c("sep", "eof") || is_op(c(";", ",")) ||
      is_kw(c("end", "else", "elseif", "case", "otherwise", "catch"))
  }

  # --- expressions ---------------------------------------------------------
  parse_expr <- function() parse_oror()
  binary_level <- function(ops, nxt) {
    function() {
      a <- nxt()
      while (is_op(ops)) {
        op <- peek_v()
        pos <<- pos + 1L
        a <- list(t = "bin", op = op, a = a, b = nxt())
      }
      a
    }
  }
  parse_range <- function() {
    a <- parse_add()
    if (is_op(":") && !(peek_t(1L) == "op" && peek_v(1L) %in% c(")", ",", "}"))) {
      pos <<- pos + 1L
      b <- parse_add()
      if (is_op(":")) {
        pos <<- pos + 1L
        c3 <- parse_add()
        return(list(t = "range", a = a, s = b, b = c3))
      }
      return(list(t = "range", a = a, s = NULL, b = b))
    }
    a
  }
  parse_unary <- function() {
    if (is_op(c("-", "+", "~", "!"))) {
      op <- peek_v()
      pos <<- pos + 1L
      return(list(t = "un", op = op, a = parse_unary()))
    }
    parse_power()
  }
  parse_pow_operand <- function() {
    if (is_op(c("-", "+", "~", "!"))) {
      op <- peek_v()
      pos <<- pos + 1L
      return(list(t = "un", op = op, a = parse_pow_operand()))
    }
    parse_postfix()
  }
  parse_power <- function() {
    a <- parse_postfix()
    while (is_op(c("^", ".^"))) {
      op <- peek_v()
      pos <<- pos + 1L
      a <- list(t = "bin", op = op, a = a, b = parse_pow_operand())
    }
    a
  }
  parse_args <- function(close) {
    args <- list()
    if (is_op(close)) {
      pos <<- pos + 1L
      return(args)
    }
    repeat {
      if (is_op(":") && peek_t(1L) == "op" && peek_v(1L) %in% c(close, ",")) {
        pos <<- pos + 1L
        args[[length(args) + 1L]] <- list(t = "colon")
      } else {
        args[[length(args) + 1L]] <- parse_expr()
      }
      if (is_op(",")) {
        pos <<- pos + 1L
        next
      }
      expect_op(close)
      break
    }
    args
  }
  parse_postfix <- function() {
    a <- parse_primary()
    repeat {
      if (is_op("(") && !(tk$sp[pos] && in_matrix > 0L)) {
        pos <<- pos + 1L
        a <- list(t = "idx", obj = a, args = parse_args(")"), kind = "(")
      } else if (is_op("{") && !(tk$sp[pos] && in_matrix > 0L)) {
        pos <<- pos + 1L
        a <- list(t = "idx", obj = a, args = parse_args("}"), kind = "{")
      } else if (is_op(".") && peek_t(1L) == "id") {
        a <- list(t = "field", obj = a, name = peek_v(1L))
        pos <<- pos + 2L
      } else if (is_op(".") && peek_t(1L) == "op" && peek_v(1L) == "(") {
        pos <<- pos + 2L
        e <- parse_expr()
        expect_op(")")
        a <- list(t = "dynfield", obj = a, name = e)
      } else if (is_op(c("'", ".'"))) {
        pos <<- pos + 1L
        a <- list(t = "un", op = "'", a = a)
      } else {
        break
      }
    }
    a
  }
  in_matrix <- 0L
  parse_rows <- function(close) {
    in_matrix <<- in_matrix + 1L
    on.exit(in_matrix <<- in_matrix - 1L)
    rows <- list()
    row <- list()
    repeat {
      if (is_op(close)) {
        pos <<- pos + 1L
        break
      }
      if (peek_t() == "eof") fail("Unterminated matrix")
      if (is_op(";") || peek_t() == "sep") {
        pos <<- pos + 1L
        if (length(row)) rows[[length(rows) + 1L]] <- row
        row <- list()
        next
      }
      if (is_op(",")) {
        pos <<- pos + 1L
        next
      }
      row[[length(row) + 1L]] <- parse_expr()
    }
    if (length(row)) rows[[length(rows) + 1L]] <- row
    rows
  }
  parse_primary <- function() {
    t <- peek_t()
    v <- peek_v()
    if (t == "num") {
      pos <<- pos + 1L
      return(list(t = "num", v = as.numeric(v)))
    }
    if (t == "str") {
      pos <<- pos + 1L
      return(list(t = "str", v = v))
    }
    if (t == "id") {
      pos <<- pos + 1L
      return(list(t = "id", name = v))
    }
    if (t == "endidx") {
      pos <<- pos + 1L
      return(list(t = "end"))
    }
    if (t == "op" && v == "(") {
      pos <<- pos + 1L
      saved <- in_matrix
      in_matrix <<- 0L
      e <- parse_expr()
      in_matrix <<- saved
      expect_op(")")
      return(list(t = "paren", a = e))
    }
    if (t == "op" && v == "[") {
      pos <<- pos + 1L
      return(list(t = "matrix", rows = parse_rows("]")))
    }
    if (t == "op" && v == "{") {
      pos <<- pos + 1L
      return(list(t = "cell", rows = parse_rows("}")))
    }
    if (t == "op" && v == "@") {
      pos <<- pos + 1L
      if (is_op("(")) {
        pos <<- pos + 1L
        prm <- character(0)
        while (!is_op(")")) {
          if (peek_t() == "id") prm <- c(prm, peek_v())
          else if (!is_op(c(",", "~"))) fail("Bad anonymous function")
          if (is_op("~")) prm <- c(prm, "~")
          pos <<- pos + 1L
        }
        pos <<- pos + 1L
        saved <- in_matrix
        in_matrix <<- 0L
        body <- parse_expr()
        in_matrix <<- saved
        return(list(t = "anon", params = prm, body = body))
      }
      if (peek_t() == "id") {
        nm <- peek_v()
        pos <<- pos + 1L
        return(list(t = "fhandle", name = nm))
      }
      fail("Bad function handle")
    }
    if (t == "op" && v == ":") {
      pos <<- pos + 1L
      return(list(t = "colon"))
    }
    fail("Unexpected token")
  }
  parse_mul <- binary_level(c("*", "/", "\\", ".*", "./", ".\\"), parse_unary)
  parse_add <- binary_level(c("+", "-"), parse_mul)
  parse_cmp <- binary_level(c("==", "~=", "!=", "<", "<=", ">", ">="),
                            parse_range)
  parse_and <- binary_level("&", parse_cmp)
  parse_or <- binary_level("|", parse_and)
  parse_andand <- binary_level("&&", parse_or)
  parse_oror <- binary_level("||", parse_andand)

  # --- statements ----------------------------------------------------------
  command_words <- c("hold", "format", "close", "clear", "clc", "warning",
                     "more", "figure", "disp", "echo", "diary", "pkg",
                     "clearvars", "drawnow", "grid", "box", "axis", "shg",
                     "dbstop", "beep")
  parse_block <- function(stops) {
    body <- list()
    repeat {
      skip_seps()
      if (peek_t() == "eof" || is_kw(stops)) break
      body[[length(body) + 1L]] <- parse_stmt()
    }
    body
  }
  parse_stmt <- function() {
    t <- peek_t()
    v <- peek_v()
    if (t == "kw") {
      pos <<- pos + 1L
      if (v == "if") {
        conds <- list(parse_expr())
        blocks <- list(parse_block(c("elseif", "else", "end")))
        els <- NULL
        repeat {
          if (is_kw("elseif")) {
            pos <<- pos + 1L
            conds[[length(conds) + 1L]] <- parse_expr()
            blocks[[length(blocks) + 1L]] <- parse_block(c("elseif", "else",
                                                           "end"))
          } else if (is_kw("else")) {
            pos <<- pos + 1L
            els <- parse_block("end")
          } else {
            break
          }
        }
        if (!is_kw("end")) fail("'if' without 'end'")
        pos <<- pos + 1L
        return(list(t = "if", conds = conds, blocks = blocks, els = els))
      }
      if (v %in% c("for", "parfor")) {
        paren <- is_op("(")
        if (paren) pos <<- pos + 1L
        if (peek_t() != "id") fail("Bad for loop")
        var <- peek_v()
        pos <<- pos + 1L
        expect_op("=")
        e <- parse_expr()
        if (paren) {
          if (is_op(",")) {
            pos <<- pos + 1L
            parse_expr()
          }
          expect_op(")")
        }
        body <- parse_block("end")
        if (!is_kw("end")) fail("'for' without 'end'")
        pos <<- pos + 1L
        return(list(t = "for", var = var, e = e, body = body))
      }
      if (v == "while") {
        cond <- parse_expr()
        body <- parse_block("end")
        if (!is_kw("end")) fail("'while' without 'end'")
        pos <<- pos + 1L
        return(list(t = "while", cond = cond, body = body))
      }
      if (v == "switch") {
        e <- parse_expr()
        skip_seps()
        cases <- list()
        other <- NULL
        repeat {
          skip_seps()
          if (is_kw("case")) {
            pos <<- pos + 1L
            ce <- parse_expr()
            cases[[length(cases) + 1L]] <- list(e = ce, body = parse_block(
              c("case", "otherwise", "end")))
          } else if (is_kw("otherwise")) {
            pos <<- pos + 1L
            other <- parse_block("end")
          } else {
            break
          }
        }
        if (!is_kw("end")) fail("'switch' without 'end'")
        pos <<- pos + 1L
        return(list(t = "switch", e = e, cases = cases, other = other))
      }
      if (v == "try") {
        body <- parse_block(c("catch", "end"))
        cvar <- NULL
        cbody <- list()
        if (is_kw("catch")) {
          pos <<- pos + 1L
          if (peek_t() == "id" && peek_t(1L) %in% c("sep", "eof")) {
            cvar <- peek_v()
            pos <<- pos + 1L
          }
          cbody <- parse_block("end")
        }
        if (!is_kw("end")) fail("'try' without 'end'")
        pos <<- pos + 1L
        return(list(t = "try", body = body, cvar = cvar, cbody = cbody))
      }
      if (v %in% c("break", "continue", "return")) return(list(t = v))
      if (v %in% c("global", "persistent")) {
        nms <- character(0)
        while (peek_t() == "id") {
          nms <- c(nms, peek_v())
          pos <<- pos + 1L
        }
        return(list(t = "global", names = nms, persistent = v == "persistent"))
      }
      fail(paste0("Unexpected '", v, "'"))
    }
    # [a, b] = f(...)
    if (t == "op" && v == "[") {
      depth <- 0L
      k <- 0L
      repeat {
        if (pos + k > n) break
        if (tk$type[pos + k] == "op" && tk$val[pos + k] %in% c("[", "(", "{"))
          depth <- depth + 1L
        if (tk$type[pos + k] == "op" && tk$val[pos + k] %in% c("]", ")", "}"))
          depth <- depth - 1L
        if (depth == 0L) break
        k <- k + 1L
      }
      if (pos + k + 1L <= n && tk$type[pos + k + 1L] == "op" &&
          tk$val[pos + k + 1L] == "=") {
        pos <<- pos + 1L
        lhs <- list()
        in_matrix <<- in_matrix + 1L
        while (!is_op("]")) {
          if (is_op(",")) {
            pos <<- pos + 1L
            next
          }
          if (is_op("~")) {
            pos <<- pos + 1L
            lhs[[length(lhs) + 1L]] <- list(t = "ignore")
            next
          }
          lhs[[length(lhs) + 1L]] <- parse_postfix()
        }
        in_matrix <<- in_matrix - 1L
        pos <<- pos + 2L
        rhs <- parse_expr()
        return(list(t = "assign", lhs = lhs, rhs = rhs))
      }
    }
    # command syntax that does something: load file var1 var2
    if (t == "id" && v %in% c("load") && peek_t(1L) %in% c("id", "num") &&
        tk$sp[pos + 1L] && !(peek_t(2L) == "op" && peek_v(2L) %in% c("=", "(", "."))) {
      pos <<- pos + 1L
      words <- character(0)
      cur <- ""
      while (!peek_t() %in% c("sep", "eof") && !is_op(c(";", ","))) {
        if (tk$sp[pos] && nzchar(cur)) {
          words <- c(words, cur)
          cur <- ""
        }
        cur <- paste0(cur, peek_v())
        pos <<- pos + 1L
      }
      if (nzchar(cur)) words <- c(words, cur)
      return(list(t = "expr", e = list(t = "idx", obj = list(t = "id", name = v),
                                       args = lapply(words, function(w) list(t = "str", v = w)),
                                       kind = "(")))
    }
    if (t == "id" && v %in% command_words && peek_t(1L) %in% c("id", "num") &&
        tk$sp[pos + 1L] && !(peek_t(2L) == "op" && peek_v(2L) == "=")) {
      while (!peek_t() %in% c("sep", "eof") && !is_op(c(";", ","))) {
        pos <<- pos + 1L
      }
      return(list(t = "nop"))
    }
    e <- parse_expr()
    if (is_op("=")) {
      pos <<- pos + 1L
      rhs <- parse_expr()
      return(list(t = "assign", lhs = list(e), rhs = rhs))
    }
    if (!at_stmt_end()) fail("Unexpected token after expression")
    list(t = "expr", e = e)
  }

  # --- program -------------------------------------------------------------
  script <- list()
  funs <- list()
  repeat {
    skip_seps()
    if (peek_t() == "eof") break
    if (is_kw("function")) {
      pos <- pos + 1L
      outs <- character(0)
      # function [a, b] = name(args) | function a = name(args) | function name
      if (is_op("[")) {
        pos <- pos + 1L
        while (!is_op("]")) {
          if (peek_t() == "id") outs <- c(outs, peek_v())
          pos <- pos + 1L
        }
        pos <- pos + 1L
        expect_op("=")
      } else if (peek_t() == "id" && peek_t(1L) == "op" && peek_v(1L) == "=") {
        outs <- peek_v()
        pos <- pos + 2L
      }
      name <- peek_v()
      pos <- pos + 1L
      ins <- character(0)
      if (is_op("(")) {
        pos <- pos + 1L
        while (!is_op(")")) {
          if (peek_t() == "id") ins <- c(ins, peek_v())
          if (is_op("~")) ins <- c(ins, "~")
          pos <- pos + 1L
        }
        pos <- pos + 1L
      }
      body <- parse_block(c("end", "function"))
      if (is_kw("end")) pos <- pos + 1L
      funs[[name]] <- list(name = name, ins = ins, outs = outs, body = body)
    } else {
      script[[length(script) + 1L]] <- parse_stmt()
    }
  }
  list(script = script, funs = funs)
}

mat_parse_cache <- new.env(parent = emptyenv())

#' Parse with a cache (eval() strings are parsed once)
#' @noRd
mat_parse_cached <- function(src) {
  key <- src
  if (nchar(key) > 2000L) return(mat_parse(src))
  hit <- mat_parse_cache[[key]]
  if (!is.null(hit)) return(hit)
  out <- mat_parse(src)
  if (length(ls(mat_parse_cache)) > 5000L) {
    rm(list = ls(mat_parse_cache), envir = mat_parse_cache)
  }
  assign(key, out, envir = mat_parse_cache)
  out
}

# ---------------------------------------------------------------------------
# Values
# ---------------------------------------------------------------------------

mat_struct <- function(...) structure(list(...), class = "mat_struct")
mat_cell <- function(items, nrow = 1L, ncol = length(items)) {
  structure(matrix(items, nrow, ncol), class = "mat_cell")
}
mat_handle <- function(fn, label = "") structure(list(fn = fn, label = label),
                                                 class = "mat_handle")
mat_multi <- function(...) structure(list(...), class = "mat_multi")

mat_is_cell <- function(x) inherits(x, "mat_cell")
mat_is_struct <- function(x) inherits(x, "mat_struct")
mat_is_handle <- function(x) inherits(x, "mat_handle")
mat_is_char <- function(x) is.character(x) && !mat_is_cell(x)

#' Numeric (or logical) value as a matrix
#' @noRd
mat_m <- function(x) {
  if (is.null(x)) return(matrix(numeric(0), 0L, 0L))
  if (mat_is_char(x)) {
    x <- if (length(x) == 1L && nzchar(x)) utf8ToInt(x) else numeric(0)
    return(matrix(as.numeric(x), nrow = if (length(x)) 1L else 0L))
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("Expected a numeric value", call. = FALSE)
  }
  if (is.null(dim(x))) return(matrix(x, nrow = if (length(x)) 1L else 0L))
  x
}
mat_scalar <- function(x) {
  x <- mat_m(x)
  if (length(x) < 1L) stop("Expected a scalar, got an empty value",
                            call. = FALSE)
  as.numeric(x[1L])
}
mat_str <- function(x) {
  if (mat_is_char(x)) return(paste(x, collapse = ""))
  if (mat_is_cell(x) && length(x) == 1L) return(mat_str(x[[1L]]))
  if (is.numeric(x) || is.logical(x)) return(intToUtf8(as.integer(x)))
  stop("Expected a string", call. = FALSE)
}
mat_true <- function(x) {
  if (mat_is_char(x)) return(nzchar(x))
  x <- mat_m(x)
  length(x) > 0L && all(!is.na(x) & x != 0)
}
mat_size <- function(x) {
  if (mat_is_char(x)) return(c(if (nzchar(x)) 1L else 0L, nchar(x)))
  if (mat_is_struct(x) || mat_is_handle(x)) return(c(1L, 1L))
  if (is.null(dim(x))) return(c(if (length(x)) 1L else 0L, length(x)))
  dim(x)
}
mat_numel <- function(x) prod(mat_size(x))

# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------

mat_new_ctx <- function(funs = list(), globals = new.env(parent = emptyenv()),
                        parent = NULL) {
  ctx <- new.env(parent = emptyenv())
  ctx$vars <- new.env(parent = emptyenv())
  ctx$funs <- funs
  ctx$globals <- globals
  ctx$global_names <- character(0)
  ctx$endval <- NULL
  ctx$output <- character(0)
  ctx$depth <- if (is.null(parent)) 0L else parent$depth + 1L
  ctx$path <- if (is.null(parent)) character(0) else parent$path
  ctx
}

mat_signal <- function(kind) {
  structure(class = c(paste0("mat_", kind), "condition"),
            list(message = kind, call = NULL))
}

mat_run_block <- function(stmts, ctx) {
  for (s in stmts) mat_exec(s, ctx)
  invisible(NULL)
}

mat_exec <- function(s, ctx) {
  switch(s$t,
    assign = {
      nout <- length(s$lhs)
      val <- if (nout > 1L) mat_eval(s$rhs, ctx, nargout = nout)
             else mat_eval(s$rhs, ctx)
      if (nout > 1L) {
        if (!inherits(val, "mat_multi")) val <- mat_multi(val)
        for (k in seq_len(nout)) {
          if (s$lhs[[k]]$t == "ignore") next
          if (k > length(val)) stop("Too many output arguments", call. = FALSE)
          mat_assign(s$lhs[[k]], val[[k]], ctx)
        }
      } else {
        if (inherits(val, "mat_multi")) val <- val[[1L]]
        mat_assign(s$lhs[[1L]], val, ctx)
      }
    },
    expr = {
      e <- s$e
      # A bare identifier naming a function is a call with no outputs
      if (e$t == "id" && !exists(e$name, envir = ctx$vars, inherits = FALSE)) {
        mat_call(e$name, list(), ctx, nargout = 0L)
      } else if (e$t == "idx" && e$obj$t == "id" && e$kind == "(" &&
                 !exists(e$obj$name, envir = ctx$vars, inherits = FALSE)) {
        mat_call(e$obj$name, mat_eval_args(e$args, ctx), ctx, nargout = 0L)
      } else {
        v <- mat_eval(e, ctx)
        if (!inherits(v, "mat_multi")) assign("ans", v, envir = ctx$vars)
      }
    },
    `if` = {
      done <- FALSE
      for (k in seq_along(s$conds)) {
        if (mat_true(mat_eval(s$conds[[k]], ctx))) {
          mat_run_block(s$blocks[[k]], ctx)
          done <- TRUE
          break
        }
      }
      if (!done && !is.null(s$els)) mat_run_block(s$els, ctx)
    },
    `for` = {
      v <- mat_eval(s$e, ctx)
      if (mat_is_cell(v)) {
        cols <- lapply(seq_len(ncol(v)), function(j) {
          mat_cell(list(v[[1L, j]]), 1L, 1L)
        })
      } else if (mat_is_char(v)) {
        cols <- as.list(strsplit(v, "")[[1]])
      } else {
        v <- mat_m(v)
        cols <- lapply(seq_len(ncol(v)), function(j) v[, j, drop = FALSE])
      }
      for (cv in cols) {
        assign(s$var, cv, envir = ctx$vars)
        r <- tryCatch({
          mat_run_block(s$body, ctx)
          "ok"
        }, mat_break = function(e) "break", mat_continue = function(e) "cont")
        if (r == "break") break
      }
    },
    `while` = {
      iter <- 0L
      while (mat_true(mat_eval(s$cond, ctx))) {
        iter <- iter + 1L
        if (iter > 1e6) stop("while loop did not terminate", call. = FALSE)
        r <- tryCatch({
          mat_run_block(s$body, ctx)
          "ok"
        }, mat_break = function(e) "break", mat_continue = function(e) "cont")
        if (r == "break") break
      }
    },
    `switch` = {
      v <- mat_eval(s$e, ctx)
      matched <- FALSE
      for (cs in s$cases) {
        cv <- mat_eval(cs$e, ctx)
        opts <- if (mat_is_cell(cv)) as.list(cv) else list(cv)
        hit <- any(vapply(opts, function(o) {
          if (mat_is_char(v) || mat_is_char(o)) {
            mat_is_char(v) && mat_is_char(o) && identical(mat_str(v), mat_str(o))
          } else {
            isTRUE(mat_scalar(v) == mat_scalar(o))
          }
        }, logical(1)))
        if (hit) {
          mat_run_block(cs$body, ctx)
          matched <- TRUE
          break
        }
      }
      if (!matched && !is.null(s$other)) mat_run_block(s$other, ctx)
    },
    `try` = {
      tryCatch(mat_run_block(s$body, ctx), error = function(e) {
        if (!is.null(s$cvar)) {
          assign(s$cvar, mat_struct(message = conditionMessage(e),
                                    identifier = ""), envir = ctx$vars)
        }
        mat_run_block(s$cbody, ctx)
      })
    },
    `break` = signalCondition(mat_signal("break")),
    `continue` = signalCondition(mat_signal("continue")),
    `return` = signalCondition(mat_signal("return")),
    global = {
      for (nm in s$names) {
        ctx$global_names <- union(ctx$global_names, nm)
        if (exists(nm, envir = ctx$globals, inherits = FALSE)) {
          assign(nm, get(nm, envir = ctx$globals), envir = ctx$vars)
        } else if (!exists(nm, envir = ctx$vars, inherits = FALSE)) {
          assign(nm, matrix(numeric(0), 0L, 0L), envir = ctx$vars)
        }
      }
    },
    nop = NULL,
    stop("Unsupported statement", call. = FALSE)
  )
  invisible(NULL)
}

#' Evaluate call/index arguments (with `end` bound to the object's size)
#' @noRd
mat_eval_args <- function(args, ctx, obj = NULL) {
  nargs <- length(args)
  out <- vector("list", nargs)
  for (k in seq_len(nargs)) {
    a <- args[[k]]
    if (a$t == "colon") {
      out[k] <- list(structure(list(), class = "mat_colon"))
      next
    }
    old <- ctx$endval
    if (!is.null(obj)) {
      sz <- mat_size(obj)
      ctx$endval <- if (nargs == 1L) prod(sz)
                    else if (k <= length(sz)) {
                      if (k == nargs) prod(sz[k:length(sz)]) else sz[k]
                    } else 1L
    }
    v <- mat_eval(a, ctx)
    ctx$endval <- old
    if (inherits(v, "mat_csl")) {
      out[[k]] <- NULL
      out <- c(out[seq_len(k - 1L)], unclass(v), out[-seq_len(k)])
    } else {
      out[k] <- list(v)
    }
  }
  out
}

mat_eval <- function(e, ctx, nargout = 1L) {
  switch(e$t,
    num = matrix(e$v, 1L, 1L),
    str = e$v,
    paren = {
      v <- mat_eval(e$a, ctx)
      if (inherits(v, "mat_multi")) v[[1L]] else v
    },
    colon = structure(list(), class = "mat_colon"),
    end = {
      if (is.null(ctx$endval)) stop("'end' outside an index", call. = FALSE)
      matrix(ctx$endval, 1L, 1L)
    },
    id = {
      nm <- e$name
      if (exists(nm, envir = ctx$vars, inherits = FALSE)) {
        get(nm, envir = ctx$vars)
      } else {
        mat_first(mat_call(nm, list(), ctx, nargout = nargout), nargout)
      }
    },
    field = {
      obj <- mat_eval(e$obj, ctx)
      mat_get_field(obj, e$name)
    },
    dynfield = {
      obj <- mat_eval(e$obj, ctx)
      mat_get_field(obj, mat_str(mat_eval(e$name, ctx)))
    },
    idx = {
      if (e$obj$t == "id" &&
          !exists(e$obj$name, envir = ctx$vars, inherits = FALSE)) {
        if (e$kind == "{") stop("Undefined cell '", e$obj$name, "'",
                                call. = FALSE)
        args <- mat_eval_args(e$args, ctx)
        return(mat_first(mat_call(e$obj$name, args, ctx, nargout = nargout),
                         nargout))
      }
      obj <- mat_eval(e$obj, ctx)
      if (mat_is_handle(obj) && e$kind == "(") {
        args <- mat_eval_args(e$args, ctx)
        return(mat_first(obj$fn(args, nargout), nargout))
      }
      args <- mat_eval_args(e$args, ctx, obj)
      if (e$kind == "{") {
        items <- mat_index(obj, args, brace = TRUE)
        if (length(items) == 1L) return(items[[1L]])
        return(structure(items, class = "mat_csl"))
      }
      mat_index(obj, args)
    },
    matrix = mat_build_matrix(e$rows, ctx),
    cell = {
      rows <- lapply(e$rows, function(r) {
        lapply(r, function(x) {
          v <- mat_eval(x, ctx)
          if (inherits(v, "mat_multi")) v[[1L]] else v
        })
      })
      if (length(rows) == 0L) return(mat_cell(list(), 0L, 0L))
      # Concatenating cells inside {} nests them, as in MATLAB
      nc <- length(rows[[1L]])
      items <- unlist(lapply(rows, function(r) {
        if (length(r) != nc) stop("Inconsistent cell row lengths",
                                  call. = FALSE)
        r
      }), recursive = FALSE)
      mat_cell(items[as.vector(t(matrix(seq_along(items), nc)))],
               length(rows), nc)
    },
    range = {
      a <- mat_scalar(mat_eval(e$a, ctx))
      b <- mat_scalar(mat_eval(e$b, ctx))
      st <- if (is.null(e$s)) 1 else mat_scalar(mat_eval(e$s, ctx))
      if (st == 0 || (st > 0 && a > b) || (st < 0 && a < b)) {
        return(matrix(numeric(0), 1L, 0L))
      }
      k <- floor((b - a) / st + 1e-10)
      matrix(a + st * (0:k), nrow = 1L)
    },
    un = {
      v <- mat_eval(e$a, ctx)
      if (inherits(v, "mat_multi")) v <- v[[1L]]
      switch(e$op,
        `-` = -mat_m(v),
        `+` = mat_m(v),
        `~` = , `!` = {
          m <- mat_m(v)
          array(!(m != 0), dim(m))
        },
        `'` = {
          if (mat_is_cell(v)) mat_cell(as.list(t(unclass(v))), ncol(v), nrow(v))
          else t(mat_m(v))
        })
    },
    bin = mat_binary(e, ctx),
    anon = {
      captured <- as.list(ctx$vars, all.names = TRUE)
      funs <- ctx$funs
      globals <- ctx$globals
      prm <- e$params
      body <- e$body
      mat_handle(function(args, nargout = 1L) {
        c2 <- mat_new_ctx(funs, globals, ctx)
        for (nm in names(captured)) assign(nm, captured[[nm]], envir = c2$vars)
        for (k in seq_along(prm)) {
          if (prm[k] != "~" && k <= length(args)) {
            assign(prm[k], args[[k]], envir = c2$vars)
          }
        }
        mat_eval(body, c2, nargout = nargout)
      }, "anonymous")
    },
    fhandle = {
      nm <- e$name
      funs <- ctx$funs
      globals <- ctx$globals
      mat_handle(function(args, nargout = 1L) {
        c2 <- mat_new_ctx(funs, globals, ctx)
        mat_call(nm, args, c2, nargout = nargout)
      }, nm)
    },
    stop("Cannot evaluate expression of type ", e$t, call. = FALSE)
  )
}

mat_first <- function(v, nargout) {
  if (inherits(v, "mat_multi") && nargout <= 1L) {
    if (length(v) == 0L) stop("Function returned no value", call. = FALSE)
    return(v[[1L]])
  }
  v
}

mat_get_field <- function(obj, name) {
  if (!mat_is_struct(obj)) stop("Field access '", name, "' on a non-struct",
                                call. = FALSE)
  if (!name %in% names(obj)) stop("Reference to non-existent field '", name,
                                  "'", call. = FALSE)
  obj[[name]]
}

mat_build_matrix <- function(rows, ctx) {
  if (length(rows) == 0L) return(matrix(numeric(0), 0L, 0L))
  built <- lapply(rows, function(r) {
    vals <- list()
    for (x in r) {
      v <- mat_eval(x, ctx)
      if (inherits(v, "mat_multi")) v <- v[[1L]]
      if (inherits(v, "mat_csl")) vals <- c(vals, unclass(v))
      else vals[length(vals) + 1L] <- list(v)
    }
    mat_hcat(vals)
  })
  mat_vcat(built)
}

mat_hcat <- function(vals) {
  vals <- Filter(function(v) mat_is_char(v) || mat_is_struct(v) ||
                   mat_is_handle(v) || mat_numel(v) > 0L, vals)
  if (length(vals) == 0L) return(matrix(numeric(0), 0L, 0L))
  if (any(vapply(vals, mat_is_cell, logical(1)))) {
    items <- list()
    nr <- NA
    for (v in vals) {
      if (!mat_is_cell(v)) v <- mat_cell(list(v), 1L, 1L)
      items <- c(items, list(v))
    }
    nr <- nrow(items[[1L]])
    return(mat_cell(unlist(lapply(items, function(m) as.list(unclass(m))),
                           recursive = FALSE), nr,
                    sum(vapply(items, ncol, integer(1)))))
  }
  if (any(vapply(vals, mat_is_char, logical(1)))) {
    return(paste(vapply(vals, function(v) {
      if (mat_is_char(v)) v else intToUtf8(as.integer(round(mat_m(v))))
    }, character(1)), collapse = ""))
  }
  if (length(vals) == 1L) return(mat_m(vals[[1L]]))
  ms <- lapply(vals, mat_m)
  nr <- vapply(ms, nrow, integer(1))
  if (length(unique(nr)) != 1L) {
    stop("Dimensions of arrays being concatenated are not consistent",
         call. = FALSE)
  }
  do.call(cbind, ms)
}

mat_vcat <- function(vals) {
  vals <- Filter(function(v) mat_is_char(v) || mat_numel(v) > 0L, vals)
  if (length(vals) == 0L) return(matrix(numeric(0), 0L, 0L))
  if (length(vals) == 1L) return(vals[[1L]])
  if (any(vapply(vals, mat_is_cell, logical(1)))) {
    items <- lapply(vals, function(v) {
      if (mat_is_cell(v)) v else mat_cell(list(v), 1L, 1L)
    })
    nc <- ncol(items[[1L]])
    m <- do.call(rbind, lapply(items, unclass))
    return(structure(m, class = "mat_cell"))
  }
  if (any(vapply(vals, mat_is_char, logical(1)))) {
    # char matrix rows: keep as a cell column (enough for name lists)
    return(mat_cell(as.list(vapply(vals, mat_str, "")), length(vals), 1L))
  }
  ms <- lapply(vals, mat_m)
  nc <- vapply(ms, ncol, integer(1))
  if (length(unique(nc)) != 1L) {
    stop("Dimensions of arrays being concatenated are not consistent",
         call. = FALSE)
  }
  do.call(rbind, ms)
}

#' Element-wise operation with MATLAB's implicit expansion
#' @noRd
mat_elementwise <- function(a, b, f) {
  a <- mat_m(a)
  b <- mat_m(b)
  da <- dim(a)
  db <- dim(b)
  if (identical(da, db)) {
    out <- f(as.vector(a), as.vector(b))
    return(array(out, da))
  }
  if (length(a) == 1L) return(array(f(as.vector(a), as.vector(b)), db))
  if (length(b) == 1L) return(array(f(as.vector(a), as.vector(b)), da))
  nr <- max(da[1], db[1])
  nc <- max(da[2], db[2])
  ok <- function(d) (d[1] %in% c(1, nr)) && (d[2] %in% c(1, nc))
  if (!ok(da) || !ok(db)) {
    stop("Matrix dimensions must agree", call. = FALSE)
  }
  ex <- function(m) m[rep(seq_len(nrow(m)), length.out = nr),
                      rep(seq_len(ncol(m)), length.out = nc), drop = FALSE]
  array(f(as.vector(ex(a)), as.vector(ex(b))), c(nr, nc))
}

mat_binary <- function(e, ctx) {
  op <- e$op
  if (op %in% c("&&", "||")) {
    a <- mat_true(mat_eval(e$a, ctx))
    if (op == "&&" && !a) return(matrix(FALSE))
    if (op == "||" && a) return(matrix(TRUE))
    return(matrix(mat_true(mat_eval(e$b, ctx))))
  }
  a <- mat_eval(e$a, ctx)
  b <- mat_eval(e$b, ctx)
  if (inherits(a, "mat_multi")) a <- a[[1L]]
  if (inherits(b, "mat_multi")) b <- b[[1L]]
  if (op %in% c("==", "~=", "!=") && mat_is_char(a) && mat_is_char(b) &&
      nchar(a) != nchar(b)) {
    stop("Matrix dimensions must agree", call. = FALSE)
  }
  switch(op,
    `+` = mat_elementwise(a, b, `+`),
    `-` = mat_elementwise(a, b, `-`),
    `.*` = mat_elementwise(a, b, `*`),
    `./` = mat_elementwise(a, b, `/`),
    `.\\` = mat_elementwise(a, b, function(x, y) y / x),
    `.^` = mat_elementwise(a, b, mat_pow),
    `*` = {
      a <- mat_m(a)
      b <- mat_m(b)
      if (length(a) == 1L || length(b) == 1L) mat_elementwise(a, b, `*`)
      else a %*% b
    },
    `/` = {
      a <- mat_m(a)
      b <- mat_m(b)
      if (length(b) == 1L) mat_elementwise(a, b, `/`)
      else t(solve(t(b), t(a)))
    },
    `\\` = {
      a <- mat_m(a)
      b <- mat_m(b)
      if (length(a) == 1L) mat_elementwise(b, a, `/`)
      else if (nrow(a) == ncol(a)) solve(a, b)
      else qr.solve(a, b)
    },
    `^` = {
      a <- mat_m(a)
      b <- mat_m(b)
      if (length(a) == 1L && length(b) == 1L) return(mat_elementwise(a, b, mat_pow))
      if (length(b) == 1L && nrow(a) == ncol(a) && b == round(b)) {
        k <- as.integer(b)
        if (k < 0L) {
          a <- solve(a)
          k <- -k
        }
        out <- diag(nrow(a))
        while (k > 0L) {
          out <- out %*% a
          k <- k - 1L
        }
        return(out)
      }
      stop("Unsupported matrix power", call. = FALSE)
    },
    `==` = mat_elementwise(a, b, `==`),
    `~=` = , `!=` = mat_elementwise(a, b, `!=`),
    `<` = mat_elementwise(a, b, `<`),
    `<=` = mat_elementwise(a, b, `<=`),
    `>` = mat_elementwise(a, b, `>`),
    `>=` = mat_elementwise(a, b, `>=`),
    `&` = mat_elementwise(a, b, function(x, y) x != 0 & y != 0),
    `|` = mat_elementwise(a, b, function(x, y) x != 0 | y != 0),
    stop("Unsupported operator ", op, call. = FALSE)
  )
}

mat_pow <- function(x, y) {
  out <- x^y
  if (any(is.nan(out) & x < 0 & !is.nan(x) & !is.nan(y))) {
    stop("Negative number raised to a non-integer power (complex result)",
         call. = FALSE)
  }
  out
}

#' Resolve index arguments to integer positions for a dimension of size n
#' @noRd
mat_index_pos <- function(a, n) {
  if (inherits(a, "mat_colon")) return(seq_len(n))
  if (is.logical(a)) return(which(as.vector(a)))
  if (mat_is_char(a)) {
    if (a == ":") return(seq_len(n))
    stop("Invalid index", call. = FALSE)
  }
  v <- as.vector(mat_m(a))
  if (is.logical(v)) return(which(v))
  if (any(v < 1 | v != round(v))) {
    stop("Index must be a positive integer", call. = FALSE)
  }
  as.integer(v)
}

mat_index <- function(obj, args, brace = FALSE) {
  if (mat_is_struct(obj) || mat_is_handle(obj)) {
    ok <- all(vapply(args, function(a) inherits(a, "mat_colon") ||
                       isTRUE(all(mat_m(a) == 1)), logical(1)))
    if (!ok) stop("Index exceeds array bounds", call. = FALSE)
    return(obj)
  }
  if (mat_is_char(obj)) {
    chars <- strsplit(obj, "")[[1]]
    i <- mat_index_pos(args[[length(args)]], length(chars))
    if (any(i > length(chars))) stop("Index exceeds string length",
                                     call. = FALSE)
    return(paste(chars[i], collapse = ""))
  }
  cell <- mat_is_cell(obj)
  m <- if (cell) unclass(obj) else mat_m(obj)
  if (length(args) == 0L) return(obj)
  if (length(args) == 1L) {
    a <- args[[1L]]
    i <- mat_index_pos(a, length(m))
    if (any(i > length(m))) stop("Index exceeds array bounds (", max(i),
                                 " > ", length(m), ")", call. = FALSE)
    vals <- m[i]
    if (brace) return(as.list(vals))
    shape <- if (inherits(a, "mat_colon")) c(length(i), 1L)
             else if (nrow(m) == 1L) c(1L, length(i))
             else if (ncol(m) == 1L) c(length(i), 1L)
             else if (!is.null(dim(a)) && length(i) == length(a)) dim(a)
             else c(1L, length(i))
    if (cell) return(structure(matrix(vals, shape[1], shape[2]),
                               class = "mat_cell"))
    return(array(vals, shape))
  }
  if (length(args) > 2L) {
    extra <- args[-(1:2)]
    if (!all(vapply(extra, function(a) inherits(a, "mat_colon") ||
                      isTRUE(all(mat_m(a) == 1)), logical(1)))) {
      stop("Arrays with more than two dimensions are not supported",
           call. = FALSE)
    }
  }
  i <- mat_index_pos(args[[1L]], nrow(m))
  j <- mat_index_pos(args[[2L]], ncol(m))
  if (any(i > nrow(m)) || any(j > ncol(m))) {
    stop("Index exceeds matrix dimensions", call. = FALSE)
  }
  out <- m[i, j, drop = FALSE]
  if (brace) return(as.list(out))
  if (cell) return(structure(out, class = "mat_cell"))
  out
}

mat_assign <- function(lhs, value, ctx) {
  switch(lhs$t,
    id = {
      assign(lhs$name, value, envir = ctx$vars)
      if (lhs$name %in% ctx$global_names) {
        assign(lhs$name, value, envir = ctx$globals)
      }
    },
    field = , dynfield = {
      name <- if (lhs$t == "field") lhs$name
              else mat_str(mat_eval(lhs$name, ctx))
      base <- mat_lvalue_current(lhs$obj, ctx)
      if (is.null(base) || (!mat_is_struct(base) && mat_numel(base) == 0L)) {
        base <- mat_struct()
      }
      if (!mat_is_struct(base)) stop("Field assignment to a non-struct",
                                     call. = FALSE)
      base[name] <- list(value)
      mat_assign(lhs$obj, base, ctx)
    },
    idx = {
      base <- mat_lvalue_current(lhs$obj, ctx)
      args <- mat_eval_args(lhs$args, ctx, base)
      mat_assign(lhs$obj, mat_assign_index(base, args, value,
                                           brace = lhs$kind == "{"), ctx)
    },
    paren = mat_assign(lhs$a, value, ctx),
    stop("Invalid assignment target", call. = FALSE)
  )
}

mat_lvalue_current <- function(node, ctx) {
  switch(node$t,
    id = if (exists(node$name, envir = ctx$vars, inherits = FALSE)) {
      get(node$name, envir = ctx$vars)
    } else NULL,
    field = , dynfield = {
      base <- mat_lvalue_current(node$obj, ctx)
      name <- if (node$t == "field") node$name
              else mat_str(mat_eval(node$name, ctx))
      if (mat_is_struct(base) && name %in% names(base)) base[[name]] else NULL
    },
    idx = {
      base <- mat_lvalue_current(node$obj, ctx)
      if (is.null(base)) return(NULL)
      args <- mat_eval_args(node$args, ctx, base)
      tryCatch({
        if (node$kind == "{") mat_index(base, args, brace = TRUE)[[1L]]
        else mat_index(base, args)
      }, error = function(e) NULL)
    },
    paren = mat_lvalue_current(node$a, ctx),
    stop("Invalid assignment target", call. = FALSE)
  )
}

mat_assign_index <- function(base, args, value, brace = FALSE) {
  cell <- brace || mat_is_cell(base) || mat_is_cell(value) && !brace &&
    is.null(base)
  if (brace || mat_is_cell(base)) {
    m <- if (is.null(base) || !mat_is_cell(base)) matrix(list(), 0L, 0L)
         else unclass(base)
    vals <- if (brace) list(value)
            else if (mat_is_cell(value)) as.list(unclass(value))
            else list(value)
    m <- mat_assign_positions(m, args, vals, fill = list(NULL))
    return(structure(m, class = "mat_cell"))
  }
  if (mat_is_char(base) && mat_is_char(value)) {
    chars <- strsplit(base, "")[[1]]
    i <- mat_index_pos(args[[length(args)]], length(chars))
    chars[i] <- strsplit(value, "")[[1]]
    return(paste(chars, collapse = ""))
  }
  m <- if (is.null(base)) matrix(numeric(0), 0L, 0L) else mat_m(base)
  v <- mat_m(value)
  if (length(v) == 0L && !is.null(dim(v)) && all(dim(v) == 0L)) {
    # deletion: x(i) = []
    if (length(args) == 1L) {
      i <- mat_index_pos(args[[1L]], length(m))
      keep <- setdiff(seq_along(m), i)
      return(if (ncol(m) == 1L) matrix(m[keep], ncol = 1L)
             else matrix(m[keep], nrow = 1L))
    }
    if (inherits(args[[1L]], "mat_colon")) {
      return(m[, -mat_index_pos(args[[2L]], ncol(m)), drop = FALSE])
    }
    return(m[-mat_index_pos(args[[1L]], nrow(m)), , drop = FALSE])
  }
  mat_assign_positions(m, args, as.vector(v), fill = 0)
}

mat_assign_positions <- function(m, args, vals, fill) {
  if (length(args) == 1L) {
    i <- mat_index_pos(args[[1L]], length(m))
    need <- if (length(i)) max(i) else 0L
    if (need > length(m)) {
      if (nrow(m) <= 1L) {
        grown <- if (is.list(m)) matrix(rep(list(NULL), need), 1L, need)
                 else matrix(fill, 1L, need)
        if (length(m)) grown[seq_along(m)] <- m
      } else if (ncol(m) == 1L) {
        grown <- if (is.list(m)) matrix(rep(list(NULL), need), need, 1L)
                 else matrix(fill, need, 1L)
        grown[seq_along(m)] <- m
      } else {
        stop("Cannot grow a matrix with a linear index", call. = FALSE)
      }
      m <- grown
    }
    if (length(vals) != 1L && length(vals) != length(i)) {
      stop("Assignment dimension mismatch", call. = FALSE)
    }
    m[i] <- if (is.list(m)) vals else rep(vals, length.out = length(i))
    return(m)
  }
  i <- mat_index_pos(args[[1L]], nrow(m))
  j <- mat_index_pos(args[[2L]], ncol(m))
  nr <- max(nrow(m), if (length(i)) max(i) else 0L)
  nc <- max(ncol(m), if (length(j)) max(j) else 0L)
  if (nr > nrow(m) || nc > ncol(m)) {
    grown <- if (is.list(m)) matrix(rep(list(NULL), nr * nc), nr, nc)
             else matrix(fill, nr, nc)
    if (length(m)) grown[seq_len(nrow(m)), seq_len(ncol(m))] <- m
    m <- grown
  }
  if (inherits(args[[1L]], "mat_colon") && nrow(m) == 0L) {
    i <- seq_len(max(1L, length(vals) %/% max(1L, length(j))))
  }
  m[i, j] <- if (is.list(m)) vals else rep(vals, length.out = length(i) * length(j))
  m
}

# ---------------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------------

mat_call <- function(name, args, ctx, nargout = 1L) {
  f <- ctx$funs[[name]]
  if (!is.null(f)) return(mat_call_user(f, args, ctx, nargout))
  b <- mat_builtins[[name]]
  if (!is.null(b)) {
    return(b(args, nargout, ctx))
  }
  if (name %in% mat_ignored) return(mat_multi())
  # a function in its own file next to the model
  for (d in ctx$path) {
    fp <- file.path(d, paste0(name, ".m"))
    if (file.exists(fp)) {
      prog <- mat_parse(paste(readLines(fp, warn = FALSE), collapse = "\n"))
      if (length(prog$funs) == 0L) break
      main <- prog$funs[[1L]]
      main$name <- name
      prog$funs[[1L]] <- main
      names(prog$funs)[1L] <- name
      for (nm in names(prog$funs)) {
        if (is.null(ctx$funs[[nm]])) ctx$funs[[nm]] <- prog$funs[[nm]]
      }
      return(mat_call_user(main, args, ctx, nargout))
    }
  }
  stop("Undefined function or variable '", name, "'", call. = FALSE)
}

mat_call_user <- function(f, args, ctx, nargout = 1L) {
  if (ctx$depth > 200L) stop("Maximum recursion depth", call. = FALSE)
  c2 <- mat_new_ctx(ctx$funs, ctx$globals, ctx)
  ins <- f$ins
  var_in <- length(ins) && ins[length(ins)] == "varargin"
  fixed <- if (var_in) ins[-length(ins)] else ins
  if (length(args) > length(fixed) && !var_in) {
    stop("Too many input arguments to ", f$name, call. = FALSE)
  }
  for (k in seq_along(fixed)) {
    if (k <= length(args) && fixed[k] != "~") {
      assign(fixed[k], args[[k]], envir = c2$vars)
    }
  }
  if (var_in) {
    extra <- if (length(args) > length(fixed)) args[-seq_along(fixed)] else list()
    assign("varargin", mat_cell(extra, 1L, length(extra)), envir = c2$vars)
  }
  assign("nargin", matrix(length(args)), envir = c2$vars)
  assign("nargout", matrix(nargout), envir = c2$vars)
  tryCatch(mat_run_block(f$body, c2), mat_return = function(e) NULL)
  outs <- f$outs
  var_out <- length(outs) && outs[length(outs)] == "varargout"
  res <- list()
  for (o in outs) {
    if (o == "varargout") {
      if (exists("varargout", envir = c2$vars, inherits = FALSE)) {
        res <- c(res, as.list(unclass(get("varargout", envir = c2$vars))))
      }
    } else if (exists(o, envir = c2$vars, inherits = FALSE)) {
      res[length(res) + 1L] <- list(get(o, envir = c2$vars))
    } else {
      if (length(res) < max(1L, nargout)) {
        stop("Output '", o, "' of ", f$name, " not assigned", call. = FALSE)
      }
      break
    }
  }
  do.call(mat_multi, res)
}

mat_ignored <- c(
  "figure", "subplot", "plot", "hold", "title", "xlabel", "ylabel", "zlabel",
  "legend", "axis", "grid", "close", "clc", "print", "saveas", "drawnow",
  "pause", "set", "box", "text", "line", "bar", "area", "fill", "semilogy",
  "semilogx", "loglog", "hist", "histogram", "surf", "mesh", "colorbar",
  "colormap", "annotation", "sgtitle", "suptitle", "xlim", "ylim", "zlim",
  "linkaxes", "tiledlayout", "nexttile", "orient", "axes", "diary", "more",
  "beep", "format", "tic", "keyboard", "dbstop", "rng", "addpath",
  "rmpath", "set_dynare_seed", "fclose", "fflush", "shg", "pkg",
  "warning_off", "clearvars", "clear")

mat_ret <- function(...) do.call(mat_multi, list(...))
mat_num_arg <- function(args, k, default = NULL) {
  if (length(args) < k) return(default)
  mat_m(args[[k]])
}

mat_reduce <- function(x, f, dim = NULL) {
  x <- mat_m(x)
  if (is.null(dim)) dim <- if (nrow(x) == 1L) 2L else 1L
  if (length(x) == 0L) return(matrix(f(numeric(0)), 1L, 1L))
  if (dim == 1L) matrix(apply(x, 2L, f), nrow = 1L)
  else matrix(apply(x, 1L, f), ncol = 1L)
}

mat_cumulate <- function(x, f) {
  x <- mat_m(x)
  if (nrow(x) == 1L || ncol(x) == 1L) return(array(f(as.vector(x)), dim(x)))
  apply(x, 2L, f)
}

mat_minmax <- function(args, nargout, which_fn, cmp) {
  if (length(args) >= 2L && mat_numel(args[[2L]]) > 0L) {
    return(mat_elementwise(args[[1L]], args[[2L]], cmp))
  }
  x <- mat_m(args[[1L]])
  if (length(x) == 0L) return(mat_ret(matrix(numeric(0), 0L, 0L),
                                      matrix(numeric(0), 0L, 0L)))
  if (nrow(x) == 1L || ncol(x) == 1L) {
    v <- as.vector(x)
    k <- which_fn(v)
    if (length(k) == 0L) k <- 1L
    return(mat_ret(matrix(v[k]), matrix(k)))
  }
  k <- apply(x, 2L, function(col) {
    r <- which_fn(col)
    if (length(r)) r else 1L
  })
  mat_ret(matrix(x[cbind(k, seq_len(ncol(x)))], nrow = 1L), matrix(k, nrow = 1L))
}

#' MATLAB sprintf: the format is applied repeatedly to all arguments
#' @noRd
mat_sprintf <- function(fmt, args) {
  fmt <- gsub("\\\\n", "\n", fmt)
  fmt <- gsub("\\\\t", "\t", fmt)
  fmt <- gsub("\\\\\\\\", "\\\\", fmt)
  vals <- list()
  for (a in args) {
    if (mat_is_char(a)) vals[[length(vals) + 1L]] <- a
    else if (mat_is_cell(a)) vals <- c(vals, as.list(unclass(a)))
    else vals <- c(vals, as.list(as.vector(mat_m(a))))
  }
  specs <- regmatches(fmt, gregexpr("%[-+ 0#]*[0-9]*(\\.[0-9]+)?[diouxXfeEgGcs]",
                                    fmt))[[1]]
  if (length(specs) == 0L || length(vals) == 0L) {
    return(gsub("%%", "%", fmt, fixed = TRUE))
  }
  pieces <- strsplit(fmt, "%[-+ 0#]*[0-9]*(\\.[0-9]+)?[diouxXfeEgGcs]",
                     perl = TRUE)[[1]]
  out <- character(0)
  k <- 1L
  repeat {
    for (s in seq_along(specs)) {
      if (s <= length(pieces)) out <- c(out, gsub("%%", "%", pieces[s]))
      if (k > length(vals)) break
      v <- vals[[k]]
      k <- k + 1L
      conv <- substr(specs[s], nchar(specs[s]), nchar(specs[s]))
      sp <- specs[s]
      if (conv %in% c("d", "i", "u")) {
        if (is.numeric(v) && v == round(v)) {
          out <- c(out, sprintf(sub("[diu]$", "d", sp), as.integer(v)))
        } else if (is.numeric(v)) {
          out <- c(out, sprintf(sub("[diu]$", "e", sp), v))
        } else {
          out <- c(out, v)
        }
      } else if (conv == "s") {
        out <- c(out, sprintf(sp, if (is.character(v)) v else format(v)))
      } else if (conv == "c") {
        out <- c(out, if (is.character(v)) v else intToUtf8(v))
      } else {
        out <- c(out, sprintf(sp, as.numeric(v)))
      }
    }
    if (k > length(vals)) break
  }
  if (length(pieces) > length(specs) && k > length(vals)) {
    out <- c(out, gsub("%%", "%", pieces[length(pieces)]))
  }
  paste(out, collapse = "")
}

mat_num2str <- function(x, fmt = NULL) {
  if (mat_is_char(x)) return(x)
  if (!is.null(fmt)) {
    if (mat_is_char(fmt)) return(mat_sprintf(fmt, list(x)))
    return(paste(formatC(as.vector(mat_m(x)), digits = mat_scalar(fmt),
                         format = "g"), collapse = "  "))
  }
  v <- as.vector(mat_m(x))
  paste(vapply(v, function(z) {
    if (is.na(z)) "NaN" else if (z == round(z) && abs(z) < 1e15) {
      format(z, scientific = FALSE)
    } else trimws(formatC(z, digits = 5L, format = "g"))
  }, ""), collapse = "  ")
}

mat_optimset <- function(args) {
  s <- mat_struct()
  k <- 1L
  if (length(args) %% 2L == 1L) k <- 2L
  while (k < length(args)) {
    s[[mat_str(args[[k]])]] <- args[[k + 1L]]
    k <- k + 2L
  }
  s
}

mat_opt <- function(opts, names, default) {
  if (!mat_is_struct(opts)) return(default)
  for (nm in names) {
    if (nm %in% names(opts) && mat_numel(opts[[nm]]) > 0L &&
        !mat_is_char(opts[[nm]])) {
      return(mat_scalar(opts[[nm]]))
    }
  }
  default
}

#' Newton's method with a numerical Jacobian, falling back to
#' Levenberg-Marquardt steps when a full step does not reduce the residual
#' @noRd
mat_newton <- function(fn, x0, tol = 1e-10, maxit = 500L) {
  x <- as.vector(x0)
  n <- length(x)
  f <- fn(x)
  if (any(!is.finite(f))) return(list(x = x, f = f, ok = FALSE))
  ssr <- sum(f^2)
  mu <- 1e-6
  for (it in seq_len(maxit)) {
    if (max(abs(f)) < tol) return(list(x = x, f = f, ok = TRUE))
    J <- matrix(0, length(f), n)
    for (j in seq_len(n)) {
      h <- 1e-7 * max(1, abs(x[j]))
      xp <- x
      xp[j] <- x[j] + h
      J[, j] <- (fn(xp) - f) / h
    }
    step <- tryCatch(-qr.solve(J, f), error = function(e) NULL)
    improved <- FALSE
    if (!is.null(step) && all(is.finite(step))) {
      lam <- 1
      for (ls in 1:30) {
        xn <- x + lam * step
        fn_x <- tryCatch(fn(xn), error = function(e) rep(NaN, length(f)))
        if (all(is.finite(fn_x)) && sum(fn_x^2) < ssr) {
          x <- xn
          f <- fn_x
          ssr <- sum(f^2)
          improved <- TRUE
          break
        }
        lam <- lam / 2
      }
    }
    if (!improved) {
      g <- crossprod(J, f)
      H <- crossprod(J)
      for (k in 1:40) {
        st <- tryCatch(-solve(H + mu * diag(max(1, max(diag(H))), n), g),
                       error = function(e) NULL)
        if (!is.null(st)) {
          xn <- x + as.vector(st)
          fn_x <- tryCatch(fn(xn), error = function(e) rep(NaN, length(f)))
          if (all(is.finite(fn_x)) && sum(fn_x^2) < ssr) {
            x <- xn
            f <- fn_x
            ssr <- sum(f^2)
            mu <- mu / 10
            improved <- TRUE
            break
          }
        }
        mu <- mu * 10
      }
    }
    if (!improved) break
  }
  list(x = x, f = f, ok = max(abs(f)) < max(tol, 1e-8))
}

mat_call_handle <- function(h, x, extra = list(), shape = NULL) {
  arg <- if (is.null(shape)) matrix(x, ncol = 1L) else array(x, shape)
  v <- if (mat_is_handle(h)) h$fn(c(list(arg), extra), 1L)
       else stop("Expected a function handle", call. = FALSE)
  if (inherits(v, "mat_multi")) v <- v[[1L]]
  as.vector(mat_m(v))
}

mat_handle_arg <- function(h, ctx) {
  if (mat_is_handle(h)) return(h)
  if (mat_is_char(h)) {
    nm <- h
    return(mat_handle(function(args, nargout = 1L) {
      mat_call(nm, args, ctx, nargout)
    }, nm))
  }
  stop("Expected a function handle", call. = FALSE)
}

mat_builtins <- list(
  pi = function(a, n, ctx) matrix(pi),
  eps = function(a, n, ctx) {
    if (length(a)) matrix(.Machine$double.eps * 2^floor(log2(abs(mat_scalar(a[[1]])))))
    else matrix(.Machine$double.eps)
  },
  `Inf` = function(a, n, ctx) mat_fill(Inf, a),
  inf = function(a, n, ctx) mat_fill(Inf, a),
  `NaN` = function(a, n, ctx) mat_fill(NaN, a),
  nan = function(a, n, ctx) mat_fill(NaN, a),
  `NA` = function(a, n, ctx) mat_fill(NaN, a),
  zeros = function(a, n, ctx) mat_fill(0, a),
  ones = function(a, n, ctx) mat_fill(1, a),
  true = function(a, n, ctx) array(TRUE, dim(mat_fill(1, a))),
  false = function(a, n, ctx) array(FALSE, dim(mat_fill(0, a))),
  realmax = function(a, n, ctx) matrix(.Machine$double.xmax),
  realmin = function(a, n, ctx) matrix(.Machine$double.xmin),
  eye = function(a, n, ctx) {
    d <- dim(mat_fill(0, a))
    m <- matrix(0, d[1], d[2])
    if (min(d) > 0) diag(m) <- 1
    m
  },
  abs = function(a, n, ctx) abs(mat_m(a[[1]])),
  sign = function(a, n, ctx) sign(mat_m(a[[1]])),
  sqrt = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    if (any(x < 0, na.rm = TRUE)) stop("sqrt of a negative number",
                                      call. = FALSE)
    sqrt(x)
  },
  exp = function(a, n, ctx) exp(mat_m(a[[1]])),
  log = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    if (any(x < 0, na.rm = TRUE)) stop("log of a negative number",
                                      call. = FALSE)
    log(x)
  },
  log10 = function(a, n, ctx) log10(mat_m(a[[1]])),
  log2 = function(a, n, ctx) log2(mat_m(a[[1]])),
  log1p = function(a, n, ctx) log1p(mat_m(a[[1]])),
  expm1 = function(a, n, ctx) expm1(mat_m(a[[1]])),
  sin = function(a, n, ctx) sin(mat_m(a[[1]])),
  cos = function(a, n, ctx) cos(mat_m(a[[1]])),
  tan = function(a, n, ctx) tan(mat_m(a[[1]])),
  asin = function(a, n, ctx) asin(mat_m(a[[1]])),
  acos = function(a, n, ctx) acos(mat_m(a[[1]])),
  atan = function(a, n, ctx) atan(mat_m(a[[1]])),
  atan2 = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], atan2),
  sinh = function(a, n, ctx) sinh(mat_m(a[[1]])),
  cosh = function(a, n, ctx) cosh(mat_m(a[[1]])),
  tanh = function(a, n, ctx) tanh(mat_m(a[[1]])),
  gamma = function(a, n, ctx) gamma(mat_m(a[[1]])),
  gammaln = function(a, n, ctx) lgamma(mat_m(a[[1]])),
  beta = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], beta),
  betaln = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], lbeta),
  erf = function(a, n, ctx) 2 * stats::pnorm(mat_m(a[[1]]) * sqrt(2)) - 1,
  erfc = function(a, n, ctx) 2 * stats::pnorm(-mat_m(a[[1]]) * sqrt(2)),
  erfinv = function(a, n, ctx) stats::qnorm((mat_m(a[[1]]) + 1) / 2) / sqrt(2),
  normcdf = function(a, n, ctx) {
    mu <- if (length(a) > 1L) mat_m(a[[2]]) else 0
    s <- if (length(a) > 2L) mat_m(a[[3]]) else 1
    array(stats::pnorm(mat_m(a[[1]]), mu, s), dim(mat_m(a[[1]])))
  },
  normpdf = function(a, n, ctx) {
    mu <- if (length(a) > 1L) mat_m(a[[2]]) else 0
    s <- if (length(a) > 2L) mat_m(a[[3]]) else 1
    array(stats::dnorm(mat_m(a[[1]]), mu, s), dim(mat_m(a[[1]])))
  },
  norminv = function(a, n, ctx) {
    mu <- if (length(a) > 1L) mat_m(a[[2]]) else 0
    s <- if (length(a) > 2L) mat_m(a[[3]]) else 1
    array(stats::qnorm(mat_m(a[[1]]), mu, s), dim(mat_m(a[[1]])))
  },
  floor = function(a, n, ctx) floor(mat_m(a[[1]])),
  ceil = function(a, n, ctx) ceiling(mat_m(a[[1]])),
  round = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    if (length(a) > 1L) {
      d <- mat_scalar(a[[2]])
      return(sign(x) * floor(abs(x) * 10^d + 0.5) / 10^d)
    }
    sign(x) * floor(abs(x) + 0.5)
  },
  fix = function(a, n, ctx) trunc(mat_m(a[[1]])),
  mod = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], function(x, y) {
    ifelse(y == 0, x, x - floor(x / y) * y)
  }),
  rem = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], function(x, y) {
    ifelse(y == 0, x, x - trunc(x / y) * y)
  }),
  hypot = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], function(x, y) sqrt(x^2 + y^2)),
  power = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], mat_pow),
  nthroot = function(a, n, ctx) mat_elementwise(a[[1]], a[[2]], function(x, y) {
    sign(x) * abs(x)^(1 / y)
  }),
  real = function(a, n, ctx) mat_m(a[[1]]),
  conj = function(a, n, ctx) mat_m(a[[1]]),
  imag = function(a, n, ctx) 0 * mat_m(a[[1]]),
  isreal = function(a, n, ctx) matrix(TRUE),
  isnan = function(a, n, ctx) is.nan(mat_m(a[[1]])) | is.na(mat_m(a[[1]])),
  isinf = function(a, n, ctx) is.infinite(mat_m(a[[1]])),
  isfinite = function(a, n, ctx) is.finite(mat_m(a[[1]])),
  isempty = function(a, n, ctx) matrix(mat_numel(a[[1]]) == 0L),
  isnumeric = function(a, n, ctx) matrix(is.numeric(a[[1]]) && !mat_is_cell(a[[1]])),
  islogical = function(a, n, ctx) matrix(is.logical(a[[1]])),
  ischar = function(a, n, ctx) matrix(mat_is_char(a[[1]])),
  iscell = function(a, n, ctx) matrix(mat_is_cell(a[[1]])),
  isstruct = function(a, n, ctx) matrix(mat_is_struct(a[[1]])),
  isscalar = function(a, n, ctx) matrix(mat_numel(a[[1]]) == 1L),
  isvector = function(a, n, ctx) matrix(min(mat_size(a[[1]])) == 1L),
  iscellstr = function(a, n, ctx) matrix(mat_is_cell(a[[1]]) &&
    all(vapply(unclass(a[[1]]), mat_is_char, logical(1)))),
  isa = function(a, n, ctx) {
    cl <- mat_str(a[[2]])
    x <- a[[1]]
    matrix(switch(cl, double = , numeric = , float = is.numeric(x) &&
                    !mat_is_cell(x),
                  char = mat_is_char(x), cell = mat_is_cell(x),
                  struct = mat_is_struct(x), logical = is.logical(x),
                  function_handle = mat_is_handle(x), FALSE))
  },
  class = function(a, n, ctx) {
    x <- a[[1]]
    if (mat_is_char(x)) "char" else if (mat_is_cell(x)) "cell"
    else if (mat_is_struct(x)) "struct" else if (mat_is_handle(x))
      "function_handle" else if (is.logical(x)) "logical" else "double"
  },
  numel = function(a, n, ctx) matrix(mat_numel(a[[1]])),
  length = function(a, n, ctx) {
    s <- mat_size(a[[1]])
    matrix(if (prod(s) == 0) 0 else max(s))
  },
  ndims = function(a, n, ctx) matrix(2),
  size = function(a, n, ctx) {
    s <- mat_size(a[[1]])
    if (length(a) > 1L) {
      k <- mat_scalar(a[[2]])
      return(matrix(if (k <= 2) s[k] else 1))
    }
    if (n <= 1L) return(mat_ret(matrix(s, nrow = 1L)))
    do.call(mat_multi, c(lapply(s, function(z) matrix(z)),
                         rep(list(matrix(1)), max(0L, n - 2L))))
  },
  sum = function(a, n, ctx) {
    d <- if (length(a) > 1L) mat_scalar(a[[2]]) else NULL
    mat_reduce(a[[1]], sum, d)
  },
  prod = function(a, n, ctx) mat_reduce(a[[1]], prod,
                                        if (length(a) > 1L) mat_scalar(a[[2]])),
  mean = function(a, n, ctx) mat_reduce(a[[1]], mean,
                                        if (length(a) > 1L) mat_scalar(a[[2]])),
  median = function(a, n, ctx) mat_reduce(a[[1]], stats::median),
  std = function(a, n, ctx) mat_reduce(a[[1]], stats::sd),
  var = function(a, n, ctx) mat_reduce(a[[1]], stats::var),
  any = function(a, n, ctx) mat_reduce(a[[1]], function(v) any(v != 0 & !is.na(v))),
  all = function(a, n, ctx) mat_reduce(a[[1]], function(v) all(v != 0)),
  nnz = function(a, n, ctx) matrix(sum(mat_m(a[[1]]) != 0)),
  cumsum = function(a, n, ctx) mat_cumulate(a[[1]], cumsum),
  cumprod = function(a, n, ctx) mat_cumulate(a[[1]], cumprod),
  diff = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    if (nrow(x) == 1L) matrix(diff(as.vector(x)), nrow = 1L)
    else apply(x, 2L, diff)
  },
  max = function(a, n, ctx) mat_minmax(a, n, which.max, pmax),
  min = function(a, n, ctx) mat_minmax(a, n, which.min, pmin),
  find = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    k <- which(as.vector(x) != 0)
    if (length(a) > 1L) k <- utils::head(k, mat_scalar(a[[2]]))
    if (nrow(x) == 1L) matrix(k, nrow = 1L) else matrix(k, ncol = 1L)
  },
  linspace = function(a, n, ctx) {
    k <- if (length(a) > 2L) mat_scalar(a[[3]]) else 100
    matrix(seq(mat_scalar(a[[1]]), mat_scalar(a[[2]]), length.out = k), nrow = 1L)
  },
  repmat = function(a, n, ctx) {
    x <- a[[1]]
    r <- mat_m(a[[2]])
    c2 <- if (length(a) > 2L) mat_scalar(a[[3]]) else if (length(r) > 1L) r[2] else r[1]
    r <- r[1]
    if (mat_is_cell(x)) {
      m <- unclass(x)
      return(structure(m[rep(seq_len(nrow(m)), r), rep(seq_len(ncol(m)), c2),
                         drop = FALSE], class = "mat_cell"))
    }
    x <- mat_m(x)
    x[rep(seq_len(nrow(x)), r), rep(seq_len(ncol(x)), c2), drop = FALSE]
  },
  reshape = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    d <- if (length(a) == 2L) as.vector(mat_m(a[[2]]))
         else vapply(a[-1], function(z) if (mat_numel(z)) mat_scalar(z) else NA, 0)
    if (anyNA(d)) d[is.na(d)] <- length(x) / prod(d[!is.na(d)])
    array(as.vector(x), d[1:2])
  },
  diag = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    if (nrow(x) == 1L || ncol(x) == 1L) {
      m <- matrix(0, length(x), length(x))
      diag(m) <- as.vector(x)
      return(m)
    }
    matrix(diag(x), ncol = 1L)
  },
  inv = function(a, n, ctx) solve(mat_m(a[[1]])),
  pinv = function(a, n, ctx) MASS_ginv(mat_m(a[[1]])),
  det = function(a, n, ctx) matrix(det(mat_m(a[[1]]))),
  trace = function(a, n, ctx) matrix(sum(diag(mat_m(a[[1]])))),
  transpose = function(a, n, ctx) t(mat_m(a[[1]])),
  kron = function(a, n, ctx) kronecker(mat_m(a[[1]]), mat_m(a[[2]])),
  chol = function(a, n, ctx) chol(mat_m(a[[1]])),
  norm = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    p <- if (length(a) > 1L) a[[2]] else 2
    if (nrow(x) == 1L || ncol(x) == 1L) {
      if (mat_is_char(p) && p == "inf" || is.numeric(p) && is.infinite(mat_scalar(p))) {
        return(matrix(max(abs(x))))
      }
      return(matrix(sum(abs(x)^mat_scalar(p))^(1 / mat_scalar(p))))
    }
    matrix(norm(x, "2"))
  },
  sort = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    desc <- length(a) > 1L && mat_is_char(a[[2]]) && tolower(a[[2]]) == "descend"
    o <- order(as.vector(x), decreasing = desc)
    mat_ret(array(as.vector(x)[o], dim(x)), array(o, dim(x)))
  },
  unique = function(a, n, ctx) {
    x <- a[[1]]
    if (mat_is_cell(x)) {
      v <- sort(unique(vapply(unclass(x), mat_str, "")))
      return(mat_cell(as.list(v), length(v), 1L))
    }
    matrix(sort(unique(as.vector(mat_m(x)))), nrow = 1L)
  },
  fliplr = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    x[, rev(seq_len(ncol(x))), drop = FALSE]
  },
  flipud = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    x[rev(seq_len(nrow(x))), , drop = FALSE]
  },
  horzcat = function(a, n, ctx) mat_hcat(a),
  vertcat = function(a, n, ctx) mat_vcat(a),
  cat = function(a, n, ctx) {
    if (mat_scalar(a[[1]]) == 1) mat_vcat(a[-1]) else mat_hcat(a[-1])
  },
  logical = function(a, n, ctx) mat_m(a[[1]]) != 0,
  double = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    storage.mode(x) <- "double"
    x
  },
  char = function(a, n, ctx) {
    x <- a[[1]]
    if (mat_is_char(x)) x
    else if (mat_is_cell(x)) mat_vcat(as.list(unclass(x)))
    else intToUtf8(as.integer(mat_m(x)))
  },
  num2str = function(a, n, ctx) mat_num2str(a[[1]], if (length(a) > 1L) a[[2]]),
  int2str = function(a, n, ctx) mat_num2str(round(mat_m(a[[1]]))),
  mat2str = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    if (length(x) == 1L) return(mat_num2str(x))
    rows <- apply(x, 1L, function(r) paste(vapply(r, mat_num2str, ""), collapse = " "))
    paste0("[", paste(rows, collapse = ";"), "]")
  },
  str2num = function(a, n, ctx) matrix(as.numeric(mat_str(a[[1]]))),
  str2double = function(a, n, ctx) matrix(suppressWarnings(as.numeric(mat_str(a[[1]])))),
  sprintf = function(a, n, ctx) mat_sprintf(mat_str(a[[1]]), a[-1]),
  fprintf = function(a, n, ctx) {
    if (length(a) >= 1L && !mat_is_char(a[[1]])) a <- a[-1]
    if (length(a)) ctx$output <- c(ctx$output, mat_sprintf(mat_str(a[[1]]), a[-1]))
    mat_multi()
  },
  disp = function(a, n, ctx) mat_multi(),
  display = function(a, n, ctx) mat_multi(),
  error = function(a, n, ctx) {
    if (length(a) == 0L) return(mat_multi())
    msg <- mat_str(a[[1]])
    rest <- a[-1]
    if (length(rest) && grepl(":", msg) && !grepl("\\s", msg) && mat_is_char(rest[[1]])) {
      msg <- mat_str(rest[[1]])
      rest <- rest[-1]
    }
    stop(if (length(rest)) mat_sprintf(msg, rest) else gsub("\\\\n", "\n", msg),
         call. = FALSE)
  },
  warning = function(a, n, ctx) mat_multi(),
  assert = function(a, n, ctx) {
    if (!mat_true(a[[1]])) stop(if (length(a) > 1L) mat_str(a[[2]]) else "Assertion failed",
                               call. = FALSE)
    mat_multi()
  },
  strcmp = function(a, n, ctx) mat_strcmp(a[[1]], a[[2]], identity),
  strcmpi = function(a, n, ctx) mat_strcmp(a[[1]], a[[2]], tolower),
  strncmp = function(a, n, ctx) {
    k <- mat_scalar(a[[3]])
    mat_strcmp(a[[1]], a[[2]], function(s) substr(s, 1L, k))
  },
  strrep = function(a, n, ctx) gsub(mat_str(a[[2]]), mat_str(a[[3]]), mat_str(a[[1]]),
                                    fixed = TRUE),
  strtrim = function(a, n, ctx) mat_map_str(a[[1]], trimws),
  deblank = function(a, n, ctx) mat_map_str(a[[1]], function(s) sub("\\s+$", "", s)),
  upper = function(a, n, ctx) mat_map_str(a[[1]], toupper),
  lower = function(a, n, ctx) mat_map_str(a[[1]], tolower),
  strcat = function(a, n, ctx) paste(vapply(a, function(x) sub("\\s+$", "", mat_str(x)), ""),
                                     collapse = ""),
  regexprep = function(a, n, ctx) gsub(mat_str(a[[2]]), mat_str(a[[3]]), mat_str(a[[1]]),
                                       perl = TRUE),
  strmatch = function(a, n, ctx) {
    s <- mat_str(a[[1]])
    lst <- a[[2]]
    names_v <- if (mat_is_cell(lst)) vapply(unclass(lst), mat_str, "") else mat_str(lst)
    exact <- length(a) > 2L
    k <- if (exact) which(sub("\\s+$", "", names_v) == s) else which(startsWith(names_v, s))
    matrix(k, ncol = 1L)
  },
  ismember = function(a, n, ctx) {
    x <- a[[1]]
    s <- a[[2]]
    if (mat_is_char(x) || mat_is_cell(x)) {
      xs <- if (mat_is_cell(x)) vapply(unclass(x), mat_str, "") else mat_str(x)
      ss <- if (mat_is_cell(s)) vapply(unclass(s), mat_str, "") else mat_str(s)
      return(mat_ret(matrix(xs %in% ss, nrow = 1L),
                     matrix(match(xs, ss, nomatch = 0L), nrow = 1L)))
    }
    xv <- mat_m(x)
    sv <- as.vector(mat_m(s))
    mat_ret(array(as.vector(xv) %in% sv, dim(xv)),
            array(match(as.vector(xv), sv, nomatch = 0L), dim(xv)))
  },
  isfield = function(a, n, ctx) {
    s <- a[[1]]
    f <- a[[2]]
    fs <- if (mat_is_cell(f)) vapply(unclass(f), mat_str, "") else mat_str(f)
    matrix(mat_is_struct(s) & fs %in% names(s), nrow = 1L)
  },
  fieldnames = function(a, n, ctx) {
    nm <- names(a[[1]])
    mat_cell(as.list(nm), length(nm), 1L)
  },
  rmfield = function(a, n, ctx) {
    s <- a[[1]]
    s[[mat_str(a[[2]])]] <- NULL
    s
  },
  struct = function(a, n, ctx) {
    s <- mat_struct()
    k <- 1L
    while (k < length(a)) {
      v <- a[[k + 1L]]
      if (mat_is_cell(v) && length(v) == 1L) v <- v[[1L]]
      s[mat_str(a[[k]])] <- list(v)
      k <- k + 2L
    }
    s
  },
  getfield = function(a, n, ctx) mat_get_field(a[[1]], mat_str(a[[2]])),
  setfield = function(a, n, ctx) {
    s <- a[[1]]
    s[mat_str(a[[2]])] <- list(a[[3]])
    s
  },
  cell = function(a, n, ctx) {
    d <- dim(mat_fill(0, a))
    mat_cell(rep(list(matrix(numeric(0), 0L, 0L)), prod(d)), d[1], d[2])
  },
  num2cell = function(a, n, ctx) {
    x <- mat_m(a[[1]])
    mat_cell(lapply(as.vector(x), function(z) matrix(z)), nrow(x), ncol(x))
  },
  cellfun = function(a, n, ctx) {
    f <- a[[1]]
    cl <- a[[2]]
    uniform <- TRUE
    if (length(a) >= 4L && tolower(mat_str(a[[3]])) == "uniformoutput") {
      uniform <- mat_true(a[[4]])
    }
    fn <- if (mat_is_char(f)) {
      switch(f, isempty = function(x) matrix(mat_numel(x) == 0L),
             length = function(x) matrix(max(mat_size(x))),
             function(x) mat_first(mat_call(f, list(x), ctx), 1L))
    } else function(x) mat_first(f$fn(list(x), 1L), 1L)
    res <- lapply(unclass(cl), fn)
    if (uniform) array(vapply(res, mat_scalar, 0), dim(unclass(cl)))
    else structure(array(res, dim(unclass(cl))), class = "mat_cell")
  },
  arrayfun = function(a, n, ctx) {
    f <- a[[1]]
    x <- mat_m(a[[2]])
    array(vapply(as.vector(x), function(z) {
      mat_scalar(mat_first(f$fn(list(matrix(z)), 1L), 1L))
    }, 0), dim(x))
  },
  deal = function(a, n, ctx) {
    if (length(a) == 1L) return(do.call(mat_multi, rep(a, max(1L, n))))
    do.call(mat_multi, a)
  },
  feval = function(a, n, ctx) {
    f <- mat_handle_arg(a[[1]], ctx)
    f$fn(a[-1], n)
  },
  func2str = function(a, n, ctx) a[[1]]$label,
  str2func = function(a, n, ctx) mat_handle_arg(mat_str(a[[1]]), ctx),
  exist = function(a, n, ctx) {
    nm <- mat_str(a[[1]])
    matrix(if (exists(nm, envir = ctx$vars, inherits = FALSE)) 1
           else if (!is.null(ctx$funs[[nm]])) 2
           else if (!is.null(mat_builtins[[nm]])) 5 else 0)
  },
  isoctave = function(a, n, ctx) matrix(FALSE),
  is_octave = function(a, n, ctx) matrix(FALSE),
  user_has_matlab_license = function(a, n, ctx) matrix(FALSE),
  user_has_octave_forge_package = function(a, n, ctx) matrix(FALSE),
  verLessThan = function(a, n, ctx) matrix(FALSE),
  matlab_ver_less_than = function(a, n, ctx) matrix(FALSE),
  octave_ver_less_than = function(a, n, ctx) matrix(FALSE),
  dynare_version = function(a, n, ctx) "6.0",
  dyn_ver = function(a, n, ctx) "6.0",
  toc = function(a, n, ctx) matrix(0),
  optimset = function(a, n, ctx) mat_optimset(a),
  optimoptions = function(a, n, ctx) mat_optimset(a),
  fsolve = function(a, n, ctx) {
    h <- mat_handle_arg(a[[1]], ctx)
    x0 <- mat_m(a[[2]])
    tol <- mat_opt(if (length(a) > 2L) a[[3]], c("TolFun", "FunctionTolerance"),
                   1e-10)
    r <- mat_newton(function(x) mat_call_handle(h, x, shape = dim(x0)),
                    as.vector(x0), tol = min(tol, 1e-10))
    mat_ret(array(r$x, dim(x0)), matrix(r$f, ncol = 1L),
            matrix(if (r$ok) 1 else -2), mat_struct(iterations = matrix(0)))
  },
  csolve = function(a, n, ctx) {
    h <- mat_handle_arg(a[[1]], ctx)
    x0 <- mat_m(a[[2]])
    crit <- if (length(a) > 3L && mat_numel(a[[4]])) mat_scalar(a[[4]]) else 1e-10
    extra <- if (length(a) > 5L) a[-(1:5)] else list()
    r <- mat_newton(function(x) mat_call_handle(h, x, extra, shape = dim(x0)),
                    as.vector(x0), tol = min(crit, 1e-10))
    mat_ret(array(r$x, dim(x0)), matrix(if (r$ok) 0 else 4))
  },
  fzero = function(a, n, ctx) {
    h <- mat_handle_arg(a[[1]], ctx)
    x0 <- as.vector(mat_m(a[[2]]))
    f <- function(x) mat_call_handle(h, x)
    if (length(x0) >= 2L) {
      lo <- x0[1]
      hi <- x0[2]
    } else {
      dx <- if (x0 == 0) 1 / 50 else abs(x0) / 50
      lo <- hi <- x0
      f0 <- f(x0)
      if (f0 == 0) return(mat_ret(matrix(x0), matrix(0), matrix(1)))
      found <- FALSE
      for (k in 1:200) {
        lo <- x0 - dx
        hi <- x0 + dx
        flo <- tryCatch(f(lo), error = function(e) NaN)
        fhi <- tryCatch(f(hi), error = function(e) NaN)
        if (is.finite(flo) && sign(flo) != sign(f0)) {
          hi <- x0
          found <- TRUE
          break
        }
        if (is.finite(fhi) && sign(fhi) != sign(f0)) {
          lo <- x0
          found <- TRUE
          break
        }
        dx <- dx * sqrt(2)
      }
      if (!found) return(mat_ret(matrix(NaN), matrix(NaN), matrix(-6)))
    }
    r <- stats::uniroot(f, c(lo, hi), tol = 1e-15, maxiter = 1000L)
    mat_ret(matrix(r$root), matrix(r$f.root), matrix(1))
  },
  fminsearch = function(a, n, ctx) {
    h <- mat_handle_arg(a[[1]], ctx)
    x0 <- mat_m(a[[2]])
    f <- function(x) mat_call_handle(h, x, shape = dim(x0))[1]
    r <- if (length(x0) == 1L) {
      stats::optim(as.vector(x0), f, method = "BFGS")
    } else {
      stats::optim(as.vector(x0), f, method = "Nelder-Mead",
                   control = list(maxit = 20000, reltol = 1e-14))
    }
    mat_ret(array(r$par, dim(x0)), matrix(r$value), matrix(1))
  },
  fminbnd = function(a, n, ctx) {
    h <- mat_handle_arg(a[[1]], ctx)
    r <- stats::optimize(function(x) mat_call_handle(h, x)[1],
                         c(mat_scalar(a[[2]]), mat_scalar(a[[3]])), tol = 1e-12)
    mat_ret(matrix(r$minimum), matrix(r$objective), matrix(1))
  },
  eval = function(a, n, ctx) {
    code <- mat_str(a[[1]])
    run <- function(code) {
      prog <- mat_parse_cached(code)
      if (n >= 1L && length(prog$script) == 1L && prog$script[[1]]$t == "expr") {
        return(mat_eval(prog$script[[1]]$e, ctx, nargout = n))
      }
      mat_run_block(prog$script, ctx)
      mat_multi()
    }
    if (length(a) > 1L) {
      return(tryCatch(run(code), error = function(e) run(mat_str(a[[2]]))))
    }
    run(code)
  },
  get_param_by_name = function(a, n, ctx) {
    M <- get0("M_", envir = ctx$globals, inherits = FALSE)
    k <- match(mat_str(a[[1]]), vapply(unclass(M$param_names), mat_str, ""))
    matrix(M$params[k])
  },
  set_param_value = function(a, n, ctx) {
    M <- get0("M_", envir = ctx$globals, inherits = FALSE)
    k <- match(mat_str(a[[1]]), vapply(unclass(M$param_names), mat_str, ""))
    M$params[k] <- mat_scalar(a[[2]])
    assign("M_", M, envir = ctx$globals)
    if (exists("M_", envir = ctx$vars, inherits = FALSE)) {
      assign("M_", M, envir = ctx$vars)
    }
    mat_multi()
  }
)

MASS_ginv <- function(X) {
  s <- svd(X)
  pos <- s$d > max(1e-12 * s$d[1], 0)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

mat_fill <- function(v, a) {
  d <- vapply(a, function(z) {
    if (mat_is_char(z)) NA_real_ else as.numeric(mat_m(z))[1]
  }, 0)
  if (length(a) == 1L && !mat_is_char(a[[1]]) && mat_numel(a[[1]]) > 1L) {
    d <- as.numeric(mat_m(a[[1]]))[1:2]
  }
  d <- d[!is.na(d)]
  if (length(d) == 0L) d <- c(1, 1)
  if (length(d) == 1L) d <- c(d, d)
  matrix(v, max(0, d[1]), max(0, d[2]))
}

mat_strcmp <- function(a, b, f) {
  one <- function(x) if (mat_is_char(x)) f(x) else NA_character_
  if (mat_is_cell(a) || mat_is_cell(b)) {
    if (mat_is_cell(a) && !mat_is_cell(b)) {
      tmp <- a
      a <- b
      b <- tmp
    }
    bs <- vapply(unclass(b), function(x) {
      r <- one(x)
      if (is.na(r)) "\001" else r
    }, "")
    as_ <- if (mat_is_cell(a)) vapply(unclass(a), one, "") else one(a)
    return(array(!is.na(as_) & as_ == bs, dim(unclass(b))))
  }
  matrix(mat_is_char(a) && mat_is_char(b) && identical(f(a), f(b)))
}

mat_map_str <- function(x, f) {
  if (mat_is_cell(x)) {
    return(structure(array(lapply(unclass(x), function(s) f(mat_str(s))),
                           dim(unclass(x))), class = "mat_cell"))
  }
  f(mat_str(x))
}

# ---------------------------------------------------------------------------
# Dynare interface
# ---------------------------------------------------------------------------

#' Run MATLAB script code in a context (errors are returned, not thrown)
#' @noRd
mat_run_script <- function(code, ctx) {
  prog <- mat_parse(code)
  for (nm in names(prog$funs)) ctx$funs[[nm]] <- prog$funs[[nm]]
  tryCatch(mat_run_block(prog$script, ctx), mat_return = function(e) NULL)
  invisible(ctx)
}

#' Dynare's M_ structure for the steady-state file
#' @noRd
mat_dynare_M <- function(fname, param_names, params, endo, exo) {
  col_cell <- function(x) mat_cell(as.list(x), length(x), 1L)
  mat_struct(
    fname = fname,
    param_names = col_cell(param_names),
    param_names_long = col_cell(param_names),
    param_names_tex = col_cell(param_names),
    params = matrix(as.numeric(params), ncol = 1L),
    param_nbr = matrix(length(param_names)),
    endo_names = col_cell(endo),
    endo_names_long = col_cell(endo),
    endo_names_tex = col_cell(endo),
    endo_nbr = matrix(length(endo)),
    orig_endo_nbr = matrix(length(endo)),
    exo_names = col_cell(exo),
    exo_names_long = col_cell(exo),
    exo_nbr = matrix(length(exo)),
    exo_det_nbr = matrix(0),
    maximum_lag = matrix(1),
    maximum_lead = matrix(1),
    Sigma_e = diag(1, max(1L, length(exo)))
  )
}

mat_dynare_options <- function() {
  mat_struct(
    order = matrix(1), linear = matrix(0), debug = matrix(0),
    qz_criterium = matrix(1 + 1e-6), steadystate_flag = matrix(1),
    solve_tolf = matrix(.Machine$double.eps^(1 / 3)),
    solve_tolx = matrix(.Machine$double.eps^(2 / 3)),
    solve_algo = matrix(4), jacobian_flag = matrix(1), gstep = matrix(1e-2),
    noprint = matrix(1), nograph = matrix(1), periods = matrix(0),
    irf = matrix(40), replic = matrix(50), loglinear = matrix(0),
    steady = mat_struct(maxit = matrix(50)), dataset = mat_struct(),
    steadystate = mat_struct(nocheck = matrix(0)), homotopy_mode = matrix(0),
    bytecode = matrix(0), block = matrix(0), use_dll = matrix(0),
    TeX = matrix(0), ramsey_policy = matrix(0), nocorr = matrix(1),
    nomoments = matrix(1), simul = mat_struct(maxit = matrix(50))
  )
}

#' A steady-state function from a Dynare `<model>_steadystate.m` file
#'
#' Returns function(params, ys0, exo) giving list(ys, params, check); params
#' are in declaration order and include any the file recalibrates.
#' @noRd
dyn_matlab_steady_state <- function(path, param_names, endo, exo) {
  code <- paste(readLines(path, warn = FALSE), collapse = "\n")
  prog <- tryCatch(mat_parse(code), error = function(e) {
    stop("Cannot parse the steady-state file ", basename(path), ": ",
         conditionMessage(e), call. = FALSE)
  })
  if (length(prog$funs) == 0L) {
    stop("The steady-state file ", basename(path), " defines no function.",
         call. = FALSE)
  }
  main <- prog$funs[[1L]]
  fname <- sub("_steadystate$", "", main$name)
  function(params, ys0, exo_vals) {
    globals <- new.env(parent = emptyenv())
    M <- mat_dynare_M(fname, param_names, params, endo, exo)
    opts <- mat_dynare_options()
    assign("M_", M, envir = globals)
    assign("options_", opts, envir = globals)
    assign("oo_", mat_struct(steady_state = matrix(ys0, ncol = 1L),
                             exo_steady_state = matrix(exo_vals, ncol = 1L)),
           envir = globals)
    ctx <- mat_new_ctx(prog$funs, globals)
    ctx$path <- dirname(path)
    args <- list(matrix(as.numeric(ys0), ncol = 1L),
                 matrix(as.numeric(exo_vals), ncol = 1L), M, opts)
    out <- tryCatch(
      mat_call_user(main, args[seq_len(min(length(args), length(main$ins)))],
                    ctx, nargout = length(main$outs)),
      error = function(e) {
        stop("Steady-state file ", basename(path), ": ", conditionMessage(e),
             call. = FALSE)
      })
    res <- stats::setNames(as.list(out), main$outs[seq_along(out)])
    ys <- as.numeric(mat_m(res[[1L]]))
    p_out <- if ("params" %in% names(res)) as.numeric(mat_m(res$params))
             else as.numeric(mat_m(get("M_", envir = globals)$params))
    check <- if ("check" %in% names(res)) mat_scalar(res$check)
             else if (length(res) >= 2L && !"params" %in% names(res)[2]) {
               mat_scalar(res[[2L]])
             } else 0
    if (length(ys) < length(endo)) ys <- c(ys, rep(0, length(endo) - length(ys)))
    list(ys = stats::setNames(ys[seq_along(endo)], endo),
         params = stats::setNames(p_out, param_names), check = check)
  }
}

# ---------------------------------------------------------------------------
# Data files, optimisation and statistics
# ---------------------------------------------------------------------------

#' Locate a data file next to the model (or in the working directory)
#' @noRd
mat_find_file <- function(name, ctx, exts = character(0)) {
  dirs <- unique(c(ctx$path, getwd()))
  cands <- name
  if (!grepl("\\.[A-Za-z0-9]+$", name)) cands <- c(paste0(name, exts), name)
  for (d in dirs) for (f in cands) {
    fp <- if (grepl("^(/|[A-Za-z]:)", f)) f else file.path(d, f)
    if (file.exists(fp) && !dir.exists(fp)) return(fp)
  }
  stop("File not found: ", name, call. = FALSE)
}

#' Convert values read by R.matlab::readMat to interpreter values
#' @noRd
mat_from_readmat <- function(x) {
  if (is.character(x)) return(paste(x, collapse = ""))
  if (is.numeric(x) || is.logical(x)) {
    if (is.null(dim(x))) return(matrix(as.numeric(x), nrow = 1L))
    return(array(as.numeric(x), dim(x)[1:2]))
  }
  if (is.list(x)) {
    d <- dim(x)
    fields <- if (!is.null(dimnames(x))) dimnames(x)[[1]] else names(x)
    if (!is.null(fields) && length(fields) == length(x)) {
      s <- mat_struct()
      for (k in seq_along(fields)) {
        s[fields[k]] <- list(mat_from_readmat(x[[k]]))
      }
      return(s)
    }
    items <- lapply(x, mat_from_readmat)
    return(mat_cell(items, 1L, length(items)))
  }
  x
}

#' Read Octave's text format (the default of Octave's save)
#' @noRd
mat_read_octave_text <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  i <- 1L
  hdr <- function(key) {
    while (i <= length(lines) && !startsWith(lines[i], paste0("# ", key, ":"))) {
      i <<- i + 1L
    }
    v <- trimws(sub("^#[^:]*:", "", lines[i]))
    i <<- i + 1L
    v
  }
  while (i <= length(lines)) {
    if (!startsWith(lines[i], "# name:")) {
      i <- i + 1L
      next
    }
    name <- hdr("name")
    type <- hdr("type")
    if (type %in% c("scalar", "bool")) {
      out[[name]] <- matrix(as.numeric(lines[i]))
      i <- i + 1L
    } else if (type %in% c("matrix", "bool matrix")) {
      r <- as.integer(hdr("rows"))
      cc <- as.integer(hdr("columns"))
      vals <- scan(text = lines[i:(i + r - 1L)], quiet = TRUE)
      out[[name]] <- matrix(vals, r, cc, byrow = TRUE)
      i <- i + r
    } else if (type %in% c("string", "sq_string")) {
      n_el <- as.integer(hdr("elements"))
      txt <- character(0)
      for (k in seq_len(n_el)) {
        hdr("length")
        txt <- c(txt, lines[i])
        i <- i + 1L
      }
      out[[name]] <- paste(txt, collapse = "")
    } else {
      stop("Octave data type '", type, "' in ", basename(path),
           " is not supported.", call. = FALSE)
    }
  }
  out
}

#' Variables stored in a data file (.mat, Octave text, or plain numbers)
#' @noRd
mat_load_file <- function(path) {
  con <- file(path, "rb")
  head <- readBin(con, "raw", 128L)
  close(con)
  htxt <- rawToChar(head[head != as.raw(0)])
  if (grepl("^# (Created by Octave|name:)", htxt)) {
    return(mat_read_octave_text(path))
  }
  if (startsWith(htxt, "MATLAB")) {
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
      stop("Reading the MATLAB file ", basename(path), " needs the ",
           "R.matlab package: install.packages(\"R.matlab\").", call. = FALSE)
    }
    if (grepl("MATLAB 7.3", htxt)) {
      stop(basename(path), " is a MATLAB v7.3 (HDF5) file; save it with ",
           "-v7 or earlier.", call. = FALSE)
    }
    raw <- R.matlab::readMat(path, fixNames = FALSE)
    return(lapply(raw, mat_from_readmat))
  }
  if (grepl("^Octave-1-", htxt)) {
    stop(basename(path), " is in Octave's binary format; save it with ",
         "save -text or -v7.", call. = FALSE)
  }
  # plain numbers (ASCII)
  m <- as.matrix(utils::read.table(path, header = FALSE))
  storage.mode(m) <- "double"
  stats::setNames(list(unname(m)),
                  gsub("[^A-Za-z0-9_]", "_",
                       tools::file_path_sans_ext(basename(path))))
}

mat_numeric_table <- function(df) {
  m <- suppressWarnings(apply(as.matrix(df), 2L, as.numeric))
  m <- matrix(m, nrow(df))
  keep_r <- which(rowSums(!is.na(m)) > 0L)
  keep_c <- which(colSums(!is.na(m)) > 0L)
  if (!length(keep_r) || !length(keep_c)) return(matrix(numeric(0), 0L, 0L))
  m <- m[min(keep_r):max(keep_r), min(keep_c):max(keep_c), drop = FALSE]
  m[is.na(m)] <- NaN
  m
}

#' Minimise f subject to bounds, linear and nonlinear constraints
#' (fmincon), by an augmented Lagrangian around L-BFGS-B
#' @noRd
mat_fmincon <- function(f, x0, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                        lb = NULL, ub = NULL, nonlcon = NULL) {
  n <- length(x0)
  lb <- if (length(lb)) rep_len(lb, n) else rep(-Inf, n)
  ub <- if (length(ub)) rep_len(ub, n) else rep(Inf, n)
  lb[!is.finite(lb)] <- -Inf
  ub[!is.finite(ub)] <- Inf
  cons <- function(x) {
    g <- numeric(0)
    h <- numeric(0)
    if (length(A)) g <- c(g, as.numeric(A %*% x - b))
    if (length(Aeq)) h <- c(h, as.numeric(Aeq %*% x - beq))
    if (!is.null(nonlcon)) {
      r <- nonlcon(x)
      g <- c(g, r$c)
      h <- c(h, r$ceq)
    }
    list(g = g, h = h)
  }
  c0 <- cons(pmin(pmax(x0, lb), ub))
  lam_g <- numeric(length(c0$g))
  lam_h <- numeric(length(c0$h))
  rho <- 10
  x <- pmin(pmax(x0, lb), ub)
  for (outer in seq_len(60L)) {
    x_prev <- x
    L <- function(z) {
      v <- f(z)
      if (!is.finite(v)) return(1e20)
      cc <- cons(z)
      if (length(cc$h)) v <- v + sum(lam_h * cc$h) + rho / 2 * sum(cc$h^2)
      if (length(cc$g)) {
        v <- v + sum(pmax(0, lam_g + rho * cc$g)^2 - lam_g^2) / (2 * rho)
      }
      v
    }
    opt <- stats::optim(x, L, method = "L-BFGS-B", lower = lb, upper = ub,
                        control = list(maxit = 2000, factr = 10))
    x <- opt$par
    cc <- cons(x)
    viol <- max(c(0, abs(cc$h), pmax(cc$g, 0)))
    if (length(cc$h)) lam_h <- lam_h + rho * cc$h
    if (length(cc$g)) lam_g <- pmax(0, lam_g + rho * cc$g)
    if (viol < 1e-9 && outer > 1L && max(abs(x - x_prev)) < 1e-10) break
    rho <- min(rho * 4, 1e8)
  }
  cc <- cons(x)
  viol <- max(c(0, abs(cc$h), pmax(cc$g, 0)))
  list(x = x, fval = f(x), exitflag = if (viol < 1e-6) 1 else -2)
}

mat_fun_of <- function(h, shape) {
  function(x) mat_scalar(mat_first(h$fn(list(array(x, shape)), 1L), 1L))
}

mat_opt_arg <- function(a, k) if (length(a) >= k && mat_numel(a[[k]]) > 0L) a[[k]] else NULL

mat_builtins$load <- function(a, n, ctx) {
  if (!length(a)) stop("load needs a file name", call. = FALSE)
  args <- vapply(a, mat_str, "")
  args <- args[!startsWith(args, "-")]
  path <- mat_find_file(args[1], ctx, c(".mat", ".txt", ".dat", ".csv"))
  vars <- mat_load_file(path)
  if (length(args) > 1L) vars <- vars[intersect(args[-1], names(vars))]
  if (n >= 1L) return(do.call(mat_struct, vars))
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = ctx$vars)
  mat_multi()
}
mat_builtins$xlsread <- function(a, n, ctx) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("xlsread needs the readxl package: install.packages(\"readxl\").",
         call. = FALSE)
  }
  path <- mat_find_file(mat_str(a[[1]]), ctx, c(".xlsx", ".xls"))
  sheet <- if (length(a) > 1L && mat_numel(a[[2]])) {
    if (mat_is_char(a[[2]])) mat_str(a[[2]]) else mat_scalar(a[[2]])
  } else 1L
  range <- if (length(a) > 2L && mat_is_char(a[[3]]) && nzchar(a[[3]])) mat_str(a[[3]]) else NULL
  df <- suppressMessages(readxl::read_excel(path, sheet = sheet, range = range,
                                            col_names = FALSE))
  num <- mat_numeric_table(df)
  if (n <= 1L) return(num)
  txt <- as.matrix(df)
  txt[!is.na(suppressWarnings(as.numeric(txt)))] <- ""
  txt[is.na(txt)] <- ""
  mat_ret(num, structure(array(as.list(txt), dim(txt)), class = "mat_cell"))
}
mat_builtins$csvread <- function(a, n, ctx) {
  path <- mat_find_file(mat_str(a[[1]]), ctx)
  r0 <- if (length(a) > 1L) mat_scalar(a[[2]]) else 0
  c0 <- if (length(a) > 2L) mat_scalar(a[[3]]) else 0
  m <- as.matrix(utils::read.csv(path, header = FALSE, skip = r0))
  storage.mode(m) <- "double"
  unname(m[, (c0 + 1):ncol(m), drop = FALSE])
}
mat_builtins$dlmread <- function(a, n, ctx) {
  path <- mat_find_file(mat_str(a[[1]]), ctx)
  sep <- if (length(a) > 1L && mat_is_char(a[[2]])) mat_str(a[[2]]) else ""
  sep <- gsub("\\\\t", "\t", sep)
  r0 <- if (length(a) > 2L) mat_scalar(a[[3]]) else 0
  c0 <- if (length(a) > 3L) mat_scalar(a[[4]]) else 0
  m <- as.matrix(utils::read.table(path, header = FALSE, sep = sep, skip = r0))
  storage.mode(m) <- "double"
  unname(m[, (c0 + 1):ncol(m), drop = FALSE])
}
mat_builtins$readmatrix <- function(a, n, ctx) {
  path <- mat_find_file(mat_str(a[[1]]), ctx)
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    return(mat_builtins$xlsread(list(path), 1L, ctx))
  }
  mat_numeric_table(utils::read.table(path, header = FALSE, fill = TRUE,
                                      sep = if (grepl("\\.csv$", path)) "," else ""))
}
mat_builtins$fmincon <- function(a, n, ctx) {
  h <- mat_handle_arg(a[[1]], ctx)
  x0 <- mat_m(a[[2]])
  get <- function(k) { v <- mat_opt_arg(a, k); if (is.null(v)) NULL else mat_m(v) }
  nl <- mat_opt_arg(a, 9)
  nonlcon <- if (is.null(nl)) NULL else {
    nh <- mat_handle_arg(nl, ctx)
    function(x) {
      r <- nh$fn(list(array(x, dim(x0))), 2L)
      if (!inherits(r, "mat_multi")) r <- mat_multi(r)
      list(c = as.numeric(mat_m(r[[1]])),
           ceq = if (length(r) > 1L) as.numeric(mat_m(r[[2]])) else numeric(0))
    }
  }
  r <- mat_fmincon(mat_fun_of(h, dim(x0)), as.vector(x0), get(3), as.vector(get(4)),
                   get(5), as.vector(get(6)), as.vector(get(7)),
                   as.vector(get(8)), nonlcon)
  mat_ret(array(r$x, dim(x0)), matrix(r$fval), matrix(r$exitflag))
}
mat_builtins$fminunc <- function(a, n, ctx) {
  h <- mat_handle_arg(a[[1]], ctx)
  x0 <- mat_m(a[[2]])
  f <- mat_fun_of(h, dim(x0))
  r <- stats::optim(as.vector(x0), f, method = "BFGS",
                    control = list(maxit = 5000, reltol = 1e-14))
  mat_ret(array(r$par, dim(x0)), matrix(r$value),
          matrix(if (r$convergence == 0) 1 else 0))
}
mat_builtins$lsqnonlin <- function(a, n, ctx) {
  h <- mat_handle_arg(a[[1]], ctx)
  x0 <- mat_m(a[[2]])
  lb <- if (length(a) > 2L && mat_numel(a[[3]])) as.vector(mat_m(a[[3]])) else -Inf
  ub <- if (length(a) > 3L && mat_numel(a[[4]])) as.vector(mat_m(a[[4]])) else Inf
  res <- function(x) mat_call_handle(h, x, shape = dim(x0))
  if (all(!is.finite(c(lb, ub)))) {
    r <- mat_newton(res, as.vector(x0))
    x <- r$x
  } else {
    x <- stats::optim(as.vector(x0), function(z) sum(res(z)^2), method = "L-BFGS-B",
                      lower = lb, upper = ub, control = list(factr = 10))$par
  }
  rr <- res(x)
  mat_ret(array(x, dim(x0)), matrix(sum(rr^2)), matrix(rr, ncol = 1L), matrix(1))
}
mat_hp <- function(y, lambda) {
  y <- mat_m(y)
  col <- nrow(y) == 1L
  if (col) y <- t(y)
  T <- nrow(y)
  D <- diff(diag(T), differences = 2L)
  trend <- solve(diag(T) + lambda * crossprod(D), y)
  if (col) list(trend = t(trend), cycle = t(y - trend)) else
    list(trend = trend, cycle = y - trend)
}
mat_builtins$hpfilter <- function(a, n, ctx) {
  lambda <- if (length(a) > 1L) mat_scalar(a[[2]]) else 1600
  r <- mat_hp(a[[1]], lambda)
  mat_ret(r$trend, r$cycle)
}
mat_builtins$sample_hp_filter <- function(a, n, ctx) {
  r <- mat_hp(a[[1]], mat_scalar(a[[2]]))
  mat_ret(r$trend, r$cycle)
}
mat_builtins$ksdensity <- function(a, n, ctx) {
  x <- as.vector(mat_m(a[[1]]))
  if (length(a) > 1L && is.numeric(a[[2]])) {
    pts <- as.vector(mat_m(a[[2]]))
    d <- stats::density(x, n = 512L)
    f <- stats::approx(d$x, d$y, xout = pts, rule = 2)$y
  } else {
    d <- stats::density(x, n = 100L)
    pts <- d$x
    f <- d$y
  }
  mat_ret(matrix(f, nrow = 1L), matrix(pts, nrow = 1L))
}
mat_builtins$interp1 <- function(a, n, ctx) {
  x <- as.vector(mat_m(a[[1]]))
  y <- as.vector(mat_m(a[[2]]))
  xi <- mat_m(a[[3]])
  method <- if (length(a) > 3L && mat_is_char(a[[4]])) tolower(mat_str(a[[4]])) else "linear"
  extrap <- any(vapply(a, function(v) mat_is_char(v) && tolower(v) == "extrap", TRUE))
  out <- if (method %in% c("spline", "pchip", "cubic")) {
    stats::spline(x, y, xout = as.vector(xi),
                  method = if (method == "spline") "fmm" else "hyman")$y
  } else {
    stats::approx(x, y, xout = as.vector(xi), rule = if (extrap) 2 else 1,
                  method = if (method %in% c("nearest", "previous")) "constant" else "linear")$y
  }
  out[is.na(out)] <- NaN
  array(out, dim(xi))
}
mat_builtins$polyfit <- function(a, n, ctx) {
  x <- as.vector(mat_m(a[[1]]))
  y <- as.vector(mat_m(a[[2]]))
  k <- mat_scalar(a[[3]])
  X <- outer(x, k:0, `^`)
  matrix(qr.solve(X, y), nrow = 1L)
}
mat_builtins$polyval <- function(a, n, ctx) {
  p <- as.vector(mat_m(a[[1]]))
  x <- mat_m(a[[2]])
  out <- 0 * x
  for (cf in p) out <- out * x + cf
  out
}
mat_builtins$prctile <- function(a, n, ctx) {
  x <- as.vector(mat_m(a[[1]]))
  p <- as.vector(mat_m(a[[2]]))
  matrix(stats::quantile(x, p / 100, names = FALSE, type = 5), nrow = 1L)
}
mat_builtins$quantile <- function(a, n, ctx) {
  x <- as.vector(mat_m(a[[1]]))
  p <- as.vector(mat_m(a[[2]]))
  matrix(stats::quantile(x, p, names = FALSE, type = 5), nrow = 1L)
}
mat_builtins$cov <- function(a, n, ctx) {
  x <- mat_m(a[[1]])
  if (length(a) > 1L) x <- cbind(as.vector(x), as.vector(mat_m(a[[2]])))
  if (nrow(x) == 1L) x <- t(x)
  v <- stats::cov(x)
  if (length(v) == 1L) matrix(v) else v
}
mat_builtins$corrcoef <- function(a, n, ctx) {
  x <- mat_m(a[[1]])
  if (length(a) > 1L) x <- cbind(as.vector(x), as.vector(mat_m(a[[2]])))
  stats::cor(x)
}
mat_builtins$corr <- mat_builtins$corrcoef
mat_builtins$trapz <- function(a, n, ctx) {
  if (length(a) == 1L) {
    y <- as.vector(mat_m(a[[1]]))
    return(matrix(sum((y[-1] + y[-length(y)]) / 2)))
  }
  x <- as.vector(mat_m(a[[1]]))
  y <- as.vector(mat_m(a[[2]]))
  matrix(sum(diff(x) * (y[-1] + y[-length(y)]) / 2))
}
mat_ignored <- c(mat_ignored, "save", "datatomfile", "dynasave", "dynatype",
                 "rplot", "print_info")
