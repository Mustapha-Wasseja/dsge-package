# Locate the companion repository of Paccagnini, *Macroeconometrics*
# (De Gruyter), https://github.com/AlessiaPaccagnini/Macroeconometrics_book
#
# The book's files are not copied into this repository. Set MACRO_BOOK_DIR
# to an existing clone, or the repository is cloned into a temporary
# directory on first use.

book_dir <- function() {
  d <- Sys.getenv("MACRO_BOOK_DIR", "")
  if (nzchar(d)) {
    if (!dir.exists(file.path(d, "Chapter6")))
      stop("MACRO_BOOK_DIR does not contain Chapter6/: ", d, call. = FALSE)
    return(normalizePath(d))
  }
  d <- file.path(tempdir(), "Macroeconometrics_book")
  if (!dir.exists(d)) {
    status <- system2("git", c("clone", "--depth", "1", "--quiet",
      "https://github.com/AlessiaPaccagnini/Macroeconometrics_book.git", d))
    if (status != 0L) stop("Could not clone the book repository.", call. = FALSE)
  }
  normalizePath(d)
}

# Evaluate only the function definitions of an R script (so that the book's
# own code can be called for comparison without running its top-level code).
source_functions <- function(file, envir = new.env(parent = globalenv())) {
  for (e in parse(file)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      eval(e, envir)
    }
  }
  envir
}

# Directory for figures and tables, next to the scripts.
out_dir <- function() {
  d <- file.path(script_dir(), "output")
  dir.create(d, showWarnings = FALSE)
  d
}

script_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
