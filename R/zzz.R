# Register the autoplot() methods with ggplot2's generic when ggplot2 is
# available (it is only suggested), now or whenever it is loaded later.

.onLoad <- function(libname, pkgname) {
  classes <- c("dsge_irf", "dsge_irf_2nd", "dsge_forecast",
               "dsge_variance_decomposition", "dsge_decomposition",
               "dsge_smoothed", "dsge_perfect_foresight", "dsge_occbin",
               "dsge_bayes")
  for (cl in classes) {
    .dsge_s3_register("ggplot2::autoplot", cl,
                      get(paste0("autoplot.", cl), mode = "function"))
  }
  invisible()
}

# Register an S3 method for a generic in a suggested package: immediately
# if the package is loaded, otherwise when it gets loaded.
.dsge_s3_register <- function(generic, class, method) {
  # evaluate now: the hook may run later, after the caller's loop moved on
  force(class)
  force(method)
  pieces <- strsplit(generic, "::", fixed = TRUE)[[1]]
  pkg <- pieces[[1]]
  gen <- pieces[[2]]
  register <- function(...) {
    envir <- asNamespace(pkg)
    if (exists(gen, envir = envir)) {
      registerS3method(gen, class, method, envir = envir)
    }
  }
  setHook(packageEvent(pkg, "onLoad"), function(...) register())
  if (isNamespaceLoaded(pkg)) register()
  invisible()
}
