## Submission summary

This is an update of the dsge package from 1.2.0 (the version on CRAN)
to 1.2.1. It improves speed, plotting and accuracy and fixes bugs; no
exported function has been removed. See NEWS.md for details.
Highlights:

- Speed: the Kalman filter, the first-order solver (cyclic reduction),
  and the tensor algebra of the second- and third-order solutions now
  run in C++ via Rcpp and RcppArmadillo. This is the first version with
  compiled code: Rcpp is a new import, and Rcpp and RcppArmadillo are
  new LinkingTo dependencies. The steady-state solver and the model
  linearisation use exact symbolic derivatives.
- Accuracy: `irf()` now iterates the state vector instead of forming
  matrix powers, which lost accuracy for models with linearly dependent
  states; imported Dynare models now match Dynare's impulse responses in
  all 53 models of the public DSGE_mod collection that both can solve.
- Bug fixes: `irf_2nd_order()` scaled shocks twice; `solve_dsge()` failed
  on models without stochastic shocks; `steady_state()` did not accept
  models imported with `read_dynare()`.
- Plots: a new, colour-vision-deficiency-safe style for all plot methods,
  and optional ggplot2 versions through `autoplot()` (ggplot2 is in
  Suggests; the methods are registered only when it is installed).
- Documentation: examples for all exported functions, and an updated
  vignette.
- `coda`, which was not used, has been removed from Suggests.

## Test environments

<!-- update before submission -->
- local: Ubuntu 24.04, R 4.3.3
- GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)
- win-builder: R-devel, R-release, R-oldrelease (to be run before
  submission)

## R CMD check results

<!-- update with the win-builder results before submission -->
`R CMD check --as-cran` locally:

0 errors | 0 warnings | 3 notes

- "checking installed package size ... NOTE: installed size is 7.1Mb;
  sub-directories of 1Mb or more: R 1.2Mb, libs 5.0Mb".
  The compiled code is C++ built with RcppArmadillo, whose templates
  generate a lot of debugging information under the default `-g`
  compiler flag. Of the 5.0 MB shared library, about 93% is debugging
  information: with `strip --strip-debug` the same library is 0.34 MB.
  The package does not ship any large data or prebuilt files, and the
  size of the installed library depends on each platform's default
  compiler flags.
- "checking for future file timestamps ... NOTE: unable to verify
  current time". This comes from the check machine being unable to reach
  an external time server, not from the package.
- "checking compilation flags used ... NOTE: Compilation used the
  following non-portable flag(s): '-mno-omit-leaf-frame-pointer'". This
  flag comes from the local R installation's own configuration (Ubuntu's
  R build), not from the package: its Makevars only define
  `ARMA_WARN_LEVEL` and link LAPACK and BLAS.

## Reverse dependencies

<!-- re-check before submission -->
There are currently no reverse dependencies: no CRAN package depends on,
imports, links to or suggests dsge
(`tools::package_dependencies("dsge", reverse = TRUE, which = "all")`,
checked 2026-09-24).
