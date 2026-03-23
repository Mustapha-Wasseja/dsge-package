# Release Checklist — dsge v0.4.0

## Pre-submission steps

### 1. Verify version consistency
- [ ] `DESCRIPTION` version: 0.4.0
- [ ] `NEWS.md` has v0.4.0 section at top
- [ ] `README.md` references v0.4.0
- [ ] `cran-comments.md` references v0.4.0

### 2. Regenerate documentation
```r
roxygen2::roxygenise()
```
Verify no new warnings. Check that `man/bayes_dsge.Rd` says
"Estimate a DSGE Model" (not "Linear DSGE Model").

### 3. Run local tests
```r
devtools::test()
```
Expected: 167 tests pass, 0 failures, 1 skip (empty test placeholder).

### 4. Build the package tarball
```bash
cd "C:\Users\musta\Dropbox\Econometrics Papers\DSGE in R"
"C:\Program Files\R\R-4.5.2\bin\R.exe" CMD build dsge
```
This produces `dsge_0.4.0.tar.gz`.

### 5. Run R CMD check --as-cran
```bash
"C:\Program Files\R\R-4.5.2\bin\R.exe" CMD check --as-cran dsge_0.4.0.tar.gz
```
Expected results:
- 0 errors
- 0 warnings (if pandoc is available for vignette build)
- 0 notes (or 1 note about "New submission" if first CRAN upload)

### 6. Check on win-builder (optional but recommended)
Upload `dsge_0.4.0.tar.gz` to:
- https://win-builder.r-project.org/upload.aspx (R-release and R-devel)

Wait for email results. Fix any issues found.

### 7. Check on R-hub (optional)
```r
rhub::check_for_cran()
```

## CRAN submission

### 8. Submit to CRAN
```r
devtools::release()
```
Or upload manually at: https://cran.r-project.org/submit.html

### 9. CRAN comments to include
The file `cran-comments.md` is ready. Key points:
- First submission (or update from previous version)
- 0 errors | 0 warnings | 0 notes
- No compiled code, minimal dependencies
- Bayesian nonlinear uses first-order approximation (not fully nonlinear)

## Post-submission

### 10. Tag the release
```bash
git tag -a v0.4.0 -m "v0.4.0: Bayesian nonlinear DSGE estimation"
git push origin v0.4.0
```

### 11. Create GitHub release (if applicable)
- Title: "dsge v0.4.0 — Bayesian Nonlinear DSGE"
- Attach `dsge_0.4.0.tar.gz`
- Include NEWS.md highlights

## Manual verification before submission

- [ ] Run `bayesian_rbc_nonlinear.R` example end-to-end
- [ ] Run `bayesian_nk_nonlinear.R` example end-to-end
- [ ] Run `bayesian_nk_fred.R` example end-to-end (requires internet)
- [ ] Verify all IRF signs are economically sensible
- [ ] Verify posterior recovers true parameters in simulation examples
- [ ] Check that the package installs cleanly from the tarball:
      `install.packages("dsge_0.4.0.tar.gz", repos = NULL)`
