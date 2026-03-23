# Release Risk Review — dsge v0.4.0

## Timing Risks

### Tests
- All Bayesian tests have `skip_on_cran()` — they will NOT run on CRAN.
- Non-Bayesian tests run in <10 seconds total. No timing risk.

### Examples
- `bayes_dsge()` example is wrapped in `\donttest{}` — will NOT run
  during `R CMD check --as-cran`.
- `estimate()` and `irf()` examples with `\donttest{}` are fast (<2 sec)
  but also guarded.
- All other examples (model specification, solving, steady state) run
  in <1 second each.

### Vignette
- Bayesian code chunks are `eval=FALSE` — no runtime cost.
- All evaluated chunks are ML estimation or model solving — fast (<5 sec total).
- Vignette requires pandoc + rmarkdown to build. On systems without pandoc,
  the vignette won't build but this does not cause an ERROR.

## Fragility Risks

### Bayesian Nonlinear Estimation
- **Steady-state solver convergence:** The Newton-Raphson solver may fail
  for extreme parameter draws. This is handled gracefully (returns -Inf,
  proposal rejected). In testing, <0.1% of draws fail.
- **Linearization singularity:** If the Jacobian is singular at a candidate
  steady state, the draw is rejected. This is rare in practice.
- **BK conditions:** Extreme parameters can violate Blanchard-Kahn saddle-path
  conditions. Again handled by returning -Inf.
- **Overall:** The sampler is robust — it never crashes. The failure rate
  was 9/10000 (0.09%) on the NK example and 0/10000 on the RBC example.

### Numerical Stability
- The Klein solver uses direct matrix inversion. For large models (>10 states),
  this could become numerically unstable. Current examples use 2-5 states.
- The Kalman filter uses standard recursions. No known stability issues
  for the model sizes tested.

### Platform Dependence
- No compiled code (pure R). Should work identically on all platforms.
- FRED data download examples require internet access — they are in
  `inst/examples/`, not in tests or vignettes.
- `numDeriv` is used for Jacobians and Hessians. This is a well-established
  CRAN package with no known issues.

## Submission Recommendation

**Recommended path: GitHub release first, then CRAN.**

Rationale:
1. The Bayesian nonlinear feature is new and has been tested on only 2
   nonlinear model specifications. A GitHub-only release period allows
   early adopters to find edge cases before CRAN exposure.
2. v0.3.0 (Bayesian linear) has not yet been submitted to CRAN either.
   Submitting v0.4.0 directly skips v0.3.0, which is fine but means
   CRAN reviewers see two major feature additions at once.
3. The vignette pandoc dependency is a known source of CRAN check
   warnings. A pre-built vignette (via `R CMD build`) resolves this,
   but requires pandoc to be installed locally.
4. All code-level checks pass cleanly (0 errors, 0 notes on code).

**Alternative:** If you want to go directly to CRAN, the package is
technically ready. The only manual step is building the vignette
(requires pandoc) and verifying on win-builder.

## Items That Could Trigger CRAN Reviewer Questions

1. **Package title** — "Dynamic Stochastic General Equilibrium Models" is
   descriptive and appropriate. No issues expected.
2. **Description field** — mentions Klein (2000) with DOI. Good.
3. **`\donttest{}` examples** — CRAN allows this for slow examples.
   The examples are correct and will pass if run manually.
4. **No reverse dependencies** — this is a new package. No breakage risk.
5. **`coda` in Suggests** — only used for optional `as.mcmc()` conversion.
   This is a standard pattern.
