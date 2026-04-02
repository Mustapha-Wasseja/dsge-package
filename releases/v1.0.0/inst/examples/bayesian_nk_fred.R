# =============================================================================
# Bayesian New Keynesian Model — US Macroeconomic Data
# =============================================================================
#
# This example estimates a standard 3-equation New Keynesian model using
# Bayesian methods on real US macroeconomic data from FRED.
#
# DATA SOURCE:
#   Federal Reserve Economic Data (FRED), St. Louis Fed
#   - GDPDEF: GDP Implicit Price Deflator (quarterly, index)
#   - FEDFUNDS: Effective Federal Funds Rate (monthly -> quarterly avg)
#
# SAMPLE PERIOD:
#   1955Q1 to 2015Q4 (244 observations)
#
# VARIABLE CONSTRUCTION:
#   p = 400 * diff(log(GDPDEF))  — annualized quarterly inflation
#   r = quarterly average of FEDFUNDS — federal funds rate
#   Both series are demeaned internally by bayes_dsge().
#
# MODEL (3-equation New Keynesian):
#   Phillips curve:  p  = beta * E[p'] + kappa * x
#   IS curve:        x  = E[x'] - (r - E[p'] - g)
#   Taylor rule:     r  = psi * p + u
#   Monetary shock:  u' = rhou * u + e_u
#   Demand shock:    g' = rhog * g + e_g
#
# The model features forward-looking inflation expectations, an output gap
# (x) that is unobserved and filtered from the data, a Taylor rule linking
# the interest rate to inflation, and two AR(1) structural shocks.
#
# PRIOR SPECIFICATION:
#   Parameter  Distribution             Prior Mean  Prior SD
#   ---------  -----------------------  ----------  --------
#   beta       beta(95, 5)              0.950       0.022
#   kappa      beta(30, 70)             0.300       0.046
#   psi        gamma(184, 122.7)        1.500       0.111
#   rhou       beta(70, 20)             0.778       0.043
#   rhog       beta(70, 20)             0.778       0.043
#   sd(e.u)    inv_gamma(0.01, 0.01)    (diffuse)
#   sd(e.g)    inv_gamma(0.01, 0.01)    (diffuse)
#
# REQUIREMENTS:
#   - Internet connection (to download from FRED)
#   - dsge package (>= 0.3.0)
#
# TO RUN:
#   library(dsge)
#   source(system.file("examples", "bayesian_nk_fred.R", package = "dsge"))
#
# EXPECTED RUNTIME: ~10-15 minutes depending on hardware
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Bayesian NK Model — US Macroeconomic Data\n")
cat("============================================================\n\n")

# =============================================================================
# Step 1: Download and prepare FRED data
# =============================================================================
cat("--- Step 1: Download and prepare FRED data ---\n\n")

# GDP Implicit Price Deflator (quarterly, seasonally adjusted)
gdpdef_url <- paste0(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?",
  "id=GDPDEF&cosd=1947-01-01&coed=2015-12-31&fq=Quarterly&fam=avg"
)
gdpdef_file <- tempfile(fileext = ".csv")
download.file(gdpdef_url, gdpdef_file, quiet = TRUE)
gdpdef <- read.csv(gdpdef_file, stringsAsFactors = FALSE)
colnames(gdpdef) <- c("date", "gdpdef")
gdpdef$date <- as.Date(gdpdef$date)
gdpdef$gdpdef <- as.numeric(gdpdef$gdpdef)
gdpdef <- gdpdef[!is.na(gdpdef$gdpdef), ]

# Effective Federal Funds Rate (quarterly average)
ffr_url <- paste0(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?",
  "id=FEDFUNDS&cosd=1954-07-01&coed=2015-12-31&fq=Quarterly&fam=avg"
)
ffr_file <- tempfile(fileext = ".csv")
download.file(ffr_url, ffr_file, quiet = TRUE)
ffr <- read.csv(ffr_file, stringsAsFactors = FALSE)
colnames(ffr) <- c("date", "ffr")
ffr$date <- as.Date(ffr$date)
ffr$ffr <- as.numeric(ffr$ffr)

# Compute annualized quarterly inflation: p = 400 * diff(log(GDPDEF))
inflation <- data.frame(
  date = gdpdef$date[-1],
  p = 400 * diff(log(gdpdef$gdpdef))
)

# Merge and trim to sample period
macro <- merge(inflation, ffr, by = "date")
macro <- macro[macro$date >= as.Date("1955-01-01") &
               macro$date <= as.Date("2015-12-31"), ]
colnames(macro)[colnames(macro) == "ffr"] <- "r"

cat(sprintf("  Sample: %d observations (%s to %s)\n",
            nrow(macro), min(macro$date), max(macro$date)))
cat(sprintf("  Inflation (p): mean = %.2f%%, sd = %.2f%%\n",
            mean(macro$p), sd(macro$p)))
cat(sprintf("  Fed Funds (r): mean = %.2f%%, sd = %.2f%%\n",
            mean(macro$r), sd(macro$r)))

dat <- data.frame(p = macro$p, r = macro$r)

# =============================================================================
# Step 2: Define the New Keynesian model
# =============================================================================
cat("\n--- Step 2: Define NK model ---\n\n")

nk <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x),       # Phillips curve
  unobs(x ~ lead(x) - (r - lead(p) - g)),      # IS curve
  obs(r   ~ psi * p + u),                       # Taylor rule
  state(u ~ rhou * u),                          # Monetary shock
  state(g ~ rhog * g),                          # Demand shock
  start = list(beta = 0.95, kappa = 0.15, psi = 1.5, rhou = 0.7, rhog = 0.9)
)
print(nk)

# =============================================================================
# Step 3: Specify priors
# =============================================================================
cat("\n--- Step 3: Prior specification ---\n\n")

# Informative priors reflecting standard macroeconomic beliefs:
#   beta:  near 1 (agents are patient; quarterly discount factor ~0.95)
#   kappa: small and positive (moderate price stickiness)
#   psi:   greater than 1 (Taylor principle: central bank raises rates
#          more than one-for-one with inflation)
#   rhou, rhog: persistent but stationary shock processes

my_priors <- list(
  beta  = prior("beta", shape1 = 95, shape2 = 5),
  kappa = prior("beta", shape1 = 30, shape2 = 70),
  psi   = prior("gamma", shape = 184, rate = 122.7),
  rhou  = prior("beta", shape1 = 70, shape2 = 20),
  rhog  = prior("beta", shape1 = 70, shape2 = 20)
  # Shock SDs default to inv_gamma(0.01, 0.01)
)

cat("  Parameter  Prior                    Mean   SD\n")
cat("  ---------  -----------------------  -----  -----\n")
cat("  beta       beta(95, 5)              0.950  0.022\n")
cat("  kappa      beta(30, 70)             0.300  0.046\n")
cat("  psi        gamma(184, 122.7)        1.500  0.111\n")
cat("  rhou       beta(70, 20)             0.778  0.043\n")
cat("  rhog       beta(70, 20)             0.778  0.043\n")
cat("  sd(e.u)    inv_gamma(0.01, 0.01)    ---    ---\n")
cat("  sd(e.g)    inv_gamma(0.01, 0.01)    ---    ---\n")

# =============================================================================
# Step 4: Bayesian estimation
# =============================================================================
cat("\n--- Step 4: MCMC estimation ---\n\n")
cat("  Chains: 2, Iterations: 20000, Warmup: 10000\n\n")

t0 <- proc.time()
fit <- bayes_dsge(nk, data = dat,
                  priors = my_priors,
                  chains = 2,
                  iter = 20000,
                  warmup = 10000,
                  seed = 42)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Elapsed: %.0f seconds (%.1f minutes)\n\n", elapsed, elapsed / 60))

# =============================================================================
# Step 5: Posterior summary
# =============================================================================
cat("--- Step 5: Posterior summary ---\n\n")
summary(fit)

# =============================================================================
# Step 6: Diagnostics and interpretation
# =============================================================================
cat("\n--- Step 6: Diagnostics and interpretation ---\n\n")

est <- coef(fit)

cat("Convergence diagnostics:\n")
cat(sprintf("  Acceptance rates:  %s\n",
            paste(round(fit$acceptance_rates, 3), collapse = ", ")))
cat(sprintf("  ESS range:         %.0f - %.0f  (higher is better)\n",
            min(fit$diagnostics$ess), max(fit$diagnostics$ess)))
cat(sprintf("  Max R-hat:         %.4f  (< 1.1 indicates convergence)\n",
            max(fit$diagnostics$rhat, na.rm = TRUE)))

cat("\nEconomic interpretation:\n")
cat(sprintf("  Discount factor (beta):    %.4f  [quarterly; ~%.1f%% annual rate]\n",
            est["beta"], (1/est["beta"] - 1) * 400))
cat(sprintf("  Phillips slope (kappa):    %.4f  [price flexibility]\n",
            est["kappa"]))
cat(sprintf("  Taylor rule coeff (psi):   %.4f  %s\n",
            est["psi"],
            ifelse(est["psi"] > 1,
                   "[> 1: Taylor principle satisfied]",
                   "[< 1: Taylor principle violated]")))
cat(sprintf("  Monetary persistence:      %.4f\n", est["rhou"]))
cat(sprintf("  Demand persistence:        %.4f\n", est["rhog"]))

# =============================================================================
# Step 7: Impulse-response functions
# =============================================================================
cat("\n--- Step 7: Posterior IRFs ---\n\n")

irfs <- irf(fit, periods = 8, n_draws = 500)
print(irfs)

# Check that IRF signs are economically sensible
irf_data <- irfs$data

cat("\nMonetary policy shock (u) at impact — contractionary:\n")
cat("  A positive monetary shock raises the interest rate.\n")
cat("  Economic theory predicts: inflation falls, output gap falls.\n")
for (v in c("p", "x", "r")) {
  val <- irf_data[irf_data$period == 0 & irf_data$impulse == "u" &
                  irf_data$response == v, "value"]
  cat(sprintf("  %s = %+7.4f\n", v, val))
}

cat("\nDemand shock (g) at impact — expansionary:\n")
cat("  A positive demand shock boosts output and inflation.\n")
cat("  The central bank raises rates in response.\n")
for (v in c("p", "x", "r")) {
  val <- irf_data[irf_data$period == 0 & irf_data$impulse == "g" &
                  irf_data$response == v, "value"]
  cat(sprintf("  %s = %+7.4f\n", v, val))
}

cat("\n============================================================\n")
cat(" Bayesian estimation complete.\n")
cat("============================================================\n")
