# =============================================================================
# Extended NK Model — 3 Observables from FRED
# =============================================================================
#
# This example extends the standard 2-observable NK model by adding the
# output gap as a third observable and a cost-push shock for identification.
#
# DATA SOURCE:
#   Federal Reserve Economic Data (FRED), St. Louis Fed
#   - GDPDEF: GDP Implicit Price Deflator (quarterly, index)
#   - FEDFUNDS: Effective Federal Funds Rate (monthly -> quarterly avg)
#   - GDPC1: Real GDP (quarterly, billions of chained 2017 dollars)
#
# SAMPLE PERIOD:
#   1955Q1 to 2015Q4 (244 observations)
#
# VARIABLE CONSTRUCTION:
#   p = 400 * diff(log(GDPDEF))  — annualized quarterly inflation
#   r = quarterly average of FEDFUNDS — federal funds rate
#   x = 100 * HP-filter cycle of log(GDPC1) — output gap (% deviation)
#       HP filter with lambda = 1600 (standard for quarterly data)
#   All series are demeaned internally by bayes_dsge().
#
# MODEL (Extended New Keynesian, 3 shocks):
#   Phillips curve:  p  = beta * E[p'] + kappa * x + cp
#   IS curve:        x  = E[x'] - (r - E[p'] - g)
#   Taylor rule:     r  = psi_p * p + psi_x * x + u
#   Monetary shock:  u'  = rhou * u + e_u
#   Demand shock:    g'  = rhog * g + e_g
#   Cost-push shock: cp' = rhocp * cp + e_cp
#
# The cost-push shock provides the third structural disturbance needed
# for identification with 3 observables. It captures supply-side pressures
# (energy prices, supply chain disruptions) that shift the Phillips curve.
#
# The Taylor rule includes both inflation and output gap responses,
# consistent with a dual-mandate central bank.
#
# REQUIREMENTS:
#   - Internet connection (to download from FRED)
#   - dsge package (>= 0.3.0)
#
# TO RUN:
#   library(dsge)
#   source(system.file("examples", "bayesian_nk_extended.R", package = "dsge"))
#
# EXPECTED RUNTIME: ~10-15 minutes depending on hardware
# =============================================================================

library(dsge)

cat("============================================================\n")
cat(" Extended NK Model — 3 Observables from FRED\n")
cat("============================================================\n\n")

# =============================================================================
# Step 1: Download and prepare FRED data
# =============================================================================
cat("--- Step 1: Download FRED data (3 series) ---\n\n")

# GDP Deflator
gdpdef_url <- paste0(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?",
  "id=GDPDEF&cosd=1947-01-01&coed=2015-12-31&fq=Quarterly&fam=avg"
)
f1 <- tempfile(fileext = ".csv")
download.file(gdpdef_url, f1, quiet = TRUE)
gdpdef <- read.csv(f1, stringsAsFactors = FALSE)
colnames(gdpdef) <- c("date", "gdpdef")
gdpdef$date <- as.Date(gdpdef$date)
gdpdef$gdpdef <- as.numeric(gdpdef$gdpdef)
gdpdef <- gdpdef[!is.na(gdpdef$gdpdef), ]

# Federal Funds Rate
ffr_url <- paste0(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?",
  "id=FEDFUNDS&cosd=1954-07-01&coed=2015-12-31&fq=Quarterly&fam=avg"
)
f2 <- tempfile(fileext = ".csv")
download.file(ffr_url, f2, quiet = TRUE)
ffr <- read.csv(f2, stringsAsFactors = FALSE)
colnames(ffr) <- c("date", "ffr")
ffr$date <- as.Date(ffr$date)
ffr$ffr <- as.numeric(ffr$ffr)

# Real GDP
rgdp_url <- paste0(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?",
  "id=GDPC1&cosd=1947-01-01&coed=2015-12-31&fq=Quarterly&fam=avg"
)
f3 <- tempfile(fileext = ".csv")
download.file(rgdp_url, f3, quiet = TRUE)
rgdp <- read.csv(f3, stringsAsFactors = FALSE)
colnames(rgdp) <- c("date", "rgdp")
rgdp$date <- as.Date(rgdp$date)
rgdp$rgdp <- as.numeric(rgdp$rgdp)
rgdp <- rgdp[!is.na(rgdp$rgdp), ]

# Compute inflation
inflation <- data.frame(
  date = gdpdef$date[-1],
  p = 400 * diff(log(gdpdef$gdpdef))
)

# Compute output gap via HP filter (lambda = 1600)
hp_filter <- function(y, lambda = 1600) {
  n <- length(y)
  D <- matrix(0, n - 2, n)
  for (i in 1:(n - 2)) {
    D[i, i] <- 1
    D[i, i + 1] <- -2
    D[i, i + 2] <- 1
  }
  tau <- solve(diag(n) + lambda * crossprod(D), y)
  list(trend = tau, cycle = y - tau)
}

hp <- hp_filter(log(rgdp$rgdp))
output_gap <- data.frame(
  date = rgdp$date,
  x = 100 * hp$cycle  # percentage deviation from trend
)

# Merge all three and trim to sample period
macro <- merge(inflation, ffr, by = "date")
macro <- merge(macro, output_gap, by = "date")
macro <- macro[macro$date >= as.Date("1955-01-01") &
               macro$date <= as.Date("2015-12-31"), ]
colnames(macro)[colnames(macro) == "ffr"] <- "r"

cat(sprintf("  Sample: %d observations (%s to %s)\n",
            nrow(macro), min(macro$date), max(macro$date)))
cat(sprintf("  Inflation (p):   mean = %6.2f, sd = %5.2f\n",
            mean(macro$p), sd(macro$p)))
cat(sprintf("  Fed Funds (r):   mean = %6.2f, sd = %5.2f\n",
            mean(macro$r), sd(macro$r)))
cat(sprintf("  Output gap (x):  mean = %6.2f, sd = %5.2f\n",
            mean(macro$x), sd(macro$x)))

dat <- data.frame(p = macro$p, r = macro$r, x = macro$x)

# =============================================================================
# Step 2: Define the extended NK model
# =============================================================================
cat("\n--- Step 2: Define extended NK model ---\n\n")

nk3 <- dsge_model(
  obs(p   ~ beta * lead(p) + kappa * x + cp),   # Phillips + cost-push
  obs(x   ~ lead(x) - (r - lead(p) - g)),       # IS curve
  obs(r   ~ psi_p * p + psi_x * x + u),         # Taylor rule (dual mandate)
  state(u  ~ rhou * u),                          # Monetary shock
  state(g  ~ rhog * g),                          # Demand shock
  state(cp ~ rhocp * cp),                        # Cost-push shock
  start = list(beta = 0.95, kappa = 0.15, psi_p = 1.5, psi_x = 0.5,
               rhou = 0.5, rhog = 0.7, rhocp = 0.5)
)
print(nk3)

# =============================================================================
# Step 3: Specify priors
# =============================================================================
cat("\n--- Step 3: Prior specification ---\n\n")

# Informative priors:
#   beta:    near 1 (quarterly discount factor)
#   kappa:   small positive (price stickiness)
#   psi_p:   > 1 for Taylor principle
#   psi_x:   positive (output gap stabilization)
#   rho's:   centered at 0.5 with moderate uncertainty

priors3 <- list(
  beta   = prior("beta", shape1 = 95, shape2 = 5),
  kappa  = prior("beta", shape1 = 30, shape2 = 70),
  psi_p  = prior("gamma", shape = 184, rate = 122.7),
  psi_x  = prior("gamma", shape = 4, rate = 8),
  rhou   = prior("beta", shape1 = 10, shape2 = 10),
  rhog   = prior("beta", shape1 = 10, shape2 = 10),
  rhocp  = prior("beta", shape1 = 10, shape2 = 10)
)

cat("  Parameter  Prior                    Mean   SD\n")
cat("  ---------  -----------------------  -----  -----\n")
cat("  beta       beta(95, 5)              0.950  0.022\n")
cat("  kappa      beta(30, 70)             0.300  0.046\n")
cat("  psi_p      gamma(184, 122.7)        1.500  0.111\n")
cat("  psi_x      gamma(4, 8)              0.500  0.250\n")
cat("  rhou       beta(10, 10)             0.500  0.109\n")
cat("  rhog       beta(10, 10)             0.500  0.109\n")
cat("  rhocp      beta(10, 10)             0.500  0.109\n")
cat("  shock SDs  inv_gamma(0.01, 0.01)    ---    ---\n")

# =============================================================================
# Step 4: Bayesian estimation
# =============================================================================
cat("\n--- Step 4: MCMC estimation ---\n\n")
cat("  Chains: 2, Iterations: 15000, Warmup: 7500\n\n")

t0 <- proc.time()
fit3 <- bayes_dsge(nk3, data = dat,
                   priors = priors3,
                   chains = 2,
                   iter = 15000,
                   warmup = 7500,
                   seed = 123)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Elapsed: %.0f seconds (%.1f minutes)\n\n", elapsed, elapsed / 60))

# =============================================================================
# Step 5: Posterior summary
# =============================================================================
cat("--- Step 5: Posterior summary ---\n\n")
summary(fit3)

# =============================================================================
# Step 6: Diagnostics and interpretation
# =============================================================================
cat("\n--- Step 6: Diagnostics and interpretation ---\n\n")

est <- coef(fit3)

cat("Convergence diagnostics:\n")
cat(sprintf("  Acceptance rates:  %s\n",
            paste(round(fit3$acceptance_rates, 3), collapse = ", ")))
cat(sprintf("  ESS range:         %.0f - %.0f\n",
            min(fit3$diagnostics$ess), max(fit3$diagnostics$ess)))
cat(sprintf("  Max R-hat:         %.4f\n",
            max(fit3$diagnostics$rhat, na.rm = TRUE)))

cat("\nEconomic interpretation:\n")
cat(sprintf("  Discount factor (beta):          %.4f\n", est["beta"]))
cat(sprintf("  Phillips curve slope (kappa):    %.4f\n", est["kappa"]))
cat(sprintf("  Taylor rule on inflation (psi_p): %.4f  %s\n",
            est["psi_p"],
            ifelse(est["psi_p"] > 1, "[Taylor principle satisfied]", "")))
cat(sprintf("  Taylor rule on output (psi_x):   %.4f\n", est["psi_x"]))
cat(sprintf("  Monetary shock persistence:      %.4f\n", est["rhou"]))
cat(sprintf("  Demand shock persistence:        %.4f\n", est["rhog"]))
cat(sprintf("  Cost-push shock persistence:     %.4f\n", est["rhocp"]))

# =============================================================================
# Step 7: Impulse-response functions
# =============================================================================
cat("\n--- Step 7: Posterior IRFs ---\n\n")

irfs <- irf(fit3, periods = 12, n_draws = 300)
print(irfs)

cat("\nMonetary shock (u) at impact — contractionary:\n")
cat("  A contractionary monetary shock raises rates, reducing output and inflation.\n")
for (v in c("p", "x", "r")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "u" &
                   irfs$data$response == v, "value"]
  if (length(val) > 0) cat(sprintf("  %s = %+7.4f\n", v, val))
}

cat("\nDemand shock (g) at impact — expansionary:\n")
cat("  A positive demand shock boosts output and inflation; the central bank responds.\n")
for (v in c("p", "x", "r")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "g" &
                   irfs$data$response == v, "value"]
  if (length(val) > 0) cat(sprintf("  %s = %+7.4f\n", v, val))
}

cat("\nCost-push shock (cp) at impact — stagflationary:\n")
cat("  A supply-side shock raises inflation while depressing output.\n")
for (v in c("p", "x", "r")) {
  val <- irfs$data[irfs$data$period == 0 & irfs$data$impulse == "cp" &
                   irfs$data$response == v, "value"]
  if (length(val) > 0) cat(sprintf("  %s = %+7.4f\n", v, val))
}

cat("\n============================================================\n")
cat(" Extended NK estimation complete.\n")
cat(" 3 observables (inflation, FFR, output gap) from FRED.\n")
cat("============================================================\n")
