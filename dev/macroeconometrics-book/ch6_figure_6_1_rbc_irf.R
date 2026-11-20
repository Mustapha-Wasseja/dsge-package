# Paccagnini, *Macroeconometrics*, Chapter 6 -- Figure 6.1
# RBC impulse responses to a one percent technology shock, with dsge.
#
# The book (Chapter6/Codes_Figures/chapter6_dsge_companion.R, Section 2)
# log-linearises the model by hand and solves it with a Blanchard-Kahn
# eigenvalue decomposition. Here the nonlinear model is written down once
# and dsge linearises and solves it (Klein / QZ). The book's own functions
# are then run on the same calibration and the two sets of responses
# compared.
#
# Run from the package root:
#   Rscript dev/macroeconometrics-book/ch6_figure_6_1_rbc_irf.R

suppressPackageStartupMessages(library(dsge))
local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
  source(file.path(here, "book_path.R"))
})

# ---------------------------------------------------------------------------
# 1. The model (book calibration: beta = 0.99, alpha = 0.33, delta = 0.025,
#    rho_A = 0.95, sigma = 1, eta = 1)
#
#    max E sum beta^t [ log C_t - N_t^(1+eta) / (1+eta) ]
#    Y = A K^alpha N^(1-alpha),  Y = C + I,  K' = (1-delta) K + I,
#    log A' = rho log A + e.
#
#    All variables are in logs, so responses are log deviations; with a
#    shock of 0.01 they are multiplied by 100 to give percent deviations.
#    k is capital at the start of the period (predetermined).
# ---------------------------------------------------------------------------

rbc <- dsgenl_model(
  "exp(-c) = beta * exp(-c(+1)) * (alpha * exp(y(+1) - k(+1)) + 1 - delta)",
  "exp(eta * n + c) = (1 - alpha) * exp(y - n)",
  "exp(y) = exp(a + alpha * k + (1 - alpha) * n)",
  "exp(y) = exp(c) + exp(i)",
  "k(+1) = log((1 - delta) * exp(k) + exp(i))",
  "a(+1) = rho * a",
  observed   = "y",
  unobserved = c("c", "i", "n"),
  endo_state = "k",
  exo_state  = "a",
  fixed = list(alpha = 0.33, beta = 0.99, delta = 0.025, eta = 1, rho = 0.95),
  ss_guess = c(y = 0, c = -0.3, i = -1.5, n = -1, k = 2, a = 0)
)

params <- c(alpha = 0.33, beta = 0.99, delta = 0.025, eta = 1, rho = 0.95)
ss  <- steady_state(rbc, params = params)
sol <- solve_dsge(rbc, params = params, shock_sd = c(a = 0.01))

cat("Steady state (levels):\n")
lv <- exp(ss$values[c("y", "c", "i", "k", "n")])
cat(sprintf("  Y = %.4f, C = %.4f, I = %.4f, K = %.4f, N = %.4f\n",
            lv["y"], lv["c"], lv["i"], lv["k"], lv["n"]))
cat(sprintf("  C/Y = %.3f, I/Y = %.3f, K/Y = %.2f\n\n",
            lv["c"] / lv["y"], lv["i"] / lv["y"], lv["k"] / lv["y"]))
print(sol)

# ---------------------------------------------------------------------------
# 2. Impulse responses, 40 quarters
#    The book plots end-of-period capital K_{t+1}; dsge's k is K_t, so the
#    book's "Capital" at quarter h is dsge's k at h + 1.
# ---------------------------------------------------------------------------

periods <- 40
ir <- irf(sol, periods = periods, se = FALSE)$data
get <- function(v) 100 * ir$value[ir$response == v][order(ir$period[ir$response == v])]
dsge_irf <- data.frame(
  quarter = 0:(periods - 1),
  Y = get("y")[1:periods], C = get("c")[1:periods], I = get("i")[1:periods],
  N = get("n")[1:periods], K = get("k")[2:(periods + 1)],
  A = get("a")[1:periods])

# ---------------------------------------------------------------------------
# 3. The book's hand-coded Blanchard-Kahn solution, same calibration
# ---------------------------------------------------------------------------

book_file <- file.path(book_dir(), "Chapter6", "Codes_Figures",
                       "chapter6_dsge_companion.R")
bk <- source_functions(book_file)
book_ss  <- bk$rbc_steady_state()
book_sol <- capture.output(s <- bk$rbc_solve_bk(book_ss))
book_irf <- bk$rbc_irfs(s, book_ss)

vars <- c("Y", "C", "I", "N", "K", "A")
cmp <- data.frame(
  variable    = vars,
  dsge_impact = sapply(vars, function(v) dsge_irf[[v]][1]),
  book_impact = sapply(vars, function(v) book_irf[[v]][1]),
  dsge_peak   = sapply(vars, function(v) max(dsge_irf[[v]])),
  peak_q      = sapply(vars, function(v) which.max(dsge_irf[[v]]) - 1),
  max_abs_diff_40q = sapply(vars, function(v) max(abs(dsge_irf[[v]] - book_irf[[v]]))),
  row.names = NULL)
cat("\nFigure 6.1: dsge vs the book's Blanchard-Kahn solution (percent)\n")
print(format(cmp, digits = 5), row.names = FALSE)
cat(sprintf("\nLargest difference over all variables and 40 quarters: %.2e\n",
            max(cmp$max_abs_diff_40q)))

utils::write.csv(dsge_irf, file.path(out_dir(), "figure_6_1_rbc_irf_dsge.csv"),
                 row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Figure 6.1, drawn from the dsge solution (the book's values as markers)
# ---------------------------------------------------------------------------

col_dsge <- "#2a78d6"; col_book <- "#eb6834"; col_grid <- "#e4e3df"
titles <- c(Y = "Output (Y)", C = "Consumption (C)", I = "Investment (I)",
            N = "Hours (N)", K = "Capital (K, end of period)",
            A = "Technology (A)")

png(file.path(out_dir(), "figure_6_1_rbc_irf_dsge.png"),
    width = 1500, height = 1000, res = 150)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2.5, 0),
    las = 1, bty = "l", col.axis = "#52514e", col.lab = "#52514e")
for (v in vars) {
  rng <- range(0, dsge_irf[[v]], book_irf[[v]])
  plot(dsge_irf$quarter, dsge_irf[[v]], type = "n", ylim = rng,
       xlab = "Quarters", ylab = "% deviation", main = titles[[v]],
       font.main = 1, cex.main = 1.05)
  abline(h = pretty(rng), col = col_grid, lwd = 0.6)
  abline(h = 0, col = "#8a8984", lwd = 0.8)
  points(book_irf$quarter[seq(1, periods, 3)], book_irf[[v]][seq(1, periods, 3)],
         pch = 1, cex = 0.9, col = col_book, lwd = 1.2)
  lines(dsge_irf$quarter, dsge_irf[[v]], col = col_dsge, lwd = 2)
  if (v == "Y")
    legend("topright", c("dsge (Klein/QZ)", "Book (Blanchard-Kahn)"),
           col = c(col_dsge, col_book), lwd = c(2, NA), pch = c(NA, 1),
           bty = "n", cex = 0.9)
}
mtext("Figure 6.1  Impulse responses to a one percent technology shock (dsge)",
      outer = TRUE, cex = 1.1, font = 2)
invisible(dev.off())
cat("\nSaved", file.path(out_dir(), "figure_6_1_rbc_irf_dsge.png"), "\n")
