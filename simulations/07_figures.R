# Step 7: report figures (vector SVG, Okabe-Ito colours, no in-image titles)

source("fmax_exact.R")
set.seed(2421974)
fig_dir <- "../04_report/fig"; dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
# FIG_PDF=1 -> vector PDF for the LaTeX report instead of SVG for the HTML one
if (nzchar(Sys.getenv("FIG_PDF")))
  svg <- function(file, width, height) cairo_pdf(sub("\\.svg$", ".pdf", file), width, height)
col_sim <- "#BFD7EA"; col_line <- "#D55E00"; col_tab <- "#0072B2"; col_grey <- "grey40"
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

dfmax <- function(c, nu) {                      # density by central difference
  h <- 1e-4 * c
  (pfmax(c + h, rep(nu, length(c))) - pfmax(c - h, rep(nu, length(c)))) / (2 * h)
}
panel <- function(s, nu, lab, xlab, sq = FALSE) {  # histogram of log(stat) + law
  z <- log(if (sq) s^2 else s)
  br <- seq(0, quantile(z, .995), length.out = 45)
  z2 <- z[z < max(br)]
  hist(z2, breaks = br, freq = FALSE, col = col_sim, border = "white",
       main = "", xlab = xlab, ylab = "Density", las = 1)
  u <- seq(0.002, max(br), length.out = 160)
  c <- exp(u)
  lines(u, dfmax(c, nu) * c * length(z) / length(z2), col = col_line, lwd = 2.2)
  ks <- suppressWarnings(ks.test(if (sq) s^2 else s, function(q) pfmax(q, rep(nu, length(q))))$p.value)
  legend("topright", bty = "n", cex = .85,
         legend = c(lab, sprintf("KS p = %.2f", ks)), text.col = c("black", col_grey))
}

## Figure 1: KaH-III simulated under H0 vs exact Hartley Fmax(3, nu)
kah3 <- function(n, p) {
  x <- scale(matrix(runif(n * p), n, p)); x <- x[order(x[, 1]), , drop = FALSE]
  y <- drop(1 + x %*% rep(1, p) + rnorm(n)); m <- floor(n / 3)
  idx <- list(1:m, (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
  ms <- sapply(idx, function(i) { f <- lm.fit(cbind(1, x[i, , drop = FALSE]), y[i]); sum(f$residuals^2) / (length(i) - p - 1) })
  max(ms) / min(ms)
}
cases1 <- list(c(15, 1), c(30, 2), c(60, 1), c(150, 3))
svg(file.path(fig_dir, "fig1_kah3_fit.svg"), width = 9, height = 6.2)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1, 1), family = "sans")
for (np in cases1) {
  n <- np[1]; p <- np[2]; nu <- floor(n / 3) - p - 1
  s <- replicate(20000, kah3(n, p))
  panel(s, nu, sprintf("n = %d, p = %d, v = %d", n, p, nu), "log(KaH-III statistic)")
}
dev.off()

## Figure 2: thesis simulated cut-offs vs the formula, 5% level
rd <- function(f) { t <- read.table(file.path("../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر", f),
  header = TRUE, check.names = FALSE); names(t) <- c("p", "n", "a005", "a01", "a025", "a05")
  t[] <- lapply(t, function(v) suppressWarnings(as.numeric(v))); t[complete.cases(t), ] }
tab1 <- rd("table KaH-III.txt"); tab2 <- rd("table KaH®V2-III.txt")
dfs <- read.csv("lts_effective_df.csv")
nu_star <- function(m, p) { d <- dfs[dfs$p == p, ]; d <- d[order(d$m), ]
  exp(approx(log(d$m), log(d$nu_star), log(m), rule = 2)$y) }
svg(file.path(fig_dir, "fig2_tables_vs_formula.svg"), width = 9, height = 3.9)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 1, 1))
for (w in 1:2) {
  tb <- if (w == 1) tab1 else tab2
  tb <- tb[tb$n <= 300, ]
  pc <- okabe[c(1, 4, 5)]; pp <- c(16, 17, 15)
  plot(NA, xlim = c(15, 300), ylim = c(1, if (w == 1) 60 else 15), log = "xy", las = 1,
       xlab = "Sample size n (log scale)", ylab = if (w == 1) "KaH-III 5% critical value" else "KaH-V2-III 5% critical value")
  for (p in 1:3) {
    t <- tb[tb$p == p, ]
    nn <- seq(min(t$n), 300, by = 3); m <- nn / 3
    cv <- if (w == 1) qfmax(rep(.95, length(nn)), m - p - 1)
          else sqrt(qfmax(rep(.95, length(nn)), sapply(m, nu_star, p = p)))
    lines(nn, cv, col = pc[p], lwd = 1.6)
    keep <- t$n <= 60 | t$n %% 15 == 0                # thin the dense large-n points
    points(t$n[keep], t$a05[keep], pch = pp[p], cex = .6, col = pc[p])
  }
  legend("topright", bty = "n", cex = .8, lwd = 1.6, pch = pp, col = pc,
         legend = c("p = 1", "p = 2", "p = 3"), title = "line = formula, points = table")
}
dev.off()

## Figure 3: effective degrees of freedom of the LTS scale
svg(file.path(fig_dir, "fig3_nu_star.svg"), width = 6.5, height = 3.9)
par(mar = c(4.2, 4.4, 1, 1))
plot(NA, xlim = c(5, 200), ylim = c(.5, .9), log = "x", las = 1,
     xlab = "Observations in each third (m)", ylab = "v* / v  (LTS d.f. / OLS d.f.)")
abline(h = .785, lty = 2, col = col_grey)
for (p in 1:3) { d <- dfs[dfs$p == p, ]; d <- d[order(d$m), ]
  lines(d$m, d$ratio, col = okabe[p + 1], lwd = 1.6); points(d$m, d$ratio, pch = 16, cex = .7, col = okabe[p + 1]) }
text(160, .775, "0.785", col = col_grey, cex = .8)
legend("bottomright", bty = "n", cex = .85, lwd = 1.6, col = okabe[2:4], legend = paste("p =", 1:3))
dev.off()

## Figure 4: KaH-V2-III simulated under H0 vs approximation Fmax(3, v*)
svg(file.path(fig_dir, "fig4_v2_fit.svg"), width = 9, height = 6.2)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1, 1))
for (f in c("n015_p1", "n030_p2", "n060_p1", "n300_p1")) {
  s <- readRDS(file.path("robust_full", paste0(f, ".rds")))
  n <- as.integer(substr(f, 2, 4)); p <- as.integer(substr(f, 7, 7)); ns <- nu_star(floor(n / 3), p)
  panel(s, ns, sprintf("n = %d, p = %d, v* = %.1f", n, p, ns), "log(KaH-V2-III statistic squared)", sq = TRUE)
}
dev.off()

## Figure 5: size and power at 5%
sp <- read.csv("size_power.csv")
sc_lab <- c(H0_clean = "H0\nclean", H0_t5 = "H0\nt5 errors", H0_outliers = "H0\n10% outliers",
            H1_clean = "H1\nclean", H1_outliers = "H1\n10% outliers")
svg(file.path(fig_dir, "fig5_size_power.svg"), width = 9, height = 4.3)
par(mfrow = c(1, 2), mar = c(4.6, 4.2, 1, 1))
for (nn in c(60, 150)) {
  d <- sp[sp$n == nn & sp$p == 1, ]; d <- d[match(names(sc_lab), d$scenario), ]
  M <- rbind(d$KaH3_dist, d$V2_dist, d$BP_Koenker)
  bp <- barplot(M, beside = TRUE, ylim = c(0, 1.08), las = 1, col = okabe[c(1, 4, 3)], border = NA,
                ylab = sprintf("Rejection rate at 5%%  (n = %d, p = 1)", nn), names.arg = rep("", 5))
  mtext(sc_lab, side = 1, at = colMeans(bp), line = 1.6, cex = .68)   # axis() drops labels that overlap
  abline(h = .05, lty = 2, col = col_line)
  abline(v = mean(bp[, 3:4]), col = "grey80")
  if (nn == 60) legend("topleft", bty = "n", cex = .8, fill = okabe[c(1, 4, 3)], border = NA,
                       legend = c("KaH-III (exact)", "KaH-V2-III (approx.)", "Koenker-BP"))
}
dev.off()
cat("figures written\n")
