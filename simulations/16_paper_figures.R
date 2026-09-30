# Step 16: figures for the Statistical Papers manuscript.
# Springer rules: vector PDF, sans-serif lettering 8-10 pt at final size, no titles in the
# artwork, legible in the review column (sn-jnl text width ~ 131 mm = 5.16 in). Writes Fig1 and Fig2;
# Fig3 (size heatmap) is written by step 22.

source("fmax_exact.R")
out <- "../paper_stat_papers/figures"; dir.create(out, showWarnings = FALSE)
W <- 5.16                                   # review-column width, inches
okabe <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73", blue = "#0072B2",
           verm = "#D55E00", pink = "#CC79A7", grey = "grey45")
setup <- function(file, h, mfrow = c(1, 1)) {
  cairo_pdf(file.path(out, file), width = W, height = h, family = "Arial", pointsize = 9)
  par(mfrow = mfrow, mar = c(3.4, 3.6, 0.8, 0.6), mgp = c(2.1, 0.6, 0), las = 1, tcl = -0.3)
}
dens_panel <- function(s, nu, xlab, lab) {
  z <- log(s); br <- seq(0, quantile(z, 0.995), length.out = 40); z2 <- z[z < max(br)]
  hist(z2, breaks = br, freq = FALSE, col = "grey85", border = "white", main = "", xlab = xlab, ylab = "Density")
  u <- seq(0.003, max(br), length.out = 200); c0 <- exp(u); h <- 1e-4 * c0
  d <- (pfmax(c0 + h, rep(nu, 200)) - pfmax(c0 - h, rep(nu, 200))) / (2 * h)
  lines(u, d * c0 * length(z) / length(z2), lwd = 1.8, col = okabe["verm"])
  legend("topright", legend = lab, bty = "n", cex = 0.9)
}

## Fig 1: null fits (a) KaH-III exact, (b) KaH-robust approximation
set.seed(2421974)
kah3 <- function(n, p) {
  x <- scale(matrix(runif(n * p), n, p)); x <- x[order(x[, 1]), , drop = FALSE]
  y <- drop(1 + x %*% rep(1, p) + rnorm(n)); m <- floor(n / 3)
  idx <- list(1:m, (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
  ms <- sapply(idx, function(i) { f <- lm.fit(cbind(1, x[i, , drop = FALSE]), y[i]); sum(f$residuals^2) / (length(i) - p - 1) })
  max(ms) / min(ms)
}
s1 <- replicate(20000, kah3(30, 2))
s2 <- readRDS("robust_a75/n060_p1.rds")^2          # step 37: alpha = 0.75, package LTS fits
load("../KOTORY/R/sysdata.rda")
nu2 <- with(.nu_table[.nu_table$alpha == 0.75 & .nu_table$p == 1 & .nu_table$m == 20, ], nu_star)
setup("Fig2.pdf", 2.5, c(1, 2))
dens_panel(s1, 7, expression(log~T[1]), c("(a) n = 30, p = 2", expression(F[max](3*","~7))))
dens_panel(s2, nu2, expression(log~T[2]^2), c("(b) n = 60, p = 1", bquote(F[max](3*","~.(round(nu2, 1))))))
dev.off()

## Fig 2: effective d.f. ratio against part size, with the asymptotic limit r(alpha)
r_theory <- function(a) {
  xi <- qchisq(a, 1); f <- function(y) dchisq(y, 1)
  T <- integrate(function(y) y * f(y), 0, xi)$value / a
  m1 <- integrate(function(y) (xi - y) * f(y), 0, xi)$value
  m2 <- integrate(function(y) (xi - y)^2 * f(y), 0, xi)$value
  2 * T^2 * a^2 / (m2 - m1^2)
}
setup("Fig1.pdf", 2.9)
cols <- c(okabe["blue"], okabe["green"], okabe["verm"]); pchs <- c(16, 17, 15)
plot(NA, xlim = c(5, 200), ylim = c(0.2, 0.95), log = "x", xlab = "Observations per part, m",
     ylab = expression(nu^"*" / nu))
for (k in 1:3) {
  a <- c(0.5, 0.75, 0.9)[k]
  t <- aggregate(ratio ~ m, data = .nu_table[.nu_table$alpha == a, ], FUN = mean)
  lines(t$m, t$ratio, col = cols[k], lwd = 1.2)
  points(t$m, t$ratio, col = cols[k], pch = pchs[k], cex = 0.8)
  abline(h = r_theory(a), col = cols[k], lty = 2)
}
legend("top", bty = "n", cex = 0.9, col = cols, pch = pchs, lwd = 1.2, horiz = TRUE,
       legend = c(expression(alpha == 0.50), expression(alpha == 0.75), expression(alpha == 0.90)))
legend("bottomleft", bty = "n", cex = 0.85, lty = c(1, 2), col = "grey30",
       legend = c("simulated (mean over p = 1,...,5)", expression("limit " * r(alpha))))
dev.off()

cat("figures written to", normalizePath(out), "\n")
