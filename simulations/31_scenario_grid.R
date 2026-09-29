# Step 31: 2 x 4 grid of the eight n = 90 settings of the factorial study in which KaH-robust has
# the largest margin over EVERY competing test (from step 32; all eight are contaminated settings).
# Each panel shows one example data set (fixed seed, not selected): the errors against the sorting
# variable, the three parts shaded, outliers circled. The printed numbers are size-adjusted powers
# from the factorial study (1000 replications), not results for the plotted data set.
s <- readRDS("factorial_summary.rds"); res <- s$res; adj <- s$adj
rivals <- setdiff(s$tests, "KaHrobust")
rlab <- c(KaH3 = "KaH-III", GQ = "Goldfeld-Quandt", BAMSET = "BAMSET", BP = "Breusch-Pagan", White = "White",
          WhiteLTS = "White, LTS screen", MGQ = "MGQ", EvansKing = "Evans-King", WK = "Wilcox-Keselman")
out <- "../paper_stat_papers/figures"
panels <- list(
  list(shape = "U",     order = "time", p = 1, title = "U-shape over time"),
  list(shape = "U",     order = "time", p = 4, title = "U-shape over time, p = 4"),
  list(shape = "bulge", order = "time", p = 1, title = "Bulge over time"),
  list(shape = "bulge", order = "time", p = 4, title = "Bulge over time, p = 4"),
  list(shape = "bulge", order = "x1",   p = 4, title = "Bulge along x1, p = 4"),
  list(shape = "U",     order = "x1",   p = 4, title = "U-shape along x1, p = 4"),
  list(shape = "mono",  order = "time", p = 1, title = "Monotone over time"),
  list(shape = "mono",  order = "x1",   p = 4, title = "Monotone along x1, p = 4"))
sdfun <- function(shape, z) switch(shape, mono = exp(0.35 * z), U = 1 + z^2, bulge = 1 + 3 * exp(-z^2))
n <- 90; m <- n / 3
cairo_pdf(file.path(out, "Fig4.pdf"), width = 5.16, height = 3.9, family = "Arial", pointsize = 7.5)
par(mfrow = c(2, 4), mar = c(2.2, 1.9, 3.7, 0.4), mgp = c(1.1, 0.3, 0), tcl = -0.2, las = 1)
for (j in seq_along(panels)) {
  pd <- panels[[j]]
  k <- res$shape == pd$shape & res$order == pd$order & res$p == pd$p & res$cont & res$n == n
  ours <- 100 * adj$KaHrobust[k]; rv <- unlist(adj[k, rivals]); rv <- rv[!is.na(rv)]   # WK undefined for p = 4
  best <- names(which.max(rv))
  set.seed(400 + j)
  x1 <- sort(scale(runif(n))[, 1])
  z <- if (pd$order == "time") as.numeric(scale(seq_len(n))) else x1
  e <- rnorm(n, 0, sdfun(pd$shape, z))
  oi <- sample(n, round(0.1 * n)); e[oi] <- rnorm(length(oi), 0, 7)
  ylim <- c(-13, 13)
  plot(NA, xlim = range(z), ylim = ylim, xlab = if (pd$order == "time") "time" else expression(x[1]),
       ylab = "", xaxt = "n", yaxt = "n")
  axis(1, labels = FALSE); axis(2, at = c(-10, 0, 10), cex.axis = 0.9)
  br <- (z[c(m, 2 * m)] + z[c(m, 2 * m) + 1]) / 2
  rect(c(par("usr")[1], br[2]), ylim[1] - 5, c(br[1], par("usr")[2]), ylim[2] + 5, col = "grey94", border = NA)
  abline(v = br, lty = 3, col = "grey50"); abline(h = 0, col = "grey75")
  ec <- pmax(pmin(e, 12.5), -12.5)
  points(z[-oi], ec[-oi], pch = 16, cex = 0.5, col = "grey25")
  points(z[oi], ec[oi], pch = 16, cex = 0.5, col = "#D55E00")
  points(z[oi], ec[oi], pch = 1, cex = 1.6, lwd = 1.2, col = "#D55E00")
  box()
  mtext(pd$title, side = 3, line = 2.3, font = 2, cex = 0.9)
  mtext(sprintf("KaH-robust %.0f%%", ours), side = 3, line = 1.15, cex = 0.82, col = "#D55E00", font = 2)
  mtext(sprintf("best rival %.0f%% (%s)", 100 * max(rv), rlab[best]), side = 3, line = 0.2, cex = 0.72, col = "grey25")
}
dev.off()
cat("Fig4.pdf (winning-settings grid) written\n")
