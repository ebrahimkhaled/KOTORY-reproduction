# Step 23: power-curve figure (Fig7) from step 21. Rejection rate against the strength gamma of
# heteroscedasticity; the gamma = 0 column is the size. KaH-robust is the highlighted series.
a <- read.csv("power_curves.csv")
out <- "../paper_stat_papers/figures"
ser <- list(
  KaHrobust = list(lab = "KaH-robust",            col = "#D55E00", lty = 1, lwd = 2.8, pch = 16, cex = 0.95),
  KaH3      = list(lab = "KaH-III",               col = "#0072B2", lty = 2, lwd = 1.4, pch = 1,  cex = 0.75),
  GQ        = list(lab = "Goldfeld-Quandt",       col = "#56B4E9", lty = 3, lwd = 1.4, pch = 2,  cex = 0.75),
  MGQ       = list(lab = "MGQ",                   col = "#CC79A7", lty = 4, lwd = 1.4, pch = 5,  cex = 0.75),
  WhiteLTS  = list(lab = "White after LTS screen", col = "#E69F00", lty = 1, lwd = 1.4, pch = 0,  cex = 0.75),
  WK        = list(lab = "Wilcox-Keselman",       col = "grey35",  lty = 5, lwd = 1.4, pch = 6,  cex = 0.75),
  BP        = list(lab = "Breusch-Pagan",         col = "#009E73", lty = 1, lwd = 1.4, pch = 3,  cex = 0.75),
  White     = list(lab = "White",                 col = "#009E73", lty = 2, lwd = 1.4, pch = 4,  cex = 0.75))
draw_order <- c("WK", "White", "BP", "WhiteLTS", "MGQ", "GQ", "KaH3", "KaHrobust")

cairo_pdf(file.path(out, "Fig7.pdf"), width = 5.16, height = 5.6, family = "Arial", pointsize = 9)
layout(rbind(1:2, 3:4, c(5, 5)), heights = c(1, 1, 0.36))
par(mgp = c(1.9, 0.5, 0), tcl = -0.25, las = 1, cex = 1)   # layout() shrinks text to 83%; undo it
for (n in c(60, 150)) for (cont in c(FALSE, TRUE)) {
  par(mar = c(3.0, if (!cont) 3.6 else 1.2, 1.6, 0.6))
  d <- a[a$n == n & a$cont == cont, ]; d <- d[order(d$g), ]
  plot(NA, xlim = c(0, 0.7), ylim = c(0, 1), xlab = expression(gamma ~ "(strength of heteroscedasticity)"),
       ylab = if (!cont) "Rejection rate (%)" else "", yaxt = "n")
  axis(2, at = seq(0, 1, 0.25), labels = if (!cont) c("0", "25", "50", "75", "100") else FALSE)
  rect(-0.03, -0.05, 0.03, 1.05, col = "grey93", border = NA)
  text(0, 1.0, "size", cex = 0.95, col = "grey40")
  abline(h = 0.05, lty = 2, col = "grey55")
  abline(h = seq(0.25, 1, 0.25), col = "grey92", lwd = 0.6)
  for (s in draw_order) {
    z <- ser[[s]]
    lines(d$g, d[[s]], col = z$col, lty = z$lty, lwd = z$lwd)
    points(d$g, d[[s]], col = z$col, pch = z$pch, cex = z$cex)
  }
  mtext(sprintf("n = %d, %s", n, if (cont) "10% outliers" else "normal errors"), side = 3, line = 0.25, cex = 1,
        font = 2)
  if (cont && n == 60) {                         # read the contaminated panel for the reader
    text(0.005, 0.87, "false alarms at γ = 0", cex = 1, font = 3, col = "#0072B2", adj = 0)
    text(0.36, 0.29, "no power under outliers", cex = 1, font = 3, col = "#009E73")
  }
  box()
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", ncol = 3, bty = "n", cex = 1,
       legend = sapply(ser, `[[`, "lab"), col = sapply(ser, `[[`, "col"), lty = sapply(ser, `[[`, "lty"),
       lwd = sapply(ser, `[[`, "lwd"), pch = sapply(ser, `[[`, "pch"), seg.len = 2.4)
dev.off()
cat("Fig7.pdf (power curves) written\n")
