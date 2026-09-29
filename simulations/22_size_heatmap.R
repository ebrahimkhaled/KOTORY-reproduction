# Step 22: size heatmap (replaces the size table). Cell = empirical size (%) at the nominal 5%
# level; colour = distance from 5% (blue conservative, pale near nominal, orange-red inflated).
a <- read.csv("all_tests_rejection.csv", check.names = FALSE)
out <- "../paper_stat_papers/figures"
tests <- c(V2_a75 = "KaH-robust", V2_a75_bs = "KaH-robust, bootstrap", KaH3 = "KaH-III",
           WK_2006 = "Wilcox-Keselman", BRW_White = "White after LTS screen", BP_Koenker = "Breusch-Pagan",
           White = "White", Zhou_2015 = "Zhou et al.", BF_3parts = "Brown-Forsythe, 3 parts",
           HMC = "Harrison-McCabe", EvansKing = "Evans-King", GQ = "Goldfeld-Quandt",
           BAMSET_3 = "BAMSET", MGQ_Rana2008 = "MGQ", V2_a90 = "KaH-robust, α = 0.90")
scen <- list(H0_clean = "Normal errors", H0_t5 = expression(italic(t)[5] * " errors"), H0_outliers = "10% outliers")
cells <- list(c(30, 1), c(60, 1), c(150, 1), c(90, 2)); clab <- c("30", "60", "150", "90*")
M <- do.call(cbind, lapply(names(scen), function(s) sapply(cells, function(np)
  unlist(a[a$scenario == s & a$n == np[1] & a$p == np[2], names(tests)]))))
M <- 100 * M                                     # rows: tests; cols: scenario x cell
col_of <- function(v) {                          # diverging palette around 5%
  br <- c(-Inf, 2, 3.5, 6.5, 8, 10, 15, 25, 50, Inf)
  pal <- c("#2166AC", "#92C5DE", "#F7F7F7", "#FDDBC7", "#F4A582", "#E08050", "#D6604D", "#B2182B", "#67001F")
  pal[findInterval(v, br, left.open = TRUE)]
}
cairo_pdf(file.path(out, "Fig3.pdf"), width = 5.16, height = 4.3, family = "Arial", pointsize = 8)
par(mar = c(2.2, 10.2, 2.6, 0.4), las = 1, xpd = NA)
nr <- nrow(M); nc <- ncol(M)
plot(NA, xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "")
for (i in 1:nr) for (j in 1:nc) {
  v <- M[i, j]
  rect(j - 1, nr - i, j, nr - i + 1, col = col_of(v), border = "white", lwd = 0.8)
  text(j - 0.5, nr - i + 0.5, formatC(v, format = "f", digits = 1), cex = 0.78,
       col = if (v > 15 || v < 2) "white" else "grey10")
}
axis(2, at = nr:1 - 0.5, labels = tests, tick = FALSE, line = -0.6, cex.axis = 0.95,
     font = 1)
axis(1, at = 1:nc - 0.5, labels = rep(clab, 3), tick = FALSE, line = -0.9, cex.axis = 0.85)
mtext("n (* p = 2)", side = 1, line = 1.1, cex = 0.8, adj = 0)
for (k in 1:3) {
  segments(4 * (k - 1) + 0.1, nr + 0.35, 4 * k - 0.1, nr + 0.35, lwd = 0.8)
  text(4 * (k - 1) + 2, nr + 0.8, scen[[k]], cex = 0.95)
}
abline(v = c(4, 8), col = "grey20", lwd = 1.2, xpd = FALSE)
rect(0, nr - 2, nc, nr, border = "grey10", lwd = 1.4, xpd = FALSE)   # frame the proposed tests
dev.off()
cat("Fig3.pdf (size heatmap) written\n")
