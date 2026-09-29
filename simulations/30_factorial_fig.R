# Step 30: headline figure of the factorial study (Fig5). (a) average size-adjusted power over all
# heteroscedastic cells, ranked, KaH-robust highlighted; (b) the same average by level of each factor.
s <- readRDS("factorial_summary.rds"); adj <- s$adj; res <- s$res
lab <- c(KaHrobust = "KaH-robust", KaH3 = "KaH-III", GQ = "Goldfeld-Quandt", BAMSET = "BAMSET",
         BP = "Breusch-Pagan", White = "White", WhiteLTS = "White after LTS screen", MGQ = "MGQ",
         EvansKing = "Evans-King", WK = "Wilcox-Keselman*")
out <- "../paper_stat_papers/figures"
avg <- sort(sapply(names(lab), function(t) mean(adj[[t]], na.rm = TRUE)))
fac <- list(shape = c(mono = "Monotone", U = "U-shape", bulge = "Bulge"),
            order = c(x1 = "Sorted by regressor", time = "Sorted by time"),
            p = c(`1` = "p = 1", `4` = "p = 4"),
            cont = c(`FALSE` = "No outliers", `TRUE` = "10% outliers"))
grp <- do.call(rbind, lapply(names(fac), function(f) do.call(rbind, lapply(names(fac[[f]]), function(l) {
  k <- as.character(res[[f]]) == l
  data.frame(level = fac[[f]][[l]], t(sapply(names(lab), function(t) mean(adj[[t]][k], na.rm = TRUE))), check.names = FALSE)
}))))
cairo_pdf(file.path(out, "Fig5.pdf"), width = 5.16, height = 4.1, family = "Arial", pointsize = 8.5)
layout(rbind(1:2, c(3, 3)), widths = c(1, 1.25), heights = c(1, 0.16))
par(mar = c(3.2, 9.2, 1.4, 0.6), las = 1, mgp = c(1.9, 0.5, 0), tcl = -0.25)
cols <- ifelse(names(avg) == "KaHrobust", "#D55E00", ifelse(names(avg) == "KaH3", "#0072B2", "grey70"))
b <- barplot(100 * avg, horiz = TRUE, names.arg = lab[names(avg)], col = cols, border = NA, xlim = c(0, 100),
             xlab = "Mean size-adjusted power (%)", cex.names = 0.95)
text(100 * avg + 1.5, b, sprintf("%.0f", 100 * avg), adj = 0, cex = 0.9,
     font = ifelse(names(avg) == "KaHrobust", 2, 1))
mtext("(a)", side = 3, line = 0.2, adj = -0.55, font = 2)
par(mar = c(3.2, 8.2, 1.4, 0.6))
show <- c("KaHrobust", "KaH3", "WhiteLTS", "White", "BP", "MGQ")
pc <- c(KaHrobust = "#D55E00", KaH3 = "#0072B2", WhiteLTS = "#E69F00", White = "#009E73", BP = "#56B4E9", MGQ = "#CC79A7")
pp <- c(KaHrobust = 16, KaH3 = 1, WhiteLTS = 0, White = 4, BP = 3, MGQ = 5)
nl <- nrow(grp)
plot(NA, xlim = c(0, 100), ylim = c(0.5, nl + 0.5), yaxt = "n", xlab = "Mean size-adjusted power (%)", ylab = "")
axis(2, at = nl:1, labels = grp$level, tick = FALSE, cex.axis = 0.95)
abline(h = c(nl - 2.5, nl - 4.5, nl - 6.5), col = "grey75")
for (t in show) points(100 * grp[[t]], nl:1, pch = pp[t], col = pc[t], cex = if (t == "KaHrobust") 1.4 else 1)
mtext("(b)", side = 3, line = 0.2, adj = -0.5, font = 2)
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", ncol = 3, bty = "n", cex = 0.9, pch = pp[show], col = pc[show], pt.cex = c(1.4, rep(1, 5)),
       legend = lab[show])
dev.off()
cat("Fig5.pdf written\n")
