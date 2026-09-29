# Step 36: Fig6. Every heteroscedastic cell of the factorial study (step 29): size-adjusted power of the
# KaH test chosen by the paper's usage rule (KaH-III on clean data, KaH-robust when outliers are
# possible) against the best of the eight other tests in that cell, a different test from cell to cell.
s <- readRDS("factorial_summary.rds"); adj <- s$adj; res <- s$res
out <- Sys.getenv("FIG_OUT", "../paper_stat_papers/figures")
rivals <- setdiff(s$tests, c("KaHrobust", "KaH3"))
ours <- 100 * ifelse(res$cont, adj$KaHrobust, adj$KaH3)
best <- 100 * apply(adj[, rivals], 1, max, na.rm = TRUE)
win <- ours > best
cat(sprintf("cells %d, rule ahead in %d (outliers %d/%d, clean %d/%d); means %.1f vs %.1f\n", length(ours), sum(win),
            sum(win & res$cont), sum(res$cont), sum(win & !res$cont), sum(!res$cont), mean(ours), mean(best)))

cairo_pdf(file.path(out, "Fig6.pdf"), width = 5.16, height = 3.2, family = "Arial", pointsize = 8.5)
layout(matrix(1:2, 1), widths = c(1, 0.78))
par(mar = c(3.2, 3.4, 0.6, 0.6), mgp = c(1.9, 0.5, 0), las = 1, tcl = -0.25, pty = "s")
plot(NA, xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i",
     xlab = "Best of the eight other tests (%)", ylab = "KaH test chosen by the rule (%)")
polygon(c(0, 100, 0), c(0, 100, 100), col = "#D55E0014", border = NA)
abline(0, 1, col = "grey40", lty = 2)
points(best[!res$cont], ours[!res$cont], pch = 1, col = "#0072B2")
points(best[res$cont], ours[res$cont], pch = 16, col = "#D55E00")
points(mean(best), mean(ours), pch = 23, bg = "black", col = "white", cex = 1.6)
box()
par(mar = c(3.2, 0.3, 0.6, 0.3), pty = "m"); plot.new()
legend("left", bty = "n", pch = c(16, 1, 23), col = c("#D55E00", "#0072B2", "white"),
       pt.bg = c(NA, NA, "black"), pt.cex = c(1, 1, 1.4), y.intersp = 1.35,
       legend = c(sprintf("10%% outliers, KaH-robust:\nahead in %d of %d settings", sum(win & res$cont), sum(res$cont)),
                  sprintf("no outliers, KaH-III:\nahead in %d of %d settings", sum(win & !res$cont), sum(!res$cont)),
                  sprintf("mean over all %d settings:\n%.0f%% against %.0f%%", length(ours), mean(ours), mean(best))),
       title = "Each point is one setting", title.adj = 0)
dev.off()
cat("Fig6.pdf (usage rule against the best other test) written\n")
