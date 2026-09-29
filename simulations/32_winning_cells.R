# Step 32: in which n = 90 settings of the factorial study does KaH-robust beat EVERY rival on
# size-adjusted power? (WK only defined for p = 1.)
s <- readRDS("factorial_summary.rds"); res <- s$res; adj <- s$adj
rivals <- setdiff(s$tests, "KaHrobust")
k <- which(res$n == 90)
tab <- do.call(rbind, lapply(k, function(i) {
  r <- unlist(adj[i, rivals]); b <- names(which.max(r))
  data.frame(shape = res$shape[i], order = res$order[i], p = res$p[i], cont = res$cont[i],
             KaHrobust = round(100 * adj$KaHrobust[i]), best_rival = b, rival = round(100 * max(r, na.rm = TRUE)),
             margin = round(100 * (adj$KaHrobust[i] - max(r, na.rm = TRUE))))
}))
tab <- tab[order(-tab$margin), ]
print(tab, row.names = FALSE)
cat("\nKaH-robust best in", sum(tab$margin > 0), "of", nrow(tab), "n = 90 settings\n")
