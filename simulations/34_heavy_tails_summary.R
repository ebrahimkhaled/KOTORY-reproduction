# Step 34: summary of the t3 study (step 33). Size at the nominal 5% level in the 12 null cells,
# and size-adjusted power over the 36 heteroscedastic cells (threshold = 5% quantile of the
# test's own p-values in the matching null cell, as in step 29). Writes the table for the paper.
a <- readRDS("heavy_tails_all.rds")
tests <- c("KaHrobust", "KaH3", "GQ", "BAMSET", "BP", "White", "WhiteLTS", "MGQ", "EvansKing", "WK")
lab <- c(KaHrobust = "KaH-robust", KaH3 = "KaH-III", GQ = "Goldfeld--Quandt", BAMSET = "BAMSET",
         BP = "Breusch--Pagan", White = "White", WhiteLTS = "White, LTS screen", MGQ = "MGQ",
         EvansKing = "Evans--King", WK = "Wilcox--Keselman")
key <- function(d) paste(d$order, d$p, d$n)
null <- a[a$shape == "none", ]; alt <- a[a$shape != "none", ]
cat("replications per null cell:", range(table(key(null))), " per alt cell:", range(table(paste(alt$shape, key(alt)))), "\n")

# size (%) per null cell, then range and mean over the 12 cells
sz <- do.call(cbind, lapply(c(tests, "KaHrobustBoot"), function(t)
  tapply(null[[t]], key(null), function(p) 100 * mean(p <= .05, na.rm = TRUE))))
colnames(sz) <- c(tests, "KaHrobustBoot")
na_boot <- sum(is.na(null$KaHrobustBoot))
cat("bootstrap NA:", na_boot, "\n")

# size-adjusted power per heteroscedastic cell
cells <- unique(alt[, c("shape", "order", "p", "n")])
adj <- t(sapply(seq_len(nrow(cells)), function(i) {
  c0 <- cells[i, ]; h <- alt[alt$shape == c0$shape & key(alt) == key(c0), ]; n0 <- null[key(null) == key(c0), ]
  sapply(tests, function(t) { if (all(is.na(h[[t]]))) return(NA)
    thr <- quantile(n0[[t]], .05, na.rm = TRUE, type = 1); 100 * mean(h[[t]] <= thr, na.rm = TRUE) })
}))
# KaH-III: size-adjust on the statistic (step 38), because its p-values underflow under heavy tails
t1 <- readRDS("heavy_kah3_T1.rds")
t1key <- function(d) paste(d$shape, d$order, d$p, d$n, round(d$fingerprint, 6))
a$KaH3T <- t1$T1[match(t1key(a), t1key(t1))]
stopifnot(!anyNA(a$KaH3T))
null$KaH3T <- a$KaH3T[a$shape == "none"]; alt$KaH3T <- a$KaH3T[a$shape != "none"]
adj[, "KaH3"] <- sapply(seq_len(nrow(cells)), function(i) {
  c0 <- cells[i, ]; T0 <- null$KaH3T[key(null) == key(c0)]
  T1 <- alt$KaH3T[alt$shape == c0$shape & key(alt) == key(c0)]
  100 * mean(T1 >= sort(T0, decreasing = TRUE)[ceiling(0.05 * length(T0))]) })
adj <- cbind(cells, adj)
mean_adj <- colMeans(adj[, tests, drop = FALSE], na.rm = TRUE)
by_shape <- sapply(tests, function(t) tapply(adj[[t]], adj$shape, mean, na.rm = TRUE))
wins <- table(factor(apply(adj[, tests, drop = FALSE], 1, function(r) names(which.max(r))), levels = tests))
out <- data.frame(test = lab[tests], size_min = apply(sz[, tests, drop = FALSE], 2, min, na.rm = TRUE), size_max = apply(sz[, tests, drop = FALSE], 2, max, na.rm = TRUE),
                  size_mean = colMeans(sz[, tests, drop = FALSE], na.rm = TRUE), adj_mean = mean_adj,
                  t(by_shape[c("mono", "U", "bulge"), ]), wins = as.integer(wins))
out <- out[order(-out$adj_mean), ]
print(format(out, digits = 3), row.names = FALSE)
cat("\nKaH-robust bootstrap size: min", min(sz[, "KaHrobustBoot"]), "max", max(sz[, "KaHrobustBoot"]),
    "mean", mean(sz[, "KaHrobustBoot"]), "\n")
print(round(sz[, c("KaHrobust", "KaHrobustBoot", "KaH3"), drop = FALSE], 1))
saveRDS(list(size = sz, adj = adj, summary = out), "heavy_tails_summary.rds")

# LaTeX table (booktabs; narrowed with font size and tabcolsep, never resizebox)
f1 <- function(x) ifelse(is.na(x), "--", sprintf("%.0f", floor(x + 0.5)))   # round half up (sprintf rounds 52.5 to 52)
rows <- sprintf("%s & %s--%s & %s & %s & %s & %s & %s & %d \\\\", out$test, f1(out$size_min), f1(out$size_max),
                f1(out$size_mean), f1(out$adj_mean), f1(out$mono), f1(out$U), f1(out$bulge), out$wins)
bt <- sz[, "KaHrobustBoot"]
tex <- c("\\begin{table}[!htbp]", "\\centering",
  sprintf("\\caption{Heavy-tailed errors: $t_3$ errors scaled to unit variance and no planted outliers, with the factorial design of Section~\\ref{sec:factorial} otherwise unchanged (two ordering variables, $p\\in\\{1,4\\}$, $n\\in\\{45,90,150\\}$); 1000 replications per cell. Size: rejection rate (\\%%) at the nominal 5\\%% level over the 12 null cells. Power: mean size-adjusted power (\\%%) over the 36 heteroscedastic cells, overall and by variance shape. Wins: cells in which the test had the highest size-adjusted power. Rows are sorted by mean size-adjusted power; the Wilcox--Keselman test, whose quantile-regression bootstrap is slow, was run for $p=1$ only. The bootstrap reference of \\KaHR{} ($B=199$) had size %.1f to %.1f\\%% (mean %.1f\\%%) in the same null cells}\\label{tab:heavy}",
          min(bt), max(bt), mean(bt)),
  "\\small\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{@{}lrrrrrrr@{}}", "\\toprule",
  " & \\multicolumn{2}{c}{Size} & \\multicolumn{4}{c}{Size-adjusted power} & \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-7}",
  "Test & range & mean & all & monotone & U & bulge & Wins \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, "../paper_stat_papers/tables/tab_heavy.tex")
cat("tab_heavy.tex written\n")
