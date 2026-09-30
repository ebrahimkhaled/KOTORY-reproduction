# Step 29: summarise the factorial study (step 28).
#  - size: rejection rate at 5% in every H0 cell (shape = none)
#  - size-adjusted power: in every H1 cell, reject below the 5% quantile of the test's p-values in
#    the matching H0 cell (same order, p, contamination, n)
#  - averages over all H1 cells and by factor; share of H1 cells in which each test is best
all <- readRDS("factorial_all.rds")
tests <- c("KaHrobust", "KaH3", "GQ", "BAMSET", "BP", "White", "WhiteLTS", "MGQ", "EvansKing", "WK")
key <- function(d) paste(d$order, d$p, d$cont, d$n)
h0 <- all[all$shape == "none", ]; h1 <- all[all$shape != "none", ]
# size per H0 cell
size <- aggregate(h0[tests] < .05, by = h0[c("order", "p", "cont", "n")], FUN = function(v) mean(v, na.rm = TRUE))
# thresholds per H0 cell
thr <- lapply(split(h0, key(h0)), function(d) sapply(tests, function(t)
  if (all(is.na(d[[t]]))) NA_real_ else unname(quantile(d[[t]], 0.05, type = 1, na.rm = TRUE))))
cellsH1 <- unique(h1[c("shape", "order", "p", "cont", "n")])
res <- do.call(rbind, lapply(seq_len(nrow(cellsH1)), function(i) {
  out <- cellsH1[i, ]
  d <- h1[h1$shape == out$shape & h1$order == out$order & h1$p == out$p & h1$cont == out$cont & h1$n == out$n, ]
  th <- thr[[key(out)]]
  for (t in tests) {
    out[[paste0(t, ".nom")]] <- mean(d[[t]] < .05, na.rm = TRUE)
    tt <- if (is.null(th)) NA_real_ else th[[t]]        # NULL if the matching H0 cell is missing
    out[[paste0(t, ".adj")]] <- if (is.na(tt)) NA_real_ else mean(d[[t]] <= tt, na.rm = TRUE)
  }
  out
}))
# KaH-III: its stored p-values underflow to 0 in some contaminated null samples, so size-adjust it on
# the statistic itself (step 38): reject when T1 is at least the ceiling(0.05 N)-th largest null value
t1 <- readRDS("factorial_kah3_T1.rds")
t1key <- function(d) paste(d$shape, d$order, d$p, d$cont, d$n, round(d$fingerprint, 6))
all$KaH3T <- t1$T1[match(t1key(all), t1key(t1))]
stopifnot(!anyNA(all$KaH3T))
h0T <- all[all$shape == "none", ]; h1T <- all[all$shape != "none", ]
for (i in seq_len(nrow(res))) {
  c0 <- res[i, ]
  T0 <- h0T$KaH3T[key(h0T) == key(c0)]
  T1 <- h1T$KaH3T[h1T$shape == c0$shape & key(h1T) == key(c0)]
  res$KaH3.adj[i] <- mean(T1 >= sort(T0, decreasing = TRUE)[ceiling(0.05 * length(T0))])
}
adj <- res[paste0(tests, ".adj")]; names(adj) <- tests
nom <- res[paste0(tests, ".nom")]; names(nom) <- tests
cat("\nSIZE (%), range over the 24 H0 cells:\n")
print(round(100 * t(sapply(tests, function(t) c(min = min(size[[t]], na.rm = TRUE), median = median(size[[t]], na.rm = TRUE),
                                                   max = max(size[[t]], na.rm = TRUE)))), 1))
p1 <- res$p == 1
cat("\nAVERAGE SIZE-ADJUSTED POWER (%) over all 72 H1 cells (WK: p = 1 cells only):\n")
avg <- sapply(tests, function(t) mean(adj[[t]], na.rm = TRUE)); print(round(100 * sort(avg, decreasing = TRUE), 1))
cat("\nAVERAGE SIZE-ADJUSTED POWER (%) over the 36 p = 1 cells (all tests comparable):\n")
avg1 <- sapply(tests, function(t) mean(adj[[t]][p1], na.rm = TRUE)); print(round(100 * sort(avg1, decreasing = TRUE), 1))
best <- apply(adj[p1, ], 1, function(v) names(v)[which.max(v)])
cat("\nShare of the 36 p = 1 cells in which each test has the highest size-adjusted power:\n")
print(round(100 * prop.table(table(factor(best, levels = tests))), 1))
bestall <- apply(adj[, setdiff(tests, "WK")], 1, function(v) names(v)[which.max(v)])
cat("\nShare of all 72 cells (WK excluded) won:\n")
print(round(100 * prop.table(table(factor(bestall, levels = setdiff(tests, "WK")))), 1))
cat("\nAverage size-adjusted power (%) by factor:\n")
for (fct in c("shape", "order", "p", "cont", "n")) {
  a <- aggregate(adj, by = res[fct], FUN = function(v) round(100 * mean(v, na.rm = TRUE), 1))
  print(a, row.names = FALSE)
}
saveRDS(list(size = size, res = res, adj = adj, nom = nom, tests = tests), "factorial_summary.rds")
