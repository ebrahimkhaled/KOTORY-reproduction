# Step 20: size-adjusted power under contamination. For each test and (n, p), the rejection
# threshold is the 5% quantile of its p-values under H0 with 10% outliers (so its size under
# contamination is exactly 5%); power is then the rate below that threshold under H1 with 10%
# outliers. Uses the per-replication p-values stored by step 11 (same data for every test).
files <- list.files("all_tests", "_ch\\d+\\.rds$", full.names = TRUE)
all <- do.call(rbind, lapply(files, readRDS))
tests <- c("KaH3", "V2_a75", "V2_a75_bs", "V2_a90", "GQ", "BP_Koenker", "White", "MGQ_Rana2008",
           "BRW_White", "WK_2006", "Zhou_2015", "EvansKing", "HMC", "BF_3parts")
cells <- unique(all[c("n", "p")])
res <- do.call(rbind, lapply(seq_len(nrow(cells)), function(k) {
  n <- cells$n[k]; p <- cells$p[k]
  h0 <- all[all$n == n & all$p == p & all$scenario == "H0_outliers", ]
  h1 <- all[all$n == n & all$p == p & all$scenario == "H1_outliers", ]
  data.frame(n = n, p = p, test = tests,
             raw = sapply(tests, function(t) mean(h1[[t]] < 0.05)),
             adjusted = sapply(tests, function(t) mean(h1[[t]] <= quantile(h0[[t]], 0.05, type = 1))))
}))
res$cell <- paste0("n", res$n, "_p", res$p)
out <- reshape(res[, c("test", "cell", "adjusted")], idvar = "test", timevar = "cell", direction = "wide")
names(out) <- sub("adjusted\\.", "", names(out))
out <- out[, c("test", "n30_p1", "n60_p1", "n150_p1", "n90_p2")]
cat("Size-adjusted power (%) under H1 with 10% outliers (threshold = 5% quantile under H0 with outliers):\n")
print(cbind(out[1], round(100 * out[-1], 1)), row.names = FALSE)
write.csv(res, "size_adjusted_power.csv", row.names = FALSE)
