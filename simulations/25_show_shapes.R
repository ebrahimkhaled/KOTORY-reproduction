# Print the shape study (step 24) as % tables: one block per shape, rows = tests, columns = outliers x n
a <- read.csv("shapes_rejection.csv")
meth <- c("KaHrobust", "KaH3", "GQ", "BAMSET", "BP", "White", "EvansKing", "MGQ", "WhiteLTS", "WK")
for (sh in c("none", "U", "bulge", "mono")) {
  d <- a[a$shape == sh, ]; d <- d[order(match(d$outl, c("none", "centre", "random")), d$n), ]
  t <- t(round(100 * as.matrix(d[meth]), 1)); colnames(t) <- paste0(substr(d$outl, 1, 4), d$n)
  cat("\n== shape:", sh, if (sh == "none") "(SIZE, % rejected, should be 5)" else "(POWER, %)", "\n")
  print(t)
}
