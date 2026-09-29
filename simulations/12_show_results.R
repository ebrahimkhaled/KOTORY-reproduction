# Print the step-11 rejection rates (% at the 5% level), one table per scenario
a <- read.csv("all_tests_rejection.csv", check.names = FALSE)
m <- setdiff(names(a), c("n", "p", "scenario"))
for (sc in unique(a$scenario)) {
  b <- a[a$scenario == sc, ]
  t <- t(round(100 * as.matrix(b[m]), 1))
  colnames(t) <- paste0("n", b$n, "_p", b$p)
  cat("\n==", sc, "(% rejected at 5%)\n")
  print(t)
}
