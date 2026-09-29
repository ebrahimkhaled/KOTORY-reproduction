# Step 2: (a) sanity-check the exact Fmax code against a fresh Monte Carlo of
# the KaH-III statistic itself (thesis design: x ~ U(0,1) standardised, sort
# by x1, 3 equal groups, OLS in each), (b) compare with the thesis table.

source("fmax_exact.R")
set.seed(2421974)

# (a) Hartley's published 5%/1% values, k = 3 (Pearson & Hartley Table 31):
#     df=4: 15.5 / 37 ; df=10: 4.85 / 7.4 ; df=30: 2.40 / 3.0 (approx.)
cat("Exact k=3, 95%/99% quantiles at df = 4, 10, 30:\n")
print(round(rbind(q95 = qfmax(rep(.95, 3), c(4, 10, 30)),
                  q99 = qfmax(rep(.99, 3), c(4, 10, 30))), 3))

kah3 <- function(n, p) {
  x <- scale(matrix(runif(n * p), n, p))
  x <- x[order(x[, 1]), , drop = FALSE]
  y <- drop(1 + x %*% rep(1, p) + rnorm(n))
  m <- floor(n / 3)
  idx <- list(1:m, (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
  ms <- sapply(idx, function(i) {
    f <- lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])
    sum(f$residuals^2) / (length(i) - p - 1)
  })
  max(ms) / min(ms)
}

R <- if (nzchar(Sys.getenv("SKIP_MC"))) 0 else 1e5
if (R > 0) {
cat("\nMonte Carlo (R =", R, ") vs exact, KaH-III:\n")
res <- do.call(rbind, lapply(list(c(15, 1), c(30, 2), c(60, 3), c(31, 1), c(100, 2)), function(np) {
  n <- np[1]; p <- np[2]; nu <- floor(n / 3) - p - 1
  s <- replicate(R, kah3(n, p))
  data.frame(n = n, p = p, nu = nu,
             mc_q95 = quantile(s, .95), ex_q95 = qfmax(.95, nu),
             mc_q99 = quantile(s, .99), ex_q99 = qfmax(.99, nu),
             ks_p = suppressWarnings(ks.test(s, function(q) pfmax(q, nu))$p.value))
}))
print(res, digits = 4, row.names = FALSE)
}
cat("nu range in thesis table:", range(floor(read.table(file.path("../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر", "table KaH-III.txt"), header = TRUE)$n / 3) - 2), "\n")

# (b) thesis table vs exact
src <- "../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر"
tab <- read.table(file.path(src, "table KaH-III.txt"), header = TRUE, check.names = FALSE)
names(tab) <- c("p", "n", "a005", "a01", "a025", "a05")
tab[] <- lapply(tab, function(v) suppressWarnings(as.numeric(v)))
tab <- tab[complete.cases(tab), ]
tab$nu <- floor(tab$n / 3) - (tab$p + 1)
al <- c(a005 = .005, a01 = .01, a025 = .025, a05 = .05)
for (a in names(al)) {
  tab[[paste0("ex_", a)]] <- qfmax(rep(1 - al[a], nrow(tab)), tab$nu)
  tab[[paste0("size_", a)]] <- 1 - pfmax(tab[[a]], tab$nu)   # true size of thesis cut-off
}
cat("\nTrue size of the thesis cut-offs under the exact law (", nrow(tab), "rows):\n")
print(sapply(names(al), function(a) summary(tab[[paste0("size_", a)]])), digits = 3)
cat("\nWorst rows (|size/nominal - 1| largest) at 1%:\n")
tab$dev01 <- tab$size_a01 / .01 - 1
print(head(tab[order(-abs(tab$dev01)), c("p", "n", "nu", "a01", "ex_a01", "size_a01")], 8), digits = 4, row.names = FALSE)
write.csv(tab, "kah3_thesis_vs_exact.csv", row.names = FALSE)
