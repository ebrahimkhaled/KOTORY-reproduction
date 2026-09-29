# KaH-III statistic = max/min of the residual mean squares of 3 equal-size
# OLS sub-regressions. Under H0 (normal, homoscedastic errors) the three
# mean squares are independent sigma^2 chi2_nu / nu, nu = m - (p+1),
# m = floor(n/3)  ->  statistic ~ Hartley's Fmax(k = 3, nu) exactly.
# Compare the exact quantiles with the thesis's simulated table.

library(SuppDists)

src <- "../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر"
tab <- read.table(file.path(src, "table KaH-III.txt"), header = TRUE, check.names = FALSE)
names(tab) <- c("p", "n", "a005", "a01", "a025", "a05")
bad <- tab[!complete.cases(suppressWarnings(sapply(tab, as.numeric))), ]
if (nrow(bad)) { cat("Non-numeric rows dropped:\n"); print(bad) }
tab[] <- lapply(tab, function(v) suppressWarnings(as.numeric(v)))
tab <- tab[complete.cases(tab), ]

tab$nu <- floor(tab$n / 3) - (tab$p + 1)
tab <- tab[tab$nu >= 1, ]

alphas <- c(a005 = .005, a01 = .01, a025 = .025, a05 = .05)
for (a in names(alphas)) {
  tab[[paste0("H_", a)]] <- qmaxFratio(1 - alphas[a], df = tab$nu, k = 3)
  tab[[paste0("rel_", a)]] <- tab[[a]] / tab[[paste0("H_", a)]] - 1
}

cat("Rows compared:", nrow(tab), "\n\n")
print(head(tab[, c("p", "n", "nu", "a05", "H_a05", "a01", "H_a01", "a005", "H_a005")], 12), digits = 4)

cat("\nRelative difference (thesis / exact - 1), summary by alpha:\n")
print(sapply(names(alphas), function(a) {
  r <- tab[[paste0("rel_", a)]]
  c(median = median(r), mean_abs = mean(abs(r)), max_abs = max(abs(r)))
}), digits = 3)

# Implied size of the thesis cut-offs under the exact law
cat("\nActual size of the thesis 5% cut-off under exact Fmax (summary):\n")
size05 <- 1 - pmaxFratio(tab$a05, df = tab$nu, k = 3)
print(summary(size05))
cat("\nActual size of the thesis 1% cut-off:\n")
print(summary(1 - pmaxFratio(tab$a01, df = tab$nu, k = 3)))

write.csv(tab, "kah3_thesis_vs_hartley.csv", row.names = FALSE)
