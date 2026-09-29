# Quick look at whatever full-statistic points are already on disk
source("fmax_exact.R")
dfs <- read.csv("lts_effective_df.csv")
nu_star <- function(m, p) {
  d <- dfs[dfs$p == p, ]; d <- d[order(d$m), ]
  exp(approx(log(d$m), log(d$nu_star), log(m), rule = 2)$y)
}
thesis <- read.table("../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر/table KaH®V2-III.txt",
                     header = TRUE, check.names = FALSE)
names(thesis) <- c("p", "n", "a005", "a01", "a025", "a05")
thesis[] <- lapply(thesis, function(v) suppressWarnings(as.numeric(v)))

rows <- lapply(list.files("robust_full", full.names = TRUE), function(f) {
  n <- as.integer(sub(".*n(\\d+)_p.*", "\\1", f)); p <- as.integer(sub(".*_p(\\d)\\.rds", "\\1", f))
  s <- readRDS(f); m <- floor(n / 3); ns <- nu_star(m, p)
  th <- thesis[thesis$p == p & thesis$n == 3 * m, ]
  q <- function(a) sqrt(qfmax(1 - a, ns))
  data.frame(n, p, nu = m - p - 1, nu_star = round(ns, 2),
             mc_q95 = quantile(s, .95), fmax_q95 = q(.05), thesis_q95 = th$a05,
             mc_q99 = quantile(s, .99), fmax_q99 = q(.01), thesis_q99 = th$a01,
             size05 = mean(s > q(.05)), size01 = mean(s > q(.01)),
             ks_p = suppressWarnings(ks.test(s^2, function(x) pfmax(x, ns))$p.value))
})
print(do.call(rbind, rows), digits = 4, row.names = FALSE)
cat("MC SE of a 5% size:", round(sqrt(.05 * .95 / 2e4), 4), " of a 1% size:", round(sqrt(.01 * .99 / 2e4), 4), "\n")
