# Step 4: validate  KaH(R)V2-III^2  ~  Fmax(k = 3, nu*)  on the FULL statistic,
# computed exactly as the thesis code does (sort by x1, 3 parts, ltsReg
# alpha = .9, mcd = TRUE, raw.scale, max/min). Per-point results saved to disk.

source("fmax_exact.R")
library(parallel)

dfs <- read.csv("lts_effective_df.csv")
nu_star <- function(m, p) {                    # interpolate on log scale in m
  d <- dfs[dfs$p == p, ]; d <- d[order(d$m), ]
  exp(approx(log(d$m), log(d$nu_star), log(m), rule = 2)$y)
}

pts <- data.frame(n = c(15, 21, 30, 45, 60, 90, 150, 300, 24, 36, 75, 120),
                  p = c(1, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2))
R <- 20000
out_dir <- "robust_full"; dir.create(out_dir, showWarnings = FALSE)

sim_point <- function(k) {
  suppressMessages(library(robustbase))
  n <- pts$n[k]; p <- pts$p[k]
  f <- file.path(out_dir, sprintf("n%03d_p%d.rds", n, p))
  if (file.exists(f)) return(readRDS(f))
  set.seed(2421974 + 7 * n + p)
  stat <- numeric(R)
  for (r in 1:R) {
    x <- scale(matrix(runif(n * p), n, p))
    x <- x[order(x[, 1]), , drop = FALSE]
    d <- data.frame(y = drop(1 + x %*% rep(1, p) + rnorm(n)), x)
    idx <- list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
    s <- sapply(idx, function(i) ltsReg(y ~ ., data = d[i, ], alpha = .9, mcd = TRUE)$raw.scale)
    stat[r] <- max(s) / min(s)
  }
  saveRDS(stat, f)
  stat
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("pts", "R", "out_dir"))
sims <- parLapplyLB(cl, seq_len(nrow(pts)), sim_point)
stopCluster(cl)

thesis <- read.table("../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر/table KaH®V2-III.txt",
                     header = TRUE, check.names = FALSE)
names(thesis) <- c("p", "n", "a005", "a01", "a025", "a05")
thesis[] <- lapply(thesis, function(v) suppressWarnings(as.numeric(v)))

rows <- lapply(seq_len(nrow(pts)), function(k) {
  n <- pts$n[k]; p <- pts$p[k]; m <- floor(n / 3)
  s <- sims[[k]]; ns <- nu_star(m, p)
  th <- thesis[thesis$p == p & thesis$n == 3 * m, ]
  approx_q <- function(a) sqrt(qfmax(1 - a, ns))
  data.frame(n = n, p = p, nu = m - p - 1, nu_star = round(ns, 2),
             mc_05 = quantile(s, .95), fm_05 = approx_q(.05), th_05 = if (nrow(th)) th$a05 else NA,
             mc_01 = quantile(s, .99), fm_01 = approx_q(.01), th_01 = if (nrow(th)) th$a01 else NA,
             size05_fm = mean(s > approx_q(.05)), size01_fm = mean(s > approx_q(.01)),
             size05_th = if (nrow(th)) mean(s > th$a05) else NA,
             ks_p = suppressWarnings(ks.test(s^2, function(q) pfmax(q, ns))$p.value))
})
res <- do.call(rbind, rows)
print(res, digits = 4, row.names = FALSE)
cat("\nMC standard error of a 5% size with R =", R, ":", round(sqrt(.05 * .95 / R), 4),
    "; of a 1% size:", round(sqrt(.01 * .99 / R), 4), "\n")
write.csv(res, "robust_validation.csv", row.names = FALSE)
