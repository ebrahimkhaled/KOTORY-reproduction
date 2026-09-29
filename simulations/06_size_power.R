# Step 6: size and power at the 5% level, table vs distribution decisions,
# with the Koenker (studentized) Breusch-Pagan test as a reference.
# Scenarios: H0 clean normal, H0 + 10% outliers, H0 with t5 errors,
# H1 (sd = exp(g * x1)) clean, H1 + 10% outliers. Each cell saved to disk.

source("fmax_exact.R")
library(parallel)

src <- "../01_extracted/sim_zip/الاختبارات وجداولها للاستخدام المباشر"
rd <- function(f) { t <- read.table(file.path(src, f), header = TRUE, check.names = FALSE)
  names(t) <- c("p", "n", "a005", "a01", "a025", "a05"); t[] <- lapply(t, function(v) suppressWarnings(as.numeric(v))); t }
tab1 <- rd("table KaH-III.txt"); tab2 <- rd("table KaH®V2-III.txt")
dfs <- read.csv("lts_effective_df.csv")
nu_star <- function(m, p) { d <- dfs[dfs$p == p, ]; d <- d[order(d$m), ]
  exp(approx(log(d$m), log(d$nu_star), log(m), rule = 2)$y) }

settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
grid <- expand.grid(s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
for (k in seq_len(nrow(settings))) {             # critical values once per setting
  n <- settings$n[k]; p <- settings$p[k]; m <- floor(n / 3)
  settings$cv1_dis[k] <- qfmax(.95, m - p - 1)
  settings$cv2_dis[k] <- sqrt(qfmax(.95, nu_star(m, p)))
  settings$cv1_tab[k] <- tab1$a05[tab1$p == p & tab1$n == 3 * m]
  settings$cv2_tab[k] <- tab2$a05[tab2$p == p & tab2$n == 3 * m]
}
R <- 4000; gam <- 0.35; out_dir <- "size_power"; dir.create(out_dir, showWarnings = FALSE)

cell <- function(g) {
  suppressMessages(library(robustbase))
  k <- grid$s[g]; sc <- grid$sc[g]; st <- settings[k, ]; n <- st$n; p <- st$p
  f <- file.path(out_dir, sprintf("n%03d_p%d_%s.rds", n, p, sc))
  if (file.exists(f)) return(readRDS(f))
  set.seed(99 * g + n)
  rej <- matrix(0, R, 5, dimnames = list(NULL, c("KaH3_dist", "KaH3_table", "V2_dist", "V2_table", "BP_Koenker")))
  for (r in 1:R) {
    x <- scale(matrix(runif(n * p), n, p)); x <- x[order(x[, 1]), , drop = FALSE]
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    y <- drop(1 + x %*% rep(1, p) + e)
    d <- data.frame(y = y, x)
    idx <- list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
    ms <- sapply(idx, function(i) { fi <- lm.fit(cbind(1, x[i, , drop = FALSE]), y[i]); sum(fi$residuals^2) / (length(i) - p - 1) })
    sc2 <- sapply(idx, function(i) ltsReg(y ~ ., data = d[i, ], alpha = .9, mcd = TRUE)$raw.scale)
    k1 <- max(ms) / min(ms); k2 <- max(sc2) / min(sc2)
    res <- lm.fit(cbind(1, x), y)$residuals^2                      # Koenker studentized BP
    bp <- n * summary(lm(res ~ x))$r.squared
    rej[r, ] <- c(k1 > st$cv1_dis, k1 > st$cv1_tab, k2 > st$cv2_dis, k2 > st$cv2_tab, bp > qchisq(.95, p))
  }
  out <- data.frame(n = n, p = p, scenario = sc, t(colMeans(rej)))
  saveRDS(out, f); out
}

cl <- makeCluster(min(nrow(grid), max(1, detectCores() - 2)))
clusterExport(cl, c("grid", "settings", "R", "gam", "out_dir", "cell"))
res <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(grid)), cell))
stopCluster(cl)
res <- res[order(res$scenario, res$p, res$n), ]
print(res, digits = 3, row.names = FALSE)
cat("MC SE at 5%:", round(sqrt(.05 * .95 / R), 4), "\n")
write.csv(res, "size_power.csv", row.names = FALSE)
