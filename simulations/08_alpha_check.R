# Step 8: does stronger LTS trimming fix the size under outliers / heavy tails?
# For alpha in {0.9, 0.75, 0.5}: (a) re-estimate nu* for the part sizes used,
# (b) rerun the step-6 size/power scenarios with the formula critical values.
# Every cell saved to disk (resume-safe).

source("fmax_exact.R")
library(parallel)
out_dir <- "alpha_check"; dir.create(out_dir, showWarnings = FALSE)

alphas   <- c(0.9, 0.75, 0.5)
settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen     <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
R_df <- 10000; R_sp <- 4000; gam <- 0.35

gen_x <- function(n, p) { x <- scale(matrix(runif(n * p), n, p)); x[order(x[, 1]), , drop = FALSE] }

## (a) effective d.f. per (alpha, m, p)
df_grid <- unique(data.frame(a = rep(alphas, each = nrow(settings)),
                             m = floor(rep(settings$n, 3) / 3), p = rep(settings$p, 3)))
df_cell <- function(k) {
  suppressMessages(library(robustbase))
  a <- df_grid$a[k]; m <- df_grid$m[k]; p <- df_grid$p[k]
  f <- file.path(out_dir, sprintf("df_a%.2f_m%03d_p%d.rds", a, m, p))
  if (file.exists(f)) return(readRDS(f))
  set.seed(k * 31 + m)
  s2 <- replicate(R_df, {
    x <- gen_x(3 * m, p)[1:m, , drop = FALSE]
    d <- data.frame(y = drop(1 + x %*% rep(1, p) + rnorm(m)), x)
    ltsReg(y ~ ., data = d, alpha = a, mcd = TRUE)$raw.scale^2
  })
  v <- var(log(s2))
  res <- data.frame(a = a, m = m, p = p,
                    nu_star = 2 * uniroot(function(z) trigamma(z) - v, c(1e-3, 1e5))$root)
  saveRDS(res, f); res
}

## (b) size / power per (alpha, setting, scenario)
sp_grid <- expand.grid(a = alphas, s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
sp_cell <- function(g, cvtab) {
  suppressMessages(library(robustbase))
  a <- sp_grid$a[g]; st <- settings[sp_grid$s[g], ]; sc <- sp_grid$sc[g]; n <- st$n; p <- st$p
  f <- file.path(out_dir, sprintf("sp_a%.2f_n%03d_p%d_%s.rds", a, n, p, sc))
  if (file.exists(f)) return(readRDS(f))
  cv <- cvtab$cv[cvtab$a == a & cvtab$m == floor(n / 3) & cvtab$p == p]
  set.seed(77 * g + n)
  rej <- replicate(R_sp, {
    x <- gen_x(n, p)
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    d <- data.frame(y = drop(1 + x %*% rep(1, p) + e), x)
    idx <- list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
    s <- sapply(idx, function(i) ltsReg(y ~ ., data = d[i, ], alpha = a, mcd = TRUE)$raw.scale)
    max(s) / min(s) > cv
  })
  res <- data.frame(alpha = a, n = n, p = p, scenario = sc, rate = mean(rej))
  saveRDS(res, f); res
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("df_grid", "sp_grid", "settings", "R_df", "R_sp", "gam", "out_dir", "gen_x", "df_cell", "sp_cell"))
dfs <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(df_grid)), function(k) df_cell(k)))
dfs$nu <- dfs$m - dfs$p - 1
dfs$ratio <- dfs$nu_star / dfs$nu
dfs$cv <- sqrt(qfmax(rep(.95, nrow(dfs)), dfs$nu_star))
cat("Effective d.f. and 5% critical values:\n"); print(dfs, digits = 3, row.names = FALSE)
clusterExport(cl, "dfs")
sp <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(sp_grid)), function(g) sp_cell(g, dfs)))
stopCluster(cl)

wide <- reshape(sp, idvar = c("n", "p", "scenario"), timevar = "alpha", direction = "wide")
wide <- wide[order(match(wide$scenario, scen), wide$p, wide$n), ]
cat("\nRejection rate at 5% by LTS alpha (MC SE 0.0034):\n"); print(wide, digits = 3, row.names = FALSE)
write.csv(dfs, "alpha_check_df.csv", row.names = FALSE)
write.csv(wide, "alpha_check_size_power.csv", row.names = FALSE)
