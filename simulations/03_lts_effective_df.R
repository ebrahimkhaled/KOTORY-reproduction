# Step 3: effective degrees of freedom of the LTS raw scale (alpha = .9, mcd)
# in ONE sub-group of size m with p regressors + intercept.
# If s^2 ~ sigma^2 chi2_nu*/nu*, then Var(log s^2) = trigamma(nu*/2).
# Invert that to get nu*. Each design point is saved to disk (resume-safe).

library(parallel)
out_dir <- "lts_df"
dir.create(out_dir, showWarnings = FALSE)

R  <- 20000
ms <- c(5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 200)
grid <- expand.grid(m = ms, p = 1:3)
grid <- grid[grid$m - grid$p - 1 >= 2, ]

one_point <- function(k) {
  suppressMessages(library(robustbase))
  m <- grid$m[k]; p <- grid$p[k]
  f <- file.path(out_dir, sprintf("m%03d_p%d.rds", m, p))
  if (file.exists(f)) return(readRDS(f))
  set.seed(1000 * m + p)
  s2 <- ols <- numeric(R)
  for (r in 1:R) {
    x <- scale(matrix(runif(3 * m * p), 3 * m, p))[1:m, , drop = FALSE]  # same design recipe as thesis
    x <- x[order(x[, 1]), , drop = FALSE]
    d <- data.frame(y = drop(1 + x %*% rep(1, p) + rnorm(m)), x)
    fit <- ltsReg(y ~ ., data = d, alpha = .9, mcd = TRUE)
    s2[r] <- fit$raw.scale^2
    ols[r] <- sum(lm.fit(cbind(1, x), d$y)$residuals^2) / (m - p - 1)
  }
  inv_trigamma <- function(v) uniroot(function(a) trigamma(a) - v, c(1e-3, 1e5))$root
  ok <- s2 > 0
  res <- list(m = m, p = p, nu = m - p - 1,
              zero_frac = mean(!ok),
              nu_star = 2 * inv_trigamma(var(log(s2[ok]))),
              nu_ols_check = 2 * inv_trigamma(var(log(ols))),   # should recover nu
              s2 = s2)
  saveRDS(res, f)
  res
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("grid", "R", "out_dir", "one_point"))
t0 <- Sys.time()
# LTS is undefined when m <= 2(p+1) (same NA cells as the thesis table) -> skip
res <- parLapplyLB(cl, seq_len(nrow(grid)), function(k) tryCatch(one_point(k), error = function(e) NULL))
stopCluster(cl)
res <- Filter(Negate(is.null), res)
cat("elapsed:", format(Sys.time() - t0), "\n")

tab <- do.call(rbind, lapply(res, function(r)
  data.frame(m = r$m, p = r$p, nu = r$nu, nu_ols_check = r$nu_ols_check,
             nu_star = r$nu_star, ratio = r$nu_star / r$nu, zero_frac = r$zero_frac)))
print(tab[order(tab$p, tab$m), ], digits = 4, row.names = FALSE)
write.csv(tab, "lts_effective_df.csv", row.names = FALSE)
