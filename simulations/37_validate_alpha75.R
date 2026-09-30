# Step 37: calibration of the Fmax(3, nu*) reference at the package default alpha = 0.75, using the
# package's own LTS fits (KOTORY:::.lts3) and its shipped nu* (kah.nu.star). Same twelve settings as
# step 4 (which used alpha = 0.9 and an earlier nu* table); normal errors, H0, observations sorted by x1.
# Every setting's statistics are saved to disk.
library(parallel)
pts <- data.frame(n = c(15, 21, 30, 45, 60, 90, 150, 300, 24, 36, 75, 120),
                  p = c(1, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2))
R <- as.integer(Sys.getenv("VAL_R", "10000"))
out_dir <- "robust_a75"; dir.create(out_dir, showWarnings = FALSE)
kotory <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/ElKotory_Thesis_to_Paper/KOTORY"

sim_point <- function(k) {
  suppressMessages(pkgload::load_all(kotory, quiet = TRUE, export_all = TRUE))
  n <- pts$n[k]; p <- pts$p[k]
  f <- file.path(out_dir, sprintf("n%03d_p%d.rds", n, p))
  if (file.exists(f)) return(readRDS(f))
  set.seed(7519 + 7 * n + p)
  stat <- vapply(seq_len(R), function(r) {
    x <- scale(matrix(runif(n * p), n, p)); x <- x[order(x[, 1]), , drop = FALSE]
    y <- drop(1 + x %*% rep(1, p) + rnorm(n))
    s <- tryCatch(.lts3(x, y, 0.75)$scale, error = function(e) rep(NA_real_, 3))
    max(s) / min(s)
  }, numeric(1))
  saveRDS(stat, f)
  stat
}
cl <- makeCluster(as.integer(Sys.getenv("VAL_CORES", "8")))
clusterExport(cl, c("pts", "R", "out_dir", "kotory"))
sims <- parLapplyLB(cl, seq_len(nrow(pts)), sim_point)
stopCluster(cl)

suppressMessages(pkgload::load_all(kotory, quiet = TRUE))
res <- do.call(rbind, lapply(seq_len(nrow(pts)), function(k) {
  n <- pts$n[k]; p <- pts$p[k]; m <- floor(n / 3); s <- sims[[k]]; s <- s[is.finite(s)]
  ns <- kah.nu.star(m, p, 0.75)
  pv <- pfmax(s^2, ns, 3, lower.tail = FALSE)
  data.frame(n = n, p = p, m = m, nu = m - p - 1, nu_star = round(ns, 2), reps = length(s),
             size05 = mean(pv < .05), size01 = mean(pv < .01),
             ks_p = suppressWarnings(ks.test(s^2, function(q) pfmax(q, ns, 3))$p.value))
}))
print(res, digits = 3, row.names = FALSE)
cat("\nMC standard error with R =", R, ": 5% size", round(sqrt(.05 * .95 / R), 4), "; 1% size", round(sqrt(.01 * .99 / R), 4), "\n")
write.csv(res, "robust_validation_a75.csv", row.names = FALSE)
