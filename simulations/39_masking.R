# Step 39: concentrated contamination and masking. One regressor, observations sorted by x1,
# variance monotone in x1 under H1 (sd = exp(0.35 x1)), n in {45, 90, 150}, 1000 replications.
#   block "mask": N(0, 7^2) outliers only in the first (low-variance) part, k = 0, 10, 20, 30% of it
#   block "rate": N(0, 7^2) outliers at random positions, 5, 10, 20, 30% of n
#   block "onesided": 10% of N(10, 1) outliers at random positions (asymmetric contamination)
# Every cell is run under H0 and H1 on independent data; H0 gives the size and the thresholds for
# size-adjusted power. KaH statistics are stored as statistics (size adjustment on the statistic),
# White after an LTS screen as a p-value. Foreign calls RNG-isolated; fingerprints must be unique.
library(parallel)
out_dir <- Sys.getenv("MASK_DIR", "masking"); dir.create(out_dir, showWarnings = FALSE)
kotory <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/ElKotory_Thesis_to_Paper/KOTORY"
cells <- rbind(
  expand.grid(block = "mask", level = c(0, 0.1, 0.2, 0.3), stringsAsFactors = FALSE),
  expand.grid(block = "rate", level = c(0.05, 0.1, 0.2, 0.3), stringsAsFactors = FALSE),
  expand.grid(block = "onesided", level = 0.1, stringsAsFactors = FALSE))
cells <- merge(cells, expand.grid(hyp = c("H0", "H1"), n = c(45, 90, 150), ch = 1:4, stringsAsFactors = FALSE))
R <- as.integer(Sys.getenv("MASK_R", "250"))
task <- function(k) {
  suppressMessages({ library(robustbase); pkgload::load_all(kotory, quiet = TRUE, export_all = TRUE) })
  keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
  white_p <- function(x, e) { W <- cbind(x, x^2); u <- e^2
    r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
    pchisq(length(u) * r2, ncol(W), lower.tail = FALSE) }
  cl <- cells[k, ]; n <- cl$n; m <- floor(n / 3); p <- 1
  f <- file.path(out_dir, sprintf("%s_%03.0f_%s_n%03d_ch%d.rds", cl$block, 100 * cl$level, cl$hyp, n, cl$ch))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(4e6 + k)
  P <- t(sapply(seq_len(R), function(r) {
    x <- sort(scale(runif(n))[, 1]); X <- matrix(x, ncol = 1)
    sdv <- if (cl$hyp == "H1") exp(0.35 * x) else rep(1, n)
    e <- rnorm(n, 0, sdv)
    if (cl$block == "mask" && cl$level > 0) {
      i <- sample(m, round(cl$level * m)); e[i] <- rnorm(length(i), 0, 7)      # first (low-variance) part
    } else if (cl$block == "rate") {
      i <- sample(n, round(cl$level * n)); e[i] <- rnorm(length(i), 0, 7)
    } else if (cl$block == "onesided") {
      i <- sample(n, round(cl$level * n)); e[i] <- rnorm(length(i), 10, 1)
    }
    y <- 1 + x + e
    ms <- sapply(list(1:m, (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n),
                 function(j) sum(lm.fit(cbind(1, x[j]), y[j])$residuals^2) / (length(j) - 2))
    s75 <- keep_rng(.lts3(X, y, 0.75)$scale); s50 <- keep_rng(.lts3(X, y, 0.5)$scale)
    wl <- keep_rng(tryCatch({ keep <- ltsReg(x = X, y = y, mcd = FALSE)$lts.wt == 1
      white_p(X[keep, , drop = FALSE], lm.fit(cbind(1, X[keep, , drop = FALSE]), y[keep])$residuals) },
      error = function(e) NA_real_))
    c(fingerprint = sum(y * seq_along(y)), T_KaH3 = max(ms) / min(ms), T_R75 = (max(s75) / min(s75))^2,
      T_R50 = (max(s50) / min(s50))^2, p_WhiteLTS = wl)
  }))
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", k)
  saveRDS(data.frame(cl[c("block", "level", "hyp", "n")], P, row.names = NULL), f)
  invisible(NULL)
}
todo <- if (nzchar(Sys.getenv("MASK_ONLY"))) as.integer(strsplit(Sys.getenv("MASK_ONLY"), ",")[[1]]) else seq_len(nrow(cells))
cl <- makeCluster(min(length(todo), 16))
clusterExport(cl, c("cells", "R", "out_dir", "kotory", "task"))
invisible(parLapplyLB(cl, todo, task))
stopCluster(cl)
all <- do.call(rbind, lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS))
stopifnot(!anyDuplicated(paste(all$block, all$level, all$hyp, all$n, round(all$fingerprint, 10))))
saveRDS(all, "masking_all.rds")
cat("rows:", nrow(all), "\n")
