# Step 33: heavy-tailed errors. The factorial design of step 28 (shape x order x p x n, no planted
# outliers) with the normal errors replaced by t3 errors scaled to unit variance, so the tails
# themselves act as scattered outliers. Same tests, same p-value code.
# KaH-robust is also run with its residual bootstrap reference (B = 199) in the null cells, to
# measure how far the normal-theory Fmax reference drifts and whether the bootstrap corrects it.
# Per-replication p-values saved; foreign calls RNG-isolated; data fingerprints must be unique.
source("fmax_exact.R")
library(parallel)
out_dir <- Sys.getenv("HT_DIR", "heavy_tails"); dir.create(out_dir, showWarnings = FALSE)
load("../KOTORY/R/sysdata.rda")
nu_star <- function(m, p, a = .75) { t <- .nu_table[.nu_table$alpha == a & .nu_table$p == p, ]; t <- t[order(t$m), ]
  (if (m > max(t$m)) mean(t$ratio[t$m >= 50]) else approx(log(t$m), t$ratio, log(m))$y) * (m - p - 1) }
cells <- expand.grid(shape = c("none", "mono", "U", "bulge"), order = c("x1", "time"), p = c(1, 4),
                     n = c(45, 90, 150), ch = 1:4, stringsAsFactors = FALSE)
R <- as.integer(Sys.getenv("HT_R", "250"))
B <- 199
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
white_p <- function(x, e) {
  W <- cbind(x, x^2); if (ncol(x) > 1) for (i in 1:(ncol(x) - 1)) for (j in (i + 1):ncol(x)) W <- cbind(W, x[, i] * x[, j])
  u <- e^2; r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  pchisq(length(u) * r2, ncol(W), lower.tail = FALSE) }
mgq <- function(x, y) { n <- length(y); m <- floor(n / 3); p <- ncol(x)
  f <- robustbase::ltsReg(x = x, y = y, mcd = FALSE); clean <- f$lts.wt == 1
  X <- cbind(1, x); Xc <- X[clean, , drop = FALSE]
  r <- drop(y - X %*% lm.fit(Xc, y[clean])$coefficients)
  r[clean] <- r[clean] / (1 - rowSums((Xc %*% solve(crossprod(Xc))) * Xc))
  s <- median(r[(n - m + 1):n]^2) / median(r[1:m]^2); df <- m - p - 1
  2 * min(pf(s, df, df), pf(s, df, df, lower.tail = FALSE)) }
bamset <- function(m0) { e <- skedastic::blus(m0, omit = "last", keepNA = FALSE)
  bartlett.test(e, cut(seq_along(e), 3, labels = FALSE))$p.value }
sdfun <- function(shape, z) switch(shape, none = rep(1, length(z)), mono = exp(0.35 * z), U = 1 + z^2,
                                   bulge = 1 + 3 * exp(-z^2))
task <- function(k) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic)
    pkgload::load_all("../KOTORY", quiet = TRUE, export_all = TRUE) })
  cl <- cells[k, ]; n <- cl$n; p <- cl$p; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("%s_%s_p%d_t3_n%03d_ch%d.rds", cl$shape, cl$order, p, n, cl$ch))
  if (file.exists(f)) return(invisible(NULL))
  nu <- nu_star(m, p)
  set.seed(3e6 + k)
  P <- t(sapply(seq_len(R), function(r) {
    x <- scale(matrix(runif(n * p), n, p)); colnames(x) <- paste0("x", 1:p)
    z <- if (cl$order == "time") as.numeric(scale(seq_len(n))) else x[, 1]
    o <- order(z); x <- x[o, , drop = FALSE]; z <- z[o]
    e <- sdfun(cl$shape, z) * rt(n, 3) / sqrt(3)
    y <- drop(1 + x %*% rep(1, p) + e)
    d <- data.frame(y = y, x); m0 <- lm(y ~ ., d)
    ms <- sapply(parts(n), function(i) sum(lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])$residuals^2) / (length(i) - p - 1))
    sL <- keep_rng(sapply(parts(n), function(i) ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = .75, mcd = FALSE)$raw.scale))
    c(fingerprint = sum(y * seq_along(y)),
      KaHrobust = 1 - pfmax((max(sL) / min(sL))^2, nu),
      KaHrobustBoot = if (cl$shape == "none") keep_rng(tryCatch(
        KOTORY:::.kah_robust_core(list(y = y, Z = x, ord = z, n = n, p = p), .75, "bootstrap", B, r)$p.value,
        error = function(e) NA)) else NA,
      KaH3 = 1 - pfmax(max(ms) / min(ms), m - p - 1),
      GQ = keep_rng(gqtest(m0, order.by = z, fraction = n - 2 * m, alternative = "two.sided")$p.value),
      BAMSET = keep_rng(bamset(m0)),
      BP = unname(bptest(m0)$p.value),
      White = white_p(x, m0$residuals),
      WhiteLTS = keep_rng({ keep <- ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
        white_p(x[keep, , drop = FALSE], lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals) }),
      MGQ = keep_rng(mgq(x, y)),
      EvansKing = keep_rng(tryCatch(evans_king(m0, deflator = "x1", method = "GLS")$p.value, error = function(e) NA)),
      WK = if (p == 1) keep_rng(tryCatch(wilcox_keselman(m0, B = 200L)$p.value, error = function(e) NA)) else NA)
  }))
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", k)
  saveRDS(data.frame(shape = cl$shape, order = cl$order, p = p, n = n, P), f)
  invisible(NULL)
}
todo <- if (nzchar(Sys.getenv("HT_ONLY"))) as.integer(strsplit(Sys.getenv("HT_ONLY"), ",")[[1]]) else seq_len(nrow(cells))
cl <- makeCluster(max(1, min(length(todo), detectCores() - 2)))
clusterExport(cl, c("cells", "R", "B", "out_dir", ".nu_table", "nu_star", "keep_rng", "parts", "white_p", "mgq",
                    "bamset", "sdfun", "task", "pfmax"))
invisible(parLapplyLB(cl, todo, task))
stopCluster(cl)
all <- do.call(rbind, lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS))
stopifnot(!anyDuplicated(paste(all$shape, all$order, all$p, all$n, round(all$fingerprint, 10))))
saveRDS(all, "heavy_tails_all.rds")
cat("rows:", nrow(all), "\n")
