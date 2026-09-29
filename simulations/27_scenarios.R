# Step 27: settings in which White's test is structurally weak (chosen for a stated reason,
# reported with their H0 versions and all rivals):
#   B  several regressors, p = 4, variance driven by x1 (White's auxiliary regression has 14 terms)
#   C  variance driven by an outside ordering variable z (e.g. time), not by the regressors
#   D  5% bad leverage points (extreme x1 with y off the line), which land in the last third
# Each with H0 and H1 (sd = exp(0.35 * sort variable)); B and C also with 10% random N(0, 7^2)
# vertical outliers. n = 45, 90, 150. All tests on the same data; RNG-isolated calls; fingerprints.
source("fmax_exact.R")
library(parallel)
out_dir <- "scenarios"; dir.create(out_dir, showWarnings = FALSE)
load("../KOTORY/R/sysdata.rda")                       # .nu_table (effective d.f.)
nu_star <- function(m, p, a = .75) { t <- .nu_table[.nu_table$alpha == a & .nu_table$p == p, ]; t <- t[order(t$m), ]
  (if (m > max(t$m)) mean(t$ratio[t$m >= 50]) else approx(log(t$m), t$ratio, log(m))$y) * (m - p - 1) }
sc <- expand.grid(set = c("B", "C", "D"), h1 = c(FALSE, TRUE), outl = c(FALSE, TRUE), n = c(45, 90, 150), ch = 1:4,
                  stringsAsFactors = FALSE)
sc <- sc[!(sc$set == "D" & sc$outl), ]                # D is itself a contamination setting
R <- 250
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
white_p <- function(x, e) {
  W <- cbind(x, x^2); if (ncol(x) > 1) for (i in 1:(ncol(x) - 1)) for (j in (i + 1):ncol(x)) W <- cbind(W, x[, i] * x[, j])
  u <- e^2; r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  pchisq(length(u) * r2, ncol(W), lower.tail = FALSE) }

task <- function(k) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic) })
  s <- sc$set[k]; h1 <- sc$h1[k]; ol <- sc$outl[k]; n <- sc$n[k]
  f <- file.path(out_dir, sprintf("%s_%s_%s_n%03d_ch%d.rds", s, if (h1) "H1" else "H0", if (ol) "out" else "clean", n, sc$ch[k]))
  if (file.exists(f)) return(invisible(NULL))
  p <- if (s == "B") 4 else 1; m <- floor(n / 3); nu <- nu_star(m, p)
  set.seed(1.1e6 + k)
  P <- t(sapply(seq_len(R), function(r) {
    x <- scale(matrix(runif(n * p), n, p)); colnames(x) <- paste0("x", 1:p)
    z <- if (s == "C") as.numeric(scale(seq_len(n))) else x[, 1]     # C: time order, independent of x
    o <- order(z); x <- x[o, , drop = FALSE]; z <- z[o]
    e <- rnorm(n, 0, if (h1) exp(0.35 * z) else 1)
    if (ol) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    y <- drop(1 + x %*% rep(1, p) + e)
    if (s == "D") {                                   # 5% bad leverage points at the far right
      i <- sample(n, round(.05 * n)); x[i, 1] <- 3 + runif(length(i), 0, 0.5); y[i] <- 1 - x[i, 1] + rnorm(length(i))
      o <- order(x[, 1]); x <- x[o, , drop = FALSE]; y <- y[o]; z <- x[, 1]
    }
    d <- data.frame(y = y, x); m0 <- lm(y ~ ., d)
    ms <- sapply(parts(n), function(i) sum(lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])$residuals^2) / (length(i) - p - 1))
    sL <- keep_rng(sapply(parts(n), function(i) ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = .75, mcd = FALSE)$raw.scale))
    c(fingerprint = sum(y * seq_along(y)),
      KaHrobust = 1 - pfmax((max(sL) / min(sL))^2, nu),
      KaH3 = 1 - pfmax(max(ms) / min(ms), m - p - 1),
      White = white_p(x, m0$residuals),
      WhiteLTS = keep_rng({ keep <- ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
        white_p(x[keep, , drop = FALSE], lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals) }),
      BP = unname(bptest(m0)$p.value),
      WK = if (p == 1) keep_rng(tryCatch(wilcox_keselman(m0, B = 200L)$p.value, error = function(e) NA)) else NA)
  }))
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", k)
  saveRDS(data.frame(set = s, h1 = h1, outl = ol, n = n, P), f)
  invisible(NULL)
}
cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("sc", "R", "out_dir", ".nu_table", "nu_star", "keep_rng", "parts", "white_p", "task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(sc)), task))
stopCluster(cl)
all <- do.call(rbind, lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS))
stopifnot(!anyDuplicated(paste(all$set, all$h1, all$outl, all$n, round(all$fingerprint, 10))))
meth <- c("KaHrobust", "KaH3", "White", "WhiteLTS", "BP", "WK")
agg <- aggregate(all[meth] < .05, by = all[c("set", "h1", "outl", "n")], FUN = function(v) mean(v, na.rm = TRUE))
agg <- agg[order(agg$set, agg$outl, agg$h1, agg$n), ]
agg[meth] <- round(100 * agg[meth], 1)
print(agg, row.names = FALSE)
write.csv(agg, "scenarios_rejection.csv", row.names = FALSE)
