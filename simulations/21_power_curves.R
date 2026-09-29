# Step 21: power curves. Rejection rate at the nominal 5% level against the strength gamma of
# heteroscedasticity (error sd = exp(gamma * x1)); gamma = 0 gives the size. Clean errors and
# 10% outliers; n = 60 and 150, p = 1. All tests on the same data sets; per-rep p-values saved;
# every call to another package isolated from the RNG; data fingerprints must be unique.
source("fmax_exact.R")
library(parallel)
out_dir <- "power_curves"; dir.create(out_dir, showWarnings = FALSE)
dfa <- read.csv("alpha_check_df.csv")
gammas <- c(0, 0.1, 0.2, 0.3, 0.4, 0.55, 0.7)
grid <- expand.grid(g = gammas, n = c(60, 150), cont = c(FALSE, TRUE), ch = 1:4)
R <- 250                                           # per task; 4 chunks -> 1000 per cell
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
white_p <- function(x, e) { W <- cbind(x, x^2); u <- e^2
  r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  pchisq(length(u) * r2, ncol(W), lower.tail = FALSE) }
mgq <- function(x, y) { n <- length(y); m <- floor(n / 3)
  f <- robustbase::ltsReg(x = x, y = y, mcd = FALSE); clean <- f$lts.wt == 1
  X <- cbind(1, x); Xc <- X[clean, , drop = FALSE]
  r <- drop(y - X %*% lm.fit(Xc, y[clean])$coefficients)
  r[clean] <- r[clean] / (1 - rowSums((Xc %*% solve(crossprod(Xc))) * Xc))
  s <- median(r[(n - m + 1):n]^2) / median(r[1:m]^2); df <- m - 2
  2 * min(pf(s, df, df), pf(s, df, df, lower.tail = FALSE)) }

task <- function(k) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic) })
  g <- grid$g[k]; n <- grid$n[k]; cont <- grid$cont[k]; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("g%.2f_n%03d_%s_ch%d.rds", g, n, if (cont) "out" else "clean", grid$ch[k]))
  if (file.exists(f)) return(invisible(NULL))
  nu <- dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == 1]
  set.seed(5e5 + k)
  P <- t(sapply(seq_len(R), function(r) {
    x <- matrix(sort(scale(runif(n))[, 1]), n, 1, dimnames = list(NULL, "x1"))
    e <- rnorm(n, 0, exp(g * x[, 1]))
    if (cont) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    y <- drop(1 + x[, 1] + e)
    d <- data.frame(y = y, x); m0 <- lm(y ~ x1, d)
    ms <- sapply(parts(n), function(i) sum(lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])$residuals^2) / (length(i) - 2))
    s <- keep_rng(sapply(parts(n), function(i) ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = .75, mcd = FALSE)$raw.scale))
    c(fingerprint = sum(y * seq_along(y)),
      KaH3 = 1 - pfmax(max(ms) / min(ms), m - 2),
      KaHrobust = 1 - pfmax((max(s) / min(s))^2, nu),
      GQ = keep_rng(gqtest(m0, order.by = x[, 1], fraction = n - 2 * m, alternative = "two.sided")$p.value),
      BP = unname(bptest(m0)$p.value),
      White = white_p(x, m0$residuals),
      MGQ = keep_rng(mgq(x, y)),
      WhiteLTS = keep_rng({ keep <- ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
        white_p(x[keep, , drop = FALSE], lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals) }),
      WK = keep_rng(tryCatch(wilcox_keselman(m0, B = 200L)$p.value, error = function(e) NA)))
  }))
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", k)
  saveRDS(data.frame(g = g, n = n, cont = cont, P), f)
  invisible(NULL)
}
cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("grid", "R", "out_dir", "dfa", "keep_rng", "parts", "white_p", "mgq", "task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(grid)), task))
stopCluster(cl)
all <- do.call(rbind, lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS))
names(all)[names(all) == "BP.BP"] <- "BP"            # files written before the unname() fix
stopifnot(!anyDuplicated(paste(all$g, all$n, all$cont, round(all$fingerprint, 10))))
meth <- c("KaH3", "KaHrobust", "GQ", "BP", "White", "MGQ", "WhiteLTS", "WK")
agg <- aggregate(all[meth] < .05, by = all[c("g", "n", "cont")], FUN = function(v) mean(v, na.rm = TRUE))
agg <- agg[order(agg$cont, agg$n, agg$g), ]
print(format(agg, digits = 2), row.names = FALSE)
write.csv(agg, "power_curves.csv", row.names = FALSE)
