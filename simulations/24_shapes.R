# Step 24: variance SHAPES, following the thesis design (thesis power code, folder
# "3...", files "test the power of the new test (p=1,s=1,pos=C/R)"): error sd = 1 + x^2 with x
# standardized (U-shape: high-low-high), vertical outliers of mean shift 4 placed in the CENTRE
# (the low-variance region, which masks the U-shape) or at random positions.
#   shapes : none (H0) | U (sd = 1 + x^2) | bulge (sd = 1 + 3 exp(-x^2), low-high-low) | mono (sd = exp(0.35 x))
#   outl   : none | centre 10% | random 10%       n : 30, 60, 150 (p = 1)
# Every call to another package isolated from the RNG; data fingerprints must be unique.
source("fmax_exact.R")
library(parallel)
out_dir <- "shapes"; dir.create(out_dir, showWarnings = FALSE)
dfa <- read.csv("alpha_check_df.csv")
grid <- expand.grid(shape = c("none", "U", "bulge", "mono"), outl = c("none", "centre", "random"),
                    n = c(30, 60, 150), ch = 1:4, stringsAsFactors = FALSE)
R <- 250                                          # per task; 4 chunks -> 1000 per cell
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
bamset <- function(m) { e <- skedastic::blus(m, omit = "last", keepNA = FALSE)
  bartlett.test(e, cut(seq_along(e), 3, labels = FALSE))$p.value }
sdfun <- function(shape, x) switch(shape, none = rep(1, length(x)), U = 1 + x^2,
                                   bulge = 1 + 3 * exp(-x^2), mono = exp(0.35 * x))
# (bulge was first 4 - x^2, which turns negative because a standardized uniform sample of
#  size 30 can reach |x| = 2.7; 1 + 3 exp(-x^2) is 4 at the centre, about 1.15 at |x| = 1.73)

task <- function(k) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic) })
  sh <- grid$shape[k]; ol <- grid$outl[k]; n <- grid$n[k]; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("%s_%s_n%03d_ch%d.rds", sh, ol, n, grid$ch[k]))
  if (file.exists(f)) return(invisible(NULL))
  nu <- dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == 1]
  set.seed(7e5 + k)
  P <- t(sapply(seq_len(R), function(r) {
    x <- matrix(sort(scale(runif(n))[, 1]), n, 1, dimnames = list(NULL, "x1"))
    mu <- rep(0, n); no <- round(0.1 * n)
    if (ol == "centre") { st <- floor((n - no) / 2); mu[(st + 1):(st + no)] <- 4 }   # thesis pos = C
    if (ol == "random") mu[sample(n, no)] <- 4                                       # thesis pos = R
    y <- drop(1 + x[, 1] + rnorm(n, mu, sdfun(sh, x[, 1])))
    d <- data.frame(y = y, x); m0 <- lm(y ~ x1, d)
    ms <- sapply(parts(n), function(i) sum(lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])$residuals^2) / (length(i) - 2))
    s <- keep_rng(sapply(parts(n), function(i) ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = .75, mcd = FALSE)$raw.scale))
    c(fingerprint = sum(y * seq_along(y)),
      KaH3 = 1 - pfmax(max(ms) / min(ms), m - 2),
      KaHrobust = 1 - pfmax((max(s) / min(s))^2, nu),
      GQ = keep_rng(gqtest(m0, order.by = x[, 1], fraction = n - 2 * m, alternative = "two.sided")$p.value),
      BAMSET = keep_rng(bamset(m0)),
      BP = unname(bptest(m0)$p.value),
      White = white_p(x, m0$residuals),
      EvansKing = keep_rng(tryCatch(evans_king(m0, deflator = "x1", method = "GLS")$p.value, error = function(e) NA)),
      MGQ = keep_rng(mgq(x, y)),
      WhiteLTS = keep_rng({ keep <- ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
        white_p(x[keep, , drop = FALSE], lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals) }),
      WK = keep_rng(tryCatch(wilcox_keselman(m0, B = 200L)$p.value, error = function(e) NA)))
  }))
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", k)
  saveRDS(data.frame(shape = sh, outl = ol, n = n, P), f)
  invisible(NULL)
}
cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("grid", "R", "out_dir", "dfa", "keep_rng", "parts", "white_p", "mgq", "bamset",
                    "sdfun", "task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(grid)), task))
stopCluster(cl)
all <- do.call(rbind, lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS))
stopifnot(!anyDuplicated(paste(all$shape, all$outl, all$n, round(all$fingerprint, 10))))
meth <- c("KaH3", "KaHrobust", "GQ", "BAMSET", "BP", "White", "EvansKing", "MGQ", "WhiteLTS", "WK")
agg <- aggregate(all[meth] < .05, by = all[c("shape", "outl", "n")], FUN = function(v) mean(v, na.rm = TRUE))
agg <- agg[order(match(agg$shape, c("none", "U", "bulge", "mono")), match(agg$outl, c("none", "centre", "random")), agg$n), ]
print(format(agg, digits = 2), row.names = FALSE)
write.csv(agg, "shapes_rejection.csv", row.names = FALSE)
