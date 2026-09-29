# Step 26: KaH-robust with alpha = 0.5 (and 0.75 on the same data) in the shape study, because
# centre-clustered outliers (10% of n = 30% of the middle third) exceed the per-part breakdown of
# alpha = 0.75 (Section 2.4). Same data generation as step 24; fingerprints checked.
source("fmax_exact.R")
library(parallel)
out_dir <- "shapes_a50"; dir.create(out_dir, showWarnings = FALSE)
dfa <- read.csv("alpha_check_df.csv")
grid <- expand.grid(shape = c("none", "U", "bulge", "mono"), outl = c("none", "centre", "random"),
                    n = c(30, 60, 150), ch = 1:4, stringsAsFactors = FALSE)
R <- 250
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
sdfun <- function(shape, x) switch(shape, none = rep(1, length(x)), U = 1 + x^2,
                                   bulge = 1 + 3 * exp(-x^2), mono = exp(0.35 * x))
task <- function(k) {
  suppressMessages(library(robustbase))
  sh <- grid$shape[k]; ol <- grid$outl[k]; n <- grid$n[k]; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("%s_%s_n%03d_ch%d.rds", sh, ol, n, grid$ch[k]))
  if (file.exists(f)) return(invisible(NULL))
  nu50 <- dfa$nu_star[dfa$a == .5 & dfa$m == m & dfa$p == 1]
  nu75 <- dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == 1]
  set.seed(9e5 + k)
  P <- t(sapply(seq_len(R), function(r) {
    x <- matrix(sort(scale(runif(n))[, 1]), n, 1)
    mu <- rep(0, n); no <- round(0.1 * n)
    if (ol == "centre") { st <- floor((n - no) / 2); mu[(st + 1):(st + no)] <- 4 }
    if (ol == "random") mu[sample(n, no)] <- 4
    y <- drop(1 + x[, 1] + rnorm(n, mu, sdfun(sh, x[, 1])))
    sc <- function(a) keep_rng(sapply(parts(n), function(i) ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = a, mcd = FALSE)$raw.scale))
    s5 <- sc(0.5); s7 <- sc(0.75)
    c(fingerprint = sum(y * seq_along(y)),
      KaHrobust50 = 1 - pfmax((max(s5) / min(s5))^2, nu50),
      KaHrobust75 = 1 - pfmax((max(s7) / min(s7))^2, nu75))
  }))
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", k)
  saveRDS(data.frame(shape = sh, outl = ol, n = n, P), f)
  invisible(NULL)
}
cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("grid", "R", "out_dir", "dfa", "keep_rng", "parts", "sdfun", "task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(grid)), task))
stopCluster(cl)
all <- do.call(rbind, lapply(list.files(out_dir, "\\.rds$", full.names = TRUE), readRDS))
agg <- aggregate(all[c("KaHrobust50", "KaHrobust75")] < .05, by = all[c("shape", "outl", "n")], FUN = mean)
agg <- agg[order(match(agg$shape, c("none", "U", "bulge", "mono")), match(agg$outl, c("none", "centre", "random")), agg$n), ]
agg[4:5] <- round(100 * agg[4:5], 1)
print(agg, row.names = FALSE)
write.csv(agg, "shapes_alpha50.csv", row.names = FALSE)
