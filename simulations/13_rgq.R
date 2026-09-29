# Step 13: Alih & Ong (2015) RGQ on the step-11 settings and scenarios, with KaH3 and
# KaH-robust (alpha .75, Fmax) re-run on the SAME data for a direct comparison.
#   RGQ_D      paper as written (rank by robust distance D, literal Step 4(b)), two-sided p
#   RGQ_D_1s   same statistic, one-sided upper-tail p (what the paper's examples appear to use)
#   RGQ_x      charitable variant: rank by x1 (classical GQ ordering), two-sided p
# The paper's worked example (Table 2) could NOT be reproduced from its description
# (see rgq_alih_ong.R); this is a good-faith implementation, reported as such.
# Guards as in step 11: keep_rng() around every foreign call; data fingerprints must be unique.

source("fmax_exact.R")
source("rgq_alih_ong.R")
library(parallel)
out_dir <- "rgq"; dir.create(out_dir, showWarnings = FALSE)

settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
R <- 1000; chunk <- 250; gam <- 0.35
dfa <- read.csv("alpha_check_df.csv")

keep_rng <- function(expr) {
  old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv()))
  expr
}
gen_x <- function(n, p) { x <- scale(matrix(runif(n * p), n, p)); x[order(x[, 1]), , drop = FALSE] }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)

one_rep <- function(x, y, st, nu75) {
  n <- st$n; p <- st$p
  pv <- c()
  ms <- sapply(parts(n), function(i) sum(lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])$residuals^2) / (length(i) - p - 1))
  pv["KaH3"] <- 1 - pfmax(max(ms) / min(ms), floor(n / 3) - p - 1)
  s <- keep_rng(sapply(parts(n), function(i) robustbase::ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = .75, mcd = FALSE)$raw.scale))
  pv["V2_a75"] <- 1 - pfmax((max(s) / min(s))^2, nu75)
  a <- tryCatch(rgq_test(x, y, "D"), error = function(e) NULL)
  pv["RGQ_D"] <- if (is.null(a)) NA else a$p.value
  pv["RGQ_D_1s"] <- if (is.null(a)) NA else pf(a$statistic, a$df, a$df, lower.tail = FALSE)
  b <- tryCatch(rgq_test(x, y, "x"), error = function(e) NULL)
  pv["RGQ_x"] <- if (is.null(b)) NA else b$p.value
  pv["RGQ_n_out"] <- if (is.null(a)) NA else a$n.outliers
  pv
}

tasks <- expand.grid(ch = seq_len(R / chunk), s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
run_task <- function(t) {
  suppressMessages(library(robustbase))
  st <- settings[tasks$s[t], ]; sc <- tasks$sc[t]; n <- st$n; p <- st$p; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("n%03d_p%d_%s_ch%d.rds", n, p, sc, tasks$ch[t]))
  if (file.exists(f)) return(invisible(NULL))
  nu75 <- dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == p]
  set.seed(9e4 + 1e3 * t + n)
  out <- vector("list", chunk)
  for (r in seq_len(chunk)) {
    x <- gen_x(n, p)
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    y <- drop(1 + x %*% rep(1, p) + e)
    out[[r]] <- c(fingerprint = sum(y * seq_along(y)), one_rep(x, y, st, nu75))
  }
  P <- do.call(rbind, out)
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", t)
  saveRDS(data.frame(n = n, p = p, scenario = sc, P, check.names = FALSE), f)
  invisible(NULL)
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("settings", "tasks", "R", "chunk", "gam", "dfa", "out_dir", "keep_rng", "gen_x", "parts",
                    "one_rep", "run_task", "pfmax", "rgq_outliers", "rgq_test"))
invisible(parLapplyLB(cl, seq_len(nrow(tasks)), run_task))
stopCluster(cl)

all <- do.call(rbind, lapply(list.files(out_dir, "_ch\\d+\\.rds$", full.names = TRUE), readRDS))
stopifnot(!anyDuplicated(paste(all$n, all$scenario, round(all$fingerprint, 10))))
meth <- c("KaH3", "V2_a75", "RGQ_D", "RGQ_D_1s", "RGQ_x")
agg <- aggregate(all[meth] < .05, by = all[c("n", "p", "scenario")], FUN = function(v) mean(v, na.rm = TRUE))
agg$RGQ_mean_outliers <- aggregate(all["RGQ_n_out"], by = all[c("n", "p", "scenario")], FUN = mean, na.rm = TRUE)$RGQ_n_out
agg$RGQ_failed <- aggregate(is.na(all["RGQ_D"]), by = all[c("n", "p", "scenario")], FUN = mean)$RGQ_D
agg <- agg[order(match(agg$scenario, scen), agg$p, agg$n), ]
print(agg, digits = 3, row.names = FALSE)
write.csv(agg, "rgq_rejection.csv", row.names = FALSE)
