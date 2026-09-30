# Step 11b: step 11 rerun with identical seeds and code, adding BAMSET computed from its
# definition (Bartlett's test on three ordered groups of BLUS residuals, as in step 18b), so that
# BAMSET is evaluated on the same data sets as every other test. Every other column must
# reproduce step 11 exactly; 11b_compare below checks this.
# Step 11: ALL tests on the SAME data (replaces steps 9 and 10, whose results were
# invalid: skedastic::bamset() and skedastic::zhou_etal() reset the global RNG to a
# fixed seed, so every rep after the first re-used the same data set; found 2026-09-28,
# old outputs kept in _invalid_rng_reset/).
# Guards: every call to another package runs inside keep_rng(), which restores the RNG
# state afterwards; each data set's fingerprint is stored and duplicates stop the run.
#
# KaH versions:  KaH3, V2_a90, V2_a75, V2_a75_kc (kurtosis-corrected nu*), V2_a75_bs (bootstrap)
# Rivals:        GQ, BP_Koenker, White, MGQ_Rana2008, WK_2006, Zhou_2015, LiYao_2019,
#                BAMSET_3 (Ramsey 1969), EvansKing (1988), HMC (Harrison-McCabe 1979, lmtest),
#                BRW_White (Berenguer-Rico & Wilms 2021), BF_3parts (Brown-Forsythe, our parts)

source("fmax_exact.R")
library(parallel)
out_dir <- "all_tests_v2"; dir.create(out_dir, showWarnings = FALSE)

settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
R <- 1000; chunk <- 50; B <- 199; gam <- 0.35
dfa <- read.csv("alpha_check_df.csv")
k0_tab <- readRDS("competitors/k0.rds")       # computed with LTS only (valid)

keep_rng <- function(expr) {                   # evaluate, then restore the global RNG state
  old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv()))
  expr
}
gen_x <- function(n, p) { x <- scale(matrix(runif(n * p), n, p)); x[order(x[, 1]), , drop = FALSE] }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
kurt <- function(r) mean(r^4) / mean(r^2)^2

lts3 <- function(x, y, a) {
  out <- lapply(parts(length(y)), function(i) {
    f <- robustbase::ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = a, mcd = FALSE)
    r <- drop(y[i] - cbind(1, x[i, , drop = FALSE]) %*% f$raw.coefficients) / f$raw.scale
    list(s = f$raw.scale, r_all = r, r_best = r[f$best])
  })
  list(s = sapply(out, `[[`, "s"), r_all = unlist(lapply(out, `[[`, "r_all")),
       r_best = unlist(lapply(out, `[[`, "r_best")))
}
white_p <- function(x, e) {
  W <- if (ncol(x) == 1) cbind(x, x^2) else cbind(x, x^2, x[, 1] * x[, 2])
  u <- e^2; r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  pchisq(length(u) * r2, ncol(W), lower.tail = FALSE)
}
mgq <- function(x, y, p) {                     # Rana, Midi & Imon (2008)
  n <- length(y); m <- floor(n / 3)
  f <- robustbase::ltsReg(x = x, y = y, mcd = FALSE); clean <- f$lts.wt == 1
  X <- cbind(1, x); Xc <- X[clean, , drop = FALSE]
  r <- drop(y - X %*% lm.fit(Xc, y[clean])$coefficients)
  r[clean] <- r[clean] / (1 - rowSums((Xc %*% solve(crossprod(Xc))) * Xc))
  s <- median(r[(n - m + 1):n]^2) / median(r[1:m]^2); df <- m - p - 1
  2 * min(pf(s, df, df), pf(s, df, df, lower.tail = FALSE))
}

one_rep <- function(x, y, st, nu) {
  n <- st$n; p <- st$p
  colnames(x) <- paste0("x", 1:p)
  d <- data.frame(y = y, x); m0 <- lm(y ~ ., data = d)
  safe <- function(expr) tryCatch(keep_rng(expr), error = function(e) NA_real_)
  pv <- c()
  fits <- lapply(parts(n), function(i) lm.fit(cbind(1, x[i, , drop = FALSE]), y[i]))
  ms <- sapply(seq_along(fits), function(j) sum(fits[[j]]$residuals^2) / (length(parts(n)[[j]]) - p - 1))
  pv["KaH3"] <- 1 - pfmax(max(ms) / min(ms), floor(n / 3) - p - 1)
  s9 <- lts3(x, y, .9)$s
  pv["V2_a90"] <- 1 - pfmax((max(s9) / min(s9))^2, nu$a90)
  L <- lts3(x, y, .75); T75 <- max(L$s) / min(L$s)
  pv["V2_a75"] <- 1 - pfmax(T75^2, nu$a75)
  f_inf <- (kurt(L$r_best) - 1) / (nu$K0 - 1)
  pv["V2_a75_kc"] <- 1 - pfmax(T75^2, nu$a75 / max(f_inf, .25))
  r <- L$r_all
  Tb <- replicate(B, { s <- lts3(x, sample(r, n, replace = TRUE), .75)$s; max(s) / min(s) })
  pv["V2_a75_bs"] <- (1 + sum(Tb >= T75)) / (B + 1)
  # rivals
  pv["GQ"] <- safe(lmtest::gqtest(m0, order.by = x[, 1], fraction = n - 2 * floor(n / 3), alternative = "two.sided")$p.value)
  pv["BP_Koenker"] <- safe(lmtest::bptest(m0)$p.value)
  pv["White"] <- white_p(x, m0$residuals)
  pv["MGQ_Rana2008"] <- safe(mgq(x, y, p))
  pv["WK_2006"] <- safe(skedastic::wilcox_keselman(m0, B = 200L)$p.value)
  pv["Zhou_2015"] <- safe(skedastic::zhou_etal(m0, method = "pooled", Bperturbed = 200L)$p.value)
  pv["LiYao_2019"] <- safe(skedastic::li_yao(m0, method = "cvt")$p.value)
  pv["BAMSET_3"] <- safe(skedastic::bamset(m0, k = 3, deflator = "x1")$p.value)
  pv["EvansKing"] <- safe(skedastic::evans_king(m0, deflator = "x1", method = "GLS")$p.value)
  pv["HMC"] <- safe(lmtest::hmctest(m0, order.by = ~ x1, data = d, point = 0.5, sim = 1000)$p.value)
  pv["BRW_White"] <- safe({
    keep <- robustbase::ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
    white_p(x[keep, , drop = FALSE], lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals)
  })
  pv["BAMSET_own"] <- safe({ e <- skedastic::blus(m0, omit = "last", keepNA = FALSE)
    bartlett.test(e, cut(seq_along(e), 3, labels = FALSE))$p.value })
  res <- unlist(lapply(fits, function(f) f$residuals))
  g <- factor(rep(1:3, sapply(parts(n), length)))
  pv["BF_3parts"] <- anova(lm(abs(res - ave(res, g, FUN = median)) ~ g))[1, "Pr(>F)"]
  pv
}

tasks <- expand.grid(ch = seq_len(R / chunk), s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
run_task <- function(t) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic) })
  st <- settings[tasks$s[t], ]; sc <- tasks$sc[t]; n <- st$n; p <- st$p; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("n%03d_p%d_%s_ch%02d.rds", n, p, sc, tasks$ch[t]))
  if (file.exists(f)) return(invisible(NULL))
  nu <- list(a90 = dfa$nu_star[dfa$a == .9 & dfa$m == m & dfa$p == p],
             a75 = dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == p],
             K0 = k0_tab$K0[k0_tab$n == n & k0_tab$p == p])
  set.seed(7e4 + 1e3 * t + n)
  out <- vector("list", chunk)
  for (r in seq_len(chunk)) {
    x <- gen_x(n, p)
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    y <- drop(1 + x %*% rep(1, p) + e)
    out[[r]] <- c(fingerprint = sum(y * seq_along(y)), one_rep(x, y, st, nu))
  }
  P <- do.call(rbind, out)
  if (anyDuplicated(round(P[, "fingerprint"], 10))) stop("duplicated data sets in task ", t)
  saveRDS(data.frame(n = n, p = p, scenario = sc, P, check.names = FALSE), f)
  invisible(NULL)
}

if (nzchar(Sys.getenv("SMOKE"))) {             # 3 reps per scenario, timed, duplicate check active
  chunk <- 3; B <- 19
  for (sc in scen) { tasks <- expand.grid(ch = 1, s = 3, sc = sc, stringsAsFactors = FALSE)
    out_dir <- tempdir(); unlink(file.path(out_dir, "*.rds"))
    tm <- system.time(run_task(1))[3]
    a <- readRDS(list.files(out_dir, "\\.rds$", full.names = TRUE)[1])
    cat(sc, sprintf("%.1fs for 3 reps\n", tm)); print(signif(a[, -(1:3)], 3)) }
  quit(save = "no")
}

cl <- makeCluster(min(12, max(1, detectCores() - 2)))
clusterExport(cl, c("settings", "tasks", "R", "chunk", "B", "gam", "dfa", "k0_tab", "out_dir",
                    "keep_rng", "gen_x", "parts", "kurt", "lts3", "white_p", "mgq", "one_rep", "run_task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(tasks)), run_task))
stopCluster(cl)

all <- do.call(rbind, lapply(list.files(out_dir, "_ch\\d+\\.rds$", full.names = TRUE), readRDS))
stopifnot(!anyDuplicated(paste(all$n, all$scenario, round(all$fingerprint, 10))))
meth <- setdiff(names(all), c("n", "p", "scenario", "fingerprint"))
agg <- aggregate(all[meth] < .05, by = all[c("n", "p", "scenario")], FUN = function(v) mean(v, na.rm = TRUE))
agg <- agg[order(match(agg$scenario, scen), agg$p, agg$n), ]
print(agg, digits = 3, row.names = FALSE)
cat("NA counts:\n"); print(colSums(is.na(all[meth])))
cat("MC SE at 5% with R =", R, ":", round(sqrt(.05 * .95 / R), 4), "\n")
write.csv(agg, "all_tests_v2_rejection.csv", row.names = FALSE)
