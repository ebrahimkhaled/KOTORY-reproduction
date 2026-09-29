# Step 9: size-fix trials for the robust KaH test + comparison with competitors.
# Every rep stores all p-values (per-task .rds, resume-safe).
#
# KaH versions:  KaH3        OLS, exact Fmax(3, v)
#                V2_a90      LTS alpha .90, Fmax(3, v*)            (thesis test)
#                V2_a75      LTS alpha .75, Fmax(3, v*)
#                V2_a75_kc   LTS alpha .75, v* / f, f = kurtosis inflation of retained residuals
#                V2_a75_bs   LTS alpha .75, residual-bootstrap p-value (B = 199)
# Competitors:   GQ (1965, middle third dropped, two-sided), BP_Koenker (1981), White (1980),
#                MGQ Rana-Midi-Imon (2008), Wilcox-Keselman (2006), Zhou-Song-Thompson (2015),
#                Li-Yao (2019, cvt)

source("fmax_exact.R")
library(parallel)
out_dir <- "competitors"; dir.create(out_dir, showWarnings = FALSE)

settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
R <- 1000; chunk <- 125; B <- 199; gam <- 0.35
dfa <- read.csv("alpha_check_df.csv")        # nu* for alpha .9 / .75 at these part sizes

gen_x <- function(n, p) { x <- scale(matrix(runif(n * p), n, p)); x[order(x[, 1]), , drop = FALSE] }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)

# LTS scales of the three parts, plus pooled standardized raw residuals of the h-subsets
lts3 <- function(x, y, a, keep = FALSE) {
  out <- lapply(parts(length(y)), function(i) {
    d <- data.frame(y = y[i], x[i, , drop = FALSE])
    f <- robustbase::ltsReg(y ~ ., data = d, alpha = a, mcd = TRUE)
    r <- drop(y[i] - cbind(1, x[i, , drop = FALSE]) %*% f$raw.coefficients) / f$raw.scale
    list(s = f$raw.scale, r_all = r, r_best = r[f$best])
  })
  list(s = sapply(out, `[[`, "s"),
       r_all = unlist(lapply(out, `[[`, "r_all")),
       r_best = unlist(lapply(out, `[[`, "r_best")))
}
kurt <- function(r) mean(r^4) / mean(r^2)^2

## K0: expected kurtosis of the retained standardized residuals under normal errors (alpha .75)
k0_tab <- NULL
k0_file <- file.path(out_dir, "k0.rds")
if (file.exists(k0_file)) k0_tab <- readRDS(k0_file) else {
  set.seed(11)
  k0_tab <- do.call(rbind, lapply(seq_len(nrow(settings)), function(k) {
    n <- settings$n[k]; p <- settings$p[k]
    kk <- replicate(2000, { x <- gen_x(n, p); kurt(lts3(x, rnorm(n), .75)$r_best) })
    data.frame(n = n, p = p, K0 = mean(kk))
  }))
  saveRDS(k0_tab, k0_file)
}
print(k0_tab)

mgq <- function(x, y, p) {                     # Rana, Midi & Imon (2008)
  n <- length(y); c0 <- n - 2 * floor(n / 3)
  d <- data.frame(y = y, x)
  f <- robustbase::ltsReg(y ~ ., data = d)
  clean <- f$lts.wt == 1
  X <- cbind(1, x)
  fit <- lm.fit(X[clean, , drop = FALSE], y[clean])
  b <- fit$coefficients
  e <- drop(y - X %*% b)
  H <- X[clean, , drop = FALSE] %*% solve(crossprod(X[clean, , drop = FALSE])) %*% t(X[clean, , drop = FALSE])
  e[clean] <- e[clean] / (1 - diag(H))          # deletion residuals of the clean points
  g1 <- 1:floor(n / 3); g2 <- (n - floor(n / 3) + 1):n
  r <- median(e[g2]^2) / median(e[g1]^2)
  df <- floor(n / 3) - (p + 1)
  2 * min(pf(r, df, df), 1 - pf(r, df, df))     # two-sided
}

one_rep <- function(x, y, st, nu) {
  n <- st$n; p <- st$p
  d <- data.frame(y = y, x); m0 <- lm(y ~ ., data = d)
  pv <- c()
  # KaH-III exact
  ms <- sapply(parts(n), function(i) { f <- lm.fit(cbind(1, x[i, , drop = FALSE]), y[i]); sum(f$residuals^2) / (length(i) - p - 1) })
  pv["KaH3"] <- 1 - pfmax(max(ms) / min(ms), floor(n / 3) - p - 1)
  # KaH-V2 alpha .9
  s9 <- lts3(x, y, .9)$s
  pv["V2_a90"] <- 1 - pfmax((max(s9) / min(s9))^2, nu$a90)
  # KaH alpha .75: formula, kurtosis-corrected, bootstrap
  L <- lts3(x, y, .75); T75 <- max(L$s) / min(L$s)
  pv["V2_a75"] <- 1 - pfmax(T75^2, nu$a75)
  f_inf <- (kurt(L$r_best) - 1) / (nu$K0 - 1)
  pv["V2_a75_kc"] <- 1 - pfmax(T75^2, nu$a75 / max(f_inf, .25))
  r <- L$r_all
  Tb <- replicate(B, { s <- lts3(x, sample(r, n, replace = TRUE), .75)$s; max(s) / min(s) })
  pv["V2_a75_bs"] <- (1 + sum(Tb >= T75)) / (B + 1)
  # competitors
  pv["GQ"] <- lmtest::gqtest(m0, order.by = x[, 1], fraction = 1 / 3, alternative = "two.sided")$p.value
  pv["BP_Koenker"] <- lmtest::bptest(m0)$p.value
  W <- if (p == 1) cbind(x, x^2) else cbind(x, x^2, x[, 1] * x[, 2])
  pv["White"] <- lmtest::bptest(m0, ~ W)$p.value
  pv["MGQ_Rana2008"] <- tryCatch(mgq(x, y, p), error = function(e) NA)
  pv["WK_2006"] <- tryCatch(skedastic::wilcox_keselman(m0, B = 200L)$p.value, error = function(e) NA)
  pv["Zhou_2015"] <- tryCatch(skedastic::zhou_etal(m0, method = "pooled", Bperturbed = 200L)$p.value, error = function(e) NA)
  pv["LiYao_2019"] <- tryCatch(skedastic::li_yao(m0, method = "cvt")$p.value, error = function(e) NA)
  pv
}

tasks <- expand.grid(ch = seq_len(R / chunk), s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
run_task <- function(t) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic) })
  st <- settings[tasks$s[t], ]; sc <- tasks$sc[t]; n <- st$n; p <- st$p; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("n%03d_p%d_%s_ch%d.rds", n, p, sc, tasks$ch[t]))
  if (file.exists(f)) return(invisible(NULL))
  nu <- list(a90 = dfa$nu_star[dfa$a == .9 & dfa$m == m & dfa$p == p],
             a75 = dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == p],
             K0 = k0_tab$K0[k0_tab$n == n & k0_tab$p == p])
  set.seed(1e4 * t + n)
  P <- t(replicate(chunk, {
    x <- gen_x(n, p)
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    y <- drop(1 + x %*% rep(1, p) + e)
    one_rep(x, y, st, nu)
  }))
  saveRDS(data.frame(n = n, p = p, scenario = sc, P, check.names = FALSE), f)
  invisible(NULL)
}

if (nzchar(Sys.getenv("SMOKE"))) {             # quick smoke test of every branch
  chunk <- 2; B <- 9
  for (sc in scen) { tasks <- expand.grid(ch = 1, s = 1, sc = sc, stringsAsFactors = FALSE)
    out_dir <- tempdir(); run_task(1)
    print(readRDS(file.path(out_dir, sprintf("n030_p1_%s_ch1.rds", sc)))) }
  quit(save = "no")
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("settings", "tasks", "R", "chunk", "B", "gam", "dfa", "k0_tab", "out_dir",
                    "gen_x", "parts", "lts3", "kurt", "mgq", "one_rep", "run_task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(tasks)), run_task))
stopCluster(cl)

all <- do.call(rbind, lapply(list.files(out_dir, "_ch\\d+\\.rds$", full.names = TRUE), readRDS))
meth <- setdiff(names(all), c("n", "p", "scenario"))
agg <- aggregate(all[meth] < .05, by = all[c("n", "p", "scenario")], FUN = function(v) mean(v, na.rm = TRUE))
agg <- agg[order(match(agg$scenario, scen), agg$p, agg$n), ]
print(agg, digits = 3, row.names = FALSE)
cat("NA counts:\n"); print(colSums(is.na(all[meth])))
cat("MC SE at 5% with R =", R, ":", round(sqrt(.05 * .95 / R), 4), "\n")
write.csv(agg, "competitors_rejection.csv", row.names = FALSE)
