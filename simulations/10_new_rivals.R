# Step 10: the rivals found by the literature search (06_literature/rivals_2015_2026.md),
# on the same settings and scenarios as step 9. KaH3 and KaH-robust (alpha .75,
# Fmax reference) are re-run on the SAME data so the rows are directly comparable.
#   BAMSET_3     Ramsey (1969), k = 3 ordered groups, skedastic::bamset
#   EvansKing    Evans & King (1988), GLS version, skedastic::evans_king
#   HMC          Harrison & McCabe (1979), lmtest::hmctest (simulated p, 1000 draws;
#                skedastic::harrison_mccabe gives a wrong exact p-value, checked 2026-09-28)
#   BRW_White    Berenguer-Rico & Wilms (2021): LTS outlier screen, then White's test on the rest
#   BF_3parts    Brown-Forsythe (median Levene) on the residuals of our three separate fits
# Every rep stores all p-values (per-task .rds, resume-safe).

source("fmax_exact.R")
library(parallel)
out_dir <- "new_rivals"; dir.create(out_dir, showWarnings = FALSE)

settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
R <- 1000; chunk <- 250; gam <- 0.35
dfa <- read.csv("alpha_check_df.csv")

gen_x <- function(n, p) { x <- scale(matrix(runif(n * p), n, p)); x[order(x[, 1]), , drop = FALSE] }
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)

white_p <- function(x, e) {
  W <- if (ncol(x) == 1) cbind(x, x^2) else cbind(x, x^2, x[, 1] * x[, 2])
  u <- e^2; r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  pchisq(length(u) * r2, ncol(W), lower.tail = FALSE)
}

one_rep <- function(x, y, st, nu75) {
  n <- st$n; p <- st$p
  colnames(x) <- paste0("x", 1:p)
  d <- data.frame(y = y, x); m0 <- lm(y ~ ., data = d)
  pv <- c()
  fits <- lapply(parts(n), function(i) lm.fit(cbind(1, x[i, , drop = FALSE]), y[i]))
  ms <- sapply(seq_along(fits), function(j) sum(fits[[j]]$residuals^2) / (length(parts(n)[[j]]) - p - 1))
  pv["KaH3"] <- 1 - pfmax(max(ms) / min(ms), floor(n / 3) - p - 1)
  s <- sapply(parts(n), function(i) robustbase::ltsReg(x = x[i, , drop = FALSE], y = y[i], alpha = .75, mcd = FALSE)$raw.scale)
  pv["V2_a75"] <- 1 - pfmax((max(s) / min(s))^2, nu75)
  pv["BAMSET_3"] <- tryCatch(skedastic::bamset(m0, k = 3, deflator = "x1")$p.value, error = function(e) NA)
  pv["EvansKing"] <- tryCatch(skedastic::evans_king(m0, deflator = "x1", method = "GLS")$p.value, error = function(e) NA)
  pv["HMC"] <- tryCatch(lmtest::hmctest(m0, order.by = ~ x1, data = d, point = 0.5, sim = 1000)$p.value, error = function(e) NA)
  pv["BRW_White"] <- tryCatch({
    keep <- robustbase::ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
    e <- lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals
    white_p(x[keep, , drop = FALSE], e)
  }, error = function(e) NA)
  res <- unlist(lapply(fits, function(f) f$residuals))
  g <- factor(rep(1:3, sapply(parts(n), length)))
  z <- abs(res - ave(res, g, FUN = median))
  pv["BF_3parts"] <- anova(lm(z ~ g))[1, "Pr(>F)"]
  pv
}

tasks <- expand.grid(ch = seq_len(R / chunk), s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
run_task <- function(t) {
  suppressMessages({ library(robustbase); library(lmtest); library(skedastic) })
  st <- settings[tasks$s[t], ]; sc <- tasks$sc[t]; n <- st$n; p <- st$p; m <- floor(n / 3)
  f <- file.path(out_dir, sprintf("n%03d_p%d_%s_ch%d.rds", n, p, sc, tasks$ch[t]))
  if (file.exists(f)) return(invisible(NULL))
  nu75 <- dfa$nu_star[dfa$a == .75 & dfa$m == m & dfa$p == p]
  set.seed(5e4 + 1e3 * t + n)
  P <- t(replicate(chunk, {
    x <- gen_x(n, p)
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    one_rep(x, drop(1 + x %*% rep(1, p) + e), st, nu75)
  }))
  saveRDS(data.frame(n = n, p = p, scenario = sc, P, check.names = FALSE), f)
  invisible(NULL)
}

if (nzchar(Sys.getenv("SMOKE"))) {
  chunk <- 2
  for (sc in scen) { tasks <- expand.grid(ch = 1, s = 4, sc = sc, stringsAsFactors = FALSE)
    out_dir <- tempdir(); run_task(1)
    print(readRDS(file.path(out_dir, sprintf("n090_p2_%s_ch1.rds", sc)))) }
  quit(save = "no")
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("settings", "tasks", "R", "chunk", "gam", "dfa", "out_dir",
                    "gen_x", "parts", "white_p", "one_rep", "run_task", "pfmax"))
invisible(parLapplyLB(cl, seq_len(nrow(tasks)), run_task))
stopCluster(cl)

all <- do.call(rbind, lapply(list.files(out_dir, "_ch\\d+\\.rds$", full.names = TRUE), readRDS))
meth <- setdiff(names(all), c("n", "p", "scenario"))
agg <- aggregate(all[meth] < .05, by = all[c("n", "p", "scenario")], FUN = function(v) mean(v, na.rm = TRUE))
agg <- agg[order(match(agg$scenario, scen), agg$p, agg$n), ]
print(agg, digits = 3, row.names = FALSE)
cat("NA counts:\n"); print(colSums(is.na(all[meth])))
write.csv(agg, "new_rivals_rejection.csv", row.names = FALSE)
