# Step 19: BAMSET (Ramsey 1969) from its definition -- Bartlett's test on three ordered groups
# of BLUS residuals -- on the step-11 settings and scenarios. Replaces the skedastic::bamset()
# column, which is miscalibrated (21% size under a clean normal null; step 18b).
suppressMessages(library(skedastic))
library(parallel)
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
own_bamset <- function(m) {
  e <- keep_rng(skedastic::blus(m, omit = "last", keepNA = FALSE))
  bartlett.test(e, cut(seq_along(e), 3, labels = FALSE))$p.value
}
settings <- data.frame(n = c(30, 60, 150, 90), p = c(1, 1, 1, 2))
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
grid <- expand.grid(s = seq_len(nrow(settings)), sc = scen, stringsAsFactors = FALSE)
R <- 1000; gam <- 0.35
cell <- function(g) {
  suppressMessages(library(skedastic))
  n <- settings$n[grid$s[g]]; p <- settings$p[grid$s[g]]; sc <- grid$sc[g]
  set.seed(3e5 + 100 * g)
  pv <- replicate(R, {
    x <- scale(matrix(runif(n * p), n, p)); x <- x[order(x[, 1]), , drop = FALSE]
    sdv <- if (startsWith(sc, "H1")) exp(gam * x[, 1]) else rep(1, n)
    e <- if (sc == "H0_t5") rt(n, 5) * sqrt(3 / 5) else rnorm(n, 0, sdv)
    if (grepl("outliers", sc)) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
    colnames(x) <- paste0("x", 1:p)
    own_bamset(lm(y ~ ., data = data.frame(y = drop(1 + x %*% rep(1, p) + e), x)))
  })
  data.frame(n = n, p = p, scenario = sc, BAMSET_own = mean(pv < 0.05))
}
cl <- makeCluster(min(20, detectCores() - 2))
clusterExport(cl, c("settings", "grid", "R", "gam", "keep_rng", "own_bamset"))
res <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(grid)), cell))
stopCluster(cl)
print(res, row.names = FALSE)
a <- read.csv("all_tests_rejection.csv", check.names = FALSE)
a$BAMSET_3 <- NULL
a <- merge(a, res, by = c("n", "p", "scenario"), sort = FALSE)
names(a)[names(a) == "BAMSET_own"] <- "BAMSET_3"
write.csv(a, "all_tests_rejection.csv", row.names = FALSE)
cat("all_tests_rejection.csv: BAMSET_3 column replaced by the definition-based version\n")
