# Step 18: is skedastic::bamset() oversized under a clean normal null, or is it our call?
# Checks deflator handling, the correction factor and the margin-omission option.
suppressMessages(library(skedastic))
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
set.seed(11)
R <- 1000; n <- 60
res <- replicate(R, {
  x <- sort(scale(runif(n))[, 1]); y <- 1 + x + rnorm(n)
  m <- lm(y ~ x1, data = data.frame(y = y, x1 = x))
  c(default   = keep_rng(bamset(m, k = 3, deflator = "x1")$p.value),
    nocorrect = keep_rng(bamset(m, k = 3, deflator = "x1", correct = FALSE)$p.value),
    col2      = keep_rng(bamset(m, k = 3, deflator = 2)$p.value))
})
cat("Rejection rate at 5% under H0 (normal errors, n = 60), R =", R, "(MC SE 0.7 pp)\n")
print(round(100 * rowMeans(res < 0.05), 1))
cat("\nQuantiles of the default p-values (uniform under a valid test):\n")
print(round(quantile(res["default", ], c(.05, .25, .5, .75, .95)), 3))
