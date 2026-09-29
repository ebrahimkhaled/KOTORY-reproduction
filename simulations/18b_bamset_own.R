# Step 18b: Ramsey's BAMSET built from its definition -- Bartlett's M test on k = 3 ordered
# groups of BLUS residuals -- using skedastic's own BLUS residuals and stats::bartlett.test.
# If this holds 5% while skedastic::bamset() gives ~21%, the package's statistic is at fault.
suppressMessages(library(skedastic))
keep_rng <- function(expr) { old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv())); expr }
own_bamset <- function(m, k = 3) {
  e <- keep_rng(skedastic::blus(m, omit = "last", keepNA = FALSE))
  g <- cut(seq_along(e), k, labels = FALSE)            # ordered, equal groups
  bartlett.test(e, g)$p.value
}
set.seed(12)
R <- 1000; n <- 60
res <- replicate(R, {
  x <- sort(scale(runif(n))[, 1]); y <- 1 + x + rnorm(n)
  m <- lm(y ~ x1, data = data.frame(y = y, x1 = x))
  c(skedastic = keep_rng(bamset(m, k = 3, deflator = "x1")$p.value), own = own_bamset(m))
})
cat("Rejection rate at 5% under H0 (normal, n = 60), R =", R, "(MC SE 0.7 pp):\n")
print(round(100 * rowMeans(res < 0.05), 1))
set.seed(13)
m <- lm(y ~ x1, data = data.frame(x1 = sort(runif(60)), y = rnorm(60)))
b <- bamset(m, k = 3, deflator = "x1"); cat("\nskedastic bamset statistic:", b$statistic, " df:", b$parameter, "\n")
e <- skedastic::blus(m, omit = "last", keepNA = FALSE)
cat("own Bartlett statistic on BLUS thirds:", bartlett.test(e, cut(seq_along(e), 3, labels = FALSE))$statistic, "\n")
