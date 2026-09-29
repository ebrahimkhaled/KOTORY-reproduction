# Robust Goldfeld-Quandt (RGQ) test of Alih & Ong (2015), J. Appl. Stat. 42(8), 1617-1634,
# coded from Sections 3.1-3.2 of the paper (no public code exists).
#
# Algorithm 1 (outlier identification, forward search on Z = [X, y], k = p + 1 columns):
#   h = floor((n + k + 1) / 2); m = coordinatewise median; A = (n-1)^-1 sum (Z_i - m)(Z_i - m)'.
#   Step 1: d_i(m, A); mean/cov of the h smallest -> d_i(Zbar_h, S_h); basic subset = k + 1 smallest.
#   Steps 2-3: grow the basic subset one point at a time (distances recomputed from its mean/cov)
#   until it holds h points.  Step 4: with r points, if the squared (r+1)-th smallest distance
#   >= c_f * chi2_{1-level, k}, flag every point with squared distance >= that bound; else grow.
#   c_f = 1 + (k+1)/(n-k) + 1/(n-h-k).
# Algorithm 2 (RGQ score): rank by the final distance D_i ("in the same way as the classical GQ
#   test, rank the observations according to the value of D_i"), split into two halves, fit OLS
#   to the non-outlying points of each half, PRESS = sum of squared leave-one-out residuals,
#   RGQ = larger PRESS / smaller PRESS ~ F((n-2k)/2, (n-2k)/2).
#
# Reading choices (the paper leaves them open), documented:
#   * the chi-square bound uses level = 0.05 unless stated;
#   * distances are compared on the squared scale on both sides of Step 4(b);
#   * PRESS is summed over the non-outlying points of each half;
#   * the p-value is two-sided, 2 * P(F > larger/smaller), capped at 1.
#   * rank = "D" follows the paper literally; rank = "x" ranks by the first regressor
#     (the classical GQ ordering), reported as a charitable variant.

rgq_outliers <- function(Z, level = 0.05, literal = TRUE) {
  # literal = TRUE follows Step 4(b) word for word: stop when the SQUARED (r+1)-th distance
  # reaches the bound, then flag points whose (unsquared) distance H reaches the same bound.
  Z <- as.matrix(Z); n <- nrow(Z); k <- ncol(Z)
  h <- floor((n + k + 1) / 2)
  cf <- 1 + (k + 1) / (n - k) + 1 / (n - h - k)
  bound <- cf * qchisq(1 - level, k)
  md <- function(S) {                        # squared distances from the mean/cov of subset S
    mu <- colMeans(Z[S, , drop = FALSE]); V <- cov(Z[S, , drop = FALSE])
    mahalanobis(Z, mu, V)
  }
  med <- apply(Z, 2, median)
  A <- crossprod(sweep(Z, 2, med)) / (n - 1)
  d0 <- mahalanobis(Z, med, A)
  S <- order(d0)[1:h]
  d <- md(S)
  S <- order(d)[1:(k + 1)]
  repeat {                                   # Steps 2-3: grow to h points
    d <- md(S)
    if (length(S) >= h) break
    S <- order(d)[1:(length(S) + 1)]
  }
  repeat {                                   # Steps 4-5
    r <- length(S)
    d <- md(S)
    if (r + 1 > n) return(list(outlier = rep(FALSE, n), D = sqrt(d)))
    d_next <- sort(d)[r + 1]
    if (d_next >= bound) return(list(outlier = (if (literal) sqrt(d) else d) >= bound, D = sqrt(d)))
    S <- order(d)[1:(r + 1)]
  }
}

rgq_test <- function(X, y, rank = c("D", "x"), level = 0.05, literal = TRUE) {
  rank <- match.arg(rank)
  X <- as.matrix(X); n <- length(y); k <- ncol(X) + 1
  o <- rgq_outliers(cbind(X, y), level, literal)
  ord <- if (rank == "D") order(o$D) else order(X[, 1])
  half <- floor(n / 2)
  groups <- list(ord[1:half], ord[(n - half + 1):n])
  press <- sapply(groups, function(g) {
    g <- g[!o$outlier[g]]
    Xg <- cbind(1, X[g, , drop = FALSE])
    f <- lm.fit(Xg, y[g])
    hat <- rowSums((Xg %*% solve(crossprod(Xg))) * Xg)
    sum((f$residuals / (1 - hat))^2)
  })
  stat <- max(press) / min(press)
  df <- (n - 2 * k) / 2
  list(statistic = stat, df = df, p.value = min(1, 2 * pf(stat, df, df, lower.tail = FALSE)),
       n.outliers = sum(o$outlier))
}

if (identical(Sys.getenv("RGQ_SELFTEST"), "1")) {                     # check against the paper's Table 1-2 (savings data)
  x <- c(8777, 9210, 9954, 10508, 10979, 11912, 12747, 13499, 14269, 15522, 16730, 17663, 18575,
         19635, 21163, 22880, 24127, 25604, 26500, 27670, 28300, 27430, 29560, 28150, 32100, 32500,
         35250, 33500, 36000, 36200, 38200)
  y <- c(264, 105, 90, 131, 122, 107, 406, 503, 431, 588, 898, 950, 779, 819, 1222, 1702, 1578, 1654,
         1400, 1829, 2200, 2017, 2105, 1600, 2250, 2420, 2570, 1720, 1900, 2100, 2300)
  yo <- y; yo[c(1, 2, 30, 31)] <- c(2644, 1050, 2.1, 2.3)
  for (lit in c(TRUE, FALSE)) for (lev in c(0.05, 0.025, 0.01)) for (rk in c("D", "x")) {
    a <- tryCatch(rgq_test(x, y, rk, lev, lit), error = function(e) NULL)
    b <- tryCatch(rgq_test(x, yo, rk, lev, lit), error = function(e) NULL)
    f <- function(r) if (is.null(r)) "   (singular fit)          " else
      sprintf("RGQ %7.4f p %.4f (out %2d)", r$statistic, r$p.value, r$n.outliers)
    cat(sprintf("literal %-5s level %.3f rank %s | clean %s | outliers %s\n", lit, lev, rk, f(a), f(b)))
  }
  cat("Paper Table 2: clean RGQ 4.9820 p 0.0021 | outliers RGQ 5.0689 p 0.0051\n")
}
