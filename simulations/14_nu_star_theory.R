# Step 14: closed-form effective degrees of freedom of the raw LTS scale.
#
# Asymptotically the raw LTS scale (coverage alpha) behaves like the lower-trimmed mean of
# Y = e^2:  T = (1/alpha) * integral_0^alpha G^{-1}(u) du,  G = cdf of Y (chi2_1 under normal errors).
# Its influence function is  IF(y) = xi - T - (xi - y)_+ / alpha,  xi = G^{-1}(alpha),
# so  n * Var(T_hat) -> Var{(xi - Y)_+} / alpha^2.
# The consistency factor is a constant and cancels in the ratio of the three parts.
# Matching Var(log T_hat) to Var(log chi2_nu / nu) ~ 2 / nu gives the effective d.f. per observation
#     r(alpha) = nu* / m  ->  2 T^2 alpha^2 / Var{(xi - Y)_+}.
# Compare with the simulated nu*/nu (large-m averages) from the KOTORY table.

r_theory <- function(alpha) {
  xi <- qchisq(alpha, 1)
  f <- function(y) dchisq(y, 1)
  T  <- integrate(function(y) y * f(y), 0, xi)$value / alpha
  m1 <- integrate(function(y) (xi - y) * f(y), 0, xi)$value
  m2 <- integrate(function(y) (xi - y)^2 * f(y), 0, xi)$value
  v  <- m2 - m1^2
  2 * T^2 * alpha^2 / v
}

load("../KOTORY/R/sysdata.rda")
sim <- aggregate(ratio ~ alpha, data = .nu_table[.nu_table$m >= 50, ], FUN = mean)
out <- data.frame(alpha = sim$alpha, simulated = sim$ratio, theory = sapply(sim$alpha, r_theory))
out$rel_diff <- out$theory / out$simulated - 1
print(out, digits = 4)
cat("sanity: alpha = 1 (no trimming) gives", r_theory(0.999999), "(should be 1)\n")

# finite-m check: theory with the regression d.f. correction nu = m - p - 1
t <- .nu_table
t$theory <- sapply(t$alpha, r_theory) * t$nu
t$rel <- t$theory / t$nu_star - 1
cat("\nRelative error of theory*nu vs simulated nu*, by m band:\n")
print(aggregate(rel ~ alpha + cut(m, c(0, 10, 20, 50, 100, 200)), data = t, FUN = function(v) round(mean(v), 3)))
