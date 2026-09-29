# Exact Hartley Fmax distribution by numerical integration.
# S_1..S_k iid chi2_nu (scale cancels). Fmax = max/min.
# P(Fmax <= c) = k * int_0^Inf f(x) [G(c x) - G(x)]^(k-1) dx
# (condition on the smallest one being at x; all others in [x, c x]).
# Integrated on the probability scale u = G(x), which is bounded and smooth
# for any nu > 0 (the x-scale form blows up at 0 when nu < 2):
#   P(Fmax <= c) = k * int_0^1 [G(c * G^-1(u)) - u]^(k-1) du

pfmax <- function(c, nu, k = 3) {
  mapply(function(cc, nu) {
    if (cc <= 1) return(0)
    g <- function(u) k * pmax(pchisq(cc * qchisq(u, nu), nu) - u, 0)^(k - 1)
    integrate(g, 0, 1, rel.tol = 1e-10, subdivisions = 2000L)$value
  }, c, nu)
}

qfmax <- function(prob, nu, k = 3) {
  mapply(function(pr, v) {
    f <- function(lc) pfmax(exp(lc), v, k) - pr
    exp(uniroot(f, c(1e-9, log(1e12)), tol = 1e-12)$root)
  }, prob, nu)
}
