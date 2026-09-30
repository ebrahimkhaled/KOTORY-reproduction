# Step 40: summary of the concentrated-contamination study (step 39). Size = rejection rate at the
# nominal 5% level under H0; power = size-adjusted, with the threshold from the H0 cell of the same
# block, level and n (KaH tests on the statistic, White after an LTS screen on its p-value).
suppressMessages(pkgload::load_all("../KOTORY", quiet = TRUE))
a <- readRDS("masking_all.rds")
cat("replications per cell:", range(table(paste(a$block, a$level, a$hyp, a$n))), "\n")
nu <- function(n, alpha) kah.nu.star(floor(n / 3), 1, alpha)
a$p_KaH3 <- pfmax(a$T_KaH3, floor(a$n / 3) - 2, 3, lower.tail = FALSE)
a$p_R75 <- mapply(function(t, n) pfmax(t, nu(n, 0.75), 3, lower.tail = FALSE), a$T_R75, a$n)
a$p_R50 <- mapply(function(t, n) pfmax(t, nu(n, 0.5), 3, lower.tail = FALSE), a$T_R50, a$n)
tests <- c(KaH3 = "T_KaH3", R75 = "T_R75", R50 = "T_R50", WhiteLTS = "p_WhiteLTS")
key <- function(d) paste(d$block, d$level, d$n)
out <- do.call(rbind, lapply(split(a, key(a)), function(d) {
  h0 <- d[d$hyp == "H0", ]; h1 <- d[d$hyp == "H1", ]
  row <- data.frame(block = d$block[1], level = d$level[1], n = d$n[1])
  for (t in names(tests)) {
    pcol <- switch(t, KaH3 = "p_KaH3", R75 = "p_R75", R50 = "p_R50", WhiteLTS = "p_WhiteLTS")
    row[[paste0(t, ".size")]] <- 100 * mean(h0[[pcol]] < 0.05, na.rm = TRUE)
    s <- tests[[t]]
    if (t == "WhiteLTS") { thr <- quantile(h0[[s]], 0.05, type = 1, na.rm = TRUE); pw <- mean(h1[[s]] <= thr, na.rm = TRUE)
    } else { T0 <- sort(h0[[s]], decreasing = TRUE); thr <- T0[ceiling(0.05 * length(T0))]; pw <- mean(h1[[s]] >= thr) }
    row[[paste0(t, ".adj")]] <- 100 * pw
  }
  row
}))
out <- out[order(out$block, out$level, out$n), ]
print(format(out, digits = 3), row.names = FALSE)
saveRDS(out, "masking_summary.rds")

# LaTeX table: size and size-adjusted power, n = 90 and 150 (n = 45 in the deposit)
f0 <- function(x) sprintf("%.0f", floor(x + 0.5))
lab <- function(b, l) switch(b, mask = sprintf("%.0f\\%% of the part", 100 * l),
                             rate = sprintf("%.0f\\%% of $n$", 100 * l), onesided = "10\\% of $n$, one-sided")
sel <- out[out$n %in% c(90, 150), ]
blocks <- list(mask = "Outliers in the low-variance part", rate = "Outliers at random positions", onesided = "Asymmetric outliers")
rows <- c()
for (b in names(blocks)) {
  rows <- c(rows, sprintf("\\multicolumn{10}{@{}l}{\\emph{%s}} \\\\", blocks[[b]]))
  for (l in sort(unique(sel$level[sel$block == b]))) for (nn in c(90, 150)) {
    r <- sel[sel$block == b & sel$level == l & sel$n == nn, ]
    rows <- c(rows, sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s \\\\", if (nn == 90) lab(b, l) else "", nn,
      f0(r$R75.size), f0(r$R50.size), f0(r$KaH3.size), f0(r$WhiteLTS.size),
      f0(r$R75.adj), f0(r$R50.adj), f0(r$KaH3.adj), f0(r$WhiteLTS.adj)))
  }
}
tex <- c("\\begin{table}[!htbp]", "\\centering",
  "\\caption{Concentrated and asymmetric contamination: rejection rate (\\%) at the nominal 5\\% level under $H_0$ and size-adjusted power (\\%) under a variance monotone in $x_1$ (standard deviation $\\exp(0.35x_1)$), $p=1$, 1000 replications per cell. Outliers are $N(0,7^2)$ errors, or $N(10,1)$ for the asymmetric case. R75 and R50: \\KaHR{} with $\\alpha=0.75$ and 0.5; White-LTS: White's test after an LTS screen. Results for $n=45$ are in the deposit}\\label{tab:masking}",
  "\\scriptsize\\setlength{\\tabcolsep}{2.5pt}",
  "\\begin{tabular}{@{}lrrrrrrrrr@{}}", "\\toprule",
  " & & \\multicolumn{4}{c}{Size, $H_0$} & \\multicolumn{4}{c}{Power, size-adjusted} \\\\",
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
  "Contamination & $n$ & R75 & R50 & KaH-III & White-LTS & R75 & R50 & KaH-III & White-LTS \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, "../paper_stat_papers/tables/tab_masking.tex")
cat("tab_masking.tex written\n")

# check for Sect. 3.5: t3 study, variance monotone in x1
h <- readRDS("heavy_tails_summary.rds")$adj
k <- h$shape == "mono" & h$order == "x1"
cat("\nt3, monotone along x1:\n"); print(round(colMeans(h[k, c("KaHrobust", "EvansKing", "WK", "MGQ", "GQ")], na.rm = TRUE), 1))
