# Step 17: LaTeX tables for the manuscript, generated from the result files so that every
# printed number traces to a released file (no hand-typed values).

out <- "../paper_stat_papers/tables"; dir.create(out, showWarnings = FALSE)
pct <- function(v) formatC(100 * v, format = "f", digits = 1)
a <- read.csv("all_tests_rejection.csv", check.names = FALSE)
rg <- read.csv("rgq_rejection.csv", check.names = FALSE)

tests <- c(KaH3 = "KaH-III", V2_a75 = "KaH-robust, 0.75", V2_a75_bs = "KaH-robust, bootstrap",
           V2_a90 = "KaH-robust, 0.90",
           GQ = "Goldfeld--Quandt", BP_Koenker = "Breusch--Pagan", White = "White",
           MGQ_Rana2008 = "MGQ", BRW_White = "White, LTS screen", WK_2006 = "Wilcox--Keselman",
           BAMSET_3 = "BAMSET", EvansKing = "Evans--King", HMC = "Harrison--McCabe",
           Zhou_2015 = "Zhou et al.", BF_3parts = "Brown--Forsythe, parts")
cells <- list(c(30, 1), c(60, 1), c(150, 1), c(90, 2))
cell_lab <- c("30,1", "60,1", "150,1", "90,2")

block <- function(scen, cols = names(tests)) {
  sapply(cells, function(np) { d <- a[a$n == np[1] & a$p == np[2] & a$scenario == scen, ]; pct(unlist(d[cols])) })
}
write_tab <- function(file, caption, label, scen_list, scen_names) {
  hdr <- paste0("& \\multicolumn{4}{c}{", scen_names, "}", collapse = " ")
  cm <- paste0("\\cmidrule(lr){", 2 + 4 * (seq_along(scen_list) - 1), "-", 1 + 4 * seq_along(scen_list), "}", collapse = "")
  sub <- paste(rep(paste0("& ", cell_lab, collapse = " "), length(scen_list)), collapse = " ")
  body <- do.call(cbind, lapply(scen_list, block))
  rows <- paste0(tests, " & ", apply(body, 1, paste, collapse = " & "), " \\\\")
  txt <- c("\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", caption, "}\\label{", label, "}"),
           "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
           paste0("\\begin{tabular}{l", strrep("r", 4 * length(scen_list)), "}"), "\\toprule",
           paste0(hdr, " \\\\"), cm, paste0("Test ", sub, " \\\\"), "\\midrule", rows, "\\bottomrule",
           "\\end{tabular}", "\\end{table}")
  writeLines(txt, file.path(out, file))
}
write_tab("tab_size.tex",
  paste("Empirical size (\\%) at the nominal 5\\% level under $H_0$; columns give $(n,p)$.",
        "1000 replications per cell, Monte Carlo standard error about 0.7 percentage points"),
  "tab:size", c("H0_clean", "H0_t5", "H0_outliers"),
  c("Normal errors", "$t_5$ errors", "10\\% outliers"))
write_tab("tab_power.tex",
  paste("Empirical power (\\%) at the nominal 5\\% level; error standard deviation $\\exp(0.35\\,x_1)$;",
        "columns give $(n,p)$. Rates for tests whose size is not controlled in",
        "Fig.~\\ref{fig:size} are not valid power"),
  "tab:power", c("H1_clean", "H1_outliers"), c("Normal errors", "10\\% outliers"))

# RGQ table (same data as KaH in that run)
rgt <- c(KaH3 = "KaH-III", V2_a75 = "KaH-robust", RGQ_D = "RGQ, as described",
         RGQ_x = "RGQ by $x_1$")
sc <- c(H0_clean = "$H_0$ normal", H0_outliers = "$H_0$ outliers", H0_t5 = "$H_0$ $t_5$",
        H1_clean = "$H_1$ normal", H1_outliers = "$H_1$ outliers")
rows <- unlist(lapply(names(sc), function(s) {
  d <- rg[rg$scenario == s, ]; d <- d[order(d$p, d$n), ]
  v <- sapply(names(rgt), function(tk) paste(sprintf("%.0f", 100 * d[[tk]]), collapse = "/"))   # whole % keeps the cells readable
  paste0(sc[s], " & ", paste(v, collapse = " & "), " \\\\")
}))
writeLines(c("\\begin{table}[!htbp]", "\\centering",
  "\\caption{Rejection rates (\\%) of the robust Goldfeld--Quandt test of \\citet{AlihOng2015} and of the KaH tests on the same data sets, an independent set of 1000 replications of the design of Section~\\ref{sec:design} (hence small differences from Fig.~\\ref{fig:size} and Tables~\\ref{tab:power} and \\ref{tab:sizeadj}). Each entry lists $(n,p)=(30,1)/(60,1)/(150,1)/(90,2)$, rounded to whole percentages}\\label{tab:rgq}",
  "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Scenario & ", paste(rgt, collapse = " & "), " \\\\"), "\\midrule", rows, "\\bottomrule",
  "\\end{tabular}", "\\end{table}"), file.path(out, "tab_rgq.tex"))

# real-data p-values
rd <- read.csv("real_data_pvalues.csv", check.names = FALSE)
keep <- c("KaH-III", "KaH robust (alpha = 0.75)", "KaH robust bootstrap", "Goldfeld-Quandt",
          "Breusch-Pagan (Koenker)", "White", "MGQ (Rana et al. 2008)", "White after LTS screen (BRW)",
          "RGQ (ranked by x)")
labs <- c("KaH-III", "KaH-robust", "KaH-robust, bootstrap", "Goldfeld--Quandt", "Breusch--Pagan",
          "White", "MGQ", "White, LTS screen", "RGQ by $x_1$")
sets <- c("Housing, clean", "Housing, outliers", "Savings, clean", "Savings, outliers",
          "Restaurant, clean", "Restaurant, outliers", "Consumption, clean", "Consumption, outliers")
fp <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3)))
rows <- sapply(seq_along(keep), function(i) {
  v <- sapply(sets, function(s) fp(rd$p[rd$data == s & rd$test == keep[i]]))
  paste0(labs[i], " & ", paste(v, collapse = " & "), " \\\\")
})
writeLines(c("\\begin{table}[!htbp]", "\\centering",
  "\\caption{$p$-values on four benchmark data sets, without and with the outliers planted by the authors of the competing robust tests ($n=20$, 31, 30 and 30). The bootstrap reference used $B=999$. A dash marks a test that could not be computed: in the housing data with outliers a matrix in the outlier screen of RGQ is singular, because the regressor takes few distinct values}\\label{tab:real}",
  "\\scriptsize", "\\setlength{\\tabcolsep}{2.5pt}", "\\begin{tabular}{lcccccccc}", "\\toprule",
  "& \\multicolumn{2}{c}{Housing} & \\multicolumn{2}{c}{Savings} & \\multicolumn{2}{c}{Restaurant} & \\multicolumn{2}{c}{Consumption} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Test & clean & outliers & clean & outliers & clean & outliers & clean & outliers \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file.path(out, "tab_real.tex"))

# power under contamination: nominal and size-adjusted, for the tests usable under outliers
sa <- read.csv("size_adjusted_power.csv")
satests <- c(V2_a75 = "KaH-robust", V2_a75_bs = "KaH-robust, bootstrap",
             BRW_White = "White, LTS screen", WK_2006 = "Wilcox--Keselman", MGQ_Rana2008 = "MGQ",
             BP_Koenker = "Breusch--Pagan", White = "White", GQ = "Goldfeld--Quandt", BAMSET_3 = "BAMSET")
satests <- satests[names(satests) %in% sa$test]
ord <- list(c(30, 1), c(60, 1), c(150, 1), c(90, 2))
h0o <- a[a$scenario == "H0_outliers", ]
# nominal-level power under outliers is in Table tab:power (same data sets), so only size and
# size-adjusted power are shown here
rows <- sapply(names(satests), function(t) {
  sz <- sapply(ord, function(np) pct(h0o[h0o$n == np[1] & h0o$p == np[2], t]))
  ad <- sapply(ord, function(np) pct(sa$adjusted[sa$test == t & sa$n == np[1] & sa$p == np[2]]))
  paste0(satests[t], " & ", paste(c(sz, ad), collapse = " & "), " \\\\")
})
writeLines(c("\\begin{table}[!htbp]", "\\centering",
  paste0("\\caption{Tests under 10\\% outliers: size under $H_0$ and size-adjusted power under $H_1$ ",
         "(rejection threshold set to the 5\\% quantile of the test's own $p$-values under $H_0$ with outliers); ",
         "power at the nominal level is in Table~\\ref{tab:power}; columns give $(n,p)$; 1000 replications on the same ",
         "data sets for every test}\\label{tab:sizeadj}"),
  "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{lrrrrrrrr}", "\\toprule",
  "& \\multicolumn{4}{c}{Size, $H_0$} & \\multicolumn{4}{c}{Power, size-adjusted} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  paste0("Test ", paste(rep(paste0("& ", cell_lab, collapse = " "), 2), collapse = " "), " \\\\"), "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file.path(out, "tab_sizeadj.tex"))

# effective d.f.: simulated vs limit
load("../KOTORY/R/sysdata.rda")
r_theory <- function(al) {
  xi <- qchisq(al, 1); f <- function(y) dchisq(y, 1)
  T <- integrate(function(y) y * f(y), 0, xi)$value / al
  m1 <- integrate(function(y) (xi - y) * f(y), 0, xi)$value
  m2 <- integrate(function(y) (xi - y)^2 * f(y), 0, xi)$value
  2 * T^2 * al^2 / (m2 - m1^2)
}
ms <- c(10, 20, 50, 100, 200)
rows <- sapply(c(0.5, 0.75, 0.9), function(al) {
  v <- sapply(ms, function(m) { t <- .nu_table[.nu_table$alpha == al & .nu_table$m == m, ]; formatC(mean(t$ratio), format = "f", digits = 3) })
  paste0(formatC(al, format = "f", digits = 2), " & ", paste(v, collapse = " & "), " & ",
         formatC(r_theory(al), format = "f", digits = 3), " \\\\")
})
writeLines(c("\\begin{table}[!htbp]", "\\centering",
  "\\caption{Simulated ratio $\\nu^*/\\nu$ (mean over the feasible $p\\le5$, all five for $m\\ge15$; 20{,}000 replications per cell) and its limit $r(\\alpha)$ from Proposition~\\ref{prop:nustar}}\\label{tab:nustar}",
  "\\small", "\\begin{tabular}{lcccccc}", "\\toprule",
  paste0("$\\alpha$ & ", paste0("$m=", ms, "$", collapse = " & "), " & $r(\\alpha)$ \\\\"), "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file.path(out, "tab_nustar.tex"))
# narrow the wide tables: sn-jnl wraps table bodies internally, so \resizebox cannot be used;
# use a smaller font and tighter column separation instead
for (f in c("tab_size.tex", "tab_power.tex", "tab_sizeadj.tex", "tab_real.tex", "tab_rgq.tex")) {
  x <- readLines(file.path(out, f))
  x <- sub("^\\\\setlength\\{\\\\tabcolsep\\}\\{[0-9.]+pt\\}$", "\\\\setlength{\\\\tabcolsep}{2pt}", x)
  writeLines(x, file.path(out, f))
}
cat("tables written:", paste(list.files(out), collapse = ", "), "\n")
