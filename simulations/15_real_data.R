# Step 15: real-data examples for the paper. Every test on the same data, clean and with the
# outliers planted by the authors of the competing tests (their own home ground), plus one
# data set with a genuine (not planted) outlier.
#   housing     Pindyck & Rubinfeld; Rana, Midi & Imon (2008) planted y1 = 4.9, y20 = 2.0
#   consumption Gujarati; Alih & Ong (2015) Table 7, planted cases 1, 2, 30 (Rana et al. 2008 use other values)
#   savings     Koutsoyiannis; Alih & Ong (2015) Table 1, planted cases 1, 2, 30, 31
#   restaurant  Montgomery et al.; Alih & Ong (2015) Table 3, planted cases 1, 26, 30
#   education   robustbase::education (Chatterjee & Price), 50 US states, 3 regressors, Alaska

suppressMessages({ devtools::load_all("../KOTORY", quiet = TRUE); library(robustbase) })
source("rgq_alih_ong.R")

white_p <- function(x, e) {
  W <- cbind(x, x^2)
  if (ncol(x) > 1) for (i in 1:(ncol(x) - 1)) for (j in (i + 1):ncol(x)) W <- cbind(W, x[, i] * x[, j])
  u <- e^2; r2 <- 1 - sum(lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  pchisq(length(u) * r2, ncol(W), lower.tail = FALSE)
}
brw_white <- function(x, y) {
  keep <- ltsReg(x = x, y = y, mcd = FALSE)$lts.wt == 1
  white_p(x[keep, , drop = FALSE], lm.fit(cbind(1, x[keep, , drop = FALSE]), y[keep])$residuals)
}

run_all <- function(label, x, y) {
  x <- as.matrix(x); colnames(x) <- paste0("x", seq_len(ncol(x)))
  d <- data.frame(y = y, x); f <- reformulate(colnames(x), "y")
  set.seed(1)
  b <- run.all.het(f, data = d, include.optional = FALSE)
  pv <- setNames(b$p_value, b$Test)
  pv["KaH robust bootstrap"] <- kah.robust.test(f, data = d, method = "bootstrap", B = 999)$p.value
  pv["RGQ (Alih-Ong, as written)"] <- tryCatch(rgq_test(x, y, "D")$p.value, error = function(e) NA)
  pv["RGQ (ranked by x)"] <- tryCatch(rgq_test(x, y, "x")$p.value, error = function(e) NA)
  pv["White after LTS screen (BRW)"] <- brw_white(x, y)
  data.frame(data = label, n = length(y), test = names(pv), p = unname(pv))
}

sets <- list()
x <- c(5,5,5,5,5,10,10,10,10,10,15,15,15,15,15,20,20,20,20,20)
y <- c(1.8,2,2,2,2.1,3.1,3.2,3.5,3.5,3.6,4.2,4.2,4.5,4.8,5,4.8,5,5.7,6,6.2)
sets$`Housing, clean` <- list(x, y); yo <- y; yo[c(1, 20)] <- c(4.9, 2.0); sets$`Housing, outliers` <- list(x, yo)

x <- c(80,100,85,110,120,115,130,140,125,90,105,160,150,165,145,180,225,200,240,185,220,210,245,260,190,205,265,270,230,250)
y <- c(55,65,70,80,79,84,98,95,90,75,74,110,113,125,108,115,140,120,145,130,152,144,175,180,135,140,178,191,137,189)
sets$`Consumption, clean` <- list(x, y)
sets$`Consumption, outliers` <- list(x, c(12,11,70,80,79,84,98,95,90,75,74,110,113,125,108,115,140,120,145,130,152,144,175,180,135,140,178,191,137,89))

x <- c(8777,9210,9954,10508,10979,11912,12747,13499,14269,15522,16730,17663,18575,19635,21163,22880,24127,25604,26500,27670,28300,27430,29560,28150,32100,32500,35250,33500,36000,36200,38200)
y <- c(264,105,90,131,122,107,406,503,431,588,898,950,779,819,1222,1702,1578,1654,1400,1829,2200,2017,2105,1600,2250,2420,2570,1720,1900,2100,2300)
sets$`Savings, clean` <- list(x, y); yo <- y; yo[c(1, 2, 30, 31)] <- c(2644, 1050, 2.1, 2.3); sets$`Savings, outliers` <- list(x, yo)

x <- c(3000,3150,3085,5225,5350,6090,8925,9015,8885,8950,9000,11345,12275,12400,12525,12310,13700,15000,15175,14995,15050,15200,15150,16800,16500,17830,19500,19200,19000,19350)
y <- c(81464,72661,72344,90743,98588,96507,126574,114133,115814,123181,131434,140564,151352,146926,130963,146630,147041,179021,166200,180732,178187,185304,155931,172579,188851,192424,203112,192482,218715,214317)
sets$`Restaurant, clean` <- list(x, y); yo <- y; yo[c(1, 26, 30)] <- c(814644, 392424, 21431); sets$`Restaurant, outliers` <- list(x, yo)

data(education, package = "robustbase")
ed <- education[order(education$X1), ]
sets$`Education (Alaska kept)` <- list(ed[, c("X1", "X2", "X3")], ed$Y)
noAK <- ed[ed$State != "AK", ]
sets$`Education (Alaska removed)` <- list(noAK[, c("X1", "X2", "X3")], noAK$Y)

res <- do.call(rbind, lapply(names(sets), function(s)
  withCallingHandlers(run_all(s, sets[[s]][[1]], sets[[s]][[2]]),
                      error = function(e) message("FAILED on data set: ", s, " -- ", conditionMessage(e),
                                                  "\n  in: ", paste(deparse(conditionCall(e)), collapse = " ")))))
w <- reshape(res[, c("data", "test", "p")], idvar = "test", timevar = "data", direction = "wide")
names(w) <- sub("^p\\.", "", names(w))
print(format(w, digits = 2), row.names = FALSE)
write.csv(res, "real_data_pvalues.csv", row.names = FALSE)
