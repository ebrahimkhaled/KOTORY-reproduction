# Step 38: the KaH-III statistic T1 for every replication of the factorial study (step 28) and of the
# heavy-tail study (step 33). Their stored p-values, 1 - pfmax(T1), underflow to 0 for about 7% of
# contaminated null samples, which ties them and distorts the size-adjusted threshold. The data sets
# are regenerated with the same seeds (every foreign call in steps 28 and 33 restores the RNG, so each
# replication draws only its data) and each is checked against the stored fingerprint.
parts <- function(n) list(1:floor(n / 3), (round(n / 3) + 1):floor(2 * n / 3), ceiling(2 * n / 3 + 1):n)
sdfun <- function(shape, z) switch(shape, none = rep(1, length(z)), mono = exp(0.35 * z), U = 1 + z^2,
                                   bulge = 1 + 3 * exp(-z^2))
t1 <- function(x, y) { n <- length(y); p <- ncol(x)
  ms <- sapply(parts(n), function(i) sum(lm.fit(cbind(1, x[i, , drop = FALSE]), y[i])$residuals^2) / (length(i) - p - 1))
  max(ms) / min(ms) }

run <- function(study) {
  if (study == "factorial") {
    cells <- expand.grid(shape = c("none", "mono", "U", "bulge"), order = c("x1", "time"), p = c(1, 4),
                         cont = c(FALSE, TRUE), n = c(45, 90, 150), ch = 1:4, stringsAsFactors = FALSE)
    file_of <- function(cl) file.path("factorial", sprintf("%s_%s_p%d_%s_n%03d_ch%d.rds", cl$shape, cl$order, cl$p,
                                                           if (cl$cont) "out" else "clean", cl$n, cl$ch))
    seed0 <- 2e6
  } else {
    cells <- expand.grid(shape = c("none", "mono", "U", "bulge"), order = c("x1", "time"), p = c(1, 4),
                         n = c(45, 90, 150), ch = 1:4, stringsAsFactors = FALSE)
    cells$cont <- FALSE
    file_of <- function(cl) file.path("heavy_tails", sprintf("%s_%s_p%d_t3_n%03d_ch%d.rds", cl$shape, cl$order, cl$p, cl$n, cl$ch))
    seed0 <- 3e6
  }
  do.call(rbind, lapply(seq_len(nrow(cells)), function(k) {
    cl <- cells[k, ]; n <- cl$n; p <- cl$p
    stored <- readRDS(file_of(cl))
    set.seed(seed0 + k)
    out <- t(sapply(seq_len(nrow(stored)), function(r) {
      x <- scale(matrix(runif(n * p), n, p)); colnames(x) <- paste0("x", 1:p)
      z <- if (cl$order == "time") as.numeric(scale(seq_len(n))) else x[, 1]
      o <- order(z); x <- x[o, , drop = FALSE]; z <- z[o]
      if (study == "factorial") {
        e <- rnorm(n, 0, sdfun(cl$shape, z))
        if (cl$cont) { i <- sample(n, round(.1 * n)); e[i] <- rnorm(length(i), 0, 7) }
      } else e <- sdfun(cl$shape, z) * rt(n, 3) / sqrt(3)
      y <- drop(1 + x %*% rep(1, p) + e)
      c(fingerprint = sum(y * seq_along(y)), T1 = t1(x, y))
    }))
    if (!isTRUE(all.equal(unname(out[, "fingerprint"]), stored$fingerprint, tolerance = 1e-10)))
      stop("fingerprint mismatch in ", study, " cell ", k)
    data.frame(shape = cl$shape, order = cl$order, p = p, cont = cl$cont, n = n, ch = cl$ch,
               fingerprint = out[, "fingerprint"], T1 = out[, "T1"])
  }))
}
for (study in c("factorial", "heavy")) {
  s <- run(study)
  cat(study, ": ", nrow(s), " replications regenerated, every fingerprint matches\n", sep = "")
  saveRDS(s, sprintf("%s_kah3_T1.rds", study))
}
