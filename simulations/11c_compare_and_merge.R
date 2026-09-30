# Step 11c: check that step 11b reproduced step 11 exactly (same data sets, same p-values for every
# test that both runs computed), then write the rejection-rate file used by the tables and figures,
# with BAMSET from its definition on the SAME data sets as every other test (column BAMSET_own of
# step 11b replaces the separate-data BAMSET rates of step 19).
rd <- function(dir) do.call(rbind, lapply(list.files(dir, "_ch\\d+\\.rds$", full.names = TRUE), readRDS))
old <- rd("all_tests"); new <- rd("all_tests_v2")
stopifnot(nrow(old) == nrow(new))
key <- function(d) paste(d$n, d$p, d$scenario, round(d$fingerprint, 8))
new <- new[match(key(old), key(new)), ]
stopifnot(!anyNA(new$fingerprint))
common <- setdiff(intersect(names(old), names(new)), c("n", "p", "scenario", "fingerprint"))
same <- sapply(common, function(t) isTRUE(all.equal(old[[t]], new[[t]], tolerance = 1e-12)))
print(same)
if (!all(same)) stop("step 11b did not reproduce step 11 for: ", paste(common[!same], collapse = ", "))
cat("all", length(common), "common columns identical over", nrow(old), "replications\n")

if (!file.exists("all_tests_rejection_step19.csv")) file.copy("all_tests_rejection.csv", "all_tests_rejection_step19.csv")
meth <- setdiff(names(new), c("n", "p", "scenario", "fingerprint"))
agg <- aggregate(new[meth] < .05, by = new[c("n", "p", "scenario")], FUN = function(v) mean(v, na.rm = TRUE))
agg$BAMSET_skedastic <- agg$BAMSET_3          # the miscalibrated skedastic::bamset(), kept for the record
agg$BAMSET_3 <- agg$BAMSET_own                 # the column name the table and figure scripts read
scen <- c("H0_clean", "H0_outliers", "H0_t5", "H1_clean", "H1_outliers")
agg <- agg[order(match(agg$scenario, scen), agg$p, agg$n), ]
write.csv(agg, "all_tests_rejection.csv", row.names = FALSE)
b19 <- read.csv("all_tests_rejection_step19.csv", check.names = FALSE)
cmp <- merge(b19[c("n", "p", "scenario", "BAMSET_3")], agg[c("n", "p", "scenario", "BAMSET_3")],
             by = c("n", "p", "scenario"), suffixes = c(".step19", ".same_data"))
print(transform(cmp, BAMSET_3.step19 = round(100 * BAMSET_3.step19, 1), BAMSET_3.same_data = round(100 * BAMSET_3.same_data, 1)),
      row.names = FALSE)
