# collect whatever factorial blocks exist (smoke test of step 29 on partial results)
all <- do.call(rbind, lapply(list.files("factorial", "[.]rds$", full.names = TRUE), readRDS))
saveRDS(all, "factorial_all.rds")
cat(nrow(all), "rows; shapes:", paste(unique(all$shape), collapse = ","), "\n")
