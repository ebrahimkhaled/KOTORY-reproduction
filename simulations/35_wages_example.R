# Step 35: the larger example of Section 4. SLID wages (carData), complete cases, wages regressed on
# education and age, observations sorted by education. Runs the KaH tests and the battery, and
# reports the residual standard deviation of each third.
suppressMessages({ pkgload::load_all("../KOTORY", quiet = TRUE); library(robustbase) })
d <- na.omit(carData::SLID[, c("wages", "education", "age")])
cat("persons:", nrow(d), "\n")
print(run.all.het(wages ~ education + age, data = d, order.by = "education"))
k3 <- kah3.test(wages ~ education + age, data = d, order.by = "education")
cat("KaH-III error scale by third of education:", sprintf("%.1f", sqrt(k3$estimate)), "\n")
