####### Bayesian Additive Regression Trees #######
# Updated: 2026-03-25
rm(list = ls()); source("load.r")
library(BART);library(dplyr); rds.path <- file.path(getwd(), "data/rds")

data <- readRDS(file.path(rds.path, "X_y.rds"))
X <- data$X; y <- data$y

# ============================================================================
set.seed(123)

mis.count <- colSums(is.na(X))
if (any(mis.count > 0)) {
  print(data.frame(
    var = names(mis.count[mis.count > 0]),
    count = mis.count[mis.count > 0]
  ))
} else {
  cat("No missing data!\n")
}

# ==============================================================================
# FIT BART MODEL
# ==============================================================================
# lbart is optimized for binary outcomes (logistic BART)

cat("\n=== Fitting BART model ===\n")
cat("Outcome: icu_mortality\n")
cat("N =", nrow(X), ", P =", ncol(X), "\n")
cat("Event rate:", round(mean(as.numeric(as.character(y)), na.rm = TRUE) * 100, 1), "%\n")

bart.model <- lbart(
  x.train = X,
  y.train = y,
  ndpost = 2000,    # Number of posterior draws
  nskip = 200,      # Burn-in iterations
  keepevery = 1,    # Keep every draw
  printevery = 100  # Progress updates
)

# ==============================================================================
# EXAMINE RESULTS
# ==============================================================================

cat("\n=== BART Model Summary ===\n")

# model free variable selection
var.importance <- bart.model$varcount / rowSums(bart.model$varcount)
mean.var.importance <- colMeans(var.importance)
hist(mean.var.importance)

top.vars <- sort(mean.var.importance, decreasing = TRUE)[1:20]
par(mar = c(15, 4, 4, 2))  # Adjust margins: bottom=15, left=4, top=4, right=2
barplot(top.vars, 
        las = 2, 
        main="BART Variable Inclusion Freq.",
        ylab="Inclusion Proportion", 
        col="steelblue")


# extracting patient phoenotypes (heterogeneity discovery)
patient.risk.means <- bart.model$prob.train.mean

# To find "phenotypes", we can cluster patients based on their full posterior risk profiles 
# or simply stratify them by their individualized conditional expectation (CATE).
# Here, we use k-means on the posterior distributions to group patients with similar risk shapes

risk.matrix <- t(bart.model$prob.train) # Transpose to have patients as rows
phoenotype.class <- kmeans(risk.matrix, centers = 4) # Adjust centers as needed
X$phoenotype.class <- kmeans(risk.matrix, centers = 4)$cluster # (assuming 4 phenotypes ?????)


# profiling discovered phenotypes
str(X[, c("age", "is_ventilated", "comorbidity_level")])

aggregate(
  cbind(age, is_ventilated, comorbidity_level) ~ phoenotype.class, 
  data = X, FUN = mean
)


pred_prob <- bart.model$prob.train.mean
cat("\n=== Predicted probability summary ===\n")
print(summary(pred_prob))

# Classification performance (at 0.5 threshold)
pred_class <- ifelse(pred_prob > 0.5, 1, 0)
conf_matrix <- table(Predicted = pred_class, Actual = y)



