####### Bayesian Additive Regression Trees #######
# Updated: 2026-03-25
rm(list = ls()); source("load.r")
library(BART);library(dplyr); library(dirichletprocess); library(mclust)
rds.path <- file.path(getwd(), "data/rds")

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

# extract M x N matrix of predicted probabilities (M = post. draws, N = patients)
prob.mat <- bart.model$prob.train # M x N
risk.means <- colMeans(prob.mat) # N-length vector of mean predicted probabilities
risk.sds <- apply(prob.mat, 2, sd) # N-length vector of uncertainty (SD) for each patient


# ==============================================================================
# Approach 1: 1D Dirichlet Process on Posterior Means (risk stratification)
# ==============================================================================
cat("\n=== Fitting 1D DP on Posterior Means ===\n")
risk.means.sc <- scale(risk.means)[, 1] # Standardize for DP clustering

dp.1d <- DirichletProcessGaussian(y = risk.means.sc)
dp.1d <- Fit(dp.1d, its = 1000)

X$dp.1d.phenotype <- as.factor(dp.1d$clusterLabels)
table(X$dp.1d.phenotype)
cat("Phoenotypes discovered by 1D DP (risk only):", length(unique(X$dp.1d.phenotype)), "\n")


# ==============================================================================
# Approach 2: MVN DP on Full Posterior Summary (uncertainty profiling)
# ==============================================================================
cat("\n=== Fitting MVN DP (Mean & SD) ===\n")

post.sum.matrix <- cbind(scale(risk.means)[, 1], scale(risk.sds)[, 1]) # Standardize mean & SD for DP clustering

# fit MVN DP mixture
# dp.2d <- DirichletProcessMvnormal(y = post.sum.matrix)
# dp.2d <- Fit(dp.2d, its = 1000)

# X$dp.2d.phenotype <- as.factor(dp.2d$clusterLabels)
# table(X$dp.2d.phoenotype)
# cat("Phoenotypes discovered by MVN DP (mean & uncertainty):", length(unique(X$dp.2d.phoenotype)), "\n")

gmm.2d <- Mclust(data=post.sum.matrix) # Optimal number of clusters based on BIC
cat("GMM BIC-optimal clusters:", gmm.2d$G, "\n")

X$gmm.2d.phoenotype <- as.factor(gmm.2d$classification)
table(X$gmm.2d.phoenotype)
cat("Phoenotypes discovered by GMM (mean & uncertainty):", length(unique(X$gmm.2d.phoenotype)), "\n")



# ==============================================================================
# EXAMINE RESULTS
# ==============================================================================
library(ggplot2); library(gridExtra)
table(Risk = X$dp.1d.phenotype, Uncertainty = X$gmm.2d.phoenotype)

X.full <- X
X.full$obs.mortality <- y
X.full$pred.risk <- risk.means
X.full$risk.uncertainty <- risk.sds


# COMPARISON 1: VISUALIZING THE PHENOTYPE GEOMETRY
p1 <- ggplot(X.full, aes(x = pred.risk, y = risk.uncertainty, color = dp.1d.phoenotype)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(title = "1D DP (Risk Only):", x = "BART Expected Mortality Risk", y = "BART Posterior SD") +
  theme(legend.position = "right")
p2 <- ggplot(X.full, aes(x = pred.risk, y = risk.uncertainty, color = gmm.2d.phoenotype)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(title = "2D DP (Risk + Uncertainty):", x = "BART Expected Mortality Risk", y = "BART Posterior SD") +
  theme(legend.position = "right")

grid.arrange(p1, p2, ncol = 1, nrow=2)





# Plotting the geometry
p1 <- ggplot(X.full, aes(x = pred.risk, y = risk.uncertainty, color = dp.1d.phenotype)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "1D DP (Blind to Uncertainty)", subtitle = "Merges Stable & Volatile Low Risk Groups",
       x = "Expected Mortality Risk", y = "Posterior SD")

p2 <- ggplot(X.full, aes(x = pred.risk, y = risk.uncertainty, color = gmm.2d.phoenotype)) +  # FIX: Update color variable
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "2D MVN GMM (Risk + Uncertainty)", subtitle = "Successfully Recovers All 3 Latent Groups",
       x = "Expected Mortality Risk", y = "Posterior SD")

grid.arrange(p1, p2, ncol = 1, nrow = 2)







# COMPARISON 2: CLINICAL PROFILING & CALIBRATION
gen.profile <- function(data, cluster_col) {
  data %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(
      N = n(),
      Obs_Mortality_Rate = round(mean(obs.mortality) * 100, 1),
      Mean_BART_Risk = round(mean(pred.risk) * 100, 1),
      Mean_Age = round(mean(age), 1),
      Mean_GCS = round(mean(gcs_mean_val, na.rm = TRUE), 1),
      Pct_Sepsis = round(mean(as.numeric(as.character(icd_sepsis))) * 100, 1),
      Pct_Ventilated = round(mean(as.numeric(as.character(is_ventilated))) * 100, 1),
      Mean_Vasopressor_Hrs = round(mean(vasopressor_duration_hrs, na.rm = TRUE), 1),
      Mean_DRG_Severity = round(mean(drg_severity_score, na.rm = TRUE), 2)
    ) %>%
    arrange(Mean_BART_Risk)
}

cat("\n--- Clinical Profile: 1D DP (10 Clusters) ---\n")
profile.1d <- gen.profile(X.full, "dp.1d.phoenotype")
print(profile.1d, width = Inf) # splitting based on artifcacts of N=1 or 2 people creating "distinct profile"

cat("\n--- Clinical Profile: 2D DP (2 Clusters) ---\n")
profile.2d <- gen.profile(X.full, "dp.2d.phoenotype")
print(profile.2d, width = Inf) # 2 major clusters with


longitudinal_df %>%
  filter(hadm_id == "114867") %>%
  write.csv(file.path(getwd(), "long-data-1.csv"), row.names = FALSE)
longitudinal_df %>%
  filter(hadm_id == "130870") %>%
  write.csv(file.path(getwd(), "long-data-2.csv"), row.names = FALSE)



# ==============================================================================
# SIMULATED DEMONSTRATION 
# ==============================================================================
rm(list = ls()); source("load.r")
library(BART);library(dplyr); library(dirichletprocess)
rds.path <- file.path(getwd(), "data/rds")
data <- readRDS(file.path(rds.path, "X_y.rds"))
X <- data$X; y <- data$y

cat("\n=== Simulating Large Scale BART Posterior Matrix ===\n")

N.per.group <- 1000
M.draws <- 1000

sim.post <- function(n, mean.risk, sd.risk, draws){
  alpha <- ( (1-mean.risk) / sd.risk^2 - 1/mean.risk ) * mean.risk^2
  beta <- alpha * (1/mean.risk - 1)
  # replicate
  replicate(n, rbeta(n=draws,shape1=alpha, shape2=beta))
}

# Group 1: Standard Low Risk
prob.g1 <- sim.post(n=N.per.group, mean.risk=0.15, sd.risk=0.02, draws=M.draws)
# Group 2: Volatile Low Risk (Same mean, high uncertainty)
prob.g2 <- sim.post(n=N.per.group, mean.risk = 0.15, sd.risk = 0.08, draws = M.draws)
# Group 3: High Risk
prob.g3 <- sim.post(n=N.per.group, mean.risk = 0.6, sd.risk = 0.1, draws = M.draws)

prob.mat <- cbind(prob.g1, prob.g2, prob.g3) # M x N matrix (1000 x 3000)
risk.means <- colMeans(prob.mat)
risk.sds <- apply(prob.mat, 2, sd)

sim.df <- data.frame(
  true.latent.class = factor(rep(c("low", "volatile.low", "high.risk"), each = N.per.group)),
  pred.risk = risk.means,
  risk.uncertainty = risk.sds
)


# ==============================================================================
# 2. FIT 1D DP (Risk Only)
# ==============================================================================
cat("\n=== Fitting 1D DP on Posterior Means ===\n")
risk.means.sc <- scale(risk.means)[, 1] 

dp.1d <- DirichletProcessGaussian(y = risk.means.sc)
dp.1d <- Fit(dp.1d, its = 500) # Reduced iterations for simulation speed

sim.df$dp.1d.phenotype <- as.factor(dp.1d$clusterLabels)

# ==============================================================================
# 3. FIT 2D DP (Risk + Uncertainty)
# ==============================================================================
cat("\n=== Fitting MVN DP (Mean & SD) ===\n")
post.sum.matrix <- cbind(scale(risk.means)[, 1], scale(risk.sds)[, 1]) 

# dp.2d <- DirichletProcessMvnormal(y = post.sum.matrix) # very slow
# dp.2d <- Fit(dp.2d, its = 500)

gmm.2d <- Mclust(data=post.sum.matrix) # Optimal number of clusters based on BIC
cat("GMM BIC-optimal clusters:", gmm.2d$G, "\n")

sim.df$gmm.2d.phoenotype <- as.factor(gmm.2d$classification)
table(sim.df$gmm.2d.phoenotype)
cat("Phoenotypes discovered by GMM (mean & uncertainty):", length(unique(sim.df$gmm.2d.phoenotype)), "\n")

sim.df$dp.2d.phenotype <- as.factor(dp.2d$clusterLabels)

# ==============================================================================
# 4. EXAMINE & VISUALIZE RESULTS
# ==============================================================================

cat("\n--- 1D DP Cluster Mapping ---\n")
print(table(True_Class = sim.df$true.latent.class, Discovered_1D = sim.df$dp.1d.phenotype))

cat("\n--- 2D GMM Cluster Mapping ---\n")
# FIX: Reference the GMM phenotype here!
print(table(True_Class = sim.df$true.latent.class, Discovered_2D = sim.df$gmm.2d.phoenotype))

# Plotting the geometry
p1 <- ggplot(sim.df, aes(x = pred.risk, y = risk.uncertainty, color = dp.1d.phenotype)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "1D DP (Blind to Uncertainty)", subtitle = "Merges Stable & Volatile Low Risk Groups",
       x = "Expected Mortality Risk", y = "Posterior SD")

p2 <- ggplot(sim.df, aes(x = pred.risk, y = risk.uncertainty, color = gmm.2d.phoenotype)) +  # FIX: Update color variable
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "2D MVN GMM (Risk + Uncertainty)", subtitle = "Successfully Recovers All 3 Latent Groups",
       x = "Expected Mortality Risk", y = "Posterior SD")

grid.arrange(p1, p2, ncol = 1, nrow = 2)