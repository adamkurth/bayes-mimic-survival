#----------------------------------------------------------------
# mimic_plots.R
# Comprehensive plotting script — loads ONLY from checkpoints
#----------------------------------------------------------------
# This script lives in output/mimic/ alongside the checkpoint files
# and produces ALL publication-quality figures for the Bayesian
# discrete-time survival analysis paper.
#
# It loads:
#   1. mimic_analysis_checkpoint.RData  — MAR draws, calibration, survival
#   2. mimic_mnar_results.RData         — MNAR posterior draws + gamma
#
# NO Stan compilation, NO sampling, NO ADVI.
#
# Outputs (saved to plots/mimic/ under project root):
#   Section A — Primary model diagnostics:
#     - mimic_baseline_hazard.pdf       (smooth line + 95% CrI)
#     - mimic_calibration.pdf
#     - mimic_survival.pdf              (spaghetti posterior draws: KM + MAR + MNAR)
#
#   Section B — MNAR sensitivity:
#     - mimic_mnar_gamma_posterior.pdf
#     - mimic_mnar_3panel_posteriors.pdf
#     - mimic_mnar_3panel_all.pdf
#     - mimic_mnar_emergency.pdf
#
# EXCLUDED by design:
#     - mimic_mnar_shift.pdf            (uninformative)
#     - mimic_overview.pdf              (old fixed-sensitivity 6-panel)
#     - mimic_sensitivity_all.pdf       (fixed-gamma grid, replaced by MNAR)
#     - mimic_forest_modelA.pdf         (redundant with 3-panel posteriors)
#     - mimic_AvsB_scatter.pdf          (not essential)
#     - mimic_mnar_forest.pdf           (redundant with 3-panel posteriors)
#     - mimic_mnar_scatter.pdf          (redundant with 3-panel posteriors)
#     - mimic_mnar_prior_sensitivity.pdf (uninformative)
#----------------------------------------------------------------

rm(list = ls())

library(ggplot2)
library(posterior)   # for subset_draws on cmdstanr draws objects

# --- Paths ---
# This script lives in output/mimic/.  Detect project root from that.
script_dir   <- "~/Documents/RStudio/bayes-mimic-survival/output/mimic"
project_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)

out_dir  <- script_dir                                       # checkpoint RData lives here
plot_dir <- file.path(project_root, "plots", "mimic")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cat(sprintf("  script_dir : %s\n", out_dir))
cat(sprintf("  plot_dir   : %s\n", plot_dir))


# =============================================================================
# 1. LOAD CHECKPOINTS
# =============================================================================
cat("\n=== Loading saved checkpoints ===\n")

# --- MNAR results (small, ~871K) ---
mnar_path <- file.path(out_dir, "mimic_mnar_results.RData")
stopifnot(file.exists(mnar_path))
load(mnar_path)  # -> mnar_results, beta_names, beta_A
cat(sprintf("  MNAR results: %d sigma_gamma fits (%s)\n",
            length(mnar_results), paste(names(mnar_results), collapse = ", ")))

# --- Main checkpoint (large, ~5.4GB) ---
chk_path <- file.path(out_dir, "mimic_analysis_checkpoint.RData")
stopifnot(file.exists(chk_path))
cat("  Loading main checkpoint (may take a moment for 5.4GB)...\n")
load(chk_path)
# Available: dat, pp_data, patients_df, draws_A, draws_B,
#            cmp_tbl, cal_summary,
#            [surv_merged, alpha_hat — if checkpoint was saved with recent code]
#            STAN_CHAINS, STAN_ITER, STAN_WARMUP
cat("  Checkpoint loaded.\n")

# Reconstruct beta_A / beta_B from comparison table (safe fallback)
beta_A <- cmp_tbl$Model_A
beta_B <- cmp_tbl$Model_B
beta_names <- colnames(dat$X)
names(beta_A) <- names(beta_B) <- beta_names
P <- dat$P
K <- dat$K


# =============================================================================
# 2. EXTRACT POSTERIOR DRAWS FROM draws_A
# =============================================================================
cat("\n[draws] Extracting posterior draws from draws_A...\n")

# draws_A is a cmdstanr draws_array.  Use posterior::as_draws_matrix
# then subset by variable name.  The "beta" prefix matches beta[1]..beta[P].
draws_mat <- as_draws_matrix(draws_A)

# --- Beta draws (P covariates) ---
beta_col_idx <- grep("^beta\\[", colnames(draws_mat))
beta_draws_mar <- as.matrix(draws_mat[, beta_col_idx])
# Sort columns numerically: beta[1], beta[2], ..., beta[P]
col_nums <- as.integer(gsub("beta\\[(\\d+)\\]", "\\1", colnames(beta_draws_mar)))
beta_draws_mar <- beta_draws_mar[, order(col_nums)]
colnames(beta_draws_mar) <- beta_names
cat(sprintf("  Beta draws: %d samples x %d covariates\n",
            nrow(beta_draws_mar), ncol(beta_draws_mar)))
stopifnot(ncol(beta_draws_mar) == P)

# --- Alpha draws (K intervals) --- for baseline hazard
alpha_col_idx <- grep("^alpha\\[", colnames(draws_mat))
alpha_draws_mat <- as.matrix(draws_mat[, alpha_col_idx])
col_nums_a <- as.integer(gsub("alpha\\[(\\d+)\\]", "\\1", colnames(alpha_draws_mat)))
alpha_draws_mat <- alpha_draws_mat[, order(col_nums_a)]
cat(sprintf("  Alpha draws: %d samples x %d intervals\n",
            nrow(alpha_draws_mat), ncol(alpha_draws_mat)))

# Reconstruct alpha_hat (posterior mean of alpha) if not in checkpoint
if (!exists("alpha_hat") || is.null(alpha_hat)) {
  cat("  Reconstructing alpha_hat from draws...\n")
  alpha_hat <- colMeans(alpha_draws_mat)
}

# Reconstruct surv_merged if not in checkpoint
if (!exists("surv_merged") || is.null(surv_merged)) {
  cat("  Reconstructing survival curves from draws...\n")
  library(survival)
  
  beta_hat <- beta_A
  eta_mat  <- outer(rep(1, dat$N), alpha_hat) + as.vector(dat$X %*% beta_hat)
  h_mat    <- plogis(eta_mat)
  log_surv_mat <- t(apply(log(1 - h_mat), 1, cumsum))
  surv_mat     <- exp(log_surv_mat)
  
  surv_mean <- colMeans(surv_mat)
  surv_lo   <- apply(surv_mat, 2, quantile, 0.025)
  surv_hi   <- apply(surv_mat, 2, quantile, 0.975)
  
  surv_df <- data.frame(
    time = (1:K) * dat$delta_width,
    model_surv = surv_mean, model_lo = surv_lo, model_hi = surv_hi
  )
  
  km_fit <- survfit(
    Surv(y_obs * dat$delta_width, delta) ~ 1,
    data = data.frame(y_obs = dat$y_obs, delta = dat$delta)
  )
  km_times   <- (1:K) * dat$delta_width
  km_summary <- summary(km_fit, times = km_times)
  km_df <- data.frame(time = km_summary$time, km_surv = km_summary$surv,
                      km_lo = km_summary$lower, km_hi = km_summary$upper)
  surv_merged <- merge(surv_df, km_df, by = "time", all.x = TRUE)
  
  rm(eta_mat, h_mat, log_surv_mat, surv_mat, surv_df, km_df, km_fit)
  cat("  Survival curves reconstructed.\n")
}

# Reconstruct cal_summary if missing
if (!exists("cal_summary") || is.null(cal_summary)) {
  cat("  Reconstructing calibration data from draws...\n")
  beta_hat <- beta_A
  eta_mat  <- outer(rep(1, dat$N), alpha_hat) + as.vector(dat$X %*% beta_hat)
  h_mat    <- plogis(eta_mat)
  log_surv_mat <- t(apply(log(1 - h_mat), 1, cumsum))
  surv_mat     <- exp(log_surv_mat)
  
  Ki_vec    <- pmin(dat$y_obs, K)
  pred_mort <- 1 - surv_mat[cbind(1:dat$N, Ki_vec)]
  obs_mort  <- dat$delta
  
  risk_decile <- cut(pred_mort,
                     breaks = quantile(pred_mort, probs = seq(0, 1, 0.1)),
                     include.lowest = TRUE, labels = 1:10)
  cal_df <- data.frame(decile = as.numeric(risk_decile),
                       predicted = pred_mort, observed = obs_mort)
  cal_summary <- aggregate(cbind(predicted, observed) ~ decile, data = cal_df, FUN = mean)
  rm(eta_mat, h_mat, log_surv_mat, surv_mat, cal_df)
  cat("  Calibration data reconstructed.\n")
}

# Free the massive draws objects to reclaim RAM
rm(draws_A, draws_B, draws_mat, pp_data, patients_df)
gc()
cat("  Freed checkpoint RAM.\n")

# Common settings
SIGMA_GAMMA_GRID <- sort(as.numeric(names(mnar_results)))
cols3 <- c("#2171b5", "#cb181d", "#238b45")
mar_mean <- colMeans(beta_draws_mar)
names(mar_mean) <- beta_names


# #############################################################################
#                       SECTION A - PRIMARY MODEL PLOTS
# #############################################################################


# =============================================================================
# A1. BASELINE HAZARD - h(t) with 95% CrI
# =============================================================================
cat("\n[A1] Baseline hazard plot...\n")

h_hat <- plogis(alpha_hat)
h_draws <- plogis(alpha_draws_mat)
h_lo <- apply(h_draws, 2, quantile, 0.025)
h_hi <- apply(h_draws, 2, quantile, 0.975)

hazard_df <- data.frame(day = 1:K, mean = h_hat, lo = h_lo, hi = h_hi)

p_hazard <- ggplot(hazard_df, aes(x = day, y = mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  labs(x = "ICU Day",
       y = expression("Baseline Hazard " * h[0](t)),
       title = "Estimated Baseline Hazard (Model A)",
       subtitle = "Posterior mean with 95% credible interval") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_baseline_hazard.pdf"), p_hazard, width = 8, height = 5)
cat("  Saved: mimic_baseline_hazard.pdf\n")

# Free alpha/hazard draws — A3 spaghetti uses alpha_hat (vector) not draws
rm(alpha_draws_mat, h_draws)


# =============================================================================
# A2. RISK-DECILE CALIBRATION (Model A)
# =============================================================================
cat("[A2] Calibration plot...\n")

p_cal <- ggplot(cal_summary, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.7) +
  geom_point(size = 3.5, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_text(aes(label = decile), vjust = -1, size = 3, color = "grey30") +
  labs(x = "Predicted Mortality (Decile Mean)",
       y = "Observed Mortality (Decile Mean)",
       title = "Risk-Decile Calibration (Model A)",
       subtitle = "Numbers indicate risk decile (1 = lowest, 10 = highest)") +
  coord_equal(xlim = c(0, max(cal_summary$predicted, cal_summary$observed) * 1.1),
              ylim = c(0, max(cal_summary$predicted, cal_summary$observed) * 1.1)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_calibration.pdf"), p_cal, width = 6, height = 6)
cat("  Saved: mimic_calibration.pdf\n")


# =============================================================================
# A3. POPULATION-AVERAGED SURVIVAL — SPAGHETTI POSTERIOR DRAWS
# =============================================================================
# Following hw3 Problem 3(d): overlay N_DRAWS transparent posterior survival
# curves behind a bold posterior-mean line for each model. This conveys
# posterior uncertainty and makes each MNAR specification visually distinct
# (a cloud of draws can't hide behind another line).
cat("[A3] Spaghetti survival curves (MAR + MNAR overlays)...\n")

N_DRAWS <- 200   # posterior samples to draw (cf. hw3 used 300)
set.seed(42)

t_axis <- surv_merged$time

# --- Subsample patients for performance ---
# dat$X is the patient-level design matrix (N x P).  pop_surv computes an
# N x K matrix per draw, so 200 draws x 4 models x 57K patients is slow.
# Subsample N_SUB patients for the spaghetti draws.
N_SUB <- 5000
set.seed(123)
sub_ids <- sample(dat$N, min(N_SUB, dat$N))
X_sub <- dat$X[sub_ids, , drop = FALSE]
cat(sprintf("  Subsampled %d of %d patients for spaghetti draws\n",
            nrow(X_sub), dat$N))

# CRITICAL: Sort mnar_keys for deterministic color assignment
mnar_keys <- as.character(sort(as.numeric(names(mnar_results))))
mnar_model_names <- paste0("MNAR (sg=", mnar_keys, ")")

cat(sprintf("  mnar_keys: %s\n", paste(mnar_keys, collapse = ", ")))
cat(sprintf("  N_DRAWS = %d posterior samples per model\n", N_DRAWS))

# ---- Helper: compute population-averaged survival for one beta vector ----
# Returns length-K vector of S(t) = (1/N) sum_i prod_{s=1}^{t} [1-h_i(s)]
pop_surv <- function(beta_vec, X_sub = NULL) {
  Xmat <- if (!is.null(X_sub)) X_sub else dat$X
  Xb  <- as.numeric(Xmat %*% as.numeric(beta_vec))
  eta <- outer(Xb, alpha_hat, FUN = "+")
  h   <- plogis(eta)
  log_S <- t(apply(log(1 - h), 1, cumsum))
  colMeans(exp(log_S))
}

# ---- 1. MAR spaghetti draws ----
cat("  Computing MAR spaghetti draws...\n")
n_total_mar <- nrow(beta_draws_mar)
draw_idx_mar <- sample(n_total_mar, min(N_DRAWS, n_total_mar))

mar_spaghetti_list <- list()
for (d in seq_along(draw_idx_mar)) {
  beta_d <- beta_draws_mar[draw_idx_mar[d], ]
  S_d <- pop_surv(beta_d, X_sub)
  mar_spaghetti_list[[d]] <- data.frame(
    time = t_axis, surv = S_d, draw = d,
    stringsAsFactors = FALSE
  )
}
mar_spaghetti <- do.call(rbind, mar_spaghetti_list)

# MAR posterior-mean curve
mar_mean_surv <- pop_surv(mar_mean, X_sub)

cat(sprintf("  MAR mean S(1)=%.4f, S(10)=%.4f, S(30)=%.4f\n",
            mar_mean_surv[1], mar_mean_surv[min(10,K)], mar_mean_surv[K]))

# ---- 2. MNAR spaghetti draws per sigma_gamma ----
mnar_spaghetti_list <- list()
mnar_mean_surv_list <- list()

for (i in seq_along(mnar_keys)) {
  sg <- mnar_keys[i]
  beta_mat <- mnar_results[[sg]]$beta   # n_draws x P
  n_avail  <- nrow(beta_mat)
  draw_idx <- sample(n_avail, min(N_DRAWS, n_avail))
  
  cat(sprintf("  Computing MNAR (sg=%s) spaghetti: %d draws from %d available...\n",
              sg, length(draw_idx), n_avail))
  
  sg_list <- list()
  for (d in seq_along(draw_idx)) {
    beta_d <- beta_mat[draw_idx[d], ]
    S_d <- pop_surv(beta_d, X_sub)
    sg_list[[d]] <- data.frame(
      time = t_axis, surv = S_d, draw = d,
      stringsAsFactors = FALSE
    )
  }
  mnar_spaghetti_list[[sg]] <- do.call(rbind, sg_list)
  
  # Posterior-mean curve
  b_mean <- colMeans(beta_mat)
  mnar_mean_surv_list[[sg]] <- pop_surv(b_mean, X_sub)
  
  cat(sprintf("    Mean S(1)=%.4f, S(10)=%.4f, S(30)=%.4f\n",
              mnar_mean_surv_list[[sg]][1],
              mnar_mean_surv_list[[sg]][min(10,K)],
              mnar_mean_surv_list[[sg]][K]))
}

# ---- 3. Build the plot layer by layer (hw3 matlines approach in ggplot) ----
cat("  Building spaghetti plot...\n")

# Legend labels — store as parseable strings, then parse in scale_*_manual
surv_labels <- c(
  "KM"  = "'Kaplan-Meier'",
  "MAR" = "'MAR (Model A)'"
)
for (sg in mnar_keys) {
  key <- paste0("MNAR_", sg)
  surv_labels[key] <- paste0("'MNAR (' * sigma[gamma] == ", sg, " * ')'")
}

p_surv <- ggplot() +
  
  # --- Layer 0: KM confidence band (background) ---
  geom_ribbon(data = data.frame(time = t_axis,
                                lo = surv_merged$km_lo,
                                hi = surv_merged$km_hi),
              aes(x = time, ymin = lo, ymax = hi),
              fill = "grey70", alpha = 0.15) +
  
  # --- Layer 1: MAR spaghetti draws (grey, very transparent) ---
  geom_line(data = mar_spaghetti,
            aes(x = time, y = surv, group = draw),
            color = "grey50", alpha = 0.04, linewidth = 0.3)

# --- Layer 2: MNAR spaghetti draws (colored, very transparent) ---
for (i in seq_along(mnar_keys)) {
  sg <- mnar_keys[i]
  p_surv <- p_surv +
    geom_line(data = mnar_spaghetti_list[[sg]],
              aes(x = time, y = surv, group = draw),
              color = cols3[i], alpha = 0.06, linewidth = 0.3)
}

# --- Layer 3: Bold posterior-mean lines ---
# KM (dashed black)
p_surv <- p_surv +
  geom_line(data = data.frame(time = t_axis, surv = surv_merged$km_surv),
            aes(x = time, y = surv, color = "KM", linetype = "KM"),
            linewidth = 0.9) +
  # MAR (grey solid)
  geom_line(data = data.frame(time = t_axis, surv = mar_mean_surv),
            aes(x = time, y = surv, color = "MAR", linetype = "MAR"),
            linewidth = 1.1)

# MNAR posterior-mean lines (colored, solid, bold)
for (i in seq_along(mnar_keys)) {
  sg <- mnar_keys[i]
  key <- paste0("MNAR_", sg)
  p_surv <- p_surv +
    geom_line(data = data.frame(time = t_axis, surv = mnar_mean_surv_list[[sg]]),
              aes(x = time, y = surv, color = !!key, linetype = !!key),
              linewidth = 1.1)
}

# --- Scales ---
color_vals <- c("KM" = "black", "MAR" = "grey40")
lty_vals   <- c("KM" = "dashed", "MAR" = "solid")
for (i in seq_along(mnar_keys)) {
  key <- paste0("MNAR_", mnar_keys[i])
  color_vals[key] <- cols3[i]
  lty_vals[key]   <- "solid"
}

p_surv <- p_surv +
  scale_color_manual(values = color_vals,
                     labels = parse(text = surv_labels)) +
  scale_linetype_manual(values = lty_vals,
                        labels = parse(text = surv_labels)) +
  labs(x = "Days Since ICU Admission",
       y = "S(t)",
       title = "Population-Averaged Survival Under MAR and MNAR",
       subtitle = paste0("Transparent bands = ", N_DRAWS,
                         " posterior draws per model; bold lines = posterior mean."),
       color = "", linetype = "") +
  coord_cartesian(ylim = c(
    min(c(mar_spaghetti$surv,
          sapply(mnar_spaghetti_list, function(d) min(d$surv))),
        na.rm = TRUE) * 0.98,
    1.0)) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_survival.pdf"), p_surv, width = 9, height = 6)
cat("  Saved: mimic_survival.pdf\n")

# Free large objects from A3
rm(mar_spaghetti, mar_spaghetti_list,
   mnar_spaghetti_list, mnar_mean_surv_list, X_sub)


# #############################################################################
#                       SECTION B - MNAR SENSITIVITY PLOTS
# #############################################################################

# Reference MNAR fit (sigma_gamma = 1.0)
ref_key <- "1"
beta_mnar_ref <- colMeans(mnar_results[[ref_key]]$beta)
shifts <- beta_mnar_ref - beta_A
names(shifts) <- beta_names

cat("\n=== Beta shifts (MAR -> MNAR, sigma_gamma = 1.0) ===\n")
shift_order <- order(abs(shifts), decreasing = TRUE)
for (j in shift_order) {
  cat(sprintf("  %-20s  %+.3f  (%.3f -> %.3f)\n",
              beta_names[j], shifts[j], beta_A[j], beta_mnar_ref[j]))
}

top_covs <- beta_names[shift_order[1:6]]
cat(sprintf("\n  Top 6: %s\n", paste(top_covs, collapse = ", ")))

# Check that MNAR beta columns match beta_names
mnar_P <- ncol(mnar_results[[ref_key]]$beta)
cat(sprintf("  MNAR beta columns: %d | beta_names: %d\n", mnar_P, length(beta_names)))
stopifnot(mnar_P == length(beta_names))


# =============================================================================
# B1. GAMMA POSTERIOR DENSITIES WITH PRIOR OVERLAYS
# =============================================================================
cat("\n[B1] Gamma posterior densities...\n")

# --- 1. Build Data Frame for Posterior Draws ---
gamma_df <- do.call(rbind, lapply(names(mnar_results), function(sg) {
  data.frame(
    sigma_gamma = paste0("sigma[gamma] == ", sg),
    gamma       = as.numeric(mnar_results[[sg]]$gamma_draws),
    stringsAsFactors = FALSE
  )
}))
gamma_df$sigma_gamma <- factor(gamma_df$sigma_gamma,
                               levels = paste0("sigma[gamma] == ", SIGMA_GAMMA_GRID))

# --- 2. Subsample for Posterior Draw Ticks (Rug Plot) ---
set.seed(42)
rug_gamma_df <- do.call(rbind, lapply(split(gamma_df, gamma_df$sigma_gamma), function(d) {
  d[sample(nrow(d), min(1000, nrow(d))), ]
}))

# --- 3. Build Data Frame for Theoretical Priors ---
x_grid <- seq(min(gamma_df$gamma) - 1, 2, length.out = 500)
prior_df <- do.call(rbind, lapply(SIGMA_GAMMA_GRID, function(sg) {
  data.frame(
    sigma_gamma = paste0("sigma[gamma] == ", sg),
    x           = x_grid,
    density     = dnorm(x_grid, 0, sg),
    stringsAsFactors = FALSE
  )
}))
prior_df$sigma_gamma <- factor(prior_df$sigma_gamma,
                               levels = levels(gamma_df$sigma_gamma))

# --- 4. Main Density Plot ---
p_gamma <- ggplot(gamma_df, aes(x = gamma)) +
  # Individual posterior draws as low-opacity ticks on the x-axis
  geom_rug(data = rug_gamma_df, aes(color = sigma_gamma),
           alpha = 0.05, length = unit(0.04, "npc"), show.legend = FALSE) +
  
  # Posterior density curves
  geom_density(aes(fill = sigma_gamma, color = sigma_gamma),
               alpha = 0.20, linewidth = 1.1) +
  
  # Prior density curves (dashed)
  geom_line(data = prior_df, aes(x = x, y = density, color = sigma_gamma),
            linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
  
  # Zero reference line for MAR
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey40", linewidth = 0.8) +
  annotate("text", x = 0.2, y = Inf, label = "MAR (gamma == 0)", vjust = 2.0, parse = TRUE,
           color = "grey30", size = 3.5, fontface = "italic") +
  
  # Apply colors and parse LaTeX labels
  scale_fill_manual(name = "Prior Width", values = cols3,
                    labels = parse(text = levels(gamma_df$sigma_gamma))) +
  scale_color_manual(name = "Prior Width", values = cols3,
                     labels = parse(text = levels(gamma_df$sigma_gamma))) +
  
  # FORCE A SINGLE LEGEND: Hide color legend, use fill legend with overridden aesthetics
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.4, color = cols3))) +
  
  labs(x = expression(gamma ~ "(selection parameter)"),
       y = "Density",
       title = expression("Posterior of " * gamma * " Under Three Prior Widths"),
       subtitle = "Solid = posterior; Dashed = prior. Bottom ticks = individual posterior draws.") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold"),
    legend.position   = c(0.15, 0.75),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80", linewidth = 0.5),
    legend.title      = element_text(face = "bold", size = 10),
    legend.key.size   = unit(0.5, "cm"),
    legend.margin     = margin(t = 4, r = 6, b = 4, l = 6),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(plot_dir, "mimic_mnar_gamma_posterior.pdf"), p_gamma, width = 8, height = 5)
cat("  Saved: mimic_mnar_gamma_posterior.pdf\n")

# =============================================================================
# B2. 3-PANEL POSTERIOR DENSITIES - MAR vs MNAR (Top 6 Covariates)
# =============================================================================
cat("[B2] 3-panel posterior density comparison (top 6)...\n")

# Helper: get column index by covariate name (safe lookup)
cov_idx <- function(cov) which(beta_names == cov)

mar_long <- do.call(rbind, lapply(top_covs, function(cov) {
  j <- cov_idx(cov)
  data.frame(parameter = cov, value = as.numeric(beta_draws_mar[, j]),
             model = "MAR (Model A)", stringsAsFactors = FALSE)
}))

panel_list <- lapply(names(mnar_results), function(sg) {
  b_draws <- mnar_results[[sg]]$beta
  mnar_long <- do.call(rbind, lapply(top_covs, function(cov) {
    j <- cov_idx(cov)
    data.frame(parameter = cov, value = as.numeric(b_draws[, j]),
               model = "MNAR", stringsAsFactors = FALSE)
  }))
  panel_df <- rbind(mar_long, mnar_long)
  panel_df$sigma_panel <- paste0("sigma[gamma] == ", sg)
  return(panel_df)
})

full_df <- do.call(rbind, panel_list)

cov_order <- top_covs[order(abs(mar_mean[top_covs]), decreasing = TRUE)]
full_df$parameter <- factor(full_df$parameter, levels = rev(cov_order))
full_df$sigma_panel <- factor(full_df$sigma_panel,
                              levels = paste0("sigma[gamma] == ", SIGMA_GAMMA_GRID))

p_3panel <- ggplot(full_df, aes(x = value, y = parameter,
                                fill = model, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_violin(position = position_identity(),
              alpha = 0.3, scale = "width",
              draw_quantiles = c(0.025, 0.5, 0.975),
              linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 2.5,
               position = position_dodge(width = 0), show.legend = FALSE) +
  facet_wrap(~ sigma_panel, ncol = 3, labeller = label_parsed) +
  scale_fill_manual(values = c("MAR (Model A)" = "#2171b5",
                               "MNAR"          = "#cb181d")) +
  scale_color_manual(values = c("MAR (Model A)" = "#08519c",
                                "MNAR"          = "#a50f15")) +
  labs(x = "Posterior Value (Log-Odds Scale)",
       y = "",
       title = "Posterior Densities: MAR vs. MNAR for Most Affected Covariates",
       subtitle = "Shaded violins = full posterior. Inner lines = 2.5%, 50%, 97.5% quantiles. Points = mean.",
       fill = "Model", color = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "bottom",
        plot.title       = element_text(face = "bold"),
        strip.text       = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor   = element_blank(),
        axis.text.y = element_text(face = "bold", size = 10))

ggsave(file.path(plot_dir, "mimic_mnar_3panel_posteriors.pdf"),
       p_3panel, width = 14, height = 7)
cat("  Saved: mimic_mnar_3panel_posteriors.pdf\n")


# =============================================================================
# B3. 3-PANEL POSTERIOR DENSITIES - ALL 17 COVARIATES
# =============================================================================
cat("[B3] Full 17-covariate 3-panel density...\n")

mar_full_long <- do.call(rbind, lapply(seq_along(beta_names), function(j) {
  data.frame(parameter = beta_names[j], value = as.numeric(beta_draws_mar[, j]),
             model = "MAR (Model A)", stringsAsFactors = FALSE)
}))

panel_full_list <- lapply(names(mnar_results), function(sg) {
  b_draws <- mnar_results[[sg]]$beta
  mnar_long <- do.call(rbind, lapply(seq_along(beta_names), function(j) {
    data.frame(parameter = beta_names[j], value = as.numeric(b_draws[, j]),
               model = "MNAR", stringsAsFactors = FALSE)
  }))
  panel_df <- rbind(mar_full_long, mnar_long)
  panel_df$sigma_panel <- paste0("sigma[gamma] == ", sg)
  return(panel_df)
})

full_all_df <- do.call(rbind, panel_full_list)

all_order <- beta_names[order(mar_mean[beta_names])]
full_all_df$parameter <- factor(full_all_df$parameter, levels = all_order)
full_all_df$sigma_panel <- factor(full_all_df$sigma_panel,
                                  levels = paste0("sigma[gamma] == ", SIGMA_GAMMA_GRID))

p_3panel_all <- ggplot(full_all_df,
                       aes(x = value, y = parameter,
                           fill = model, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_violin(position = position_identity(),
              alpha = 0.25, scale = "width",
              draw_quantiles = c(0.5), linewidth = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 1.8,
               position = position_dodge(width = 0), show.legend = FALSE) +
  facet_wrap(~ sigma_panel, ncol = 3, labeller = label_parsed) +
  scale_fill_manual(values = c("MAR (Model A)" = "#2171b5",
                               "MNAR"          = "#cb181d")) +
  scale_color_manual(values = c("MAR (Model A)" = "#08519c",
                                "MNAR"          = "#a50f15")) +
  labs(x = "Posterior Value (Log-Odds)",
       y = "",
       title = "All 17 Covariate Posteriors: MAR vs. MNAR",
       subtitle = "Vertical median lines and mean points. Wider priors allow greater MNAR departure.",
       fill = "Model", color = "Model") +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        plot.title       = element_text(face = "bold"),
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major.y = element_line(color = "grey95"),
        panel.grid.minor   = element_blank(),
        axis.text.y = element_text(size = 8))

ggsave(file.path(plot_dir, "mimic_mnar_3panel_all.pdf"),
       p_3panel_all, width = 16, height = 10)
cat("  Saved: mimic_mnar_3panel_all.pdf\n")


# =============================================================================
# B4. EMERGENCY COVARIATE DEEP DIVE — MAR vs MNAR across sigma_gamma
# =============================================================================
cat("[B4] Emergency covariate deep dive...\n")

emerg_idx <- which(beta_names == "emergency")

# --- Build density data: MAR + one MNAR per sigma_gamma ---
emerg_mar <- data.frame(
  value = as.numeric(beta_draws_mar[, emerg_idx]),
  model = "MAR (Model A)",
  stringsAsFactors = FALSE
)

emerg_mnar_list <- lapply(names(mnar_results), function(sg) {
  data.frame(
    value = as.numeric(mnar_results[[sg]]$beta[, emerg_idx]),
    model = paste0("MNAR (sigma_gamma = ", sg, ")"),
    stringsAsFactors = FALSE
  )
})

emerg_df <- do.call(rbind, c(list(emerg_mar), emerg_mnar_list))

# Keep simple text for factor levels to ensure color mapping works
emerg_model_levels <- c("MAR (Model A)", paste0("MNAR (sigma_gamma = ", SIGMA_GAMMA_GRID, ")"))
emerg_df$model <- factor(emerg_df$model, levels = emerg_model_levels)

# --- Subsample for Posterior Draw Ticks (Rug Plot) ---
set.seed(42)
rug_df <- do.call(rbind, lapply(split(emerg_df, emerg_df$model), function(d) {
  d[sample(nrow(d), min(1000, nrow(d))), ]
}))

# --- Compute summary stats for annotation ---
emerg_stats <- data.frame(
  model     = levels(emerg_df$model),
  post_mean = NA_real_,
  q025      = NA_real_,
  q975      = NA_real_,
  stringsAsFactors = FALSE
)
# MAR
emerg_stats$post_mean[1] <- mean(emerg_mar$value)
emerg_stats$q025[1]      <- quantile(emerg_mar$value, 0.025)
emerg_stats$q975[1]      <- quantile(emerg_mar$value, 0.975)

# MNAR per sigma_gamma
for (i in seq_along(SIGMA_GAMMA_GRID)) {
  sg_key <- as.character(SIGMA_GAMMA_GRID[i])
  draws_i <- as.numeric(mnar_results[[sg_key]]$beta[, emerg_idx])
  emerg_stats$post_mean[i + 1] <- mean(draws_i)
  emerg_stats$q025[i + 1]      <- quantile(draws_i, 0.025)
  emerg_stats$q975[i + 1]      <- quantile(draws_i, 0.975)
}

# --- MATCH COLORS TO B1 ---
emerg_cols <- c(
  "MAR (Model A)" = "grey30",
  setNames(
    cols3,
    paste0("MNAR (sigma_gamma = ", SIGMA_GAMMA_GRID, ")")
  )
)

# Parse labels for proper LaTeX sigma rendering in the legend
parsed_labels <- parse(text = c(
  "'MAR (Model A)'",
  paste0("'MNAR (' * sigma[gamma] == ", SIGMA_GAMMA_GRID, " * ')'")
))

# --- Main density plot ---
p_emerg <- ggplot(emerg_df, aes(x = value, fill = model, color = model)) +
  # 1. Plot individual posterior draws as low-opacity ticks on the x-axis
  geom_rug(data = rug_df, alpha = 0.05, length = unit(0.04, "npc"), show.legend = FALSE) +
  
  # 2. Main density curves
  geom_density(alpha = 0.15, linewidth = 1.0) +
  
  # 3. Vertical lines at posterior means
  geom_vline(data = emerg_stats, aes(xintercept = post_mean, color = model),
             linetype = "solid", linewidth = 0.7, show.legend = FALSE) +
  
  # 4. Zero reference line
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey60", linewidth = 0.6) +
  
  # Annotate the MAR mean
  annotate("text",
           x = emerg_stats$post_mean[1], y = Inf,
           label = sprintf("MAR: %.3f", emerg_stats$post_mean[1]),
           vjust = 2.0, hjust = -0.1, size = 3.3, color = "grey30", fontface = "bold") +
  
  # Annotate the widest MNAR mean
  annotate("text",
           x = emerg_stats$post_mean[nrow(emerg_stats)], y = Inf,
           label = sprintf("MNAR (2.0): %.3f", emerg_stats$post_mean[nrow(emerg_stats)]),
           vjust = 3.5, hjust = -0.1, size = 3.3, color = cols3[3], fontface = "bold") +
  
  # Annotate the shift (Delta Beta)
  annotate("segment",
           x = emerg_stats$post_mean[1], xend = emerg_stats$post_mean[nrow(emerg_stats)],
           y = 0.08, yend = 0.08,
           arrow = arrow(length = unit(0.15, "cm"), ends = "both"),
           color = "grey30", linewidth = 0.5) +
  annotate("text",
           x = mean(c(emerg_stats$post_mean[1], emerg_stats$post_mean[nrow(emerg_stats)])),
           y = 0.08,
           label = sprintf("Delta*beta == %.3f",
                           emerg_stats$post_mean[nrow(emerg_stats)] - emerg_stats$post_mean[1]),
           parse = TRUE, vjust = -0.5, size = 3.2, color = "grey20") +
  
  # Apply updated colors and parsed labels
  scale_fill_manual(values = emerg_cols, labels = parsed_labels) +
  scale_color_manual(values = emerg_cols, labels = parsed_labels) +
  
  labs(x = expression(beta[emergency] ~ "(log-odds scale)"),
       y = "Posterior Density",
       title = expression("Emergency Admission Effect: MAR vs. MNAR Across " * sigma[gamma]),
       subtitle = paste0("Informative censoring attenuates the emergency effect from ",
                         sprintf("%.3f to %.3f", emerg_stats$post_mean[1],
                                 emerg_stats$post_mean[nrow(emerg_stats)]),
                         ".\nVertical bottom ticks represent individual MCMC posterior draws."),
       fill = "Model", color = "Model") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold"),
    legend.position   = c(0.15, 0.75),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80", linewidth = 0.5),
    legend.title      = element_blank(),
    legend.key.size   = unit(0.4, "cm"),
    legend.margin     = margin(t = 3, r = 5, b = 3, l = 5),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(plot_dir, "mimic_mnar_emergency.pdf"), p_emerg, width = 9, height = 6)
cat("  Saved: mimic_mnar_emergency.pdf\n")

# --- Supplementary: CrI comparison table printed to console ---
cat("\n  Emergency posterior summary:\n")
cat(sprintf("  %-30s  Mean     [95%% CrI]\n", "Model"))
for (r in seq_len(nrow(emerg_stats))) {
  cat(sprintf("  %-30s  %+.3f   [%.3f, %.3f]\n",
              emerg_stats$model[r],
              emerg_stats$post_mean[r],
              emerg_stats$q025[r], emerg_stats$q975[r]))
}
cat("\n")


# =============================================================================
# DONE
# =============================================================================
cat("\n=== All plots complete ===\n")
cat("  No models were refit. All data loaded from checkpoints.\n")
cat("\n  Section A (primary model):\n")
cat("    mimic_baseline_hazard.pdf       (smooth line + 95% CrI)\n")
cat("    mimic_calibration.pdf\n")
cat("    mimic_survival.pdf              (spaghetti posterior draws)\n")
cat("\n  Section B (MNAR sensitivity):\n")
cat("    mimic_mnar_gamma_posterior.pdf\n")
cat("    mimic_mnar_3panel_posteriors.pdf\n")
cat("    mimic_mnar_3panel_all.pdf\n")
cat("    mimic_mnar_emergency.pdf\n")