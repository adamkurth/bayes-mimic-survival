# ============================================================================
# mimic_plots.R — All figures for the paper
# ============================================================================
# Self-contained: loads checkpoints, extracts draws, computes everything,
# and produces the 5 publication figures referenced in final-project.tex:
#
#   Figure 1  mimic_followup.pdf
#   Figure 2  mimic_baseline_hazard.pdf
#   Figure 3  mimic_calibration.pdf
#   Figure 4  mimic_survival.pdf
#   Figure 5  mimic_mnar_gamma_posterior.pdf
#   Figure 6  mimic_mnar_emergency.pdf
#
# Inputs (from mimic_analysis.R):
#   output/mimic/mimic_analysis_checkpoint.RData
#   output/mimic/mimic_mnar_results.RData
# ============================================================================

rm(list = ls())

library(ggplot2)
library(posterior)

# ---- Paths -----------------------------------------------------------------
project_root <- "~/Documents/RStudio/bayes-mimic-survival"
out_dir      <- file.path(project_root, "output/mimic")
plot_dir     <- file.path(project_root, "plots/mimic")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# ============================================================================
# 1. LOAD CHECKPOINTS
# ============================================================================

load(file.path(out_dir, "mimic_mnar_results.RData"))   # → mnar_results, beta_names, beta_A
load(file.path(out_dir, "mimic_analysis_checkpoint.RData"))
# → dat, draws_A, draws_B, cmp_tbl, cal_summary, surv_merged, alpha_hat, ...

# Reconstruct beta vectors from checkpoint
if (!exists("beta_A")) { beta_A <- cmp_tbl$Model_A; names(beta_A) <- colnames(dat$X) }
beta_names <- colnames(dat$X)
P <- dat$P;  K <- dat$K

cat(sprintf("[loaded] N=%d  K=%d  P=%d  MNAR fits=%d\n",
            dat$N, K, P, length(mnar_results)))


# ============================================================================
# 2. EXTRACT POSTERIOR DRAWS
# ============================================================================

draws_mat <- as_draws_matrix(draws_A)

# Beta draws (N_samples x P)
beta_idx       <- grep("^beta\\[", colnames(draws_mat))
beta_draws_mar <- as.matrix(draws_mat[, beta_idx])
col_order      <- order(as.integer(gsub("beta\\[(\\d+)\\]", "\\1", colnames(beta_draws_mar))))
beta_draws_mar <- beta_draws_mar[, col_order]
colnames(beta_draws_mar) <- beta_names

# Alpha draws (N_samples x K)
alpha_idx       <- grep("^alpha\\[", colnames(draws_mat))
alpha_draws_mat <- as.matrix(draws_mat[, alpha_idx])
col_order_a     <- order(as.integer(gsub("alpha\\[(\\d+)\\]", "\\1", colnames(alpha_draws_mat))))
alpha_draws_mat <- alpha_draws_mat[, col_order_a]

if (!exists("alpha_hat") || is.null(alpha_hat)) alpha_hat <- colMeans(alpha_draws_mat)

# Reconstruct cal_summary / surv_merged from draws if missing from checkpoint
if (!exists("cal_summary") || is.null(cal_summary) ||
    !exists("surv_merged") || is.null(surv_merged)) {
  cat("[reconstruct] Computing calibration + survival from draws...\n")
  library(survival)
  eta_mat      <- outer(rep(1, dat$N), alpha_hat) + as.vector(dat$X %*% beta_A)
  h_mat        <- plogis(eta_mat)
  log_surv_mat <- t(apply(log(1 - h_mat), 1, cumsum))
  surv_mat     <- exp(log_surv_mat)
  
  # Calibration
  Ki_vec      <- pmin(dat$y_obs, K)
  pred_mort   <- 1 - surv_mat[cbind(1:dat$N, Ki_vec)]
  risk_decile <- cut(pred_mort, breaks = quantile(pred_mort, seq(0,1,0.1)),
                     include.lowest = TRUE, labels = 1:10)
  cal_summary <- aggregate(cbind(predicted = pred_mort, observed = dat$delta)
                           ~ as.numeric(risk_decile), FUN = mean)
  names(cal_summary)[1] <- "decile"
  
  # Survival + KM
  surv_df    <- data.frame(time = (1:K)*dat$delta_width,
                           model_surv = colMeans(surv_mat),
                           model_lo = apply(surv_mat,2,quantile,0.025),
                           model_hi = apply(surv_mat,2,quantile,0.975))
  km_fit     <- survfit(Surv(y_obs*dat$delta_width, delta) ~ 1,
                        data = data.frame(y_obs=dat$y_obs, delta=dat$delta))
  km_s       <- summary(km_fit, times = surv_df$time)
  km_df      <- data.frame(time=km_s$time, km_surv=km_s$surv,
                           km_lo=km_s$lower, km_hi=km_s$upper)
  surv_merged <- merge(surv_df, km_df, by="time", all.x=TRUE)
  rm(eta_mat, h_mat, log_surv_mat, surv_mat, surv_df, km_df, km_fit)
}

# Free large objects
rm(draws_A, draws_B, draws_mat, pp_data, patients_df); gc()

# Common settings
SIGMA_GAMMA_GRID <- sort(as.numeric(names(mnar_results)))
cols3    <- c("#2171b5", "#cb181d", "#238b45")
mar_mean <- colMeans(beta_draws_mar)


# ============================================================================
# FIGURE: BASELINE HAZARD (Model A)
# ============================================================================
cat("[plot] Baseline hazard...\n")

h_hat    <- plogis(alpha_hat)
h_draws  <- plogis(alpha_draws_mat)
hazard_df <- data.frame(day  = 1:K,
                        mean = h_hat,
                        lo   = apply(h_draws, 2, quantile, 0.025),
                        hi   = apply(h_draws, 2, quantile, 0.975))

p_hazard <- ggplot(hazard_df, aes(day, mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  labs(x = "ICU Day", y = expression("Baseline Hazard " * h[0](t)),
       title = "Estimated Baseline Hazard (Model A)",
       subtitle = "Posterior mean with 95% credible interval") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_baseline_hazard.pdf"), p_hazard, width = 8, height = 5)
rm(alpha_draws_mat, h_draws)


# ============================================================================
# FIGURE: CALIBRATION (Model A)
# ============================================================================
cat("[plot] Calibration...\n")

p_cal <- ggplot(cal_summary, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_text(aes(label = decile), vjust = -1, size = 3, color = "grey30") +
  labs(x = "Predicted Mortality (Decile Mean)", y = "Observed Mortality (Decile Mean)",
       title = "Risk-Decile Calibration (Model A)",
       subtitle = "Numbers indicate risk decile (1 = lowest, 10 = highest)") +
  coord_equal(xlim = c(0, max(cal_summary$predicted, cal_summary$observed) * 1.1),
              ylim = c(0, max(cal_summary$predicted, cal_summary$observed) * 1.1)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_calibration.pdf"), p_cal, width = 6, height = 6)


# ============================================================================
# FIGURE: SPAGHETTI SURVIVAL (MAR + MNAR)
# ============================================================================
cat("[plot] Spaghetti survival curves...\n")

N_DRAWS <- 200
set.seed(42)
t_axis <- surv_merged$time

# Subsample patients for speed
N_SUB   <- 5000
set.seed(123)
sub_ids <- sample(dat$N, min(N_SUB, dat$N))
X_sub   <- dat$X[sub_ids, , drop = FALSE]

mnar_keys <- as.character(sort(as.numeric(names(mnar_results))))

# Population-averaged survival from one beta vector
pop_surv <- function(beta_vec) {
  Xb    <- as.numeric(X_sub %*% as.numeric(beta_vec))
  eta   <- outer(Xb, alpha_hat, FUN = "+")
  log_S <- t(apply(log(1 - plogis(eta)), 1, cumsum))
  colMeans(exp(log_S))
}

# MAR draws
n_mar        <- nrow(beta_draws_mar)
idx_mar      <- sample(n_mar, min(N_DRAWS, n_mar))
mar_spaghetti <- do.call(rbind, lapply(idx_mar, function(d)
  data.frame(time = t_axis, surv = pop_surv(beta_draws_mar[d, ]), draw = d)))
mar_mean_surv <- pop_surv(mar_mean)

# MNAR draws per sigma_gamma
mnar_spaghetti <- list();  mnar_mean_surv <- list()
for (sg in mnar_keys) {
  bmat <- mnar_results[[sg]]$beta
  idx  <- sample(nrow(bmat), min(N_DRAWS, nrow(bmat)))
  mnar_spaghetti[[sg]] <- do.call(rbind, lapply(idx, function(d)
    data.frame(time = t_axis, surv = pop_surv(bmat[d, ]), draw = d)))
  mnar_mean_surv[[sg]] <- pop_surv(colMeans(bmat))
}

# Build plot
surv_labels <- c(KM = "'Kaplan-Meier'", MAR = "'MAR (Model A)'")
for (sg in mnar_keys)
  surv_labels[paste0("MNAR_", sg)] <- paste0("'MNAR (' * sigma[gamma] == ", sg, " * ')'")

p_surv <- ggplot() +
  geom_ribbon(data = data.frame(time = t_axis, lo = surv_merged$km_lo, hi = surv_merged$km_hi),
              aes(time, ymin = lo, ymax = hi), fill = "grey70", alpha = 0.15) +
  geom_line(data = mar_spaghetti, aes(time, surv, group = draw),
            color = "grey50", alpha = 0.04, linewidth = 0.3)

for (i in seq_along(mnar_keys))
  p_surv <- p_surv +
  geom_line(data = mnar_spaghetti[[mnar_keys[i]]],
            aes(time, surv, group = draw),
            color = cols3[i], alpha = 0.06, linewidth = 0.3)

# Bold mean lines
p_surv <- p_surv +
  geom_line(data = data.frame(time = t_axis, surv = surv_merged$km_surv),
            aes(time, surv, color = "KM", linetype = "KM"), linewidth = 0.9) +
  geom_line(data = data.frame(time = t_axis, surv = mar_mean_surv),
            aes(time, surv, color = "MAR", linetype = "MAR"), linewidth = 1.1)

for (i in seq_along(mnar_keys)) {
  key <- paste0("MNAR_", mnar_keys[i])
  p_surv <- p_surv +
    geom_line(data = data.frame(time = t_axis, surv = mnar_mean_surv[[mnar_keys[i]]]),
              aes(time, surv, color = !!key, linetype = !!key), linewidth = 1.1)
}

color_vals <- c(KM = "black", MAR = "grey40")
lty_vals   <- c(KM = "dashed", MAR = "solid")
for (i in seq_along(mnar_keys)) {
  k <- paste0("MNAR_", mnar_keys[i])
  color_vals[k] <- cols3[i];  lty_vals[k] <- "solid"
}

p_surv <- p_surv +
  scale_color_manual(values = color_vals, labels = parse(text = surv_labels)) +
  scale_linetype_manual(values = lty_vals, labels = parse(text = surv_labels)) +
  labs(x = "Days Since ICU Admission", y = "S(t)",
       title = "Population-Averaged Survival Under MAR and MNAR",
       subtitle = paste0(N_DRAWS, " posterior draws per model; bold = posterior mean"),
       color = "", linetype = "") +
  coord_cartesian(ylim = c(
    min(c(mar_spaghetti$surv,
          sapply(mnar_spaghetti, function(d) min(d$surv))), na.rm = TRUE) * 0.98, 1)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"), panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_survival.pdf"), p_surv, width = 9, height = 6)
rm(mar_spaghetti, mnar_spaghetti, mnar_mean_surv, X_sub)


# ============================================================================
# FIGURE: GAMMA POSTERIOR DENSITIES
# ============================================================================
cat("[plot] Gamma posterior...\n")

gamma_df <- do.call(rbind, lapply(names(mnar_results), function(sg)
  data.frame(sigma_gamma = paste0("sigma[gamma] == ", sg),
             gamma = as.numeric(mnar_results[[sg]]$gamma_draws))))
gamma_df$sigma_gamma <- factor(gamma_df$sigma_gamma,
                               levels = paste0("sigma[gamma] == ", SIGMA_GAMMA_GRID))

# Subsample for rug
set.seed(42)
rug_df <- do.call(rbind, lapply(split(gamma_df, gamma_df$sigma_gamma),
                                function(d) d[sample(nrow(d), min(1000, nrow(d))), ]))

# Prior curves
x_grid   <- seq(min(gamma_df$gamma) - 1, 2, length.out = 500)
prior_df <- do.call(rbind, lapply(SIGMA_GAMMA_GRID, function(sg)
  data.frame(sigma_gamma = paste0("sigma[gamma] == ", sg),
             x = x_grid, density = dnorm(x_grid, 0, sg))))
prior_df$sigma_gamma <- factor(prior_df$sigma_gamma, levels = levels(gamma_df$sigma_gamma))

p_gamma <- ggplot(gamma_df, aes(gamma)) +
  geom_rug(data = rug_df, aes(color = sigma_gamma),
           alpha = 0.05, length = unit(0.04, "npc"), show.legend = FALSE) +
  geom_density(aes(fill = sigma_gamma, color = sigma_gamma), alpha = 0.20, linewidth = 1.1) +
  geom_line(data = prior_df, aes(x, density, color = sigma_gamma),
            linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey40", linewidth = 0.8) +
  annotate("text", x = 0.2, y = Inf, label = "MAR (gamma == 0)", vjust = 2,
           parse = TRUE, color = "grey30", size = 3.5, fontface = "italic") +
  scale_fill_manual(name = "Prior Width", values = cols3,
                    labels = parse(text = levels(gamma_df$sigma_gamma))) +
  scale_color_manual(name = "Prior Width", values = cols3,
                     labels = parse(text = levels(gamma_df$sigma_gamma))) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.4, color = cols3))) +
  labs(x = expression(gamma ~ "(selection parameter)"), y = "Density",
       title = expression("Posterior of " * gamma * " Under Three Prior Widths"),
       subtitle = "Solid = posterior; Dashed = prior. Bottom ticks = individual posterior draws.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.15, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80"),
        legend.title = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_mnar_gamma_posterior.pdf"), p_gamma, width = 8, height = 5)


# ============================================================================
# FIGURE: EMERGENCY COVARIATE DEEP DIVE
# ============================================================================
cat("[plot] Emergency covariate...\n")

emerg_idx <- which(beta_names == "emergency")

emerg_df <- do.call(rbind, c(
  list(data.frame(value = as.numeric(beta_draws_mar[, emerg_idx]),
                  model = "MAR (Model A)")),
  lapply(names(mnar_results), function(sg)
    data.frame(value = as.numeric(mnar_results[[sg]]$beta[, emerg_idx]),
               model = paste0("MNAR (sigma_gamma = ", sg, ")")))
))

model_levels <- c("MAR (Model A)",
                  paste0("MNAR (sigma_gamma = ", SIGMA_GAMMA_GRID, ")"))
emerg_df$model <- factor(emerg_df$model, levels = model_levels)

# Rug subsample
set.seed(42)
rug_e <- do.call(rbind, lapply(split(emerg_df, emerg_df$model),
                               function(d) d[sample(nrow(d), min(1000, nrow(d))), ]))

# Summary stats
emerg_stats <- do.call(rbind, lapply(levels(emerg_df$model), function(m) {
  v <- emerg_df$value[emerg_df$model == m]
  data.frame(model = m, post_mean = mean(v),
             q025 = quantile(v, 0.025), q975 = quantile(v, 0.975))
}))

emerg_cols <- c("MAR (Model A)" = "grey30",
                setNames(cols3, paste0("MNAR (sigma_gamma = ", SIGMA_GAMMA_GRID, ")")))

parsed_labels <- parse(text = c("'MAR (Model A)'",
                                paste0("'MNAR (' * sigma[gamma] == ", SIGMA_GAMMA_GRID, " * ')'")))

p_emerg <- ggplot(emerg_df, aes(value, fill = model, color = model)) +
  geom_rug(data = rug_e, alpha = 0.05, length = unit(0.04, "npc"), show.legend = FALSE) +
  geom_density(alpha = 0.15, linewidth = 1.0) +
  geom_vline(data = emerg_stats, aes(xintercept = post_mean, color = model),
             linetype = "solid", linewidth = 0.7, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey60", linewidth = 0.6) +
  annotate("text", x = emerg_stats$post_mean[1], y = Inf,
           label = sprintf("MAR: %.3f", emerg_stats$post_mean[1]),
           vjust = 2, hjust = -0.1, size = 3.3, color = "grey30", fontface = "bold") +
  annotate("text", x = emerg_stats$post_mean[nrow(emerg_stats)], y = Inf,
           label = sprintf("MNAR (2.0): %.3f", emerg_stats$post_mean[nrow(emerg_stats)]),
           vjust = 3.5, hjust = -0.1, size = 3.3, color = cols3[3], fontface = "bold") +
  annotate("segment",
           x = emerg_stats$post_mean[1], xend = emerg_stats$post_mean[nrow(emerg_stats)],
           y = 0.08, yend = 0.08,
           arrow = arrow(length = unit(0.15, "cm"), ends = "both"), color = "grey30") +
  annotate("text",
           x = mean(c(emerg_stats$post_mean[1], emerg_stats$post_mean[nrow(emerg_stats)])),
           y = 0.08,
           label = sprintf("Delta*beta == %.3f",
                           emerg_stats$post_mean[nrow(emerg_stats)] - emerg_stats$post_mean[1]),
           parse = TRUE, vjust = -0.5, size = 3.2, color = "grey20") +
  scale_fill_manual(values = emerg_cols, labels = parsed_labels) +
  scale_color_manual(values = emerg_cols, labels = parsed_labels) +
  labs(x = expression(beta[emergency] ~ "(log-odds scale)"),
       y = "Posterior Density",
       title = expression("Emergency Admission Effect: MAR vs. MNAR Across " * sigma[gamma]),
       fill = "Model", color = "Model") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.15, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_mnar_emergency.pdf"), p_emerg, width = 9, height = 6)


# ============================================================================
# DONE
# ============================================================================
cat("\n=== All paper figures saved to", plot_dir, "===\n")
cat("  mimic_baseline_hazard.pdf\n")
cat("  mimic_calibration.pdf\n")
cat("  mimic_survival.pdf\n")
cat("  mimic_mnar_gamma_posterior.pdf\n")
cat("  mimic_mnar_emergency.pdf\n")
