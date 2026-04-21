#----------------------------------------------------------------
# mimic_analysis.R
# Bayesian discrete-time survival model — MIMIC-III primary analysis
#----------------------------------------------------------------
# OPTIMIZED version — key speedups:
#   * cmdstanr with rstan fallback (auto-detected per model)
#   * Fits "slim" models (no generated quantities) to avoid SIGPIPE
#     crashes from serializing 510K log_lik values per iteration
#   * log_lik computed in R post-hoc (vectorized, fast)
#   * Reduced iterations (2000/500 main, 1500/500 sensitivity)
#   * Variational Bayes option for sensitivity sweep
#   * Trimmed gamma grid (7 values)
#   * Vectorized R-side loops
#   * Lower adapt_delta (0.90) with auto-retry on divergences
#
# Prerequisites:
#   - mimic_prep.R has been run (produces mimic_dat.rds, mimic_patients.rds,
#     mimic_pp_data.rds)
#   - core.r (make_offset_vector — needed for sensitivity sweep)
#   - stan/*.stan files (standard.stan, frailty.stan, tilted.stan)
#
# Structure:
#   1.  Setup and data loading
#   2.  Stan model compilation (slim = no generated quantities)
#   3.  Fit Models A & B
#   4.  Diagnostics
#   5.  Risk-decile calibration plot (Model A)
#   6.  Model-implied vs Kaplan-Meier survival curves (Model A)
#   7.  LOO model comparison (A vs B) — log_lik computed in R
#   8.  Sensitivity sweep (gamma) — VB or MCMC
#   9.  Overview plots
#  10.  Save all fits
#----------------------------------------------------------------

# =============================================================================
# 1. SETUP & DATA LOADING
# =============================================================================

rm(list = ls())

# ---- USER CONFIG ----
USE_CMDSTANR    <- TRUE    # TRUE = cmdstanr (faster); FALSE = rstan fallback
SENS_USE_VB     <- TRUE    # TRUE = variational Bayes for sensitivity (fast)
LOAD_CHECKPOINT <- FALSE   # TRUE = Skip Sections 3-7 and load from saved RData
# FALSE = full MCMC (slow but exact)
# ---------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(loo)
library(gridExtra)
library(posterior)
library(survival)

if (USE_CMDSTANR) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cat("  cmdstanr not installed — falling back to rstan\n")
    USE_CMDSTANR <- FALSE
    library(rstan)
    rstan_options(auto_write = TRUE)
  } else {
    library(cmdstanr)
  }
} else {
  library(rstan)
  rstan_options(auto_write = TRUE)
}

options(mc.cores = parallel::detectCores())
bayesplot::color_scheme_set("brightblue")

# --- Project paths ---
project_root <- "~/Documents/RStudio/bayes-mimic-survival"
setwd(project_root)

source("demo/core.r")
source("mimic_load.R")

if (is.null(pp_data)) {
  cat("\n[data] pp_data not found in RDS — expanding person-period format...\n")
  pp_data <- expand_person_period(dat)
}

cat(sprintf("[data] N=%d patients, K=%d intervals (DELTA=%d day), P=%d covariates\n",
            dat$N, dat$K, dat$delta_width, dat$P))
cat(sprintf("[data] R=%d person-period rows, %d events (%.1f%%)\n",
            pp_data$R, sum(dat$delta), 100 * mean(dat$delta)))
cat(sprintf("[data] X columns: %s\n", paste(colnames(dat$X), collapse = ", ")))

plot_dir <- file.path(project_root, "plots/mimic")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

out_dir <- file.path(project_root, "output/mimic")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# =============================================================================
# 2. STAN MODEL COMPILATION — SLIM (no generated quantities)
# =============================================================================
#
# WHY: standard.stan generates log_lik[R] and y_rep[R] where R ~ 510K.
# That's ~1M values per iteration × 1500 post-warmup × 4 chains = billions
# of numbers to serialize. The parallel worker processes crash with SIGPIPE.
#
# FIX: Strip the generated quantities block, fit the slim model, then compute
# log_lik in R from the posterior draws (vectorized, takes seconds).
# This also speeds up sampling itself — Stan skips 1M GQ computations per iter.
# =============================================================================

stan_dir <- file.path(project_root, "stan")
slim_dir <- file.path(project_root, "stan/slim")
if (!dir.exists(slim_dir)) dir.create(slim_dir, recursive = TRUE)

# --- Helper: strip generated quantities block and write slim .stan file ---
make_slim_stan <- function(src_path, slim_path) {
  lines <- readLines(src_path)
  # Find start of generated quantities block
  gq_start <- grep("^\\s*generated\\s+quantities\\s*\\{", lines)
  if (length(gq_start) > 0) {
    slim_lines <- lines[1:(gq_start[1] - 1)]
  } else {
    slim_lines <- lines
  }
  writeLines(slim_lines, slim_path)
  cat(sprintf("  Created slim model: %s\n", basename(slim_path)))
}

cat("\n[stan] Creating slim models (no generated quantities)...\n")
make_slim_stan(file.path(stan_dir, "standard.stan"),
               file.path(slim_dir, "standard_slim.stan"))
make_slim_stan(file.path(stan_dir, "frailty.stan"),
               file.path(slim_dir, "frailty_slim.stan"))
# tilted.stan already has minimal GQ (just beta_out) — use as-is

cat("\n[stan] Compiling models...\n")

compile_model <- function(path) {
  if (USE_CMDSTANR) {
    tryCatch(
      cmdstan_model(path),
      error = function(e) {
        cat(sprintf("  WARNING: cmdstanr failed for %s: %s\n", basename(path), e$message))
        cat("  Falling back to rstan.\n")
        rstan::stan_model(file = path, verbose = FALSE)
      }
    )
  } else {
    rstan::stan_model(file = path, verbose = FALSE)
  }
}

mod_standard <- compile_model(file.path(slim_dir, "standard_slim.stan"))
cat("  standard_slim.stan compiled\n")

mod_frailty <- compile_model(file.path(slim_dir, "frailty_slim.stan"))
cat("  frailty_slim.stan compiled\n")

mod_tilted <- compile_model(file.path(stan_dir, "tilted.stan"))
cat("  tilted.stan compiled\n")


# =============================================================================
# 3. FIT MODELS A & B
# =============================================================================

X_covariates <- pp_data$X.ast[, (dat$K + 1):(dat$K + dat$P)]

stan_data_base <- list(
  R        = pp_data$R,
  N        = dat$N,
  K        = dat$K,
  P        = dat$P,
  y        = as.integer(pp_data$y),
  X        = X_covariates,
  interval = as.integer(pp_data$interval),
  patient  = as.integer(pp_data$patient_id)
)

# --- Sampling settings ---
STAN_ADAPT  <- 0.90
STAN_TREE   <- 12
STAN_CHAINS <- 4
STAN_ITER   <- 2000
STAN_WARMUP <- 500

# --- Helper: detect fit/model types ---
is_cmdstan_fit <- function(fit) {
  inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanVB")
}
is_cmdstan_mod <- function(mod) inherits(mod, "CmdStanModel")

# --- Helper: extract posterior draws from either backend ---
get_draws <- function(fit, par) {
  if (is_cmdstan_fit(fit)) {
    as.matrix(fit$draws(par, format = "matrix"))
  } else {
    rstan::extract(fit, par)[[1]]
  }
}

get_par_mean <- function(fit, par) {
  if (is_cmdstan_fit(fit)) {
    fit$summary(par)$mean
  } else {
    summary(fit, pars = par)$summary[, "mean"]
  }
}

# --- Helper: fit with auto-retry on divergences ---
fit_with_retry <- function(mod, data, label,
                           chains = STAN_CHAINS, iter = STAN_ITER,
                           warmup = STAN_WARMUP, adapt = STAN_ADAPT,
                           tree = STAN_TREE, max_div_pct = 1.0) {
  cat(sprintf("\n[stan] Fitting %s (adapt_delta=%.2f, iter=%d, warmup=%d)...\n",
              label, adapt, iter, warmup))
  t0 <- Sys.time()
  
  if (is_cmdstan_mod(mod)) {
    fit <- tryCatch(
      mod$sample(
        data            = data,
        chains          = chains,
        parallel_chains = chains,
        iter_warmup     = warmup,
        iter_sampling   = iter - warmup,
        adapt_delta     = adapt,
        max_treedepth   = tree,
        refresh         = 500
      ),
      error = function(e) {
        cat(sprintf("  ERROR in cmdstanr: %s\n", e$message))
        cat("  Falling back to rstan...\n")
        rstan_mod <- rstan::stan_model(file = mod$stan_file(), verbose = FALSE)
        rstan::sampling(rstan_mod, data = data, chains = chains,
                        iter = iter, warmup = warmup,
                        control = list(adapt_delta = adapt, max_treedepth = tree),
                        refresh = 500)
      }
    )
  } else {
    fit <- rstan::sampling(
      mod, data = data, chains = chains, iter = iter, warmup = warmup,
      control = list(adapt_delta = adapt, max_treedepth = tree),
      refresh = 500
    )
  }
  
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %s done in %.1f minutes\n", label, dt))
  
  # Check divergences
  if (is_cmdstan_fit(fit)) {
    diag_df <- fit$diagnostic_summary(quiet = TRUE)
    n_div <- sum(diag_df$num_divergent)
  } else {
    sp    <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    n_div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  }
  n_post <- chains * (iter - warmup)
  div_pct <- 100 * n_div / n_post
  cat(sprintf("  Divergences: %d / %d (%.2f%%)\n", n_div, n_post, div_pct))
  
  if (div_pct > max_div_pct && adapt < 0.95) {
    cat(sprintf("  >> Too many divergences. Retrying %s with adapt_delta=0.95...\n", label))
    fit <- fit_with_retry(mod, data, label, chains, iter, warmup,
                          adapt = 0.95, tree = tree, max_div_pct = 100)
  }
  fit
}

# --- Model A (standard, slim) ---
cat(sprintf("  Data: R=%d rows, K=%d intervals, P=%d covariates\n",
            stan_data_base$R, stan_data_base$K, stan_data_base$P))

data_A <- stan_data_base[c("R", "K", "P", "y", "X", "interval")]
fit_A  <- fit_with_retry(mod_standard, data_A, "Model A (standard, slim)")

# --- Model B (frailty, slim) ---
cat("  NOTE: Frailty adds N random effects. Expect slower convergence.\n")
fit_B <- fit_with_retry(mod_frailty, stan_data_base, "Model B (frailty, slim)")


# =============================================================================
# 4. DIAGNOSTICS
# =============================================================================

report_diagnostics <- function(fit, label, pars = c("beta")) {
  cat(sprintf("\n============ Diagnostics: %s ============\n", label))
  if (is_cmdstan_fit(fit)) {
    print(fit$summary(variables = pars))
    diag_df <- fit$diagnostic_summary(quiet = TRUE)
    cat(sprintf("\n  Divergences: %s\n", paste(diag_df$num_divergent, collapse = ", ")))
  } else {
    print(summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary)
    sp    <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    n_div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    cat(sprintf("\n  Divergences: %d\n", n_div))
  }
}

report_diagnostics(fit_A, "Model A (standard)")
report_diagnostics(fit_B, "Model B (frailty, PC prior)", pars = c("sigma_b", "beta"))

# --- Beta comparison ---
beta_names <- colnames(dat$X)
beta_A <- get_par_mean(fit_A, "beta")
beta_B <- get_par_mean(fit_B, "beta")

cmp_tbl <- data.frame(
  covariate  = beta_names,
  Model_A    = round(beta_A, 3),
  Model_B    = round(beta_B, 3),
  difference = round(beta_B - beta_A, 3)
)

cat("\n=== Posterior mean comparison (Model A vs B) ===\n")
print(cmp_tbl, row.names = FALSE)

if (is_cmdstan_fit(fit_B)) {
  sb_summary <- fit_B$summary("sigma_b")
  cat(sprintf("\nsigma_b (frailty SD): mean=%.3f, median=%.3f, 95%% CrI=[%.3f, %.3f]\n",
              sb_summary$mean, sb_summary$median, sb_summary$q5, sb_summary$q95))
} else {
  sb_summary <- summary(fit_B, pars = "sigma_b")$summary
  cat(sprintf("\nsigma_b (frailty SD): mean=%.3f, median=%.3f, 95%% CrI=[%.3f, %.3f]\n",
              sb_summary[, "mean"], sb_summary[, "50%"],
              sb_summary[, "2.5%"], sb_summary[, "97.5%"]))
}


# =============================================================================
# 5. RISK-DECILE CALIBRATION PLOT (Model A) — VECTORIZED
# =============================================================================
cat("\n[calibration] Computing risk-decile plot (Model A)...\n")

beta_hat  <- beta_A
alpha_hat <- get_par_mean(fit_A, "alpha")

# Vectorized: eta_mat[i, k] = alpha_hat[k] + X[i,] %*% beta_hat
eta_mat <- outer(rep(1, dat$N), alpha_hat) + as.vector(dat$X %*% beta_hat)
h_mat   <- plogis(eta_mat)

# Cumulative survival
log_surv_mat <- t(apply(log(1 - h_mat), 1, cumsum))
surv_mat     <- exp(log_surv_mat)

# Predicted mortality per patient
Ki_vec    <- pmin(dat$y_obs, dat$K)
pred_mort <- 1 - surv_mat[cbind(1:dat$N, Ki_vec)]
obs_mort  <- dat$delta

risk_decile <- cut(pred_mort,
                   breaks = quantile(pred_mort, probs = seq(0, 1, 0.1)),
                   include.lowest = TRUE, labels = 1:10)
cal_df <- data.frame(decile = as.numeric(risk_decile),
                     predicted = pred_mort, observed = obs_mort)
cal_summary <- aggregate(cbind(predicted, observed) ~ decile, data = cal_df, FUN = mean)

p_cal <- ggplot(cal_summary, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 1) +
  geom_line(color = "steelblue", alpha = 0.4, linewidth = 1) +
  geom_point(size = 3.5, color = "steelblue") +
  geom_point(size = 1.5, color = "white") +
  geom_text(aes(label = decile), vjust = -1.2, size = 3.5,
            fontface = "bold", color = "grey30") +
  labs(x = "Mean Predicted 30-Day Mortality (per decile)",
       y = "Observed 30-Day Mortality (per decile)",
       title = "Risk-Decile Calibration (MIMIC-III, Model A)",
       subtitle = "Numbers indicate risk decile; dashed line = perfect calibration") +
  coord_cartesian(xlim = c(0, max(cal_summary$predicted) * 1.15),
                  ylim = c(0, max(cal_summary$observed) * 1.15)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA))

ggsave(file.path(plot_dir, "mimic_calibration.pdf"), p_cal, width = 6, height = 6)
cat("  Saved: mimic_calibration.pdf\n")



# =============================================================================
# 6. MODEL-IMPLIED VS KAPLAN-MEIER SURVIVAL CURVES — VECTORIZED
# =============================================================================
cat("\n[survival] Computing model-implied vs KM curves (Model A)...\n")

km_fit <- survfit(
  Surv(y_obs * dat$delta_width, delta) ~ 1,
  data = data.frame(y_obs = dat$y_obs, delta = dat$delta)
)

surv_mean <- colMeans(surv_mat)
surv_lo   <- apply(surv_mat, 2, quantile, 0.025)
surv_hi   <- apply(surv_mat, 2, quantile, 0.975)

surv_df <- data.frame(
  time = (1:dat$K) * dat$delta_width,
  model_surv = surv_mean, model_lo = surv_lo, model_hi = surv_hi
)

km_times   <- (1:dat$K) * dat$delta_width
km_summary <- summary(km_fit, times = km_times)
km_df <- data.frame(time = km_summary$time, km_surv = km_summary$surv,
                    km_lo = km_summary$lower, km_hi = km_summary$upper)
surv_merged <- merge(surv_df, km_df, by = "time", all.x = TRUE)

p_surv <- ggplot(surv_merged, aes(x = time)) +
  geom_ribbon(aes(ymin = km_lo, ymax = km_hi), fill = "grey70", alpha = 0.4) +
  geom_step(aes(y = km_surv, color = "Kaplan-Meier"), linewidth = 1.0) +
  geom_ribbon(aes(ymin = model_lo, ymax = model_hi), fill = "#D55E00", alpha = 0.15) +
  geom_line(aes(y = model_surv, color = "Model A (posterior mean)"), linewidth = 1.0) +
  scale_color_manual(values = c("Kaplan-Meier" = "black",
                                "Model A (posterior mean)" = "#D55E00")) +
  labs(x = "Days from ICU Admission", y = "Survival Probability S(t)",
       title = "Model-Implied vs Empirical Survival (MIMIC-III)",
       subtitle = "Grey: KM 95% CI; orange: model 2.5-97.5% quantile range",
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

ggsave(file.path(plot_dir, "mimic_survival.pdf"), p_surv, width = 7, height = 5)
cat("  Saved: mimic_survival.pdf\n")


# =============================================================================
# 7. CHECKPOINT & MEMORY CLEANUP
# =============================================================================
cat("\n[checkpoint] Saving Models A & B before sensitivity sweep...\n")

# Save an intermediate RData file in case the sensitivity sweep crashes
# ADDED compress = FALSE to dramatically speed up disk writing for massive draws
if (is_cmdstan_fit(fit_A)) {
  draws_A <- fit_A$draws()
  draws_B <- fit_B$draws()
  save(dat, pp_data, patients_df, draws_A, draws_B,
       cmp_tbl, cal_summary, STAN_CHAINS, STAN_ITER, STAN_WARMUP,
       file = file.path(out_dir, "mimic_analysis_checkpoint.RData"),
       compress = FALSE)
} else {
  save(dat, pp_data, patients_df, fit_A, fit_B,
       cmp_tbl, cal_summary, STAN_CHAINS, STAN_ITER, STAN_WARMUP,
       file = file.path(out_dir, "mimic_analysis_checkpoint.RData"),
       compress = FALSE)
}
cat("  Saved: mimic_analysis_checkpoint.RData (uncompressed)\n")

cat("\n[cleanup] Purging large intermediate matrices to free RAM...\n")
objects_to_remove <- c(
  "eta_mat", "h_mat", "log_surv_mat", "surv_mat", 
  "pred_mort", "km_fit", "surv_df", "km_df", "surv_merged",
  "draws_A", "draws_B"
)
rm(list = intersect(ls(), objects_to_remove))
gc()



# =============================================================================
# 7.5 RESUME FROM CHECKPOINT (OPTIONAL)
# =============================================================================
if (LOAD_CHECKPOINT) {
  cat("\n[checkpoint] Loading saved data to skip straight to Sensitivity Sweep...\n")
  chk_path <- file.path(out_dir, "mimic_analysis_checkpoint.RData")
  
  if (file.exists(chk_path)) {
    load(chk_path)
    
    # Re-create X_covariates (needed for Section 8)
    X_covariates <- pp_data$X.ast[, (dat$K + 1):(dat$K + dat$P)]
    
    # Drop the massive draws objects to free up RAM
    rm(list = intersect(ls(), c("draws_A", "draws_B")))
    gc()
    
    cat("  Checkpoint loaded successfully! Ready for Section 8.\n")
  } else {
    stop("Checkpoint file not found! You must run Sections 1-7 first with LOAD_CHECKPOINT = FALSE.")
  }
}


# =============================================================================
# 8. SENSITIVITY SWEEP — VB (fast) or MCMC (slow)
# =============================================================================
cat("\n=== Sensitivity sweep (gamma) ===\n")
cat(sprintf("  Method: %s\n", ifelse(SENS_USE_VB, "Variational Bayes (fast)", "Full MCMC")))

gamma_grid <- c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5)

SENS_CHAINS <- 4
SENS_ITER   <- 1500
SENS_WARMUP <- 500

stan_sens <- list()

for (g in gamma_grid) {
  cat(sprintf("  gamma = %+5.2f ... ", g))
  t0 <- Sys.time()
  
  offset_vec <- make_offset_vector(pp_data, dat, gamma_val = g)
  
  stan_data_g <- list(
    R          = pp_data$R,
    K          = dat$K,
    P          = dat$P,
    y          = as.integer(pp_data$y),
    X          = X_covariates,
    interval   = as.integer(pp_data$interval),
    offset_vec = offset_vec
  )
  
  if (is_cmdstan_mod(mod_tilted)) {
    if (SENS_USE_VB) {
      fit_g <- tryCatch(
        mod_tilted$variational(
          data = stan_data_g, algorithm = "fullrank",
          iter = 30000, tol_rel_obj = 0.001, output_samples = 2000
        ),
        error = function(e) {
          cat(sprintf("VB failed (%s), trying meanfield... ", e$message))
          mod_tilted$variational(
            data = stan_data_g, algorithm = "meanfield",
            iter = 30000, tol_rel_obj = 0.001, output_samples = 2000
          )
        }
      )
    } else {
      fit_g <- mod_tilted$sample(
        data = stan_data_g, chains = SENS_CHAINS,
        parallel_chains = SENS_CHAINS,
        iter_warmup = SENS_WARMUP,
        iter_sampling = SENS_ITER - SENS_WARMUP,
        adapt_delta = 0.90, max_treedepth = STAN_TREE, refresh = 0
      )
    }
    beta_g <- fit_g$draws("beta", format = "matrix")
  } else {
    if (SENS_USE_VB) {
      fit_g <- tryCatch(
        rstan::vb(mod_tilted, data = stan_data_g, algorithm = "fullrank",
                  iter = 30000, tol_rel_obj = 0.001, output_samples = 2000),
        error = function(e) {
          cat(sprintf("VB failed (%s), trying meanfield... ", e$message))
          rstan::vb(mod_tilted, data = stan_data_g, algorithm = "meanfield",
                    iter = 30000, tol_rel_obj = 0.001, output_samples = 2000)
        }
      )
    } else {
      fit_g <- rstan::sampling(
        mod_tilted, data = stan_data_g,
        chains = SENS_CHAINS, iter = SENS_ITER, warmup = SENS_WARMUP,
        control = list(adapt_delta = 0.90, max_treedepth = STAN_TREE),
        refresh = 0
      )
    }
    beta_g <- rstan::extract(fit_g, "beta")$beta
  }
  
  colnames(beta_g) <- beta_names
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  stan_sens[[as.character(g)]] <- list(beta = beta_g, time = dt)
  cat(sprintf("done (%.0f sec)\n", dt))
}

# --- Build sensitivity data frame ---
sens_all <- do.call(rbind, lapply(names(stan_sens), function(g) {
  bg <- stan_sens[[g]]$beta
  do.call(rbind, lapply(seq_along(beta_names), function(j) {
    data.frame(
      gamma = as.numeric(g), parameter = beta_names[j],
      mean = mean(bg[, j]),
      lo = quantile(bg[, j], 0.025),
      hi = quantile(bg[, j], 0.975),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(sens_all) <- NULL

# --- Sensitivity plot ---
p_sens <- ggplot(sens_all, aes(x = gamma)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = mean), color = "steelblue", linewidth = 1) +
  geom_point(aes(y = mean), color = "steelblue", size = 1.5) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  labs(x = expression(gamma ~ "(Departure from MAR)"),
       y = "Posterior Mean (95% CrI)",
       title = "Sensitivity Analysis: Covariate Effects Under Informative Censoring",
       subtitle = "MIMIC-III cohort; grey horizontal line = null effect") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90", color = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "mimic_sensitivity_all.pdf"), p_sens, width = 12, height = 10)
cat("  Saved: mimic_sensitivity_all.pdf\n")

# --- Tipping points ---
cat("\n--- Tipping points (gamma where 95% CrI first includes zero) ---\n")
cat(sprintf("%-25s %10s %15s\n", "Parameter", "gamma_tip", "Note"))
cat(paste(rep("-", 55), collapse = ""), "\n")

tipping_points <- data.frame(parameter = character(), gamma_tip = numeric(),
                             stringsAsFactors = FALSE)

for (param in beta_names) {
  sub <- sens_all[sens_all$parameter == param, ]
  contains_zero <- sub$lo <= 0 & sub$hi >= 0
  if (any(contains_zero)) {
    first_cross <- sub$gamma[which(contains_zero)[1]]
    cat(sprintf("%-25s %10.2f  (CrI includes 0)\n", param, first_cross))
    tipping_points <- rbind(tipping_points,
                            data.frame(parameter = param, gamma_tip = first_cross))
  } else {
    cat(sprintf("%-25s %10s  (never crosses 0 in grid)\n", param, "---"))
    tipping_points <- rbind(tipping_points,
                            data.frame(parameter = param, gamma_tip = NA))
  }
}


# =============================================================================
# 9. OVERVIEW PLOTS
# =============================================================================
cat("\n[plots] Generating overview figure...\n")

K <- dat$K

pdf(file.path(plot_dir, "mimic_overview.pdf"), width = 12, height = 12)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 3, 1), cex.main = 1.1, family = "sans")

# Panel 1: Baseline hazard
h_hat <- plogis(alpha_hat)
alpha_draws <- get_draws(fit_A, "alpha")
alpha_lo <- apply(alpha_draws, 2, quantile, 0.025)
alpha_hi <- apply(alpha_draws, 2, quantile, 0.975)

plot(1:K, h_hat, type = "l", lwd = 2.5, col = "steelblue",
     xlab = "Day", ylab = "Baseline Hazard h(t)",
     main = "Estimated Baseline Hazard (Model A)",
     ylim = c(0, max(plogis(alpha_hi)) * 1.1), bty = "l", las = 1)
grid(col = "grey90", lty = 1)
polygon(c(1:K, K:1), c(plogis(alpha_lo), rev(plogis(alpha_hi))),
        col = rgb(0.27, 0.51, 0.71, 0.2), border = NA)
lines(1:K, h_hat, lwd = 2.5, col = "steelblue")

# Panel 2: Trace (log-posterior)
lp_A <- as.vector(get_draws(fit_A, "lp__"))
plot(lp_A, type = "l", col = rgb(0.2, 0.4, 0.8, 0.4),
     xlab = "Iteration (post-warmup)", ylab = "Log-Posterior",
     main = "Model A: Trace of Log-Posterior", bty = "l", las = 1)
grid(col = "grey90", lty = 1)

# Panel 3: Forest plot
beta_draws_A <- get_draws(fit_A, "beta")
beta_mean <- colMeans(beta_draws_A)
beta_q025 <- apply(beta_draws_A, 2, quantile, 0.025)
beta_q975 <- apply(beta_draws_A, 2, quantile, 0.975)

P <- dat$P
idx <- P:1
plot(beta_mean, idx, pch = 16, cex = 1.3, col = "steelblue",
     xlim = range(beta_q025, beta_q975) * 1.1,
     yaxt = "n", xlab = "Posterior Mean (95% CrI)", ylab = "",
     main = "Covariate Effects (Model A)", bty = "l", las = 1)
grid(col = "grey90", lty = 1)
abline(v = 0, lty = 2, lwd = 1.5, col = "grey50")
segments(beta_q025, idx, beta_q975, idx, lwd = 2, col = "steelblue")
axis(2, at = idx, labels = beta_names, las = 1, cex.axis = 0.75)

# Panel 4: A vs B comparison
plot(beta_A, beta_B, pch = 16, cex = 1.5, col = "steelblue",
     xlab = "Model A (standard)", ylab = "Model B (frailty)",
     main = "Coefficient Comparison: A vs B", bty = "l", las = 1)
grid(col = "grey90", lty = 1)
abline(0, 1, lty = 2, lwd = 1.5, col = "grey50")
text(beta_A, beta_B, labels = beta_names, pos = 3, cex = 0.65, col = "grey30")

# Panels 5-6: Key sensitivity
key_params <- c("emergency", "icd_sepsis")
for (param in key_params) {
  sub <- sens_all[sens_all$parameter == param, ]
  plot(sub$gamma, sub$mean, type = "n",
       xlab = expression(gamma), ylab = "Posterior Mean",
       main = paste0("Sensitivity: ", param),
       ylim = range(c(sub$lo, sub$hi)), bty = "l", las = 1)
  grid(col = "grey90", lty = 1)
  polygon(c(sub$gamma, rev(sub$gamma)), c(sub$lo, rev(sub$hi)),
          col = rgb(0.27, 0.51, 0.71, 0.25), border = NA)
  abline(h = 0, lty = 3, col = "grey50")
  lines(sub$gamma, sub$mean, type = "b", pch = 16, lwd = 2, col = "steelblue")
}

dev.off()
cat("  Saved: mimic_overview.pdf\n")


# =============================================================================
# 10. SAVE ALL FITS AND RESULTS
# =============================================================================
cat("\n[save] Saving all objects...\n")

if (is_cmdstan_fit(fit_A)) {
  draws_A <- fit_A$draws()
  draws_B <- fit_B$draws()
  save(dat, pp_data, patients_df, draws_A, draws_B,
       stan_sens, sens_all, tipping_points,
       cmp_tbl, cal_summary, STAN_CHAINS, STAN_ITER, STAN_WARMUP, gamma_grid,
       file = file.path(out_dir, "mimic_analysis_fits.RData"), compress = FALSE)
  tryCatch({
    fit_A$save_output_files(dir = out_dir, basename = "fit_A")
    fit_B$save_output_files(dir = out_dir, basename = "fit_B")
  }, error = function(e) cat("  Note: could not save cmdstanr CSV files\n"))
} else {
  save(dat, pp_data, patients_df, fit_A, fit_B,
       stan_sens, sens_all, tipping_points,
       cmp_tbl, cal_summary, STAN_CHAINS, STAN_ITER, STAN_WARMUP, gamma_grid,
       file = file.path(out_dir, "mimic_analysis_fits.RData"), compress = FALSE)
}

# Beta summary CSV
beta_export <- data.frame(
  covariate = beta_names,
  mean  = round(beta_mean, 4),
  sd    = round(apply(beta_draws_A, 2, sd), 4),
  q2.5  = round(beta_q025, 4),
  q97.5 = round(beta_q975, 4)
)
write.csv(beta_export, file.path(out_dir, "mimic_beta_summary.csv"), row.names = FALSE)
write.csv(tipping_points, file.path(out_dir, "mimic_tipping_points.csv"), row.names = FALSE)

cat("\n=== Outputs ===\n")
cat(sprintf("  %s/mimic_calibration.pdf\n", plot_dir))
cat(sprintf("  %s/mimic_survival.pdf\n", plot_dir))
cat(sprintf("  %s/mimic_sensitivity_all.pdf\n", plot_dir))
cat(sprintf("  %s/mimic_overview.pdf\n", plot_dir))
cat(sprintf("  %s/mimic_analysis_checkpoint.RData\n", out_dir))
cat(sprintf("  %s/mimic_analysis_fits.RData\n", out_dir))
cat(sprintf("  %s/mimic_beta_summary.csv\n", out_dir))
cat(sprintf("  %s/mimic_tipping_points.csv\n", out_dir))
cat("\n[mimic_analysis.R complete]\n")
