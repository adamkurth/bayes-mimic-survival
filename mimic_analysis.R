# ============================================================================
# mimic_analysis.R — Bayesian discrete-time survival (MIMIC-III)
# ============================================================================
# Fits three models via CmdStanR:
#   Model A  — standard logistic hazard (slim, no GQ)
#   Model B  — frailty extension        (slim, no GQ)
#   MNAR     — selection-model sensitivity (ADVI, 3 prior widths)
#
# Saves two checkpoint files consumed by mimic_plots.R:
#   mimic_analysis_checkpoint.RData  (dat, draws, calibration, survival)
#   mimic_mnar_results.RData         (MNAR posterior draws + gamma)
#
# Prerequisites:
#   mimic_prep.R:   mimic_dat.rds, mimic_patients.rds, mimic_pp_data.rds
#   mimic_load.R:   loads the above into dat, pp_data, patients_df
#   stan/*.stan:    standard.stan, frailty.stan, mnar_sensitivity.stan
# ============================================================================

rm(list = ls())

# ---- Config ----------------------------------------------------------------
LOAD_CHECKPOINT <- TRUE   # TRUE = skip fitting, load saved checkpoint
LOAD_MNAR       <- TRUE   # TRUE = skip MNAR, load saved results

STAN_CHAINS <- 4
STAN_WARMUP <- 500
STAN_ITER   <- 2000        # sampling iterations (post-warmup)
STAN_ADAPT  <- 0.90
STAN_TREE   <- 12
# ----------------------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(data.table)
library(survival)

options(mc.cores = parallel::detectCores())

# ---- Paths -----------------------------------------------------------------
project_root <- "~/Documents/RStudio/bayes-mimic-survival"
setwd(project_root)

out_dir  <- file.path(project_root, "output/mimic")
stan_dir <- file.path(project_root, "stan")
slim_dir <- file.path(stan_dir, "slim")
for (d in c(out_dir, slim_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)


# ============================================================================
# 1. DATA
# ============================================================================
source("mimic_load.R")   # → dat, pp_data, patients_df

if (is.null(pp_data)) pp_data <- expand_person_period(dat)

cat(sprintf("[data] N=%d  K=%d  P=%d  R=%d  events=%d (%.1f%%)\n",
            dat$N, dat$K, dat$P, pp_data$R, sum(dat$delta), 100*mean(dat$delta)))

beta_names   <- colnames(dat$X)
X_covariates <- pp_data$X.ast[, (dat$K + 1):(dat$K + dat$P)]

# Censored-terminal indicator for MNAR model:
#   1 at the terminal row of each censored patient, 0 elsewhere
terminal_pos       <- cumsum(dat$y_obs)
censored_terminal  <- rep(0L, pp_data$R)
censored_terminal[terminal_pos[dat$delta == 0]] <- 1L


# ============================================================================
# 2. COMPILE STAN MODELS
# ============================================================================

# Strip generated-quantities block → slim model (avoids serializing 510K
# log_lik values per iteration, which causes SIGPIPE crashes)
make_slim <- function(src, dst) {
  lines <- readLines(src)
  gq <- grep("^\\s*generated\\s+quantities\\s*\\{", lines)
  writeLines(if (length(gq)) lines[1:(gq[1]-1)] else lines, dst)
}

make_slim(file.path(stan_dir, "standard.stan"), file.path(slim_dir, "standard_slim.stan"))
make_slim(file.path(stan_dir, "frailty.stan"),  file.path(slim_dir, "frailty_slim.stan"))

cat("[stan] Compiling models...\n")
mod_A    <- cmdstan_model(file.path(slim_dir, "standard_slim.stan"))
mod_B    <- cmdstan_model(file.path(slim_dir, "frailty_slim.stan"))
mod_mnar <- cmdstan_model(file.path(stan_dir, "mnar_sensitivity.stan"))
cat("  All models compiled.\n")


# ============================================================================
# 3. FIT MODELS A & B
# ============================================================================

if (!LOAD_CHECKPOINT) {
  
  data_A <- list(R = pp_data$R, K = dat$K, P = dat$P,
                 y = as.integer(pp_data$y), X = X_covariates,
                 interval = as.integer(pp_data$interval))
  
  data_B <- c(data_A, list(N = dat$N, patient = as.integer(pp_data$patient_id)))
  
  cat("\n[fit] Model A (standard)...\n")
  fit_A <- mod_A$sample(data = data_A, chains = STAN_CHAINS,
                        parallel_chains = STAN_CHAINS,
                        iter_warmup = STAN_WARMUP,
                        iter_sampling = STAN_ITER,
                        adapt_delta = STAN_ADAPT,
                        max_treedepth = STAN_TREE, refresh = 500)
  
  cat("\n[fit] Model B (frailty)...\n")
  fit_B <- mod_B$sample(data = data_B, chains = STAN_CHAINS,
                        parallel_chains = STAN_CHAINS,
                        iter_warmup = STAN_WARMUP,
                        iter_sampling = STAN_ITER,
                        adapt_delta = STAN_ADAPT,
                        max_treedepth = STAN_TREE, refresh = 500)
  
  # ---- Diagnostics ----------------------------------------------------------
  cat("\n[diag] Model A\n"); print(fit_A$summary("beta"))
  cat("\n[diag] Model B\n"); print(fit_B$summary(c("sigma_b", "beta")))
  
  sb <- fit_B$summary("sigma_b")
  cat(sprintf("\nsigma_b: mean=%.3f  95%% CrI=[%.3f, %.3f]\n", sb$mean, sb$q5, sb$q95))
  
  # ---- Extract key quantities -----------------------------------------------
  draws_A   <- fit_A$draws()
  draws_B   <- fit_B$draws()
  alpha_hat <- fit_A$summary("alpha")$mean
  beta_A    <- fit_A$summary("beta")$mean
  beta_B    <- fit_B$summary("beta")$mean
  names(beta_A) <- names(beta_B) <- beta_names
  
  cmp_tbl <- data.frame(covariate = beta_names,
                        Model_A = round(beta_A, 3),
                        Model_B = round(beta_B, 3),
                        difference = round(beta_B - beta_A, 3))
  cat("\n[results] Beta comparison:\n"); print(cmp_tbl, row.names = FALSE)
  
  # ---- Calibration (Model A) ------------------------------------------------
  eta_mat      <- outer(rep(1, dat$N), alpha_hat) + as.vector(dat$X %*% beta_A)
  h_mat        <- plogis(eta_mat)
  log_surv_mat <- t(apply(log(1 - h_mat), 1, cumsum))
  surv_mat     <- exp(log_surv_mat)
  
  Ki_vec       <- pmin(dat$y_obs, dat$K)
  pred_mort    <- 1 - surv_mat[cbind(1:dat$N, Ki_vec)]
  risk_decile  <- cut(pred_mort,
                      breaks = quantile(pred_mort, seq(0, 1, 0.1)),
                      include.lowest = TRUE, labels = 1:10)
  cal_summary  <- aggregate(cbind(predicted = pred_mort, observed = dat$delta)
                            ~ as.numeric(risk_decile), FUN = mean)
  names(cal_summary)[1] <- "decile"
  
  # ---- KM + model survival --------------------------------------------------
  surv_df <- data.frame(time = (1:dat$K) * dat$delta_width,
                        model_surv = colMeans(surv_mat),
                        model_lo   = apply(surv_mat, 2, quantile, 0.025),
                        model_hi   = apply(surv_mat, 2, quantile, 0.975))
  
  km_fit     <- survfit(Surv(y_obs * dat$delta_width, delta) ~ 1,
                        data = data.frame(y_obs = dat$y_obs, delta = dat$delta))
  km_summary <- summary(km_fit, times = surv_df$time)
  km_df      <- data.frame(time = km_summary$time, km_surv = km_summary$surv,
                           km_lo = km_summary$lower, km_hi = km_summary$upper)
  surv_merged <- merge(surv_df, km_df, by = "time", all.x = TRUE)
  
  # ---- Save checkpoint ------------------------------------------------------
  save(dat, pp_data, patients_df, draws_A, draws_B,
       cmp_tbl, cal_summary, surv_merged, alpha_hat, beta_A, beta_B, beta_names,
       STAN_CHAINS, STAN_ITER, STAN_WARMUP,
       file = file.path(out_dir, "mimic_analysis_checkpoint.RData"),
       compress = FALSE)
  cat("[checkpoint] Saved: mimic_analysis_checkpoint.RData\n")
  
  rm(eta_mat, h_mat, log_surv_mat, surv_mat, pred_mort, surv_df, km_df, km_fit)
  gc()
  
} else {
  # ---- Resume from checkpoint -----------------------------------------------
  cat("[checkpoint] Loading mimic_analysis_checkpoint.RData...\n")
  load(file.path(out_dir, "mimic_analysis_checkpoint.RData"))
  X_covariates      <- pp_data$X.ast[, (dat$K + 1):(dat$K + dat$P)]
  beta_names        <- colnames(dat$X)
  if (!exists("beta_A")) { beta_A <- cmp_tbl$Model_A; names(beta_A) <- beta_names }
  terminal_pos      <- cumsum(dat$y_obs)
  censored_terminal <- rep(0L, pp_data$R)
  censored_terminal[terminal_pos[dat$delta == 0]] <- 1L
  rm(list = intersect(ls(), c("draws_A", "draws_B"))); gc()
  cat("  Checkpoint loaded.\n")
}


# ============================================================================
# 4. MNAR SENSITIVITY (ADVI)
# ============================================================================

if (!LOAD_MNAR) {
  
  cat("\n[mnar] Fitting MNAR model via ADVI (fullrank)...\n")
  SIGMA_GAMMA_GRID <- c(0.5, 1.0, 2.0)
  
  data_mnar_base <- list(R = pp_data$R, K = dat$K, P = dat$P,
                         y = as.integer(pp_data$y), X = X_covariates,
                         interval = as.integer(pp_data$interval),
                         censored_terminal = censored_terminal)
  
  mnar_results <- list()
  
  for (sg in SIGMA_GAMMA_GRID) {
    cat(sprintf("  sigma_gamma = %.1f ... ", sg))
    t0 <- Sys.time()
    
    fit_mnar <- mod_mnar$variational(
      data = c(data_mnar_base, list(sigma_gamma = sg)),
      algorithm = "fullrank", iter = 30000,
      tol_rel_obj = 0.001, output_samples = 2000
    )
    
    gamma_draws <- as.vector(fit_mnar$draws("gamma_sens", format = "matrix"))
    beta_mnar   <- fit_mnar$draws("beta", format = "matrix")
    colnames(beta_mnar) <- beta_names
    
    mnar_results[[as.character(sg)]] <- list(
      sigma_gamma   = sg,
      gamma_draws   = gamma_draws,
      gamma_summary = c(mean = mean(gamma_draws), sd = sd(gamma_draws),
                        q025 = unname(quantile(gamma_draws, 0.025)),
                        q975 = unname(quantile(gamma_draws, 0.975))),
      beta = beta_mnar
    )
    
    dt <- round(difftime(Sys.time(), t0, units = "secs"))
    gs <- mnar_results[[as.character(sg)]]$gamma_summary
    cat(sprintf("done (%ss)  gamma=%.3f [%.3f, %.3f]\n",
                dt, gs["mean"], gs["q025"], gs["q975"]))
  }
  
  save(mnar_results, beta_names, beta_A,
       file = file.path(out_dir, "mimic_mnar_results.RData"), compress = FALSE)
  cat("[checkpoint] Saved: mimic_mnar_results.RData\n")
  
} else {
  cat("[checkpoint] Loading mimic_mnar_results.RData...\n")
  load(file.path(out_dir, "mimic_mnar_results.RData"))
  cat(sprintf("  Loaded %d MNAR fits.\n", length(mnar_results)))
}

# ---- Print summary ----------------------------------------------------------
cat("\n=== MNAR Summary ===\n")
for (sg_key in names(mnar_results)) {
  gs <- mnar_results[[sg_key]]$gamma_summary
  cat(sprintf("  sg=%-3s  gamma=%.3f (sd=%.3f)  95%% CrI=[%.3f, %.3f]\n",
              sg_key, gs["mean"], gs["sd"], gs["q025"], gs["q975"]))
}

if ("1" %in% names(mnar_results)) {
  shifts <- colMeans(mnar_results[["1"]]$beta) - beta_A
  cat("\n=== Beta shifts (MAR -> MNAR, sg=1.0) ===\n")
  for (j in order(abs(shifts), decreasing = TRUE))
    cat(sprintf("  %-20s  %+.3f\n", beta_names[j], shifts[j]))
}

cat("\n[mimic_analysis.R complete]\n")

