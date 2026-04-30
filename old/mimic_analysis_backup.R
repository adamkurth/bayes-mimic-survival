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
#   8.  [OPTIONAL] Fixed sensitivity sweep (gamma grid, point-mass priors)
#   8.5 MNAR model (gamma estimated with proper prior) — Selection factorization
#   8.6 MNAR tables (LaTeX)
#   9.  Save results (all plotting handled by mimic_plots.R)
#----------------------------------------------------------------

# =============================================================================
# 1. SETUP & DATA LOADING
# =============================================================================

rm(list = ls())

# ---- USER CONFIG ----
USE_CMDSTANR       <- TRUE    # TRUE = cmdstanr (faster); FALSE = rstan fallback
SENS_USE_VB        <- TRUE    # TRUE = variational Bayes for sensitivity (fast)
LOAD_CHECKPOINT    <- TRUE    # TRUE = Skip Sections 3-7 and load from saved RData
LOAD_MNAR          <- TRUE    # TRUE = Skip Section 8.5 and load saved MNAR results
RUN_FIXED_SWEEP    <- FALSE   # TRUE = run fixed-gamma grid (Section 8); FALSE = skip
# (MNAR model in Section 8.5 replaces this)
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

mod_mnar <- compile_model(file.path(stan_dir, "mnar_sensitivity.stan"))
cat("  mnar_sensitivity.stan compiled\n")


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

# --- Censored-terminal indicator (for MNAR model, Section 8.5) ---
# Binary vector: 1 at the terminal row of each censored patient, 0 elsewhere.
# This is the discrete-time analog of the missingness indicator M_i from
# the selection factorization (Note Set 13): gamma * censored_terminal[r]
# shifts the hazard at the last observed interval for censored patients.
terminal_positions <- cumsum(dat$y_obs)
censored_mask      <- dat$delta == 0
censored_terminal  <- rep(0, pp_data$R)
censored_terminal[terminal_positions[censored_mask]] <- 1
cat(sprintf("[data] censored_terminal: %d of %d rows flagged (%.1f%%)\n",
            sum(censored_terminal), pp_data$R, 100 * mean(censored_terminal)))

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

cat(sprintf("  Calibration computed: %d deciles\n", nrow(cal_summary)))
# NOTE: Plotting moved to mimic_plots.R



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

cat("  Survival curves computed.\n")
# NOTE: Plotting moved to mimic_plots.R


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
       cmp_tbl, cal_summary, surv_merged, alpha_hat,
       STAN_CHAINS, STAN_ITER, STAN_WARMUP,
       file = file.path(out_dir, "mimic_analysis_checkpoint.RData"),
       compress = FALSE)
} else {
  save(dat, pp_data, patients_df, fit_A, fit_B,
       cmp_tbl, cal_summary, surv_merged, alpha_hat,
       STAN_CHAINS, STAN_ITER, STAN_WARMUP,
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
    
    # Re-create objects needed for Sections 8+
    X_covariates <- pp_data$X.ast[, (dat$K + 1):(dat$K + dat$P)]
    beta_names   <- colnames(dat$X)
    
    # Reconstruct beta_A and beta_B from the saved comparison table
    beta_A <- cmp_tbl$Model_A
    beta_B <- cmp_tbl$Model_B
    names(beta_A) <- names(beta_B) <- beta_names
    
    # Re-create censored_terminal indicator (needed for MNAR model, Section 8.5)
    terminal_positions <- cumsum(dat$y_obs)
    censored_mask      <- dat$delta == 0
    censored_terminal  <- rep(0, pp_data$R)
    censored_terminal[terminal_positions[censored_mask]] <- 1
    
    # Drop the massive draws objects to free up RAM
    rm(list = intersect(ls(), c("draws_A", "draws_B")))
    gc()
    
    cat("  Checkpoint loaded successfully! Ready for Section 8.\n")
    
    # --- Optionally load saved MNAR results (skip Section 8.5) ---
    if (LOAD_MNAR) {
      mnar_chk_path <- file.path(out_dir, "mimic_mnar_results.RData")
      if (file.exists(mnar_chk_path)) {
        cat("[checkpoint] Loading saved MNAR results...\n")
        load(mnar_chk_path)   # loads: mnar_results, beta_names, beta_A
        cat(sprintf("  Loaded %d MNAR fits from mimic_mnar_results.RData\n",
                    length(mnar_results)))
      } else {
        cat("  No MNAR checkpoint found — will run Section 8.5 from scratch.\n")
        LOAD_MNAR <- FALSE
      }
    }
    
  } else {
    stop("Checkpoint file not found! You must run Sections 1-7 first with LOAD_CHECKPOINT = FALSE.")
  }
}


# =============================================================================
# 8. SENSITIVITY SWEEP — Fixed gamma grid (OPTIONAL)
# =============================================================================
# NOTE: The MNAR model in Section 8.5 subsumes this analysis by estimating
# gamma with a proper prior. Set RUN_FIXED_SWEEP = TRUE to also run the
# point-mass grid sweep (tilted.stan). Otherwise, skip to Section 8.5.
# =============================================================================

# Ensure beta_names exists (may be missing if running from checkpoint)
if (!exists("beta_names")) beta_names <- colnames(dat$X)

if (RUN_FIXED_SWEEP) {
  
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
  
  # NOTE: Sensitivity plotting moved to mimic_plots.R
  
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
  
} else {
  cat("\n[skip] Section 8: Fixed sweep disabled (RUN_FIXED_SWEEP = FALSE).\n")
  cat("  Proceeding directly to MNAR model (Section 8.5).\n")
  # Initialize empty objects so downstream save doesn't fail
  stan_sens <- list()
  sens_all <- data.frame()
  tipping_points <- data.frame(parameter = character(), gamma_tip = numeric(),
                               stringsAsFactors = FALSE)
}  # end if (RUN_FIXED_SWEEP)


# =============================================================================
# 8.5 MNAR MODEL — JOINTLY ESTIMATE GAMMA (Selection Factorization)
# =============================================================================
# This section implements the MNAR approach:
#
#   f(m, y, x | omega) = f(m | y, x; xi) * f(y | x; eta) * f(x; theta)
#
# Instead of fixing gamma at each grid point (point-mass prior, Section 8),
# we place a proper prior gamma ~ N(0, sigma_gamma) and let Stan estimate
# the posterior f(gamma | D^o). This parallels mnar_mod.stan from the lecture
# where eps2 ~ N(0, 1) encodes prior uncertainty about MNAR departure.
#
# sigma_gamma controls prior beliefs:
#   - sigma_gamma = 1: moderate uncertainty around MAR
#   - sigma_gamma = 0.5: tighter prior favoring MAR
#   - sigma_gamma = 2: diffuse prior allowing large departures
#
# The posterior of gamma tells us what the data + prior jointly support
# about the degree of informative censoring.
# =============================================================================

if (LOAD_MNAR && exists("mnar_results") && length(mnar_results) > 0) {
  cat("\n=== MNAR Model: Loaded from checkpoint (skipping VB fitting) ===\n")
  cat(sprintf("  %d sigma_gamma values: %s\n", length(mnar_results),
              paste(names(mnar_results), collapse = ", ")))
} else {
  
  cat("\n=== MNAR Model: Estimating gamma with N(0, sigma_gamma) prior ===\n")
  
  SIGMA_GAMMA_GRID <- c(0.5, 1.0, 2.0)
  
  stan_data_mnar_base <- list(
    R                 = pp_data$R,
    K                 = dat$K,
    P                 = dat$P,
    y                 = as.integer(pp_data$y),
    X                 = X_covariates,
    interval          = as.integer(pp_data$interval),
    censored_terminal = censored_terminal
  )
  
  mnar_results <- list()
  
  for (sg in SIGMA_GAMMA_GRID) {
    cat(sprintf("  sigma_gamma = %.1f ... ", sg))
    t0 <- Sys.time()
    
    stan_data_mnar <- c(stan_data_mnar_base, list(sigma_gamma = sg))
    
    if (is_cmdstan_mod(mod_mnar)) {
      if (SENS_USE_VB) {
        fit_mnar <- tryCatch(
          mod_mnar$variational(
            data = stan_data_mnar, algorithm = "fullrank",
            iter = 30000, tol_rel_obj = 0.001, output_samples = 2000
          ),
          error = function(e) {
            cat(sprintf("VB failed (%s), trying MCMC... ", e$message))
            mod_mnar$sample(
              data = stan_data_mnar, chains = STAN_CHAINS,
              parallel_chains = STAN_CHAINS,
              iter_warmup = STAN_WARMUP,
              iter_sampling = STAN_ITER - STAN_WARMUP,
              adapt_delta = 0.90, max_treedepth = STAN_TREE, refresh = 0
            )
          }
        )
      } else {
        fit_mnar <- mod_mnar$sample(
          data = stan_data_mnar, chains = STAN_CHAINS,
          parallel_chains = STAN_CHAINS,
          iter_warmup = STAN_WARMUP,
          iter_sampling = STAN_ITER - STAN_WARMUP,
          adapt_delta = 0.90, max_treedepth = STAN_TREE, refresh = 500
        )
      }
      gamma_draws  <- as.vector(fit_mnar$draws("gamma_sens", format = "matrix"))
      beta_mnar    <- fit_mnar$draws("beta", format = "matrix")
    } else {
      if (SENS_USE_VB) {
        fit_mnar <- tryCatch(
          rstan::vb(mod_mnar, data = stan_data_mnar, algorithm = "fullrank",
                    iter = 30000, tol_rel_obj = 0.001, output_samples = 2000),
          error = function(e) {
            cat(sprintf("VB failed (%s), trying MCMC... ", e$message))
            rstan::sampling(mod_mnar, data = stan_data_mnar,
                            chains = STAN_CHAINS, iter = STAN_ITER, warmup = STAN_WARMUP,
                            control = list(adapt_delta = 0.90, max_treedepth = STAN_TREE),
                            refresh = 0)
          }
        )
      } else {
        fit_mnar <- rstan::sampling(
          mod_mnar, data = stan_data_mnar,
          chains = STAN_CHAINS, iter = STAN_ITER, warmup = STAN_WARMUP,
          control = list(adapt_delta = 0.90, max_treedepth = STAN_TREE),
          refresh = 500
        )
      }
      gamma_draws <- rstan::extract(fit_mnar, "gamma_sens")$gamma_sens
      beta_mnar   <- rstan::extract(fit_mnar, "beta")$beta
    }
    
    colnames(beta_mnar) <- beta_names
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    gamma_summary <- c(
      mean   = mean(gamma_draws),
      median = median(gamma_draws),
      sd     = sd(gamma_draws),
      q025   = unname(quantile(gamma_draws, 0.025)),
      q975   = unname(quantile(gamma_draws, 0.975))
    )
    
    cat(sprintf("done (%.0f sec)\n", dt))
    cat(sprintf("    gamma posterior: mean=%.3f, median=%.3f, 95%% CrI=[%.3f, %.3f]\n",
                gamma_summary["mean"], gamma_summary["median"],
                gamma_summary["q025"], gamma_summary["q975"]))
    
    mnar_results[[as.character(sg)]] <- list(
      sigma_gamma   = sg,
      gamma_draws   = gamma_draws,
      gamma_summary = gamma_summary,
      beta          = beta_mnar,
      time          = dt
    )
  }
  
  # --- Save MNAR checkpoint ---
  mnar_chk_path <- file.path(out_dir, "mimic_mnar_results.RData")
  cat(sprintf("\n[checkpoint] Saving MNAR results to %s ...\n", mnar_chk_path))
  save(mnar_results, beta_names, beta_A,
       file = mnar_chk_path, compress = FALSE)
  cat("  Saved: mimic_mnar_results.RData\n")
  
}  # end if/else LOAD_MNAR guard for Section 8.5

# --- Report MNAR results ---
cat("\n=== MNAR Model Summary ===\n")
cat(sprintf("%-15s %10s %10s %15s\n",
            "sigma_gamma", "gamma_mean", "gamma_sd", "95% CrI"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (sg_key in names(mnar_results)) {
  gs <- mnar_results[[sg_key]]$gamma_summary
  cat(sprintf("%-15s %10.3f %10.3f    [%.3f, %.3f]\n",
              sg_key, gs["mean"], gs["sd"], gs["q025"], gs["q975"]))
}

# --- Comparison: beta under MNAR vs MAR (Model A) ---
cat("\n=== Beta comparison: MAR (Model A) vs MNAR (sigma_gamma=1.0) ===\n")
if ("1" %in% names(mnar_results)) {
  beta_mnar_mean <- colMeans(mnar_results[["1"]]$beta)
  mnar_cmp <- data.frame(
    covariate = beta_names,
    MAR       = round(beta_A, 3),
    MNAR      = round(beta_mnar_mean, 3),
    shift     = round(beta_mnar_mean - beta_A, 3)
  )
  print(mnar_cmp, row.names = FALSE)
}


# =============================================================================
# 8.6 MNAR TABLES (plots moved to mimic_plots.R)
# =============================================================================
cat("\n[tables] All plotting is handled by mimic_plots.R\n")

# ---------- MNAR summary table (LaTeX) ----------
cat("\n[tables] Writing MNAR LaTeX tables...\n")
table_dir <- file.path(project_root, "output", "tables")
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

# Table: gamma posterior summary
gamma_tex <- file.path(table_dir, "mimic_mnar_gamma.tex")
cat("\\begin{table}[t]\n\\centering\n", file = gamma_tex)
cat("\\caption{Posterior summary of $\\gamma$ under three prior widths $\\sigma_\\gamma$.}\n", file = gamma_tex, append = TRUE)
cat("\\label{tab:mnar-gamma}\n", file = gamma_tex, append = TRUE)
cat("\\begin{tabular}{ccccc}\n\\toprule\n", file = gamma_tex, append = TRUE)
cat("$\\sigma_\\gamma$ & Mean & SD & 95\\% CrI \\\\\n\\midrule\n", file = gamma_tex, append = TRUE)

for (sg_key in names(mnar_results)) {
  gs <- mnar_results[[sg_key]]$gamma_summary
  cat(sprintf("%.1f & %.3f & %.3f & [%.3f, %.3f] \\\\\n",
              as.numeric(sg_key), gs["mean"], gs["sd"], gs["q025"], gs["q975"]),
      file = gamma_tex, append = TRUE)
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n", file = gamma_tex, append = TRUE)
cat(sprintf("  Saved: %s\n", gamma_tex))

# Table: beta comparison MAR vs MNAR
if ("1" %in% names(mnar_results)) {
  beta_mnar_mean <- colMeans(mnar_results[["1"]]$beta)
  beta_mnar_lo   <- apply(mnar_results[["1"]]$beta, 2, quantile, 0.025)
  beta_mnar_hi   <- apply(mnar_results[["1"]]$beta, 2, quantile, 0.975)
  
  mnar_beta_tex <- file.path(table_dir, "mimic_mnar_beta.tex")
  cat("\\begin{table}[t]\n\\centering\n\\small\n", file = mnar_beta_tex)
  cat("\\caption{Covariate effects under MAR (Model~A) vs.~MNAR ($\\sigma_\\gamma = 1$).}\n",
      file = mnar_beta_tex, append = TRUE)
  cat("\\label{tab:mnar-beta}\n", file = mnar_beta_tex, append = TRUE)
  cat("\\begin{tabular}{lcccc}\n\\toprule\n", file = mnar_beta_tex, append = TRUE)
  cat("Covariate & MAR $\\hat\\beta$ & MNAR $\\hat\\beta$ & 95\\% CrI (MNAR) & $\\Delta\\beta$ \\\\\n\\midrule\n",
      file = mnar_beta_tex, append = TRUE)
  
  for (j in seq_along(beta_names)) {
    shade <- if (j %% 2 == 0) "\\rowcolor{rowgray}" else ""
    cat(sprintf("%s %s & %.3f & %.3f & [%.3f, %.3f] & %s%.3f \\\\\n",
                shade, gsub("_", "\\\\_", beta_names[j]),
                beta_A[j], beta_mnar_mean[j],
                beta_mnar_lo[j], beta_mnar_hi[j],
                ifelse(beta_mnar_mean[j] - beta_A[j] < 0, "$-$", ""),
                abs(beta_mnar_mean[j] - beta_A[j])),
        file = mnar_beta_tex, append = TRUE)
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n", file = mnar_beta_tex, append = TRUE)
  cat(sprintf("  Saved: %s\n", mnar_beta_tex))
}


# =============================================================================
# 9. SAVE RESULTS (all plotting handled by mimic_plots.R)
# =============================================================================
cat("\n[save] Saving final results...\n")

# Save comprehensive fits file (works from checkpoint or full run)
save_objs <- c("dat", "pp_data", "cmp_tbl", "beta_names", "beta_A",
               "cal_summary", "surv_merged", "alpha_hat")
if (exists("patients_df"))   save_objs <- c(save_objs, "patients_df")
if (exists("draws_A"))       save_objs <- c(save_objs, "draws_A")
if (exists("draws_B"))       save_objs <- c(save_objs, "draws_B")
if (exists("mnar_results"))  save_objs <- c(save_objs, "mnar_results")
if (exists("stan_sens"))     save_objs <- c(save_objs, "stan_sens")
if (exists("sens_all"))      save_objs <- c(save_objs, "sens_all")
if (exists("tipping_points")) save_objs <- c(save_objs, "tipping_points")

save(list = save_objs,
     file = file.path(out_dir, "mimic_analysis_fits.RData"),
     compress = FALSE)
cat(sprintf("  Saved: %s/mimic_analysis_fits.RData\n", out_dir))

# Beta summary CSV (from checkpoint or full run)
if (exists("cmp_tbl")) {
  beta_export <- data.frame(
    covariate = beta_names,
    Model_A   = round(beta_A, 4),
    stringsAsFactors = FALSE
  )
  if ("Model_B" %in% names(cmp_tbl)) {
    beta_export$Model_B <- round(cmp_tbl$Model_B, 4)
  }
  write.csv(beta_export, file.path(out_dir, "mimic_beta_summary.csv"), row.names = FALSE)
  cat(sprintf("  Saved: %s/mimic_beta_summary.csv\n", out_dir))
}

cat("\n=== Outputs ===\n")
cat(sprintf("  %s/mimic_analysis_checkpoint.RData\n", out_dir))
cat(sprintf("  %s/mimic_mnar_results.RData\n", out_dir))
cat(sprintf("  %s/mimic_analysis_fits.RData\n", out_dir))
cat(sprintf("  %s/mimic_beta_summary.csv\n", out_dir))
cat("  NOTE: All plots are generated by mimic_plots.R\n")
cat("\n[mimic_analysis.R complete]\n")