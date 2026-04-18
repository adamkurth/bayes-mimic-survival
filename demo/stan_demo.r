#----------------------------------------------------------------
# stan_demo.r
# Bayesian discrete-time survival model with Stan (rstan) -- demo
#----------------------------------------------------------------
# Provides robust and fully specified Stan implementation that complements the PG Gibbs sampler in `core.r`.
# In particular, this fixes the Neal's-funnel mixing problem in the frailty variance by using non-centered
# parameterization and HMC and adds full diagnostics / model comparison

# inherit from core.r
#   This script SOURCES core.R to reuse:
#     - generate_synthetic_icu()
#     - expand_person_period()
#     - pg_gibbs_standard(), pg_gibbs_frailty(), pg_gibbs_tilted()
#     - make_offset_vector()
#

# Structure:
#   1.  Setup and data generation
#   2.  PG Gibbs fits (for comparison)
#   3.  Stan Model A: standard discrete-time hazard (no frailty)
#   4.  Stan Model B: shared frailty (non-centered, PC prior)
#   5.  Stan Model A-tilted: standard + exponential tilting (sensitivity)
#   6.  Stan Model D: 2-class latent mixture
#   7.  Compile and fit Models A, B, D
#   8.  Diagnostics: Rhat, ESS, divergences
#   9.  Risk-decile calibration plot (Model A)
#  10.  Model-implied vs empirical survival curves (Model A)
#  11.  LOO / WAIC model comparison (A vs B vs D)
#  12.  Sensitivity sweep (gamma) via A-tilted — all betas, tipping points
#  13.  Known-DGP validation (independent vs informative censoring)
#  14.  Discretization sensitivity study (DELTA = 1, 3, 7, 14)
#  15.  Combined comparison table
#  16.  Overview plots
#  17.  Save all fits

# 1. setup ----------------------------------------------------------------
rm(list = ls())
library(BayesLogit); library(coda); library(ggplot2); library(bayesplot); library(loo); library(gridExtra); 
library(posterior); library(rstan); library(survival)

rstan_options(auto_write = TRUE)         # cache compiled models
options(mc.cores = parallel::detectCores())
bayesplot::color_scheme_set("brightblue")
set.seed(2530)

# source
setwd("~/Documents/RStudio/bayes-mimic-survival/")
source("demo/core.r") # inherit 


# 2. data generation & gibbs fits -------------------------------------------
N_PATIENTS <- 300   # Reduced from 3000 for testing; increase after debugging
DELTA      <- 1
K_MAX      <- ceiling(90 / DELTA)   # 180-day (~ 6 month) observation window

dat     <- generate_synthetic_icu(N = N_PATIENTS, K = K_MAX, delta = DELTA)
pp_data <- expand_person_period(dat)

cat(sprintf("\n[data] N=%d, K=%d, R=%d person-periods, %d events (%.1f%%)\n",
            dat$N, dat$K, pp_data$R, sum(dat$delta), 100*mean(dat$delta)))

# PG Gibbs fits for side-by-side comparison
fit_pg_std     <- pg_gibbs_standard(pp_data, n_iter = 1000, burn_in = 300)
fit_pg_frailty <- pg_gibbs_frailty(pp_data, dat, n_iter = 1000, burn_in = 300)

# 3. Stan Model A: discrete-time hazard, PRIMARY MODEL ------------------------

stan_code_standard <- "
data {
  int<lower=1> R;
  int<lower=1> K;
  int<lower=1> P;
  array[R] int<lower=0, upper=1> y;
  matrix[R, P] X;
  array[R] int<lower=1, upper=K> interval;
}
parameters {
  vector[K] alpha;
  vector[P] beta;
  real<lower=0> sigma_alpha;
}
model {
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);
  beta ~ normal(0, sqrt(10));

  {
    vector[R] eta;
    for (r in 1:R) eta[r] = alpha[interval[r]] + X[r] * beta;
    y ~ bernoulli_logit(eta);
  }
}
generated quantities {
  vector[R] log_lik;
  array[R] int y_rep;
  for (r in 1:R) {
    real eta_r = alpha[interval[r]] + X[r] * beta;
    log_lik[r] = bernoulli_logit_lpmf(y[r] | eta_r);
    y_rep[r]   = bernoulli_logit_rng(eta_r);
  }
}
"


# 4. Stan Model B: shared frailty (non-centered, PC prior) ------------------------
# Retained as simulation-study evidence. Demonstrates that continuous
# patient-level frailty is weakly identified from sparse binary data.
# Key features:
#   - Non-centered: b_raw ~ N(0,1), b = sigma_b * b_raw (flattens funnel)
#   - Penalized complexity prior: sigma_b ~ exponential(1)
#     (Simpson et al. 2017 — shrinks toward the null of no heterogeneity)
stan_code_frailty <- "
data {
  int<lower=1> R;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  array[R] int<lower=0, upper=1> y;
  matrix[R, P] X;
  array[R] int<lower=1, upper=K> interval;
  array[R] int<lower=1, upper=N> patient;
}
parameters {
  vector[K] alpha;
  vector[P] beta;
  vector[N] b_raw;
  real<lower=0> sigma_b;
  real<lower=0> sigma_alpha;
}
transformed parameters {
  vector[N] b = sigma_b * b_raw;   // NON-CENTERED
}
model {
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);

  beta    ~ normal(0, sqrt(10));
  b_raw   ~ std_normal();
  sigma_b ~ exponential(1);        // PC prior: shrinks toward sigma_b = 0

  {
    vector[R] eta;
    for (r in 1:R) eta[r] = alpha[interval[r]] + X[r] * beta + b[patient[r]];
    y ~ bernoulli_logit(eta);
  }
}
generated quantities {
  vector[R] log_lik;
  array[R] int y_rep;
  for (r in 1:R) {
    real eta_r = alpha[interval[r]] + X[r] * beta + b[patient[r]];
    log_lik[r] = bernoulli_logit_lpmf(y[r] | eta_r);
    y_rep[r]   = bernoulli_logit_rng(eta_r);
  }
}
"

# Also keep the half-t prior version for prior sensitivity comparison
stan_code_frailty_halft <- gsub(
  "sigma_b ~ exponential(1);        // PC prior: shrinks toward sigma_b = 0",
  "sigma_b ~ student_t(3, 0, 1);    // half-t prior (weaker regularization)",
  stan_code_frailty, fixed = TRUE
)


# 5. Stan Model C: sensitivity with STANDARD MODEL + exponential tilting ------------------------
# Same as frailty model but with a data-supplied offset vector o[r] that
# encodes gamma * I(C_i = t, delta_i = 0) at the terminal censored interval.

stan_code_tilted_std <- "
data {
  int<lower=1> R;
  int<lower=1> K;
  int<lower=1> P;
  array[R] int<lower=0, upper=1> y;
  matrix[R, P] X;
  array[R] int<lower=1, upper=K> interval;
  vector[R] offset_vec;
}
parameters {
  vector[K] alpha;
  vector[P] beta;
  real<lower=0> sigma_alpha;
}
model {
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);
  beta ~ normal(0, sqrt(10));

  {
    vector[R] eta;
    for (r in 1:R) eta[r] = alpha[interval[r]] + X[r] * beta + offset_vec[r];
    y ~ bernoulli_logit(eta);
  }
}
generated quantities {
  vector[P] beta_out = beta;
}
"

# 6. 2D latent mixture (continuous frailty)------------------------------------------------ 
# Parsimonious alternative to continuous frailty. Captures unobserved patient
# heterogeneity with 2 class-specific intercepts + 1 mixing proportion
# (3 extra parameters total, vs N for continuous frailty).
# Uses ordered class effects for identifiability and log_sum_exp for
# numerically stable marginalization over classes.
stan_code_latent_class <- "
data {
  int<lower=1> R;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  array[R] int<lower=0, upper=1> y;
  matrix[R, P] X;
  array[R] int<lower=1, upper=K> interval;
  array[R] int<lower=1, upper=N> patient;
  // Per-patient indexing for efficient marginalization
  array[N] int<lower=1, upper=R> pat_start;
  array[N] int<lower=1> pat_len;
}
parameters {
  vector[K] alpha;
  vector[P] beta;
  ordered[2] class_effect;      // ordered for identifiability
  real<lower=0, upper=1> pi_1;  // P(class 1) — low-risk class
  real<lower=0> sigma_alpha;
}
model {
  // Priors
  alpha[1] ~ normal(0, 3);
  for (t in 2:K) alpha[t] ~ normal(alpha[t-1], sigma_alpha);
  sigma_alpha ~ student_t(3, 0, 1);
  beta ~ normal(0, sqrt(10));
  class_effect ~ normal(0, 2);
  pi_1 ~ beta(2, 2);           // weakly informative, centered at 0.5

  // Pre-compute base linear predictor (no class effect)
  vector[R] eta_base;
  for (r in 1:R)
    eta_base[r] = alpha[interval[r]] + X[r] * beta;

  // Marginalize over 2 classes for each patient
  for (i in 1:N) {
    real lp1 = log(pi_1);
    real lp2 = log1m(pi_1);
    for (j in 0:(pat_len[i] - 1)) {
      int r = pat_start[i] + j;
      lp1 += bernoulli_logit_lpmf(y[r] | eta_base[r] + class_effect[1]);
      lp2 += bernoulli_logit_lpmf(y[r] | eta_base[r] + class_effect[2]);
    }
    target += log_sum_exp(lp1, lp2);
  }
}
generated quantities {
  // Observation-level marginal log-likelihood (marginalized over classes)
  vector[R] log_lik;
  {
    vector[R] eta_base;
    for (r in 1:R)
      eta_base[r] = alpha[interval[r]] + X[r] * beta;

    for (r in 1:R) {
      real lp1 = log(pi_1)    + bernoulli_logit_lpmf(y[r] | eta_base[r] + class_effect[1]);
      real lp2 = log1m(pi_1)  + bernoulli_logit_lpmf(y[r] | eta_base[r] + class_effect[2]);
      log_lik[r] = log_sum_exp(lp1, lp2);
    }
  }

  // Posterior class probabilities for each patient (useful for diagnostics)
  vector[N] prob_class2;
  {
    vector[R] eta_base2;
    for (r in 1:R)
      eta_base2[r] = alpha[interval[r]] + X[r] * beta;

    for (i in 1:N) {
      real lp1 = log(pi_1);
      real lp2 = log1m(pi_1);
      for (j in 0:(pat_len[i] - 1)) {
        int r = pat_start[i] + j;
        lp1 += bernoulli_logit_lpmf(y[r] | eta_base2[r] + class_effect[1]);
        lp2 += bernoulli_logit_lpmf(y[r] | eta_base2[r] + class_effect[2]);
      }
      prob_class2[i] = inv_logit(lp2 - lp1);  // P(class 2 | data)
    }
  }
}
"




# 7.  COMPILE AND FIT MODELS A, B, D ------------------------------------------------

# --- Stan data ---
stan_data_base <- list(
  R        = pp_data$R,
  N        = dat$N,
  K        = dat$K,
  P        = dat$P,
  y        = as.integer(pp_data$y),
  X        = pp_data$X.ast[, (dat$K+1):(dat$K+dat$P)],
  interval = as.integer(pp_data$interval),
  patient  = as.integer(pp_data$patient_id)
)

# Per-patient indexing for latent class model
pat_len   <- as.integer(table(pp_data$patient_id))
pat_start <- as.integer(c(1, cumsum(pat_len[-length(pat_len)]) + 1))
stan_data_lc <- c(stan_data_base, list(pat_start = pat_start, pat_len = pat_len))

# Shared sampling settings — robust configuration
STAN_CTRL   <- list(adapt_delta = 0.99, max_treedepth = 12)
STAN_CHAINS <- 4
STAN_ITER   <- 4000
STAN_WARMUP <- 1500

# --- Model A (standard) ---
cat("\n[stan] Compiling and fitting Model A (standard)...\n")
mod_standard <- stan_model(model_code = stan_code_standard, verbose = FALSE)
fit_stan_std <- sampling(mod_standard,
                         data    = stan_data_base[c("R","K","P","y","X","interval")],
                         chains  = STAN_CHAINS,
                         iter    = STAN_ITER,
                         warmup  = STAN_WARMUP,
                         control = STAN_CTRL)

# --- Model B (frailty, PC prior) ---
cat("\n[stan] Compiling and fitting Model B (frailty, PC prior)...\n")
mod_frailty <- stan_model(model_code = stan_code_frailty, verbose = FALSE)
fit_stan_frailty <- sampling(mod_frailty,
                             data    = stan_data_base,
                             chains  = STAN_CHAINS,
                             iter    = STAN_ITER,
                             warmup  = STAN_WARMUP,
                             control = STAN_CTRL)

# --- Model B (frailty, half-t prior) for prior sensitivity ---
cat("\n[stan] Fitting Model B with half-t prior (prior sensitivity)...\n")
mod_frailty_ht <- stan_model(model_code = stan_code_frailty_halft, verbose = FALSE)
fit_stan_frailty_ht <- sampling(mod_frailty_ht,
                                data    = stan_data_base,
                                chains  = STAN_CHAINS,
                                iter    = STAN_ITER,
                                warmup  = STAN_WARMUP,
                                control = STAN_CTRL)

# --- Model D (2-class latent mixture) ---
cat("\n[stan] Compiling and fitting Model D (latent class)...\n")
mod_latent_class <- stan_model(model_code = stan_code_latent_class, verbose = FALSE)
fit_stan_lc <- sampling(mod_latent_class,
                        data    = stan_data_lc,
                        chains  = STAN_CHAINS,
                        iter    = STAN_ITER,
                        warmup  = STAN_WARMUP,
                        control = STAN_CTRL)



# 8. diagnostics and comparison (Rhat, ESS, divergences) ----------------------------------

report_diagnostics <- function(fit, label, pars = c("sigma_b","sigma_alpha","beta")) {
  cat(sprintf("\n=== Diagnostics: %s ===\n", label))
  print(summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary)

  sp    <- get_sampler_params(fit, inc_warmup = FALSE)
  n_div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  n_tot <- STAN_CHAINS * (STAN_ITER - STAN_WARMUP)
  cat(sprintf("\n  Divergences: %d / %d post-warmup\n", n_div, n_tot))

  if ("sigma_b" %in% pars) {
    ess  <- summary(fit, pars = "sigma_b")$summary[, "n_eff"]
    rhat <- summary(fit, pars = "sigma_b")$summary[, "Rhat"]
    cat(sprintf("  sigma_b: ESS = %.1f, Rhat = %.4f\n", ess, rhat))
  }
}

# Model A
report_diagnostics(fit_stan_std, "Model A (standard)",
                   pars = c("sigma_alpha", "beta"))

# Model B (both priors)
report_diagnostics(fit_stan_frailty, "Model B (frailty, PC prior)")
report_diagnostics(fit_stan_frailty_ht, "Model B (frailty, half-t prior)")

# Model D
report_diagnostics(fit_stan_lc, "Model D (latent class)",
                   pars = c("sigma_alpha", "beta", "class_effect", "pi_1"))

# PG vs Stan ESS comparison for sigma_b
pg_ess      <- coda::effectiveSize(as.mcmc(fit_pg_frailty$sigma_b2))
stan_ess_pc <- summary(fit_stan_frailty, pars = "sigma_b")$summary[, "n_eff"]
stan_ess_ht <- summary(fit_stan_frailty_ht, pars = "sigma_b")$summary[, "n_eff"]

cat("\n--- sigma_b ESS comparison ---\n")
cat(sprintf("  PG Gibbs (centered, InvGamma)    ESS = %7.1f  (of 2000)\n", pg_ess))
cat(sprintf("  Stan NUTS (non-centered, PC)     ESS = %7.1f  (of %d)\n",
            stan_ess_pc, STAN_CHAINS * (STAN_ITER - STAN_WARMUP)))
cat(sprintf("  Stan NUTS (non-centered, half-t) ESS = %7.1f  (of %d)\n",
            stan_ess_ht, STAN_CHAINS * (STAN_ITER - STAN_WARMUP)))

# Prior sensitivity for sigma_b
sb_pc <- summary(fit_stan_frailty, pars = "sigma_b")$summary
sb_ht <- summary(fit_stan_frailty_ht, pars = "sigma_b")$summary
cat("\n--- sigma_b prior sensitivity ---\n")
cat(sprintf("  PC prior (exp(1)):  mean = %.3f [%.3f, %.3f]\n",
            sb_pc[,"mean"], sb_pc[,"2.5%"], sb_pc[,"97.5%"]))
cat(sprintf("  half-t(3,0,1):      mean = %.3f [%.3f, %.3f]\n",
            sb_ht[,"mean"], sb_ht[,"2.5%"], sb_ht[,"97.5%"]))
cat(sprintf("  Truth:              sigma_b = %.3f\n", dat$sigma_b_true))

# Latent class summary
lc_class <- summary(fit_stan_lc, pars = c("class_effect", "pi_1"))$summary
cat("\n--- Latent class model ---\n")
cat(sprintf("  Class 1 (low-risk)  effect = %.3f [%.3f, %.3f]\n",
            lc_class["class_effect[1]","mean"],
            lc_class["class_effect[1]","2.5%"],
            lc_class["class_effect[1]","97.5%"]))
cat(sprintf("  Class 2 (high-risk) effect = %.3f [%.3f, %.3f]\n",
            lc_class["class_effect[2]","mean"],
            lc_class["class_effect[2]","2.5%"],
            lc_class["class_effect[2]","97.5%"]))
cat(sprintf("  P(class 1) = %.3f [%.3f, %.3f]\n",
            lc_class["pi_1","mean"],
            lc_class["pi_1","2.5%"],
            lc_class["pi_1","97.5%"]))
cat(sprintf("  Implied sigma_b (approx) = %.3f\n",
            abs(lc_class["class_effect[2]","mean"] - lc_class["class_effect[1]","mean"]) / 2))


# 9. risk-decile calibration plot  (MODEL A) -------------------------------------------------------
# Using the STANDARD model (well-identified) for calibration.
# This validates the event model without making claims about censoring.

cat("\n[calibration] Computing risk-decile plot (Model A)...\n")

# Get posterior mean beta and alpha from standard model
beta_hat  <- summary(fit_stan_std, pars = "beta")$summary[, "mean"]
alpha_hat <- summary(fit_stan_std, pars = "alpha")$summary[, "mean"]

# Patient-level predicted mortality: P(event before K) = 1 - prod(1 - h_t)
pred_mort <- numeric(dat$N)
for (i in 1:dat$N) {
  eta_i  <- alpha_hat + as.numeric(dat$X[i, ] %*% beta_hat)
  h_i    <- plogis(eta_i)
  pred_mort[i] <- 1 - prod(1 - h_i[1:min(dat$y_obs[i], dat$K)])
}

obs_mort <- dat$delta

# Create deciles
risk_decile <- cut(pred_mort,
                   breaks = quantile(pred_mort, probs = seq(0, 1, 0.1)),
                   include.lowest = TRUE, labels = 1:10)
cal_df <- data.frame(decile    = as.numeric(risk_decile),
                     predicted = pred_mort,
                     observed  = obs_mort)

cal_summary <- aggregate(cbind(predicted, observed) ~ decile, data = cal_df, FUN = mean)

p_calibration <- ggplot(cal_summary, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 1) +
  geom_line(color = "steelblue", alpha = 0.4, linewidth = 1) +
  geom_point(size = 3.5, color = "steelblue") +
  geom_point(size = 1.5, color = "white") + # Creates a nice "donut" effect
  geom_text(aes(label = decile), vjust = -1.2, size = 3.5, fontface = "bold", color = "grey30") +
  labs(x = "Mean predicted mortality (per decile)",
       y = "Observed mortality (per decile)",
       title = "Risk-Decile Calibration (Model A)",
       subtitle = "Numbers indicate risk decile; dashed line represents perfect calibration") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

ggsave("stan_demo_calibration.pdf", p_calibration, width = 6, height = 6)




# 10. model-implied vs KM curves (MODEL A)-------------------------------------------------------


cat("\n[survival] Computing model-implied vs KM curves (Model A)...\n")

# Empirical Kaplan-Meier
km_fit <- survfit(Surv(y_obs, delta) ~ 1,
                  data = data.frame(y_obs = dat$y_obs, delta = dat$delta))

# Model-implied marginal survival: S(t) = E_x[ prod_{s=1}^{t} (1 - h_s(x)) ]
surv_model <- matrix(0, nrow = dat$N, ncol = dat$K)
for (i in 1:dat$N) {
  eta_i <- alpha_hat + as.numeric(dat$X[i, ] %*% beta_hat)
  h_i   <- plogis(eta_i)
  surv_model[i, ] <- cumprod(1 - h_i)
}
surv_mean <- colMeans(surv_model)

# Build comparison
surv_df <- data.frame(time = 1:dat$K, model_surv = surv_mean)
km_summary <- summary(km_fit, times = 1:dat$K)
km_df <- data.frame(time = km_summary$time, km_surv = km_summary$surv,
                    km_lo = km_summary$lower, km_hi = km_summary$upper)
surv_merged <- merge(surv_df, km_df, by = "time", all.x = TRUE)

p_survival <- ggplot(surv_merged, aes(x = time)) +
  geom_ribbon(aes(ymin = km_lo, ymax = km_hi), fill = "grey70", alpha = 0.4) +
  geom_step(aes(y = km_surv, color = "Kaplan-Meier"), linewidth = 1.2) +
  geom_line(aes(y = model_surv, color = "Stan Model A"), linewidth = 1.2) +
  scale_color_manual(values = c("Kaplan-Meier" = "black",
                                "Stan Model A" = "#D55E00")) +
  labs(x = "Interval (days)", y = "Survival Probability S(t)",
       title = "Model-Implied vs Empirical Survival (Model A)",
       subtitle = "Grey ribbon represents 95% CI for Kaplan-Meier estimates",
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = c(0.85, 0.85), # Moves legend inside plot
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

ggsave("stan_demo_survival.pdf", p_survival, width = 7, height = 5)




# 11. LOO / WAIC MODEL COMPARISON (A vs B vs D) ----------------------------------

# what is loo and waic? how to compute from log_lik? how to interpret?
# loo: leave-one-out cross-validation estimate of out-of-sample predictive fit, 
#   using Pareto-smoothed importance sampling to approximate the LOO posteriors without refitting
# waic: widely applicable information criterion, an approximation to out-of-sample predictive fit
#    based on the log pointwise predictive density, with a penalty for effective number of parameters

cat("\n=== LOO model comparison ===\n")

extract_loglik <- function(fit, par = "log_lik") rstan::extract(fit, par)[[1]]

ll_std     <- extract_loglik(fit_stan_std)
ll_frailty <- extract_loglik(fit_stan_frailty)
ll_lc      <- extract_loglik(fit_stan_lc)

loo_std     <- loo::loo(ll_std, cores = 2)
loo_frailty <- loo::loo(ll_frailty, cores = 2)
loo_lc      <- loo::loo(ll_lc, cores = 2)

cat("\n-- Model A (standard) --\n"); print(loo_std)
cat("\n-- Model B (frailty, PC prior) --\n"); print(loo_frailty)
cat("\n-- Model D (latent class) --\n"); print(loo_lc)

cat("\n-- Pairwise comparison: A vs B --\n")
print(loo::loo_compare(loo_std, loo_frailty))
cat("\n-- Pairwise comparison: A vs D --\n")
print(loo::loo_compare(loo_std, loo_lc))
cat("\n-- Three-way comparison --\n")
print(loo::loo_compare(loo_std, loo_frailty, loo_lc))

waic_std     <- loo::waic(ll_std)
waic_frailty <- loo::waic(ll_frailty)
waic_lc      <- loo::waic(ll_lc)
cat(sprintf("\nWAIC Model A = %.2f (SE %.2f)\n",
            waic_std$estimates["waic","Estimate"], waic_std$estimates["waic","SE"]))
cat(sprintf("WAIC Model B = %.2f (SE %.2f)\n",
            waic_frailty$estimates["waic","Estimate"], waic_frailty$estimates["waic","SE"]))
cat(sprintf("WAIC Model D = %.2f (SE %.2f)\n",
            waic_lc$estimates["waic","Estimate"], waic_lc$estimates["waic","SE"]))


# # 12. SENSITIVITY SWEEP (gamma) — STANDARD MODEL + TILTING in Stan --------------------------------------------------
# runs off Model A (no frailty). This is the centerpiece analysis.
# The tilting parameter gamma encodes departure from independent censoring
# (Scharfstein, Rotnitzky & Robins 1999). gamma = 0 corresponds to MAR.

cat("\n=== Sensitivity sweep (gamma) via Model A-tilted ===\n")

mod_tilted_std <- stan_model(model_code = stan_code_tilted_std, verbose = FALSE)

# Grid: finer near zero, coarser at extremes
gamma_grid <- c(-2.0, -1.5, -1.0, -0.5, -0.25, 0, 0.25, 0.5, 1.0, 1.5, 2.0)

# Sensitivity sweep settings — well-identified model, converges fast
SENS_CHAINS <- 4
SENS_ITER   <- 2500
SENS_WARMUP <- 1000

stan_sens <- list()

for (g in gamma_grid) {
  cat(sprintf("  gamma = %+5.2f ... ", g))
  offset_vec  <- make_offset_vector(pp_data, dat, gamma_val = g)

  # Data for tilted standard model (no patient ID needed)
  stan_data_g <- list(
    R          = pp_data$R,
    K          = dat$K,
    P          = dat$P,
    y          = as.integer(pp_data$y),
    X          = pp_data$X.ast[, (dat$K+1):(dat$K+dat$P)],
    interval   = as.integer(pp_data$interval),
    offset_vec = offset_vec
  )

  fit_g <- sampling(mod_tilted_std, data = stan_data_g,
                    chains  = SENS_CHAINS,
                    iter    = SENS_ITER,
                    warmup  = SENS_WARMUP,
                    control = STAN_CTRL,
                    refresh = 0)

  beta_g <- rstan::extract(fit_g, "beta")$beta
  colnames(beta_g) <- colnames(pp_data$X.ast)[(dat$K+1):(dat$K+dat$P)]

  sp_g    <- get_sampler_params(fit_g, inc_warmup = FALSE)
  n_div_g <- sum(sapply(sp_g, function(x) sum(x[, "divergent__"])))

  # ESS for beta (no sigma_b to worry about)
  beta_ess <- summary(fit_g, pars = "beta")$summary[, "n_eff"]

  stan_sens[[as.character(g)]] <- list(
    beta    = beta_g,
    n_div   = n_div_g,
    ess_beta_min = min(beta_ess)
  )
  cat(sprintf("done (div=%d, min_ESS_beta=%.0f)\n", n_div_g, min(beta_ess)))
}

# --- Build sensitivity data frame for ALL betas ---
beta_names <- colnames(pp_data$X.ast)[(dat$K+1):(dat$K+dat$P)]
sens_all <- do.call(rbind, lapply(names(stan_sens), function(g) {
  bg <- stan_sens[[g]]$beta
  do.call(rbind, lapply(seq_along(beta_names), function(j) {
    data.frame(
      gamma     = as.numeric(g),
      parameter = beta_names[j],
      mean      = mean(bg[, j]),
      lo        = quantile(bg[, j], 0.025),
      hi        = quantile(bg[, j], 0.975),
      truth     = dat$beta_true[j],
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(sens_all) <- NULL

# --- Sensitivity plot: faceted by parameter ---
p_sens_all <- ggplot(sens_all, aes(x = gamma)) +
  geom_hline(aes(yintercept = 0, linetype = "Null Effect"), color = "grey60", linewidth = 0.8) +
  geom_hline(aes(yintercept = truth, color = "Truth"), linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "95% CrI"), alpha = 0.2) +
  geom_line(aes(y = mean, color = "Stan Model A (Tilted)"), linewidth = 1) +
  geom_point(aes(y = mean, color = "Stan Model A (Tilted)"), size = 1.5) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  scale_color_manual(name = "Estimate", 
                     values = c("Truth" = "#D55E00", "Stan Model A (Tilted)" = "steelblue")) +
  scale_fill_manual(name = "Uncertainty", values = c("95% CrI" = "steelblue")) +
  scale_linetype_manual(name = "Reference", values = c("Null Effect" = "solid")) +
  labs(x = expression(gamma ~ "(Departure from MAR)"), 
       y = "Posterior Mean",
       title = "Sensitivity Analysis: Covariate Effects Under Informative Censoring") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("stan_demo_sensitivity_all.pdf", p_sens_all, width = 10, height = 8)


# --- Tipping points: gamma at which 95% CrI first includes 0 ---
cat("\n--- Tipping points (gamma where 95% CrI crosses zero) ---\n")
cat(sprintf("%-15s %10s %15s\n", "Parameter", "gamma_tip", "Note"))

for (param in beta_names) {
  sub <- sens_all[sens_all$parameter == param, ]
  contains_zero <- sub$lo <= 0 & sub$hi >= 0
  if (any(contains_zero)) {
    first_cross <- sub$gamma[which(contains_zero)[1]]
    cat(sprintf("%-15s %10.2f  (interval includes 0)\n", param, first_cross))
  } else {
    cat(sprintf("%-15s %10s  (never crosses 0 in grid)\n", param, "---"))
  }
}


# 13. known-DGP validation (indep vs. informative) --------------------------------------------------
# Now that core.r is FIXED (cens_frailty_coef is actually used), this
# comparison should show a meaningful contrast.

cat("\n=== Known-DGP validation ===\n")

set.seed(2530)  # same seed for comparable cohort structure

cat("[DGP] Generating data with INDEPENDENT censoring (cens_frailty_coef = 0)...\n")
dat_indep <- generate_synthetic_icu(N = N_PATIENTS, K = K_MAX, delta = DELTA,
                                    cens_frailty_coef = 0)
pp_indep  <- expand_person_period(dat_indep)

cat(sprintf("  Independent: N=%d, R=%d, events=%d (%.1f%%)\n",
            dat_indep$N, pp_indep$R, sum(dat_indep$delta), 100*mean(dat_indep$delta)))
cat(sprintf("  Informative: N=%d, R=%d, events=%d (%.1f%%)\n",
            dat$N, pp_data$R, sum(dat$delta), 100*mean(dat$delta)))

stan_data_indep <- list(
  R = pp_indep$R, N = dat_indep$N, K = dat_indep$K, P = dat_indep$P,
  y = as.integer(pp_indep$y),
  X = pp_indep$X.ast[, (dat_indep$K+1):(dat_indep$K+dat_indep$P)],
  interval = as.integer(pp_indep$interval),
  patient  = as.integer(pp_indep$patient_id)
)

cat("[DGP] Fitting standard + frailty under independent censoring...\n")
fit_indep_std <- sampling(mod_standard,
                          data    = stan_data_indep[c("R","K","P","y","X","interval")],
                          chains  = STAN_CHAINS, iter = STAN_ITER,
                          warmup  = STAN_WARMUP, control = STAN_CTRL, refresh = 200)
fit_indep_fra <- sampling(mod_frailty,
                          data    = stan_data_indep,
                          chains  = STAN_CHAINS, iter = STAN_ITER,
                          warmup  = STAN_WARMUP, control = STAN_CTRL, refresh = 200)

# sigma_b recovery comparison
sb_indep  <- summary(fit_indep_fra, pars = "sigma_b")$summary
sb_inform <- summary(fit_stan_frailty, pars = "sigma_b")$summary

cat("\n--- sigma_b recovery: independent vs informative censoring ---\n")
cat(sprintf("  Independent: sigma_b = %.3f [%.3f, %.3f]  (truth = %.3f)\n",
            sb_indep[,"mean"], sb_indep[,"2.5%"], sb_indep[,"97.5%"], dat_indep$sigma_b_true))
cat(sprintf("  Informative: sigma_b = %.3f [%.3f, %.3f]  (truth = %.3f)\n",
            sb_inform[,"mean"], sb_inform[,"2.5%"], sb_inform[,"97.5%"], dat$sigma_b_true))

# Beta comparison: independent censoring should recover truth better
beta_indep_std  <- summary(fit_indep_std, pars = "beta")$summary[, "mean"]
beta_indep_fra  <- summary(fit_indep_fra, pars = "beta")$summary[, "mean"]
beta_inform_std <- summary(fit_stan_std, pars = "beta")$summary[, "mean"]
beta_inform_fra <- summary(fit_stan_frailty, pars = "beta")$summary[, "mean"]

dgp_tbl <- data.frame(
  parameter       = names(dat$beta_true),
  truth           = round(as.numeric(dat$beta_true), 3),
  indep_std       = round(beta_indep_std, 3),
  indep_frailty   = round(beta_indep_fra, 3),
  inform_std      = round(beta_inform_std, 3),
  inform_frailty  = round(beta_inform_fra, 3)
)
cat("\n--- Beta recovery: independent vs informative censoring ---\n")
print(dgp_tbl, row.names = FALSE)

# LOO comparison under each censoring regime
ll_indep_std <- extract_loglik(fit_indep_std)
ll_indep_fra <- extract_loglik(fit_indep_fra)
loo_indep_std <- loo::loo(ll_indep_std, cores = 2)
loo_indep_fra <- loo::loo(ll_indep_fra, cores = 2)

cat("\n--- LOO under independent censoring ---\n")
print(loo::loo_compare(loo_indep_std, loo_indep_fra))

cat("\n--- LOO under informative censoring ---\n")
print(loo::loo_compare(loo_std, loo_frailty))


# 14. DISCRETIZATION SENSITIVITY STUDY -------------------------------------------------------
# Focus on BETA STABILITY from Model A (which is well-identified).
# Secondarily shows sigma_b instability from Model B as a cautionary note.

cat("\n=== Discretization sensitivity study ===\n")

delta_grid   <- c(1, 3, 7, 14)
disc_results <- list()

for (d in delta_grid) {
  cat(sprintf("\n[disc] DELTA = %d ...\n", d))

  set.seed(2530)
  K_d   <- ceiling(180 / d)   # match the main analysis window
  dat_d <- generate_synthetic_icu(N = N_PATIENTS, K = K_d, delta = d)
  pp_d  <- expand_person_period(dat_d)

  sd_d <- list(
    R = pp_d$R, N = dat_d$N, K = dat_d$K, P = dat_d$P,
    y = as.integer(pp_d$y),
    X = pp_d$X.ast[, (dat_d$K+1):(dat_d$K+dat_d$P)],
    interval = as.integer(pp_d$interval),
    patient  = as.integer(pp_d$patient_id)
  )

  fit_d_std <- sampling(mod_standard,
                        data = sd_d[c("R","K","P","y","X","interval")],
                        chains = 4, iter = 3000, warmup = 1000,
                        control = STAN_CTRL, refresh = 0)

  fit_d_fra <- sampling(mod_frailty,
                        data = sd_d,
                        chains = 4, iter = 3000, warmup = 1000,
                        control = STAN_CTRL, refresh = 0)

  ll_d_std <- extract_loglik(fit_d_std)
  ll_d_fra <- extract_loglik(fit_d_fra)
  loo_d_std <- loo::loo(ll_d_std, cores = 2)
  loo_d_fra <- loo::loo(ll_d_fra, cores = 2)

  sp_d <- get_sampler_params(fit_d_fra, inc_warmup = FALSE)

  # Standard model beta CIs for stability assessment
  beta_std_summary <- summary(fit_d_std, pars = "beta")$summary
  beta_std_lo <- beta_std_summary[, "2.5%"]
  beta_std_hi <- beta_std_summary[, "97.5%"]

  disc_results[[as.character(d)]] <- list(
    delta       = d,
    K           = K_d,
    R           = pp_d$R,
    n_events    = sum(dat_d$delta),
    beta_std    = beta_std_summary[, "mean"],
    beta_std_lo = beta_std_lo,
    beta_std_hi = beta_std_hi,
    beta_fra    = summary(fit_d_fra, pars = "beta")$summary[, "mean"],
    sigma_b     = summary(fit_d_fra, pars = "sigma_b")$summary[, "mean"],
    ess_sigma_b = summary(fit_d_fra, pars = "sigma_b")$summary[, "n_eff"],
    rhat_sigma_b = summary(fit_d_fra, pars = "sigma_b")$summary[, "Rhat"],
    n_div       = sum(sapply(sp_d, function(x) sum(x[, "divergent__"]))),
    loo_std     = loo_d_std$estimates["looic", "Estimate"],
    loo_fra     = loo_d_fra$estimates["looic", "Estimate"],
    p_loo_fra   = loo_d_fra$estimates["p_loo", "Estimate"]
  )
}

# --- Discretization summary table ---
cat("\n--- Discretization sensitivity summary ---\n")
cat(sprintf("%-6s %4s %7s %5s %8s %7s %7s %5s %9s %9s %9s\n",
            "DELTA", "K", "R", "evts", "sigma_b", "ESS_sb", "Rhat",
            "div", "LOOIC_std", "LOOIC_fra", "p_loo_fra"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (d in as.character(delta_grid)) {
  r <- disc_results[[d]]
  cat(sprintf("%-6d %4d %7d %5d %8.3f %7.1f %7.4f %5d %9.1f %9.1f %9.1f\n",
              r$delta, r$K, r$R, r$n_events, r$sigma_b,
              r$ess_sigma_b, r$rhat_sigma_b, r$n_div,
              r$loo_std, r$loo_fra, r$p_loo_fra))
}

# --- Beta stability across discretizations (standard model) ---
disc_beta_df <- do.call(rbind, lapply(names(disc_results), function(d) {
  r <- disc_results[[d]]
  data.frame(
    delta     = r$delta,
    parameter = names(dat$beta_true),
    truth     = as.numeric(dat$beta_true),
    std       = r$beta_std,
    std_lo    = r$beta_std_lo,
    std_hi    = r$beta_std_hi,
    frailty   = r$beta_fra,
    stringsAsFactors = FALSE
  )
}))
rownames(disc_beta_df) <- NULL

# --- Beta stability across discretizations (standard model) ---
p_disc <- ggplot(disc_beta_df, aes(x = factor(delta))) +
  geom_hline(aes(yintercept = truth, linetype = "Truth"), color = "darkred", linewidth = 0.8) +
  geom_pointrange(aes(y = std, ymin = std_lo, ymax = std_hi, 
                      color = "Standard", shape = "Standard"),
                  linewidth = 0.6, position = position_dodge(width = 0.4)) +
  geom_point(aes(y = frailty, color = "Frailty", shape = "Frailty"), 
             size = 2.5, position = position_dodge(width = 0.4)) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  scale_color_manual(name = "Model", values = c("Standard" = "grey40", "Frailty" = "steelblue")) +
  scale_shape_manual(name = "Model", values = c("Standard" = 16, "Frailty" = 17)) +
  scale_linetype_manual(name = "Reference", values = c("Truth" = "dashed")) +
  labs(x = expression(Delta ~ "(Interval width in days)"), 
       y = "Posterior Mean",
       title = "Discretization Sensitivity: Coefficient Stability") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("stan_demo_discretization.pdf", p_disc, width = 10, height = 8)

# 15. combined model comparison table -------------------------------------------------------

stan_beta_summary <- function(fit, P) {
  as.numeric(summary(fit, pars = "beta")$summary[, "mean"])
}

cmp_tbl <- data.frame(
  parameter    = names(dat$beta_true),
  truth        = round(as.numeric(dat$beta_true), 3),
  pg_standard  = round(colMeans(fit_pg_std$beta), 3),
  pg_frailty   = round(colMeans(fit_pg_frailty$beta), 3),
  stan_std     = round(stan_beta_summary(fit_stan_std, dat$P), 3),
  stan_fra_pc  = round(stan_beta_summary(fit_stan_frailty, dat$P), 3),
  stan_fra_ht  = round(stan_beta_summary(fit_stan_frailty_ht, dat$P), 3),
  stan_lc      = round(stan_beta_summary(fit_stan_lc, dat$P), 3)
)

cat("\n=== Combined posterior-mean comparison table ===\n")
print(cmp_tbl, row.names = FALSE)

pg_sigma_b <- sqrt(mean(fit_pg_frailty$sigma_b2))
stan_sb_pc <- summary(fit_stan_frailty, pars = "sigma_b")$summary[, "mean"]
stan_sb_ht <- summary(fit_stan_frailty_ht, pars = "sigma_b")$summary[, "mean"]

cat(sprintf("\nsigma_b: truth=%.3f | PG=%.3f | Stan(PC)=%.3f | Stan(half-t)=%.3f\n",
            dat$sigma_b_true, pg_sigma_b, stan_sb_pc, stan_sb_ht))
cat(sprintf("ESS(sigma_b): PG=%.0f | Stan(PC)=%.0f | Stan(half-t)=%.0f\n",
            pg_ess, stan_ess_pc, stan_ess_ht))

write.csv(cmp_tbl, "stan_vs_pg_comparison.csv", row.names = FALSE)

# 12. final plots -----------------------------------------------------------

pdf("stan_demo_overview.pdf", width = 12, height = 16)
par(mfrow = c(4, 2), mar = c(4.5, 4.5, 3, 1), cex.main = 1.2, family = "sans")

# --- Panel 1: Baseline hazard ---
K <- dat$K
alpha_stan <- summary(fit_stan_std, pars = "alpha")$summary[, "mean"]

plot(1:K, plogis(dat$alpha_true), type = "l", lwd = 2.5, col = "black",
     xlab = "Interval (days)", ylab = "h(t)", main = "Baseline Hazard",
     ylim = range(plogis(dat$alpha_true), plogis(alpha_stan)) * c(0.5, 1.3),
     bty = "l", las = 1)
grid(col = "grey90", lty = 1)
lines(1:K, plogis(alpha_stan), col = "steelblue", lwd = 2.5, lty = 2)
legend("topleft", c("True", "Stan Model A"), col = c("black","steelblue"), 
       lty = 1:2, lwd = 2.5, cex = 0.85, bty = "n")

# --- Panel 2: Stan Model A Trace (Log-Posterior) ---
lp_stan_A <- rstan::extract(fit_stan_std, "lp__")$lp__
plot(lp_stan_A, type = "l", col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "Iteration", ylab = "Log-Posterior",
     main = "Stan Model A Trace: Log-Posterior (lp__)",
     bty = "l", las = 1)
grid(col = "grey90", lty = 1)

# --- Panel 3: Stan Model A Trace (Beta 1) ---
beta1_stan_A <- rstan::extract(fit_stan_std, "beta")$beta[, 1]
plot(beta1_stan_A, type = "l", col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "Iteration", ylab = names(dat$beta_true)[1],
     main = paste("Stan Model A Trace: Coefficient (", names(dat$beta_true)[1], ")", sep=""),
     bty = "l", las = 1)
grid(col = "grey90", lty = 1)
abline(h = dat$beta_true[1], lty = 2, lwd = 2, col = "#D55E00")

# --- Panel 4: Beta comparison (Profile Plot instead of Jitter) ---
beta_pg   <- colMeans(fit_pg_frailty$beta)
beta_stan <- stan_beta_summary(fit_stan_std, dat$P)
beta_lc   <- stan_beta_summary(fit_stan_lc, dat$P)
idx <- seq_along(dat$beta_true)

plot(idx, dat$beta_true, type = "b", pch = 16, cex = 1.5, lwd = 2, col = "black",
     xaxt = "n", xlab = "", ylab = "Coefficient Value",
     main = "Covariate Effects: Truth vs Models (Profile)",
     ylim = range(c(dat$beta_true, beta_pg, beta_stan, beta_lc)) * 1.2,
     bty = "l", las = 1)
grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
abline(h = 0, lty = 3, col = "grey50")

# Plot each model connecting the dots. Much easier to read than jitter.
lines(idx, beta_pg,   type = "b", pch = 17, col = "#D55E00",  cex = 1.3, lwd = 1.5, lty = 2)
lines(idx, beta_stan, type = "b", pch = 15, col = "steelblue", cex = 1.3, lwd = 1.5, lty = 3)
lines(idx, beta_lc,   type = "b", pch = 18, col = "#009E73",  cex = 1.5, lwd = 1.5, lty = 4)

axis(1, idx, labels = names(dat$beta_true), las = 2, cex.axis = 0.8)
legend("bottomleft", c("True", "PG Frailty", "Stan A (Std)", "Stan D (LC)"),
       pch = c(16, 17, 15, 18), col = c("black", "#D55E00", "steelblue", "#009E73"), 
       lty = c(1, 2, 3, 4), lwd = 1.5, cex = 0.85, bty = "n")

# --- Panels 5-7: Sensitivity sweep (smooth polygons) ---
key_params <- c("emergency", "icd_sepsis", "vasopressor")
for (param in key_params) {
  sub <- sens_all[sens_all$parameter == param, ]
  plot(sub$gamma, sub$mean, type = "n", 
       xlab = expression(gamma), ylab = "Posterior Mean",
       main = paste0("Sensitivity: ", param),
       ylim = range(c(sub$lo, sub$hi, sub$truth[1])),
       bty = "l", las = 1)
  grid(col = "grey90", lty = 1)
  
  polygon(c(sub$gamma, rev(sub$gamma)), c(sub$lo, rev(sub$hi)), 
          col = rgb(0.27, 0.51, 0.71, 0.25), border = NA)
  
  abline(h = sub$truth[1], lty = 2, lwd = 1.5, col = "#D55E00")
  abline(h = 0, lty = 3, col = "grey50")
  lines(sub$gamma, sub$mean, type = "b", pch = 16, lwd = 2, col = "steelblue")
}

# --- Panel 8: Discretization sigma_b ---
disc_sb <- sapply(disc_results, function(r) r$sigma_b)
disc_d  <- sapply(disc_results, function(r) r$delta)

plot(disc_d, disc_sb, type = "b", pch = 16, lwd = 2, col = "steelblue",
     xlab = expression(Delta ~ "(Interval width in days)"), ylab = expression(hat(sigma)[b]),
     main = expression("Discretization:" ~ hat(sigma)[b] ~ "vs Width"),
     log = "x", bty = "l", las = 1)
grid(col = "grey90", lty = 1, equilogs = FALSE)
abline(h = dat$sigma_b_true, lty = 2, lwd = 1.5, col = "#D55E00")
text(disc_d, disc_sb, labels = paste0("K=", sapply(disc_results, function(r) r$K)),
     pos = 3, cex = 0.85, col = "grey30")

dev.off()


# save fits for reuse (e.g. in the paper) --------------------------------------------
save(dat, pp_data, dat_indep, pp_indep,
     fit_pg_std, fit_pg_frailty,
     fit_stan_std, fit_stan_frailty, fit_stan_frailty_ht, fit_stan_lc,
     fit_indep_std, fit_indep_fra,
     stan_sens, disc_results,
     loo_std, loo_frailty, loo_lc, loo_indep_std, loo_indep_fra,
     cmp_tbl, dgp_tbl, sens_all, disc_beta_df,
     file = "stan_demo_fits.RData")

cat("\n=== Outputs ===\n")
cat("  stan_demo_calibration.pdf       risk-decile calibration (Model A)\n")
cat("  stan_demo_survival.pdf          model-implied vs KM survival (Model A)\n")
cat("  stan_demo_sensitivity_all.pdf   sensitivity sweep, all betas (A-tilted)\n")
cat("  stan_demo_discretization.pdf    beta stability across DELTA\n")
cat("  stan_demo_overview.pdf          8-panel overview\n")
cat("  stan_vs_pg_comparison.csv       posterior-mean table (all models)\n")
cat("  stan_demo_fits.RData            all fits for reuse\n")
cat("\n[stan_demo.r complete]\n")
