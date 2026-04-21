# mimic_summary.R
# ==============================================================================
# MIMIC-III Cohort Summary Statistics for Paper Appendix
# ==============================================================================
# Produces:
#   1. Table 1: Cohort characteristics stratified by outcome
#   2. Table 2: Observation time & data structure
#   3. Table 3: Covariate matrix description (page-safe)
#   4. Table 4: Censoring pattern by interval
#   5. LaTeX code for all tables → mimic_tables.tex
#   6. PDF figure: follow-up time distribution
#
# OUTCOME DECOMPOSITION:
#   Within the K_MAX observation window (K_MAX * DELTA days):
#     - Died:       event_status == 1 AND time_to_event <= window
#     - Discharged:  event_status == 0 AND time_to_event < window
#                   (clinically censored — the informative mechanism)
#     - Survived:    time_to_event >= window
#                   (administratively censored at K_MAX — non-informative)
#
# PREREQUISITES:
#   source("mimic_load.R")   # loads dat, patients_df
# ==============================================================================

if (!exists("dat") || !exists("patients_df")) {
  stop("Run source('mimic_load.R') first to load dat and patients_df.")
}

library(dplyr)
library(tidyr)

DELTA  <- dat$delta_width
K      <- dat$K
N      <- dat$N
P      <- dat$P
X      <- dat$X
WINDOW <- K * DELTA   # observation window in days

out.path <- file.path(
  "/Users/adamkurth/Documents/RStudio/bayes-mimic-survival", "output"
)
if (!dir.exists(out.path)) dir.create(out.path, recursive = TRUE)

# ==============================================================================
# OUTCOME DECOMPOSITION
# ==============================================================================
# This is the key distinction for the paper:
#   1. Died within window       — the events
#   2. Discharged alive < window — clinically censored (potentially informative)
#   3. Survived past window     — administratively censored (non-informative)

patients_df <- patients_df %>%
  mutate(
    outcome_group = case_when(
      delta == 1                                      ~ "Died",
      event_status == 0 & time_to_event_days < WINDOW ~ "Discharged",
      TRUE                                            ~ "Survived"
    ),
    outcome_group = factor(outcome_group,
                           levels = c("Died", "Discharged", "Survived"))
  )

grp <- patients_df %>%
  group_by(outcome_group) %>%
  summarize(n = n(), .groups = "drop") %>%
  mutate(pct = 100 * n / N)

n_died        <- grp$n[grp$outcome_group == "Died"]
n_discharged  <- grp$n[grp$outcome_group == "Discharged"]
n_survived    <- grp$n[grp$outcome_group == "Survived"]

# Subsets for stratified summaries
died_df       <- patients_df %>% filter(outcome_group == "Died")
discharged_df <- patients_df %>% filter(outcome_group == "Discharged")
survived_df   <- patients_df %>% filter(outcome_group == "Survived")

# For the model, "censored" = discharged + survived (delta == 0)
n_censored_model <- n_discharged + n_survived

message(sprintf("\n  Outcome decomposition within %d-day window:", WINDOW))
message(sprintf("    Died within window:          %s (%.1f%%)",
                formatC(n_died, big.mark = ","), 100 * n_died / N))
message(sprintf("    Discharged alive (< %d days): %s (%.1f%%)",
                WINDOW, formatC(n_discharged, big.mark = ","), 100 * n_discharged / N))
message(sprintf("    Survived past %d days:        %s (%.1f%%)",
                WINDOW, formatC(n_survived, big.mark = ","), 100 * n_survived / N))
message(sprintf("    Total censored (model):      %s (%.1f%%)",
                formatC(n_censored_model, big.mark = ","), 100 * n_censored_model / N))

# ==============================================================================
# HELPERS
# ==============================================================================
fmt_n_pct <- function(n, total) {
  sprintf("%s (%.1f\\%%)", formatC(n, format = "d", big.mark = ","), 100 * n / total)
}
fmt_median_iqr <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}
fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f $\\pm$ %.1f", mean(x), sd(x))
}

# ==============================================================================
# TABLE 1: COHORT CHARACTERISTICS BY OUTCOME (3-way split)
# ==============================================================================
message("\n==== Generating Table 1 ====")

# Helper to pull a binary sum from a subset
bin_sum <- function(df, var) sum(df[[var]], na.rm = TRUE)

tex1 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  sprintf("\\caption{Baseline characteristics of the MIMIC-III analytic cohort within the %d-day observation window, stratified by outcome. \\emph{Discharged} denotes patients censored alive before the window closes (the clinically informative censoring mechanism); \\emph{Survived} denotes patients still alive at %d days (administratively censored). Values are $n$ (\\%%) for categorical variables and mean $\\pm$ SD or median [IQR] for continuous variables.}", WINDOW, WINDOW),
  "\\label{tab:cohort}",
  "\\begin{tabular}{l C{2.4cm} C{2.0cm} C{2.2cm} C{2.0cm}}",
  "\\toprule",
  sprintf(" & \\textbf{Overall} & \\textbf{Died} & \\textbf{Discharged} & \\textbf{Survived} \\\\"),
  sprintf(" & (N=%s) & (n=%s) & (n=%s) & (n=%s) \\\\",
          formatC(N, big.mark = ","),
          formatC(n_died, big.mark = ","),
          formatC(n_discharged, big.mark = ","),
          formatC(n_survived, big.mark = ",")),
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Demographics}} \\\\[2pt]"
)

# --- Demographics ---
tex1 <- c(tex1, sprintf("Age, years (mean $\\pm$ SD) & %s & %s & %s & %s \\\\",
                        fmt_mean_sd(patients_df$age), fmt_mean_sd(died_df$age),
                        fmt_mean_sd(discharged_df$age), fmt_mean_sd(survived_df$age)))

tex1 <- c(tex1, sprintf("Male sex & %s & %s & %s & %s \\\\",
                        fmt_n_pct(bin_sum(patients_df, "male"), N),
                        fmt_n_pct(bin_sum(died_df, "male"), n_died),
                        fmt_n_pct(bin_sum(discharged_df, "male"), n_discharged),
                        fmt_n_pct(bin_sum(survived_df, "male"), n_survived)))

tex1 <- c(tex1, sprintf("Emergency admission & %s & %s & %s & %s \\\\",
                        fmt_n_pct(bin_sum(patients_df, "emergency"), N),
                        fmt_n_pct(bin_sum(died_df, "emergency"), n_died),
                        fmt_n_pct(bin_sum(discharged_df, "emergency"), n_discharged),
                        fmt_n_pct(bin_sum(survived_df, "emergency"), n_survived)))

tex1 <- c(tex1, sprintf("Follow-up, days (median [IQR]) & %s & %s & %s & %s \\\\",
                        fmt_median_iqr(patients_df$time_to_event_days),
                        fmt_median_iqr(died_df$time_to_event_days),
                        fmt_median_iqr(discharged_df$time_to_event_days),
                        fmt_median_iqr(survived_df$time_to_event_days)))

# --- Comorbidities ---
tex1 <- c(tex1,
          "\\midrule",
          "\\multicolumn{5}{l}{\\textit{Comorbidities (ICD-9)}} \\\\[2pt]")

icd_map <- c(
  icd_sepsis     = "Sepsis",
  icd_hf         = "Heart failure",
  icd_aki        = "Acute kidney injury",
  icd_resp       = "Respiratory failure",
  icd_pneumonia  = "Pneumonia",
  icd_diabetes   = "Diabetes",
  icd_liver      = "Liver disease",
  icd_malignancy = "Malignancy",
  icd_copd       = "COPD"
)

for (var in names(icd_map)) {
  tex1 <- c(tex1, sprintf("\\quad %s & %s & %s & %s & %s \\\\",
                          icd_map[var],
                          fmt_n_pct(bin_sum(patients_df, var), N),
                          fmt_n_pct(bin_sum(died_df, var), n_died),
                          fmt_n_pct(bin_sum(discharged_df, var), n_discharged),
                          fmt_n_pct(bin_sum(survived_df, var), n_survived)))
}

# --- Clinical indicators ---
tex1 <- c(tex1,
          "\\midrule",
          "\\multicolumn{5}{l}{\\textit{Clinical Indicators}} \\\\[2pt]")

clin_map <- c(
  gcs_impaired      = "GCS impaired ($\\leq 12$)",
  is_ventilated     = "Mechanical ventilation",
  on_vasopressors   = "Vasopressor use",
  received_dialysis = "Dialysis",
  positive_blood_cx = "Positive blood culture"
)

for (var in names(clin_map)) {
  tex1 <- c(tex1, sprintf("\\quad %s & %s & %s & %s & %s \\\\",
                          clin_map[var],
                          fmt_n_pct(bin_sum(patients_df, var), N),
                          fmt_n_pct(bin_sum(died_df, var), n_died),
                          fmt_n_pct(bin_sum(discharged_df, var), n_discharged),
                          fmt_n_pct(bin_sum(survived_df, var), n_survived)))
}

tex1 <- c(tex1,
          "\\bottomrule",
          "\\end{tabular}",
          "\\end{table}")


# ==============================================================================
# TABLE 2: DATA STRUCTURE & OBSERVATION WINDOW
# ==============================================================================
message("==== Generating Table 2 ====")

pp_rows <- sum(dat$y_obs)

tex2 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  sprintf("\\caption{Summary of the discrete-time data structure within the %d-day observation window ($\\Delta = %d$ days, $K = %d$ intervals). Clinically censored patients are those discharged alive before %d days; administratively censored patients survived the full window.}",
          WINDOW, DELTA, K, WINDOW),
  "\\label{tab:data-structure}",
  "\\begin{tabular}{l r}",
  "\\toprule",
  sprintf("Unique admissions ($N$) & %s \\\\", formatC(N, big.mark = ",")),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Outcome within observation window}} \\\\[2pt]",
  sprintf("\\quad Died within %d days & %s (%.1f\\%%) \\\\",
          WINDOW, formatC(n_died, big.mark = ","), 100 * n_died / N),
  sprintf("\\quad Discharged alive $<$ %d days & %s (%.1f\\%%) \\\\",
          WINDOW, formatC(n_discharged, big.mark = ","), 100 * n_discharged / N),
  sprintf("\\quad Survived $\\geq$ %d days & %s (%.1f\\%%) \\\\",
          WINDOW, formatC(n_survived, big.mark = ","), 100 * n_survived / N),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Discrete-time structure}} \\\\[2pt]",
  sprintf("\\quad Discretization width ($\\Delta$) & %d days \\\\", DELTA),
  sprintf("\\quad Maximum intervals ($K$) & %d \\\\", K),
  sprintf("\\quad Observation window & %d days \\\\", WINDOW),
  sprintf("\\quad Static covariates ($P$) & %d \\\\", P),
  sprintf("\\quad Person-period rows ($\\mathcal{R}$) & %s \\\\",
          formatC(pp_rows, big.mark = ",")),
  sprintf("\\quad Event rows ($y_{i,t}=1$) & %s (%.3f\\%% of $\\mathcal{R}$) \\\\",
          formatC(sum(dat$delta), big.mark = ","), 100 * sum(dat$delta) / pp_rows),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Follow-up time, days}} \\\\[2pt]",
  sprintf("\\quad Overall: median [IQR] & %s \\\\",
          fmt_median_iqr(patients_df$time_to_event_days)),
  sprintf("\\quad Decedents: median [IQR] & %s \\\\",
          fmt_median_iqr(died_df$time_to_event_days)),
  sprintf("\\quad Discharged: median [IQR] & %s \\\\",
          fmt_median_iqr(discharged_df$time_to_event_days)),
  sprintf("\\quad Survived $\\geq$ %d days & %s patients \\\\",
          WINDOW, formatC(n_survived, big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}")


# ==============================================================================
# TABLE 3: COVARIATE DESCRIPTIONS (fits on page — two-column layout)
# ==============================================================================
message("==== Generating Table 3 ====")

# Use a compact tabularx layout: variable, description, prevalence
# No "Type" column (all binary except age_std, obvious from context)

cov_rows <- data.frame(
  var = colnames(X),
  stringsAsFactors = FALSE
) %>% mutate(
  label = case_when(
    var == "age_std"          ~ "Age (standardized to $\\mu=0$, $\\sigma=1$)",
    var == "male"             ~ "Male sex",
    var == "emergency"        ~ "Emergency admission",
    var == "icd_sepsis"       ~ "Sepsis (995.91/92, 038.x)",
    var == "icd_hf"           ~ "Heart failure (428.x)",
    var == "icd_resp"         ~ "Respiratory failure (518.8x)",
    var == "icd_aki"          ~ "Acute kidney injury (584.x)",
    var == "gcs_impaired"     ~ "GCS $\\leq 12$",
    var == "vasopressor"      ~ "Any vasopressor",
    var == "icd_pneumonia"    ~ "Pneumonia (480--486)",
    var == "icd_diabetes"     ~ "Diabetes (250.x)",
    var == "icd_liver"        ~ "Liver disease (571.x)",
    var == "icd_malignancy"   ~ "Malignancy (140--208)",
    var == "icd_copd"         ~ "COPD (490--492, 494, 496)",
    var == "is_ventilated"    ~ "Mechanical ventilation (proc 96.7x)",
    var == "received_dialysis"~ "Hemodialysis (proc 39.95)",
    var == "positive_blood_cx"~ "Positive blood culture",
    TRUE ~ var
  ),
  prev = sapply(seq_len(P), function(j) {
    if (colnames(X)[j] == "age_std") {
      sprintf("%.1f $\\pm$ %.1f yr",
              mean(patients_df$age, na.rm = TRUE),
              sd(patients_df$age, na.rm = TRUE))
    } else {
      sprintf("%.1f\\%%", 100 * mean(X[, j]))
    }
  })
)

tex3 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  sprintf("\\caption{Covariates in the design matrix $\\mathbf{X}$ ($P = %d$). The first nine covariates (above the rule) correspond to the simulation study; the remaining eight are additional clinical predictors. All variables except age are binary (0/1). ICD-9 codes shown in parentheses.}", P),
  "\\label{tab:covariates}",
  "\\begin{tabular}{r l l r}",
  "\\toprule",
  " & \\textbf{Variable} & \\textbf{Description} & \\textbf{Prevalence} \\\\",
  "\\midrule")

for (i in 1:nrow(cov_rows)) {
  if (i == 10) tex3 <- c(tex3, "\\midrule")
  tex3 <- c(tex3, sprintf("%d & \\texttt{%s} & %s & %s \\\\",
                          i,
                          gsub("_", "\\\\_", cov_rows$var[i]),
                          cov_rows$label[i],
                          cov_rows$prev[i]))
}

tex3 <- c(tex3,
          "\\bottomrule",
          "\\end{tabular}",
          "\\end{table}")


# ==============================================================================
# TABLE 4: CENSORING PATTERN BY INTERVAL
# ==============================================================================
message("==== Generating Table 4 ====")

interval_summary <- data.frame(
  interval  = 1:K,
  day_start = (0:(K-1)) * DELTA + 1,
  day_end   = (1:K) * DELTA
) %>% mutate(
  n_at_risk  = sapply(interval, function(t) sum(dat$y_obs >= t)),
  n_events   = sapply(interval, function(t) sum(dat$y_obs == t & dat$delta == 1)),
  n_censored = sapply(interval, function(t) sum(dat$y_obs == t & dat$delta == 0)),
  n_exits    = n_events + n_censored,
  cens_pct   = ifelse(n_exits > 0, 100 * n_censored / n_exits, 0),
  cum_events = cumsum(n_events),
  cum_cens   = cumsum(n_censored)
)

# Show all intervals (K=30 fits on a page)
tex4 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  sprintf("\\caption{Risk set composition across discrete-time intervals ($\\Delta = %d$ days). \\emph{Censored} here denotes patients discharged alive during each interval. The censoring proportion among exits quantifies the rate at which the informative censoring mechanism depletes the risk set---the empirical pattern interrogated by the sensitivity parameter $\\gamma$ in Section~\\ref{sc:sensitivity}.}", DELTA),
  "\\label{tab:censoring-pattern}",
  "\\begin{tabular}{r r r r r r}",
  "\\toprule",
  "Interval & Days & At Risk & Events & Discharged & \\% Disch. \\\\",
  "\\midrule")

# Show selected intervals to keep it compact
show_idx <- c(1:5, 10, 15, 20, 25, K)
show_idx <- unique(show_idx[show_idx <= K])

for (idx in show_idx) {
  r <- interval_summary[idx, ]
  # Add visual separator before the last block
  if (idx == show_idx[length(show_idx)] && idx > show_idx[length(show_idx) - 1] + 1) {
    tex4 <- c(tex4, "\\addlinespace")
  }
  tex4 <- c(tex4, sprintf("%d & %d--%d & %s & %s & %s & %.1f\\%% \\\\",
                          r$interval, r$day_start, r$day_end,
                          formatC(r$n_at_risk, big.mark = ","),
                          formatC(r$n_events, big.mark = ","),
                          formatC(r$n_censored, big.mark = ","),
                          r$cens_pct))
}

# Totals row
tex4 <- c(tex4,
          "\\midrule",
          sprintf("\\textbf{Total} & 1--%d & --- & %s & %s & --- \\\\",
                  WINDOW,
                  formatC(sum(interval_summary$n_events), big.mark = ","),
                  formatC(sum(interval_summary$n_censored), big.mark = ",")),
          "\\bottomrule",
          "\\end{tabular}",
          "\\end{table}")


# ==============================================================================
# SAVE ALL LATEX
# ==============================================================================
message("\n==== Saving LaTeX ====")

all_tex <- c(
  "% ==============================================================================",
  "% MIMIC-III Cohort Summary Tables (auto-generated by mimic_summary.R)",
  "% Include in appendix with: \\input{mimic_tables}",
  "% ==============================================================================",
  "",
  "% --- Table 1: Cohort Characteristics ---",
  tex1, "",
  "% --- Table 2: Data Structure ---",
  tex2, "",
  "% --- Table 3: Covariate Descriptions ---",
  tex3, "",
  "% --- Table 4: Censoring Pattern ---",
  tex4)

tex_file <- file.path(out.path, "mimic_tables.tex")
writeLines(all_tex, tex_file)
message(sprintf("  Saved: %s", tex_file))


# ==============================================================================
# FIGURE: FOLLOW-UP DISTRIBUTION (3-panel)
# ==============================================================================
message("==== Generating figure ====")

fig_file <- file.path(out.path, "mimic_followup.pdf")

pdf(fig_file, width = 10, height = 3.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), cex.lab = 0.9, cex.main = 1.0)

# --- Panel A: Follow-up time by outcome ---
# Stacked histograms
brks <- seq(0, max(patients_df$time_to_event_days, na.rm = TRUE) + 5, by = 3)
h_died <- hist(died_df$time_to_event_days, breaks = brks, plot = FALSE)
h_disc <- hist(discharged_df$time_to_event_days, breaks = brks, plot = FALSE)

ymax <- max(h_disc$counts, h_died$counts) * 1.1
plot(h_disc, col = adjustcolor("steelblue", 0.5), border = "steelblue",
     main = "(a) Follow-Up Time Distribution",
     xlab = "Days from ICU Admission", ylab = "Admissions",
     ylim = c(0, ymax), xlim = c(0, WINDOW * 1.5))
plot(h_died, col = adjustcolor("firebrick", 0.6), border = "firebrick", add = TRUE)
abline(v = WINDOW, col = "gray30", lty = 2, lwd = 1.5)
text(WINDOW, ymax * 0.9, sprintf("%d-day\nwindow", WINDOW),
     pos = 4, col = "gray30", cex = 0.75)
legend("bottomright",
       legend = c("Discharged alive", "Died"),
       fill = c(adjustcolor("steelblue", 0.5), adjustcolor("firebrick", 0.6)),
       border = c("steelblue", "firebrick"),
       cex = 0.7, bg = "white")

# --- Panel B: Risk set over intervals ---
plot(interval_summary$interval, interval_summary$n_at_risk,
     type = "s", lwd = 2, col = "black",
     main = "(b) Risk Set Over Intervals",
     xlab = sprintf("Interval (DELTA = %d days)", DELTA),
     ylab = "Patients at Risk",
     ylim = c(0, N * 1.05))

# Shade the discharged portion
polygon(c(1, interval_summary$interval, K),
        c(N, N - interval_summary$cum_cens, N - tail(interval_summary$cum_cens, 1)),
        col = adjustcolor("steelblue", 0.2), border = NA)
polygon(c(1, interval_summary$interval, K),
        c(N, N - interval_summary$cum_events, N - tail(interval_summary$cum_events, 1)),
        col = adjustcolor("firebrick", 0.15), border = NA)

lines(interval_summary$interval, interval_summary$n_at_risk,
      type = "s", lwd = 2, col = "black")

legend("topright",
       legend = c("At risk", "Cumul. discharged", "Cumul. deaths"),
       fill = c(NA, adjustcolor("steelblue", 0.3), adjustcolor("firebrick", 0.2)),
       border = c(NA, "steelblue", "firebrick"),
       lwd = c(2, NA, NA), col = c("black", NA, NA),
       cex = 0.65, bg = "white")

dev.off()
message(sprintf("  Saved: %s", fig_file))


# ==============================================================================
# CONSOLE SUMMARY
# ==============================================================================
message("\n============================================================")
message(sprintf("MIMIC-III COHORT SUMMARY (%d-day window)", WINDOW))
message("============================================================")
message(sprintf("  N = %s admissions", formatC(N, big.mark = ",")))
message(sprintf("  Died within %d days:           %s (%.1f%%)",
                WINDOW, formatC(n_died, big.mark = ","), 100 * n_died / N))
message(sprintf("  Discharged alive < %d days:    %s (%.1f%%)",
                WINDOW, formatC(n_discharged, big.mark = ","), 100 * n_discharged / N))
message(sprintf("  Survived >= %d days:           %s (%.1f%%)",
                WINDOW, formatC(n_survived, big.mark = ","), 100 * n_survived / N))
message("")
message(sprintf("  Follow-up (overall):  median %.1f days [%.1f, %.1f]",
                median(patients_df$time_to_event_days),
                quantile(patients_df$time_to_event_days, 0.25),
                quantile(patients_df$time_to_event_days, 0.75)))
message(sprintf("  Follow-up (died):     median %.1f days",
                median(died_df$time_to_event_days)))
message(sprintf("  Follow-up (discharged): median %.1f days",
                median(discharged_df$time_to_event_days)))
message(sprintf("  Person-period rows:   %s", formatC(sum(dat$y_obs), big.mark = ",")))
message(sprintf("  DELTA = %d days, K = %d, window = %d days", DELTA, K, WINDOW))
message("============================================================")
message(sprintf("\nOutput:"))
message(sprintf("  %s", tex_file))
message(sprintf("  %s", fig_file))
message("============================================================")


# ==============================================================================
# TABLE 5: MODEL DIAGNOSTICS & BETA COMPARISON (from saved analysis fits)
# ==============================================================================
# Loads mimic_analysis_fits.RData (or checkpoint) and generates:
#   - Tab 5: Beta comparison table (Model A mean, Model B mean, difference)
#            with alternating row shading — ready for main text
#   - Tab 6: Full diagnostics table (Rhat, ESS_bulk, ESS_tail) for both models
#            side-by-side — ready for appendix
#
# Also writes standalone .tex files for flexible \input{} inclusion.
# ==============================================================================

message("\n==== Generating Table 5 & 6: Model Diagnostics & Beta Comparison ====")

library(posterior)

# mimic_analysis.R saves to output/mimic/ under the project root
mimic_out <- file.path(out.path, "mimic")

fits_path <- file.path(mimic_out, "mimic_analysis_fits.RData")
if (!file.exists(fits_path)) {
  fits_path <- file.path(mimic_out, "mimic_analysis_checkpoint.RData")
}

message(sprintf("  Looking for fits at: %s", fits_path))

if (!file.exists(fits_path)) {
  message("  WARNING: No analysis RData found. Skipping Tables 5 & 6.")
  message("  Run mimic_analysis.R first, then re-run this script.")
} else {
  message(sprintf("  Loading: %s", fits_path))
  
  # Load into a temporary env so we don't clobber the summary workspace
  analysis_env <- new.env()
  load(fits_path, envir = analysis_env)
  
  # ---- Extract draws ----
  # mimic_analysis.R saves draws_A / draws_B (posterior::draws objects) when
  # cmdstanr was used, or fit_A / fit_B (stanfit objects) when rstan was used.
  
  if (exists("draws_A", envir = analysis_env)) {
    # cmdstanr path — draws are posterior::draws_array objects
    drA <- analysis_env$draws_A
    drB <- analysis_env$draws_B
  } else if (exists("fit_A", envir = analysis_env)) {
    # rstan path — convert to posterior draws
    drA <- as_draws(analysis_env$fit_A)
    drB <- as_draws(analysis_env$fit_B)
  } else {
    stop("Cannot find draws_A/draws_B or fit_A/fit_B in the loaded RData.")
  }
  
  beta_names <- colnames(dat$X)
  P_cov <- length(beta_names)
  
  # ---- Diagnostics for Model A (beta only) ----
  summ_A <- summarise_draws(
    subset_draws(drA, variable = paste0("beta[", 1:P_cov, "]")),
    mean, median, sd, mad,
    ~quantile(.x, probs = 0.05),
    ~quantile(.x, probs = 0.95),
    rhat, ess_bulk, ess_tail
  )
  summ_A$covariate <- beta_names
  
  # ---- Diagnostics for Model B (sigma_b + beta) ----
  summ_B_beta <- summarise_draws(
    subset_draws(drB, variable = paste0("beta[", 1:P_cov, "]")),
    mean, median, sd, mad,
    ~quantile(.x, probs = 0.05),
    ~quantile(.x, probs = 0.95),
    rhat, ess_bulk, ess_tail
  )
  summ_B_beta$covariate <- beta_names
  
  summ_B_sigma <- summarise_draws(
    subset_draws(drB, variable = "sigma_b"),
    mean, median, sd, mad,
    ~quantile(.x, probs = 0.05),
    ~quantile(.x, probs = 0.95),
    rhat, ess_bulk, ess_tail
  )
  summ_B_sigma$covariate <- "sigma_b"
  
  # ---- Console report (mirrors report_diagnostics output) ----
  message("\n============ Diagnostics: Model A (standard) ============")
  print(summ_A[, c("variable", "mean", "median", "sd", "mad",
                   "5%", "95%", "rhat", "ess_bulk", "ess_tail")])
  
  message("\n============ Diagnostics: Model B (frailty, PC prior) ============")
  print(rbind(
    summ_B_sigma[, c("variable", "mean", "median", "sd", "mad",
                     "5%", "95%", "rhat", "ess_bulk", "ess_tail")],
    summ_B_beta[, c("variable", "mean", "median", "sd", "mad",
                    "5%", "95%", "rhat", "ess_bulk", "ess_tail")]
  ))
  
  # ==== LaTeX: Beta comparison table (main text) ====
  # Programmatically generated from posterior draws — replaces hardcoded table
  mean_A  <- summ_A$mean
  mean_B  <- summ_B_beta$mean
  diff_AB <- mean_B - mean_A
  
  # Covariate display names (LaTeX-safe, matches paper order)
  cov_display <- c(
    "age (std)", "male", "emergency", "icd\\_sepsis", "icd\\_hf",
    "icd\\_resp", "icd\\_aki", "gcs\\_impaired", "vasopressor",
    "icd\\_pneumonia", "icd\\_diabetes", "icd\\_liver", "icd\\_malignancy",
    "icd\\_copd", "is\\_ventilated", "received\\_dialysis", "positive\\_blood\\_cx"
  )
  
  tex5 <- c(
    "\\begin{table}[t]",
    "\\centering",
    "\\small",
    "\\caption{Posterior mean covariate effects (log-odds ratio for daily hazard) under Model~A and Model~B. Difference column shows the attenuation correction from including frailty.}",
    "\\label{tab:beta-comparison}",
    "\\begin{tabular}{lrrr}",
    "\\toprule",
    "\\textbf{Covariate} & \\textbf{Model A} & \\textbf{Model B} & \\textbf{Difference} \\\\",
    "\\midrule"
  )
  
  for (j in 1:P_cov) {
    shade <- if (j %% 2 == 1) "\\rowcolor{rowgray} " else ""
    # Format values: negative numbers get $-$ for proper LaTeX minus sign
    fmt_val <- function(x) {
      if (x < 0) sprintf("$-$%.3f", abs(x)) else sprintf("%.3f", x)
    }
    d_str <- ifelse(diff_AB[j] < 0,
                    sprintf("$-$%.3f", abs(diff_AB[j])),
                    sprintf("+%.3f", diff_AB[j]))
    tex5 <- c(tex5, sprintf("%s%s & %s & %s & %s \\\\",
                            shade, cov_display[j],
                            fmt_val(mean_A[j]),
                            fmt_val(mean_B[j]),
                            d_str))
  }
  
  tex5 <- c(tex5,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}"
  )
  
  # ==== LaTeX: Full diagnostics table (appendix) ====
  # Model A and Model B side by side: Mean, SD, Rhat, ESS_bulk, ESS_tail
  
  tex6 <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\footnotesize",
    "\\caption{MCMC diagnostics for covariate effects ($\\hat{\\beta}$) under both models. $\\hat{R}$ values below 1.01 and ESS above 400 indicate reliable convergence. Model~B additionally estimates $\\sigma_b$ (frailty SD). All chains: 0 divergent transitions.}",
    "\\label{tab:diagnostics}",
    "\\begin{tabular}{l rrrrr rrrrr}",
    "\\toprule",
    " & \\multicolumn{5}{c}{\\textbf{Model A (Standard)}} & \\multicolumn{5}{c}{\\textbf{Model B (Frailty)}} \\\\",
    "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}",
    "\\textbf{Parameter} & Mean & SD & $\\hat{R}$ & ESS$_b$ & ESS$_t$ & Mean & SD & $\\hat{R}$ & ESS$_b$ & ESS$_t$ \\\\",
    "\\midrule"
  )
  
  # sigma_b row (Model B only)
  sb <- summ_B_sigma
  tex6 <- c(tex6, sprintf(
    "$\\sigma_b$ & --- & --- & --- & --- & --- & %.3f & %.3f & %.3f & %.0f & %.0f \\\\",
    sb$mean, sb$sd, sb$rhat, sb$ess_bulk, sb$ess_tail
  ))
  tex6 <- c(tex6, "\\midrule")
  
  for (j in 1:P_cov) {
    shade <- if (j %% 2 == 1) "\\rowcolor{rowgray} " else ""
    a <- summ_A[j, ]
    b <- summ_B_beta[j, ]
    tex6 <- c(tex6, sprintf(
      "%s%s & %.3f & %.4f & %.3f & %.0f & %.0f & %.3f & %.4f & %.3f & %.0f & %.0f \\\\",
      shade, cov_display[j],
      a$mean, a$sd, a$rhat, a$ess_bulk, a$ess_tail,
      b$mean, b$sd, b$rhat, b$ess_bulk, b$ess_tail
    ))
  }
  
  tex6 <- c(tex6,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}"
  )
  
  # ---- Append to the combined LaTeX output ----
  all_tex <- c(
    all_tex, "",
    "% --- Table 5: Beta Comparison (main text) ---",
    tex5, "",
    "% --- Table 6: Full MCMC Diagnostics (appendix) ---",
    tex6
  )
  
  # Re-write the combined mimic_tables.tex
  writeLines(all_tex, tex_file)
  message(sprintf("  Updated: %s (now includes Tables 5 & 6)", tex_file))
  
  # ---- Also write standalone files for flexible \input{} ----
  beta_cmp_file <- file.path(out.path, "mimic_beta_comparison.tex")
  writeLines(tex5, beta_cmp_file)
  message(sprintf("  Saved standalone: %s", beta_cmp_file))
  
  diag_file <- file.path(out.path, "mimic_diagnostics.tex")
  writeLines(tex6, diag_file)
  message(sprintf("  Saved standalone: %s", diag_file))
  
  # Cleanup
  rm(analysis_env, drA, drB, summ_A, summ_B_beta, summ_B_sigma)
  gc()
}
