# mimic_prep.R
# ==============================================================================
# MIMIC-III Data Preparation for Bayesian Discrete-Time Survival Analysis
# ==============================================================================
#
# PURPOSE:
#   Extract clinical variables from MIMIC-III and produce objects that plug
#   directly into core.r's expand_person_period() and the Stan analysis pipeline.
#
# OUTPUT (saved to out.path):
#   mimic_patients.rds  — patient-level data frame (one row per admission)
#   mimic_dat.rds       — list matching generate_synthetic_icu() interface:
#                          $N, $K, $P, $X, $y_obs, $delta
#                          so that expand_person_period(dat) works directly.
#
# USAGE:
#   1. Set DATA_SOURCE and paths below.
#   2. source("mimic_prep.R")          — runs once, saves .rds files
#   3. source("mimic_load.R")          — loads saved objects into workspace
#
# NOTE ON COLUMN NAMES:
#   The full MIMIC-III database (physionet.org/files/mimiciii/1.4) uses
#   UPPERCASE column names (HADM_ID, SUBJECT_ID, INTIME, etc.) while the
#   demo database uses lowercase. We normalize to lowercase immediately
#   after each fread() call so all downstream code works for either source.
#
#   Key tables and their columns used here:
#     ADMISSIONS:    hadm_id, subject_id, admittime, dischtime, deathtime,
#                    admission_type, hospital_expire_flag
#     PATIENTS:      subject_id, gender, dob
#     ICUSTAYS:      hadm_id, icustay_id, intime, outtime, first_careunit, los
#     CHARTEVENTS:   hadm_id, itemid, valuenum
#     LABEVENTS:     hadm_id, itemid, valuenum
#     DIAGNOSES_ICD: hadm_id, icd9_code
#     PROCEDURES_ICD: hadm_id, icd9_code
#     INPUTEVENTS_MV: hadm_id, itemid
#     MICROBIOLOGYEVENTS: hadm_id, spec_type_desc, org_itemid
#
# ==============================================================================

rm(list = ls())
library(data.table)
library(dplyr)
library(lubridate)
library(tidyr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DATA_SOURCE <- "full"     # "demo" or "full"
DELTA       <- 1          # discretization width in days (daily)
                          # simulation showed beta stable across Delta in {1,3,7,14}
                          # daily gives finest resolution where data is richest
K_MAX       <- 30         # max discrete intervals (30-day mortality window)
                          # standard ICU outcome metric; risk set still ~3K at day 30

base.path <- "/Users/adamkurth/Documents/RStudio/bayes-mimic-survival"

if (DATA_SOURCE == "full") {
  data.path <- file.path(base.path, "data/physionet.org/files/mimiciii/1.4")
} else {
  data.path <- file.path(base.path, "data/mimic-iii-clinical-database-demo-1.4")
}

out.path <- file.path(base.path, "data/rds")
if (!dir.exists(out.path)) dir.create(out.path, recursive = TRUE)
setwd(data.path)

# ==============================================================================
# HELPER: normalize column names to lowercase after fread
# ==============================================================================
# Full MIMIC uses UPPERCASE (HADM_ID), demo uses lowercase (hadm_id).
# This makes all downstream code work for either source.
read_mimic <- function(filename, ...) {
  dt <- fread(filename, ...)
  setnames(dt, tolower(names(dt)))
  dt
}

# ==============================================================================
# ITEM ID DICTIONARIES
# ==============================================================================
# These IDs map to D_ITEMS.csv and D_LABITEMS.csv.
# MIMIC-III has two charting systems (CareVue: older IDs, MetaVision: 220xxx+).

vital_ids <- list(
  heart_rate = c(211, 220045),      # bpm
  resp_rate  = c(618, 220210),      # breaths/min
  spo2       = c(646, 220277),      # %
  temp_c     = c(676, 223762),      # Celsius
  temp_f     = c(678, 223761),      # Fahrenheit (will convert)
  gcs        = c(198, 220739)       # 3-15 scale
)
all_vital_ids <- unlist(vital_ids)

lab_ids <- list(
  creatinine = 50912,     # mg/dL — renal function
  platelets  = 51265,     # K/uL — coagulation / SOFA
  wbc        = 51301,     # K/uL — infection marker
  lactate    = 50813      # mmol/L — tissue perfusion / sepsis
)
all_lab_ids <- unlist(lab_ids)

vasopressor_ids_mv <- c(
  221289,   # Epinephrine
  221662,   # Dopamine
  221749,   # Phenylephrine
  221906,   # Norepinephrine
  222315    # Vasopressin
)

# ==============================================================================
# SECTION 1: CORE TABLES — Demographics, Admissions, ICU Stays
# ==============================================================================
message("==== SECTION 1: Loading core tables ====")

admissions <- read_mimic("ADMISSIONS.csv")
patients   <- read_mimic("PATIENTS.csv")
icustays   <- read_mimic("ICUSTAYS.csv")

message(sprintf("  ADMISSIONS: %d rows, columns: %s",
                nrow(admissions), paste(names(admissions), collapse = ", ")))
message(sprintf("  PATIENTS:   %d rows", nrow(patients)))
message(sprintf("  ICUSTAYS:   %d rows", nrow(icustays)))

# First ICU stay per hospital admission (defines Time Zero)
first_icu <- icustays %>%
  group_by(hadm_id) %>%
  arrange(intime) %>%
  slice(1) %>%
  ungroup() %>%
  select(hadm_id, icustay_id,
         icu_intime = intime, icu_outtime = outtime,
         first_careunit, los_icu = los) %>%
  as.data.table()

# Build cohort
cohort <- merge(admissions, patients, by = "subject_id", all.x = TRUE)
cohort <- merge(cohort, first_icu, by = "hadm_id", all.x = TRUE)

cohort <- cohort %>%
  mutate(
    admittime   = as.POSIXct(admittime),
    dischtime   = as.POSIXct(dischtime),
    deathtime   = as.POSIXct(deathtime),
    dob         = as.POSIXct(dob),
    icu_intime  = as.POSIXct(icu_intime),
    icu_outtime = as.POSIXct(icu_outtime),
    # Age: MIMIC shifts DOB for patients >89 (privacy), cap at 90
    age = as.numeric(difftime(admittime, dob, units = "days")) / 365.25,
    age = ifelse(age > 150, 90, age),
    hospital_expire_flag = as.integer(hospital_expire_flag)
  ) %>%
  as.data.table()

message(sprintf("  Cohort: %d admissions, %d unique patients",
                nrow(cohort), length(unique(cohort$subject_id))))

# ==============================================================================
# SECTION 2: VITAL SIGN SUMMARIES (per-admission means)
# ==============================================================================
message("==== SECTION 2: Loading vital signs from CHARTEVENTS ====")

# For full MIMIC, CHARTEVENTS is ~330M rows. Filter on itemid during read.
chartevents <- read_mimic("CHARTEVENTS.csv")[itemid %in% all_vital_ids]

message(sprintf("  CHARTEVENTS (filtered): %d rows", nrow(chartevents)))

vitals_long <- chartevents %>%
  filter(!is.na(valuenum)) %>%
  mutate(
    vital_type = case_when(
      itemid %in% vital_ids$heart_rate ~ "hr",
      itemid %in% vital_ids$resp_rate  ~ "resp_rate",
      itemid %in% vital_ids$spo2       ~ "spo2",
      itemid %in% vital_ids$temp_c     ~ "temp_c",
      itemid %in% vital_ids$temp_f     ~ "temp_f",
      itemid %in% vital_ids$gcs        ~ "gcs",
      TRUE ~ "other"
    ),
    # Convert Fahrenheit to Celsius
    valuenum = ifelse(vital_type == "temp_f", (valuenum - 32) * 5 / 9, valuenum),
    vital_type = ifelse(vital_type == "temp_f", "temp_c", vital_type)
  ) %>%
  select(hadm_id, vital_type, valuenum)

vital_means <- vitals_long %>%
  group_by(hadm_id, vital_type) %>%
  summarize(mean_val = mean(valuenum, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = vital_type, values_from = mean_val,
              names_prefix = "vital_") %>%
  as.data.table()

message(sprintf("  Vital means: %d admissions", nrow(vital_means)))

# Free memory — chartevents is huge
rm(chartevents, vitals_long); gc()

# ==============================================================================
# SECTION 3: LAB SUMMARIES (per-admission means)
# ==============================================================================
message("==== SECTION 3: Loading labs from LABEVENTS ====")

labevents <- read_mimic("LABEVENTS.csv")[itemid %in% all_lab_ids]

message(sprintf("  LABEVENTS (filtered): %d rows", nrow(labevents)))

lab_means <- labevents %>%
  filter(!is.na(valuenum)) %>%
  mutate(
    lab_type = case_when(
      itemid == lab_ids$creatinine ~ "creatinine",
      itemid == lab_ids$platelets  ~ "platelets",
      itemid == lab_ids$wbc        ~ "wbc",
      itemid == lab_ids$lactate    ~ "lactate",
      TRUE ~ "other"
    )
  ) %>%
  group_by(hadm_id, lab_type) %>%
  summarize(mean_val = mean(valuenum, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = lab_type, values_from = mean_val,
              names_prefix = "lab_") %>%
  as.data.table()

message(sprintf("  Lab means: %d admissions", nrow(lab_means)))
rm(labevents); gc()

# ==============================================================================
# SECTION 4: DIAGNOSES (ICD-9 binary flags)
# ==============================================================================
message("==== SECTION 4: Loading diagnoses ====")

dx_icd <- read_mimic("DIAGNOSES_ICD.csv")

dx_flags <- dx_icd %>%
  group_by(hadm_id) %>%
  summarize(
    icd_sepsis     = as.integer(any(grepl("^99591|^99592|^038", icd9_code))),
    icd_hf         = as.integer(any(grepl("^428", icd9_code))),
    icd_aki        = as.integer(any(grepl("^584", icd9_code))),
    icd_resp       = as.integer(any(grepl("^51881|^51882|^51884", icd9_code))),
    icd_pneumonia  = as.integer(any(grepl("^48[0-6]", icd9_code))),
    icd_diabetes   = as.integer(any(grepl("^250", icd9_code))),
    icd_liver      = as.integer(any(grepl("^571", icd9_code))),
    icd_malignancy = as.integer(any(grepl("^1[4-9][0-9]|^20[0-8]", icd9_code))),
    icd_copd       = as.integer(any(grepl("^496|^49[0-2]|^494", icd9_code))),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Diagnoses: %d admissions with ICD codes", nrow(dx_flags)))
rm(dx_icd); gc()

# ==============================================================================
# SECTION 5: PROCEDURES (ventilation, dialysis)
# ==============================================================================
message("==== SECTION 5: Loading procedures ====")

proc_icd <- read_mimic("PROCEDURES_ICD.csv")

proc_flags <- proc_icd %>%
  group_by(hadm_id) %>%
  summarize(
    is_ventilated     = as.integer(any(grepl("^967", icd9_code))),
    received_dialysis = as.integer(any(icd9_code == "3995")),
    .groups = "drop"
  ) %>%
  as.data.table()

rm(proc_icd); gc()

# ==============================================================================
# SECTION 6: VASOPRESSORS (INPUTEVENTS_MV)
# ==============================================================================
message("==== SECTION 6: Loading vasopressors ====")

inputevents_mv <- read_mimic("INPUTEVENTS_MV.csv")

vaso_flags <- inputevents_mv %>%
  filter(itemid %in% vasopressor_ids_mv) %>%
  group_by(hadm_id) %>%
  summarize(on_vasopressors = 1L, .groups = "drop") %>%
  as.data.table()

message(sprintf("  Vasopressors: %d admissions", nrow(vaso_flags)))
rm(inputevents_mv); gc()

# ==============================================================================
# SECTION 7: MICROBIOLOGY (positive blood culture)
# ==============================================================================
message("==== SECTION 7: Loading microbiology ====")

micro <- read_mimic("MICROBIOLOGYEVENTS.csv")

blood_cx <- micro %>%
  filter(spec_type_desc == "BLOOD CULTURE") %>%
  group_by(hadm_id) %>%
  summarize(
    positive_blood_cx = as.integer(any(!is.na(org_itemid))),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Blood cultures: %d admissions", nrow(blood_cx)))
rm(micro); gc()

# ==============================================================================
# SECTION 8: MERGE INTO ANALYTIC DATASET
# ==============================================================================
message("\n==== SECTION 8: Merging into analytic dataset ====")

df <- cohort %>%
  select(hadm_id, subject_id, admittime, dischtime, deathtime,
         icu_intime, icu_outtime, first_careunit,
         age, gender, admission_type, hospital_expire_flag)

df <- merge(df, vital_means, by = "hadm_id", all.x = TRUE)
df <- merge(df, lab_means,   by = "hadm_id", all.x = TRUE)
df <- merge(df, dx_flags,    by = "hadm_id", all.x = TRUE)
df <- merge(df, proc_flags,  by = "hadm_id", all.x = TRUE)
df <- merge(df, vaso_flags,  by = "hadm_id", all.x = TRUE)
df <- merge(df, blood_cx,    by = "hadm_id", all.x = TRUE)

message(sprintf("  Merged dataset: %d rows x %d columns", nrow(df), ncol(df)))

# Free intermediate tables
rm(cohort, vital_means, lab_means, dx_flags, proc_flags, vaso_flags, blood_cx); gc()

# ==============================================================================
# SECTION 9: CONSTRUCT SURVIVAL VARIABLES & CLEAN COVARIATES
# ==============================================================================
message("==== SECTION 9: Constructing survival variables ====")

patients_df <- df %>%
  as_tibble() %>%
  mutate(
    # --- Survival outcomes ---
    end_time = coalesce(deathtime, dischtime),
    time_to_event_days = as.numeric(difftime(end_time, icu_intime, units = "days")),
    event_status = as.integer(!is.na(deathtime)),

    # --- Demographics (numeric) ---
    male      = as.integer(gender == "M"),
    emergency = as.integer(admission_type == "EMERGENCY"),

    # --- GCS impairment (binary: mean GCS <= 12) ---
    gcs_impaired = as.integer(!is.na(vital_gcs) & vital_gcs <= 12),

    # --- Fill NA binary flags with 0 (absence = not coded) ---
    across(starts_with("icd_"), ~ ifelse(is.na(.), 0L, .)),
    is_ventilated     = ifelse(is.na(is_ventilated), 0L, is_ventilated),
    received_dialysis = ifelse(is.na(received_dialysis), 0L, received_dialysis),
    on_vasopressors   = ifelse(is.na(on_vasopressors), 0L, on_vasopressors),
    positive_blood_cx = ifelse(is.na(positive_blood_cx), 0L, positive_blood_cx)
  ) %>%
  # --- Filter to valid survival data ---
  filter(
    !is.na(time_to_event_days),
    time_to_event_days > 0,
    !is.na(icu_intime)
  ) %>%
  # --- Discretize time ---
  mutate(
    y_obs_raw = ceiling(time_to_event_days / DELTA),
    y_obs     = pmin(y_obs_raw, K_MAX),
    # If patient survived past K_MAX, administratively censored
    delta = ifelse(y_obs_raw > K_MAX & event_status == 1, 0L, event_status)
  )

message(sprintf("  Patients after survival filter: %d", nrow(patients_df)))
message(sprintf("  Event rate: %.1f%%", 100 * mean(patients_df$event_status)))

rm(df); gc()

# ==============================================================================
# SECTION 10: BUILD NUMERIC X MATRIX
# ==============================================================================
message("==== SECTION 10: Building X matrix ====")

# Standardize age
patients_df <- patients_df %>%
  mutate(age_std = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE))

# Core covariates (match simulation in core.r):
x_cols_core <- c(
  "age_std", "male", "emergency",
  "icd_sepsis", "icd_hf", "icd_resp", "icd_aki",
  "gcs_impaired", "on_vasopressors"
)

# Extended clinical predictors:
x_cols_extended <- c(
  "icd_pneumonia", "icd_diabetes", "icd_liver",
  "icd_malignancy", "icd_copd",
  "is_ventilated", "received_dialysis", "positive_blood_cx"
)

x_cols <- c(x_cols_core, x_cols_extended)

# Check availability
available_cols <- x_cols[x_cols %in% names(patients_df)]
missing_cols   <- setdiff(x_cols, available_cols)
if (length(missing_cols) > 0) {
  message("  WARNING: columns not found: ", paste(missing_cols, collapse = ", "))
}

# Build X, drop incomplete cases
X <- patients_df %>% select(all_of(available_cols)) %>% as.matrix()

complete <- complete.cases(X)
n_dropped <- sum(!complete)
message(sprintf("  Dropping %d rows with NA covariates (%d remain)",
                n_dropped, sum(complete)))

patients_df <- patients_df[complete, ]
X <- X[complete, , drop = FALSE]

# Rename to match simulation convention
colnames(X)[colnames(X) == "on_vasopressors"] <- "vasopressor"

# ==============================================================================
# SECTION 11: COMPUTE DIMENSIONS
# ==============================================================================
K <- as.integer(max(patients_df$y_obs))
N <- nrow(patients_df)
P <- ncol(X)

message(sprintf("\n  ---- Final cohort ----"))
message(sprintf("  N = %d patients", N))
message(sprintf("  P = %d covariates", P))
message(sprintf("  K = %d intervals (DELTA = %d days, max window = %d days)",
                K, DELTA, K * DELTA))
message(sprintf("  Events:   %d (%.1f%%)", sum(patients_df$delta),
                100 * mean(patients_df$delta)))
message(sprintf("  Censored: %d (%.1f%%)", sum(1 - patients_df$delta),
                100 * mean(1 - patients_df$delta)))

# Quick covariate prevalence summary
message("\n  Covariate prevalences:")
for (j in 1:P) {
  cname <- colnames(X)[j]
  if (cname == "age_std") {
    message(sprintf("    %-20s  mean=%.2f  sd=%.2f", cname, mean(X[,j]), sd(X[,j])))
  } else {
    message(sprintf("    %-20s  %.1f%%", cname, 100 * mean(X[,j])))
  }
}

# ==============================================================================
# SECTION 12: PACKAGE INTO dat LIST
# ==============================================================================
# Matches generate_synthetic_icu() output so expand_person_period(dat) works.

dat <- list(
  N     = N,
  K     = K,
  P     = P,
  delta_width = DELTA,
  X     = X,
  y_obs = as.integer(patients_df$y_obs),
  delta = as.integer(patients_df$delta)
)

# ==============================================================================
# SECTION 13: PERSON-PERIOD EXPANSION
# ==============================================================================
# Run expand_person_period() here so mimic_load.R can load the result
# directly, avoiding the ~350K-row expansion at analysis time.
message("\n==== SECTION 13: Person-period expansion ====")

source(file.path(base.path, "demo/core.r"))
pp_data <- expand_person_period(dat)

message(sprintf("  R = %d person-period rows", pp_data$R))
message(sprintf("  X.ast dimensions: %d x %d", nrow(pp_data$X.ast), ncol(pp_data$X.ast)))

# ==============================================================================
# SECTION 14: SAVE
# ==============================================================================
message("\n==== SECTION 14: Saving ====")

saveRDS(patients_df, file = file.path(out.path, "mimic_patients.rds"))
saveRDS(dat,         file = file.path(out.path, "mimic_dat.rds"))
saveRDS(pp_data,     file = file.path(out.path, "mimic_pp_data.rds"))

message(sprintf("  Saved to: %s", out.path))
message("    mimic_patients.rds  — patient-level data frame")
message("    mimic_dat.rds       — list for expand_person_period(dat)")
message("    mimic_pp_data.rds   — person-period expanded data (ready for Stan)")

# ==============================================================================
# SECTION 15: SUMMARY
# ==============================================================================
message("\n============ MIMIC PREP COMPLETE ============")
message(sprintf("  dat$N = %d, dat$K = %d, dat$P = %d", dat$N, dat$K, dat$P))
message(sprintf("  dat$X columns: %s", paste(colnames(dat$X), collapse = ", ")))
message(sprintf("  Event rate: %.1f%%", 100 * mean(dat$delta)))
message(sprintf("  Median observed time: %.0f days", median(dat$y_obs) * DELTA))
message(sprintf("  Max observed time: %.0f days", max(dat$y_obs) * DELTA))
message(sprintf("  Person-period rows: %d", pp_data$R))
message("")
message("  Next steps:")
message("    source('mimic_load.R')   # loads dat, patients_df, pp_data")
message("    source('mimic_analysis.R')")
message("=========================================")

setwd(base.path)
