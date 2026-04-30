# data-prep.R
# ==============================================================================
# MIMIC-III Data Extraction for Bayesian Hierarchical Survival Modeling
# ==============================================================================
#
# PURPOSE:
# --------
# Extract clinical variables from MIMIC-III for Bayesian hierarchical survival
# analysis. Creates THREE output datasets with distinct structures:
#
#   1. ready_df (Static Baseline Matrix X)
#      - Dimension: N x P where N = number of admissions, P = number of covariates
#      - One row per hospital admission (unit of analysis)
#      - Contains: demographics, aggregated vitals/labs, comorbidities, procedures
#      - Represents everything known at "Time Zero" (ICU admission)
#      - Use for: standard survival analysis, binary outcome models, latent class discovery
#
#   2. longitudinal_df (Time-Varying Covariate Matrix Z)
#      - Dimension: variable (one row per measurement occasion)
#      - Structure: Z_{i,t} = vitals for patient i at time t
#      - Contains: raw vital sign measurements with timestamps
#      - Use for: joint longitudinal-survival models, trajectory analysis
#
#   3. person_period_df (Discrete-Time Survival Format)
#      - Dimension: sum of T_i periods across all patients
#      - Structure: one row per patient i per discrete period t
#      - Contains: y_{i,t} binary outcome (1 = event in period t), time-varying covariates
#      - Use for: Polya-Gamma augmented Gibbs sampling, discrete-time hazard models
#
# CONCEPTUAL FRAMEWORK:
# ---------------------
# For patient i with T_i discrete time periods:
#
#   Static covariates (X_i):
#     - Age, gender, admission type, comorbidities, etc.
#     - These are REPEATED across all periods for the same patient
#     - Captured in ready_df and carried to person_period_df
#
#   Time-varying covariates (Z_{i,t}):
#     - Heart rate, MAP, temperature, etc. at each time point
#     - Different values at each t
#     - Captured in longitudinal_df, aggregated to periods in person_period_df
#
#   Discrete-time outcome:
#     - y_{i,t} = 1 if patient dies in period t, 0 otherwise
#     - Patient who dies:     y_{i,1}=0, y_{i,2}=0, y_{i,T_i}=1
#     - Patient censored:     y_{i,1}=0, y_{i,2}=0, y_{i,T_i}=0
#
# ==============================================================================

rm(list = ls())

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
library(data.table)   # Fast CSV reading with fread()
library(dplyr)        # Data manipulation (pipes, group_by, summarize)
library(lubridate)    # Date/time parsing
library(tidyr)        # Reshaping (pivot_wider, uncount, replace_na)
library(zoo)          # na.locf() for last observation carried forward
library(forcats)      # Factor manipulation (fct_na_value_to_level)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DATA_SOURCE <- "demo"  # Options: "demo" or "full"

# Paths
base.path <- "/Users/adamkurth/Documents/RStudio/bayes-mimic-survival"
data.path <- file.path(base.path, "data/mimic-iii-clinical-database-demo-1.4/")
out.path <- file.path(base.path, "data/rds")

setwd(data.path)

# ==============================================================================
# ITEM ID DICTIONARIES
# ==============================================================================
# These IDs map to D_ITEMS.csv and D_LABITEMS.csv in the MIMIC-III database.
# MIMIC-III has two charting systems (CareVue and MetaVision) with different ID schemes.

# --- VITAL SIGNS (from CHARTEVENTS) ---
# CareVue (CV) was the older system; MetaVision (MV) is newer.
# We must capture both to avoid losing data.
vital_ids <- list(
  heart_rate = c(211, 220045),      # CV: 211, MV: 220045 (bpm)
  sbp        = c(51, 220050),       # Systolic BP: CV: 51, MV: 220050 (mmHg)
  dbp        = c(8364, 220051),     # Diastolic BP: CV: 8364, MV: 220051 (mmHg)
  map        = c(52, 220052),       # Mean Arterial Pressure: CV: 52, MV: 220052 (mmHg)
  resp_rate  = c(618, 220210),      # Respiratory Rate: CV: 618, MV: 220210 (breaths/min)
  spo2       = c(646, 220277),      # SpO2: CV: 646, MV: 220277 (%)
  temp_c     = c(676, 223762),      # Temperature Celsius: CV: 676, MV: 223762
  temp_f     = c(678, 223761),      # Temperature Fahrenheit: CV: 678, MV: 223761 (will convert to C)
  gcs        = c(198, 220739)       # GCS Total: CV: 198, MV: 220739 (3-15 scale)
)
all_vital_ids <- unlist(vital_ids)

# --- LABORATORY VALUES (from LABEVENTS) ---
# Selected based on SOFA score components and clinical relevance for ICU mortality
lab_ids <- list(
  lactate    = 50813,    # Blood lactate - tissue perfusion, sepsis marker (mmol/L)
  creatinine = 50912,    # Serum creatinine - renal function (mg/dL)
  bilirubin  = 50885,    # Total bilirubin - hepatic function, SOFA component (mg/dL)
  platelets  = 51265,    # Platelet count - coagulation, SOFA component (K/uL)
  wbc        = 51301,    # White blood cells - infection/immune response (K/uL)
  hemoglobin = 51222,    # Hemoglobin - anemia, bleeding (g/dL)
  glucose    = 50931,    # Serum glucose - metabolic status (mg/dL)
  bun        = 51006,    # Blood urea nitrogen - renal function (mg/dL)
  sodium     = 50983,    # Serum sodium - electrolyte balance (mEq/L)
  potassium  = 50971,    # Serum potassium - electrolyte, cardiac risk (mEq/L)
  pao2       = 50821,    # Arterial pO2 - oxygenation, SOFA respiratory (mmHg)
  pco2       = 50818,    # Arterial pCO2 - ventilation status (mmHg)
  ph         = 50820     # Arterial pH - acid-base status (unitless, 7.35-7.45 normal)
)
all_lab_ids <- unlist(lab_ids)

# --- VASOPRESSORS (from INPUTEVENTS_MV) ---
# These are strong indicators of hemodynamic instability/shock
vasopressor_ids_mv <- c(
  221289,  # Epinephrine - cardiac arrest, anaphylaxis
  221662,  # Dopamine - inotropic support
  221749,  # Phenylephrine - pure alpha-agonist vasoconstrictor
  221906,  # Norepinephrine - first-line for septic shock
  222315   # Vasopressin - second-line adjunct
)

# ==============================================================================
# SECTION 1: CORE TABLES - Demographics, Admissions, ICU Stays
# ==============================================================================
# These tables form the COHORT - the patient-admission level structure
# that serves as the backbone for dat (the static X matrix)

message("=== SECTION 1: Loading core tables ===")

admissions <- fread("ADMISSIONS.csv")
patients   <- fread("PATIENTS.csv")
icustays   <- fread("ICUSTAYS.csv")

# ---------------------------------------------------------------------------
# DEFINE TIME ZERO: First ICU admission within each hospital stay
# ---------------------------------------------------------------------------
# Justification: The ICU admission marks entry to highest-acuity care.
# Survival from this point is the clinically relevant outcome.
# We use the FIRST ICU stay to avoid survivorship bias.

first_icu <- icustays %>%
  group_by(hadm_id) %>%
  arrange(intime) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    hadm_id,
    icustay_id,
    icu_intime = intime,      # TIME ZERO for survival analysis
    icu_outtime = outtime,
    first_careunit,           # ICU type: MICU, SICU, CCU, CSRU, TSICU
    los_icu = los             # Length of stay in ICU (days)
  ) %>%
  as.data.table()

# Count total ICU stays per admission - measure of clinical complexity
icu_count <- icustays %>%
  group_by(hadm_id) %>%
  summarize(
    n_icu_stays = n(),
    total_icu_los = sum(los, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.table()

# ---------------------------------------------------------------------------
# BUILD COHORT: Merge patient demographics with ICU stay info
# ---------------------------------------------------------------------------
# This creates the patient-admission level dataset that becomes "ready_df"

cohort <- merge(admissions, patients, by = "subject_id", all.x = TRUE)
cohort <- merge(cohort, first_icu, by = "hadm_id", all.x = TRUE)
cohort <- merge(cohort, icu_count, by = "hadm_id", all.x = TRUE)

# Process timestamps and demographics
cohort <- cohort %>%
  mutate(
    # Parse timestamps
    admittime   = as.POSIXct(admittime),
    dischtime   = as.POSIXct(dischtime),
    deathtime   = as.POSIXct(deathtime),
    dob         = as.POSIXct(dob),
    icu_intime  = as.POSIXct(icu_intime),
    icu_outtime = as.POSIXct(icu_outtime),

    # Age calculation
    # MIMIC-III shifts DOB for patients > 89 years (privacy protection)
    # This results in apparent ages > 200; we cap at 90 per MIMIC convention
    age = as.numeric(difftime(admittime, dob, units = "days")) / 365.25,
    age = ifelse(age > 150, 90, age),

    # Convert demographics to factors for hierarchical modeling
    gender = as.factor(gender),
    admission_type = as.factor(admission_type),
    insurance = as.factor(insurance),
    ethnicity = as.factor(ethnicity),
    marital_status = as.factor(marital_status),

    # Hospital death indicator (from ADMISSIONS table)
    hospital_expire_flag = as.integer(hospital_expire_flag)
  ) %>%
  # Remove columns we don't need to avoid clutter (use any_of to avoid errors if missing)
  select(-any_of(c("row_id", "row_id.x", "row_id.y"))) %>%
  as.data.table()

message(sprintf("  Cohort: %d admissions, %d unique patients",
                nrow(cohort), length(unique(cohort$subject_id))))

# ==============================================================================
# SECTION 2: VITAL SIGNS (CHARTEVENTS)
# ==============================================================================
# Creates TWO outputs:
#   1. vitals_wide: One row per (admission, charttime) - forms basis of longitudinal_df
#   2. vital_summaries: One row per admission with aggregated stats - merges into ready_df

message("=== SECTION 2: Loading vital signs from CHARTEVENTS ===")

chartevents <- fread("CHARTEVENTS.csv")[itemid %in% all_vital_ids]

# ---------------------------------------------------------------------------
# Create long-format vitals with harmonized variable names
# ---------------------------------------------------------------------------
vitals_long <- chartevents %>%
  filter(!is.na(valuenum)) %>%
  mutate(
    # Map CareVue and MetaVision IDs to common variable names
    vital_type = case_when(
      itemid %in% vital_ids$heart_rate ~ "hr",
      itemid %in% vital_ids$sbp        ~ "sbp",
      itemid %in% vital_ids$dbp        ~ "dbp",
      itemid %in% vital_ids$map        ~ "map",
      itemid %in% vital_ids$resp_rate  ~ "resp_rate",
      itemid %in% vital_ids$spo2       ~ "spo2",
      itemid %in% vital_ids$temp_c     ~ "temp_c",
      itemid %in% vital_ids$temp_f     ~ "temp_f",  # Will convert below
      itemid %in% vital_ids$gcs        ~ "gcs",
      TRUE ~ "other"
    ),
    charttime = as.POSIXct(charttime)
  ) %>%
  # Convert Fahrenheit to Celsius for uniformity
  mutate(
    valuenum = ifelse(vital_type == "temp_f", (valuenum - 32) * 5 / 9, valuenum),
    vital_type = ifelse(vital_type == "temp_f", "temp_c", vital_type)
  ) %>%
  select(subject_id, hadm_id, icustay_id, charttime, vital_type, valuenum)

# ---------------------------------------------------------------------------
# Create wide-format vitals: one row per (admission, time point)
# This becomes the basis for longitudinal_df (Z_{i,t} matrix)
# ---------------------------------------------------------------------------
vitals_wide <- vitals_long %>%
  group_by(subject_id, hadm_id, icustay_id, charttime, vital_type) %>%
  summarize(value = mean(valuenum, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = vital_type, values_from = value) %>%
  as.data.table()

message(sprintf("  Vitals: %d measurement occasions across %d admissions",
                nrow(vitals_wide), length(unique(vitals_wide$hadm_id))))

# ---------------------------------------------------------------------------
# Per-admission vital summaries for ready_df (static X matrix)
# ---------------------------------------------------------------------------
# These aggregate statistics summarize the entire ICU stay into single values
# Justification: Different summaries capture different aspects:
#   - mean: overall level
#   - min/max: extremes (decompensation events)
#   - sd: variability (instability)
#   - n_obs: data completeness

vital_summaries <- vitals_long %>%
  group_by(hadm_id, vital_type) %>%
  summarize(
    mean_val = mean(valuenum, na.rm = TRUE),
    min_val = min(valuenum, na.rm = TRUE),
    max_val = max(valuenum, na.rm = TRUE),
    sd_val = sd(valuenum, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = vital_type,
    values_from = c(mean_val, min_val, max_val, sd_val, n_obs),
    names_glue = "{vital_type}_{.value}"
  ) %>%
  as.data.table()

# ==============================================================================
# SECTION 3: LABORATORY VALUES (LABEVENTS)
# ==============================================================================
# Creates lab_summaries: One row per admission with aggregated stats

message("=== SECTION 3: Loading laboratory values ===")

labevents <- fread("LABEVENTS.csv")[itemid %in% all_lab_ids]

labs_long <- labevents %>%
  filter(!is.na(valuenum)) %>%
  mutate(
    lab_type = case_when(
      itemid == lab_ids$lactate    ~ "lactate",
      itemid == lab_ids$creatinine ~ "creatinine",
      itemid == lab_ids$bilirubin  ~ "bilirubin",
      itemid == lab_ids$platelets  ~ "platelets",
      itemid == lab_ids$wbc        ~ "wbc",
      itemid == lab_ids$hemoglobin ~ "hemoglobin",
      itemid == lab_ids$glucose    ~ "glucose",
      itemid == lab_ids$bun        ~ "bun",
      itemid == lab_ids$sodium     ~ "sodium",
      itemid == lab_ids$potassium  ~ "potassium",
      itemid == lab_ids$pao2       ~ "pao2",
      itemid == lab_ids$pco2       ~ "pco2",
      itemid == lab_ids$ph         ~ "ph",
      TRUE ~ "other"
    ),
    charttime = as.POSIXct(charttime)
  ) %>%
  select(subject_id, hadm_id, charttime, lab_type, valuenum)

# Per-admission lab summaries
# Justification: first/last values capture trajectory direction
lab_summaries <- labs_long %>%
  group_by(hadm_id, lab_type) %>%
  summarize(
    first_val = first(valuenum),  # Baseline at admission
    last_val = last(valuenum),    # Final status
    min_val = min(valuenum, na.rm = TRUE),
    max_val = max(valuenum, na.rm = TRUE),
    mean_val = mean(valuenum, na.rm = TRUE),
    n_labs = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = lab_type,
    values_from = c(first_val, last_val, min_val, max_val, mean_val, n_labs),
    names_glue = "{lab_type}_{.value}"
  ) %>%
  as.data.table()

message(sprintf("  Labs: %d measurements across %d admissions",
                nrow(labs_long), length(unique(labs_long$hadm_id))))

# ==============================================================================
# SECTION 4: MICROBIOLOGY - Infection Markers
# ==============================================================================
message("=== SECTION 4: Loading microbiology data ===")

micro <- fread("MICROBIOLOGYEVENTS.csv")

# Positive blood cultures indicate bacteremia - strong mortality predictor
blood_culture <- micro %>%
  filter(spec_type_desc == "BLOOD CULTURE") %>%
  group_by(hadm_id) %>%
  summarize(
    n_blood_cultures = n(),
    positive_blood_culture = as.integer(any(!is.na(org_itemid))),
    n_positive_cultures = sum(!is.na(org_itemid)),
    n_organisms = n_distinct(org_itemid, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.table()

# Any positive culture (any specimen type)
any_positive_culture <- micro %>%
  filter(!is.na(org_itemid)) %>%
  group_by(hadm_id) %>%
  summarize(
    any_positive_culture = 1L,
    n_positive_specimens = n(),
    n_unique_organisms = n_distinct(org_itemid),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Microbiology: %d admissions with blood cultures",
                nrow(blood_culture)))

# ==============================================================================
# SECTION 5: DIAGNOSES (ICD-9) AND DRG CODES
# ==============================================================================
message("=== SECTION 5: Loading diagnosis codes ===")

dx_icd <- fread("DIAGNOSES_ICD.csv")

# Create comorbidity flags based on ICD-9 patterns
# These codes are assigned retrospectively for billing but provide
# standardized disease categorization
dx_summaries <- dx_icd %>%
  group_by(hadm_id) %>%
  summarize(
    n_diagnoses = n(),
    icd_sepsis = as.integer(any(grepl("^99591|^99592|^038", icd9_code))),
    icd_heart_failure = as.integer(any(grepl("^428", icd9_code))),
    icd_aki = as.integer(any(grepl("^584", icd9_code))),
    icd_ckd = as.integer(any(grepl("^585", icd9_code))),
    icd_resp_failure = as.integer(any(grepl("^51881|^51882|^51884", icd9_code))),
    icd_pneumonia = as.integer(any(grepl("^48[0-6]", icd9_code))),
    icd_diabetes = as.integer(any(grepl("^250", icd9_code))),
    icd_liver_disease = as.integer(any(grepl("^571", icd9_code))),
    icd_malignancy = as.integer(any(grepl("^1[4-9][0-9]|^20[0-8]", icd9_code))),
    icd_stroke = as.integer(any(grepl("^43[0-8]", icd9_code))),
    icd_mi = as.integer(any(grepl("^410", icd9_code))),
    icd_copd = as.integer(any(grepl("^496|^49[0-2]|^494", icd9_code))),
    .groups = "drop"
  ) %>%
  as.data.table()

# DRG codes for severity (1-4 scale, higher = worse)
drg <- fread("DRGCODES.csv")

drg_summaries <- drg %>%
  filter(!is.na(drg_mortality)) %>%
  group_by(hadm_id) %>%
  summarize(
    drg_mortality_score = max(drg_mortality, na.rm = TRUE),
    drg_severity_score = max(drg_severity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    drg_mortality_score = ifelse(is.infinite(drg_mortality_score), NA, drg_mortality_score),
    drg_severity_score = ifelse(is.infinite(drg_severity_score), NA, drg_severity_score)
  ) %>%
  as.data.table()

message(sprintf("  Diagnoses: %d admissions with ICD codes", nrow(dx_summaries)))

# ==============================================================================
# SECTION 6: PROCEDURES (ICD-9 & CPT)
# ==============================================================================
message("=== SECTION 6: Loading procedure data ===")

proc_icd <- fread("PROCEDURES_ICD.csv")
cptevents <- fread("CPTEVENTS.csv")

# Key procedure flags - these serve as both severity indicators and treatments
proc_summaries <- proc_icd %>%
  group_by(hadm_id) %>%
  summarize(
    n_procedures = n(),
    is_ventilated = as.integer(any(grepl("^967", icd9_code))),       # Mechanical ventilation
    received_dialysis = as.integer(any(icd9_code == "3995")),        # Hemodialysis
    has_central_line = as.integer(any(grepl("^3893|^3897", icd9_code))),
    received_transfusion = as.integer(any(grepl("^990", icd9_code))),
    .groups = "drop"
  ) %>%
  as.data.table()

cpt_summaries <- cptevents %>%
  group_by(hadm_id) %>%
  summarize(
    total_cpt_procedures = n(),
    n_distinct_cpt = n_distinct(cpt_cd),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Procedures: %d admissions with ICD, %d with CPT",
                nrow(proc_summaries), nrow(cpt_summaries)))

# ==============================================================================
# SECTION 7: VASOPRESSORS & INPUTS (INPUTEVENTS_MV)
# ==============================================================================
message("=== SECTION 7: Loading vasopressor data ===")

inputevents_mv <- fread("INPUTEVENTS_MV.csv")

# Vasopressor use - one of the strongest predictors of ICU mortality
vasopressor_use <- inputevents_mv %>%
  filter(itemid %in% vasopressor_ids_mv) %>%
  group_by(hadm_id) %>%
  summarize(
    received_vasopressors = 1L,
    n_vasopressor_doses = n(),
    vasopressor_duration_hrs = sum(as.numeric(difftime(
      as.POSIXct(endtime), as.POSIXct(starttime), units = "hours"
    )), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.table()

# Total IV fluid volume (non-vasopressor)
fluid_summary <- inputevents_mv %>%
  filter(!(itemid %in% vasopressor_ids_mv)) %>%
  group_by(hadm_id) %>%
  summarize(
    total_input_volume_ml = sum(amount, na.rm = TRUE),
    n_input_events = n(),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Vasopressors: %d admissions received vasopressors",
                nrow(vasopressor_use)))

# ==============================================================================
# SECTION 8: OUTPUT EVENTS (Urine Output)
# ==============================================================================
message("=== SECTION 8: Loading output events ===")

outputevents <- fread("OUTPUTEVENTS.csv")

output_summary <- outputevents %>%
  group_by(hadm_id) %>%
  summarize(
    total_output_ml = sum(as.numeric(value), na.rm = TRUE),
    n_output_events = n(),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Outputs: %d admissions with output data", nrow(output_summary)))

# ==============================================================================
# SECTION 9: SERVICES (Medical vs Surgical)
# ==============================================================================
message("=== SECTION 9: Loading service assignments ===")

services <- fread("SERVICES.csv")

service_summary <- services %>%
  group_by(hadm_id) %>%
  arrange(transfertime) %>%
  summarize(
    first_service = first(curr_service),
    n_service_changes = n() - 1,
    is_surgical = as.integer(any(grepl("SURG|ORTHO|GYN|NSURG|CSURG|VSURG|TSURG",
                                       curr_service, ignore.case = TRUE))),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Services: %d admissions with service data", nrow(service_summary)))

# ==============================================================================
# SECTION 10: TRANSFERS
# ==============================================================================
message("=== SECTION 10: Loading transfer data ===")

transfers <- fread("TRANSFERS.csv")

transfer_summary <- transfers %>%
  group_by(hadm_id) %>%
  summarize(
    n_transfers = n() - 1,
    n_unique_units = n_distinct(curr_careunit, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.table()

# ==============================================================================
# SECTION 11: PRESCRIPTIONS
# ==============================================================================
message("=== SECTION 11: Loading prescription data ===")

prescriptions <- fread("PRESCRIPTIONS.csv")

rx_summary <- prescriptions %>%
  group_by(hadm_id) %>%
  summarize(
    n_prescriptions = n(),
    n_unique_drugs = n_distinct(drug),
    n_unique_drug_classes = n_distinct(formulary_drug_cd, na.rm = TRUE),
    on_antibiotic = as.integer(any(grepl("CILLIN|CEPH|MYCIN|FLOXACIN|CYCLINE|SULFA",
                                         drug, ignore.case = TRUE))),
    on_anticoagulant = as.integer(any(grepl("HEPARIN|WARFARIN|ENOXAPARIN|RIVAROXABAN",
                                            drug, ignore.case = TRUE))),
    on_insulin = as.integer(any(grepl("INSULIN", drug, ignore.case = TRUE))),
    on_sedative = as.integer(any(grepl("PROPOFOL|MIDAZOLAM|LORAZEPAM|FENTANYL|MORPHINE",
                                       drug, ignore.case = TRUE))),
    .groups = "drop"
  ) %>%
  as.data.table()

message(sprintf("  Prescriptions: %d admissions with medication data", nrow(rx_summary)))

# ==============================================================================
# SECTION 12: CALLOUTS (Discharge Planning)
# ==============================================================================
message("=== SECTION 12: Loading callout data ===")

callout <- fread("CALLOUT.csv")

callout_summary <- callout %>%
  group_by(hadm_id) %>%
  summarize(
    n_callouts = n(),
    callout_requested = 1L,
    .groups = "drop"
  ) %>%
  as.data.table()

# ==============================================================================
# SECTION 13: MERGE INTO COMPREHENSIVE ANALYTIC DATASET
# ==============================================================================
# This builds the STATIC BASELINE MATRIX X
# All variables here are considered "known at time zero" (ICU admission)
# Even though some (like labs) may span the entire stay, they're aggregated
# into single values per admission for the purpose of this matrix.

message("\n=== SECTION 13: Merging into analytic dataset ===")

analytic_df <- cohort

# Sequentially merge all summary tables
analytic_df <- merge(analytic_df, vital_summaries, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, lab_summaries, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, blood_culture, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, any_positive_culture, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, dx_summaries, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, drg_summaries, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, proc_summaries, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, cpt_summaries, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, vasopressor_use, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, fluid_summary, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, output_summary, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, service_summary, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, transfer_summary, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, rx_summary, by = "hadm_id", all.x = TRUE)
analytic_df <- merge(analytic_df, callout_summary, by = "hadm_id", all.x = TRUE)

# ==============================================================================
# SECTION 14: CREATE ready_df (STATIC BASELINE MATRIX X)
# ==============================================================================
# This is the final patient-admission level dataset
# Structure: N rows (admissions) x P columns (covariates)
# Each row represents ONE admission at "time zero"

message("\n=== SECTION 14: Creating ready_df (static X matrix) ===")

ready_df <- analytic_df %>%
  as_tibble() %>%
  mutate(
    # -------------------------------------------------------------------------
    # SURVIVAL OUTCOMES
    # -------------------------------------------------------------------------
    # end_time: death if died, discharge if censored
    end_time = coalesce(deathtime, dischtime),

    # Time from ICU admission (Time Zero) to event/censoring
    time_to_event_hours = as.numeric(difftime(end_time, icu_intime, units = "hours")),
    time_to_event_days = time_to_event_hours / 24,

    # Event indicator: delta_i = 1 if died, 0 if censored (survived to discharge)
    event_status = as.integer(!is.na(deathtime)),

    # ICU-specific mortality (died during ICU stay, not just hospitalization)
    icu_mortality = as.integer(!is.na(deathtime) & deathtime <= icu_outtime),

    # -------------------------------------------------------------------------
    # HANDLE MISSING BINARY INDICATORS
    # -------------------------------------------------------------------------
    # Justification: Absence of a coded diagnosis/procedure/medication means
    # that condition was not present or that intervention was not performed
    across(starts_with("icd_"), ~ ifelse(is.na(.), 0L, .)),
    across(starts_with("is_"), ~ ifelse(is.na(.), 0L, .)),
    across(starts_with("received_"), ~ ifelse(is.na(.), 0L, .)),
    across(starts_with("has_"), ~ ifelse(is.na(.), 0L, .)),
    across(starts_with("on_"), ~ ifelse(is.na(.), 0L, .)),
    positive_blood_culture = ifelse(is.na(positive_blood_culture), 0L, positive_blood_culture),
    any_positive_culture = ifelse(is.na(any_positive_culture), 0L, any_positive_culture),
    callout_requested = ifelse(is.na(callout_requested), 0L, callout_requested),

    # -------------------------------------------------------------------------
    # HANDLE MISSING COUNTS
    # -------------------------------------------------------------------------
  n_diagnoses = ifelse(is.na(n_diagnoses), 0L, n_diagnoses),
    n_procedures = ifelse(is.na(n_procedures), 0L, n_procedures),
    n_prescriptions = ifelse(is.na(n_prescriptions), 0L, n_prescriptions),
    n_transfers = ifelse(is.na(n_transfers), 0L, n_transfers),
    n_callouts = ifelse(is.na(n_callouts), 0L, n_callouts),
    n_icu_stays = ifelse(is.na(n_icu_stays), 1L, n_icu_stays),
    total_cpt_procedures = ifelse(is.na(total_cpt_procedures), 0L, total_cpt_procedures),
    n_vasopressor_doses = ifelse(is.na(n_vasopressor_doses), 0L, n_vasopressor_doses),
    vasopressor_duration_hrs = ifelse(is.na(vasopressor_duration_hrs), 0, vasopressor_duration_hrs),
    # DRG scores default to 1 (lowest severity) - conservative choice
    drg_mortality_score = ifelse(is.na(drg_mortality_score), 1, drg_mortality_score),
    drg_severity_score = ifelse(is.na(drg_severity_score), 1, drg_severity_score),
    # -------------------------------------------------------------------------
    # DERIVED VARIABLES
    # -------------------------------------------------------------------------
    # Comorbidity burden: sum of all ICD flags (0-12 range)
    comorbidity_count = icd_sepsis + icd_heart_failure + icd_aki + icd_ckd +
      icd_resp_failure + icd_pneumonia + icd_diabetes + icd_liver_disease +
      icd_malignancy + icd_stroke + icd_mi + icd_copd,

    # -------------------------------------------------------------------------
    # SCALED CONTINUOUS VARIABLES
    # -------------------------------------------------------------------------
    # Justification: Standardization (mean=0, sd=1) allows consistent
    # weakly informative priors like N(0, 2.5) across all predictors
    scaled_age = scale(age)[, 1],
    scaled_los_icu = scale(los_icu)[, 1],
    scaled_n_diagnoses = scale(n_diagnoses)[, 1],
    scaled_n_procedures = scale(n_procedures)[, 1],
    scaled_n_prescriptions = scale(n_prescriptions)[, 1],
    scaled_comorbidity = scale(comorbidity_count)[, 1],

    # -------------------------------------------------------------------------
    # FACTOR VARIABLES FOR HIERARCHICAL MODELING
    # -------------------------------------------------------------------------
    subject_id = as.factor(subject_id),
    hadm_id = as.factor(hadm_id),
    first_careunit = as.factor(first_careunit),
    first_service = as.factor(first_service)
  ) %>%
  # Filter to valid survival data (positive time to event)
  filter(
    !is.na(time_to_event_hours),
    time_to_event_hours > 0
  ) %>%
  # -------------------------------------------------------------------------
  # MULTICOLLINEARITY & HIGH MISSINGNESS REDUCTION 
  # -------------------------------------------------------------------------
  select(
    # 1. Drop multicollinear aggregates (keep only mean/sd/n_obs)
    -ends_with("_min_val"),
    -ends_with("_max_val"),
    -ends_with("_first_val"),
    -ends_with("_last_val"),
    
    # NOTE THIS IS ONLY ON DEMO DATA - MAY WANT TO KEEP IN FULL DATASET
    # 2. Drop covariates with > 25% missingness to preserve cohort size
    -starts_with("dbp_"),               # ~81% missing 
    -starts_with("map_"),               # ~58% missing
    -starts_with("sbp_"),               # ~57% missing
    -starts_with("total_input_"),       # ~45% missing
    -starts_with("n_input_"),           # ~45% missing
    -starts_with("n_positive_spec"),    # ~44% missing
    -starts_with("n_unique_org"),       # ~44% missing
    -dod_hosp,                          # ~37% missing
    -starts_with("bilirubin_"),         # ~32% missing
    -starts_with("pao2_"),              # ~29% missing
    -starts_with("pco2_"),              # ~29% missing
    -edregtime, -edouttime,             # ~28% missing
    -starts_with("n_blood_cultures"),   # ~27% missing
    -starts_with("n_positive_cult"),    # ~27% missing
    -starts_with("n_organisms"),        # ~27% missing
    -starts_with("ph_"),                # ~26% missing
    # Lactate is at 21% missing; decide if you want to keep or drop it.
    # If dropping lactate entirely uncomment below:
    -starts_with("lactate_"),
    -dod_ssn                            # ~19% missing
  ) %>%
  # ===========================================================================
  # SECTION 14b: CLEAN HIGH-CARDINALITY CHARACTER VARIABLES
  # ===========================================================================
  # Remove variables problematic for BART modeling:
  # - High-cardinality text fields (many unique levels create sparse splits)
  # - POSIXct/datetime columns (should use derived numeric features instead)
  # - Redundant or leakage variables
  select(
    -admission_location,     # High-cardinality text (7+ categories with long strings)
    -discharge_location,     # High-cardinality text (15+ categories) - also LEAKAGE (known only at discharge)
    -language,               # High-cardinality character (many languages + empty strings)
    -religion,               # High-cardinality character (many religions + empty strings)
    -diagnosis,              # Free-text primary diagnosis - extremely high cardinality
    -dob,                    # POSIXct datetime - use 'age' instead
    -dod,                    # POSIXct datetime - use 'event_status' instead
    -expire_flag,            # Outcome indicator - LEAKAGE (redundant with event_status)
    -has_chartevents_data    # Binary indicator - all ICU patients have chart data
  ) %>%
  # ===========================================================================
  # SECTION 14c: DISCRETIZE CONTINUOUS VARIABLES INTO FACTORS
  # ===========================================================================
  # Following the pattern from data_rebin.R: as.factor(ifelse(...))
  # Creates binary or simple categorical factors for BART
  #
  # Clinical thresholds based on:
  # - SOFA score criteria
  # - SIRS criteria
  # - Standard clinical reference ranges
  # ===========================================================================
  mutate(
# -------------------------------------------------------------------------
    # DEMOGRAPHICS & STAY CHARACTERISTICS
    # -------------------------------------------------------------------------
    age_group = factor(case_when(
      age >= 80 ~ "Very Old",
      age >= 65 ~ "Elderly",
      TRUE      ~ "Adult"
    ), levels = c("Adult", "Elderly", "Very Old")),

    los_group = factor(case_when(
      los_icu >= 7 ~ "Extended",
      los_icu >= 3 ~ "Prolonged",
      TRUE         ~ "Short"
    ), levels = c("Short", "Prolonged", "Extended")),

    multiple_icu_stays = as.factor(ifelse(n_icu_stays > 1, 1, 0)),
    gender = as.factor(gender),
    # -------------------------------------------------------------------------
    # VITAL SIGNS
    # -------------------------------------------------------------------------
    gcs_status = factor(case_when(
      gcs_mean_val <= 8 ~ "Severe Impairment",
      gcs_mean_val < 15 ~ "Impaired",
      TRUE              ~ "Normal"
    ), levels = c("Normal", "Impaired", "Severe Impairment")),

    hr_status = factor(case_when(
      hr_mean_val > 100 ~ "Tachycardia",
      hr_mean_val < 60  ~ "Bradycardia",
      TRUE              ~ "Normal"
    ), levels = c("Normal", "Bradycardia", "Tachycardia")),

    rr_status = factor(case_when(
      resp_rate_mean_val > 30 ~ "Severe Tachypnea",
      resp_rate_mean_val > 20 ~ "Tachypnea",
      TRUE                    ~ "Normal"
    ), levels = c("Normal", "Tachypnea", "Severe Tachypnea")),

    spo2_status = factor(case_when(
      spo2_mean_val < 90 ~ "Severe Hypoxia",
      spo2_mean_val < 95 ~ "Hypoxia",
      TRUE               ~ "Normal"
    ), levels = c("Normal", "Hypoxia", "Severe Hypoxia")),

    temp_status = factor(case_when(
      temp_c_mean_val > 38 ~ "Fever",
      temp_c_mean_val < 36 ~ "Hypothermia",
      TRUE                 ~ "Normal"
    ), levels = c("Normal", "Hypothermia", "Fever")),

    # -------------------------------------------------------------------------
    # LABORATORY VALUES
    # -------------------------------------------------------------------------
    bun_status = factor(case_when(
      bun_mean_val > 40 ~ "Very High",
      bun_mean_val > 20 ~ "High",
      TRUE              ~ "Normal"
    ), levels = c("Normal", "High", "Very High")),

    creat_status = factor(case_when(
      creatinine_mean_val >= 2.0 ~ "Dysfunction",
      creatinine_mean_val >= 1.2 ~ "Elevated",
      TRUE                       ~ "Normal"
    ), levels = c("Normal", "Elevated", "Dysfunction")),

    glucose_status = factor(case_when(
      glucose_mean_val > 180 ~ "Hyperglycemia",
      glucose_mean_val < 70  ~ "Hypoglycemia",
      TRUE                   ~ "Normal"
    ), levels = c("Normal", "Hypoglycemia", "Hyperglycemia")),

    hgb_status = factor(case_when(
      hemoglobin_mean_val < 7  ~ "Severe Anemia",
      hemoglobin_mean_val < 10 ~ "Anemia",
      TRUE                     ~ "Normal"
    ), levels = c("Normal", "Anemia", "Severe Anemia")),

    plt_status = factor(case_when(
      platelets_mean_val < 50  ~ "Severe Thrombocytopenia",
      platelets_mean_val < 150 ~ "Thrombocytopenia",
      TRUE                     ~ "Normal"
    ), levels = c("Normal", "Thrombocytopenia", "Severe Thrombocytopenia")),

    k_status = factor(case_when(
      potassium_mean_val > 5.0 ~ "Hyperkalemia",
      potassium_mean_val < 3.5 ~ "Hypokalemia",
      TRUE                     ~ "Normal"
    ), levels = c("Normal", "Hypokalemia", "Hyperkalemia")),

    na_status = factor(case_when(
      sodium_mean_val > 145 ~ "Hypernatremia",
      sodium_mean_val < 135 ~ "Hyponatremia",
      TRUE                  ~ "Normal"
    ), levels = c("Normal", "Hyponatremia", "Hypernatremia")),

    wbc_status = factor(case_when(
      wbc_mean_val > 12 ~ "Leukocytosis",
      wbc_mean_val < 4  ~ "Leukopenia",
      TRUE              ~ "Normal"
    ), levels = c("Normal", "Leukopenia", "Leukocytosis")),

    # -------------------------------------------------------------------------
    # VARIABILITY & COMPLEXITY
    # -------------------------------------------------------------------------
    hr_unstable = as.factor(ifelse(hr_sd_val > 15, 1, 0)),
    temp_unstable = as.factor(ifelse(temp_c_sd_val > 0.5, 1, 0)),

    comorbidity_level = factor(case_when(
      comorbidity_count >= 5 ~ "Very High",
      comorbidity_count >= 3 ~ "High",
      TRUE                   ~ "Low/Moderate"
    ), levels = c("Low/Moderate", "High", "Very High")),

    complexity_status = factor(case_when(
      n_prescriptions > 20 | n_diagnoses > 10 | n_procedures > 5 ~ "High Complexity",
      TRUE ~ "Standard"
    ), levels = c("Standard", "High Complexity")),

    # -------------------------------------------------------------------------
    # EXISTING CATEGORICAL VARIABLES - Ensure all are factors
    # -------------------------------------------------------------------------
    admission_type = as.factor(admission_type),
    marital_status = as.factor(marital_status),
    ethnicity = as.factor(ethnicity),
    insurance = as.factor(insurance),
    first_careunit = as.factor(first_careunit),
    first_service = as.factor(first_service),

    # -------------------------------------------------------------------------
    # EXISTING BINARY INDICATORS - Ensure all are factors (already 0/1)
    # -------------------------------------------------------------------------
    icd_sepsis = as.factor(icd_sepsis),
    icd_heart_failure = as.factor(icd_heart_failure),
    icd_aki = as.factor(icd_aki),
    icd_ckd = as.factor(icd_ckd),
    icd_resp_failure = as.factor(icd_resp_failure),
    icd_pneumonia = as.factor(icd_pneumonia),
    icd_diabetes = as.factor(icd_diabetes),
    icd_liver_disease = as.factor(icd_liver_disease),
    icd_malignancy = as.factor(icd_malignancy),
    icd_stroke = as.factor(icd_stroke),
    icd_mi = as.factor(icd_mi),
    icd_copd = as.factor(icd_copd),
    is_ventilated = as.factor(is_ventilated),
    received_dialysis = as.factor(received_dialysis),
    has_central_line = as.factor(has_central_line),
    received_transfusion = as.factor(received_transfusion),
    received_vasopressors = as.factor(received_vasopressors),
    is_surgical = as.factor(is_surgical),
    on_antibiotic = as.factor(on_antibiotic),
    on_anticoagulant = as.factor(on_anticoagulant),
    on_insulin = as.factor(on_insulin),
    on_sedative = as.factor(on_sedative),
    positive_blood_culture = as.factor(positive_blood_culture),
    any_positive_culture = as.factor(any_positive_culture),
    callout_requested = as.factor(callout_requested)
  ) %>%
  # # ===========================================================================
  # # SECTION 14d: DROP MULTICOLLINEARITY & HIGH-MISSINGNESS VARIABLES
  # # Justification:
  # # - Multicollinear aggregates (mean/sd/n_obs capture the same info as min/max/first/last)
  # # - High missingness variables reduce sample size and may not add enough predictive value to justify the loss
  # # ===========================================================================
  select(
    -gcs_sd_val, -gcs_n_obs,
    -hr_mean_val, -hr_sd_val, -hr_n_obs,
    -resp_rate_mean_val, -resp_rate_sd_val, -resp_rate_n_obs,
    -spo2_mean_val, -spo2_sd_val, -spo2_n_obs,
    -temp_c_mean_val, -temp_c_sd_val, -temp_c_n_obs,
    -bun_mean_val, -bun_n_labs,
    -creatinine_mean_val, -creatinine_n_labs,
    -glucose_mean_val, -glucose_n_labs,
    -hemoglobin_mean_val, -hemoglobin_n_labs,
    -platelets_mean_val, -platelets_n_labs,
    -potassium_mean_val, -potassium_n_labs,
    -sodium_mean_val, -sodium_n_labs,
    -wbc_mean_val, -wbc_n_labs,
  ) %>% 
  select(
    -starts_with("scaled_")
  )



message(sprintf("  ready_df: %d admissions x %d variables",
                nrow(ready_df), ncol(ready_df)))

# ==============================================================================
# SECTION 15: CREATE longitudinal_df (TIME-VARYING MATRIX Z)
# ==============================================================================
# This is the time-varying covariate dataset
# Structure: One row per vital sign measurement occasion
# Z_{i,t} = covariates for patient i at time t
# Used for: joint longitudinal-survival models, trajectory analysis

message("\n=== SECTION 15: Creating longitudinal_df (time-varying Z matrix) ===")

# Select only the columns we need from cohort to avoid duplicate columns
cohort_for_merge <- cohort %>%
  select(
    subject_id, hadm_id,
    # Timestamps needed for time calculations
    icu_intime, icu_outtime, deathtime, dischtime,
    # Demographics for clustering
    first_careunit, age, gender,
    # Outcome info
    hospital_expire_flag
  )

# Merge vitals with cohort info
# Using specific columns avoids .x/.y suffix issues
longitudinal_df <- vitals_wide %>%
  left_join(cohort_for_merge, by = c("subject_id", "hadm_id")) %>%
  as_tibble() %>%
  mutate(
    # Convert timestamps
    charttime = as.POSIXct(charttime),
    icu_intime = as.POSIXct(icu_intime),

    # -------------------------------------------------------------------------
    # TIME INDEX: hours since ICU admission
    # -------------------------------------------------------------------------
    # This is the "t" in Z_{i,t}
    hours_since_icu_admit = as.numeric(difftime(charttime, icu_intime, units = "hours")),

    # -------------------------------------------------------------------------
    # SURVIVAL OUTCOME INFO (for joint modeling)
    # -------------------------------------------------------------------------
    deathtime = as.POSIXct(deathtime),
    dischtime = as.POSIXct(dischtime),
    end_time = coalesce(deathtime, dischtime),
    event_status = as.integer(!is.na(deathtime)),

    # Time to event for this patient (constant within patient)
    time_to_event_hours = as.numeric(difftime(end_time, icu_intime, units = "hours")),

    # -------------------------------------------------------------------------
    # SCALED VITALS
    # -------------------------------------------------------------------------
    scaled_hr = scale(hr)[, 1],
    scaled_sbp = scale(sbp)[, 1],
    scaled_dbp = scale(dbp)[, 1],
    scaled_map = scale(map)[, 1],
    scaled_resp_rate = scale(resp_rate)[, 1],
    scaled_spo2 = scale(spo2)[, 1],
    scaled_temp_c = scale(temp_c)[, 1],
    scaled_time = scale(hours_since_icu_admit)[, 1],

    # Factor IDs for random effects
    subject_id = as.factor(subject_id),
    hadm_id = as.factor(hadm_id),
    first_careunit = as.factor(first_careunit)
  ) %>%
  # -------------------------------------------------------------------------
  # FILTER TO VALID TIME WINDOW
  # -------------------------------------------------------------------------
  # Only keep measurements:
  #   1. After ICU admission (hours_since_icu_admit >= 0)
  #   2. Within first 7 days (168 hours) - focus on acute phase
  filter(
    hours_since_icu_admit >= 0,
    hours_since_icu_admit <= 168
  )

message(sprintf("  longitudinal_df: %d observations x %d variables",
                nrow(longitudinal_df), ncol(longitudinal_df)))

# ==============================================================================
# SECTION 16: CREATE person_period_df (DISCRETE-TIME FORMAT)
# ==============================================================================
# This is the discrete-time survival dataset for Polya-Gamma augmentation
# Structure: One row per patient i per discrete time period t
#
# DISCRETE-TIME SURVIVAL FRAMEWORK:
# ---------------------------------
# Continuous time is discretized into 24-hour periods (days)
# For patient i observed for T_i periods:
#   - y_{i,t} = 1 if event occurs in period t, 0 otherwise
#   - Patient dies in period 3:  y_{i,1}=0, y_{i,2}=0, y_{i,3}=1
#   - Patient censored in period 3: y_{i,1}=0, y_{i,2}=0, y_{i,3}=0
#
# WHY DISCRETE-TIME?
# -----------------
# 1. Transforms survival into logistic regression (binary outcomes)
# 2. Enables Polya-Gamma augmentation for conjugate Gibbs sampling
# 3. Natural for time-varying covariates (one covariate value per period)
# 4. Computational efficiency for large datasets like MIMIC-III

message("\n=== SECTION 16: Creating person_period_df (discrete-time format) ===")
# -------------------------------------------------------------------------
# STEP 1: Aggregate longitudinal data (Z) into 24-hour periods
# -------------------------------------------------------------------------
longitudinal_daily <- longitudinal_df %>%
  mutate(
    # Assign each measurement to a discrete period (day)
    # period 1 = hours 0-24, period 2 = hours 24-48, etc.
    period = ceiling(hours_since_icu_admit / 24)
  ) %>%
  # Keep first 7 days only
  filter(period >= 1 & period <= 7) %>%
  group_by(hadm_id, period) %>%
  # Compute period-level summaries (mean, min, max, count)
  summarize(
    across(
      c(hr, map, resp_rate, spo2, temp_c),
      list(
        mean = ~ mean(., na.rm = TRUE),
        min  = ~ if(all(is.na(.))) NA_real_ else min(., na.rm = TRUE),
        max  = ~ if(all(is.na(.))) NA_real_ else max(., na.rm = TRUE),
        n    = ~ sum(!is.na(.))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# STEP 2 & 3: Expand static variables (X) to person-period format
# -------------------------------------------------------------------------
person_period_df <- ready_df %>%
  mutate(
    # Total periods at risk: ceiling of ICU LOS, capped at 7 days
    max_period = pmin(ceiling(los_icu), 7),
    max_period = pmax(max_period, 1)  # At least 1 period
  ) %>%
  # Duplicate each row 'max_period' times, creating a period counter
  uncount(max_period, .id = "period")

# -------------------------------------------------------------------------
# STEP 4: Construct discrete-time outcome y_{i,t}
# -------------------------------------------------------------------------
person_period_df <- person_period_df %>%
  group_by(hadm_id) %>%
  mutate(
    max_period = max(period),
    is_last_period = (period == max_period),
    # DISCRETE-TIME OUTCOME:
    # 1 only if (a) this is the last period AND (b) the patient died
    y_it = as.integer(is_last_period & event_status == 1)
  ) %>%
  ungroup() %>%
  select(-is_last_period)

# -------------------------------------------------------------------------
# STEP 5 & 6: Merge Z_{i,t} and Handle Missing Data
# -------------------------------------------------------------------------
person_period_df <- person_period_df %>%
  left_join(longitudinal_daily, by = c("hadm_id", "period")) %>%
  group_by(hadm_id) %>%
  mutate(
    # MISSINGNESS INDICATORS (Informative sparsity)
    hr_missing = as.integer(is.na(hr_mean)),
    temp_c_missing = as.integer(is.na(temp_c_mean)),
    map_missing = as.integer(is.na(map_mean)),
    
    # LOCF IMPUTATION: carry previous day's vitals forward if unmeasured
    across(matches("^(hr|map|resp_rate|spo2|temp_c)_(mean|min|max|n)$"), 
           ~ zoo::na.locf(., na.rm = FALSE))
  ) %>%
  ungroup() %>%
  # POPULATION MEDIAN for remaining NAs (patients with no prior values)
  mutate(
    across(matches("^(hr|map|resp_rate|spo2|temp_c)_(mean|min|max|n)$"), 
           ~ replace_na(., median(., na.rm = TRUE)))
  )

# # -------------------------------------------------------------------------
# # STEP 6: Handle missing time-varying covariates
# # -------------------------------------------------------------------------
# # Strategy:
# #   1. Create missingness indicators (informative in ICU context)
# #   2. Apply LOCF (last observation carried forward) within patient
# #   3. Fill remaining NAs with population median
# 
# person_period_df <- person_period_df %>%
#   group_by(hadm_id) %>%
#   mutate(
#     # MISSINGNESS INDICATORS
#     # Justification: In ICU, missing vitals may indicate patient stability
#     # (nurse not recording as frequently) - this is informative
#     hr_missing = as.integer(is.na(hr_mean)),
#     temp_c_missing = as.integer(is.na(temp_c_mean)),
#     map_missing = as.integer(is.na(map_mean)),
# 
#     # LOCF IMPUTATION
#     # Justification: Use most recent value as best guess for current period
#     # Note: This artificially shrinks variance - consider Bayesian imputation
#     hr_mean = zoo::na.locf(hr_mean, na.rm = FALSE),
#     temp_c_mean = zoo::na.locf(temp_c_mean, na.rm = FALSE),
#     map_mean = zoo::na.locf(map_mean, na.rm = FALSE),
#     resp_rate_mean = zoo::na.locf(resp_rate_mean, na.rm = FALSE),
#     spo2_mean = zoo::na.locf(spo2_mean, na.rm = FALSE)
#   ) %>%
#   ungroup() %>%
#   # POPULATION MEDIAN for remaining NAs (patients with no prior values)
#   mutate(
#     hr_mean = replace_na(hr_mean, median(hr_mean, na.rm = TRUE)),
#     temp_c_mean = replace_na(temp_c_mean, median(temp_c_mean, na.rm = TRUE)),
#     map_mean = replace_na(map_mean, median(map_mean, na.rm = TRUE)),
#     resp_rate_mean = replace_na(resp_rate_mean, median(resp_rate_mean, na.rm = TRUE)),
#     spo2_mean = replace_na(spo2_mean, median(spo2_mean, na.rm = TRUE))
#   )
# 
# # Convert IDs to factors for hierarchical modeling
# person_period_df <- person_period_df %>%
#   mutate(
#     hadm_id = as.factor(hadm_id),
#     subject_id = as.factor(subject_id),
#     first_careunit = as.factor(first_careunit)
#   )
# 
# message(sprintf("  person_period_df: %d rows from %d admissions",
#                 nrow(person_period_df), length(unique(person_period_df$hadm_id))))


# ==============================================================================
# SECTION 17: FINAL CLEANING OF ready_df
# ==============================================================================
# Justification:
# - Remove any remaining variables that are IDs, outcomes, text, time-to-event, or high-cardinality
#   text that could cause issues for BART modeling

base.path <- "/Users/adamkurth/Documents/RStudio/bayes-mimic-survival"
data.path <- file.path(base.path, "data/mimic-iii-clinical-database-demo-1.4/")
out.path <- file.path(base.path, "data/rds")
setwd(base.path); 
# source("load.r")


exclude <- c(
  # --- PREVIOUS EXCLUSIONS ---
  "subject_id", "icustay_id",                   
  "event_status", "hospital_expire_flag",                  
  "time_to_event_hours", "time_to_event_days",             
  "admittime", "dischtime", "deathtime",                   
  "icu_intime", "icu_outtime", "end_time",                 
  "discharge_location", "expire_flag",                     
  "dob", "dod", "dod_hosp", "dod_ssn",                     
  "diagnosis", "marital_status", "religion", "language",
  
  # --- NEW: EXACT MULTICOLLINEARITY (SCALED DUPLICATES) ---
  "scaled_age", "scaled_los_icu", "scaled_n_diagnoses", 
  "scaled_n_procedures", "scaled_n_prescriptions", 
  "scaled_comorbidity",
  
  # --- NEW: STRUCTURAL REDUNDANCY ---
  "total_icu_los",          # Redundant with los_icu
  "comorbidity_count",      # Exact sum of the icd_* binary flags
  "received_vasopressors"   # Redundant with vasopressor_duration_hrs
)

# Note: 'icu_mortality' was explicitly removed from the exclude list 
# above so it can be extracted as the target 'y' below.

# Apply exclusions and drop the EHR artifact _n_obs / _n_labs columns
dat <- ready_df %>% 
  select(-any_of(exclude)) %>%
  select(-ends_with("_n_obs"), -ends_with("_n_labs")) %>% 
  select(-any_of(c("total_output_ml", "n_output_events", "n_distinct_cpt", 
                   "n_unique_drugs", "n_unique_drug_classes",
                   "edregtime", "edouttime"))) %>%
  mutate_if(is.character, as.factor)

# ethnicity
dat$ethnicity <- recode(dat$ethnicity,
  "WHITE" = "1",
  "BLACK/AFRICAN AMERICAN" = "2",
  "ASIAN" = "3",
  "AMERICAN INDIAN/ALASKA NATIVE FEDERALLY RECOGNIZED TRIBE" = "4",
  "HISPANIC OR LATINO" = "5",
  "HISPANIC/LATINO - PUERTO RICAN" = "5",
  "OTHER" = "6",
  "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" = "7",
  "UNKNOWN/NOT SPECIFIED" = "0",
  "OTHER/UNKNOWN" = "0",
  "UNABLE TO OBTAIN" = "0",
)

dat$ethnicity[is.na(dat$ethnicity)] <- "0"
dat$ethnicity <- factor(dat$ethnicity, levels = c("0", "1", "2", "3", "4", "5", "6", "7"))

cat(" Before NA omission: ", nrow(dat), " rows\n") # 122
dat <- na.omit(dat)
cat(" After NA omission: ", nrow(na.omit(dat)), " rows\n") # 122

y <- dat$icu_mortality
X <- as.data.frame(dat %>% select(-icu_mortality, -hadm_id)) # Drop ID now
rownames(X) <- valid_hadm_ids # Optional: cleanly names the rows in your matrix

# -----------------------------------------------------------------------------
# CREATE S MATRIX
# -----------------------------------------------------------------------------
S <- person_period_df %>%
  # Filter to matched cohort FIRST, before any column drops
  filter(hadm_id %in% valid_hadm_ids) %>%
  select(-any_of(exclude)) %>%
  select(-ends_with("_n_obs"), -ends_with("_n_labs")) %>% 
  select(-any_of(c("total_output_ml", "n_output_events", "n_distinct_cpt", 
                   "n_unique_drugs", "n_unique_drug_classes",
                   "edregtime", "edouttime"))) %>%
  mutate_if(is.character, as.factor)

# Recode ethnicity in S to match X exactly
S$ethnicity <- recode(S$ethnicity,
                      "WHITE" = "1", "BLACK/AFRICAN AMERICAN" = "2", "ASIAN" = "3",
                      "AMERICAN INDIAN/ALASKA NATIVE FEDERALLY RECOGNIZED TRIBE" = "4",
                      "HISPANIC OR LATINO" = "5", "HISPANIC/LATINO - PUERTO RICAN" = "5",
                      "OTHER" = "6", "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" = "7",
                      "UNKNOWN/NOT SPECIFIED" = "0", "OTHER/UNKNOWN" = "0", "UNABLE TO OBTAIN" = "0"
)

S$ethnicity[is.na(S$ethnicity)] <- "0"
S$ethnicity <- factor(S$ethnicity, levels = c("0", "1", "2", "3", "4", "5", "6", "7"))

cat(" Before NA omission: ", nrow(S), " rows\n") # 434 rows
S <- na.omit(S)
cat(" After NA omission: ", nrow(na.omit(S)), " rows\n") # 424 rows

# Extract y vector for Polya-Gamma
y_S <- S$y_it

# Remove target labels and ID from the actual S design matrix
# Note: Keeping `period` is often useful for baseline hazard approximations
S_mat <- as.data.frame(S %>% select(-y_it, -icu_mortality, -hadm_id))

# ==============================================================================
# SECTION 17: SUMMARY AND OUTPUT
# ==============================================================================
message("\n============ DATA PREPARATION COMPLETE ============")
message(sprintf("  ready_df         : %d admissions x %d variables (Static X matrix)",
                nrow(ready_df), ncol(ready_df)))
message(sprintf("  longitudinal_df  : %d observations x %d variables (Time-varying Z matrix)",
                nrow(longitudinal_df), ncol(longitudinal_df)))
message(sprintf("  person_period_df : %d rows x %d variables (Discrete-time format)",
                nrow(person_period_df), ncol(person_period_df)))
message(sprintf("  Mortality rate   : %.1f%%",
                100 * mean(ready_df$event_status, na.rm = TRUE)))
message(sprintf("  'X' matrix        : %d admissions x %d predictors",
                nrow(X), ncol(X)))
message(sprintf("  'y' vector        : %d admissions with %d events (icu_mortality)",
                length(y), sum(y)))

# ==============================================================================
# SECTION 18: SAVE OUTPUT
# ==============================================================================
message("\n=== Saving output files ===")

saveRDS(ready_df, file = file.path(out.path, "ready_df.rds"))
saveRDS(longitudinal_df, file = file.path(out.path, "longitudinal_df.rds"))
saveRDS(person_period_df, file = file.path(out.path, "person_period_df.rds"))
saveRDS (list(X = X, y = y), file = file.path(out.path, "X_y.rds"))
saveRDS(list(S = S_mat, y = y_S), file = file.path(out.path, "S_y.rds"))


message(sprintf("  Files saved to: %s", out.path))
message("\nData objects available:")
message("  - ready_df: Patient-level static baseline matrix (X)")
message("  - longitudinal_df: Time-varying vital signs (Z_{i,t})")
message("  - person_period_df: Discrete-time survival format for Polya-Gamma")

setwd(base.path) # Return to base path



