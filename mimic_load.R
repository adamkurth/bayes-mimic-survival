# mimic_load.R
# ==============================================================================
# Load MIMIC-III survival objects into workspace
# ==============================================================================
# Loads the output of mimic_prep.R and makes objects ready for analysis.
#
# Objects loaded:
#   patients_df  — patient-level data frame (N rows)
#   dat          — list matching generate_synthetic_icu() interface
#                  ($N, $K, $P, $X, $y_obs, $delta, $delta_width)
#   pp_data      — person-period expanded data (from expand_person_period)
#                  ($pp, $X.ast, $y, $patient_id, $interval, $R, $K, $P, $N)
#
# After loading, you can immediately feed pp_data into Stan — no need to
# re-run expand_person_period().
#
# ==============================================================================

# --- Detect project root ---
get_project_root <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  if (sys.nframe() > 0) {
    for (i in seq_len(sys.nframe())) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
    }
  }
  if (file.exists("data/rds/mimic_dat.rds")) return(getwd())
  return("/Users/adamkurth/Documents/RStudio/bayes-mimic-survival")
}

root <- get_project_root()
rds.dir <- file.path(root, "data/rds")

# --- Load dat and patients_df ---
dat_path      <- file.path(rds.dir, "mimic_dat.rds")
patients_path <- file.path(rds.dir, "mimic_patients.rds")
pp_path       <- file.path(rds.dir, "mimic_pp_data.rds")

if (!file.exists(dat_path)) {
  stop("mimic_dat.rds not found at ", dat_path,
       "\nRun mimic_prep.R first to generate this file.")
}

dat         <- readRDS(dat_path)
patients_df <- readRDS(patients_path)

# --- Load pp_data (person-period expansion) ---
if (file.exists(pp_path)) {
  pp_data <- readRDS(pp_path)
  message(sprintf("  pp_data loaded: R=%d person-period rows", pp_data$R))
} else {
  message("  WARNING: mimic_pp_data.rds not found. You will need to run:")
  message("    source('demo/core.r')")
  message("    pp_data <- expand_person_period(dat)")
  pp_data <- NULL
}

# --- Summary ---
message("\n========== MIMIC DATA LOADED ==========")
message(sprintf("  dat$N = %d patients", dat$N))
message(sprintf("  dat$K = %d max intervals (DELTA = %d day)", dat$K, dat$delta_width))
message(sprintf("  dat$P = %d covariates", dat$P))
message(sprintf("  dat$X columns: %s", paste(colnames(dat$X), collapse = ", ")))
message(sprintf("  Events: %d (%.1f%%)", sum(dat$delta), 100 * mean(dat$delta)))
message(sprintf("  Censored: %d (%.1f%%)", sum(1 - dat$delta), 100 * mean(1 - dat$delta)))
if (!is.null(pp_data)) {
  message(sprintf("  pp_data$R = %d person-period rows", pp_data$R))
}
message("========================================\n")
