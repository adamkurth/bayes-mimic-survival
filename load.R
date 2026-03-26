# load.R
# Robust data loading with path validation and error handling
# PHP2530 Bayesian Statistics Final Project
# Source this file to load ready_df and longitudinal_df into your environment

# ============================================================================
# CONFIGURATION - Centralized path management
# ============================================================================
# Detect project root - works in RStudio, sourced scripts, or command line
get_project_root <- function() {
  # Try rstudioapi first (interactive RStudio)
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  # Try sys.frame for sourced scripts
  if (sys.nframe() > 0) {
    for (i in seq_len(sys.nframe())) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile)) {
        return(dirname(normalizePath(ofile)))
      }
    }
  }
  # Fall back to current working directory or manual path
  if (file.exists("data/ready_df.rds")) {
    return(getwd())
  }
  # Manual fallback
  return("/Users/adamkurth/Documents/RStudio/bayes-mimic-survival")
}

root <- get_project_root()
data.rds.dir <- file.path(root, "data/rds")
# ============================================================================
# LOAD DATA WITH ERROR HANDLING
# ============================================================================
d.files <- list(
    ready_df = file.path(data.rds.dir, "ready_df.rds"),
    longitudinal_df = file.path(data.rds.dir, "longitudinal_df.rds"),
    person_period_df = file.path(data.rds.dir, "person_period_df.rds")
)
load.data <- function(filepath, name) {
    if (!file.exists(filepath)) {
        warning("File not found: ", filepath, "\n",
                        "Make sure data-prep.R has been run to generate this file.")
        return(NULL)
    }
    
    tryCatch(
        {
            data <- readRDS(filepath)
            message("✓ Loaded ", name, " (", nrow(data), " rows, ", 
                            ncol(data), " columns)")
            return(data)
        },
        error = function(e) {
            warning("Error loading ", name, ": ", e$message)
            return(NULL)
        }
    )
}

# Load all data files
ready_df <- load.data(d.files$ready_df, "ready_df")
longitudinal_df <- load.data(d.files$longitudinal_df, "longitudinal_df")
person_period_df <- load.data(d.files$person_period_df, "person_period_df")
# ============================================================================
# VALIDATION & SUMMARY
# ============================================================================

# Check that at least one dataset loaded successfully
if (is.null(ready_df) && is.null(longitudinal_df) && is.null(person_period_df)) {
    stop("Failed to load any data files. ",
             "Run data-prep.R first to generate the required .rds files.")
}

# Summary of loaded objects
message("\n========== DATA LOADED ==========")
if (!is.null(ready_df)) {
    message("ready_df: ", nrow(ready_df), " admissions")
    message("  Variables: ", paste(names(ready_df)[1:5], collapse=", "), ", ...")
}
if (!is.null(longitudinal_df)) {
    message("longitudinal_df: ", nrow(longitudinal_df), " measurements")
    message("  Variables: ", paste(names(longitudinal_df)[1:5], collapse=", "), ", ...")
}
if (!is.null(person_period_df)) {
    message("person_period_df: ", nrow(person_period_df), " person-period records")
    message("  Variables: ", paste(names(person_period_df)[1:5], collapse=", "), ", ...")
}
message("================================\n")

setwd(root)