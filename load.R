# load.R
# Robust data loading with path validation and error handling
# PHP2530 Bayesian Statistics Final Project
# Source this file to load X, y, S_mat, and y_S into your environment

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
  if (file.exists("data/rds/X_y.rds")) {
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

# --- Intermediate data frames (already made, commented out) ---
# d.files <- list(
#     ready_df = file.path(data.rds.dir, "ready_df.rds"),
#     longitudinal_df = file.path(data.rds.dir, "longitudinal_df.rds"),
#     person_period_df = file.path(data.rds.dir, "person_period_df.rds")
# )
# load.data <- function(filepath, name) {
#     if (!file.exists(filepath)) {
#         warning("File not found: ", filepath, "\n",
#                         "Make sure data-prep.R has been run to generate this file.")
#         return(NULL)
#     }
#
#     tryCatch(
#         {
#             data <- readRDS(filepath)
#             message("✓ Loaded ", name, " (", nrow(data), " rows, ",
#                             ncol(data), " columns)")
#             return(data)
#         },
#         error = function(e) {
#             warning("Error loading ", name, ": ", e$message)
#             return(NULL)
#         }
#     )
# }
#
# # Load all data files
# ready_df <- load.data(d.files$ready_df, "ready_df")
# longitudinal_df <- load.data(d.files$longitudinal_df, "longitudinal_df")
# person_period_df <- load.data(d.files$person_period_df, "person_period_df")

# --- Final processed matrices: X_y.rds and S_y.rds ---
X_y.path <- file.path(data.rds.dir, "X_y.rds")
S_y.path <- file.path(data.rds.dir, "S_y.rds")

# Load X_y.rds (contains list with X and y)
if (file.exists(X_y.path)) {
    X_y <- readRDS(X_y.path)
    X <- X_y$X
    y <- X_y$y
    message("✓ Loaded X_y.rds: X (", nrow(X), " x ", ncol(X), "), y (", length(y), " elements)")
} else {
    warning("File not found: ", X_y.path, "\n",
            "Make sure data-prep.R has been run to generate this file.")
    X <- NULL
    y <- NULL
}

# Load S_y.rds (contains list with S and y_S)
if (file.exists(S_y.path)) {
    S_y <- readRDS(S_y.path)
    S <- S_y$S
    y_S <- S_y$y
    message("✓ Loaded S_y.rds: S (", nrow(S), " x ", ncol(S), "), y_S (", length(y_S), " elements)")
} else {
    warning("File not found: ", S_y.path, "\n",
            "Make sure data-prep.R has been run to generate this file.")
    S <- NULL
    y_S <- NULL
}

# ============================================================================
# VALIDATION & SUMMARY
# ============================================================================

# Check that at least one dataset loaded successfully
if (is.null(X) && is.null(S)) {
    stop("Failed to load any data files. ",
             "Run data-prep.R first to generate the required .rds files.")
}

# Summary of loaded objects
message("\n========== DATA LOADED ==========")
if (!is.null(X)) {
    message("X: ", nrow(X), " rows x ", ncol(X), " columns (baseline covariate matrix)")
    message("  Columns: ", paste(head(colnames(X), 5), collapse=", "), ", ...")
}
if (!is.null(y)) {
    message("y: ", length(y), " elements (outcome vector)")
}
if (!is.null(S)) {
    message("S: ", nrow(S), " rows x ", ncol(S), " columns (time-varying matrix)")
    message("  Columns: ", paste(head(colnames(S), 5), collapse=", "), ", ...")
}
if (!is.null(y_S)) {
    message("y_S: ", length(y_S), " elements (time-varying outcome)")
}
message("================================\n")

setwd(root)