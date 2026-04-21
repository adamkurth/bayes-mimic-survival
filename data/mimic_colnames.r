# ==============================================================================
# MIMIC-III Column Names Reference (Fast Header-Only)
# ==============================================================================
# Purpose: Extract and catalog all column names from MIMIC-III CSV files
# for proper R variable management during analysis
#
# Usage:
#   source("mimic_colnames.r")
#
# Output:
#   - Console: formatted list of all column names by table
#   - mimic_columns_ref.R: R code with column vectors for each table
# ==============================================================================

library(data.table)

# Set working directory to MIMIC data location
mimic_dir <- "/Users/adamkurth/Documents/RStudio/bayes-mimic-survival/data/physionet.org/files/mimiciii/1.4"

# Get all CSV files
csv_files <- list.files(mimic_dir, pattern = "\\.csv$", full.names = FALSE)
csv_files <- sort(csv_files)

cat("\n")
cat("================================================================================\n")
cat("MIMIC-III Column Names Reference\n")
cat("================================================================================\n")
cat(sprintf("Directory: %s\n", mimic_dir))
cat(sprintf("Total CSV files: %d\n\n", length(csv_files)))

# Create list to store all column names
mimic_columns <- list()

# Extract column names from each CSV (header only, fast)
for (csv_file in csv_files) {
  filepath <- file.path(mimic_dir, csv_file)
  
  # Use system command to get header only (much faster than reading full file)
  header_raw <- system(paste("head -1", shQuote(filepath)), intern = TRUE)
  
  # Parse header manually
  colnames_vec <- unlist(strsplit(header_raw, ","))
  colnames_vec <- gsub('^[[:space:]]*"|"[[:space:]]*$', '', colnames_vec)
  
  # Store in list
  table_name <- gsub("\\.csv$", "", csv_file)
  mimic_columns[[table_name]] <- colnames_vec
  
  # Print formatted output
  cat(sprintf("\n%-30s (%d columns)\n", paste0("--- ", csv_file, " ---"), length(colnames_vec)))
  cat(sprintf("%-35s | %s\n", "Column Name", "Index"))
  cat(strrep("-", 80), "\n")
  
  for (i in seq_along(colnames_vec)) {
    cat(sprintf("%-35s | %2d\n", colnames_vec[i], i))
  }
}

cat("\n")
cat("================================================================================\n")
cat("Summary Statistics\n")
cat("================================================================================\n")
total_cols <- sum(sapply(mimic_columns, length))
cat(sprintf("Total columns across all tables: %d\n", total_cols))
cat(sprintf("Average columns per table: %.1f\n", mean(sapply(mimic_columns, length))))
cat(sprintf("Min columns: %d (%s)\n", 
            min(sapply(mimic_columns, length)),
            names(mimic_columns)[which.min(sapply(mimic_columns, length))]))
cat(sprintf("Max columns: %d (%s)\n\n", 
            max(sapply(mimic_columns, length)),
            names(mimic_columns)[which.max(sapply(mimic_columns, length))]))

# ============================================================================
# Generate R code for column references (easy copy-paste)
# ============================================================================
cat("================================================================================\n")
cat("R Code for Column References (copy-paste ready)\n")
cat("================================================================================\n")
cat("\n# Column lists for each table:\n\n")

r_code <- "# MIMIC-III Column References\n# Auto-generated column vectors for use in R\n\n"

for (table_name in names(mimic_columns)) {
  cols <- mimic_columns[[table_name]]
  col_str <- paste(sprintf('"%s"', cols), collapse = ", ")
  code_line <- sprintf("%s_cols <- c(%s)\n", tolower(table_name), col_str)
  cat(code_line)
  r_code <- paste0(r_code, code_line)
}

# Save R code to file
r_file <- file.path(mimic_dir, "../../../mimic_columns_ref.R")
writeLines(r_code, r_file)
cat("\n✓ Saved R reference code to: mimic_columns_ref.R\n")

cat("\n")
cat("================================================================================\n")
cat("Usage in Your Analysis\n")
cat("================================================================================\n")
cat("\nTo access in future R sessions:\n")
cat("  source('mimic_columns_ref.R')\n")
cat("  # Then use: admissions_cols, patients_cols, icustays_cols, etc.\n\n")
cat("Example data.table selection:\n")
cat("  admissions <- fread('ADMISSIONS.csv')[, ..admissions_cols]\n")
cat("  icustays <- fread('ICUSTAYS.csv')[, ..icustays_cols]\n\n")





