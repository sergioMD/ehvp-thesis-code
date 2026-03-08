# ============================================================================
# 00_RUN_ALL.R - Master Script: EHVP Full Analysis Pipeline
# ============================================================================
# Paper III of PhD thesis: Effects of the Extended Home Visiting Programme
# on child health outcomes in Sweden.
#
# Description: Runs the complete analysis from matching to final results,
#              including main analysis and dosage sensitivity analyses.
#
# Author:  Sergio Flores
# Date:    2026
# Paper:   Flores et al., "Effects of the Swedish Extended Home Visiting
#          Programme on child health outcomes: A register-based study"
#
# Required packages: here (and all packages loaded by sourced scripts)
#
# USAGE:
#   1. Set RUN_MATCHING = TRUE if you need to redo matching (slow)
#   2. Set which sensitivity analyses to run
#   3. Source this file: source("00_run_all.R")
#
# OUTPUT:
#   Results saved to: output/RESULTS_Main_optimal_YYYY-MM-DD/
#   Sensitivity results saved to: output/RESULTS_Sensitivity70_optimal_YYYY-MM-DD/
# ============================================================================

# Clear environment
rm(list = ls())
gc()

#------------------------------------------------------------------------------
# USER CONFIGURATION
#------------------------------------------------------------------------------

# Run matching? Set FALSE if matching already saved
RUN_MATCHING <- TRUE

# Run main analysis (60% threshold)?
RUN_MAIN_ANALYSIS <- TRUE

# Run sensitivity analyses? (Set to NULL to skip)
# Sensitivity thresholds (from config.R defaults, override here if needed)
SENSITIVITY_THRESHOLDS <- c(0.70, 0.80)  # Or NULL to skip sensitivity analyses

# Verbose logging?
VERBOSE <- TRUE

#------------------------------------------------------------------------------
# SETUP
#------------------------------------------------------------------------------

setwd(here::here())

# Logging function
log_msg <- function(...) {
  if (VERBOSE) {
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
    cat(msg, "\n")
  }
}

log_msg("=== EHVP Analysis Pipeline Started ===")
log_msg("Working directory: ", getwd())

# Check required files exist
required_files <- c(

  "code/config.R",
  "code/helper_functions.R",
  "code/outcomes_config.R",
  "code/load_matching_adapter.R",
  "code/02_analysis.R"
)

if (RUN_MATCHING) {
 required_files <- c(required_files, "code/01_matching.R")
}

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing required files:\n  ", paste(missing, collapse = "\n  "))
}

log_msg("All required files found.")

#------------------------------------------------------------------------------
# STEP 1: MATCHING (Optional)
#------------------------------------------------------------------------------

if (RUN_MATCHING) {
  log_msg("=== STEP 1: Running Matching ===")
  start_time <- Sys.time()

  tryCatch({
    source("code/01_matching.R")
    log_msg("Matching completed successfully.")
    log_msg("Time elapsed: ", round(difftime(Sys.time(), start_time, units = "mins"), 1), " minutes")
  }, error = function(e) {
    stop("Matching failed with error: ", e$message)
  })

} else {
  log_msg("=== STEP 1: Skipping Matching (using saved results) ===")
}

#------------------------------------------------------------------------------
# STEP 2: MAIN ANALYSIS (60% threshold)
#------------------------------------------------------------------------------

if (RUN_MAIN_ANALYSIS) {
  log_msg("=== STEP 2: Running Main Analysis (60% threshold) ===")
  start_time <- Sys.time()

  # Set sensitivity to NULL for main analysis
  SENSITIVITY_DOSAGE <- NULL

  tryCatch({
    source("code/02_analysis.R")

    log_msg("Main analysis completed successfully.")
    log_msg("Time elapsed: ", round(difftime(Sys.time(), start_time, units = "mins"), 1), " minutes")
    log_msg("Results saved to: ", ROOT_OUTPUT_DIR)

  }, error = function(e) {
    warning("Main analysis failed with error: ", e$message)
    print(traceback())
  })

  # Clean up large objects from main analysis
  rm(list = intersect(ls(), c("Stacked_data_individual_level", "SD_bal", "S",
                               "analysis_data_2018_individuals_final",
                               "analysis_data_2019_individuals_final",
                               "analysis_data_2020_individuals_final")))
  gc()
}

#------------------------------------------------------------------------------
# STEP 3: SENSITIVITY ANALYSES (70%, 80% thresholds)
#------------------------------------------------------------------------------

if (!is.null(SENSITIVITY_THRESHOLDS) && length(SENSITIVITY_THRESHOLDS) > 0) {

  for (threshold in SENSITIVITY_THRESHOLDS) {
    pct <- as.integer(threshold * 100)
    log_msg(sprintf("=== STEP 3: Running Sensitivity Analysis (%d%% threshold) ===", pct))
    start_time <- Sys.time()

    # Set sensitivity threshold
    SENSITIVITY_DOSAGE <- threshold

    tryCatch({
      source("code/02_analysis.R")

      log_msg(sprintf("Sensitivity analysis (%d%%) completed successfully.", pct))
      log_msg("Time elapsed: ", round(difftime(Sys.time(), start_time, units = "mins"), 1), " minutes")
      log_msg("Results saved to: ", ROOT_OUTPUT_DIR)

    }, error = function(e) {
      warning(sprintf("Sensitivity analysis (%d%%) failed with error: %s", pct, e$message))
      print(traceback())
    })

    # Clean up large objects between sensitivity runs
    rm(list = intersect(ls(), c("Stacked_data_individual_level", "SD_bal", "S",
                                 "analysis_data_2018_individuals_final",
                                 "analysis_data_2019_individuals_final",
                                 "analysis_data_2020_individuals_final")))
    gc()
  }
}

#------------------------------------------------------------------------------
# SUMMARY
#------------------------------------------------------------------------------

log_msg("=== Pipeline Complete ===")
log_msg("Check the output/ folder for results:")
log_msg("  - RESULTS_Main_optimal_*        : Main analysis (60%)")
if (!is.null(SENSITIVITY_THRESHOLDS)) {
  for (threshold in SENSITIVITY_THRESHOLDS) {
    pct <- as.integer(threshold * 100)
    log_msg(sprintf("  - RESULTS_Sensitivity%d_optimal_* : Sensitivity (%d%%)", pct, pct))
  }
}

log_msg("=== Done! ===")
