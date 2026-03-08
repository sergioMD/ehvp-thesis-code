# ============================================================================
# CONFIG.R - Shared Configuration for EHVP Effects Analysis
# ============================================================================
# Paper III of PhD thesis: Effects of the Extended Home Visiting Programme
# on child health outcomes in Sweden.
#
# Description: Paths, treatment thresholds, matching specifications, and
#              regional settings for the DeSO-level matched analysis.
#
# Author:  Sergio Flores
# Date:    2026
# Paper:   Flores et al., "Effects of the Swedish Extended Home Visiting
#          Programme on child health outcomes: A register-based study"
#
# Required packages: here
#
# OUTCOME DEFINITIONS are in outcomes_config.R (source separately)
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("Loading Configuration...\n")
cat("============================================================================\n\n")

# ============================================================================
# 1. PATHS
# ============================================================================

# Project root directory (uses here::here() for portability)
PROJECT_DIR <- here::here()

# DATA_DIR: Set to location of your register data files
# Example: DATA_DIR <- file.path(PROJECT_DIR, "data")
DATA_DIR <- file.path(PROJECT_DIR, "data")
RAW_PANEL_DATA_PATH <- file.path(DATA_DIR, "final_panel_data_complete.csv")
SPATIAL_DATA_PATH <- file.path(DATA_DIR, "deso_aggregated_dosage_spatial.rds")

# Matching output (stays under output/matching)
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")
MATCHING_OUTPUT_DIR <- file.path(OUTPUT_DIR, "matching")
MATCHING_OBJECTS_PATH <- file.path(MATCHING_OUTPUT_DIR, "optimal", "matching_objects.RData")
ALL_SPECS_PATH <- file.path(MATCHING_OUTPUT_DIR, "all_specs", "all_matching_specs.RData")

# ============================================================================
# 2. TREATMENT THRESHOLDS
# ============================================================================
#
# TREATMENT THRESHOLD RATIONALE:
# 60% chosen as "majority rule" -- DeSO classified as treated if majority of
# births exposed to EHVP. This balances:
#   - Sample size preservation (158 treated DeSOs at 60%)
#   - Treatment intensity (meaningful exposure to programme)
#   - Comparability with control group (clear separation from <10% controls)
#
# Sensitivity analyses at 50%, 70%, 80% bracket this choice to assess
# robustness of findings to threshold selection.
# Approximate sample sizes by threshold:
#   50%: ~200 treated DeSOs (more inclusive)
#   60%: 158 treated DeSOs (main analysis)
#   70%: ~107 treated DeSOs (higher intensity)
#   80%: ~51 treated DeSOs (highest intensity, low power)
#
TREATED_THRESHOLD <- 0.60    # DeSO is "Treated" if max dosage >= 0.60
CONTROL_THRESHOLD <- 0.10    # DeSO is "Untreated" if max dosage <= 0.10 (becomes "Matched Control" after matching)

# Sensitivity analysis thresholds (used by 00_run_all.R)
# These filter treated DeSOs by minimum dosage to test robustness
SENSITIVITY_THRESHOLDS <- c(
  0.50,  # 50% - lower threshold (more inclusive, tests if weaker treatment still shows effects)
  0.70,  # 70% - moderate intensity filter
  0.80   # 80% - high intensity filter (CAUTION: low statistical power)
)

# ============================================================================
# 3. MATCHING SPECIFICATIONS
# ============================================================================

COHORTS_TO_MATCH <- c(2018, 2019, 2020)

MATCHING_COVARIATES <- c(
  "prop_low_education",
  "prop_foreign_born",
  "prop_unemployed",
  "prop_single_parent",
  "prop_low_birth_weight",
  "avg_mother_age_at_birth"
)

# All 6 matching specifications
MATCHING_SPECS <- list(
  spec_1_region_caliper = list(
    name = "Exact Region + Caliper 0.2",
    exact = "lan_code",
    caliper = 0.2,
    ps_includes_region = FALSE,
    description = "Exact match on region with 0.2 caliper"
  ),
  spec_2_region_no_caliper = list(
    name = "Exact Region, No Caliper",
    exact = "lan_code",
    caliper = NULL,
    ps_includes_region = FALSE,
    description = "Exact region without caliper constraint"
  ),
  spec_3_urban_caliper = list(
    name = "Urbanicity + Caliper 0.2",
    exact = "urban_rural_category",
    caliper = 0.2,
    ps_includes_region = FALSE,
    description = "Match within urbanicity with 0.2 caliper"
  ),
  spec_4_urban_no_caliper = list(
    name = "Urbanicity Only, No Caliper [OPTIMAL]",
    exact = "urban_rural_category",
    caliper = NULL,
    ps_includes_region = FALSE,
    description = "OPTIMAL: Best balance with full sample retention"
  ),
  spec_5_soft_region_caliper = list(
    name = "Soft Region (in PS) + Caliper 0.2",
    exact = "urban_rural_category",
    caliper = 0.2,
    ps_includes_region = TRUE,
    description = "Region as PS covariate, not exact match"
  ),
  spec_6_soft_region_no_caliper = list(
    name = "Soft Region (in PS), No Caliper",
    exact = "urban_rural_category",
    caliper = NULL,
    ps_includes_region = TRUE,
    description = "Region as PS covariate without caliper"
  )
)

# --------------------------------------------------------------------------
# ACTIVE MATCHING SPECIFICATION
# --------------------------------------------------------------------------
ACTIVE_MATCHING_SPEC <- "spec_4_urban_no_caliper"  # OPTIMAL

# Extract parameters for backward compatibility
MATCHING_METHOD <- "nearest"
MATCHING_EXACT <- MATCHING_SPECS[[ACTIVE_MATCHING_SPEC]]$exact
MATCHING_CALIPER <- MATCHING_SPECS[[ACTIVE_MATCHING_SPEC]]$caliper
MATCHING_REPLACE <- TRUE
MATCHING_RATIO <- 5

# ============================================================================
# 4. ANALYSIS CONFIGURATION
# ============================================================================

# --- TOGGLE: Dosage Sensitivity Analysis ---
# Controls which treated DeSOs are included in analysis (post-matching filter)
#
# Options:
#   SENSITIVITY_DOSAGE <- NULL   -> Main analysis (all matched treated, 60%+)
#   SENSITIVITY_DOSAGE <- 0.70   -> Sensitivity: only 70%+ dosage DeSOs
#   SENSITIVITY_DOSAGE <- 0.80   -> Sensitivity: only 80%+ dosage DeSOs
#
# NOTE: Matching always uses 60% threshold. This filter is applied AFTER matching.
# NOTE: If SENSITIVITY_DOSAGE is already defined (e.g., by 00_run_all.R), it won't be overwritten.

if (!exists("SENSITIVITY_DOSAGE")) {
  SENSITIVITY_DOSAGE <- NULL   # Default: Main analysis (60%)
}

# Regional analysis configuration
REGIONS_MAP <- list(
  "National"        = NULL,
  "Stockholm"       = "01",
  "Skane"           = "12",
  "Vastra_Gotaland" = "14"
)

# ============================================================================
# 5. MODEL COVARIATES
# ============================================================================

MODEL_COVARIATES <- c(
  "cni_score",
  "foreign_born",
  "single_parent",
  "low_education",
  "unemployed",
  "mother_birth_year",
  "low_birth_weight"
)

# ============================================================================
# 6. VISUALIZATION SETTINGS
# ============================================================================

SMD_EXCELLENT <- 0.10
SMD_ACCEPTABLE <- 0.25

PLOT_COLORS <- list(
  treated = "#D55E00",
  control = "#0072B2",
  fuzzy = "grey85",
  significant = "#E69F00"
)

# ============================================================================
# 7. CREATE MATCHING OUTPUT DIRECTORIES
# ============================================================================

create_matching_dirs <- function() {
  dirs <- c(
    MATCHING_OUTPUT_DIR,
    file.path(MATCHING_OUTPUT_DIR, "optimal"),
    file.path(MATCHING_OUTPUT_DIR, "all_specs"),
    file.path(MATCHING_OUTPUT_DIR, "sensitivity_analysis"),
    file.path(MATCHING_OUTPUT_DIR, "balance_tables"),
    file.path(MATCHING_OUTPUT_DIR, "love_plots"),
    file.path(MATCHING_OUTPUT_DIR, "maps")
  )
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  cat("Matching output directories ready.\n")
}

# ============================================================================
# 8. VALIDATION
# ============================================================================

validate_config <- function() {
  errors <- character(0)

  if (!file.exists(RAW_PANEL_DATA_PATH)) {
    errors <- c(errors, sprintf("Raw panel data not found: %s", RAW_PANEL_DATA_PATH))
  }

  if (!file.exists(SPATIAL_DATA_PATH)) {
    cat(sprintf("  Warning: Spatial data not found: %s\n", SPATIAL_DATA_PATH))
  }

  if (!ACTIVE_MATCHING_SPEC %in% names(MATCHING_SPECS)) {
    errors <- c(errors, sprintf("Invalid ACTIVE_MATCHING_SPEC: %s", ACTIVE_MATCHING_SPEC))
  }

  if (length(errors) > 0) {
    stop("Configuration validation failed:\n", paste("  -", errors, collapse = "\n"))
  }

  cat("Configuration validated successfully.\n")
}

# ============================================================================
# Print summary
# ============================================================================

cat("Configuration loaded:\n")
cat(sprintf("  Project directory: %s\n", PROJECT_DIR))
cat(sprintf("  Cohorts: %s\n", paste(COHORTS_TO_MATCH, collapse = ", ")))
cat(sprintf("  Treatment threshold: %.0f%%\n", TREATED_THRESHOLD * 100))
cat(sprintf("  Active matching spec: %s\n", ACTIVE_MATCHING_SPEC))
if (is.null(SENSITIVITY_DOSAGE)) {
  cat("  Analysis mode: Main (60%)\n")
} else {
  cat(sprintf("  Analysis mode: Sensitivity (%.0f%%)\n", SENSITIVITY_DOSAGE * 100))
}
cat("\n")
cat("\n")
