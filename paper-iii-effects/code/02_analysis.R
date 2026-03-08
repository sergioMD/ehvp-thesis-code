#------------------------------------------------------------------------------#
# EHVP Effects Analysis — Individual-Level Event Study
#------------------------------------------------------------------------------#
#
# Description:
#   Main analysis script for Paper III (Effects of the Extended Home Visiting
#   Programme). Implements staggered difference-in-differences with multiple
#   estimators (stacked TWFE, Sun-Abraham, Callaway-Sant'Anna) on individual-
#   level register data. Includes heterogeneity analysis by neighbourhood
#   deprivation (CNI quartiles), formal pre-trends testing, and table/figure
#   generation.
#
# Paper:
#   Flores S. Effects of the Extended Home Visiting Programme on child health
#   outcomes: a quasi-experimental study using Swedish register data.
#
# Author: Sergio Flores
# Dependencies: data.table, fixest, ggplot2, patchwork, MatchIt, cobalt,
#               modelsummary, did, ggtext, openxlsx, kableExtra, dplyr
#
# Pipeline order: config.R -> 01_matching.R -> 02_analysis.R
# Assumes: config.R provides paths, thresholds, and matching spec selection
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# General setup and data preparation. Run sections 1-6 once
#------------------------------------------------------------------------------#

# ==== 1: Working directory, packages, libraries and functions ====
# SETUP - import data, install and/or load packages, load helper functions

# Load configuration (paths, thresholds, matching spec selection)
source("code/config.R")

# Load matching objects (uses ACTIVE_MATCHING_SPEC from config)
source("code/load_matching_adapter.R")


# Installing Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

library(pacman)

pacman::p_load(
  data.table, ggplot2, fixest, broom, stringr, MatchIt, cobalt,
  modelsummary, sf, patchwork, kableExtra, forcats, dplyr, did, ggtext, openxlsx,
  character.only = FALSE
)

#Helper functions
source("code/helper_functions.R")

# Outcome definitions (single source of truth)
source("code/outcomes_config.R")


cat("\n--- Loading the complete raw panel data... (This may take a moment) ---\n")

# 1. Load the very large original dataset
full_raw_data <- fread(RAW_PANEL_DATA_PATH)



# ============================================================================
# 2. DERIVED VARIABLE CREATION
# ============================================================================
cat("--- Deriving Outcome Variables --- \n")

# A. EMERGENCY VISITS (Inpatient + Outpatient) - PRIMARY OUTCOME
full_raw_data[, count_emergency :=
                fifelse(is.na(inpatient_emergency_events), 0L, inpatient_emergency_events) +
                fifelse(is.na(outpatient_emergency_events), 0L, outpatient_emergency_events)]
full_raw_data[, count_planned :=
                (fifelse(is.na(inpatient_total_events), 0L, inpatient_total_events) +
                 fifelse(is.na(outpatient_total_events), 0L, outpatient_total_events)) - count_emergency]
full_raw_data[, total_visits :=
                fifelse(is.na(inpatient_total_events), 0L, inpatient_total_events) +
                fifelse(is.na(outpatient_total_events), 0L, outpatient_total_events)]

# B. PLANNED CARE SHARE (Exploratory)
full_raw_data[, planned_care_share := fifelse(total_visits > 0, count_planned / total_visits, NA_real_)]
full_raw_data[, emergency_share_ratio := fifelse(total_visits > 0, count_emergency / total_visits, NA_real_)]

# C. SENTINEL INJURIES (Burns + Poisonings) - PRIMARY OUTCOME
full_raw_data[, sentinel_injury_events :=
                fifelse(is.na(inpatient_trauma_burns_events), 0L, inpatient_trauma_burns_events) +
                fifelse(is.na(inpatient_trauma_poisoning_icd_events), 0L, inpatient_trauma_poisoning_icd_events)]

# D. COMBINED ASTHMA EVENTS (Inpatient + Outpatient) - PRIMARY OUTCOME
full_raw_data[, asthma_combined_events :=
                fifelse(is.na(inpatient_asthma_events), 0L, inpatient_asthma_events) +
                fifelse(is.na(outpatient_asthma_events), 0L, outpatient_asthma_events)]

# E. GENDER EQUITY IN PARENTAL LEAVE - SECONDARY OUTCOME
# High NA rate expected: parental leave typically taken in year 0-1 only
full_raw_data[, fp_father_share := fifelse(
  (mother_fp_net_days_year + father_fp_net_days_year) > 0,
  father_fp_net_days_year / (mother_fp_net_days_year + father_fp_net_days_year),
  NA_real_
)]
full_raw_data[, fp_father_any := fifelse(father_fp_net_days_year > 0, 1L, 0L)]

# F. Injury alias (backward compatibility)
if ("inpatient_trauma_any_diag_events" %in% names(full_raw_data)) {
  full_raw_data[, inpatient_any_injury_events := inpatient_trauma_any_diag_events]
} else {
  full_raw_data[, inpatient_any_injury_events := 0L]
}

# G. Extensive Margin Flags
full_raw_data[, had_any_inpatient_year := fifelse(inpatient_total_events > 0, 1L, 0L)]
full_raw_data[, had_avoidable_event_year := fifelse((inpatient_avoidable_events + outpatient_avoidable_events) > 0, 1L, 0L)]

cat("--- All Derived Variables Created Successfully --- \n")
cat("   PRIMARY: sentinel_injury_events, count_emergency, asthma_combined_events\n")
cat("   DERIVED (kept for data): fp_father_share, planned_care_share (not in main analysis)\n")


# 3. Set key for speed
setkey(full_raw_data, deso_2021)

# --- Filter for Cohort 2018 ---
cat("--- Filtering raw data for 2018 matched sample ---\n")
matched_data_2018 <- match.data(match_results_list_stratified$`2018`$match_object)
matched_desos_2018 <- as.character(unique(matched_data_2018$deso_id))
analysis_data_2018_individuals <- full_raw_data[.(matched_desos_2018), nomatch = 0L]

# --- Filter for Cohort 2019 ---
cat("--- Filtering raw data for 2019 matched sample ---\n")
matched_data_2019 <- match.data(match_results_list_stratified$`2019`$match_object)
matched_desos_2019 <- as.character(unique(matched_data_2019$deso_id))
analysis_data_2019_individuals <- full_raw_data[.(matched_desos_2019), nomatch = 0L]

# --- Filter for Cohort 2020 ---
cat("--- Filtering raw data for 2020 matched sample ---\n")
matched_data_2020 <- match.data(match_results_list_stratified$`2020`$match_object)
matched_desos_2020 <- as.character(unique(matched_data_2020$deso_id))
analysis_data_2020_individuals <- full_raw_data[.(matched_desos_2020), nomatch = 0L]

# --- Clean up to save memory ---
cat("--- Cleaning up large data objects from memory ---\n")
rm(full_raw_data, matched_data_2018, matched_data_2019, matched_data_2020)
gc()


#------------------------------------------------------------------------------#
# NOTE: The analysis data was retrieved from the original individual level
# data. Treated and control units were extracted from the original individual
# level data set. This data set was very large so this was done once and the
# subsets saved directly.
# The analysis subsets are retrieved by load("analysis_data_individual_level.RData")
# and were created by the following code (example for cohort 2018):
#
#   # Create the following objects (one for every cohort).
#     matched_data_2018 <- match.data(match_results_list_stratified$`2018`$match_object)
#     setDT(matched_data_2018)
#     treated_desos_for_2018 <- matched_data_2018[treated_binary == 1, deso_id]
#     control_desos_for_2018 <- matched_data_2018[treated_binary == 0, deso_id]
#     matched_sample_2018_DeSO_id <- c(treated_desos_for_2018, control_desos_for_2018)
#
#   # Raw individual level file (very large, load only if necessary)
#     final_panel_data_individual_level <- fread("final_panel_data_complete_with_FK_MBR_Meds_20250720.csv")
#     setDT(final_panel_data_individual_level)
#     ids <- unique(as.character(matched_sample_2018_DeSO_id))
#
#   # optional: set key for speed
#     setkey(final_panel_data_individual_level, deso_2021)
#
#     analysis_data_2018_individuals <- final_panel_data_individual_level[.(ids), nomatch = 0L]
#------------------------------------------------------------------------------#

# ==== 2: Extracting matched data AND WEIGHTS ====
# When matching with replacement, MatchIt creates a 'weights' column.
# These weights must be preserved and used in outcome regressions.
# Reference: Ho et al. (2007), Stuart (2010) - weights account for control reuse

cat("\n--- Extracting matched data with weights ---\n")

# Extract the matched DeSOs WITH weights
matched_data_2018 <- match.data(match_results_list_stratified$`2018`$match_object)
matched_data_2019 <- match.data(match_results_list_stratified$`2019`$match_object)
matched_data_2020 <- match.data(match_results_list_stratified$`2020`$match_object)

setDT(matched_data_2018)
setDT(matched_data_2019)
setDT(matched_data_2020)

# Extract treated DeSOs (for treatment assignment)
treated_desos_for_2018 <- matched_data_2018[treated_binary == 1, deso_id]
treated_desos_for_2019 <- matched_data_2019[treated_binary == 1, deso_id]
treated_desos_for_2020 <- matched_data_2020[treated_binary == 1, deso_id]

# =============================================================================
# CREATE DeSO-LEVEL WEIGHTS TABLE (critical for matching with replacement)
# The 'weights' column from match.data() is at the DeSO level.
# Treated units always get weight=1; controls get weight based on reuse count.
# =============================================================================
deso_weights_2018 <- matched_data_2018[, .(deso_id, match_weight = weights)][, cohort := "2018"]
deso_weights_2019 <- matched_data_2019[, .(deso_id, match_weight = weights)][, cohort := "2019"]
deso_weights_2020 <- matched_data_2020[, .(deso_id, match_weight = weights)][, cohort := "2020"]

# Combine into master weights table
DESO_WEIGHTS <- rbindlist(list(deso_weights_2018, deso_weights_2019, deso_weights_2020))
setnames(DESO_WEIGHTS, "deso_id", "deso_2021")  # Match column name in analysis data

cat(sprintf("  Created weights table: %d DeSO-cohort combinations\n", nrow(DESO_WEIGHTS)))
cat(sprintf("  Weight range: [%.2f, %.2f]\n", min(DESO_WEIGHTS$match_weight), max(DESO_WEIGHTS$match_weight)))
cat(sprintf("  Control DeSOs with weight > 1 (used multiple times): %d\n",
            sum(DESO_WEIGHTS$match_weight > 1)))

# Cleaning (keep DESO_WEIGHTS for later use)
rm(matched_data_2018, matched_data_2019, matched_data_2020,
   deso_weights_2018, deso_weights_2019, deso_weights_2020)




# ============================================================================
#  GLOBAL FILTER FOR FIRST-BORNS
# ============================================================================
cat("\n--- APPLYING GLOBAL FILTER: FIRST-BORN CHILDREN ONLY ---\n")

# Check counts before
n_before <- nrow(analysis_data_2018_individuals) + nrow(analysis_data_2019_individuals) + nrow(analysis_data_2020_individuals)

# Apply Filter (child_ordnrmor == 1)
analysis_data_2018_individuals <- analysis_data_2018_individuals[child_ordnrmor == 1]
analysis_data_2019_individuals <- analysis_data_2019_individuals[child_ordnrmor == 1]
analysis_data_2020_individuals <- analysis_data_2020_individuals[child_ordnrmor == 1]

n_after <- nrow(analysis_data_2018_individuals) + nrow(analysis_data_2019_individuals) + nrow(analysis_data_2020_individuals)

cat(sprintf("--- Data filtered to First-Borns. Rows reduced from %s to %s ---\n",
            format(n_before, big.mark=","), format(n_after, big.mark=",")))


# ============================================================================
# CONFIGURATION & DYNAMIC FOLDER SETUP
# ============================================================================

# --- 1. ANALYSIS IDENTIFIER ---
# Include matching spec in folder name (so different specs don't overwrite)
spec_label <- if (ACTIVE_MATCHING_SPEC == "spec_4_urban_no_caliper") "optimal" else ACTIVE_MATCHING_SPEC

# --- 2. PRE-REQ: Create 'analysis_group' Variable ---
# (Must be done before filtering)
analysis_data_2018_individuals[, analysis_group := fifelse(deso_2021 %in% treated_desos_for_2018, "Treated", "Control")]
analysis_data_2019_individuals[, analysis_group := fifelse(deso_2021 %in% treated_desos_for_2019, "Treated", "Control")]
analysis_data_2020_individuals[, analysis_group := fifelse(deso_2021 %in% treated_desos_for_2020, "Treated", "Control")]

# --- 3. DEFINE DATED OUTPUT DIRECTORY ---
today_str <- format(Sys.Date(), "%Y-%m-%d")

if (!is.null(SENSITIVITY_DOSAGE)) {
  # Sensitivity Analysis (70% or 80%)
  dosage_pct <- as.integer(SENSITIVITY_DOSAGE * 100)
  folder_name <- sprintf("RESULTS_Sensitivity%d_%s_%s", dosage_pct, spec_label, today_str)
  ROOT_OUTPUT_DIR <- file.path(getwd(), "output", folder_name)

  cat(sprintf("\n!!! CONFIGURATION: SENSITIVITY ANALYSIS (%d%% DOSAGE) !!!\n", dosage_pct))
  cat(sprintf("!!! MATCHING SPEC: %s !!!\n", MATCHING_SPECS[[ACTIVE_MATCHING_SPEC]]$name))
  cat(sprintf("--- Outputs will be saved to: output/%s ---\n", folder_name))

  # Logic to filter treated units
  apply_dosage_filter <- function(dt, threshold) {
    return(dt[analysis_group == "Control" | (analysis_group == "Treated" & expected_deso_treatment_dosage >= threshold)])
  }

  # Apply Filter
  analysis_data_2018_individuals <- apply_dosage_filter(analysis_data_2018_individuals, SENSITIVITY_DOSAGE)
  analysis_data_2019_individuals <- apply_dosage_filter(analysis_data_2019_individuals, SENSITIVITY_DOSAGE)
  analysis_data_2020_individuals <- apply_dosage_filter(analysis_data_2020_individuals, SENSITIVITY_DOSAGE)

  cat(sprintf("--- Weakly treated units (<%d%%) removed from memory. ---\n", dosage_pct))

} else {
  # Main Analysis (60% - all matched treated units)
  folder_name <- paste0("RESULTS_Main_", spec_label, "_", today_str)
  ROOT_OUTPUT_DIR <- file.path(getwd(), "output", folder_name)

  cat(sprintf("\n!!! CONFIGURATION: MAIN ANALYSIS (STANDARD 60%%) !!!\n"))
  cat(sprintf("!!! MATCHING SPEC: %s !!!\n", MATCHING_SPECS[[ACTIVE_MATCHING_SPEC]]$name))
  cat(sprintf("--- Outputs will be saved to: output/%s ---\n", folder_name))
  cat("--- Full matched sample used. ---\n")
}

# --- 4. CREATE DIRECTORY STRUCTURE IMMEDIATELY ---
if (!dir.exists(ROOT_OUTPUT_DIR)) dir.create(ROOT_OUTPUT_DIR, recursive = TRUE)

# Create subfolders now so they are ready for the loop and save functions
dir.create(file.path(ROOT_OUTPUT_DIR, "plots"), showWarnings = FALSE)
dir.create(file.path(ROOT_OUTPUT_DIR, "tables"), showWarnings = FALSE)
dir.create(file.path(ROOT_OUTPUT_DIR, "plots", "heterogeneity"), showWarnings = FALSE)
dir.create(file.path(ROOT_OUTPUT_DIR, "plots", "sanity_check"), showWarnings = FALSE)

# ==== 3: Dropping redundant variables from the analysis data  ====

# analysis_group already created above in Configuration section

analysis_data_2018_individuals[, first_treat_year := fifelse(analysis_group == "Treated", 2018, Inf) ]
analysis_data_2019_individuals[, first_treat_year := fifelse(analysis_group == "Treated", 2019, Inf) ]
analysis_data_2020_individuals[, first_treat_year := fifelse(analysis_group == "Treated", 2020, Inf) ]

# Variables to keep (defined in outcomes_config.R)
k <- KEEP_VARS_FLAT
analysis_data_2018_individuals <- analysis_data_2018_individuals[,..k]
analysis_data_2019_individuals <- analysis_data_2019_individuals[,..k]
analysis_data_2020_individuals <- analysis_data_2020_individuals[,..k]

# Cleaning
rm(k)




# ==== 4: Structuring the analysis data for event study analysis ====



# Discard all observations where year is less than child_birth year
analysis_data_2018_individuals <- analysis_data_2018_individuals[year >= child_birth_year]
analysis_data_2019_individuals <- analysis_data_2019_individuals[year >= child_birth_year]
analysis_data_2020_individuals <- analysis_data_2020_individuals[year >= child_birth_year]


# Dummy variable for treated deso
analysis_data_2018_individuals[, deso_treated_binary := fifelse(deso_2021 %in% treated_desos_for_2018, 1, 0)]
analysis_data_2019_individuals[, deso_treated_binary := fifelse(deso_2021 %in% treated_desos_for_2019, 1, 0)]
analysis_data_2020_individuals[, deso_treated_binary := fifelse(deso_2021 %in% treated_desos_for_2020, 1, 0)]


# Event time is based on when the child was born
analysis_data_2018_individuals[,event_time :=  (child_birth_year - 2018)]
analysis_data_2019_individuals[,event_time :=  (child_birth_year - 2019)]
analysis_data_2020_individuals[,event_time :=  (child_birth_year - 2020)]

#-----------------------------------------------------------------------------#
# NOTE: This operation removes data from 2022 since there are no children
# born in 2022 in any of the data sets. Verify by running:
#
#       max(analysis_data_2018_individuals$child_birth_year)
#       max(analysis_data_2018_individuals$year)
#
#       max(analysis_data_2019_individuals$child_birth_year)
#       max(analysis_data_2018_individuals$year)
#
#       max(analysis_data_2018_individuals$year)
#       max(analysis_data_2020_individuals$child_birth_year)
#
#       prior to this section.
#------------------------------------------------------------------------------#

# ==== 5: Mapping outcome variables and co-variates to individual-event-time pair ====

#------------------------------------------------------------------------------#
# Collapse outcome variables onto one row per individual in the data set.
#
# A child is treated when born (the intervention targeted newborns) IF the
# birth date >= programme start year AND the child was born in a programme-
# exposed DeSO. A child born in 2017 in a DeSO that started the programme
# in 2018 is untreated (event_time = -1) but belongs to the treated DeSO.
# Pre-treatment, this child should be compared to a child born in 2017 in
# a non-exposed DeSO.
#
# A child born in 2018 in the exposed DeSO is treated at event_time = 0;
# born in 2019, event_time = +1; and so on. These individuals are compared
# to those born in non-exposed DeSOs at the corresponding event-times.
#
# Outcomes are collapsed so that outcomes in the birth year map to
# event_time = 0, outcomes one year after birth map to event_time = 1, etc.
# This creates a child-age dimension: every child has outcomes observed from
# birth onwards, enabling examination of effects at ages 0, 1, 2, 3.
#
# A child born in a treated cohort in 2018 has immediate effects that year
# (event_time = 0) and also effects in 2019, 2020, etc. (event_times 1, 2, ...).
# Keeping one row per year from birth onward would mix outcomes at
# event_time = 1 for children born in 2018 (one-year-olds) with newborns
# in 2019, confounding the comparisons.
#
# Retaining one row per year would yield something closer to a cumulative
# effect estimate rather than an age-specific effect.
#------------------------------------------------------------------------------#

# Use master list of outcome variables from outcomes_config.R
all_outcome_vars <- ALL_OUTCOME_VARS

### 2018 Cohort
analysis_data_2018_individuals_final <- collapse_outcomes_reltime(
  analysis_data_2018_individuals,
  outcomes = all_outcome_vars,
  horizon = 0:4,
  id = "lopnr",
  year = "year",
  anchor = "child_birth_year"
)

### 2019 Cohort
analysis_data_2019_individuals_final <- collapse_outcomes_reltime(
  analysis_data_2019_individuals,
  outcomes = all_outcome_vars,
  horizon = 0:4,
  id = "lopnr",
  year = "year",
  anchor = "child_birth_year"
)

### 2020 Cohort
analysis_data_2020_individuals_final <- collapse_outcomes_reltime(
  analysis_data_2020_individuals,
  outcomes = all_outcome_vars,
  horizon = 0:4,
  id = "lopnr",
  year = "year",
  anchor = "child_birth_year"
)

# Base outcomes for summary calculation
base_outcomes <- all_outcome_vars

# Summing outcome variables year 0-1 and 0-3
add_rel_summaries(analysis_data_2018_individuals_final, bases = base_outcomes)
add_rel_summaries(analysis_data_2019_individuals_final, bases = base_outcomes)
add_rel_summaries(analysis_data_2020_individuals_final, bases = base_outcomes)

# Cleaning
rm(analysis_data_2018_individuals,
   analysis_data_2019_individuals,
   analysis_data_2020_individuals)

#------------------------------------------------------------------------------#
# NOTE: What is collapsed should be all outcome variables and covariates that
# change over time. In the current setup, all covariates are fixed in time
# so they do not need to be collapsed.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# NOTE: The maximum event time is now 3. No children in the data set were born
# in 2022. Verify by running:
#
#     max(analysis_data_2018_individuals_final$event_time[analysis_data_2018_individuals_final$event_time != Inf])
#     max(analysis_data_2019_individuals_final$event_time[analysis_data_2019_individuals_final$event_time != Inf])
#     max(analysis_data_2020_individuals_final$event_time[analysis_data_2020_individuals_final$event_time != Inf])
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# NOTE:
# When summing outcome variables over the first four years, some cohorts lose
# observations in the right tail of the event-time distribution.
#
# Reason: Observations for the four-year sum exist only for children born
# between 2012 and 2019, since the data ends in 2022 (cannot sum over four
# years with only three available). This corresponds to event times -6 to 1
# for the 2018 cohort, -7 to 0 for 2019, and -8 to -1 for 2020.
# Hollow points in plots indicate incomplete-sample event-times.
#
# Verification code (un-comment and run to check):
#
# IsData_0_3 <- make_cohort_intersection(
#   dts = list(
#     `2018` = analysis_data_2018_individuals_final,
#     `2019` = analysis_data_2019_individuals_final,
#     `2020` = analysis_data_2020_individuals_final
#   ),
#   sum_suffix = "_rel_time0_3_sum$"
# )
#
# IsData_0_1 <- make_cohort_intersection(
#   dts = list(
#     `2018` = analysis_data_2018_individuals_final,
#     `2019` = analysis_data_2019_individuals_final,
#     `2020` = analysis_data_2020_individuals_final
#   ),
#   sum_suffix = "_rel_time01_sum$"
# )
#
#------------------------------------------------------------------------------#

# ==== 6: Creating the stacked data, the CSA data the weights for aggregating models ====

# Creating stacked data
setDT(analysis_data_2018_individuals_final)
setDT(analysis_data_2019_individuals_final)
setDT(analysis_data_2020_individuals_final)

analysis_data_2018_individuals_final[, cohort := "2018"]
analysis_data_2019_individuals_final[, cohort := "2019"]
analysis_data_2020_individuals_final[, cohort := "2020"]

analysis_data_2018_individuals_final[, stack_deso_id := paste0(cohort, ":", deso_2021)]
analysis_data_2019_individuals_final[, stack_deso_id := paste0(cohort, ":", deso_2021)]
analysis_data_2020_individuals_final[, stack_deso_id := paste0(cohort, ":", deso_2021)]

analysis_data_2018_individuals_final[, coh_time := interaction(cohort, year, drop = TRUE)]
analysis_data_2019_individuals_final[, coh_time := interaction(cohort, year, drop = TRUE)]
analysis_data_2020_individuals_final[, coh_time := interaction(cohort, year, drop = TRUE)]

Stacked_data_individual_level <- rbindlist(
  list(
    analysis_data_2018_individuals_final,
    analysis_data_2019_individuals_final,
    analysis_data_2020_individuals_final
  ),
  use.names = TRUE, fill = TRUE
)

# Create clean factors
setDT(Stacked_data_individual_level)
Stacked_data_individual_level[, cohort_factor := as.factor(cohort)]
Stacked_data_individual_level[, trt    := factor(deso_treated_binary, levels = c(0, 1), labels = c("ctrl","trt"))]

# Build a combined factor: cohort x treated
Stacked_data_individual_level[, coh_trt := interaction(cohort_factor, trt, drop = TRUE)]

#  Keep only the cohort x treated levels (end with ".trt")
keep_treated_lvls <- levels(Stacked_data_individual_level$coh_trt)
keep_treated_lvls <- keep_treated_lvls[endsWith(keep_treated_lvls, ".trt")]

# =============================================================================
# DYNAMIC COHORT HANDLING
# Calculate N_treated dynamically from data instead of hardcoding.
# This allows easy exclusion of cohorts (e.g., 2020 for COVID robustness).
# =============================================================================

# Calculate N_treated per cohort DYNAMICALLY from the actual data
N_by_cohort <- Stacked_data_individual_level[
  is.finite(first_treat_year),
  .(N_treated = uniqueN(deso_2021[deso_treated_binary == 1])),
  by = cohort
][order(cohort)]

cat(sprintf("\n--- Cohort Sizes (N Treated DeSOs) ---\n"))
print(N_by_cohort)

# Derive active cohorts from config (allows exclusion via COHORTS_TO_MATCH)
active_cohorts <- as.character(COHORTS_TO_MATCH)
n_active_cohorts <- length(active_cohorts)

# =============================================================================
# MERGE MATCHING WEIGHTS INTO STACKED DATA
# Each observation needs the weight from its cohort-specific matching
# =============================================================================
cat("\n--- Merging matching weights into analysis data ---\n")

# Merge weights by DeSO and cohort
Stacked_data_individual_level <- merge(
  Stacked_data_individual_level,
  DESO_WEIGHTS,
  by = c("deso_2021", "cohort"),
  all.x = TRUE
)

# Verify merge success
n_missing_weights <- sum(is.na(Stacked_data_individual_level$match_weight))
if (n_missing_weights > 0) {
  warning(sprintf("WARNING: %d observations have missing weights! Setting to 1.", n_missing_weights))
  # Set missing weights to 1 (conservative assumption - treats as unweighted)
  Stacked_data_individual_level[is.na(match_weight), match_weight := 1]
} else {
  cat("  All observations successfully merged with weights.\n")
}

cat(sprintf("  Mean weight: %.3f | SD: %.3f\n",
            mean(Stacked_data_individual_level$match_weight),
            sd(Stacked_data_individual_level$match_weight)))
cat(sprintf("  Observations with weight > 1: %d (%.1f%%)\n",
            sum(Stacked_data_individual_level$match_weight > 1),
            100 * mean(Stacked_data_individual_level$match_weight > 1)))

# Creating balanced stacked panel
# DYNAMIC: requires all active cohorts to be present
SD_bal <- Stacked_data_individual_level[
  , if (uniqueN(cohort, na.rm = TRUE) == n_active_cohorts) .SD, by = event_time
]

# Creating data for Callaway Sant'Anna
# Copying stacked data to not interfere with the other uses of the data
S <- copy(Stacked_data_individual_level)

# Calculate number of treated in each cohort
N_treated <- S[is.finite(first_treat_year),
               .(N = uniqueN(deso_2021[deso_treated_binary == 1])),
               by = cohort][order(cohort)]
N_treated[,cohort := as.integer(cohort)]
# Ensure numeric time and group (g)
S[, year := as.integer(year)]

# Cohorts for CSA loop - derived from config (not hardcoded)
cohorts <- COHORTS_TO_MATCH


#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------#
# Section 6: Regional Analysis Loop (WITH SMART DATA FILTER)
#------------------------------------------------------------------------------#

# Regional and outcome configuration (from config files)
regions_map <- REGIONS_MAP
base_outcomes <- BASE_OUTCOMES
pretty_labels <- PRETTY_LABELS
covars <- MODEL_COVARIATES

# Formulas
treat_twfe  <- "i(event_time, deso_treated_binary, ref = -1)"
treat_sunab <- "sunab(first_treat_year, year, ref.p = -1)"
treat_cohtrt <- "i(event_time, coh_trt, ref = -1, keep2 = keep_treated_lvls)"
fe_disagg  <- "deso_2021 + year"
fe_stacked <- "stack_deso_id + coh_time"

# --- START REGION LOOP ---
# Create a master list to store results for ALL regions
MASTER_REGIONAL_RESULTS <- list()
for (reg_name in names(regions_map)) {

  reg_code <- regions_map[[reg_name]]
  cat(sprintf("\n### STARTING ANALYSIS FOR REGION: %s ###\n", reg_name))

  # Function to subset data by region code
  filter_region <- function(dt, code) {
    if (is.null(code)) return(copy(dt))
    return(dt[substr(deso_2021, 1, 2) == code])
  }

  # Create region-specific datasets (Base versions)
  DT_2018_Reg <- filter_region(analysis_data_2018_individuals_final, reg_code)
  DT_2019_Reg <- filter_region(analysis_data_2019_individuals_final, reg_code)
  DT_2020_Reg <- filter_region(analysis_data_2020_individuals_final, reg_code)
  DT_Stacked_Reg <- filter_region(Stacked_data_individual_level, reg_code)
  DT_Bal_Reg    <- filter_region(SD_bal, reg_code)
  DT_S_Reg      <- filter_region(S, reg_code)

  if (nrow(DT_Stacked_Reg) == 0) { cat(sprintf("Skipping %s: No data.\n", reg_name)); next }

  # Recalculate Weights for this Region (Base) - DYNAMIC cohorts
  N_treated_reg <- DT_S_Reg[is.finite(first_treat_year), .(N = uniqueN(deso_2021[deso_treated_binary == 1])), by = cohort][order(cohort)]
  N_treated_reg[, cohort := as.character(cohort)]
  # Initialize with active cohorts (from config), not hardcoded
  N_by_cohort_reg <- data.table(cohort = active_cohorts, N_treated = 0)
  if(nrow(N_treated_reg) > 0) N_by_cohort_reg[match(N_treated_reg$cohort, cohort), N_treated := N_treated_reg$N]

  # Initialize Containers
  MODELS  <- list(); RESULTS <- list()
  if (!exists("CSA_ERROR_LOG")) CSA_ERROR_LOG <- data.table()

  # --- Run Analysis Loop ---
  for (use_covariates in c(TRUE, FALSE)) {
    covars_on <- if (use_covariates) "covars_on" else "covars_off"
    note_base <- sprintf("%s Analysis (First-Borns): %s covariates", reg_name, if (use_covariates) "with" else "NO")
    covar_str  <- if (use_covariates) paste(covars, collapse = " + ") else NULL

    # --- OUTCOME LOOP WITH SMART FILTER ---
    for (outcome_idx in seq_along(base_outcomes)) {
      y_base <- base_outcomes[outcome_idx]
      label_base <- pretty_labels[outcome_idx]

      # --- BIRTH-TIME-ONLY CHECK ---
      # Outcomes like fp_father_share are only meaningful at birth;
      # restrict analysis to horizon_k = 0 only
      if (is_birth_time_only(y_base)) {
        horizon_range <- 0
        if(use_covariates) cat(sprintf("   [Birth-Time-Only] %s: Analyzing at age 0 only.\n", label_base))
      } else {
        horizon_range <- 0:3
      }

      # --- SMART FILTER: SELECT ANALYTIC SAMPLE ---
      # Apply quality filter to outcomes with incomplete BVC/registry coverage
      # (defined in SMART_FILTER_OUTCOMES in outcomes_config.R)
      is_smart_filter_outcome <- needs_smart_filter(y_base)

      if (is_smart_filter_outcome) {

        # TIER 2: Filter to Valid Reporters

        # Identify Valid DeSOs (mean > 0 or !NA)
        # Calculated on the Stacked Region data to define the universe.
        #
        # For age-dependent outcomes like dtap_3dos_utd_by_18m, the raw column
        # at anchor year (age 0) is 0/NA because measurement happens later.
        # Check _rel_time columns where data actually exists.
        # For dtap (measured at ~18 months), data appears at rel_time2 (age 2).

        # Build list of columns to check for this outcome
        check_cols <- y_base  # Default: check raw column
        if (y_base == "dtap_3dos_utd_by_18m") {
          # Vaccination measured at ~18 months - data appears at rel_time2 (age 2)
          check_cols <- paste0(y_base, "_rel_time2")
        }
        # Add other age-dependent outcomes here if needed

        # Check if any of the specified columns have data
        has_data_check <- function(dt, cols) {
          for (col in cols) {
            if (col %in% names(dt)) {
              if (y_base %in% c("dft_score", "ds_score")) {
                if (any(!is.na(dt[[col]]))) return(TRUE)
              } else {
                if (any(dt[[col]] > 0, na.rm = TRUE)) return(TRUE)
              }
            }
          }
          return(FALSE)
        }

        valid_reporters <- DT_Stacked_Reg[, .(
          has_data = has_data_check(.SD, check_cols)
        ), by = deso_2021][has_data == TRUE, deso_2021]

        # Subset the datasets
        DT_Analytic <- DT_Stacked_Reg[deso_2021 %in% valid_reporters]
        DT_2018_Sub <- DT_2018_Reg[deso_2021 %in% valid_reporters]
        DT_2019_Sub <- DT_2019_Reg[deso_2021 %in% valid_reporters]
        DT_2020_Sub <- DT_2020_Reg[deso_2021 %in% valid_reporters]
        DT_Bal_Sub  <- DT_Bal_Reg[deso_2021 %in% valid_reporters]

        # SAFETY CHECK: Enough Treated units?
        n_trt_remaining <- DT_Analytic[deso_treated_binary == 1, uniqueN(deso_2021)]

        if (n_trt_remaining < 10) {
          if(use_covariates) cat(sprintf("   SKIPPING %s: Too few treated units (%d) after quality filter.\n", y_base, n_trt_remaining))
          next
        }

        # Recalculate Weights for this subset - DYNAMIC cohorts
        N_treated_sub <- DT_Analytic[is.finite(first_treat_year), .(N = uniqueN(deso_2021[deso_treated_binary == 1])), by = cohort]
        N_by_cohort_sub <- data.table(cohort = active_cohorts, N_treated = 0)
        if(nrow(N_treated_sub)>0) N_by_cohort_sub[match(N_treated_sub$cohort, cohort), N_treated := N_treated_sub$N]

        if(use_covariates) cat(sprintf("   [Smart Filter] %s: Analyzed on %d Valid DeSOs (Treated n=%d).\n",
                                       label_base, uniqueN(DT_Analytic$deso_2021), n_trt_remaining))

      } else {
        # TIER 1: Use Full Region Data
        DT_Analytic <- DT_Stacked_Reg
        DT_2018_Sub <- DT_2018_Reg
        DT_2019_Sub <- DT_2019_Reg
        DT_2020_Sub <- DT_2020_Reg
        DT_Bal_Sub  <- DT_Bal_Reg
        N_by_cohort_sub <- N_by_cohort_reg
      }
      # --- END SMART FILTER ---

      # Use horizon_range from birth-time-only check above
      variants <- make_outcome_variants(y_base, label_base, horizon_k = horizon_range)

      for (v in variants) {
        y <- v$y; rel_key <- v$rel_key; label <- v$label; note <- note_base

        rhs_twfe <- build_rhs(treat_twfe,  covar_str)
        rhs_sunab <- build_rhs(treat_sunab, covar_str)
        rhs_cohtrt <- build_rhs(treat_cohtrt, covar_str)

        TWFE_disagg_fml <- build_fml(y, rhs_twfe,  fe_disagg)
        fml_cohtrt      <- build_fml(y, rhs_cohtrt, fe_stacked)
        Sunab_fml       <- build_fml(y, rhs_sunab,  fe_stacked)

        # A: Separate TWFE (Safe) - WITH MATCHING WEIGHTS
        # Weights account for control DeSOs used multiple times in matching with replacement
        safe_feols <- function(fml, data) {
          tryCatch({
            if ("match_weight" %in% names(data)) {
              fixest::feols(fml, data = data, weights = ~match_weight, cluster = ~ deso_2021)
            } else {
              fixest::feols(fml, data = data, cluster = ~ deso_2021)
            }
          }, error = function(e) NULL)
        }
        TWFE_2018 <- safe_feols(TWFE_disagg_fml, DT_2018_Sub)
        TWFE_2019 <- safe_feols(TWFE_disagg_fml, DT_2019_Sub)
        TWFE_2020 <- safe_feols(TWFE_disagg_fml, DT_2020_Sub)

        list_results <- list()
        if (!is.null(TWFE_2018)) list_results[["2018"]] <- merge_fixest_pvals(extract_twfe_terms(TWFE_2018, "2018"), TWFE_2018)
        if (!is.null(TWFE_2019)) list_results[["2019"]] <- merge_fixest_pvals(extract_twfe_terms(TWFE_2019, "2019"), TWFE_2019)
        if (!is.null(TWFE_2020)) list_results[["2020"]] <- merge_fixest_pvals(extract_twfe_terms(TWFE_2020, "2020"), TWFE_2020)

        if (length(list_results) == 0) next

        TWFE_res_ALL <- rbindlist(list_results, use.names = TRUE, fill = TRUE)
        N_w_subset <- N_by_cohort_sub[cohort %in% names(list_results)]
        TWFE_Agg_W <- add_sig_from_se(agg_eventstudy(TWFE_res_ALL, N_w_subset), "estimate", "se")

        # B: Stacked
        stacked_TWFE <- safe_feols(fml_cohtrt, DT_Analytic)
        stacked_Agg_W <- if(!is.null(stacked_TWFE)) add_sig_from_se(agg_eventstudy(extract_stacked_model_terms(stacked_TWFE), N_by_cohort_sub), "estimate", "se") else NULL

        # C: Balanced
        bal_TWFE <- safe_feols(fml_cohtrt, DT_Bal_Sub)
        bal_Agg_W <- if(!is.null(bal_TWFE)) add_sig_from_se(agg_eventstudy(extract_stacked_model_terms(bal_TWFE), N_by_cohort_sub), "estimate", "se") else NULL

        # D: Sunab
        sunab_mod <- safe_feols(Sunab_fml, DT_Analytic)
        sunab_res <- if(!is.null(sunab_mod)) merge_fixest_pvals(extract_sunab_terms(sunab_mod), sunab_mod, "estimate", "std.error") else NULL

        # Save
        RESULTS[[covars_on]][[rel_key]][[y_base]] <- list(
          TWFE_weighted = TWFE_Agg_W,
          stacked_TWFE_weighted = stacked_Agg_W,
          balanced_stacked_TWFE_weighted = bal_Agg_W,
          SUNAB = sunab_res,
          SUNAB_model = sunab_mod,  # Stored for pre-trends testing
          TWFE_2018 = list_results[["2018"]],
          TWFE_2019 = list_results[["2019"]],
          TWFE_2020 = list_results[["2020"]],
          TWFE_per_cohort = TWFE_res_ALL
        )
      }
    }
  }

  # Dose-response analysis removed in favour of threshold sensitivity.
  # Threshold sensitivity is handled via SENSITIVITY_DOSAGE in config.R (60%, 70%, 80%).

  # ==========================================================================
  # PLOTTING SECTION (Adjusts for Smart Filter Drops)
  # ==========================================================================

  # 1. Identify which outcomes actually generated results
  # Check the "covars_on" / "age0" slot to see what keys exist
  successful_outcomes <- names(RESULTS[["covars_on"]][["age0"]])

  # Intersect with base list to maintain correct order
  outcomes_to_plot <- intersect(base_outcomes, successful_outcomes)
  labels_to_plot   <- pretty_labels[outcomes_to_plot]

  if (length(outcomes_to_plot) == 0) {
    cat(sprintf("WARNING: No valid results generated for region %s. Skipping plots.\n", reg_name))
    next
  }

  cat(sprintf("Generating plots for %d successful outcomes (skipped %d)...\n",
              length(outcomes_to_plot), length(base_outcomes) - length(outcomes_to_plot)))

  # 2. Create Directories
  dir_sanity_reg <- file.path(ROOT_OUTPUT_DIR, "plots", reg_name, "sanity_check")
  dir_cohort_reg <- file.path(ROOT_OUTPUT_DIR, "plots", reg_name, "cohort_plots")
  dir.create(dir_sanity_reg, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_cohort_reg, recursive = TRUE, showWarnings = FALSE)

  # 3. Generate SANITY CHECK Master Plot
  cat(sprintf("   [PLOT] Generating Master Plot for %s (%d outcomes)...\n", reg_name, length(outcomes_to_plot)))
  tryCatch({
    PLOTS_on  <- rebuild_PLOTS_from_RESULTS(RESULTS, outcomes_to_plot, labels_to_plot,
                                            covars_on = TRUE,
                                            include = c("TWFE_weighted", "stacked_TWFE_weighted", "SUNAB"),
                                            plot_name = "cs", connect = FALSE)

    if (is.null(PLOTS_on) || length(PLOTS_on) == 0) {
      cat("   [PLOT] WARNING: rebuild_PLOTS_from_RESULTS returned NULL/empty\n")
    } else {
      master_on <- make_master_grid(PLOTS_on, outcomes_to_plot, labels_to_plot,
                                    covars_on = TRUE, plot_name = "cs",
                                    grid_title = paste(reg_name, "- With Covariates"))

      if (is.null(master_on)) {
        cat("   [PLOT] WARNING: make_master_grid returned NULL\n")
      } else {
        png_path <- file.path(dir_sanity_reg, paste0("master_covars_on_", reg_name, ".png"))
        pdf_path <- file.path(dir_sanity_reg, paste0("master_covars_on_", reg_name, ".pdf"))

        # Force delete existing files to get a fresh save
        if (file.exists(png_path)) file.remove(png_path)
        if (file.exists(pdf_path)) file.remove(pdf_path)

        save_master(master_on, PLOTS_on, covars_on = TRUE,
                    file_pdf = pdf_path,
                    file_png = png_path)

        if (file.exists(png_path)) {
          cat(sprintf("   [PLOT] SUCCESS: Saved %s (%.1f KB)\n", basename(png_path), file.size(png_path)/1024))
        } else {
          cat(sprintf("   [PLOT] FAILED: File not created: %s\n", png_path))
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("   [PLOT] ERROR creating Master Plot: %s\n", e$message))
    cat(sprintf("   [PLOT] Error call: %s\n", deparse(e$call)[1]))
  })

  # 4. Generate COHORT PLOTS
  cat(sprintf("   Creating Per-Cohort Grids for %s...\n", reg_name))

  # A. With Covariates
  tryCatch({
    g_rel0_on <- make_per_cohort_twfe_grid(
      RESULTS, outcomes_to_plot, labels_to_plot,
      rel_key = "age0", covars_on = TRUE
    )

    if (!is.null(g_rel0_on)) {
      save_patchwork(
        g_rel0_on,
        file.path(dir_cohort_reg, paste0("per_cohort_age0_covars_on_", reg_name, ".pdf")),
        n_panels = length(outcomes_to_plot), ncol = 1, cell_w = 12, cell_h = 2.8
      )
      cat("   - Saved: Covariates ON grid\n")
    } else {
      cat("   - Warning: g_rel0_on was NULL (could not build grid)\n")
    }
  }, error = function(e) cat(sprintf("   Error creating Cohort Grid (On): %s\n", e$message)))

  # B. Without Covariates
  tryCatch({
    g_rel0_off <- make_per_cohort_twfe_grid(
      RESULTS, outcomes_to_plot, labels_to_plot,
      rel_key = "age0", covars_on = FALSE
    )

    if (!is.null(g_rel0_off)) {
      save_patchwork(
        g_rel0_off,
        file.path(dir_cohort_reg, paste0("per_cohort_age0_covars_off_", reg_name, ".pdf")),
        n_panels = length(outcomes_to_plot), ncol = 1, cell_w = 12, cell_h = 2.8
      )
      cat("   - Saved: Covariates OFF grid\n")
    }
  }, error = function(e) cat(sprintf("   Error creating Cohort Grid (Off): %s\n", e$message)))

  # Save this region's results into the master container before the loop restarts
  MASTER_REGIONAL_RESULTS[[reg_name]] <- RESULTS

  # Clean up memory
  rm(DT_2018_Reg, DT_2019_Reg, DT_2020_Reg, DT_Stacked_Reg, DT_Bal_Reg, DT_S_Reg); gc()
}


# ============================================================================
# POST-ANALYSIS CLEANUP: PRUNE EMPTY OUTCOMES
# ============================================================================
# The Smart Filter might have skipped some Tier 2 outcomes (e.g., Smoking).
# Remove them from 'base_outcomes' so the plotting functions do not crash.

cat("\n--- Checking which outcomes generated results ---\n")

# Check which outcomes exist in the results list
# Check ALL age keys, not just age0: some outcomes like dtap only have results at age2
if (exists("MASTER_REGIONAL_RESULTS") && length(MASTER_REGIONAL_RESULTS) > 0) {
  first_region <- names(MASTER_REGIONAL_RESULTS)[1]
  all_age_keys <- c("age0", "age1", "age2", "age3", "age01sum", "age0_3sum")
  successful_outcomes <- unique(unlist(lapply(all_age_keys, function(ak) {
    names(MASTER_REGIONAL_RESULTS[[first_region]][["covars_on"]][[ak]])
  })))
} else {
  all_age_keys <- c("age0", "age1", "age2", "age3", "age01sum", "age0_3sum")
  successful_outcomes <- unique(unlist(lapply(all_age_keys, function(ak) {
    names(RESULTS[["covars_on"]][[ak]])
  })))
}

# Identify dropped outcomes
dropped_outcomes <- setdiff(base_outcomes, successful_outcomes)

if (length(dropped_outcomes) > 0) {
  cat("WARNING: The following outcomes were skipped (Smart Filter) and will be excluded from plots/tables:\n")
  print(dropped_outcomes)

  # Update the master lists to only include successful outcomes
  base_outcomes <- intersect(base_outcomes, successful_outcomes)
  pretty_labels <- pretty_labels[base_outcomes]

} else {
  cat("Success: All outcomes in the list generated results.\n")
}

cat(sprintf("Proceeding to plots with %d active outcomes.\n", length(base_outcomes)))

# ============================================================================


# ==============================================================================#
# =================== HETEROGENEITY: QUARTILE ANALYSIS =========================#
# ==============================================================================#

cat("\n\n--- Starting Heterogeneity Analysis: Top vs Bottom CNI Quartiles ---\n")

# 1. Extract Baseline CNI (one year before treatment) - DYNAMIC cohorts
baseline_cni_map <- data.table()
for (cohort_year in active_cohorts) {  # Use active_cohorts from config
  match_obj <- match_results_list_stratified[[cohort_year]]$match_object
  cohort_data <- as.data.table(match.data(match_obj))
  baseline_data <- cohort_data[year == (as.numeric(cohort_year) - 1)]
  baseline_cni_map <- rbind(baseline_cni_map, baseline_data[, .(deso_id, cni_score)])
}
baseline_cni_map <- unique(baseline_cni_map)

# 2. Calculate Quartiles
q25 <- quantile(baseline_cni_map$cni_score, 0.25, na.rm = TRUE)
q75 <- quantile(baseline_cni_map$cni_score, 0.75, na.rm = TRUE)

cat(sprintf("--- CNI Thresholds: Bottom Quartile < %.3f | Top Quartile > %.3f ---\n", q25, q75))

# 3. Define Groups (Drop the middle 50%)
baseline_cni_map[, cni_group := fcase(
  cni_score <= q25, "Low_CNI",
  cni_score >= q75, "High_CNI",
  default = NA_character_
)]

# Filter list of DeSOs
high_cni_desos <- baseline_cni_map[cni_group == "High_CNI", deso_id]
low_cni_desos  <- baseline_cni_map[cni_group == "Low_CNI", deso_id]

# 4. Create Sub-Datasets
Stacked_data_individual_level[, cni_group := fcase(
  deso_2021 %in% high_cni_desos, "High_CNI",
  deso_2021 %in% low_cni_desos, "Low_CNI",
  default = "Middle"
)]

Stacked_data_high_cni <- Stacked_data_individual_level[cni_group == "High_CNI"]
Stacked_data_low_cni  <- Stacked_data_individual_level[cni_group == "Low_CNI"]

cat(sprintf("--- Sample Sizes: High CNI (Top 25%%): %d | Low CNI (Bottom 25%%): %d ---\n",
            nrow(Stacked_data_high_cni), nrow(Stacked_data_low_cni)))

# 5. Run Sun & Abraham Models
HET_RESULTS <- list()

for (use_covariates in c(TRUE, FALSE)) {
  covars_on <- if (use_covariates) "covars_on" else "covars_off"
  covar_str <- if (use_covariates) paste(covars, collapse = " + ") else NULL

  for (outcome_idx in seq_along(base_outcomes)) {
    y_base <- base_outcomes[outcome_idx]
    label_base <- pretty_labels[outcome_idx]
    variants <- make_outcome_variants(y_base, label_base, horizon_k = 0:3)

    for (v in variants) {
      y <- v$y; rel_key <- v$rel_key
      rhs_sunab <- build_rhs(treat_sunab, covar_str)
      Sunab_fml <- build_fml(y, rhs_sunab, fe_stacked)

      # High CNI - WITH MATCHING WEIGHTS
      tryCatch({
        mod_high <- feols(Sunab_fml, data = Stacked_data_high_cni,
                          weights = ~match_weight, cluster = ~ deso_2021)
        HET_RESULTS[[covars_on]][[rel_key]][[y_base]][["High_CNI"]] <- merge_fixest_pvals(extract_sunab_terms(mod_high), mod_high)
      }, error = function(e) NULL)

      # Low CNI - WITH MATCHING WEIGHTS
      tryCatch({
        mod_low <- feols(Sunab_fml, data = Stacked_data_low_cni,
                         weights = ~match_weight, cluster = ~ deso_2021)
        HET_RESULTS[[covars_on]][[rel_key]][[y_base]][["Low_CNI"]] <- merge_fixest_pvals(extract_sunab_terms(mod_low), mod_low)
      }, error = function(e) NULL)
    }
  }
}
cat("--- Heterogeneity analysis (Quartiles) complete. ---\n")






#------------------------------------------------------------------------------#
# Create plots after the for-loop
#------------------------------------------------------------------------------#

# ==== Sanity check master plot ====

# Choose which models to display in the master plot by supplying them to the
# "include" input of rebuild_PLOTS_from_RESULTS.
# Available model keys:
#------------------------------------------------------------------------------#
# "TWFE_weighted"                   Weighted average of disaggregated models
# "stacked_TWFE_weighted"           Stacked model
# "balanced_stacked_TWFE_weighted"  Stacked model (balanced panel)
# "SUNAB"                           Sun and Abraham estimator
# "CSA_dynamic_weighted"            Callaway Sant'Anna (IPW)
#------------------------------------------------------------------------------#

# Wrap entire post-loop plotting in tryCatch so tables still generate even if plots fail
tryCatch({
  cat("\n--- Starting Post-Loop Plotting ---\n")

  # Five models comparison.
  # Stacked TWFE is the PRIMARY estimator; Sun-Abraham and others are for robustness.
  PLOTS_on  <- rebuild_PLOTS_from_RESULTS(
    RESULTS, base_outcomes, pretty_labels,
    covars_on = TRUE,
    include   = c("stacked_TWFE_weighted", "TWFE_weighted",
                  "balanced_stacked_TWFE_weighted", "SUNAB",
                  "CSA_dynamic_weighted" ),
    plot_name = "cs",
    connect   = FALSE
  )
  PLOTS_off <- rebuild_PLOTS_from_RESULTS(
    RESULTS, base_outcomes, pretty_labels,
    covars_on = FALSE,
    include   = c("stacked_TWFE_weighted", "TWFE_weighted",
                  "balanced_stacked_TWFE_weighted", "SUNAB",
                  "CSA_dynamic_weighted" ),
    plot_name = "cs",
    connect   = FALSE
  )
  PLOTS <- PLOTS_on
  PLOTS$covars_off <- PLOTS_off$covars_off

  master_on_5models  <- make_master_grid(PLOTS, base_outcomes, pretty_labels,
                                         covars_on = TRUE,  plot_name = "cs",
                                         grid_title = "Five models — with covariates")

master_off_5models <- make_master_grid(PLOTS, base_outcomes, pretty_labels,
                                       covars_on = FALSE, plot_name = "cs",
                                       grid_title = "Five models — no covariates")

# Four models (excluding Callaway Sant'Anna)
PLOTS_on  <- rebuild_PLOTS_from_RESULTS(
  RESULTS, base_outcomes, pretty_labels,
  covars_on = TRUE,
  include   = c("TWFE_weighted","stacked_TWFE_weighted",
                "balanced_stacked_TWFE_weighted", "SUNAB"),
  plot_name = "cs",
  connect   = FALSE
)
PLOTS_off <- rebuild_PLOTS_from_RESULTS(
  RESULTS, base_outcomes, pretty_labels,
  covars_on = FALSE,
  include   = c("TWFE_weighted","stacked_TWFE_weighted",
                "balanced_stacked_TWFE_weighted", "SUNAB"),
  plot_name = "cs",
  connect   = FALSE
)
PLOTS <- PLOTS_on
PLOTS$covars_off <- PLOTS_off$covars_off

master_on_4models  <- make_master_grid(PLOTS, base_outcomes, pretty_labels,
                                       covars_on = TRUE,  plot_name = "cs",
                                       grid_title = "Four models — with covariates")

master_off_4models <- make_master_grid(PLOTS, base_outcomes, pretty_labels,
                                       covars_on = FALSE, plot_name = "cs",
                                       grid_title = "Four models — no covariates")

#----------------------------------------------------------------------------#
# NOTE: Hollow points indicate that not all cohorts contribute to that
# event-time. All three cohorts overlap in event times (-6, 1). Data from
# event times (-8, -7) and (2, 3) are cohort-incomplete. Verify:
#
# common_event_times <- Reduce(
#  intersect,
#  list(
#    unique(analysis_data_2018_individuals_final$event_time),
#    unique(analysis_data_2019_individuals_final$event_time),
#    unique(analysis_data_2020_individuals_final$event_time)
#  )
#  )
#  sort(common_event_times)
#
#----------------------------------------------------------------------------#

# ==== Comparison plots  ====

# Produces plots comparing model results with and without covariates for a
# specified outcome horizon.
# Parameters:
#   model_key: which model to compare
#   rel_key (age_key): which child age to examine
#   connect (T/F): connect point estimates with lines
#   ncol: number of columns

#------------------------------------------------------------------------------#
# Model keys:
# "TWFE_weighted"                   Weighted average of disaggregated models
# "stacked_TWFE_weighted"           Stacked model
# "balanced_stacked_TWFE_weighted"  Stacked model (balanced panel)
# "SUNAB"                           Sun and Abraham estimator
# "CSA_dynamic_weighted"            Callaway Sant'Anna (IPW)
#
# Rel keys:
# "age0"            First year
# "rel1"            Second year
# "rel2"            Third year
# "rel3"            Fourth year
# "age01sum"        First AND second year summed
# "age0_3sum"       All years summed
#------------------------------------------------------------------------------#



# Stacked TWFE, first year (age0)
g1 <- make_model_covariate_compare_grid(
  RESULTS, base_outcomes, pretty_labels,
  rel_key = "age0",
  model_key = "stacked_TWFE_weighted",
  title_size = 18
)

# Stacked TWFE, first two years (age01sum)
g2 <- make_model_covariate_compare_grid(
  RESULTS, base_outcomes, pretty_labels,
  rel_key = "age01sum",
  model_key = "stacked_TWFE_weighted",
  title_size = 18
)

# Stacked TWFE, ALL years (age0_3sum)
g3 <- make_model_covariate_compare_grid(
  RESULTS, base_outcomes, pretty_labels,
  rel_key = "age0_3sum",
  model_key = "stacked_TWFE_weighted",
  connect = FALSE,
  title_size = 18,
  hollow_ticks = c(-8, -7, 0, 1)   # custom hollowing for incomplete-sample event-times
)

# ==== Per cohort plots ====

# Per cohort plots for rel_time0
g_rel0_on <- make_per_cohort_twfe_grid(
  RESULTS, base_outcomes, pretty_labels,
  rel_key = "age0", covars_on = TRUE
)

# Same rel_time0, no covariates
g_rel2_off <- make_per_cohort_twfe_grid(
  RESULTS, base_outcomes, pretty_labels,
  rel_key = "age0", covars_on = FALSE
)

#-------------------------------------------------------------------------------



#------------------------------------------------------------------------------#
# Saving the plots
#------------------------------------------------------------------------------#

# ==== Folders ====
dir_sanity <- file.path(ROOT_OUTPUT_DIR, "plots", "sanity_check")
dir_fig13  <- file.path(ROOT_OUTPUT_DIR, "plots", "fig_1_3")
dir_cohort <- file.path(ROOT_OUTPUT_DIR, "plots", "cohort_plots")

invisible(lapply(list(dir_sanity, dir_fig13, dir_cohort),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# ==== Master grids -> plots/sanity_check ====
save_master(
  master_on_5models,  PLOTS, covars_on = TRUE,
  file_pdf = file.path(dir_sanity, "master_covars_on_5models.pdf"),
  file_png = file.path(dir_sanity, "master_covars_on_5models.png")
)

save_master(
  master_off_5models, PLOTS, covars_on = FALSE,
  file_pdf = file.path(dir_sanity, "master_covars_off_5models.pdf"),
  file_png = file.path(dir_sanity, "master_covars_off_5models.png")
)

save_master(
  master_on_4models,  PLOTS, covars_on = TRUE,
  file_pdf = file.path(dir_sanity, "master_covars_on_4models.pdf"),
  file_png = file.path(dir_sanity, "master_covars_on_4models.png")
)

save_master(
  master_off_4models, PLOTS, covars_on = FALSE,
  file_pdf = file.path(dir_sanity, "master_covars_off_4models.pdf"),
  file_png = file.path(dir_sanity, "master_covars_off_4models.png")
)

# ==== Covariate comparison grids -> plots/fig_1_3 ====
n_outcomes <- length(base_outcomes)

save_patchwork(
  g1, file.path(dir_fig13, "cov_comp_sunab_rel0.pdf"),
  n_panels = n_outcomes, ncol = 3
)

save_patchwork(
  g2, file.path(dir_fig13, "cov_comp_sunab_rel01sum.pdf"),
  n_panels = n_outcomes, ncol = 3
)

save_patchwork(
  g3, file.path(dir_fig13, "cov_comp_sunab_rel0_3sum.pdf"),
  n_panels = n_outcomes, ncol = 3
)


# ==== Per-cohort grids -> plots/cohort_plots ====
save_patchwork(
  g_rel0_on, file.path(dir_cohort, "per_cohort_age0_covars_on.pdf"),
  n_panels = n_outcomes, ncol = 1, cell_w = 12, cell_h = 2.8
)

save_patchwork(
  g_rel2_off, file.path(dir_cohort, "per_cohort_age0_covars_off.pdf"),
  n_panels = n_outcomes, ncol = 1, cell_w = 12, cell_h = 2.8
)


#-------------------------------------------------------------------------------#




# ==============================================================================#
# =================== HETEROGENEITY PLOTTING SECTION ===========================#
# ==============================================================================#

cat("\n--- Generating Heterogeneity Plots ---\n")

# 1. Create directory
dir_het    <- file.path(ROOT_OUTPUT_DIR, "plots", "heterogeneity")
dir.create(dir_het, recursive = TRUE, showWarnings = FALSE)

# 2. Plotting function for High vs Low CNI comparison
make_cni_compare_grid <- function(
    HET_RESULTS,
    base_outcomes,
    pretty_labels,
    rel_key,
    covars_on = TRUE,
    ncol = 3,
    title_size = 18,
    subtitle_fill = "#E8F1FD",
    subtitle_box_col = "grey60"
) {

  cov_key <- if (covars_on) "covars_on" else "covars_off"

  # Dynamic title based on the time horizon (rel_key)
  rel_title <- switch(rel_key,
                      "age0"      = "First Year (Age 0)",
                      "rel1"      = "Second Year (Age 1)",
                      "age01sum"  = "Sum of Years 0-1",
                      "age0_3sum" = "Cumulative Sum (Years 0-3)",
                      rel_key) # Fallback

  grid_title <- sprintf("Heterogeneity by Baseline CNI - %s", rel_title)

  plot_list <- list()

  for (i in seq_along(base_outcomes)) {
    y_base <- base_outcomes[i]
    pretty <- pretty_labels[i]

    # Safely extract the High and Low tables
    dt_high <- tryCatch(normalize_plot_dt(HET_RESULTS[[cov_key]][[rel_key]][[y_base]][["High_CNI"]]), error = function(e) NULL)
    dt_low  <- tryCatch(normalize_plot_dt(HET_RESULTS[[cov_key]][[rel_key]][[y_base]][["Low_CNI"]]), error = function(e) NULL)

    # Combine into a list for the plotter
    dts <- list()
    if (!is.null(dt_high) && nrow(dt_high) > 0) dts[["High CNI (Vulnerable)"]] <- dt_high
    if (!is.null(dt_low)  && nrow(dt_low)  > 0) dts[["Low CNI (Less Vulnerable)"]]  <- dt_low

    if (length(dts) == 0) next

    # Create the single cell plot
    p <- make_event_study_plot(
      dts = dts,
      label = pretty,
      note = NULL,
      connect = TRUE,
      legend_title = "Subgroup",
      title_prefix = "",
      palette_line = c("High CNI (Vulnerable)" = "#D55E00", "Low CNI (Less Vulnerable)" = "#0072B2"),
      palette_fill = c("High CNI (Vulnerable)" = "#D55E00", "Low CNI (Less Vulnerable)" = "#0072B2")
    )

    # Add box header styling
    p <- p +
      labs(title = NULL, subtitle = pretty) +
      theme(
        plot.subtitle = ggtext::element_textbox_simple(
          padding    = margin(4, 6, 4, 6),
          margin     = margin(b = 6),
          fill       = subtitle_fill,
          box.colour = subtitle_box_col,
          linewidth  = 0.3,
          r          = unit(2, "pt"),
          size       = 11,
          face       = "bold"
        ),
        legend.position = "bottom"
      )

    plot_list[[i]] <- p
  }

  # Remove NULLs
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) return(NULL)

  # Stitch together using Patchwork
  pw <- patchwork::wrap_plots(plot_list, ncol = ncol) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = grid_title,
      theme = theme(
        plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
        legend.position = "bottom"
      )
    )

  return(pw)
}


# 3. Execution Loop: Generate and Save Plots for Key Horizons

horizons_to_plot <- c("age0", "age01sum", "age0_3sum")

for (rk in horizons_to_plot) {

  cat(sprintf("   Generating grid for horizon: %s...\n", rk))

  g_het <- make_cni_compare_grid(
    HET_RESULTS,
    base_outcomes,
    pretty_labels,
    rel_key = rk,
    covars_on = TRUE
  )

  if (!is.null(g_het)) {
    file_name <- file.path(dir_het, paste0("heterogeneity_sunab_", rk, ".pdf"))

    save_patchwork(
      g_het,
      file_name,
      n_panels = length(base_outcomes),
      ncol = 3,
      also_png = TRUE
    )
    cat(sprintf("   Saved: %s\n", file_name))
  }
}

cat("\n--- All Heterogeneity Plots Saved to 'plots/heterogeneity' ---\n")

}, error = function(e) {
  cat(sprintf("\n--- WARNING: Post-analysis plotting failed: %s ---\n", e$message))
  cat("--- Continuing to table extraction... ---\n\n")
})






# ==============================================================================#
# TABLE EXTRACTION - COMPLETE SELF-CONTAINED BLOCK
# ==============================================================================#

cat("\n", rep("=", 70), "\n", sep = "")
cat("STARTING TABLE EXTRACTION\n")
cat(rep("=", 70), "\n\n", sep = "")

# ==============================================================================
# DEFINE ALL HELPER FUNCTIONS
# ==============================================================================

# --- Function 1: Extract all results to master table ---
extract_all_results_to_table <- function(RESULTS,
                                         model_keys = c("TWFE_weighted",
                                                        "stacked_TWFE_weighted",
                                                        "balanced_stacked_TWFE_weighted",
                                                        "SUNAB",
                                                        "TWFE_2018",
                                                        "TWFE_2019",
                                                        "TWFE_2020")) {
  all_rows <- list()
  row_idx <- 1

  for (cov_key in names(RESULTS)) {
    for (rel_key in names(RESULTS[[cov_key]])) {
      for (y_base in names(RESULTS[[cov_key]][[rel_key]])) {
        for (model_key in model_keys) {

          res <- NULL
          res <- tryCatch({ RESULTS[[cov_key]][[rel_key]][[y_base]][[model_key]] }, error = function(e) NULL)

          if (is.null(res)) next
          if (!is.data.frame(res)) next
          if (nrow(res) == 0) next

          dt <- copy(as.data.table(res))

          if ("std.error" %in% names(dt) && !"se" %in% names(dt)) {
            setnames(dt, "std.error", "se")
          }

          dt[, covariates := cov_key]
          dt[, time_horizon := rel_key]
          dt[, outcome := y_base]
          dt[, model := model_key]

          all_rows[[row_idx]] <- dt
          row_idx <- row_idx + 1
        }
      }
    }
  }

  cat(sprintf("   Collected %d result tables.\n", length(all_rows)))

  if (length(all_rows) == 0) {
    warning("No results extracted!")
    return(NULL)
  }

  master <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)

  first_cols <- c("covariates", "time_horizon", "outcome", "model", "event_time",
                  "estimate", "se", "p.value", "conf.low", "conf.high")
  existing_first <- intersect(first_cols, names(master))
  other_cols <- setdiff(names(master), existing_first)
  setcolorder(master, c(existing_first, other_cols))

  return(master)
}

# --- Function 2: Create summary table (wide format) ---
make_summary_table <- function(master_table,
                               rel_key_filter = "age0",
                               cov_filter = "covars_on",
                               event_times = c(-1, 0, 1, 2, 3)) {

  dt <- master_table[covariates == cov_filter & time_horizon == rel_key_filter]
  if (nrow(dt) == 0) return(NULL)
  dt <- dt[event_time %in% event_times]

  # Order outcomes: PRIMARY first, then SECONDARY (from outcomes_config.R)
  outcome_order <- c(OUTCOMES_PRIMARY, OUTCOMES_SECONDARY)
  dt[, outcome := factor(outcome, levels = outcome_order)]
  dt <- dt[!is.na(outcome)]  # Drop outcomes not in the ordered list

  # Add outcome tier indicator (for tiered FDR adjustment)
  dt[, outcome_tier := fifelse(as.character(outcome) %in% OUTCOMES_PRIMARY, "Primary", "Secondary")]

  # Tiered Benjamini-Hochberg FDR adjustment
  # Apply FDR correction SEPARATELY within primary and secondary outcomes.
  # This prevents secondary outcomes from penalizing primary outcome significance.
  # Reference: Benjamini & Hochberg (1995)
  dt[outcome_tier == "Primary", p.value.fdr := p.adjust(p.value, method = "BH")]
  dt[outcome_tier == "Secondary", p.value.fdr := p.adjust(p.value, method = "BH")]

  dt[, est_formatted := sprintf("%.3f%s", estimate,
                                fifelse(p.value < 0.01, "***",
                                        fifelse(p.value < 0.05, "**",
                                                fifelse(p.value < 0.1, "*", ""))))]
  dt[, se_formatted := sprintf("(%.3f)", se)]
  dt[, cell := paste0(est_formatted, "\n", se_formatted)]
  dt[, fdr_formatted := sprintf("[FDR: %.3f]", p.value.fdr)]

  wide <- dcast(dt, outcome + model ~ event_time, value.var = "cell", fill = "")
  return(wide)
}

# --- Function 2b: Create FDR summary table (wide format with p-values and FDR) ---
make_fdr_summary_table <- function(master_table,
                                   rel_key_filter = "age0",
                                   cov_filter = "covars_on",
                                   model_filter = "stacked_TWFE_weighted",
                                   event_times = c(0, 1)) {

  dt <- master_table[covariates == cov_filter &
                       time_horizon == rel_key_filter &
                       model == model_filter]
  if (nrow(dt) == 0) return(NULL)
  dt <- dt[event_time %in% event_times]

  # Order outcomes: PRIMARY first, then SECONDARY
  outcome_order <- c(OUTCOMES_PRIMARY, OUTCOMES_SECONDARY)
  dt[, outcome := factor(outcome, levels = outcome_order)]
  dt <- dt[!is.na(outcome)]

  # Add outcome tier indicator
  dt[, outcome_tier := fifelse(as.character(outcome) %in% OUTCOMES_PRIMARY, "Primary", "Secondary")]

  # Tiered Benjamini-Hochberg FDR adjustment
  dt[outcome_tier == "Primary", p.value.fdr := p.adjust(p.value, method = "BH")]
  dt[outcome_tier == "Secondary", p.value.fdr := p.adjust(p.value, method = "BH")]

  # Add significance indicators
  dt[, sig_raw := fifelse(p.value < 0.01, "***",
                          fifelse(p.value < 0.05, "**",
                                  fifelse(p.value < 0.1, "*", "")))]
  dt[, sig_fdr := fifelse(p.value.fdr < 0.01, "***",
                          fifelse(p.value.fdr < 0.05, "**",
                                  fifelse(p.value.fdr < 0.1, "*", "")))]

  # Create formatted columns
  dt[, estimate_fmt := sprintf("%.4f", estimate)]
  dt[, se_fmt := sprintf("%.4f", se)]
  dt[, p_raw_fmt := sprintf("%.4f%s", p.value, sig_raw)]
  dt[, p_fdr_fmt := sprintf("%.4f%s", p.value.fdr, sig_fdr)]

  # Create output table (long format for clarity)
  result <- dt[, .(outcome, event_time, estimate = estimate_fmt, se = se_fmt,
                   p_value = p_raw_fmt, p_value_FDR = p_fdr_fmt)]

  # Add pretty labels
  result[, outcome_label := PRETTY_LABELS_MAP[outcome]]
  setcolorder(result, c("outcome_label", "outcome", "event_time", "estimate", "se", "p_value", "p_value_FDR"))
  setorder(result, outcome, event_time)
  result[, outcome := as.character(outcome)]

  return(result)
}


# --- Function 3: Export to Excel ---
export_results_to_excel <- function(master_table,
                                    file_path = "results_tables.xlsx",
                                    base_outcomes,
                                    pretty_labels) {

  wb <- createWorkbook()

  addWorksheet(wb, "Master_All_Results")
  writeData(wb, "Master_All_Results", master_table)

  summary_rel0 <- make_summary_table(master_table, "age0", "covars_on")
  if (!is.null(summary_rel0)) {
    addWorksheet(wb, "Summary_Year0_CovOn")
    writeData(wb, "Summary_Year0_CovOn", summary_rel0)
  }

  summary_rel01 <- make_summary_table(master_table, "age01sum", "covars_on")
  if (!is.null(summary_rel01)) {
    addWorksheet(wb, "Summary_Years01_CovOn")
    writeData(wb, "Summary_Years01_CovOn", summary_rel01)
  }

  summary_rel03 <- make_summary_table(master_table, "age0_3sum", "covars_on")
  if (!is.null(summary_rel03)) {
    addWorksheet(wb, "Summary_AllYears_CovOn")
    writeData(wb, "Summary_AllYears_CovOn", summary_rel03)
  }

  saveWorkbook(wb, file_path, overwrite = TRUE)
  cat(sprintf("   Excel file saved: %s\n", file_path))
}

# --- Function 4: Extract heterogeneity results ---
extract_heterogeneity_results <- function(HET_RESULTS,
                                          subgroups = c("High_CNI", "Low_CNI")) {
  all_rows <- list()
  row_idx <- 1

  for (cov_key in names(HET_RESULTS)) {
    for (rel_key in names(HET_RESULTS[[cov_key]])) {
      for (y_base in names(HET_RESULTS[[cov_key]][[rel_key]])) {
        for (subgroup in subgroups) {

          res <- NULL
          try({ res <- HET_RESULTS[[cov_key]][[rel_key]][[y_base]][[subgroup]] }, silent = TRUE)

          if (is.null(res)) next
          if (!is.data.frame(res)) next
          if (nrow(res) == 0) next

          dt <- copy(as.data.table(res))
          if ("std.error" %in% names(dt) && !"se" %in% names(dt)) {
            setnames(dt, "std.error", "se")
          }

          dt[, covariates := cov_key]
          dt[, time_horizon := rel_key]
          dt[, outcome := y_base]
          dt[, subgroup := subgroup]

          all_rows[[row_idx]] <- dt
          row_idx <- row_idx + 1
        }
      }
    }
  }

  if (length(all_rows) == 0) return(NULL)

  het_master <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)

  first_cols <- c("covariates", "time_horizon", "outcome", "subgroup", "event_time",
                  "estimate", "se", "p.value", "conf.low", "conf.high")
  existing_first <- intersect(first_cols, names(het_master))
  other_cols <- setdiff(names(het_master), existing_first)
  setcolorder(het_master, c(existing_first, other_cols))

  return(het_master)
}

# --- Function 5: Create paper-ready table ---
make_paper_table <- function(master_table,
                             outcomes_to_show,
                             pretty_labels,
                             model_key = "stacked_TWFE_weighted",
                             cov_filter = "covars_on",
                             rel_key = "age0",
                             event_times = c(0, 1, 2, 3)) {

  dt <- master_table[model == model_key &
                       covariates == cov_filter &
                       time_horizon == rel_key &
                       event_time %in% event_times]

  if (nrow(dt) == 0) return(NULL)

  label_map <- setNames(pretty_labels, outcomes_to_show)
  dt[, outcome_label := label_map[outcome]]
  dt <- dt[!is.na(outcome_label)]

  if (nrow(dt) == 0) return(NULL)

  dt[, estimate_str := sprintf("%.4f", estimate)]
  dt[, se_str := sprintf("(%.4f)", se)]
  dt[, stars := fifelse(p.value < 0.01, "***",
                        fifelse(p.value < 0.05, "**",
                                fifelse(p.value < 0.1, "*", "")))]
  dt[, cell := paste0(estimate_str, stars)]


  # Order outcomes: PRIMARY first, then SECONDARY
  outcome_order_labels <- PRETTY_LABELS_MAP[c(OUTCOMES_PRIMARY, OUTCOMES_SECONDARY)]
  dt[, outcome_label := factor(outcome_label, levels = outcome_order_labels)]
  dt <- dt[!is.na(outcome_label)]

  wide_est <- dcast(dt, outcome_label ~ paste0("t", event_time), value.var = "cell", fill = "")
  wide_se <- dcast(dt, outcome_label ~ paste0("t", event_time), value.var = "se_str", fill = "")

  result_rows <- list()
  for (i in 1:nrow(wide_est)) {
    result_rows[[2*i - 1]] <- wide_est[i]
    se_row <- wide_se[i]
    se_row$outcome_label <- ""
    result_rows[[2*i]] <- se_row
  }

  paper_table <- rbindlist(result_rows, fill = TRUE)
  return(paper_table)
}

# --- Function 6: Print quick summary to console ---
print_quick_summary <- function(master_table,
                                outcome_var,
                                model_key = "stacked_TWFE_weighted",
                                cov_filter = "covars_on") {

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("OUTCOME:", outcome_var, "\n")
  cat("MODEL:", model_key, "| COVARIATES:", cov_filter, "\n")
  cat(rep("=", 70), "\n\n", sep = "")

  dt <- master_table[outcome == outcome_var &
                       model == model_key &
                       covariates == cov_filter]

  for (rk in c("age0", "age01sum", "age0_3sum")) {
    sub <- dt[time_horizon == rk]
    if (nrow(sub) == 0) next

    cat("Time Horizon:", rk, "\n")
    cat(rep("-", 50), "\n", sep = "")

    sub[, sig := fifelse(p.value < 0.01, "***",
                         fifelse(p.value < 0.05, "**",
                                 fifelse(p.value < 0.1, "*", "")))]

    print(sub[, .(event_time,
                  estimate = round(estimate, 4),
                  se = round(se, 4),
                  p.value = round(p.value, 4),
                  sig)], row.names = FALSE)
    cat("\n")
  }
}

# --- Function 7: Create High vs Low CNI comparison table ---
make_cni_comparison_table <- function(het_master,
                                      rel_key_filter = "age0",
                                      cov_filter = "covars_on",
                                      event_times = c(-2, -1, 0, 1, 2, 3),
                                      pretty_labels = NULL,
                                      base_outcomes = NULL) {

  # Filter data
  dt <- het_master[covariates == cov_filter &
                     time_horizon == rel_key_filter &
                     event_time %in% event_times]

  if (nrow(dt) == 0) return(NULL)

  # Format estimates with stars
  dt[, est_star := sprintf("%.4f%s", estimate,
                           fifelse(p.value < 0.01, "***",
                                   fifelse(p.value < 0.05, "**",
                                           fifelse(p.value < 0.1, "*", ""))))]
  dt[, se_fmt := sprintf("(%.4f)", se)]

  # Add pretty labels if provided
  if (!is.null(pretty_labels) && !is.null(base_outcomes)) {
    label_map <- setNames(pretty_labels, base_outcomes)
    dt[, outcome_label := label_map[outcome]]
  } else {
    dt[, outcome_label := outcome]
  }

  # Create separate tables for High and Low CNI
  high_dt <- dt[subgroup == "High_CNI"]
  low_dt <- dt[subgroup == "Low_CNI"]

  # Pivot each to wide format
  if (nrow(high_dt) > 0) {
    high_wide <- dcast(high_dt, outcome_label ~ event_time,
                       value.var = "est_star", fill = "")
    setnames(high_wide, old = names(high_wide)[-1],
             new = paste0("High_t", names(high_wide)[-1]))
  } else {
    high_wide <- NULL
  }

  if (nrow(low_dt) > 0) {
    low_wide <- dcast(low_dt, outcome_label ~ event_time,
                      value.var = "est_star", fill = "")
    setnames(low_wide, old = names(low_wide)[-1],
             new = paste0("Low_t", names(low_wide)[-1]))
  } else {
    low_wide <- NULL
  }

  # Merge High and Low side by side
  if (!is.null(high_wide) && !is.null(low_wide)) {
    comparison <- merge(high_wide, low_wide, by = "outcome_label", all = TRUE)
  } else if (!is.null(high_wide)) {
    comparison <- high_wide
  } else if (!is.null(low_wide)) {
    comparison <- low_wide
  } else {
    return(NULL)
  }

  return(comparison)
}

# --- Function 8: Create detailed CNI comparison (long format) ---
make_cni_comparison_detailed <- function(het_master,
                                         rel_key_filter = "age0",
                                         cov_filter = "covars_on",
                                         event_times = c(-1, 0, 1, 2)) {

  dt <- het_master[covariates == cov_filter &
                     time_horizon == rel_key_filter &
                     event_time %in% event_times]

  if (nrow(dt) == 0) return(NULL)

  # Create wide format: columns for each subgroup's estimate and SE
  wide <- dcast(dt,
                outcome + event_time ~ subgroup,
                value.var = c("estimate", "se", "p.value"),
                fill = NA)

  # Add significance stars
  wide[, High_CNI_star := fifelse(p.value_High_CNI < 0.01, "***",
                                  fifelse(p.value_High_CNI < 0.05, "**",
                                          fifelse(p.value_High_CNI < 0.1, "*", "")))]
  wide[, Low_CNI_star := fifelse(p.value_Low_CNI < 0.01, "***",
                                 fifelse(p.value_Low_CNI < 0.05, "**",
                                         fifelse(p.value_Low_CNI < 0.1, "*", "")))]

  # Calculate difference (High - Low)
  wide[, diff_estimate := estimate_High_CNI - estimate_Low_CNI]

  # Format for display (handle NA values)
  wide[, High_CNI_fmt := fifelse(is.na(estimate_High_CNI), "NA",
                                  sprintf("%.4f%s", estimate_High_CNI, High_CNI_star))]
  wide[, Low_CNI_fmt := fifelse(is.na(estimate_Low_CNI), "NA",
                                 sprintf("%.4f%s", estimate_Low_CNI, Low_CNI_star))]
  wide[, diff_fmt := fifelse(is.na(diff_estimate), "NA",
                              sprintf("%.4f", diff_estimate))]

  # Select final columns
  result <- wide[, .(outcome, event_time,
                     High_CNI = High_CNI_fmt,
                     High_CNI_SE = round(se_High_CNI, 4),
                     Low_CNI = Low_CNI_fmt,
                     Low_CNI_SE = round(se_Low_CNI, 4),
                     Difference = diff_fmt)]

  # Order by PRIMARY then SECONDARY outcomes
  outcome_order <- c(OUTCOMES_PRIMARY, OUTCOMES_SECONDARY)
  result[, outcome := factor(outcome, levels = outcome_order)]
  setorder(result, outcome, event_time)
  result[, outcome := as.character(outcome)]

  return(result)
}

cat("   All helper functions defined.\n\n")

# ==============================================================================
# EXECUTE TABLE EXTRACTION (FIXED FOR REGIONAL LOOP)
# ==============================================================================

dir_tables <- file.path(ROOT_OUTPUT_DIR, "tables")
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

# DEBUG: Check if MASTER_REGIONAL_RESULTS exists and has content
cat("DEBUG: Checking MASTER_REGIONAL_RESULTS...\n")
if (!exists("MASTER_REGIONAL_RESULTS")) {
  cat("ERROR: MASTER_REGIONAL_RESULTS does not exist!\n")
} else {
  cat(sprintf("DEBUG: MASTER_REGIONAL_RESULTS has %d regions: %s\n",
              length(MASTER_REGIONAL_RESULTS),
              paste(names(MASTER_REGIONAL_RESULTS), collapse = ", ")))
}

# Loop through each region stored in MASTER list
for (reg_name in names(MASTER_REGIONAL_RESULTS)) {

  cat(sprintf("\n--- Creating All Tables for Region: %s ---\n", reg_name))

  # 1. Extract this region's data
  current_results <- MASTER_REGIONAL_RESULTS[[reg_name]]

  # DEBUG: Check structure
  cat(sprintf("DEBUG: current_results has %d top-level keys: %s\n",
              length(names(current_results)),
              paste(names(current_results), collapse = ", ")))

  # 2. Convert to a single flat Master Table for this region
  master_table_reg <- extract_all_results_to_table(
    current_results,
    model_keys = c("TWFE_weighted", "stacked_TWFE_weighted", "SUNAB")
  )

  if (!is.null(master_table_reg)) {
    # --- SAVE RAW DATA ---
    fwrite(master_table_reg, file.path(dir_tables, paste0("results_master_", reg_name, ".csv")))

    # --- EXCEL WORKBOOK ---
    export_results_to_excel(
      master_table_reg,
      file_path = file.path(dir_tables, paste0("EHV_results_", reg_name, ".xlsx")),
      base_outcomes = base_outcomes,
      pretty_labels = pretty_labels
    )

    # --- STEP 3: SUMMARY TABLES (Region-Specific) ---
    # A. Year 0 (With Covariates)
    summary_y0_cov <- make_summary_table(master_table_reg, "age0", "covars_on")
    if (!is.null(summary_y0_cov)) fwrite(summary_y0_cov, file.path(dir_tables, paste0("summary_y0_cov_", reg_name, ".csv")))

    # B. Years 0-1 Sum
    summary_y01_cov <- make_summary_table(master_table_reg, "age01sum", "covars_on")
    if (!is.null(summary_y01_cov)) fwrite(summary_y01_cov, file.path(dir_tables, paste0("summary_y01_cov_", reg_name, ".csv")))

    # C. All Years Sum (includes birth-time-only outcomes from age0)
    summary_all_cov <- make_summary_table(master_table_reg, "age0_3sum", "covars_on")

    # C2. Replace birth-time-only outcomes with age0 results (not age0_3sum).
    # These outcomes (e.g., fp_father_share) are only meaningful at age 0.
    birth_time_only_summary <- make_summary_table(master_table_reg, "age0", "covars_on")
    if (!is.null(birth_time_only_summary) && !is.null(summary_all_cov)) {
      # Get birth-time-only outcomes from outcomes_config.R
      birth_time_outcomes <- if (exists("BIRTH_TIME_ONLY_OUTCOMES")) BIRTH_TIME_ONLY_OUTCOMES else c("fp_father_share")

      # Remove birth-time-only outcomes from age0_3sum results (they are artifacts)
      summary_all_cov <- summary_all_cov[!outcome %in% birth_time_outcomes]

      # Add the correct age0 results for birth-time-only outcomes
      birth_time_rows <- birth_time_only_summary[outcome %in% birth_time_outcomes]
      if (nrow(birth_time_rows) > 0) {
        summary_all_cov <- rbind(summary_all_cov, birth_time_rows, fill = TRUE)
        cat(sprintf("   Replaced %d birth-time-only outcomes with age0 results: %s\n",
                    nrow(birth_time_rows), paste(unique(birth_time_rows$outcome), collapse = ", ")))
      }

      # Also add any outcomes that are ONLY in age0 (missing from age0_3sum entirely)
      missing_outcomes <- setdiff(unique(birth_time_only_summary$outcome), unique(summary_all_cov$outcome))
      if (length(missing_outcomes) > 0) {
        extra_rows <- birth_time_only_summary[outcome %in% missing_outcomes]
        if (nrow(extra_rows) > 0) {
          summary_all_cov <- rbind(summary_all_cov, extra_rows, fill = TRUE)
          cat(sprintf("   Added %d additional age0-only outcomes to summary_all_cov: %s\n",
                      length(missing_outcomes), paste(missing_outcomes, collapse = ", ")))
        }
      }
    } else if (is.null(summary_all_cov) && !is.null(birth_time_only_summary)) {
      summary_all_cov <- birth_time_only_summary
    }
    if (!is.null(summary_all_cov)) fwrite(summary_all_cov, file.path(dir_tables, paste0("summary_all_cov_", reg_name, ".csv")))

    # D. FDR-adjusted summary table (post-treatment only, includes birth-time-only)
    fdr_summary <- make_fdr_summary_table(master_table_reg, "age0_3sum", "covars_on")

    # D2. Replace birth-time-only outcomes with age0 results
    fdr_birth_time <- make_fdr_summary_table(master_table_reg, "age0", "covars_on")
    if (!is.null(fdr_birth_time) && !is.null(fdr_summary)) {
      birth_time_outcomes <- if (exists("BIRTH_TIME_ONLY_OUTCOMES")) BIRTH_TIME_ONLY_OUTCOMES else c("fp_father_share")

      # Remove birth-time-only outcomes from age0_3sum results
      fdr_summary <- fdr_summary[!outcome %in% birth_time_outcomes]

      # Add the correct age0 results for birth-time-only outcomes
      fdr_birth_rows <- fdr_birth_time[outcome %in% birth_time_outcomes]
      if (nrow(fdr_birth_rows) > 0) {
        fdr_summary <- rbind(fdr_summary, fdr_birth_rows, fill = TRUE)
        cat(sprintf("   Replaced %d birth-time-only outcomes with age0 results: %s\n",
                    nrow(fdr_birth_rows), paste(unique(fdr_birth_rows$outcome), collapse = ", ")))
      }

      # Also add any outcomes that are ONLY in age0 (missing from age0_3sum entirely)
      missing_fdr <- setdiff(unique(fdr_birth_time$outcome), unique(fdr_summary$outcome))
      if (length(missing_fdr) > 0) {
        fdr_rows <- fdr_birth_time[outcome %in% missing_fdr]
        if (nrow(fdr_rows) > 0) {
          fdr_summary <- rbind(fdr_summary, fdr_rows, fill = TRUE)
          cat(sprintf("   Added %d additional age0-only outcomes to FDR summary: %s\n",
                      length(missing_fdr), paste(missing_fdr, collapse = ", ")))
        }
      }
    } else if (is.null(fdr_summary) && !is.null(fdr_birth_time)) {
      fdr_summary <- fdr_birth_time
    }
    if (!is.null(fdr_summary)) fwrite(fdr_summary, file.path(dir_tables, paste0("summary_FDR_", reg_name, ".csv")))

    # --- STEP 4: PAPER-READY TABLES (Stacked TWFE, With Covariates) ---
    paper_cumulative <- make_paper_table(
      master_table_reg,
      outcomes_to_show = base_outcomes,
      pretty_labels = pretty_labels,
      model_key = "stacked_TWFE_weighted",
      cov_filter = "covars_on",
      rel_key = "age0_3sum",
      event_times = c(-6, -4, -2, -1, 0, 1)
    )
    if (!is.null(paper_cumulative)) fwrite(paper_cumulative, file.path(dir_tables, paste0("paper_table_cumulative_", reg_name, ".csv")))

    cat(sprintf("   Successfully saved all CSV and Excel tables for %s.\n", reg_name))
  }
}

# --- STEP 5: HETEROGENEITY (Run once globally/National) ---
cat("\n--- Extracting National Heterogeneity (High vs Low CNI) ---\n")
if (exists("HET_RESULTS") && length(HET_RESULTS) > 0) {
  het_master <- extract_heterogeneity_results(HET_RESULTS)
  if (!is.null(het_master)) {
    fwrite(het_master, file.path(dir_tables, "heterogeneity_results_master.csv"))
    # Detail comparison
    cni_detailed <- make_cni_comparison_detailed(het_master, "age0_3sum", "covars_on")
    if (!is.null(cni_detailed)) fwrite(cni_detailed, file.path(dir_tables, "CNI_comparison_detailed_National.csv"))
  }
}

# ===========================================================================
# STEP 6: PRE-TRENDS FORMAL TESTING
# ===========================================================================
# Joint Wald test on pre-treatment coefficients for parallel trends assumption.
# Tests H0: All pre-treatment event-time coefficients = 0
# Reference: Roth (2022) "Pretest with Caution" AER: Insights
# ===========================================================================

cat("\n--- Running Formal Pre-Trends Tests (Joint Wald) ---\n")

# Run pre-trends test on National RESULTS (the primary analysis)
if (exists("RESULTS") && length(RESULTS) > 0) {

  # Test all base outcomes for which SUNAB models exist
  pre_trends_summary <- compile_pre_trends_summary(
    RESULTS = RESULTS,
    outcomes = base_outcomes,
    model_key = "SUNAB_model",
    rel_key = "age0",
    covars_on = TRUE,
    alpha = 0.10
  )

  if (!is.null(pre_trends_summary) && nrow(pre_trends_summary) > 0) {
    # Save to tables directory
    fwrite(pre_trends_summary, file.path(dir_tables, "pre_trends_test_results.csv"))

    # Print summary to console
    cat("\nPre-Trends Test Summary (SUNAB, age0, with covariates):\n")
    cat(rep("-", 60), "\n")

    n_tested <- sum(!is.na(pre_trends_summary$p_value))
    n_rejected <- sum(pre_trends_summary$reject_null == TRUE, na.rm = TRUE)

    cat(sprintf("  Outcomes tested: %d\n", n_tested))
    cat(sprintf("  Pre-trends rejected (p < 0.10): %d\n", n_rejected))

    if (n_rejected > 0) {
      concerns <- pre_trends_summary[reject_null == TRUE, outcome]
      cat("\n  ** OUTCOMES WITH PRE-TRENDS CONCERNS:\n")
      for (oc in concerns) {
        row <- pre_trends_summary[outcome == oc]
        cat(sprintf("     - %s (F=%.2f, p=%.3f)\n", oc, row$f_stat, row$p_value))
      }
      cat("\n  Consider additional robustness checks for these outcomes.\n")
      cat("  See Roth (2022) for discussion of pre-trends testing limitations.\n")
    } else {
      cat("\n  All outcomes pass pre-trends test at alpha = 0.10.\n")
    }
    cat(rep("-", 60), "\n")
  }

  # Also run without covariates for comparison
  pre_trends_nocov <- compile_pre_trends_summary(
    RESULTS = RESULTS,
    outcomes = base_outcomes,
    model_key = "SUNAB_model",
    rel_key = "age0",
    covars_on = FALSE,
    alpha = 0.10
  )

  if (!is.null(pre_trends_nocov) && nrow(pre_trends_nocov) > 0) {
    fwrite(pre_trends_nocov, file.path(dir_tables, "pre_trends_test_results_nocov.csv"))
  }

} else {
  cat("  Warning: RESULTS not available for pre-trends testing.\n")
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("TABLE EXTRACTION COMPLETE!\n")
cat(rep("=", 70), "\n")
