# ============================================================================
# OUTCOMES_CONFIG.R - Central Configuration for All Outcome Variables
# ============================================================================
# Paper III of PhD thesis: Effects of the Extended Home Visiting Programme
# on child health outcomes in Sweden.
#
# Description: Defines all outcome variables, their groupings, tier
#              classifications (primary/secondary/exploratory), directional
#              hypotheses, and helper functions for the effects analysis.
#
# Author:  Sergio Flores
# Date:    2026
# Paper:   Flores et al., "Effects of the Swedish Extended Home Visiting
#          Programme on child health outcomes: A register-based study"
#
# Required packages: tools (base R)
#
# CRITICAL TERMINOLOGY:
#   - event_time: Birth cohort relative to treatment (x-axis of event studies)
#   - rel_time / age: Child's age when outcome is measured (0 = first year of life)
#   - Reference period: event_time = -1 (last untreated birth cohort)
#
# OUTCOME SELECTION RATIONALE:
#   Final set reflects co-author consensus across multiple rounds of review.
#   Main analysis focuses on child health outcomes with clear programme theory
#   links and unambiguous directional hypotheses.
#   After final review: PRIMARY = 4 child health outcomes, SECONDARY = 1
#   (avoidable hospitalisations). FDR correction recalculates with fewer tests.
#
# OUTCOMES DROPPED during review:
#   - home_visits_this_year: fidelity measure, reported in separate paper
#   - fp_father_share: parent outcome, not child health
#   - parental_tfp_gross_days_year: parent outcome (VAB), not child health
#   - parent_had_dx_mental_health_common_year: parent outcome, not child health
#   - rapid_repeat_birth: no clear theoretical prediction
#   - mother_is_employed: too many confounders outside programme influence
#   - planned_care_share: demoted to exploratory (interpretation not straightforward)
#   - Diabetes from avoidable conditions: not amenable to prevention at this age
# ============================================================================

cat("\n--- Loading Outcome Configuration ---\n")

# ============================================================================
# SECTION 1: OUTCOME VARIABLE MAPS
# ============================================================================

# --- Child Healthcare Utilisation ---
ES_CHILD_HEALTH_UTIL_MAP <- c(
  "Avoidable Inpatient Events" = "inpatient_avoidable_events",
  "Avoidable Outpatient Visits" = "outpatient_avoidable_events",
  "Emergency Inpatient Events" = "inpatient_emergency_events",
  "Emergency Outpatient Visits" = "outpatient_emergency_events",
  "Total Inpatient Days (LOS)" = "inpatient_total_los_days"
)

# --- Child Injury (Safety Outcomes) ---
ES_CHILD_INJURY_MAP <- c(
  "Sentinel Injuries (Burns/Poisonings)" = "sentinel_injury_events",
  "Inpatient Burns Events" = "inpatient_trauma_burns_events",
  "Inpatient Poisoning Events" = "inpatient_trauma_poisoning_icd_events",
  "Inpatient Fall E-Code Events" = "inpatient_ecode_fall_events"
)

# --- Child Respiratory/Infection Outcomes ---
ES_CHILD_INFECTION_MAP <- c(
  "Inpatient Lower Respiratory Events" = "inpatient_respiratory_acute_lower_events",
  "Inpatient Gastroenteritis Events" = "inpatient_gastroenteritis_noninfectious_events",
  "Inpatient Asthma Events" = "inpatient_asthma_events",
  "Outpatient Asthma Visits" = "outpatient_asthma_events",
  "Combined Asthma Events" = "asthma_combined_events"
)

# --- BVC Process (Requires Smart Filter) ---
ES_BVC_PROCESS_MAP <- c(
  "BVC Home Visits This Year" = "home_visits_this_year",
  "Total Referrals This Year" = "total_referrals_this_year",
  "Speech Therapy Referrals" = "speech_referrals_this_year",
  "Psychologist Referrals" = "psych_referrals_this_year"
)

# --- Parental Health Diagnoses ---
# NOTE: Prescription data would be more sensitive than diagnoses for detecting
# mental health changes, but parental prescription data is not available in
# this dataset. Using diagnoses as proxy.
ES_PARENTAL_HEALTH_DX_MAP <- c(
  "Pr(Parent Dx: Common Mental Health)" = "parent_had_dx_mental_health_common_year",
  "Pr(Parent Dx: Stress/Burnout)" = "parent_had_dx_stress_burnout_year",
  "Pr(Parent Dx: Postpartum Depression)" = "parent_had_dx_postpartum_depression_year",
  "Pr(Parent Dx: Substance Use)" = "parent_had_dx_substance_use_year"
)

# --- Parental Benefits / Gender Equity ---
ES_PARENTAL_BENEFITS_MAP <- c(
  "Parental Sickness Days (SJP)" = "parental_sjp_net_days_year",
  "Care of Sick Child Days (VAB)" = "parental_tfp_gross_days_year",
  "Mother Parental Leave Days (FP)" = "mother_fp_net_days_year",
  "Father Parental Leave Days (FP)" = "father_fp_net_days_year",
  "Father Share of Parental Leave" = "fp_father_share"
)

# --- Child Oral Health ---
ES_ORAL_HEALTH_MAP <- c(
  "Dental Caries Score (DFT)" = "dft_score",
  "Dental Decay Score (DS)" = "ds_score",
  "Pr(Has Dental Exam)" = "has_dental_exam"
)

# --- Child Medication ---
ES_CHILD_MEDICATION_MAP <- c(
  "Pr(Any Systemic Antibiotic)" = "had_antibiotic_systemic_year",
  "Pr(Any Asthma/Respiratory Med)" = "had_asthma_respiratory_year",
  "Pr(Any Laxative)" = "had_laxative_year",
  "Systemic Antibiotic Prescriptions (count)" = "events_antibiotic_systemic_year",
  "Asthma/Respiratory Prescriptions (count)" = "events_asthma_respiratory_year"
)

# ============================================================================
# SECTION 2: DERIVED OUTCOMES (Created in 02_analysis.R)
# ============================================================================

# These variables are computed from raw data in the analysis script:
#
# sentinel_injury_events := inpatient_trauma_burns_events + inpatient_trauma_poisoning_icd_events
# count_emergency := inpatient_emergency_events + outpatient_emergency_events
# asthma_combined_events := inpatient_asthma_events + outpatient_asthma_events
# fp_father_share := father_fp_net_days_year / (mother_fp_net_days_year + father_fp_net_days_year)
#
# DROPPED (per co-author review):
# - planned_care_share (demoted to exploratory -- interpretation not straightforward)
# - rapid_repeat_birth (no clear hypothesis)
# - early_detection_events (mechanism unclear)

ES_DERIVED_MAP <- c(
  "Sentinel Injuries (Burns/Poisonings)" = "sentinel_injury_events",
  "Total Emergency Visit Count" = "count_emergency",
  "Combined Asthma Events (IP+OP)" = "asthma_combined_events",
  "Father Share of Parental Leave" = "fp_father_share"
)

# ============================================================================
# SECTION 3: SMART FILTER OUTCOMES
# ============================================================================
# These outcomes require filtering to DeSOs with valid BVC/registry reporting.

SMART_FILTER_OUTCOMES <- c(
  # BVC Process
  "home_visits_this_year",
  "total_referrals_this_year",
  "speech_referrals_this_year",
  "psych_referrals_this_year",

  # Dental (requires dental exam coverage) - NOTE: 0% for treatment cohorts
  "has_dental_exam",
  "dft_score",
  "ds_score"
)

# ============================================================================
# SECTION 3b: BIRTH-TIME-ONLY OUTCOMES
# ============================================================================
# These outcomes are only meaningful at child age 0 (birth year).
# mother/father_fp_net_days_year are 0% positive after age 1,
# so fp_father_share should only be analysed at event_time = 0.

BIRTH_TIME_ONLY_OUTCOMES <- c(
  # fp_father_share: DROPPED -- parent outcome, not child health
  "mother_fp_net_days_year",   # Mother FP days - 0% after age 1
  "father_fp_net_days_year"    # Father FP days - 0% after age 1
)

# ============================================================================
# SECTION 4: ACTIVE OUTCOME SETS - MAIN ANALYSIS
# ============================================================================
# Final set based on co-author consensus.
#
# DIRECTIONAL HYPOTHESES (must have clear predictions):
#   (up arrow) = Increase expected | (down arrow) = Decrease expected | ? = Ambiguous/exploratory

# -----------------------------------------------------------------------------
# PRIMARY OUTCOMES (Tier 1 - Pre-registered, confirmatory)
# -----------------------------------------------------------------------------
# These outcomes have:
#   - Clear programme theory link
#   - Unambiguous directional hypotheses
#   - Co-author approval

OUTCOMES_PRIMARY <- c(
  # DROPPED: home_visits_this_year (fidelity measure, reported in separate paper)

  # 1. SAFETY - Core prevention outcome
  "sentinel_injury_events",           # Decrease: Reduction in burns + poisonings

  # 2. HEALTHCARE UTILISATION
  "count_emergency",                  # Decrease: Reduction in emergency visits

  # 3. INFECTIONS - SES gradient outcome
 "inpatient_respiratory_acute_lower_events",  # Decrease: Reduction (unequally distributed by SES)

  # 4. ASTHMA - Avoidable condition
  "asthma_combined_events"            # Decrease: Reduction expected
)

# -----------------------------------------------------------------------------
# SECONDARY OUTCOMES (Tier 2 - Exploratory, hypothesis-generating)
# -----------------------------------------------------------------------------
# These outcomes have programme theory links but may have:
#   - Lower power (rare events)
#   - More complex interpretation
#   - Emerged from heterogeneity analysis

OUTCOMES_SECONDARY <- c(
  # 5. AVOIDABLE HOSPITALISATIONS (excluding diabetes -- not amenable to prevention)
  "inpatient_avoidable_events"        # Decrease: Reduction in preventable admissions

  # DROPPED: fp_father_share (parent outcome, not child health)
  # DROPPED: parental_tfp_gross_days_year (parent outcome, VAB)
  # DROPPED: parent_had_dx_mental_health_common_year (parent outcome)
)

# -----------------------------------------------------------------------------
# EXPLORATORY OUTCOMES (Tier 3 - Robustness checks only)
# -----------------------------------------------------------------------------
# Not part of main analysis but available for sensitivity checks

OUTCOMES_EXPLORATORY <- c(
  "planned_care_share",               # Demoted per co-author feedback
  "inpatient_gastroenteritis_noninfectious_events",
  "parental_sjp_net_days_year",

  # BVC Referrals - Require Smart Filter
  "total_referrals_this_year",        # 81% coverage - moderate power
  "speech_referrals_this_year",       # 46% coverage - LOW POWER warning

  # Parental mental health medications
  # From LMED: ATC N06A (antidepressants) and N05B/N05C (anxiolytics/sedatives)
  "parent_any_antidepressant_year",   # ? Better detection vs true change
  "parent_any_anxiolytic_year"        # ? Better detection vs true change
)

# -----------------------------------------------------------------------------
# DROPPED OUTCOMES (Do NOT analyse)
# -----------------------------------------------------------------------------
#
# DROPPED PER CO-AUTHOR CONSENSUS:
# rapid_repeat_birth      - No clear theoretical prediction
# mother_is_employed      - Too many confounders outside programme influence
# had_laxative_year       - 0% prevalence in data
# Diabetes (from avoidable) - Not amenable to prevention at this age
#
# -----------------------------------------------------------------------------
# DROPPED - PARENT OUTCOMES (final review round)
# -----------------------------------------------------------------------------
# Dropped from main analysis to focus on child health outcomes:
# home_visits_this_year              - Fidelity measure, reported in separate paper
# fp_father_share                    - Parent outcome, not child health
# parental_tfp_gross_days_year       - Parent outcome (VAB), not child health
# parent_had_dx_mental_health_common_year - Parent outcome, not child health
#
# After removal: PRIMARY = 4 outcomes, SECONDARY = 1 outcome
# FDR correction recalculates automatically with fewer tests
#
# -----------------------------------------------------------------------------
# DROPPED - VACCINATION OUTCOMES (data audit)
# -----------------------------------------------------------------------------
# ALL vaccination outcomes dropped due to Svevac registry limitations.
# Registry covers only ~17% of children; near-zero variation in matched sample.
#
# dtap_3dos_utd_by_18m     - 39.8% coverage, 1.4% positive, methodological mismatch
# mmr1_utd_by_24m          - 35.7% coverage, 2.4% positive
# mmr2_utd_by_7yr          - 6.7% coverage, 0% positive (no variation)
# dtap_booster4_utd_by_6yr - 10.6% coverage, 0% positive (no variation)
# total_doses_mmr_by_year  - 100% non-missing but only 2.9% positive (97% zeros)
# total_doses_dtap_combo_by_year - 100% non-missing but only 3.0% positive
# has_mmr_by_year          - 100% non-missing, 0% positive (no variation)
# has_rota_by_year         - 100% non-missing, 0% positive (no variation)
# has_pneumo_by_year       - 100% non-missing, 0% positive (no variation)
#
# -----------------------------------------------------------------------------
# DROPPED - SMOKING OUTCOME (data audit)
# -----------------------------------------------------------------------------
# smoking_ever_home        - DIFFERENTIAL REPORTING ARTIFACT
#                            100% technical coverage BUT 31x differential:
#                            Treated 0.24% vs Control 7.49% (-7.25 ppt)
#                            Within-region analysis shows NO differential;
#                            artifact from national control pool composition.
#                            Cannot use as outcome OR covariate.
#
# -----------------------------------------------------------------------------
# DROPPED - BREASTFEEDING OUTCOMES (data audit)
# -----------------------------------------------------------------------------
# bf_any_duration_days     - 5.2% coverage, 33x differential (T: 0.2% vs C: 6.7%)
# bf_exclusive_duration_days - Same as above; severe selection bias risk
# bf_status_initial        - 100% non-missing but only 4.4% positive, -5.4 ppt diff
# bf_status_latest         - 100% non-missing but only 1.6% positive, -2.0 ppt diff
# bf_max_status            - 100% non-missing but only 4.5% positive, -5.5 ppt diff
#
# -----------------------------------------------------------------------------
# DROPPED - DENTAL OUTCOMES (data limitation)
# -----------------------------------------------------------------------------
# dft_score                - SKaPa only has births 2012-2015; 0% for 2018-2020
# ds_score                 - Same data limitation
# has_dental_exam          - Same; data linked via personal identification number
#                            with no crosswalk available for study cohorts

# ============================================================================
# SECTION 5: COMBINED SETS AND LABELS
# ============================================================================

# Combined base outcomes for main analysis
BASE_OUTCOMES <- c(OUTCOMES_PRIMARY, OUTCOMES_SECONDARY)

# Pretty labels for plotting (order matches BASE_OUTCOMES)
PRETTY_LABELS_MAP <- c(
  # Primary (4 child health outcomes)
  "sentinel_injury_events"                    = "Sentinel Injuries (Burns/Poisonings)",
  "count_emergency"                           = "Emergency Visits (IP+OP)",
  "inpatient_respiratory_acute_lower_events"  = "Respiratory Infection Hospitalisations",
  "asthma_combined_events"                    = "Asthma Events (IP+OP)",

  # Secondary (1 child health outcome)
  "inpatient_avoidable_events"                = "Avoidable Hospitalisations"
)

# Extract labels in order of BASE_OUTCOMES
PRETTY_LABELS <- PRETTY_LABELS_MAP[BASE_OUTCOMES]

# ============================================================================
# SECTION 6: DIRECTIONAL HYPOTHESES (For interpretation)
# ============================================================================
# Used in results tables and discussion

HYPOTHESIS_DIRECTION <- c(
  # Primary - all have clear predictions
  "sentinel_injury_events"                    = "decrease",
  "count_emergency"                           = "decrease",
  "inpatient_respiratory_acute_lower_events"  = "decrease",
  "asthma_combined_events"                    = "decrease",

  # Secondary
  "inpatient_avoidable_events"                = "decrease"
)

# ============================================================================
# SECTION 7: ALL OUTCOME VARIABLES (For data subsetting / collapse)
# ============================================================================

ALL_OUTCOME_VARS <- c(
  # --- DERIVED OUTCOMES (Created in 02_analysis.R) ---
  "sentinel_injury_events",
  "count_emergency",
  "asthma_combined_events",
  "fp_father_share",
  "planned_care_share",  # Keep for exploratory

  # --- RAW OUTCOMES: Healthcare Utilisation ---
  "inpatient_avoidable_events",
  "outpatient_avoidable_events",
  "inpatient_emergency_events",
  "outpatient_emergency_events",
  "inpatient_total_los_days",

  # --- RAW OUTCOMES: Injuries ---
  "inpatient_trauma_burns_events",
  "inpatient_trauma_poisoning_icd_events",
  "inpatient_ecode_fall_events",

  # --- RAW OUTCOMES: Infections/Respiratory ---
  "inpatient_respiratory_acute_lower_events",
  "inpatient_gastroenteritis_noninfectious_events",
  "inpatient_asthma_events",
  "outpatient_asthma_events",

  # --- RAW OUTCOMES: Parental Health ---
  "parent_had_dx_mental_health_common_year",
  "parent_had_dx_stress_burnout_year",
  "parent_had_dx_postpartum_depression_year",

  # --- RAW OUTCOMES: Parental Benefits ---
  "parental_sjp_net_days_year",
  "parental_tfp_gross_days_year",
  "mother_fp_net_days_year",
  "father_fp_net_days_year",

  # --- RAW OUTCOMES: Oral Health ---
  "dft_score",
  "ds_score",
  "has_dental_exam",

  # --- RAW OUTCOMES: BVC Process ---
  "home_visits_this_year",
  "total_referrals_this_year",
  "speech_referrals_this_year",
  "psych_referrals_this_year",

  # --- RAW OUTCOMES: Vaccination ---
  # dtap_3dos_utd_by_18m: DROPPED -- see DROPPED section
  "total_doses_dtap_combo_by_year",
  "total_doses_mmr_by_year",

  # --- RAW OUTCOMES: Child Medication ---
  "had_antibiotic_systemic_year",
  "had_asthma_respiratory_year",
  "events_antibiotic_systemic_year",
  "events_asthma_respiratory_year",

  # --- RAW OUTCOMES: Parental Medication ---
  "parent_any_antidepressant_year",
  "parent_any_anxiolytic_year"
)

# ============================================================================
# SECTION 8: KEEP LIST (For data subsetting in 02_analysis.R)
# ============================================================================

KEEP_VARS <- list(
  IDs = c(
    "lopnr",
    "deso_2021",
    "year",
    "child_ordnrmor"
  ),

  Covariates = c(
    "cni_score",
    "foreign_born",
    "single_parent",
    "low_education",
    "unemployed",
    "mother_birth_year",
    "low_birth_weight"
  ),

  Model_and_Cohort_Variables = c(
    "child_birth_year",
    "first_treat_year",
    "analysis_group",
    "expected_deso_treatment_dosage"
  ),

  # For gender equity calculation
  Parental_Leave = c(
    "mother_fp_net_days_year",
    "father_fp_net_days_year"
  ),

  Outcome_Variables = ALL_OUTCOME_VARS
)

# Flatten for easy use
KEEP_VARS_FLAT <- unique(unlist(KEEP_VARS))

# ============================================================================
# SECTION 9: HELPER FUNCTIONS
# ============================================================================

#' Check if an outcome requires Smart Filter
needs_smart_filter <- function(var) {
  var %in% SMART_FILTER_OUTCOMES
}

#' Get label for a variable name
get_label <- function(var) {
  if (var %in% names(PRETTY_LABELS_MAP)) {
    return(PRETTY_LABELS_MAP[var])
  }
  # Fallback: clean up variable name
  gsub("_", " ", tools::toTitleCase(var))
}

#' Get hypothesis direction for a variable
get_hypothesis <- function(var) {
  if (var %in% names(HYPOTHESIS_DIRECTION)) {
    return(HYPOTHESIS_DIRECTION[var])
  }
  return("unknown")
}

#' Check if outcome is primary
is_primary <- function(var) {
  var %in% OUTCOMES_PRIMARY
}

#' Check if outcome is secondary
is_secondary <- function(var) {
  var %in% OUTCOMES_SECONDARY
}

#' Check if outcome is birth-time only
#' These outcomes should only be analysed at event_time = 0
is_birth_time_only <- function(var) {
  var %in% BIRTH_TIME_ONLY_OUTCOMES
}

# ============================================================================
# SECTION 10: PRINT SUMMARY
# ============================================================================

cat("\n  === OUTCOME CONFIGURATION SUMMARY ===\n")
cat(sprintf("  PRIMARY outcomes (confirmatory): %d\n", length(OUTCOMES_PRIMARY)))
for (v in OUTCOMES_PRIMARY) {
  cat(sprintf("    - %s [%s]\n", get_label(v), get_hypothesis(v)))
}
cat(sprintf("\n  SECONDARY outcomes (exploratory): %d\n", length(OUTCOMES_SECONDARY)))
for (v in OUTCOMES_SECONDARY) {
  cat(sprintf("    - %s [%s]\n", get_label(v), get_hypothesis(v)))
}
cat(sprintf("\n  EXPLORATORY outcomes (robustness): %d\n", length(OUTCOMES_EXPLORATORY)))
cat(sprintf("  Smart Filter outcomes: %d\n", length(SMART_FILTER_OUTCOMES)))
cat(sprintf("  Birth-Time-Only outcomes: %d\n", length(BIRTH_TIME_ONLY_OUTCOMES)))
cat(sprintf("  Total outcome variables tracked: %d\n", length(ALL_OUTCOME_VARS)))
cat("\n  DROPPED per co-author consensus:\n")
cat("    - rapid_repeat_birth (no clear hypothesis)\n")
cat("    - mother_is_employed (too many confounders)\n")
cat("    - Diabetes in avoidable (not amenable to prevention)\n")
cat("\n  DROPPED (final review, focus on child health):\n")
cat("    - home_visits_this_year (fidelity, separate paper)\n")
cat("    - fp_father_share (parent outcome)\n")
cat("    - parental_tfp_gross_days_year (parent outcome)\n")
cat("    - parent_had_dx_mental_health_common_year (parent outcome)\n")
cat("\n  DROPPED due to DATA LIMITATION:\n")
cat("    - dft_score, ds_score, has_dental_exam\n")
cat("      SKaPa only has births 2012-2015; linked via pseudonymised identifier (no crosswalk)\n")
cat("\n")
