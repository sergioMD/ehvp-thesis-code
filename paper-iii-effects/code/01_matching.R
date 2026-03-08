# ============================================================================
# 01_matching.R
# Propensity score matching for the EHVP effects analysis
# ============================================================================
#
# Paper:   "Effects of the Extended Home Visiting Programme on Child
#           Health Outcomes: A Staggered Difference-in-Differences Analysis"
#           (Paper III in doctoral thesis)
# Author:  Sergio Flores
# Date:    2025
#
# Description:
#   Loads DeSO-level panel data, classifies neighbourhoods into treated
#   and control groups based on programme dosage, and performs stratified
#   nearest-neighbour propensity score matching by treatment cohort.
#   Produces balance diagnostics (tables, love plots, maps) and saves
#   matched objects for use by 02_analysis.R.
#
# Matching specification:
#   - Method: nearest-neighbour propensity score
#   - Exact matching on urbanicity (Urban / Suburban / Rural)
#   - No caliper (maximises sample retention)
#   - With replacement
#   - Ratio: up to 5 controls per treated unit
#   - Regional confounding controlled via fixed effects in the outcome model
#
# Dependencies:
#   data.table, ggplot2, MatchIt, cobalt, dplyr, stringr,
#   kableExtra, sf, patchwork
#
# Usage:
#   Source config.R first (sets paths and parameters), then run this script.
#   Expected working directory: project root.
#
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("01_MATCHING.R - EHVP Optimal Matching Analysis\n")
cat("============================================================================\n")
cat(sprintf("Run started: %s\n", Sys.time()))
cat("============================================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Source configuration (sets paths, thresholds, covariate lists)
source(here::here("code", "config.R"))

# Validate configuration
validate_config()

# Create output directories
create_matching_dirs()

# Load packages
cat("\n--- Loading packages ---\n")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

pacman::p_load(
  data.table, ggplot2, MatchIt, cobalt, dplyr, stringr,
  kableExtra, sf, patchwork,
  character.only = FALSE
)

# ============================================================================
# PHASE 1: DATA PREPARATION
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("PHASE 1: DATA PREPARATION\n")
cat("============================================================================\n\n")

# --- Load raw panel data ---
cat("Loading raw panel data...\n")
cat(sprintf("  Source: %s\n", RAW_PANEL_DATA_PATH))

full_panel_data <- fread(RAW_PANEL_DATA_PATH)
cat(sprintf("  Loaded: %s rows x %s columns\n",
            format(nrow(full_panel_data), big.mark = ","),
            ncol(full_panel_data)))

# --- Extract birth DeSO for each child ---
cat("\nExtracting birth DeSO mapping...\n")

deso_cols <- grep("^deso_\\d{4}$", names(full_panel_data), value = TRUE)
deso_long <- melt(
  full_panel_data[, c("lopnr", "child_birth_year", ..deso_cols)],
  id.vars = c("lopnr", "child_birth_year"),
  measure.vars = deso_cols,
  value.name = "birth_deso_code"
)
deso_long[, deso_year := as.integer(gsub("deso_", "", variable))]
birth_deso_map <- unique(deso_long[child_birth_year == deso_year, .(lopnr, birth_deso_code)])

analysis_data <- merge(full_panel_data, birth_deso_map, by = "lopnr")

# --- Aggregate to DeSO-birth year level ---
cat("Aggregating to DeSO-birth year panel...\n")

first_year_data <- analysis_data[
  year == child_birth_year &
    !is.na(birth_deso_code) &
    birth_deso_code != "9999"
]

deso_birth_year_panel <- first_year_data[, .(
  num_children = .N,
  prop_low_birth_weight = mean(low_birth_weight, na.rm = TRUE),
  avg_mother_age_at_birth = mean(child_birth_year - mother_birth_year, na.rm = TRUE),
  prop_low_education = mean(low_education, na.rm = TRUE),
  prop_foreign_born = mean(foreign_born, na.rm = TRUE),
  prop_unemployed = mean(unemployed, na.rm = TRUE),
  prop_single_parent = mean(single_parent, na.rm = TRUE)
), by = .(deso_id = birth_deso_code, year = child_birth_year)]

# --- Add CNI and treatment dosage ---
cat("Adding CNI and treatment dosage...\n")

cni_dose_vars <- unique(full_panel_data[, .(
  deso_id = deso_2021,
  year,
  cni_score,
  expected_deso_treatment_dosage
)])

deso_birth_year_panel <- merge(
  deso_birth_year_panel,
  cni_dose_vars,
  by = c("deso_id", "year"),
  all.x = TRUE
)

# Clean and add geographic identifiers
deso_birth_year_panel <- deso_birth_year_panel[deso_id != "9999" & !is.na(deso_id)]

deso_birth_year_panel[, `:=`(
  urban_rural_category = fcase(
    substr(deso_id, 5, 5) == 'A', "Rural",
    substr(deso_id, 5, 5) == 'B', "Suburban",
    substr(deso_id, 5, 5) == 'C', "Urban",
    default = "Unknown"
  ),
  lan_code = substr(deso_id, 1, 2)
)]

setkey(deso_birth_year_panel, deso_id, year)

cat(sprintf("  Created panel: %s DeSO-years across %s unique DeSOs\n",
            format(nrow(deso_birth_year_panel), big.mark = ","),
            format(uniqueN(deso_birth_year_panel$deso_id), big.mark = ",")))

# Save intermediate data
fwrite(deso_birth_year_panel,
       file.path(MATCHING_OUTPUT_DIR, "deso_birth_year_panel.csv"))

# Clean up memory
rm(full_panel_data, analysis_data, first_year_data, deso_long, birth_deso_map, cni_dose_vars)
gc()

cat("\nPHASE 1 COMPLETE\n")

# ============================================================================
# PHASE 2: DeSO CLASSIFICATION
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("PHASE 2: DeSO CLASSIFICATION\n")
cat("============================================================================\n\n")

# Calculate maximum dosage per DeSO
deso_max_dose <- deso_birth_year_panel[, .(
  max_dose = max(expected_deso_treatment_dosage, na.rm = TRUE)
), by = deso_id]
deso_max_dose[is.infinite(max_dose), max_dose := NA]

# Calculate pre-treatment CNI (average up to 2020)
deso_avg_cni <- deso_birth_year_panel[
  year <= 2020,
  .(avg_cni_pre_treatment = mean(cni_score, na.rm = TRUE)),
  by = deso_id
]

# Merge and classify
deso_summary <- merge(deso_max_dose, deso_avg_cni, by = "deso_id", all = TRUE)

deso_summary[, analysis_group := fcase(
  max_dose >= TREATED_THRESHOLD, "Treated",
  max_dose <= CONTROL_THRESHOLD, "Control",
  default = "Fuzzy Zone"
)]

# Print classification summary
group_counts <- deso_summary[, .N, by = analysis_group]
cat("DeSO Classification:\n")
for (i in 1:nrow(group_counts)) {
  cat(sprintf("  %s: %d DeSOs\n", group_counts$analysis_group[i], group_counts$N[i]))
}

# Save classification
fwrite(deso_summary, file.path(MATCHING_OUTPUT_DIR, "deso_classification.csv"))

# --- Create exploratory plots ---
cat("\nCreating exploratory plots...\n")

# Dosage distribution plot
dosage_dist_plot <- ggplot(deso_summary[!is.na(max_dose)], aes(x = max_dose)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  geom_vline(xintercept = CONTROL_THRESHOLD, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = TREATED_THRESHOLD, color = "darkred", linetype = "dashed", linewidth = 1) +
  annotate("text", x = CONTROL_THRESHOLD - 0.02, y = Inf,
           label = sprintf("Control (<= %.0f%%)", CONTROL_THRESHOLD * 100),
           hjust = 1, vjust = 2, color = "darkblue", size = 3) +
  annotate("text", x = TREATED_THRESHOLD + 0.02, y = Inf,
           label = sprintf("Treated (>= %.0f%%)", TREATED_THRESHOLD * 100),
           hjust = 0, vjust = 2, color = "darkred", size = 3) +
  labs(
    title = "Distribution of Maximum Treatment Dosage Across All DeSOs",
    subtitle = "Vertical lines indicate thresholds for group assignment",
    x = "Maximum Treatment Dosage",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(MATCHING_OUTPUT_DIR, "dosage_distribution.png"),
       dosage_dist_plot, width = 10, height = 6, dpi = 300)

cat("\nPHASE 2 COMPLETE\n")

# ============================================================================
# PHASE 3: OPTIMAL MATCHING
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("PHASE 3: OPTIMAL MATCHING\n")
cat("============================================================================\n")
cat("\nSPECIFICATION:\n")
cat("  - Method: Nearest-neighbor propensity score matching\n")
cat("  - Exact matching: Urbanicity only (Urban/Suburban/Rural)\n")
cat("  - Caliper: NONE (maximizes sample retention)\n")
cat("  - Replacement: YES\n")
cat("  - Ratio: Up to 5 controls per treated\n")
cat("  - Regional confounding: Controlled via FE in outcome model\n")
cat(sprintf("\nCohorts: %s\n", paste(COHORTS_TO_MATCH, collapse = ", ")))
cat("============================================================================\n\n")

# --- Prepare matching panel ---
deso_panel <- merge(
  deso_birth_year_panel,
  deso_summary[, .(deso_id, analysis_group)],
  by = "deso_id",
  all.x = TRUE
)

# Calculate first treatment year for each DeSO
first_treat_year_map <- deso_panel[, .(
  first_treat_year = {
    valid_years <- year[expected_deso_treatment_dosage >= TREATED_THRESHOLD]
    if (length(valid_years) == 0) NA_real_ else as.numeric(min(valid_years))
  }
), by = deso_id]

deso_panel <- merge(deso_panel, first_treat_year_map, by = "deso_id", all.x = TRUE)
deso_panel[analysis_group == "Control", first_treat_year := Inf]

# --- Matching formula ---
matching_formula <- as.formula(
  paste("treated_binary ~", paste(MATCHING_COVARIATES, collapse = " + "))
)

# --- Main matching loop ---
# Named 'match_results_list_stratified' for compatibility with 02_analysis.R
match_results_list_stratified <- list()
all_balance_results <- list()

for (cohort in COHORTS_TO_MATCH) {

  cat(sprintf("\n--- Cohort %d ---\n", cohort))

  baseline_year <- cohort - 1

  # Identify treated DeSOs for this cohort
  treated_desos <- deso_panel[first_treat_year == cohort, unique(deso_id)]
  untreated_desos <- deso_panel[analysis_group == "Control", unique(deso_id)]

  cat(sprintf("  Treated DeSOs in this cohort: %d\n", length(treated_desos)))
  cat(sprintf("  Potential untreated DeSOs: %d\n", length(untreated_desos)))

  # Create matching dataset at baseline
  matching_data <- deso_panel[
    year == baseline_year &
      (deso_id %in% treated_desos | deso_id %in% untreated_desos)
  ]
  matching_data[, treated_binary := fifelse(deso_id %in% treated_desos, 1L, 0L)]

  # Impute missing covariate values with column means
  for (covar in MATCHING_COVARIATES) {
    if (any(!is.finite(matching_data[[covar]]))) {
      impute_val <- mean(matching_data[[covar]], na.rm = TRUE)
      matching_data[!is.finite(get(covar)), (covar) := impute_val]
    }
  }

  # Impute CNI if needed (for diagnostics)
  if (any(!is.finite(matching_data$cni_score))) {
    matching_data[!is.finite(cni_score), cni_score := mean(matching_data$cni_score, na.rm = TRUE)]
  }

  # Urbanicity distribution
  urban_dist <- matching_data[, .(
    Treated = sum(treated_binary == 1),
    Control = sum(treated_binary == 0)
  ), by = urban_rural_category]
  cat("  Urbanicity distribution:\n")
  print(urban_dist)

  # --- Perform matching ---
  tryCatch({
    match_obj <- matchit(
      matching_formula,
      data = matching_data,
      method = MATCHING_METHOD,
      exact = as.formula(paste("~", MATCHING_EXACT)),
      distance = "glm",
      ratio = MATCHING_RATIO,
      replace = MATCHING_REPLACE
      # No caliper -- key design choice to maximise sample retention
    )

    md <- match.data(match_obj)
    n_treated <- sum(md$treated_binary == 1)
    n_control <- sum(md$treated_binary == 0)

    cat(sprintf("  Matched: %d treated -> %d controls\n", n_treated, n_control))

    # Store results
    match_results_list_stratified[[as.character(cohort)]] <- list(
      match_object = match_obj,
      matched_data = as.data.table(md),
      cohort_year = cohort,
      baseline_year = baseline_year,
      n_treated = n_treated,
      n_controls = n_control
    )

    # --- Balance diagnostics ---

    # Balance table
    tryCatch({
      balance_stats <- bal.tab(
        match_obj,
        stats = c("m", "v"),
        thresholds = c(m = SMD_EXCELLENT),
        un = TRUE,
        addl = ~ cni_score
      )

      if (!is.null(balance_stats$Balance)) {
        balance_df <- as.data.table(balance_stats$Balance, keep.rownames = "Variable")
        balance_df[, Cohort := cohort]

        # Rename columns for clarity
        if ("Diff.Adj" %in% names(balance_df)) setnames(balance_df, "Diff.Adj", "SMD_Matched")
        if ("Diff.Un" %in% names(balance_df)) setnames(balance_df, "Diff.Un", "SMD_Unmatched")
        if ("M.0.Adj" %in% names(balance_df)) setnames(balance_df, "M.0.Adj", "Mean_MatchedControl")
        if ("M.1.Adj" %in% names(balance_df)) setnames(balance_df, "M.1.Adj", "Mean_Treated")
        if ("M.0.Un" %in% names(balance_df)) setnames(balance_df, "M.0.Un", "Mean_Unmatched_Control")
        if ("M.1.Un" %in% names(balance_df)) setnames(balance_df, "M.1.Un", "Mean_Unmatched_Treated")

        all_balance_results[[as.character(cohort)]] <- balance_df

        fwrite(balance_df,
               file.path(MATCHING_OUTPUT_DIR, "balance_tables",
                         sprintf("balance_cohort_%d.csv", cohort)))

        balance_df %>%
          kable("html", caption = sprintf("Covariate Balance: Cohort %d (Treated vs Matched Controls)", cohort), digits = 3) %>%
          kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
          save_kable(file = file.path(MATCHING_OUTPUT_DIR, "balance_tables",
                                      sprintf("balance_cohort_%d.html", cohort)))
      }
    }, error = function(e) {
      cat(sprintf("  Warning: Balance table error: %s\n", e$message))
    })

    # Love plot
    tryCatch({
      love_plot <- love.plot(
        match_obj,
        stats = "mean.diffs",
        binary = "std",
        abs = TRUE,
        var.order = "alphabetical",  # Fixed order for consistency across cohorts
        thresholds = c(m = SMD_EXCELLENT),
        colors = c("red", "blue"),
        shapes = c("circle filled", "triangle filled"),
        sample.names = c("Unmatched", "Matched"),
        addl = ~ cni_score,
        title = NULL
      ) +
        labs(
          title = sprintf("Covariate Balance: Cohort %d", cohort),
          subtitle = "Urbanicity Only, No Caliper, With Replacement | Absolute SMD",
          x = "Absolute Standardized Mean Difference"
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")

      ggsave(file.path(MATCHING_OUTPUT_DIR, "love_plots",
                       sprintf("love_plot_cohort_%d.png", cohort)),
             love_plot, width = 10, height = 7, dpi = 300)
    }, error = function(e) {
      cat(sprintf("  Warning: Love plot error: %s\n", e$message))
    })

    # Propensity score overlap plot
    tryCatch({
      ps_plot <- bal.plot(match_obj, var.name = "distance", which = "both",
                          type = "density", mirror = TRUE,
                          sample.names = c("Unmatched", "Matched")) +
        labs(title = sprintf("Propensity Score Overlap: Cohort %d", cohort),
             fill = "Sample") +
        theme_minimal()

      ggsave(file.path(MATCHING_OUTPUT_DIR, "love_plots",
                       sprintf("ps_overlap_cohort_%d.png", cohort)),
             ps_plot, width = 10, height = 6, dpi = 300)
    }, error = function(e) {
      cat(sprintf("  Warning: PS overlap plot error: %s\n", e$message))
    })

  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", e$message))
  })
}

# ============================================================================
# COMPILE RESULTS
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("COMPILING RESULTS\n")
cat("============================================================================\n\n")

# Compile matched DeSO IDs
matched_deso_ids <- unique(unlist(lapply(match_results_list_stratified, function(x) {
  if (!is.null(x$matched_data)) x$matched_data$deso_id else character(0)
})))

matched_treated_ids <- unique(unlist(lapply(match_results_list_stratified, function(x) {
  if (!is.null(x$matched_data)) {
    x$matched_data[treated_binary == 1, deso_id]
  } else character(0)
})))

matched_control_ids <- unique(unlist(lapply(match_results_list_stratified, function(x) {
  if (!is.null(x$matched_data)) {
    x$matched_data[treated_binary == 0, deso_id]
  } else character(0)
})))

final_matched_panel <- deso_panel[deso_id %in% matched_deso_ids]

cat(sprintf("Total matched sample: %d treated + %d controls = %d unique DeSOs\n",
            length(matched_treated_ids),
            length(matched_control_ids),
            length(matched_deso_ids)))

# --- Matching summary table ---
matching_summary <- rbindlist(lapply(names(match_results_list_stratified), function(cohort_str) {
  x <- match_results_list_stratified[[cohort_str]]
  data.table(
    Cohort = as.integer(cohort_str),
    Specification = "Urbanicity Only, No Caliper",
    N_Treated = x$n_treated,
    N_Matched_Controls = x$n_controls
  )
}))

fwrite(matching_summary, file.path(MATCHING_OUTPUT_DIR, "matching_summary.csv"))

matching_summary %>%
  kable("html", caption = "Matching Summary by Cohort") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
  save_kable(file = file.path(MATCHING_OUTPUT_DIR, "matching_summary.html"))

cat("\nMatching Summary:\n")
print(matching_summary)

# --- Combined love plot ---
cat("\nCreating combined balance plot...\n")

# Pretty labels for variables
var_labels <- c(
  "prop_low_education" = "Prop. Low Education",
  "prop_foreign_born" = "Prop. Foreign-Born",
  "prop_unemployed" = "Prop. Unemployed",
  "prop_single_parent" = "Prop. Single Parent",
  "prop_low_birth_weight" = "Prop. Low Birth Weight",
  "avg_mother_age_at_birth" = "Avg. Mother Age",
  "cni_score" = "CNI Score"
)

overall_balance <- rbindlist(lapply(names(match_results_list_stratified), function(cohort_str) {
  x <- match_results_list_stratified[[cohort_str]]
  md <- x$matched_data

  covars <- c(MATCHING_COVARIATES, "cni_score")

  rbindlist(lapply(covars, function(var) {
    if (!var %in% names(md)) return(NULL)
    mean_t <- md[treated_binary == 1, mean(get(var), na.rm = TRUE)]
    mean_c <- md[treated_binary == 0, mean(get(var), na.rm = TRUE)]
    var_t <- md[treated_binary == 1, var(get(var), na.rm = TRUE)]
    var_c <- md[treated_binary == 0, var(get(var), na.rm = TRUE)]
    sd_pool <- sqrt((var_t + var_c) / 2)
    data.table(Cohort = cohort_str, Variable = var, Abs_SMD = abs((mean_t - mean_c) / sd_pool))
  }))
}))

# Apply pretty labels
overall_balance[, Variable := fifelse(Variable %in% names(var_labels), var_labels[Variable], Variable)]

combined_plot <- ggplot(overall_balance, aes(x = Abs_SMD, y = Variable, color = Cohort, shape = Cohort)) +
  geom_vline(xintercept = SMD_EXCELLENT, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = SMD_ACCEPTABLE, linetype = "dotted", color = "grey70") +
  geom_point(size = 3.5, position = position_dodge(width = 0.4)) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Covariate Balance: Treated vs Matched Controls (All Cohorts)",
    subtitle = sprintf("Optimal Specification | Dashed = %.2f (excellent), Dotted = %.2f (acceptable)",
                       SMD_EXCELLENT, SMD_ACCEPTABLE),
    x = "Absolute Standardized Mean Difference",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(MATCHING_OUTPUT_DIR, "love_plot_all_cohorts.png"),
       combined_plot, width = 10, height = 7, dpi = 300)

# --- Balance summary ---
balance_summary <- overall_balance[, .(
  Max_SMD = max(Abs_SMD, na.rm = TRUE),
  Mean_SMD = mean(Abs_SMD, na.rm = TRUE),
  Vars_Above_0.1 = sum(Abs_SMD > SMD_EXCELLENT),
  Vars_Above_0.25 = sum(Abs_SMD > SMD_ACCEPTABLE)
), by = Cohort]

fwrite(balance_summary, file.path(MATCHING_OUTPUT_DIR, "balance_summary.csv"))

cat("\nBalance Summary (Max |SMD| per cohort):\n")
print(balance_summary)

# ============================================================================
# DIAGNOSTIC VISUALIZATIONS
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("DIAGNOSTIC VISUALIZATIONS\n")
cat("============================================================================\n\n")

# Create diagnostic plots directory
dir.create(file.path(MATCHING_OUTPUT_DIR, "diagnostics"), showWarnings = FALSE)
diag_dir <- file.path(MATCHING_OUTPUT_DIR, "diagnostics")

# -----------------------------------------------------------------------------
# 1. CNI density: treated vs non-treated (before matching)
# -----------------------------------------------------------------------------
cat("--- Creating CNI density plots ---\n")

# Pre-matching comparison
cni_pre_match <- deso_summary[analysis_group %in% c("Treated", "Control"),
                              .(deso_id, analysis_group, avg_cni_pre_treatment)]
cni_pre_match[, group_label := fifelse(analysis_group == "Treated", "Treated", "Non-Treated")]

p_cni_pre <- ggplot(cni_pre_match[!is.na(avg_cni_pre_treatment)],
                    aes(x = avg_cni_pre_treatment, fill = group_label, color = group_label)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("Treated" = "#D55E00", "Non-Treated" = "#0072B2"), name = "") +
  scale_color_manual(values = c("Treated" = "#D55E00", "Non-Treated" = "#0072B2"), name = "") +
  labs(
    title = "CNI Distribution: Treated vs Non-Treated DeSOs",
    subtitle = "Before matching - All potential units (Non-treated = All untreated DeSOs)",
    x = "Care Need Index (CNI) Score",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(diag_dir, "cni_density_pre_matching.png"),
       p_cni_pre, width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 2. CNI density: treated vs matched controls (after matching)
# -----------------------------------------------------------------------------

matched_cni <- deso_summary[deso_id %in% matched_deso_ids,
                            .(deso_id, analysis_group, avg_cni_pre_treatment)]
matched_cni[, group_label := fifelse(analysis_group == "Treated", "Treated", "Matched Control")]

p_cni_post <- ggplot(matched_cni[!is.na(avg_cni_pre_treatment)],
                     aes(x = avg_cni_pre_treatment, fill = group_label, color = group_label)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("Treated" = "#D55E00", "Matched Control" = "#0072B2"), name = "") +
  scale_color_manual(values = c("Treated" = "#D55E00", "Matched Control" = "#0072B2"), name = "") +
  labs(
    title = "CNI Distribution: Treated vs Matched Control DeSOs",
    subtitle = sprintf("After matching - %d Treated, %d Matched Controls",
                       sum(matched_cni$group_label == "Treated"),
                       sum(matched_cni$group_label == "Matched Control")),
    x = "Care Need Index (CNI) Score",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(diag_dir, "cni_density_post_matching.png"),
       p_cni_post, width = 10, height = 6, dpi = 300)

# Combined comparison panel
p_cni_combined <- p_cni_pre / p_cni_post +
  plot_annotation(title = "CNI Balance: Before and After Matching")

ggsave(file.path(diag_dir, "cni_density_comparison.png"),
       p_cni_combined, width = 10, height = 10, dpi = 300)

cat("  Saved CNI density plots.\n")

# -----------------------------------------------------------------------------
# 3. Covariate distributions: boxplots by group
# -----------------------------------------------------------------------------
cat("--- Creating covariate distribution plots ---\n")

# Matched data for covariate comparison
matched_covars <- deso_birth_year_panel[deso_id %in% matched_deso_ids & year == 2017]
matched_covars <- merge(matched_covars,
                        deso_summary[, .(deso_id, analysis_group)],
                        by = "deso_id")
matched_covars[, group_label := fifelse(analysis_group == "Treated", "Treated", "Matched Control")]

# Reshape to long format for faceted boxplot
covar_vars <- c("prop_low_education", "prop_foreign_born", "prop_unemployed",
                "prop_single_parent", "prop_low_birth_weight", "avg_mother_age_at_birth", "cni_score")

covar_long <- melt(matched_covars,
                   id.vars = c("deso_id", "group_label"),
                   measure.vars = covar_vars,
                   variable.name = "covariate",
                   value.name = "value")

covar_labels <- c(
  "prop_low_education" = "Prop. Low Education",
  "prop_foreign_born" = "Prop. Foreign Born",
  "prop_unemployed" = "Prop. Unemployed",
  "prop_single_parent" = "Prop. Single Parent",
  "prop_low_birth_weight" = "Prop. Low Birth Weight",
  "avg_mother_age_at_birth" = "Avg. Mother Age",
  "cni_score" = "CNI Score"
)
covar_long[, covariate_label := covar_labels[as.character(covariate)]]

p_boxplot <- ggplot(covar_long[!is.na(value)],
                    aes(x = group_label, y = value, fill = group_label)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ covariate_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Treated" = "#D55E00", "Matched Control" = "#0072B2")) +
  labs(
    title = "Covariate Distributions: Treated vs Matched Controls",
    subtitle = "Post-matching balance check (Baseline year 2017)",
    x = "",
    y = "Value"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(diag_dir, "covariate_boxplots_post_matching.png"),
       p_boxplot, width = 14, height = 8, dpi = 300)

cat("  Saved covariate distribution plots.\n")

# -----------------------------------------------------------------------------
# 4. Treatment rollout timeline
# -----------------------------------------------------------------------------
cat("--- Creating treatment rollout timeline ---\n")

rollout_data <- deso_birth_year_panel[deso_id %in% deso_summary[analysis_group == "Treated", deso_id]]
rollout_data <- rollout_data[expected_deso_treatment_dosage >= TREATED_THRESHOLD,
                             .(first_treated_year = min(year)), by = deso_id]

rollout_summary <- rollout_data[, .(
  newly_treated = .N
), by = first_treated_year][order(first_treated_year)]

rollout_summary[, cumulative := cumsum(newly_treated)]
rollout_summary[, pct_total := round(cumulative / max(cumulative) * 100, 1)]

p_rollout <- ggplot(rollout_summary, aes(x = first_treated_year)) +
  geom_bar(aes(y = newly_treated), stat = "identity", fill = "#D55E00", alpha = 0.7) +
  geom_line(aes(y = cumulative), color = "#0072B2", linewidth = 1.2) +
  geom_point(aes(y = cumulative), color = "#0072B2", size = 2) +
  geom_text(aes(y = cumulative, label = cumulative), vjust = -0.5, color = "#0072B2", size = 3) +
  scale_x_continuous(breaks = rollout_summary$first_treated_year) +
  labs(
    title = "EHVP Treatment Rollout Timeline",
    subtitle = "Bars = Newly treated DeSOs per year | Line = Cumulative total",
    x = "Year First Achieved >=60% Dosage",
    y = "Number of DeSOs"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(diag_dir, "treatment_rollout_timeline.png"),
       p_rollout, width = 10, height = 6, dpi = 300)

fwrite(rollout_summary, file.path(diag_dir, "treatment_rollout_summary.csv"))

cat("  Saved rollout timeline.\n")

# -----------------------------------------------------------------------------
# 5. Binscatter: treatment dose vs CNI
# -----------------------------------------------------------------------------
cat("--- Creating dose-CNI relationship plot ---\n")

dose_cni_data <- deso_summary[max_dose > 0 & !is.na(avg_cni_pre_treatment)]

# Create decile bins
dose_cni_data[, dose_decile := cut(max_dose,
                                    breaks = quantile(max_dose, probs = seq(0, 1, 0.1), na.rm = TRUE),
                                    include.lowest = TRUE)]

binned_dose_cni <- dose_cni_data[!is.na(dose_decile), .(
  mean_dose = mean(max_dose, na.rm = TRUE),
  mean_cni = mean(avg_cni_pre_treatment, na.rm = TRUE),
  n_desos = .N
), by = dose_decile]

p_binscatter <- ggplot(binned_dose_cni, aes(x = mean_dose, y = mean_cni)) +
  geom_point(aes(size = n_desos), color = "#D55E00", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", linetype = "dashed") +
  scale_size_continuous(range = c(3, 8), name = "N DeSOs") +
  labs(
    title = "Relationship: Treatment Intensity and Socioeconomic Need",
    subtitle = "Each point = average within treatment dose decile (confirms targeting)",
    x = "Mean Maximum Treatment Dose",
    y = "Mean Pre-Treatment CNI Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(diag_dir, "binscatter_dose_vs_cni.png"),
       p_binscatter, width = 10, height = 7, dpi = 300)

cat("  Saved dose-CNI binscatter.\n")

# -----------------------------------------------------------------------------
# 6. Descriptive statistics tables
# -----------------------------------------------------------------------------
cat("--- Creating descriptive statistics tables ---\n")

# Pre-matching: All treated vs all untreated
desc_pre <- deso_summary[analysis_group %in% c("Treated", "Control"), .(
  N = .N,
  Mean_CNI = round(mean(avg_cni_pre_treatment, na.rm = TRUE), 2),
  SD_CNI = round(sd(avg_cni_pre_treatment, na.rm = TRUE), 2),
  Mean_Dose = round(mean(max_dose, na.rm = TRUE), 3)
), by = .(Group = analysis_group)]

fwrite(desc_pre, file.path(diag_dir, "descriptive_stats_pre_matching.csv"))

# Post-matching: Treated vs matched controls
desc_post <- matched_cni[, .(
  N = .N,
  Mean_CNI = round(mean(avg_cni_pre_treatment, na.rm = TRUE), 2),
  SD_CNI = round(sd(avg_cni_pre_treatment, na.rm = TRUE), 2)
), by = .(Group = group_label)]

fwrite(desc_post, file.path(diag_dir, "descriptive_stats_post_matching.csv"))

cat("  Saved descriptive statistics tables.\n")

# -----------------------------------------------------------------------------
# 7. Validation tables
# -----------------------------------------------------------------------------
cat("--- Creating validation tables ---\n")

# --- TABLE 1: Sample overview ---
sample_overview <- data.table(
  Metric = c(
    "Total DeSO-birth year observations",
    "Total unique DeSOs (with valid data)",
    "Treated DeSOs (>=60% dosage)",
    "Untreated DeSOs (<=10% dosage)",
    "Fuzzy Zone DeSOs (10-60%)",
    "Matched Treated DeSOs",
    "Matched Control DeSOs",
    "Total Matched DeSOs"
  ),
  Count = c(
    nrow(deso_birth_year_panel),
    uniqueN(deso_birth_year_panel$deso_id),
    nrow(deso_summary[analysis_group == "Treated"]),
    nrow(deso_summary[analysis_group == "Control"]),
    nrow(deso_summary[analysis_group == "Fuzzy Zone"]),
    length(matched_treated_ids),
    length(matched_control_ids),
    length(matched_deso_ids)
  )
)
fwrite(sample_overview, file.path(diag_dir, "table_1_sample_overview.csv"))

# --- TABLE 2: Treatment classification with CNI ---
treatment_classification <- deso_summary[, .(
  N_DeSOs = .N,
  Mean_CNI = round(mean(avg_cni_pre_treatment, na.rm = TRUE), 1),
  SD_CNI = round(sd(avg_cni_pre_treatment, na.rm = TRUE), 1),
  Mean_Max_Dose = round(mean(max_dose, na.rm = TRUE), 3)
), by = .(Group = analysis_group)]
setorder(treatment_classification, -N_DeSOs)
fwrite(treatment_classification, file.path(diag_dir, "table_2_treatment_classification.csv"))

# --- TABLE 3: Matching variables used ---
matching_vars_table <- data.table(
  Variable = c(MATCHING_COVARIATES, "CNI Score (for diagnostics)"),
  Description = c(
    "Proportion with low education (< high school)",
    "Proportion foreign-born",
    "Proportion unemployed",
    "Proportion single-parent households",
    "Proportion low birth weight (<2500g)",
    "Average mother age at birth",
    "Care Need Index (composite; shown in balance but not in PS)"
  ),
  Used_in_PS = c(rep("Yes", length(MATCHING_COVARIATES)), "No (diagnostic only)")
)
fwrite(matching_vars_table, file.path(diag_dir, "table_3_matching_variables.csv"))

# --- TABLE 4: Detailed descriptive statistics (matched sample) ---
baseline_matched <- deso_birth_year_panel[deso_id %in% matched_deso_ids & year == 2017]
baseline_matched <- merge(baseline_matched,
                          deso_summary[, .(deso_id, analysis_group)],
                          by = "deso_id")

desc_vars <- c("cni_score", "prop_low_education", "prop_foreign_born",
               "prop_unemployed", "prop_single_parent", "prop_low_birth_weight",
               "avg_mother_age_at_birth")

covar_pretty <- c(
  "cni_score" = "CNI Score",
  "prop_low_education" = "Prop. Low Education",
  "prop_foreign_born" = "Prop. Foreign-Born",
  "prop_unemployed" = "Prop. Unemployed",
  "prop_single_parent" = "Prop. Single Parent",
  "prop_low_birth_weight" = "Prop. Low Birth Weight",
  "avg_mother_age_at_birth" = "Avg. Mother Age"
)

desc_table_rows <- lapply(desc_vars, function(v) {
  treated_vals <- baseline_matched[analysis_group == "Treated", get(v)]
  control_vals <- baseline_matched[analysis_group == "Control", get(v)]

  data.table(
    Covariate = covar_pretty[v],
    Treated_Mean = round(mean(treated_vals, na.rm = TRUE), 3),
    Treated_SD = round(sd(treated_vals, na.rm = TRUE), 3),
    Treated_Min = round(min(treated_vals, na.rm = TRUE), 3),
    Treated_Max = round(max(treated_vals, na.rm = TRUE), 3),
    Control_Mean = round(mean(control_vals, na.rm = TRUE), 3),
    Control_SD = round(sd(control_vals, na.rm = TRUE), 3),
    Control_Min = round(min(control_vals, na.rm = TRUE), 3),
    Control_Max = round(max(control_vals, na.rm = TRUE), 3)
  )
})
detailed_desc_stats <- rbindlist(desc_table_rows)

fwrite(detailed_desc_stats, file.path(diag_dir, "table_4_descriptive_statistics_matched.csv"))

# HTML version
tryCatch({
  detailed_desc_stats %>%
    kable("html", caption = sprintf("Descriptive Statistics: Matched Sample (N=%d Treated, N=%d Control)",
                                    nrow(baseline_matched[analysis_group == "Treated"]),
                                    nrow(baseline_matched[analysis_group == "Control"])),
          col.names = c("Covariate", "Mean", "SD", "Min", "Max", "Mean", "SD", "Min", "Max")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
    add_header_above(c(" " = 1, "TREATED" = 4, "MATCHED CONTROL" = 4)) %>%
    save_kable(file = file.path(diag_dir, "table_4_descriptive_statistics_matched.html"))
}, error = function(e) cat(sprintf("  Note: Could not create HTML table: %s\n", e$message)))

# --- TABLE 5: Combined balance table (all cohorts) ---
all_balance <- rbindlist(lapply(names(match_results_list_stratified), function(cohort_str) {
  res <- match_results_list_stratified[[cohort_str]]
  if (!is.null(res$balance_table)) {
    bal <- as.data.table(res$balance_table)
    bal[, cohort := cohort_str]
    return(bal)
  }
  return(NULL)
}), fill = TRUE)

if (nrow(all_balance) > 0) {
  combined_balance <- all_balance[, .(
    Mean_SMD = round(mean(abs(Diff.Adj), na.rm = TRUE), 3),
    Max_SMD = round(max(abs(Diff.Adj), na.rm = TRUE), 3),
    Balanced = fifelse(max(abs(Diff.Adj), na.rm = TRUE) < 0.25, "Yes", "No")
  ), by = .(Variable = rownames)]

  fwrite(combined_balance, file.path(diag_dir, "table_5_combined_balance.csv"))
}

# --- TABLE 6: Cohort weights for stacked analysis ---
cohort_weights <- data.table(
  Cohort = as.integer(names(match_results_list_stratified)),
  Treated_DeSOs = sapply(match_results_list_stratified, function(x) x$n_treated),
  Matched_Controls = sapply(match_results_list_stratified, function(x) x$n_controls)
)
cohort_weights[, Total := Treated_DeSOs + Matched_Controls]
cohort_weights[, Weight_Pct := round(Treated_DeSOs / sum(Treated_DeSOs) * 100, 1)]
cohort_weights[, Ratio := sprintf("1:%.1f", Matched_Controls / Treated_DeSOs)]

fwrite(cohort_weights, file.path(diag_dir, "table_6_cohort_weights.csv"))

# --- TABLE 7: Sensitivity analysis thresholds ---
sensitivity_thresholds <- data.table(
  Threshold = c("50%", "60% (base)", "70%", "80%"),
  Treated_DeSOs = c(
    nrow(deso_summary[max_dose >= 0.50]),
    nrow(deso_summary[max_dose >= 0.60]),
    nrow(deso_summary[max_dose >= 0.70]),
    nrow(deso_summary[max_dose >= 0.80])
  )
)
sensitivity_thresholds[, Reduction := c(
  NA,
  "---",
  sprintf("-%d%%", round((1 - Treated_DeSOs[3]/Treated_DeSOs[2]) * 100)),
  sprintf("-%d%%", round((1 - Treated_DeSOs[4]/Treated_DeSOs[2]) * 100))
)]
sensitivity_thresholds[, Trade_off := c(
  "More power, more dilution",
  "Current specification",
  "Cleaner treatment, less power",
  "Purest treatment definition"
)]

fwrite(sensitivity_thresholds, file.path(diag_dir, "table_7_sensitivity_thresholds.csv"))

cat("  Saved validation tables (7 tables).\n")

cat(sprintf("\nDiagnostic plots saved to: %s\n", diag_dir))

# ============================================================================
# GEOGRAPHIC MAPS
# ============================================================================

cat("\n--- Creating geographic maps ---\n")

if (file.exists(SPATIAL_DATA_PATH)) {
  tryCatch({
    deso_sf <- readRDS(SPATIAL_DATA_PATH) %>%
      st_transform(crs = 4326) %>%
      dplyr::rename(deso_id = desokod)

    map_data <- deso_sf %>%
      dplyr::mutate(
        status = dplyr::case_when(
          deso_id %in% matched_treated_ids ~ "Treated",
          deso_id %in% matched_control_ids ~ "Matched Control",
          TRUE ~ "Not in Analysis"
        )
      )

    map_colors <- c("Treated" = "#D55E00", "Matched Control" = "#0072B2", "Not in Analysis" = "grey90")
    map_data$status <- factor(map_data$status, levels = names(map_colors))

    national_map <- ggplot(data = map_data) +
      geom_sf(aes(fill = status), color = "white", lwd = 0.02) +
      scale_fill_manual(values = map_colors, name = "Analysis Sample", drop = TRUE) +
      labs(
        title = "Matched Sample: Geographic Distribution",
        subtitle = sprintf("%d Treated + %d Matched Controls",
                           length(matched_treated_ids), length(matched_control_ids))
      ) +
      theme_void(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave(file.path(MATCHING_OUTPUT_DIR, "maps", "map_matched_sample.png"),
           national_map, width = 10, height = 9, dpi = 300)

    cat("  Maps saved.\n")

  }, error = function(e) {
    cat(sprintf("  Map generation failed: %s\n", e$message))
  })
} else {
  cat("  Spatial data not found, skipping maps.\n")
}

# ============================================================================
# SAVE OUTPUTS FOR ANALYSIS SCRIPT
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("SAVING OUTPUTS\n")
cat("============================================================================\n\n")

# Save matched data by cohort
for (cohort_str in names(match_results_list_stratified)) {
  cohort_data <- match_results_list_stratified[[cohort_str]]$matched_data
  fwrite(cohort_data, file.path(MATCHING_OUTPUT_DIR,
                                 sprintf("matched_data_cohort_%s.csv", cohort_str)))
}

# Save main panel
fwrite(final_matched_panel, file.path(MATCHING_OUTPUT_DIR, "final_matched_panel.csv"))

# Save all objects needed by analysis script
save(
  # DeSO-level objects
  matched_deso_ids,
  matched_treated_ids,
  matched_control_ids,
  final_matched_panel,
  deso_panel,
  deso_summary,
  deso_birth_year_panel,

  # Matching objects
  match_results_list_stratified,
  matching_summary,
  balance_summary,

  # Parameters
  COHORTS_TO_MATCH,
  TREATED_THRESHOLD,
  CONTROL_THRESHOLD,
  MATCHING_COVARIATES,

  file = MATCHING_OBJECTS_PATH
)

cat(sprintf("Matching objects saved to: %s\n", MATCHING_OBJECTS_PATH))

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("MATCHING COMPLETE\n")
cat("============================================================================\n")
cat("\n")
cat("SPECIFICATION:\n")
cat("  - Exact match on: Urbanicity (Urban/Suburban/Rural)\n")
cat("  - Caliper: None\n")
cat("  - Replacement: Yes\n")
cat("  - Ratio: Up to 5:1\n")
cat("  - Region: NOT matched on (controlled via FE in outcome model)\n")
cat("\n")
cat("SAMPLE:\n")
for (i in 1:nrow(matching_summary)) {
  cat(sprintf("  Cohort %d: %d treated -> %d controls\n",
              matching_summary$Cohort[i],
              matching_summary$N_Treated[i],
              matching_summary$N_Matched_Controls[i]))
}
cat(sprintf("\n  TOTAL: %d treated + %d control DeSOs\n",
            length(matched_treated_ids), length(matched_control_ids)))
cat("\n")
cat("BALANCE (Max |SMD| per cohort):\n")
for (i in 1:nrow(balance_summary)) {
  status <- ifelse(balance_summary$Max_SMD[i] < SMD_EXCELLENT, "Excellent",
                   ifelse(balance_summary$Max_SMD[i] < SMD_ACCEPTABLE, "Acceptable", "Poor"))
  cat(sprintf("  Cohort %s: %.3f %s\n",
              balance_summary$Cohort[i],
              balance_summary$Max_SMD[i],
              status))
}
cat("\n")
cat("OUTPUT FILES:\n")
cat(sprintf("  %s/\n", MATCHING_OUTPUT_DIR))
cat("    - deso_birth_year_panel.csv\n")
cat("    - deso_classification.csv\n")
cat("    - final_matched_panel.csv\n")
cat("    - matched_data_cohort_YYYY.csv\n")
cat("    - matching_summary.csv/.html\n")
cat("    - balance_summary.csv\n")
cat("    - balance_tables/balance_cohort_YYYY.csv/.html\n")
cat("    - love_plots/love_plot_cohort_YYYY.png\n")
cat("    - love_plots/love_plot_all_cohorts.png\n")
cat("    - love_plots/ps_overlap_cohort_YYYY.png\n")
cat("    - maps/map_matched_sample.png\n")
cat(sprintf("    - matching_objects.RData (for 02_analysis.R)\n"))
cat("\n")
cat(sprintf("Run completed: %s\n", Sys.time()))
cat("============================================================================\n")
