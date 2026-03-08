# =============================================================================
# Title:        CNI Measurement Validity Analysis
# Description:  Tests whether the Care Need Index (CNI) predicts pediatric
#               healthcare burden at the neighbourhood level. The analysis
#               includes construct validity (saturation), comparative fit
#               against a household vulnerability index (HVI), functional-form
#               checks, age heterogeneity, out-of-sample prediction, regional
#               validity, targeting efficiency, within-neighbourhood fixed
#               effects, lagged-HVI sensitivity, and a Skane deep-dive.
# Paper:        Paper I -- Targeting (EHVP thesis)
# Author:       Sergio Flores
# Date:         2026-01-04
# Dependencies: data.table, fixest, ggplot2, modelsummary, scales, parallel
#
# USAGE:
#   Set `base_dir` and `output_root` below to match your local file layout,
#   then source this script. All outputs (tables as CSV, figures as PNG) are
#   written to a date-stamped subfolder under `output_root`.
# =============================================================================

rm(list = ls()); gc()

library(data.table)
library(fixest)
library(ggplot2)
library(modelsummary)
library(scales)
library(parallel)

setDTthreads(parallel::detectCores())
setFixest_nthreads(parallel::detectCores())
cat("Using", parallel::detectCores(), "CPU cores for parallel processing\n")

# --- CONFIGURATION ----------------------------------------------------------
# Set these two paths before running the script.
base_dir    <- file.path("data")
output_root <- file.path("output")
output_dir  <- file.path(output_root,
                         paste0("CNI_Validity_", format(Sys.Date(), "%Y%m%d")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: print a summary to console and save as CSV
save_and_print <- function(dt, name) {
  cat(paste0("\n--- RESULTS: ", name, " ---\n"))
  print(dt)
  fwrite(dt, file.path(output_dir, paste0(name, ".csv")))
}

# ---------------------------------------------------------------------------
# PHASE 1: DATA LOADING AND PREPARATION
# ---------------------------------------------------------------------------
cat("--- PHASE 1: Data Preparation ---\n")
input_file <- file.path(base_dir,
                         "final_panel_data_complete_with_FK_MBR_Meds_20260104.csv")

cols <- c("lopnr", "year", "deso_2021", "lan_res", "child_birth_year", "child_kon",
          "relative_cni", "cni_score",
          "mother_Sun2000niva", "father_Sun2000niva",
          "mother_SyssStat11", "father_SyssStat11",
          "mother_FamTypF", "father_FamTypF",
          "mother_fodelseland", "father_fodelseland",
          "parent_had_dx_mental_health_common_year",
          "inpatient_avoidable_events", "outpatient_avoidable_events",
          "inpatient_emergency_events", "events_antibiotic_systemic_year",
          "parental_tfp_gross_days_year",
          "inpatient_trauma_maltreatment_syndrome_events",
          "outpatient_trauma_maltreatment_syndrome_events",
          "inpatient_trauma_burns_events", "outpatient_trauma_burns_events",
          "inpatient_trauma_poisoning_icd_events",
          "outpatient_trauma_poisoning_icd_events")

# `lopnr` is a pseudonymised individual identifier assigned by the register
# authority; no personal identification numbers are present in the data.
dt <- fread(input_file,
            select = intersect(cols, names(fread(input_file, nrows = 1))))

# Restrict to study period 2012-2017, ages 0-5, and the four study regions
dt <- dt[year %between% c(2012, 2017) & lan_res %in% c(1, 12, 14, 18)]

dt[, child_age := year - child_birth_year]
dt <- dt[child_age %between% c(0, 5)]

cat("After year/region/age filter. N =", format(nrow(dt), big.mark = ","), "\n")

# --- Construct the Household Vulnerability Index (HVI) ---------------------
# Each flag captures one dimension of household-level disadvantage.

# Low parental education (pre-secondary or short secondary)
dt[, hvi_flag_low_ed := fcase(
  (mother_Sun2000niva %in% c(100, 200, 310, 312)) |
    (father_Sun2000niva %in% c(100, 200, 310, 312)), 1L,
  default = 0L)]

# Parental unemployment or labour-market exclusion
dt[, hvi_flag_unemp := fcase(
  (mother_SyssStat11 %in% c(5, 6, 7)) |
    (father_SyssStat11 %in% c(5, 6, 7)), 1L,
  default = 0L)]

# Single-parent household
dt[, hvi_flag_single := fcase(
  (mother_FamTypF %in% c(31, 32, 41, 42)) |
    (father_FamTypF %in% c(31, 32, 41, 42)), 1L,
  default = 0L)]

# At least one foreign-born parent (exclude empty strings from classification)
dt[, hvi_flag_foreign := fcase(
  (mother_fodelseland != "" & mother_fodelseland != "SVERIGE") |
    (father_fodelseland != "" & father_fodelseland != "SVERIGE"), 1L,
  default = 0L)]

# Parental common mental disorder diagnosis during the year
dt[, hvi_flag_mental := fcase(
  parent_had_dx_mental_health_common_year == 1, 1L,
  default = 0L)]

dt[, hvi_score := rowSums(.SD), .SDcols = patterns("hvi_flag_")]

# Standardise exposure variables (wrap in as.numeric to avoid matrix columns)
dt[, cni_std := as.numeric(scale(relative_cni))]
dt[, hvi_std := as.numeric(scale(hvi_score))]
dt[, cni_sq  := cni_std^2]

# --- Define outcome constructs ---------------------------------------------
dt[, val_burden     := inpatient_avoidable_events + outpatient_avoidable_events]
dt[, val_safety     := rowSums(.SD, na.rm = TRUE), .SDcols = patterns("trauma")]
dt[, val_efficiency := inpatient_emergency_events]
dt[, val_practice   := events_antibiotic_systemic_year]
dt[, val_capacity   := parental_tfp_gross_days_year]

outcomes <- c("val_burden", "val_safety", "val_efficiency",
              "val_practice", "val_capacity")
for (v in outcomes) dt[is.na(get(v)), (v) := 0]

# Region labels for stratified analyses
dt[, region_name := fcase(
  lan_res == 1,  "Stockholm",
  lan_res == 12, "Skane",
  lan_res == 14, "Vastra Gotaland",
  lan_res == 18, "Orebro")]

cat("Full sample N =", format(nrow(dt), big.mark = ","), "\n")
cat("Missing CNI N =", format(sum(is.na(dt$cni_std)), big.mark = ","), "\n")

# Common analytic sample: complete cases on both indices and the main outcome
dt_analytic <- dt[!is.na(cni_std) & !is.na(hvi_std) & !is.na(val_burden)]
cat("Analytic sample N (complete cases) =",
    format(nrow(dt_analytic), big.mark = ","), "\n")

# ---------------------------------------------------------------------------
# PHASE 2: CONSTRUCT VALIDITY -- Saturation / Quadratic Check
# ---------------------------------------------------------------------------
cat("\n--- PHASE 2: Construct Validity (Saturation Check) ---\n")

# A. Decile plot with quadratic overlay
dt[, cni_decile := cut(relative_cni,
                       breaks = quantile(relative_cni, probs = 0:10/10,
                                         na.rm = TRUE),
                       labels = 1:10, include.lowest = TRUE)]
grad_burden <- dt[, .(mean_burden = mean(val_burden, na.rm = TRUE)),
                  by = cni_decile]

gg_sat <- ggplot(grad_burden,
                 aes(x = as.numeric(cni_decile), y = mean_burden)) +
  geom_point(size = 4, color = "#D55E00") +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2),
              se = FALSE, color = "black", linetype = "dashed") +
  labs(title = "Saturation Check: Clinical Burden",
       subtitle = "Mean burden by CNI decile with quadratic fit overlay",
       x = "CNI Decile (1 = Low, 10 = High)",
       y = "Mean Events per Child-Year") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "Fig1_Burden_Saturation_Check.png"),
       gg_sat, width = 8, height = 6)

# B. Formal quadratic term test for all outcomes (DeSO-clustered SEs)
sat_res <- data.table()
for (out in outcomes) {
  m_quad <- fepois(as.formula(paste(out,
                    "~ cni_std + cni_sq + factor(child_age) | year")),
                   data = dt_analytic, cluster = "deso_2021")
  ct <- coeftable(m_quad)
  sat_res <- rbind(sat_res, data.table(
    Outcome    = out,
    Linear_Beta = round(ct["cni_std", 1], 4),
    Linear_SE   = round(ct["cni_std", 2], 4),
    Quad_Beta   = round(ct["cni_sq", 1], 4),
    Quad_SE     = round(ct["cni_sq", 2], 4),
    P_Quad      = round(ct["cni_sq", 4], 4)
  ))
}
save_and_print(sat_res, "Table1_Quadratic_Saturation_Tests")

# ---------------------------------------------------------------------------
# PHASE 3: VARIANCE-EXPLAINED HORSE RACE (CNI vs HVI)
# ---------------------------------------------------------------------------
cat("\n--- PHASE 3: Variance Explained Horse Race ---\n")

# All four models are estimated on the same analytic sample so that
# pseudo-R-squared values are directly comparable.
f0 <- val_burden ~ factor(child_age) | lan_res + year
f1 <- val_burden ~ cni_std + factor(child_age) | lan_res + year
f2 <- val_burden ~ hvi_std + factor(child_age) | lan_res + year
f3 <- val_burden ~ cni_std + hvi_std + factor(child_age) | lan_res + year

m0 <- fepois(f0, dt_analytic, cluster = "deso_2021")
m1 <- fepois(f1, dt_analytic, cluster = "deso_2021")
m2 <- fepois(f2, dt_analytic, cluster = "deso_2021")
m3 <- fepois(f3, dt_analytic, cluster = "deso_2021")

get_fit <- function(m, name) {
  r2_val <- r2(m, "pr2")["pr2"]
  ll     <- logLik(m)
  n      <- nobs(m)
  data.table(Model = name, N = n,
             Pseudo_R2 = as.numeric(r2_val),
             LogLik    = as.numeric(ll))
}

horse_race <- rbind(
  get_fit(m0, "M0: Demographics Only"),
  get_fit(m1, "M1: Demo + CNI"),
  get_fit(m2, "M2: Demo + HVI"),
  get_fit(m3, "M3: Demo + CNI + HVI")
)

base_r2 <- horse_race[Model == "M0: Demographics Only", Pseudo_R2]
horse_race[, Delta_R2_vs_Base := Pseudo_R2 - base_r2]
horse_race[, Delta_R2_pct := round(Delta_R2_vs_Base * 100, 4)]

# Verify that all models use the same N
cat("Sample sizes (should all match):\n")
print(horse_race[, .(Model, N)])

save_and_print(horse_race, "Table2_Variance_Explained_Horse_Race")

# Report the incremental R-squared ratio
cni_delta <- horse_race[Model == "M1: Demo + CNI", Delta_R2_vs_Base]
hvi_delta <- horse_race[Model == "M2: Demo + HVI", Delta_R2_vs_Base]
cat(sprintf("\nCNI incremental R-squared: %.5f\n", cni_delta))
cat(sprintf("HVI incremental R-squared: %.5f\n", hvi_delta))
cat(sprintf("Ratio (CNI / HVI): %.2fx\n", cni_delta / hvi_delta))

# ---------------------------------------------------------------------------
# PHASE 4: FUNCTIONAL FORM -- OLS vs PPML
# ---------------------------------------------------------------------------
cat("\n--- PHASE 4: Functional Form (OLS vs PPML) ---\n")

m_ppml <- fepois(val_burden ~ cni_std + factor(child_age) | lan_res + year,
                 data = dt_analytic, cluster = "deso_2021")

m_ols <- feols(log(1 + val_burden) ~ cni_std + factor(child_age) | lan_res + year,
               data = dt_analytic, cluster = "deso_2021")

robust_check <- data.table(
  Specification = c("PPML (Count)", "OLS (Log-Linear)"),
  Outcome       = "Clinical Burden",
  Beta_CNI  = c(coeftable(m_ppml)["cni_std", 1],
                coeftable(m_ols)["cni_std", 1]),
  SE        = c(coeftable(m_ppml)["cni_std", 2],
                coeftable(m_ols)["cni_std", 2]),
  P_Value   = c(coeftable(m_ppml)["cni_std", 4],
                coeftable(m_ols)["cni_std", 4]),
  Significant = c(coeftable(m_ppml)["cni_std", 4] < 0.05,
                  coeftable(m_ols)["cni_std", 4] < 0.05)
)
save_and_print(robust_check, "Table3_Functional_Form_Robustness")

# ---------------------------------------------------------------------------
# PHASE 5: AGE HETEROGENEITY
# ---------------------------------------------------------------------------
cat("\n--- PHASE 5: Age Heterogeneity ---\n")

run_age_split <- function(subset_data, label) {
  m <- fepois(val_burden ~ cni_std + factor(child_age) | lan_res + year,
              data = subset_data, cluster = "deso_2021")
  data.table(
    Group     = label,
    N         = nobs(m),
    Beta_CNI  = round(coeftable(m)["cni_std", 1], 4),
    SE        = round(coeftable(m)["cni_std", 2], 4),
    Pseudo_R2 = round(r2(m, "pr2")["pr2"], 5)
  )
}

age_res <- rbind(
  run_age_split(dt_analytic[child_age %in% 0:2], "Infants/Toddlers (0-2)"),
  run_age_split(dt_analytic[child_age %in% 3:5], "Preschoolers (3-5)")
)
save_and_print(age_res, "Table4_Age_Heterogeneity")

# ---------------------------------------------------------------------------
# PHASE 6: OUT-OF-SAMPLE VALIDATION
# ---------------------------------------------------------------------------
cat("\n--- PHASE 6: Out-of-Sample Validation ---\n")

set.seed(42)

train_idx <- sample(seq_len(nrow(dt_analytic)),
                    size = 0.8 * nrow(dt_analytic))
dt_train  <- dt_analytic[train_idx]
dt_test   <- dt_analytic[-train_idx]

cat(sprintf("Training N: %s, Test N: %s\n",
            format(nrow(dt_train), big.mark = ","),
            format(nrow(dt_test),  big.mark = ",")))

m_cni  <- fepois(val_burden ~ cni_std + factor(child_age) | lan_res + year,
                 data = dt_train)
m_hvi  <- fepois(val_burden ~ hvi_std + factor(child_age) | lan_res + year,
                 data = dt_train)
m_both <- fepois(val_burden ~ cni_std + hvi_std + factor(child_age) | lan_res + year,
                 data = dt_train)

pred_cni  <- predict(m_cni,  newdata = dt_test)
pred_hvi  <- predict(m_hvi,  newdata = dt_test)
pred_both <- predict(m_both, newdata = dt_test)

rmse_calc <- function(actual, pred) {
  sqrt(mean((actual - pred)^2, na.rm = TRUE))
}

oos_res <- data.table(
  Model    = c("CNI Only", "HVI Only", "Both"),
  OOS_RMSE = c(rmse_calc(dt_test$val_burden, pred_cni),
               rmse_calc(dt_test$val_burden, pred_hvi),
               rmse_calc(dt_test$val_burden, pred_both))
)
save_and_print(oos_res, "Table5_Out_of_Sample_RMSE")

# ---------------------------------------------------------------------------
# PHASE 7: REGIONAL VALIDITY
# ---------------------------------------------------------------------------
cat("\n--- PHASE 7: Regional Validity ---\n")

r2_reg_res <- data.table()
regions    <- c("Stockholm", "Skane", "Vastra Gotaland", "Orebro")

for (r in regions) {
  sub_dt <- dt_analytic[region_name == r]

  m_base <- fepois(val_burden ~ factor(child_age) | year,
                   data = sub_dt, cluster = "deso_2021")
  m_cni  <- fepois(val_burden ~ cni_std + factor(child_age) | year,
                   data = sub_dt, cluster = "deso_2021")

  r2_gain <- r2(m_cni, "pr2")["pr2"] - r2(m_base, "pr2")["pr2"]

  r2_reg_res <- rbind(r2_reg_res, data.table(
    Region         = r,
    N              = nobs(m_cni),
    Beta_CNI       = round(coeftable(m_cni)["cni_std", 1], 4),
    SE             = round(coeftable(m_cni)["cni_std", 2], 4),
    Incremental_R2 = round(as.numeric(r2_gain), 5)
  ))
}
save_and_print(r2_reg_res, "Table6_Regional_Validity_Metrics")

# Regional gradient plot
dt_analytic[, cni_decile := cut(relative_cni,
                                breaks = quantile(relative_cni,
                                                  probs = 0:10/10,
                                                  na.rm = TRUE),
                                labels = 1:10, include.lowest = TRUE)]

reg_gradient <- dt_analytic[, .(mean_burden = mean(val_burden, na.rm = TRUE),
                                N = .N),
                            by = .(region_name, cni_decile)]
reg_gradient <- reg_gradient[!is.na(cni_decile) & N > 50]

gg_reg_grad <- ggplot(reg_gradient,
                      aes(x = as.numeric(cni_decile), y = mean_burden,
                          color = region_name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~region_name, scales = "free_y") +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Regional CNI-Burden Gradients",
       subtitle = "Mean clinical burden by national CNI decile and region",
       x = "National CNI Decile (1 = Low, 10 = High)",
       y = "Mean Annual Burden (Events per Child)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "Fig3B_Regional_Gradients.png"),
       gg_reg_grad, width = 10, height = 8)

# ---------------------------------------------------------------------------
# PHASE 8: TARGETING EFFICIENCY (Concentration Curve)
# ---------------------------------------------------------------------------
cat("\n--- PHASE 8: Targeting Efficiency ---\n")

dt_rank <- dt_analytic[!is.na(relative_cni) & !is.na(val_burden)]
setorder(dt_rank, -relative_cni)

dt_rank[, cum_pop        := seq_len(.N)]
dt_rank[, cum_pop_pct    := cum_pop / .N]
dt_rank[, cum_burden     := cumsum(val_burden)]
dt_rank[, cum_burden_pct := cum_burden / sum(val_burden)]

thresholds     <- c(0.1, 0.2, 0.3, 0.5)
efficiency_tab <- data.table()

for (t in thresholds) {
  row_idx  <- which.min(abs(dt_rank$cum_pop_pct - t))
  captured <- dt_rank[row_idx, cum_burden_pct]

  efficiency_tab <- rbind(efficiency_tab, data.table(
    Targeted_Population = paste0("Top ", t * 100, "%"),
    Burden_Captured     = round(captured, 4),
    Lift                = round(captured / t, 3)
  ))
}
save_and_print(efficiency_tab, "Table7_Targeting_Efficiency")

# Concentration curve
dt_plot <- dt_rank[seq(1, nrow(dt_rank), length.out = 1000)]
gg_conc <- ggplot(dt_plot, aes(x = cum_pop_pct, y = cum_burden_pct)) +
  geom_line(linewidth = 1.5, color = "#D55E00") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Targeting Validity: Concentration Curve",
       subtitle = "Cumulative burden captured by targeting the top X% by CNI",
       x = "Cumulative % of Population (sorted by CNI)",
       y = "Cumulative % of Burden Captured") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "Fig4_Targeting_Efficiency.png"),
       gg_conc, width = 8, height = 6)

# ---------------------------------------------------------------------------
# PHASE 9: WITHIN-DeSO FIXED EFFECTS
# ---------------------------------------------------------------------------
cat("\n--- PHASE 9: Within-DeSO Fixed Effects Analysis ---\n")
cat("Testing whether within-neighbourhood CNI variation predicts burden\n")

# The DeSO fixed-effects specification absorbs all time-invariant
# neighbourhood characteristics, providing a stronger test of whether
# changes in CNI track changes in health outcomes within the same area.

m_within  <- fepois(val_burden ~ cni_std + factor(child_age) | deso_2021 + year,
                    data = dt_analytic, cluster = "deso_2021")
m_between <- fepois(val_burden ~ cni_std + factor(child_age) | lan_res + year,
                    data = dt_analytic, cluster = "deso_2021")

within_res <- data.table(
  Specification = c("Between-DeSO (Region FE)", "Within-DeSO (DeSO FE)"),
  Beta_CNI  = c(coeftable(m_between)["cni_std", 1],
                coeftable(m_within)["cni_std", 1]),
  SE        = c(coeftable(m_between)["cni_std", 2],
                coeftable(m_within)["cni_std", 2]),
  P_Value   = c(coeftable(m_between)["cni_std", 4],
                coeftable(m_within)["cni_std", 4])
)
save_and_print(within_res, "Table8_Within_DeSO_Analysis")

cat("\nIf the within-DeSO coefficient is smaller, between-area differences\n")
cat("may partly reflect sorting rather than causal neighbourhood effects.\n")

# ---------------------------------------------------------------------------
# PHASE 10: LAGGED HVI SENSITIVITY
# ---------------------------------------------------------------------------
cat("\n--- PHASE 10: Lagged HVI Sensitivity ---\n")
cat("Using t-1 HVI to address potential reverse causality\n")

# Lag the HVI score by one year within each child
setorder(dt_analytic, lopnr, year)
dt_analytic[, hvi_score_lag := shift(hvi_score, n = 1, type = "lag"),
            by = lopnr]
dt_analytic[, hvi_std_lag := as.numeric(scale(hvi_score_lag))]

dt_lag <- dt_analytic[!is.na(hvi_std_lag)]
cat("Sample with lagged HVI: N =", format(nrow(dt_lag), big.mark = ","), "\n")

m_hvi_contemp <- fepois(val_burden ~ hvi_std + factor(child_age) | lan_res + year,
                        data = dt_lag, cluster = "deso_2021")
m_hvi_lagged  <- fepois(val_burden ~ hvi_std_lag + factor(child_age) | lan_res + year,
                        data = dt_lag, cluster = "deso_2021")
m_cni_only    <- fepois(val_burden ~ cni_std + factor(child_age) | lan_res + year,
                        data = dt_lag, cluster = "deso_2021")
m_cni_lag_hvi <- fepois(val_burden ~ cni_std + hvi_std_lag + factor(child_age) | lan_res + year,
                        data = dt_lag, cluster = "deso_2021")

lag_res <- data.table(
  Model = c("HVI (Contemporaneous)", "HVI (Lagged t-1)",
            "CNI Only", "CNI + Lagged HVI"),
  N = c(nobs(m_hvi_contemp), nobs(m_hvi_lagged),
        nobs(m_cni_only), nobs(m_cni_lag_hvi)),
  Beta_Main = c(
    coeftable(m_hvi_contemp)["hvi_std", 1],
    coeftable(m_hvi_lagged)["hvi_std_lag", 1],
    coeftable(m_cni_only)["cni_std", 1],
    coeftable(m_cni_lag_hvi)["cni_std", 1]
  ),
  SE = c(
    coeftable(m_hvi_contemp)["hvi_std", 2],
    coeftable(m_hvi_lagged)["hvi_std_lag", 2],
    coeftable(m_cni_only)["cni_std", 2],
    coeftable(m_cni_lag_hvi)["cni_std", 2]
  ),
  Pseudo_R2 = c(
    r2(m_hvi_contemp, "pr2")["pr2"],
    r2(m_hvi_lagged, "pr2")["pr2"],
    r2(m_cni_only, "pr2")["pr2"],
    r2(m_cni_lag_hvi, "pr2")["pr2"]
  )
)
save_and_print(lag_res, "Table9_Lagged_HVI_Sensitivity")

cat("\nIf the lagged HVI coefficient is smaller than the contemporaneous one,\n")
cat("reverse causality (child health -> parental outcomes) may be at play.\n")

# ---------------------------------------------------------------------------
# PHASE 11: SKANE DEEP-DIVE
# ---------------------------------------------------------------------------
cat("\n--- PHASE 11: Skane Deep-Dive ---\n")
cat("Investigating the weaker CNI validity observed in Skane\n")

# Malmo (kommun 1280) is the largest city in Skane; DeSO codes starting
# with "1280" identify Malmo neighbourhoods.

dt_skane <- dt_analytic[region_name == "Skane"]
dt_skane[, is_malmo := substr(deso_2021, 1, 4) == "1280"]

cat("\nSkane sample breakdown:\n")
cat("  Malmo:          N =", format(sum(dt_skane$is_malmo),  big.mark = ","), "\n")
cat("  Rest of Skane:  N =", format(sum(!dt_skane$is_malmo), big.mark = ","), "\n")

m_malmo       <- fepois(val_burden ~ cni_std + factor(child_age) | year,
                        data = dt_skane[is_malmo == TRUE],
                        cluster = "deso_2021")
m_other_skane <- fepois(val_burden ~ cni_std + factor(child_age) | year,
                        data = dt_skane[is_malmo == FALSE],
                        cluster = "deso_2021")
m_all_skane   <- fepois(val_burden ~ cni_std + factor(child_age) | year,
                        data = dt_skane, cluster = "deso_2021")

skane_res <- data.table(
  Subregion = c("All Skane", "Malmo Only", "Rest of Skane"),
  N = c(nobs(m_all_skane), nobs(m_malmo), nobs(m_other_skane)),
  Beta_CNI = c(
    coeftable(m_all_skane)["cni_std", 1],
    coeftable(m_malmo)["cni_std", 1],
    coeftable(m_other_skane)["cni_std", 1]
  ),
  SE = c(
    coeftable(m_all_skane)["cni_std", 2],
    coeftable(m_malmo)["cni_std", 2],
    coeftable(m_other_skane)["cni_std", 2]
  ),
  P_Value = c(
    coeftable(m_all_skane)["cni_std", 4],
    coeftable(m_malmo)["cni_std", 4],
    coeftable(m_other_skane)["cni_std", 4]
  )
)
save_and_print(skane_res, "Table10_Skane_DeepDive_MalmoVsRest")

# CNI distribution by region
cat("\nCNI Distribution by Region:\n")
cni_dist <- dt_analytic[, .(
  Mean_CNI = mean(relative_cni, na.rm = TRUE),
  SD_CNI   = sd(relative_cni, na.rm = TRUE),
  Min_CNI  = min(relative_cni, na.rm = TRUE),
  Max_CNI  = max(relative_cni, na.rm = TRUE),
  N_DeSO   = uniqueN(deso_2021)
), by = region_name]
print(cni_dist)
fwrite(cni_dist, file.path(output_dir, "Table11_CNI_Distribution_by_Region.csv"))

# Burden distribution by region
cat("\nBurden Distribution by Region:\n")
burden_dist <- dt_analytic[, .(
  Mean_Burden = mean(val_burden, na.rm = TRUE),
  SD_Burden   = sd(val_burden, na.rm = TRUE),
  Pct_Zero    = mean(val_burden == 0) * 100
), by = region_name]
print(burden_dist)
fwrite(burden_dist,
       file.path(output_dir, "Table12_Burden_Distribution_by_Region.csv"))

# Malmo vs rest-of-Skane gradient plot
skane_gradient <- dt_skane[, .(mean_burden = mean(val_burden, na.rm = TRUE),
                               N = .N),
                           by = .(is_malmo, cni_decile)]
skane_gradient[, location := ifelse(is_malmo, "Malmo", "Rest of Skane")]

gg_skane <- ggplot(skane_gradient[!is.na(cni_decile)],
                   aes(x = as.numeric(cni_decile), y = mean_burden,
                       color = location)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Skane Deep-Dive: Malmo vs Rest of Region",
       subtitle = "CNI-Burden gradient comparison",
       x = "National CNI Decile",
       y = "Mean Burden (Events per Child-Year)",
       color = "Location") +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "Fig5_Skane_DeepDive.png"),
       gg_skane, width = 9, height = 6)

# ---------------------------------------------------------------------------
# APPENDIX: COMPLETE REGIONAL-OUTCOME MATRIX
# ---------------------------------------------------------------------------
cat("\n--- APPENDIX: Complete Regional-Outcome Matrix ---\n")

res_matrix <- data.table()
for (r in regions) {
  for (o in outcomes) {
    m <- fepois(as.formula(paste(o,
                  "~ cni_std + factor(child_age) | year")),
                data = dt_analytic[region_name == r],
                cluster = "deso_2021")
    res_matrix <- rbind(res_matrix, data.table(
      Region  = r,
      Outcome = o,
      Beta    = round(coeftable(m)["cni_std", 1], 4),
      SE      = round(coeftable(m)["cni_std", 2], 4),
      P_Value = round(coeftable(m)["cni_std", 4], 4)
    ))
  }
}
save_and_print(res_matrix, "TableS5_Complete_Regional_Outcome_Matrix")

# ---------------------------------------------------------------------------
# WRAP-UP
# ---------------------------------------------------------------------------
cat("\n=========================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=========================================================================\n")
cat("\nOutput directory:", output_dir, "\n")
cat("\nOutputs:\n")
cat("  Table 1:  Quadratic saturation tests\n")
cat("  Table 2:  Variance-explained horse race (CNI vs HVI)\n")
cat("  Table 3:  Functional form (OLS vs PPML)\n")
cat("  Table 4:  Age heterogeneity\n")
cat("  Table 5:  Out-of-sample RMSE\n")
cat("  Table 6:  Regional validity metrics\n")
cat("  Table 7:  Targeting efficiency\n")
cat("  Table 8:  Within-DeSO fixed effects\n")
cat("  Table 9:  Lagged HVI sensitivity\n")
cat("  Table 10: Skane -- Malmo vs rest\n")
cat("  Table 11-12: Distribution summaries\n")
cat("  Table S5: Complete regional-outcome matrix\n")
cat("  Figures 1, 3B, 4, 5\n")
