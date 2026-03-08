# ==============================================================================
# Title:        Budget Impact Analysis of the Extended Home Visiting Programme
# Description:  Estimates national-level programme costs under three delivery
#               scenarios (efficient, base case, intensive) and runs sensitivity
#               analyses on uptake rates and regional/municipal payer splits.
# Paper:        Paper IV in Flores (2026), doctoral thesis, Uppsala University
# Author:       Sergio Flores
# Date:         2024
# Dependencies: tidyverse
# ==============================================================================

# --- Setup --------------------------------------------------------------------
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# ==============================================================================
# PART 1: INPUT DATA
# ==============================================================================

# Site-level cost and enrolment data.
# Total costs in EUR (2023 prices) from programme accounting records.
# Enrolled families from site-level data collection (Brunnberg, 2024).
# Payer split reflects the proportion borne by the regional health authority
# versus the municipality; 50/50 is assumed unless site records indicate
# otherwise.

site_data <- tibble(
  Site = c("Site A", "Site B", "Site C", "Site D"),

  # Total programme cost per site (EUR, 2023 prices)
  Total_Cost_EUR = c(47152, 37335, 51936, 26683),

  # Number of enrolled families per site
  Families_Enrolled = c(49, 45, 44, 48),

  # Share of cost borne by the regional health authority
  Pct_Region = c(0.55, 0.50, 0.60, 0.50)
)

# Unit cost per enrolled family
site_data <- site_data %>%
  mutate(Cost_Per_Family = Total_Cost_EUR / Families_Enrolled)

print("--- Step 1: Unit Costs Per Site (Calculated) ---")
print(site_data)


# ==============================================================================
# PART 2: NATIONAL POPULATION (From the Effects Paper)
# ==============================================================================

# Total sample of first-born children (2012-2022) from national register data
total_sample_size <- 157844
years_of_data <- 11

# Average annual eligible cohort
N_eligible <- round(total_sample_size / years_of_data)
print(paste("Annual Eligible Births (Estimated):", N_eligible))


# ==============================================================================
# PART 3: MAIN SCENARIO ANALYSIS
# ==============================================================================

# Three scenarios reflecting different delivery models:
#   1. Efficient  -- uses the lowest observed unit cost (Site D, integrated model)
#   2. Base Case  -- average unit cost across all four sites
#   3. Intensive  -- uses the highest observed unit cost (Site C, dispersed model)

uptake_base <- 0.70  # assumed 70% uptake among eligible families

bia_scenarios <- tibble(
  Scenario = c(
    "1. Efficient (Integrated)",
    "2. Base Case (Average)",
    "3. Intensive (Dispersed)"
  ),

  Unit_Cost = c(
    site_data$Cost_Per_Family[site_data$Site == "Site D"],
    mean(site_data$Cost_Per_Family),
    site_data$Cost_Per_Family[site_data$Site == "Site C"]
  ),

  Region_Share = c(
    site_data$Pct_Region[site_data$Site == "Site D"],
    mean(site_data$Pct_Region),
    site_data$Pct_Region[site_data$Site == "Site C"]
  )
)

# Calculate total budget and payer-specific costs
bia_results <- bia_scenarios %>%
  mutate(
    Target_Families = round(N_eligible * uptake_base),

    # Total annual budget (million EUR)
    Total_Budget_M_EUR = round((Target_Families * Unit_Cost) / 1e6, 2),

    # Split by payer
    Region_M_EUR = round(Total_Budget_M_EUR * Region_Share, 2),
    Muni_M_EUR   = round(Total_Budget_M_EUR * (1 - Region_Share), 2)
  )

print("=== FINAL TABLE: MAIN BUDGET IMPACT (Million EUR) ===")
print(bia_results)


# ==============================================================================
# PART 4: SENSITIVITY ANALYSIS (Uptake and Payer Split)
# ==============================================================================

# Grid of uptake rates and region/municipality cost-sharing proportions.
# Uses the base-case (average) unit cost throughout.

sensitivity_grid <- expand_grid(
  Uptake       = c(0.60, 0.70, 0.80, 0.90),
  Region_Split = c(0.40, 0.50, 0.60)
)

sensitivity_results <- sensitivity_grid %>%
  mutate(
    Scenario = "Base Case (Average Cost)",
    Unit_Cost = mean(site_data$Cost_Per_Family),

    Total_M_EUR  = round((N_eligible * Uptake * Unit_Cost) / 1e6, 2),
    Region_M_EUR = round(Total_M_EUR * Region_Split, 2),
    Muni_M_EUR   = round(Total_M_EUR * (1 - Region_Split), 2)
  )

print("=== SENSITIVITY ANALYSIS TABLE ===")
print(sensitivity_results)
