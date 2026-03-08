# Analysis Code for Doctoral Thesis

**Reaching the Unreached: Targeting, Effectiveness, and Costs of an Extended Home Visiting Programme in Sweden**

Sergio Flores
Department of Public Health and Caring Sciences, Uppsala University
2026

[![DOI](https://zenodo.org/badge/DOI/XXXXX.svg)](https://doi.org/XXXXX)

## Overview

This repository contains the analysis code for a doctoral thesis examining the Extended Home Visiting Programme (EHVP) in Sweden. The thesis comprises four studies:

| Paper | Title | Analysis |
|-------|-------|----------|
| I | Targeting disadvantage in universal systems: a register-based study of Care Need Index validity for predicting paediatric health burden | `paper-i-targeting/` |
| II | Reaching the unreached: a qualitative study of families' and professionals' experiences of an extended home visiting programme | *(Qualitative; no statistical code)* |
| III | Effects of an extended home visiting programme on healthcare utilisation during the first year of life | `paper-iii-effects/` |
| IV | Scaling an equity-oriented extended home visiting programme in Sweden: costs and cross-sectoral fiscal distribution | `paper-iv-costing/` |

## Data availability

The analyses use Swedish population register data linked across multiple agencies (Statistics Sweden, the National Board of Health and Welfare, the Social Insurance Agency). Individual-level data cannot be shared publicly due to Swedish data protection regulations (GDPR and the Swedish Public Access to Information and Secrecy Act).

Researchers wishing to replicate these analyses may apply for data access through:
- **Statistics Sweden (SCB)**: [www.scb.se](https://www.scb.se)
- **National Board of Health and Welfare (Socialstyrelsen)**: [www.socialstyrelsen.se](https://www.socialstyrelsen.se)

Data access requires ethical approval from the Swedish Ethical Review Authority (Etikprovningsmyndigheten). The ethical approvals for this thesis are referenced in the individual papers.

## Repository structure

```
ehvp-thesis-code/
├── README.md
├── LICENSE
├── .gitignore
├── paper-i-targeting/
│   └── code/
│       └── cni_validation.R        # CNI validity analysis
├── paper-iii-effects/
│   └── code/
│       ├── config.R                # Paths, thresholds, specifications
│       ├── 00_run_all.R            # Pipeline orchestration
│       ├── outcomes_config.R       # Outcome variable definitions
│       ├── 01_matching.R           # DeSO-level propensity score matching
│       ├── 02_analysis.R           # Stacked TWFE estimation
│       └── helper_functions.R      # Plotting, extraction, pre-trends
└── paper-iv-costing/
    └── code/
        └── costing_analysis.R      # Budget impact and scaling scenarios
```

## Software requirements

All analyses were conducted in R (version 4.3+). Key packages:

| Package | Version | Purpose |
|---------|---------|---------|
| `data.table` | >= 1.14 | Data manipulation |
| `fixest` | >= 0.11 | Two-way fixed effects estimation |
| `MatchIt` | >= 4.5 | Propensity score matching |
| `cobalt` | >= 4.5 | Balance diagnostics |
| `did` | >= 2.1 | Callaway-Sant'Anna estimator |
| `ggplot2` | >= 3.4 | Visualisation |
| `sf` | >= 1.0 | Spatial data handling |
| `modelsummary` | >= 1.4 | Regression tables |
| `tidyverse` | >= 2.0 | General data manipulation |

Install all dependencies:

```r
install.packages(c(
  "data.table", "fixest", "MatchIt", "cobalt", "did",
  "ggplot2", "sf", "modelsummary", "tidyverse",
  "patchwork", "ggtext", "kableExtra", "scales",
  "broom", "stringr", "forcats", "openxlsx", "here"
))
```

## Reproducing the analyses

### Paper I (Targeting)

```r
# Set working directory to paper-i-targeting/
# Edit the data path at the top of the script
source("code/cni_validation.R")
```

### Paper III (Effects)

```r
# Set working directory to paper-iii-effects/
# 1. Edit paths in code/config.R
# 2. Run the pipeline:
source("code/00_run_all.R")
```

The pipeline runs in three stages: (1) propensity score matching, (2) main analysis at the 60% treatment threshold, and (3) sensitivity analyses at 70% and 80% thresholds.

### Paper IV (Costing)

```r
# Set working directory to paper-iv-costing/
# The costing script uses aggregate data embedded in the code
# (no register data required)
source("code/costing_analysis.R")
```

## Citation

If you use this code, please cite:

```
Flores, S. (2026). Reaching the Unreached: Targeting, Effectiveness, and Costs of
an Extended Home Visiting Programme in Sweden [Doctoral thesis, Uppsala University].
```

## Contact

Sergio Flores
Department of Public Health and Caring Sciences
Uppsala University
sergio.flores@pubcare.uu.se

## License

MIT License. See [LICENSE](LICENSE) for details.
