#------------------------------------------------------------------------------#
# helper_functions.R
#
# Description: Helper functions for difference-in-differences event-study
#              analysis of child healthcare utilisation outcomes.
#
# Paper:       "Effects of the Extended Home Visiting Programme on child
#               healthcare utilisation: A difference-in-differences study"
#               (Paper III, doctoral thesis)
#
# Author:      Sergio Flores
# Date:        2026
#
# Dependencies: data.table, ggplot2, fixest, did, patchwork, ggtext
#
# Contents:
#   1)  collapse_outcomes_reltime   -- reshape panel to wide by relative time
#   2)  make_event_study_plot       -- multi-series event-study plot
#       plot_eventstudy_per_cohort  -- faceted per-cohort event-study plot
#   3)  extract_twfe_terms          -- extract event-time coefficients (TWFE)
#       extract_stacked_model_terms -- extract event-time coefficients (stacked)
#       extract_sunab_terms         -- extract event-time coefficients (SUNAB)
#   4)  agg_eventstudy              -- cohort-weighted aggregation (delta method)
#   5)  sig_stars, add_sig_from_se, pvals_from_fixest, merge_fixest_pvals
#   6)  make_outcome_variants, build_rhs, build_fml
#   7)  .strip_for_grid, .rel_order_from, make_master_grid
#   8)  save_master
#   9)  csa_err_handler
#   10) csa_agg_dynamic_many, csa_extract_dynamic_results,
#       csa_weight_dynamic_results
#   11) add_rel_summaries
#   12) normalize_plot_dt, series_from_RESULT_cell, label_from_rel_key,
#       make_model_covariate_compare_grid, make_per_cohort_twfe_grid,
#       rebuild_PLOTS_from_RESULTS
#   13) save_patchwork
#   14) make_cohort_intersection
#   15) test_pre_trends, compile_pre_trends_summary
#------------------------------------------------------------------------------#


# ==== 1) Collapse panel data to one row per individual ====

collapse_outcomes_reltime <- function(
    DT,
    outcomes = c("inpatient_avoidable_events",
                 "outpatient_avoidable_events",
                 "inpatient_emergency_events",
                 "outpatient_emergency_events",
                 "inpatient_total_los_days"),
    horizon = 0:4,
    id = "lopnr",
    year = "year",
    anchor = "child_birth_year",
    keep_anchor_only = TRUE,
    order_rel_right = TRUE
) {
  # Work on a copy
  if (!is.data.table(DT)) DT <- as.data.table(DT) else DT <- copy(DT)

  # Normalize names (allow bare names or strings)
  to_name <- function(x) if (is.character(x)) x else deparse(substitute(x))
  id     <- to_name(id);  year <- to_name(year);  anchor <- to_name(anchor)

  # Validate and map to numeric indices
  nm <- names(DT)
  idx_id  <- match(id, nm)
  idx_y   <- match(year, nm)
  idx_a   <- match(anchor, nm)
  idx_out <- match(outcomes, nm)

  miss <- c(setNames(idx_id, id),
            setNames(idx_y, year),
            setNames(idx_a, anchor),
            setNames(idx_out, outcomes))
  miss <- names(miss[is.na(miss)])
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  # Relative time
  DT[, rel_time := .SD[[1]] - .SD[[2]], .SDcols = c(idx_y, idx_a)]

  # Build long only for desired horizon
  DT_sub <- DT[rel_time %in% horizon, c(id, "rel_time", outcomes), with = FALSE]

  long <- melt(DT_sub,
               id.vars = c(id, "rel_time"),
               measure.vars = outcomes,
               variable.name = "outcome",
               value.name = "val")

  long[, col := sprintf("%s_rel_time%s", outcome, rel_time)]

  # Wide per id
  wide <- dcast(long, as.formula(sprintf("%s ~ col", id)), value.var = "val")

  # Join back to all rows, then optionally keep only anchor row
  out <- wide[DT, on = id]
  out[, rel_time := NULL]

  if (keep_anchor_only) {
    nmo <- names(out)
    iy  <- match(year,   nmo)
    ia  <- match(anchor, nmo)
    out <- out[out[[iy]] == out[[ia]]]
  }

  # Move rel-time columns to the right
  if (order_rel_right) {
    desired_rel <- unlist(lapply(outcomes, function(o) paste0(o, "_rel_time", horizon)))
    rel_cols    <- intersect(desired_rel, names(out))
    other_cols  <- setdiff(names(out), rel_cols)
    setcolorder(out, c(other_cols, rel_cols))
  }

  out[]
}


# ==== 2) Plotting functions ====

make_event_study_plot <- function(
    dts,
    label,
    note = NULL,
    series = NULL,
    title_prefix = "Event-study:",
    hollow_n = 2,            # Ignored if hollow_ticks is not NULL
    connect = FALSE,
    line_width = 0.3,
    point_size = 2.8,
    errorbar_width = 0.15,
    x_shift_width = 0.12,
    x_offsets = NULL,
    legend_title = "Models",
    palette_fill = NULL,
    palette_line = NULL,
    base_size = 12,
    star_size = 6,
    star_pad  = 0.03,
    hollow_ticks = c(-8, -7, 2, 3)
) {
  # Accept single dataset
  if (inherits(dts, c("data.frame","data.table"))) {
    dts <- list(dts)
    if (is.null(series)) series <- "Series 1"
  }
  stopifnot(is.list(dts), length(dts) >= 1)

  # Series names
  if (is.null(series)) {
    series <- names(dts)
    if (is.null(series)) series <- paste("Series", seq_along(dts))
  } else {
    stopifnot(length(series) == length(dts))
  }

  # Bind into one DT with series id (keep stars/p.value if present)
  L <- lapply(dts, function(x) data.table::as.data.table(x))
  for (i in seq_along(L)) {
    # Coerce safely and repair names
    df <- as.data.frame(L[[i]], stringsAsFactors = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    data.table::setDT(df)

    # Required + optional columns
    required <- c("event_time","estimate","conf.low","conf.high")
    if (!all(required %in% names(df))) {
      stop(sprintf(
        "Input series %d (%s) is missing required columns.\nHas: %s\nNeeds: %s",
        i,
        ifelse(length(series) >= i && !is.null(series[i]), series[i], "<unnamed>"),
        paste(names(df), collapse = ", "),
        paste(required, collapse = ", ")
      ))
    }
    optional <- intersect(c("stars","p.value"), names(df))
    select_cols <- c(required, optional)

    # Subset using with=FALSE
    df <- df[, select_cols, with = FALSE]

    # Row sanity
    df <- df[is.finite(event_time) & is.finite(estimate)]

    # Ensure optional cols exist
    if (!("stars" %in% names(df)))   df[, stars   := NA_character_]
    if (!("p.value" %in% names(df))) df[, p.value := NA_real_]

    # Ensure a filled point at event_time = -1 if missing
    if (!(-1 %in% df$event_time)) {
      df <- data.table::rbindlist(list(
        df,
        data.table::data.table(
          event_time = -1L,
          estimate   = 0,
          conf.low   = 0,
          conf.high  = 0,
          stars      = "",
          p.value    = as.numeric(NA)
        )
      ), use.names = TRUE, fill = TRUE)
    }

    # Attach series label
    df[, series := series[i]]

    L[[i]] <- df
  }

  DT <- data.table::rbindlist(L, use.names = TRUE, fill = TRUE)

  bad_names <- which(!nzchar(names(DT)) | is.na(names(DT)))
  if (length(bad_names)) {
    stop(sprintf("Combined DT still has empty names at positions: %s",
                 paste(bad_names, collapse = ", ")))
  }

  # Union of event times for breaks + hollow set
  evs <- sort(unique(DT$event_time))
  if (!is.null(hollow_ticks)) {
    hollow_set <- intersect(hollow_ticks, evs)
  } else {
    hollow_set <- if (length(evs) >= 2 * hollow_n && hollow_n > 0)
      c(head(evs, hollow_n), tail(evs, hollow_n)) else numeric(0)
  }
  DT[, hollow := event_time %in% hollow_set]

  # Compute per-series x-offsets
  k <- length(series)
  if (is.null(x_offsets)) {
    x_offsets <- if (k == 1) 0 else seq(-x_shift_width, x_shift_width, length.out = k)
  } else {
    stopifnot(length(x_offsets) == k)
  }
  offset_map <- setNames(x_offsets, series)
  DT[, x_shift := event_time + offset_map[series]]

  # Colours
  if (is.null(palette_fill)) palette_fill <- c("#4C78A8","#F58518","#54A24B","#72B7B2","#FF9DA6","#E45756")
  if (is.null(palette_line)) palette_line <- c("#2F5D8A","#B05E12","#3A7A3F","#47807B","#B35A62","#B23A3F")
  palette_fill <- rep(palette_fill, length.out = k); names(palette_fill) <- series
  palette_line <- rep(palette_line, length.out = k); names(palette_line) <- series
  series_levels <- unique(DT$series)

  # Best-star logic (one star per event_time, no x-shift)
  norm_star <- function(s) {
    s <- ifelse(is.na(s), "", s)
    gsub("\\.", "\u00b7", s)
  }
  DT[, stars := norm_star(stars)]

  fill_idx <- (DT$stars == "" | is.na(DT$stars)) & is.finite(DT$p.value)
  if (any(fill_idx)) {
    DT$stars[fill_idx] <- ifelse(DT$p.value[fill_idx] <= 0.001, "***",
                                 ifelse(DT$p.value[fill_idx] <= 0.01,  "**",
                                        ifelse(DT$p.value[fill_idx] <= 0.05,  "*",
                                               ifelse(DT$p.value[fill_idx] <= 0.10,  "\u00b7", ""))))
  }

  levels_star <- c("", "\u00b7", "*", "**", "***")
  DT[, star_rank := match(stars, levels_star, nomatch = 1L) - 1L]

  Ymax <- DT[, .(y_top = max(conf.high, na.rm = TRUE)), by = event_time]
  rng  <- range(c(DT$conf.low, DT$conf.high), na.rm = TRUE)
  pad  <- diff(rng) * star_pad
  best_star <- DT[, .SD[which.max(star_rank)][1L], by = event_time, .SDcols = c("star_rank","stars")]
  STAR <- merge(best_star[, .(event_time, stars, star_rank)],
                Ymax[, .(event_time, y_top)], by = "event_time", all.x = TRUE)
  STAR <- STAR[is.finite(y_top) & star_rank > 0]
  STAR[, y_star := y_top + pad]

  subtitle_text <- if (!is.null(note) && !is.na(note) && nzchar(note)) note else NULL

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "red") +
    ggplot2::geom_errorbar(
      data = DT,
      ggplot2::aes(x = x_shift, ymin = conf.low, ymax = conf.high, color = series),
      width = errorbar_width
    ) +
    ggplot2::geom_point(
      data = DT[hollow == FALSE],
      ggplot2::aes(x = x_shift, y = estimate, color = series, fill = series),
      shape = 21, size = point_size, stroke = 0.9
    ) +
    ggplot2::geom_point(
      data = DT[hollow == TRUE],
      ggplot2::aes(x = x_shift, y = estimate, color = series),
      shape = 21, size = point_size, fill = "white", stroke = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = STAR,
      ggplot2::aes(x = event_time, y = y_star, label = stars),
      size = star_size, vjust = 0, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = evs) +
    ggplot2::scale_color_manual(
      values = palette_line[series_levels],
      breaks = series_levels,
      name   = legend_title
    ) +
    ggplot2::scale_fill_manual(
      values = palette_fill[series_levels],
      breaks = series_levels,
      guide  = "none"
    ) +
    ggplot2::labs(
      title    = sprintf("%s %s", title_prefix, label),
      subtitle = subtitle_text,
      x = "Birth cohort (event time)",
      y = "Estimate (95% CI)"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.border     = ggplot2::element_rect(color = "grey20", fill = NA, linewidth = 0.6),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      axis.ticks       = ggplot2::element_line(color = "grey40"),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.margin      = ggplot2::margin(t = 10, r = 5, b = 5, l = 5)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(shape = 21, fill = unname(palette_fill[series_levels]))
    )) +
    ggplot2::coord_cartesian(clip = "off")

  if (isTRUE(connect)) {
    p <- p + ggplot2::geom_line(
      data = DT[order(series, event_time)],
      ggplot2::aes(x = x_shift, y = estimate, color = series),
      linewidth = line_width
    )
  }

  p
}

plot_eventstudy_per_cohort <- function(dt, label, note) {
  x_breaks <- sort(unique(dt$event_time))
  yr  <- range(dt$conf.low, dt$conf.high, na.rm = TRUE)
  pad <- diff(yr) * 0.04
  y_lims <- c(yr[1] - pad, yr[2] + pad)

  ggplot(dt, aes(event_time, estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "red") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.18, linewidth = 0.4) +
    facet_wrap(~ cohort, nrow = 1, scales = "fixed") +
    scale_x_continuous(breaks = x_breaks, minor_breaks = NULL,
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = sprintf("Event-study: %s", label),
         subtitle = note, x = "Birth cohort (event time)", y = "Estimate (95% CI)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12, margin = margin(t = 6)),
      axis.ticks.length.x = unit(7, "pt"),
      axis.ticks.x = element_line(linewidth = 0.5),
      strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12, margin = margin(t = 6, b = 6)),
      panel.spacing.x = unit(16, "pt"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(t = 16, r = 14, b = 18, l = 14),
      panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),
      plot.title = element_text(face = "bold")
    )
}

# ==== 3) Convenience extractor functions for event-time terms ====

# Keep ONLY dynamic event-time terms from TWFE-style specs
extract_twfe_terms <- function(m, cohort_tag = NULL) {
  df <- data.frame(
    term      = names(coef(m)),
    estimate  = unname(coef(m)),
    std.error = unname(se(m)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  data.table::setDT(df)

  # Parse event time from patterns like "...::-3:..."
  df[, event_time := suppressWarnings(
    as.integer(sub("^.*::\\s*(-?\\d+)\\s*(?::.*)?$", "\\1", term, perl = TRUE))
  )]

  # Keep only rows with a parsed event_time (drops covariates/intercept)
  df <- df[!is.na(event_time)]

  # CIs after filtering
  df[, `:=`(
    conf.low  = estimate - qnorm(0.975) * std.error,
    conf.high = estimate + qnorm(0.975) * std.error
  )]

  if (!is.null(cohort_tag)) df[, cohort := cohort_tag]
  df[]
}

# Keep ONLY event-time terms from the stacked model specs
extract_stacked_model_terms <- function(model, alpha = 0.05) {
  stopifnot(inherits(model, "fixest"))
  z <- qnorm(1 - alpha/2)

  dt <- data.table::data.table(
    term      = names(coef(model)),
    estimate  = unname(coef(model)),
    std.error = unname(se(model))
  )

  # Keep only event-study terms: "event_time::k:coh_trt::<cohort>.<group>"
  dt <- dt[grepl("^event_time::", term)]

  # Parse fields
  dt[, event_time := as.integer(sub("^.*?event_time::(-?\\d+).*?$", "\\1", term))]
  dt[, cohort     := sub("^.*?coh_trt::(.*)\\..*$", "\\1", term)]

  # CIs
  dt[, `:=`(
    conf.low  = estimate - z * std.error,
    conf.high = estimate + z * std.error
  )]

  data.table::setorder(dt, cohort, event_time)
  dt[, .(cohort, event_time, estimate, std.error, conf.low, conf.high, term)]
}


# Keep ONLY dynamic event-time terms from SUNAB-style specs
extract_sunab_terms <- function(m) {
  df <- data.frame(
    term      = names(coef(m)),
    estimate  = unname(coef(m)),
    std.error = unname(se(m)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  data.table::setDT(df)

  # SUNAB terms often end with ':: <k>' -- be permissive but require an event-time integer
  df[, event_time := suppressWarnings(
    as.integer(sub("^.*::\\s*(-?\\d+)\\s*(?:\\(.*\\))?\\s*$", "\\1", term, perl = TRUE))
  )]

  # Drop everything without an event_time (covariates, static, etc.)
  df <- df[!is.na(event_time)]

  # CIs after filtering
  df[, `:=`(
    conf.low  = estimate - qnorm(0.975) * std.error,
    conf.high = estimate + qnorm(0.975) * std.error
  )]

  df[]
}

# ==== 4) Cohort-weighted aggregation of event-study estimates ====
#
# Aggregates cohort-specific event-study estimates into a single set of
# event-time coefficients using treated-cohort-size weights.
# Variance computed via delta method assuming cross-cohort independence.

agg_eventstudy <- function(es_dt, N_by_cohort) {

  dt <- data.table::copy(data.table::as.data.table(es_dt))

  # Weighted aggregate over cohorts
  data.table::setorder(dt, event_time)

  # Ensure cohort column types match before join
  N_by_cohort <- data.table::as.data.table(N_by_cohort)
  if ("cohort" %in% names(dt) && "cohort" %in% names(N_by_cohort)) {
    dt[, cohort := as.character(cohort)]
    N_by_cohort[, cohort := as.character(cohort)]
  }

  dt[N_by_cohort, N_treated := i.N_treated, on = .(cohort)]
  dt[, denom := sum(N_treated), by = event_time]
  dt[, w     := N_treated / denom]

  out <- dt[
    , .(
      estimate = sum(w * estimate),
      var_hat  = sum((w^2) * (std.error^2))  # delta-method assuming independence
    ),
    by = event_time
  ][order(event_time)]

  out[, se := sqrt(var_hat)]
  out[, `:=`(
    conf.low  = estimate - qnorm(0.975) * se,
    conf.high = estimate + qnorm(0.975) * se
  )]

  out[]
}

# ==== 5) Significance stars and p-value helpers ====

# p-values -> star symbols
sig_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01,  "**",
                       ifelse(p <= 0.05, "*",
                              ifelse(p <= 0.10, "\u00b7", "")
                       ))))
}

# Add p-values and stars from estimate + se columns.
# Uses normal approximation by default (z = beta/se). Set use_t = TRUE with df for t.
add_sig_from_se <- function(DT, est_col = "estimate", se_col = "std.error", df = Inf){
  stopifnot(est_col %in% names(DT), se_col %in% names(DT))
  z <- DT[[est_col]] / DT[[se_col]]
  p <- 2 * pnorm(-abs(z))
  DT$p.value <- p
  DT$stars   <- sig_stars(p)
  DT$sig     <- DT$p.value <= 0.05
  DT
}

# Grab clustered p-values from fixest directly (more faithful than z-approx)
pvals_from_fixest <- function(mod){
  ct <- as.data.frame(fixest::coeftable(mod))
  ct$term <- rownames(ct)
  data.table::setDT(ct)[, .(term, p.value = `Pr(>|t|)`)]
}

# Left-join p-values from fixest coeftable; fallback to z-approx for missing rows
merge_fixest_pvals <- function(DT, mod, est_col = "estimate", se_col = "std.error"){
  data.table::setDT(DT)
  pv <- pvals_from_fixest(mod)
  out <- pv[DT, on = .(term), nomatch = 0L][
    DT, on = .(term), roll = TRUE]
  # Fill NA p.values from estimate/se where terms were relabelled
  fill <- is.na(out$p.value)
  if (any(fill)){
    z  <- out[[est_col]][fill] / out[[se_col]][fill]
    out$p.value[fill] <- 2 * pnorm(-abs(z))
  }
  out[, `:=`(stars = sig_stars(p.value), sig = p.value <= 0.05)]
  out[]
}

# ==== 6) Outcome variants and formula construction ====

# Build the set of outcome variants (single-k + cumulative sums) for a given base outcome.
# Labels indicate CHILD AGE at outcome measurement (not event time / birth cohort).
make_outcome_variants <- function(y_base, label_base, horizon_k = 0:3) {

  # Age-specific labels
  age_labels <- c(
    "0" = "first year of life",
    "1" = "age 1",
    "2" = "age 2",
    "3" = "age 3",
    "4" = "age 4"
  )

  # Single-k variants (outcomes at specific child ages)
  k_list <- lapply(horizon_k, function(k) {
    age_desc <- if (as.character(k) %in% names(age_labels)) age_labels[as.character(k)] else sprintf("age %d", k)
    list(
      y       = sprintf("%s_rel_time%d", y_base, k),
      rel_key = sprintf("age%d", k),
      label   = sprintf("%s (%s)", label_base, age_desc)
    )
  })
  # sum(0,1) -- first two years of life
  sum01 <- list(
    y       = sprintf("%s_rel_time01_sum", y_base),
    rel_key = "age01sum",
    label   = sprintf("%s (cumulative: ages 0-1)", label_base)
  )
  # sum(0..3) -- first four years of life
  sum03 <- list(
    y       = sprintf("%s_rel_time0_3_sum", y_base),
    rel_key = "age0_3sum",
    label   = sprintf("%s (cumulative: ages 0-3)", label_base)
  )
  c(k_list, list(sum01, sum03))
}


build_rhs <- function(treat, covar_str) paste(c(treat, covar_str), collapse = " + ")

build_fml <- function(y, rhs, fe) as.formula(sprintf("`%s` ~ %s | %s", y, rhs, fe))

# ==== 7) Master grid plot construction ====

# Strip titles for embedding in a grid
.strip_for_grid <- function(p, collapse_shape_fill = TRUE) {
  p <- p + labs(title = NULL, subtitle = NULL)
  if (collapse_shape_fill) {
    p <- p + guides(
      shape = "none",
      alpha = "none",
      size  = "none",
      fill  = "none"
    )
  }
  p
}

# Detect rel-key ordering (supports both "age0" and legacy "rel0" formats)
.rel_order_from <- function(PLOTS_branch) {
  rel_keys <- names(PLOTS_branch)

  # Try age-based keys first
  age_keys <- intersect(rel_keys, c("age0", "age1", "age2", "age3"))
  if (length(age_keys) > 0) {
    ks <- sort(as.integer(sub("^age", "", age_keys)))
    return(list(keys = paste0("age", ks), ks = ks))
  }

  # Fall back to legacy rel-based keys
  rel_numeric <- intersect(rel_keys, c("rel0", "rel1", "rel2", "rel3"))
  if (length(rel_numeric) > 0) {
    ks <- sort(as.integer(sub("^rel", "", rel_numeric)))
    return(list(keys = paste0("rel", ks), ks = ks))
  }

  list(keys = character(0), ks = integer(0))
}

# Assemble a master grid of event-study plots
make_master_grid <- function(PLOTS, base_outcomes, pretty_labels,
                             covars_on = TRUE, plot_name = "cs",
                             grid_title = NULL,
                             col_header_size = 12,
                             main_title_size = 18,
                             legend_box_size = 0.3,
                             subtitle_fill = "#E8F1FD",
                             subtitle_box_col = "grey60",
                             subtitle_box_size = 0.3) {

  cov_key <- if (covars_on) "covars_on" else "covars_off"
  stopifnot(cov_key %in% names(PLOTS))

  ro <- .rel_order_from(PLOTS[[cov_key]])
  rel_keys <- ro$keys
  rel_nums <- ro$ks

  first_rel <- rel_keys[[1]]
  available_outcomes <- intersect(base_outcomes, names(PLOTS[[cov_key]][[first_rel]]))
  stopifnot(length(available_outcomes) > 0)

  grid_plots <- list()
  for (oi in seq_along(available_outcomes)) {
    outcome <- available_outcomes[[oi]]
    pretty  <- pretty_labels[match(outcome, base_outcomes)]

    for (rj in seq_along(rel_keys)) {
      rk <- rel_keys[[rj]]
      p <- PLOTS[[cov_key]][[rk]][[outcome]][[plot_name]]
      if (is.null(p)) p <- ggplot() + theme_void()

      p <- .strip_for_grid(p)

      # Column headers: "Outcome in year k"
      if (oi == 1) {
        p <- p +
          ggtitle(sprintf("Outcome in year %d", rel_nums[[rj]])) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold", size = col_header_size))
      }
      # Row subtitles (boxed header via ggtext)
      if (rj == 1) {
        p <- p +
          labs(subtitle = pretty) +
          theme(
            plot.subtitle = ggtext::element_textbox_simple(
              padding   = margin(4, 6, 4, 6),
              margin    = margin(b = 6),
              fill      = subtitle_fill,
              box.colour= subtitle_box_col,
              linewidth = subtitle_box_size,
              r         = unit(2, "pt"),
              size      = 11,
              face      = "bold"
            )
          )
      }

      grid_plots[[length(grid_plots) + 1L]] <- p
    }
  }

  pw <- wrap_plots(
    grid_plots,
    ncol   = length(rel_keys),
    guides = "collect",
    byrow  = TRUE
  )

  pw <- pw + plot_annotation(
    title = grid_title %||% if (covars_on) "With covariates" else "No covariates",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = main_title_size),

      legend.position       = "bottom",
      legend.direction      = "horizontal",
      legend.box            = "horizontal",
      legend.box.background = element_rect(color = "grey40", fill = NA, linewidth = legend_box_size),
      legend.background     = element_rect(fill = "white", color = NA),
      legend.key.width      = unit(16, "pt"),
      legend.key.height     = unit(10, "pt"),
      legend.spacing.x      = unit(6, "pt")
    )
  )

  pw <- pw & guides(
    color    = guide_legend(
      nrow = 1,
      override.aes = list(shape = 16, alpha = 1, size = 3, stroke = 0)
    ),
    linetype = guide_legend(nrow = 1)
  )
}

# ==== 8) Save function for master plots ====

save_master <- function(plot_obj, PLOTS, covars_on = TRUE,
                        cell_w = 3.2,      # inches per column
                        cell_h = 2.6,      # inches per row
                        extra_w = 0.6,     # padding for margins
                        extra_h = 1.4,     # padding for title + legend
                        file_pdf = "master.pdf",
                        file_png = NULL,
                        dpi = 600) {

  cov_key <- if (covars_on) "covars_on" else "covars_off"
  rel_keys <- names(PLOTS[[cov_key]])
  n_cols <- length(rel_keys)

  first_rel <- rel_keys[[1]]
  n_rows <- length(names(PLOTS[[cov_key]][[first_rel]]))

  width_in  <- n_cols * cell_w + extra_w
  height_in <- n_rows * cell_h + extra_h

  # PDF (vector)
  ggsave(file_pdf, plot = plot_obj, device = cairo_pdf,
         width = width_in, height = height_in, units = "in")

  # Optional hi-res PNG
  if (!is.null(file_png)) {
    ggsave(file_png, plot = plot_obj, dpi = dpi,
           width = width_in, height = height_in, units = "in", bg = "white")
  }

  invisible(list(width_in = width_in, height_in = height_in))
}

# ==== 9) Error handler for Callaway-Sant'Anna sub-loop ====

# Error handler factory: logs into iter_error_log in the caller env and returns NULL
csa_err_handler <- function(covars_on, rel_key, y, y_base, g, log_env, log_name = "iter_error_log") {
  force(covars_on); force(rel_key); force(y); force(y_base); force(g)
  stopifnot(is.environment(log_env))
  function(e) {
    msg <- conditionMessage(e)
    message(sprintf("[CSA SKIP] covars=%s, %s, y=%s, g=%s: %s",
                    covars_on, rel_key, y, g, msg))
    new_row <- data.table(
      covars_on    = covars_on,
      rel_key      = rel_key,
      outcome      = y,
      outcome_base = y_base,
      cohort       = as.integer(g),
      msg          = msg
    )
    log_env[[log_name]] <- rbind(log_env[[log_name]], new_row, fill = TRUE)
    NULL
  }
}

# ==== 10) Callaway-Sant'Anna extraction and aggregation ====

#------------------------------------------------------------------------------#
# csa_agg_dynamic_many
# Aggregates dynamic effects for multiple att_gt objects.
# Returns: $ok (named list of successful aggte objects), $err (error log DT)
#------------------------------------------------------------------------------#
csa_agg_dynamic_many <- function(results_ok, covars_on, rel_key, y, y_base) {
  stopifnot(is.list(results_ok))
  err_log <- data.table(
    covars_on    = character(),
    rel_key      = character(),
    outcome      = character(),
    outcome_base = character(),
    cohort       = integer(),
    msg          = character()
  )
  ok <- list()

  for (g in names(results_ok)) {
    es <- tryCatch(
      did::aggte(results_ok[[g]], type = "dynamic"),
      error = function(e) {
        msg <- paste0("aggte failed: ", conditionMessage(e))
        message(sprintf("[CSA aggte SKIP] covars=%s, %s, y=%s, g=%s: %s",
                        covars_on, rel_key, y, g, msg))
        err_log <<- rbind(err_log, data.table(
          covars_on    = covars_on,
          rel_key      = rel_key,
          outcome      = y,
          outcome_base = y_base,
          cohort       = as.integer(g),
          msg          = msg
        ), fill = TRUE)
        NULL
      }
    )
    if (!is.null(es)) ok[[g]] <- es
  }
  list(ok = ok, err = err_log)
}

#------------------------------------------------------------------------------#
# csa_extract_dynamic_results
# Input: named list of aggte(type="dynamic") objects
# Output: DT_stacked with one row per (cohort, event_time)
#------------------------------------------------------------------------------#
csa_extract_dynamic_results <- function(agg_results_ok) {
  if (!length(agg_results_ok)) return(NULL)
  stopifnot(is.list(agg_results_ok))
  requireNamespace("data.table")

  DT_by_cohort <- lapply(names(agg_results_ok), function(g) {
    es <- agg_results_ok[[g]]
    out <- data.table::data.table(
      cohort     = as.integer(g),
      event_time = es$egt,
      estimate   = es$att.egt,
      se         = es$se.egt
    )
    z <- stats::qnorm(0.975)
    out[, `:=`(
      conf.low  = estimate - z * se,
      conf.high = estimate + z * se
    )]
    out[]
  })
  names(DT_by_cohort) <- names(agg_results_ok)

  DT_stacked <- data.table::rbindlist(DT_by_cohort, use.names = TRUE, fill = TRUE)
  data.table::setorder(DT_stacked, cohort, event_time)
  DT_stacked
}

#------------------------------------------------------------------------------#
# csa_weight_dynamic_results
# Input:
#   DT_stacked: output from csa_extract_dynamic_results()
#   N_treated : data.table with columns: cohort, N  (treated counts per cohort)
#   ci_level  : 0.95 by default
#   keep_weights: if TRUE, also returns per-(event_time, cohort) weights
# Output:
#   data.table with event_time, estimate, se, conf.low, conf.high
#   (or list(agg=..., weights=...) if keep_weights=TRUE)
#
# Variance: delta-method under cross-cohort independence assumption.
#------------------------------------------------------------------------------#
csa_weight_dynamic_results <- function(DT_stacked, N_treated, ci_level = 0.95, keep_weights = FALSE) {
  if (is.null(DT_stacked) || !nrow(DT_stacked)) return(NULL)
  stopifnot(all(c("cohort", "event_time", "estimate", "se") %in% names(DT_stacked)))
  stopifnot(all(c("cohort", "N") %in% names(N_treated)))
  requireNamespace("data.table")

  DTw <- N_treated[DT_stacked, on = "cohort"]
  # Drop rows where N is missing (cohort not found in N_treated)
  if (anyNA(DTw$N)) {
    data.table::set(DTw, which(is.na(DTw$N)), "N", NA_real_)
    missing_coh <- sort(unique(DTw[is.na(N), cohort]))
    message(sprintf("[CSA weight] Dropping cohorts with missing N: %s",
                    paste(missing_coh, collapse = ", ")))
    DTw <- DTw[!is.na(N)]
  }
  if (!nrow(DTw)) return(NULL)

  data.table::setorder(DTw, event_time)
  DTw[, N_sum := sum(N), by = event_time]
  DTw <- DTw[N_sum > 0]
  if (!nrow(DTw)) return(NULL)

  DTw[, w := N / N_sum]

  # Aggregate across cohorts by event_time
  DT_agg <- DTw[, .(
    estimate = sum(w * estimate),
    var_hat  = sum((w^2) * (se^2))   # independence across cohorts assumption
  ), by = event_time][order(event_time)]

  z <- stats::qnorm(0.5 + ci_level/2)
  DT_agg[, `:=`(
    se        = sqrt(var_hat),
    conf.low  = estimate - z * sqrt(var_hat),
    conf.high = estimate + z * sqrt(var_hat)
  )][, var_hat := NULL][]

  if (keep_weights) {
    return(list(
      agg     = DT_agg,
      weights = DTw[, .(event_time, cohort, N, w)]
    ))
  }
  DT_agg
}


# ==== 11) Cumulative outcome summaries over relative time ====

# Add two summary variables per outcome (ages 0-1 and 0-3).
# Rule: define the sum ONLY if all required columns exist AND no row values are NA.
add_rel_summaries <- function(DT, bases = base_outcomes) {
  stopifnot(data.table::is.data.table(DT))
  for (base in bases) {
    cols2 <- paste0(base, "_rel_time", 0:1)
    cols4 <- paste0(base, "_rel_time", 0:3)

    # --- Sum 0-1 ---
    new2 <- paste0(base, "_rel_time01_sum")
    miss2 <- setdiff(cols2, names(DT))
    if (length(miss2)) {
      warning(sprintf("Missing columns for %s (0-1): %s", base, paste(miss2, collapse = ", ")))
      DT[, (new2) := NA_real_]
    } else {
      DT[, (new2) := rowSums(.SD, na.rm = FALSE), .SDcols = cols2]
    }

    # --- Sum 0-3 ---
    new4 <- paste0(base, "_rel_time0_3_sum")
    miss4 <- setdiff(cols4, names(DT))
    if (length(miss4)) {
      warning(sprintf("Missing columns for %s (0-3): %s", base, paste(miss4, collapse = ", ")))
      DT[, (new4) := NA_real_]
    } else {
      DT[, (new4) := rowSums(.SD, na.rm = FALSE), .SDcols = cols4]
    }
  }
  invisible(DT)
}


# ==== 12) Plot extraction and grid construction from RESULTS ====

# Friendly labels for legend
.default_label_map <- c(
  TWFE_weighted                  = "Weighted average of disaggregated models",
  stacked_TWFE_weighted          = "Stacked model",
  balanced_stacked_TWFE_weighted = "Stacked model (balanced panel)",
  SUNAB                          = "Sun and Abraham",
  CSA_dynamic_weighted           = "Callaway Sant'Anna (IPW)"
)

# Normalize naming/columns coming from different extractors
normalize_plot_dt <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  dt <- data.table::as.data.table(dt)
  if (!"se" %in% names(dt) && "std.error" %in% names(dt)) data.table::setnames(dt, "std.error", "se")
  dt
}

# Extract the series from one RESULTS "cell"
series_from_RESULT_cell <- function(cell, include = names(.default_label_map),
                                    label_map = .default_label_map) {
  out <- list()
  for (key in include) {
    if (!is.null(cell[[key]])) {
      dt <- normalize_plot_dt(cell[[key]])
      if (!is.null(dt) && nrow(dt)) out[[ label_map[[key]] ]] <- dt
    }
  }
  out
}

# Create label from rel_time/age variables
label_from_rel_key <- function(rel_key, pretty_label) {
  # Handle age-based keys (current format)
  if (rel_key == "age01sum")  return(sprintf("%s (cumulative: ages 0-1)", pretty_label))
  if (rel_key == "age0_3sum") return(sprintf("%s (cumulative: ages 0-3)", pretty_label))
  # Handle old rel-based keys (backward compatibility)
  if (rel_key == "rel01sum")  return(sprintf("%s - sum of years 0-1", pretty_label))
  if (rel_key == "rel0_3sum") return(sprintf("%s - sum of years 0-3", pretty_label))
  # age{K} or rel{K} where K is 0..3
  if (grepl("^age", rel_key)) {
    k <- suppressWarnings(as.integer(sub("^age", "", rel_key)))
    age_desc <- c("0" = "first year of life", "1" = "age 1", "2" = "age 2", "3" = "age 3")
    return(sprintf("%s (%s)", pretty_label, age_desc[as.character(k)]))
  }
  k <- suppressWarnings(as.integer(sub("^rel", "", rel_key)))
  sprintf("%s - %d years after treatment", pretty_label, k)
}

# Compare a single model's results (with vs without covariates)
# across all base outcomes, for a chosen rel_key.
# model_key: one of "TWFE_weighted", "stacked_TWFE_weighted",
#            "balanced_stacked_TWFE_weighted", "SUNAB", "CSA_dynamic_weighted"
# Requires: ggtext, patchwork
make_model_covariate_compare_grid <- function(
    RESULTS, base_outcomes, pretty_labels,
    rel_key,
    model_key = "SUNAB",
    ncol = 3,
    connect = FALSE,
    legend_title = "Specification",
    subtitle_fill = "#E8F1FD",
    subtitle_box_col = "grey60",
    subtitle_box_size = 0.3,
    title_size = 20,
    title_face = "bold",
    title_center = TRUE,
    legend_box_size = 0.3,
    hollow_ticks = NULL,
    hollow_n = 2
) {
  label_map <- c(
    TWFE_weighted                  = "Weighted average of disaggregated models",
    stacked_TWFE_weighted          = "Stacked model",
    balanced_stacked_TWFE_weighted = "Stacked model (balanced panel)",
    SUNAB                          = "Sun and Abraham",
    CSA_dynamic_weighted           = "Callaway Sant'Anna (IPW)"
  )
  if (!model_key %in% names(label_map)) {
    stop("Unknown model_key: ", model_key,
         ". Choose one of: ", paste(names(label_map), collapse = ", "))
  }
  model_label <- label_map[[model_key]]

  .rel_title <- function(rk) {
    # Age-based keys
    if (rk == "age01sum")  return("First two years of life")
    if (rk == "age0_3sum") return("First four years of life")
    if (grepl("^age[0-9]$", rk)) {
      kk <- as.integer(sub("^age","", rk))
      titles <- c("0" = "First year of life", "1" = "Age 1", "2" = "Age 2", "3" = "Age 3")
      return(titles[as.character(kk)])
    }
    # Legacy rel-based keys
    if (rk == "rel01sum")  return("Outcomes year 1 and 2")
    if (rk == "rel0_3sum") return("Outcomes first 4 years")
    kk <- suppressWarnings(as.integer(sub("^rel","", rk)))
    if (!is.na(kk)) return(sprintf("Outcome in year %d", kk))
    rk
  }
  grid_title <- sprintf("%s covariate comparison - %s", model_label, .rel_title(rel_key))

  plist <- vector("list", length(base_outcomes))

  for (i in seq_along(base_outcomes)) {
    y_base <- base_outcomes[i]
    pretty <- pretty_labels[i]

    on  <- try(normalize_plot_dt(RESULTS$covars_on [[rel_key]][[y_base]][[model_key]]),  silent = TRUE)
    off <- try(normalize_plot_dt(RESULTS$covars_off[[rel_key]][[y_base]][[model_key]]),  silent = TRUE)
    on  <- if (inherits(on,  "try-error")) NULL else on
    off <- if (inherits(off, "try-error")) NULL else off

    dts <- list()
    if (!is.null(on)  && nrow(on))  dts[["With covariates"]]    <- on
    if (!is.null(off) && nrow(off)) dts[["Without covariates"]] <- off
    if (!length(dts)) next

    p <- make_event_study_plot(
      dts = dts,
      label = pretty,
      note  = sprintf("%s estimated on individual-level data", model_label),
      connect = connect,
      legend_title = legend_title,
      title_prefix = "",
      hollow_ticks = hollow_ticks,
      hollow_n = hollow_n
    )

    # Light-blue header per cell
    p <- p +
      labs(title = NULL, subtitle = pretty) +
      theme(
        plot.subtitle = ggtext::element_textbox_simple(
          padding   = margin(4, 6, 4, 6),
          margin    = margin(b = 6),
          fill      = subtitle_fill,
          box.colour= subtitle_box_col,
          linewidth = subtitle_box_size,
          r         = unit(2, "pt"),
          size      = 11,
          face      = "bold"
        )
      )

    plist[[i]] <- p
  }

  plots <- Filter(Negate(is.null), plist)
  if (!length(plots)) return(NULL)

  pw <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = grid_title,
      theme = theme(
        plot.title = element_text(
          face = title_face,
          size = title_size,
          hjust = if (isTRUE(title_center)) 0.5 else 0
        ),
        legend.position       = "bottom",
        legend.direction      = "horizontal",
        legend.box            = "horizontal",
        legend.box.background = element_rect(color = "grey40", fill = NA, linewidth = legend_box_size),
        legend.background     = element_rect(fill = "white", color = NA),
        legend.key.width      = unit(16, "pt"),
        legend.key.height     = unit(10, "pt"),
        legend.spacing.x      = unit(6, "pt")
      )
    )

  pw <- pw & guides(
    color = guide_legend(nrow = 1, override.aes = list(shape = 16, size = 3, alpha = 1, stroke = 0)),
    linetype = guide_legend(nrow = 1),
    shape = "none",
    fill  = "none"
  )
  pw
}


# Per-cohort TWFE grid
# Requires: plot_eventstudy_per_cohort(), ggtext, patchwork
make_per_cohort_twfe_grid <- function(
    RESULTS, base_outcomes, pretty_labels,
    rel_key,
    covars_on = TRUE,
    ncol = 1,
    title_size = 20,
    title_face = "bold",
    title_center = TRUE,
    subtitle_fill = "#E8F1FD",
    subtitle_box_col = "grey60",
    subtitle_box_size = 0.3,
    subtitle_font_size = 11,
    subtitle_radius_pt = 2
) {
  cov_key <- if (covars_on) "covars_on" else "covars_off"

  rel_title <- switch(
    rel_key,
    # Age-based format
    age0 = "First year of life",
    age1 = "Age 1",
    age2 = "Age 2",
    age3 = "Age 3",
    age01sum  = "Cumulative: ages 0-1",
    age0_3sum = "Cumulative: ages 0-3",
    # Legacy rel-based format
    rel01sum  = "Sum of years 0-1",
    rel0_3sum = "Sum of years 0-3",
    {
      if (grepl("^age", rel_key)) {
        kk <- suppressWarnings(as.integer(sub("^age","", rel_key)))
        if (!is.na(kk)) sprintf("Age %d", kk) else rel_key
      } else {
        kk <- suppressWarnings(as.integer(sub("^rel","", rel_key)))
        if (!is.na(kk)) sprintf("Outcome in year %d", kk) else rel_key
      }
    }
  )
  grid_title <- sprintf("Per-cohort TWFE - %s (%s)",
                        rel_title,
                        if (covars_on) "with covariates" else "no covariates")

  plots <- vector("list", length(base_outcomes))
  idx <- 0L

  for (i in seq_along(base_outcomes)) {
    y_base <- base_outcomes[i]
    pretty <- pretty_labels[i]

    cell <- try(RESULTS[[cov_key]][[rel_key]][[y_base]], silent = TRUE)
    if (inherits(cell, "try-error") || is.null(cell) || is.null(cell$TWFE_per_cohort)) next

    dt <- data.table::as.data.table(cell$TWFE_per_cohort)
    if (!nrow(dt)) next

    p <- plot_eventstudy_per_cohort(dt, label = pretty, note = NULL) +
      labs(title = NULL, subtitle = pretty) +
      theme(
        plot.subtitle = ggtext::element_textbox_simple(
          padding    = margin(4, 6, 4, 6),
          margin     = margin(b = 6),
          fill       = subtitle_fill,
          box.colour = subtitle_box_col,
          linewidth  = subtitle_box_size,
          r          = unit(subtitle_radius_pt, "pt"),
          size       = subtitle_font_size,
          face       = "bold"
        )
      )

    idx <- idx + 1L
    plots[[idx]] <- p
  }

  plots <- plots[seq_len(idx)]
  if (!length(plots)) return(NULL)

  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_layout(heights = rep(1, length(plots))) +
    patchwork::plot_annotation(
      title = grid_title,
      theme = theme(
        plot.title = element_text(
          face  = title_face,
          size  = title_size,
          hjust = if (isTRUE(title_center)) 0.5 else 0
        )
      )
    )
}

# Rebuild PLOTS structure from RESULTS for use with make_master_grid
rebuild_PLOTS_from_RESULTS <- function(
    RESULTS, base_outcomes, pretty_labels,
    covars_on = TRUE,
    include   = c("TWFE_weighted","stacked_TWFE_weighted","balanced_stacked_TWFE_weighted","SUNAB","CSA_dynamic_weighted"),
    label_map = .default_label_map,
    connect   = FALSE,
    plot_name = "cs"
) {
  cov_key <- if (covars_on) "covars_on" else "covars_off"
  if (is.null(RESULTS[[cov_key]])) stop("No RESULTS for cov_key: ", cov_key)

  ro <- .rel_order_from(RESULTS[[cov_key]])
  rel_keys <- ro$keys
  if (!length(rel_keys)) {
    available_keys <- paste(names(RESULTS[[cov_key]]), collapse = ", ")
    stop("No age/rel keys (0..3) found under RESULTS[[", cov_key, "]]. Available keys: ", available_keys)
  }

  PLOTS <- list(); PLOTS[[cov_key]] <- list()

  for (rk in rel_keys) {
    PLOTS[[cov_key]][[rk]] <- list()
    for (i in seq_along(base_outcomes)) {
      y_base <- base_outcomes[i]
      pretty <- pretty_labels[i]
      cell   <- RESULTS[[cov_key]][[rk]][[y_base]]
      if (is.null(cell)) next

      dts <- series_from_RESULT_cell(cell, include = include, label_map = label_map)
      if (!length(dts)) next

      label <- label_from_rel_key(rk, pretty)
      note  <- sprintf("Estimated on individual-level data with %s",
                       if (covars_on) "covariates" else "NO covariates")

      tmp <- list()
      tmp[[plot_name]] <- make_event_study_plot(
        dts     = dts,
        label   = label,
        note    = note,
        connect = connect
      )
      PLOTS[[cov_key]][[rk]][[y_base]] <- tmp
    }
  }
  PLOTS
}

# ==== 13) Save patchwork grid plots ====

save_patchwork <- function(plot, file,
                           n_panels,
                           ncol = 3,
                           cell_w = 4.2,
                           cell_h = 3.0,
                           extra_w = 0.6,
                           extra_h = 1.2,
                           dpi = 600,
                           also_png = TRUE) {
  n_rows <- ceiling(n_panels / ncol)
  width  <- ncol  * cell_w + extra_w
  height <- n_rows * cell_h + extra_h
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  if (also_png) {
    ggsave(sub("\\.pdf$", ".png", file), plot, width = width, height = height,
           dpi = dpi, limitsize = FALSE)
  }
}


# ==== 14) Cohort-intersection indicator for summed variables ====

make_cohort_intersection <- function(dts,
                                     sum_suffix = "_rel_time0_3_sum$",
                                     label_prefix = "Cohort ") {
  stopifnot(is.list(dts), length(dts) >= 1L)

  # Build per-cohort indicator (0/1) table
  build_indic <- function(DT, cohort_name) {
    stopifnot(is.data.table(DT), "event_time" %in% names(DT))
    sum_cols <- grep(sum_suffix, names(DT), value = TRUE)

    # Count non-NA by event_time for selected summary columns
    if (length(sum_cols)) {
      counts <- DT[, lapply(.SD, function(x) sum(!is.na(x))), by = event_time, .SDcols = sum_cols]
      setorder(counts, event_time)
      # 1 if ALL summary variables have count > 0 at this event_time
      ind <- counts[, setNames(
        list(as.integer(rowSums(.SD > 0L) == length(.SD))),
        paste0(label_prefix, cohort_name)
      ), by = event_time, .SDcols = sum_cols]
    } else {
      # No summary columns found: keep event_time domain, set indicator NA
      ind <- unique(DT[, .(event_time)])
      ind[, (paste0(label_prefix, cohort_name)) := NA_integer_]
    }
    ind
  }

  indic_list <- Map(build_indic, dts, names(dts))

  # Full outer join over event_time
  out <- Reduce(function(x, y) merge(x, y, by = "event_time", all = TRUE), indic_list)

  # Replace NA with 0 (no info at that event_time => treat as absent)
  for (cn in setdiff(names(out), "event_time")) {
    set(out, which(is.na(out[[cn]])), cn, 0L)
  }
  setorder(out, event_time)
  out[]
}

# ==== 15) Pre-Trends Testing ====
#
# Formal joint hypothesis tests for the parallel trends assumption in
# difference-in-differences designs.
#
# Reference: Roth (2022) "Pretest with Caution: Event-Study Estimates after
#            Testing for Parallel Pre-Trends" AER: Insights
#------------------------------------------------------------------------------#

#' Test Pre-Trends Assumption via Joint Wald Test
#'
#' Performs a joint Wald test on pre-treatment coefficients (event_time < reference_period)
#' to formally test the parallel trends assumption. Uses fixest::wald for cluster-robust
#' inference when available.
#'
#' @param model A fixest model object with event-study specification
#' @param reference_period The omitted reference period (default: -1)
#' @param alpha Significance level for rejection (default: 0.10)
#' @param verbose Print detailed output (default: FALSE)
#'
#' @return A list with:
#'   - f_stat: Wald F-statistic
#'   - df1, df2: Degrees of freedom
#'   - p_value: P-value from F-distribution
#'   - reject_null: Logical, TRUE if parallel trends rejected at alpha
#'   - n_pre_coefs: Number of pre-treatment coefficients tested
#'   - significant_pre_coefs: Character vector of individually significant (p<0.10) pre-period terms
#'   - pre_coef_details: Data.table with coefficient details
test_pre_trends <- function(model, reference_period = -1, alpha = 0.10, verbose = FALSE) {
  stopifnot(inherits(model, "fixest"))

  coef_names <- names(coef(model))

  # Identify pre-treatment event-time coefficients
  # Pattern matches SUNAB and TWFE-style coefficient names
  pre_pattern <- "::\\s*(-?\\d+)\\s*(?:\\(.*\\))?\\s*$"

  pre_coefs <- character(0)
  for (i in seq_along(coef_names)) {
    match <- regmatches(coef_names[i], regexec(pre_pattern, coef_names[i], perl = TRUE))[[1]]
    if (length(match) >= 2) {
      et <- suppressWarnings(as.integer(match[2]))
      if (!is.na(et) && et < reference_period) {
        pre_coefs <- c(pre_coefs, coef_names[i])
      }
    }
  }

  # Return NA results if no pre-treatment coefficients
  if (length(pre_coefs) == 0) {
    if (verbose) message("No pre-treatment coefficients found for testing.")
    return(list(
      f_stat = NA_real_,
      df1 = NA_integer_,
      df2 = NA_integer_,
      p_value = NA_real_,
      reject_null = NA,
      n_pre_coefs = 0L,
      significant_pre_coefs = character(0),
      pre_coef_details = NULL
    ))
  }

  if (verbose) {
    cat(sprintf("Testing %d pre-treatment coefficients:\n", length(pre_coefs)))
    cat(paste("  ", pre_coefs, collapse = "\n"), "\n")
  }

  # Build coefficient details table
  ct <- as.data.frame(fixest::coeftable(model))
  ct$term <- rownames(ct)
  pre_ct <- ct[ct$term %in% pre_coefs, ]

  pre_details <- data.table::data.table(
    term = pre_ct$term,
    estimate = pre_ct$Estimate,
    std_error = pre_ct$`Std. Error`,
    t_value = pre_ct$`t value`,
    p_value = pre_ct$`Pr(>|t|)`
  )

  # Individually significant pre-period coefficients
  sig_pre <- pre_details[p_value < alpha, term]

  # Joint Wald test using fixest::wald
  # H0: All pre-treatment coefficients = 0 (parallel trends holds)
  wald_result <- tryCatch({
    hyp <- paste(sprintf("`%s` = 0", pre_coefs), collapse = ", ")
    fixest::wald(model, hyp, print = FALSE)
  }, error = function(e) {
    if (verbose) message("Wald test failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(wald_result) || (length(wald_result) == 1 && is.na(wald_result))) {
    return(list(
      f_stat = NA_real_,
      df1 = length(pre_coefs),
      df2 = NA_integer_,
      p_value = NA_real_,
      reject_null = NA,
      n_pre_coefs = length(pre_coefs),
      significant_pre_coefs = as.character(sig_pre),
      pre_coef_details = pre_details
    ))
  }

  # fixest::wald returns a named vector with: stat, p, df1, df2
  f_stat <- wald_result[["stat"]]
  p_value <- wald_result[["p"]]
  df1 <- wald_result[["df1"]]
  df2 <- wald_result[["df2"]]

  if (is.null(f_stat) || is.na(f_stat)) {
    return(list(
      f_stat = NA_real_,
      df1 = if (!is.null(df1) && !is.na(df1)) as.integer(df1) else length(pre_coefs),
      df2 = if (!is.null(df2) && !is.na(df2)) as.integer(df2) else NA_integer_,
      p_value = NA_real_,
      reject_null = NA,
      n_pre_coefs = length(pre_coefs),
      significant_pre_coefs = as.character(sig_pre),
      pre_coef_details = pre_details
    ))
  }

  if (verbose) {
    cat(sprintf("\nJoint Wald Test for Parallel Trends:\n"))
    cat(sprintf("  F(%d, %d) = %.3f, p = %.4f\n", df1, df2, f_stat, p_value))
    if (p_value < alpha) {
      cat(sprintf("  ** REJECT parallel trends at alpha = %.2f **\n", alpha))
    } else {
      cat(sprintf("  Fail to reject parallel trends at alpha = %.2f\n", alpha))
    }
    if (length(sig_pre) > 0) {
      cat(sprintf("\n  Individually significant pre-period terms (p < %.2f):\n", alpha))
      cat(paste("    ", sig_pre, collapse = "\n"), "\n")
    }
  }

  list(
    f_stat = f_stat,
    df1 = df1,
    df2 = df2,
    p_value = p_value,
    reject_null = p_value < alpha,
    n_pre_coefs = length(pre_coefs),
    significant_pre_coefs = as.character(sig_pre),
    pre_coef_details = pre_details
  )
}


#' Compile Pre-Trends Summary Table Across Outcomes
#'
#' Applies test_pre_trends() to multiple outcomes and returns a summary table.
#'
#' @param RESULTS The nested RESULTS list from the analysis script
#' @param outcomes Character vector of base outcome names to test
#' @param model_key Which model to test (default: "SUNAB_model")
#' @param rel_key Which relative time to test (default: "age0")
#' @param covars_on Logical, test covariates-on models (default: TRUE)
#' @param alpha Significance level (default: 0.10)
#'
#' @return A data.table with one row per outcome containing:
#'   - outcome, f_stat, df1, df2, p_value, reject_null
#'   - n_pre_coefs, n_significant_pre, significant_terms, interpretation
compile_pre_trends_summary <- function(
    RESULTS,
    outcomes,
    model_key = "SUNAB_model",
    rel_key = "age0",
    covars_on = TRUE,
    alpha = 0.10
) {
  cov_key <- if (covars_on) "covars_on" else "covars_off"

  results_list <- lapply(outcomes, function(y_base) {
    cell <- tryCatch({
      RESULTS[[cov_key]][[rel_key]][[y_base]]
    }, error = function(e) NULL)

    if (is.null(cell)) {
      return(data.table::data.table(
        outcome = y_base,
        f_stat = NA_real_,
        df1 = NA_integer_,
        df2 = NA_integer_,
        p_value = NA_real_,
        reject_null = NA,
        n_pre_coefs = NA_integer_,
        n_significant_pre = NA_integer_,
        significant_terms = NA_character_,
        interpretation = "Model not available"
      ))
    }

    model <- cell[[model_key]]
    if (is.null(model) || !inherits(model, "fixest")) {
      return(data.table::data.table(
        outcome = y_base,
        f_stat = NA_real_,
        df1 = NA_integer_,
        df2 = NA_integer_,
        p_value = NA_real_,
        reject_null = NA,
        n_pre_coefs = NA_integer_,
        n_significant_pre = NA_integer_,
        significant_terms = NA_character_,
        interpretation = sprintf("%s model not available", model_key)
      ))
    }

    # Run pre-trends test
    pt <- test_pre_trends(model, reference_period = -1, alpha = alpha, verbose = FALSE)

    interp <- if (is.na(pt$p_value)) {
      "Test not available"
    } else if (pt$reject_null) {
      sprintf("CONCERN: Pre-trends rejected (p=%.3f)", pt$p_value)
    } else {
      sprintf("OK: Pre-trends not rejected (p=%.3f)", pt$p_value)
    }

    data.table::data.table(
      outcome = y_base,
      f_stat = pt$f_stat,
      df1 = pt$df1,
      df2 = pt$df2,
      p_value = pt$p_value,
      reject_null = pt$reject_null,
      n_pre_coefs = pt$n_pre_coefs,
      n_significant_pre = length(pt$significant_pre_coefs),
      significant_terms = if (length(pt$significant_pre_coefs) > 0)
        paste(pt$significant_pre_coefs, collapse = "; ") else NA_character_,
      interpretation = interp
    )
  })

  data.table::rbindlist(results_list)
}
