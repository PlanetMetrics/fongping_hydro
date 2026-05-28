# =============================================================================
# hydro_plot.R
# Fengping River Hydropower Simulation — Plot Module
#
# Purpose
#   All ggplot2 chart builders for fp_hydro_main.qmd.
#   Each function accepts tidy data and returns a ggplot object.
#   Saving is handled by the calling notebook, not here.
#
# Naming convention
#   plot_*()   returns a ggplot object ready to print or save
#
# Shared theme
#   theme_fengping()   minimal theme with consistent typography
#
# Colour palette
#   .SCENARIO_COLS    named vector for baseline / SSP2-4.5 / SSP5-8.5
#   .EFLOW_COL        colour for e-flow threshold lines
#
# Dependencies
#   tidyverse, scales   loaded externally in fp_hydro_main.qmd
#
# Author  [your name]
# Date    2025
# =============================================================================

library(ggplot2)
library(scales)
library(tidyverse)

# -----------------------------------------------------------------------------
# Shared visual constants
# -----------------------------------------------------------------------------

# Scenario colour palette (colour-blind safe)
.SCENARIO_COLS <- c(
  "Baseline"   = "#888780",
  "SSP2-4.5"   = "#1D9E75",
  "SSP5-8.5"   = "#D85A30"
)

# E-flow threshold line colour
.EFLOW_COL <- "#C0392B"

# Sediment scenarios
.SED_COLS <- c(
  "Without sediment" = "#1D9E75",
  "With sediment"    = "#D85A30"
)


# -----------------------------------------------------------------------------
# theme_fengping()
#
# Shared ggplot2 theme for all plots in this project.
# -----------------------------------------------------------------------------

theme_fengping <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(size = base_size + 1,
                                     face = "plain",
                                     margin = margin(b = 4)),
      plot.subtitle   = element_text(size = base_size - 1,
                                     colour = "grey40",
                                     margin = margin(b = 8)),
      plot.caption    = element_text(size = base_size - 2,
                                     colour = "grey55",
                                     hjust = 0),
      axis.title      = element_text(size = base_size - 1),
      axis.text       = element_text(size = base_size - 2),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92"),
      legend.position   = "bottom",
      legend.title      = element_text(size = base_size - 2),
      legend.text       = element_text(size = base_size - 2),
      strip.text        = element_text(size = base_size - 1,
                                       face = "plain")
    )
}


# =============================================================================
# HYDROLOGY PLOTS
# =============================================================================

# -----------------------------------------------------------------------------
# plot_annual_flow_trend()
#
# Line chart of annual mean flow with a LOESS smoother.
# Annotates the three EIA comparison periods.
#
# Arguments
#   annual_stats   data.frame   output of compute_annual_stats()
#                               columns: year, mean_Q, max_Q, Q05
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_annual_flow_trend <- function(annual_stats) {
  
  # EIA period shading
  periods <- data.frame(
    label  = c("EIA\n1959–1995", "EIA Diff\n1969–2005", "Current\n1983–2019"),
    xmin   = c(1959, 1969, 1983),
    xmax   = c(1995, 2005, 2019),
    fill   = c("#F4A582", "#92C5DE", "#A1D99B"),
    y      = c(3, 6, 9)
  )
  
  ggplot(annual_stats, aes(x = year)) +
    # Period shading
    geom_rect(data = periods,
              aes(xmin = xmin, xmax = xmax,
                  ymin = -Inf, ymax = Inf, fill = label),
              alpha = 0.08, inherit.aes = FALSE) +
    scale_fill_manual(
      values = setNames(periods$fill, periods$label),
      name   = "Reference period"
    ) +
    # Annual mean flow
    geom_line(aes(y = mean_Q), colour = "grey50", linewidth = 0.5) +
    geom_smooth(aes(y = mean_Q), method = "loess", span = 0.4,
                colour = "#1D4E89", fill = "#AEC6E8",
                linewidth = 0.9, alpha = 0.25, se = TRUE) +
    labs(
      title    = "Annual Mean Daily Flow — Lishan Station (01T230)",
      subtitle = "LOESS smoother with 95% confidence band; shaded = EIA reference periods",
      x        = "Year",
      y        = "Annual mean flow (cms)",
      caption  = "Source: Water Resources Agency HYDROINFO; station 01T230 (1958–2025)"
    ) +
    theme_fengping()
}


# -----------------------------------------------------------------------------
# plot_rbi_trend()
#
# Scatter plot of annual Richards-Baker Flashiness Index with trend line.
#
# Arguments
#   annual_stats   data.frame   output of compute_annual_stats()
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_rbi_trend <- function(annual_stats) {
  
  ggplot(annual_stats, aes(x = year, y = rbi)) +
    geom_point(colour = "#D85A30", alpha = 0.7, size = 1.8) +
    geom_smooth(method = "lm", colour = "#1D4E89",
                fill = "#AEC6E8", alpha = 0.25,
                linewidth = 0.9, se = TRUE) +
    geom_hline(yintercept = 0.4, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    annotate("text", x = min(annual_stats$year, na.rm = TRUE) + 1,
             y = 0.42, label = "Flashy threshold (RBI = 0.4)",
             size = 3, colour = "grey40", hjust = 0) +
    labs(
      title    = "Richards-Baker Flashiness Index — Annual Trend",
      subtitle = "Higher RBI indicates more flashy (typhoon-dominated) flow regime",
      x        = "Year",
      y        = "RBI (dimensionless)",
      caption  = paste(
        "RBI = sum|q_i - q_{i-1}| / sum(q_i)",
        "\nBaker et al. (2004) JAWRA 40(2):503-522"
      )
    ) +
    theme_fengping()
}


# -----------------------------------------------------------------------------
# plot_rolling_stats()
#
# Line chart of rolling 30-year mean and CV to visualise non-stationarity.
#
# Arguments
#   rolling_stats   data.frame   output of compute_rolling_stats()
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_rolling_stats <- function(rolling_stats) {
  
  df_long <- rolling_stats |>
    select(year_centre, roll_mean_Q, roll_cv_Q) |>
    pivot_longer(-year_centre,
                 names_to  = "metric",
                 values_to = "value") |>
    mutate(
      metric_label = if_else(
        metric == "roll_mean_Q",
        "Rolling mean flow (cms)",
        "Rolling CV of flow"
      )
    )
  
  ggplot(df_long, aes(x = year_centre, y = value,
                      colour = metric_label)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
    scale_colour_manual(
      values = c("Rolling mean flow (cms)" = "#1D9E75",
                 "Rolling CV of flow"      = "#D85A30"),
      guide  = "none"
    ) +
    labs(
      title    = "Rolling 30-Year Flow Statistics",
      subtitle = "Non-stationarity assessment; each point = mean over 30-year window",
      x        = "Window centre year",
      y        = NULL
    ) +
    theme_fengping()
}


# -----------------------------------------------------------------------------
# plot_scenario_envelope()
#
# Ribbon plot of Monte Carlo P10/P50/P90 mean annual flow for three scenarios.
#
# Arguments
#   climate_out   list   output of build_climate_scenarios()
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_scenario_envelope <- function(climate_out) {
  
  # Extract annual mean Q from each simulation for each scenario
  df <- purrr::imap_dfr(climate_out, function(scen, ssp_name) {
    
    label <- .SSP_PARAMS[[ssp_name]]$label
    
    scen$summary$annual_stats |>
      mutate(scenario = label)
  })
  
  df_sum <- df |>
    group_by(scenario) |>
    summarise(
      P10 = quantile(mean_Q, 0.10),
      P25 = quantile(mean_Q, 0.25),
      P50 = quantile(mean_Q, 0.50),
      P75 = quantile(mean_Q, 0.75),
      P90 = quantile(mean_Q, 0.90),
      .groups = "drop"
    ) |>
    mutate(scenario = factor(scenario, levels = rev(
      purrr::map_chr(.SSP_PARAMS, "label")
    )))
  
  ggplot(df_sum, aes(y = scenario, colour = scenario)) +
    geom_linerange(aes(xmin = P10, xmax = P90),
                   linewidth = 2.5, alpha = 0.35) +
    geom_linerange(aes(xmin = P25, xmax = P75),
                   linewidth = 4.5, alpha = 0.45) +
    geom_point(aes(x = P50), size = 3.5) +
    scale_colour_manual(values = .SCENARIO_COLS, guide = "none") +
    labs(
      title    = "Simulated Annual Mean Flow by Climate Scenario",
      subtitle = "Inner bar = P25–P75; outer bar = P10–P90; point = P50",
      x        = "Annual mean flow (cms)",
      y        = NULL,
      caption  = paste0(
        "Monte Carlo bootstrap (n = ",
        climate_out[[1]]$n_sim,
        " simulated years); SSP parameters: IPCC AR6"
      )
    ) +
    theme_fengping()
}


# =============================================================================
# RESERVOIR / EFLOW PLOTS
# =============================================================================

# -----------------------------------------------------------------------------
# plot_eflow_satisfaction()
#
# Bar chart of e-flow satisfaction rate (%) by e-flow mode and scenario.
#
# Arguments
#   results_df   data.frame   columns: scenario, eflow_mode,
#                             satisfaction_pct (0–100)
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_eflow_satisfaction <- function(results_df) {
  
  ggplot(results_df,
         aes(x = eflow_mode, y = satisfaction_pct,
             fill = scenario)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.6, alpha = 0.85) +
    geom_hline(yintercept = 100, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = .SCENARIO_COLS, name = "Climate scenario") +
    scale_y_continuous(limits   = c(0, 105),
                       labels   = function(x) paste0(x, "%"),
                       breaks   = seq(0, 100, 20)) +
    labs(
      title    = "Ecological Flow Satisfaction Rate",
      subtitle = "Percentage of days where inflow meets e-flow requirement",
      x        = "E-flow mode",
      y        = "Satisfaction rate (%)",
      caption  = paste(
        "Committed: W1 = 0.48 cms, W2 = 0.06 cms (EIA commitment)",
        "\nRecommended: W1 = 5.0 cms, W2 = 1.0 cms (expert minimum)"
      )
    ) +
    theme_fengping()
}


# -----------------------------------------------------------------------------
# plot_mcfi_comparison()
#
# Bar chart of Multi-demand Conflict Frequency Index by scenario and e-flow.
#
# Arguments
#   mcfi_df   data.frame   columns: scenario, eflow_mode, mcfi (0–1)
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_mcfi_comparison <- function(mcfi_df) {
  
  ggplot(mcfi_df,
         aes(x = scenario, y = mcfi * 100, fill = eflow_mode)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.6, alpha = 0.85) +
    geom_text(aes(label = paste0(round(mcfi * 100, 1), "%")),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3) +
    scale_fill_brewer(palette = "Set2", name = "E-flow mode") +
    scale_y_continuous(limits = c(0, 105),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Multi-Demand Conflict Frequency Index (MCFI)",
      subtitle = "Fraction of days where total demand exceeds available flow",
      x        = "Climate scenario",
      y        = "MCFI (%)"
    ) +
    theme_fengping()
}


# =============================================================================
# SEDIMENT PLOTS
# =============================================================================

# -----------------------------------------------------------------------------
# plot_trap_efficiency()
#
# Line chart of reservoir trap efficiency and remaining storage over 35 years,
# comparing W1 and W2.
#
# Arguments
#   trap_W1   data.frame   output of estimate_trap_efficiency("W1")
#   trap_W2   data.frame   output of estimate_trap_efficiency("W2")
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_trap_efficiency <- function(trap_W1, trap_W2) {
  
  df <- bind_rows(
    trap_W1 |> mutate(weir = "Plant 1 (Lower, W1)"),
    trap_W2 |> mutate(weir = "Plant 2 (Upper, W2)")
  )
  
  p1 <- ggplot(df, aes(x = year, y = TE * 100, colour = weir)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("#1D9E75", "#D85A30"), name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(subtitle = "Trap efficiency (TE)",
         x = "Year of operation", y = "TE (%)") +
    theme_fengping() +
    theme(legend.position = "none")
  
  p2 <- ggplot(df, aes(x = year,
                       y = S_remaining_m3 / 1e6,
                       colour = weir)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("#1D9E75", "#D85A30"), name = NULL) +
    labs(subtitle = "Remaining effective storage",
         x = "Year of operation",
         y = "Storage (million m³)") +
    theme_fengping()
  
  # Stack with patchwork if available, otherwise return list
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    (p1 / p2) +
      plot_annotation(
        title   = "Weir Sedimentation Over 35-Year Project Lifetime",
        caption = paste(
          "Brune (1953) trap efficiency curve approximation.",
          "\nDeposit bulk density = 1300 kg/m³",
          "(Wang et al. 2018, Water 10(8):1034)"
        )
      )
  } else {
    list(trap_efficiency = p1, storage = p2)
  }
}


# =============================================================================
# FINANCE PLOTS
# =============================================================================

# -----------------------------------------------------------------------------
# plot_npv_feasibility_frontier()
#
# Heatmap of NPV across interest rate × LTV grid, with NPV = 0 contour
# showing the feasibility frontier. Optionally overlays with-sediment and
# without-sediment scenarios side by side.
#
# Arguments
#   sensitivity_df   data.frame   output of run_sensitivity()
#                                 columns: scenario, r_loan_pct, ltv_pct,
#                                 npv_b_ntd, viable
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_npv_feasibility_frontier <- function(sensitivity_df) {
  
  ggplot(sensitivity_df,
         aes(x = r_loan_pct, y = ltv_pct, fill = npv_b_ntd)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_contour(aes(z = npv_b_ntd), breaks = 0,
                 colour = "white", linewidth = 1.2,
                 linetype = "dashed") +
    scale_fill_gradient2(
      low      = "#D85A30",
      mid      = "white",
      high     = "#1D9E75",
      midpoint = 0,
      name     = "NPV\n(NTD billion)",
      labels   = function(x) sprintf("%.1f", x)
    ) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~scenario, ncol = 2) +
    labs(
      title    = "NPV Feasibility Frontier",
      subtitle = "Dashed white line = NPV breakeven (= 0); green = viable, red = loss",
      x        = "Loan interest rate (%)",
      y        = "Loan-to-value ratio (%)",
      caption  = paste(
        "CAPEX = NTD 9.7 billion (Yongwei investor conference 2025);",
        "FIT = NTD 2.8599/kWh (114年度, MOEA Energy Administration);",
        "Equity discount rate = 8%; Project lifetime = 35 years"
      )
    ) +
    theme_fengping() +
    theme(
      legend.position = "right",
      panel.grid      = element_blank()
    )
}


# -----------------------------------------------------------------------------
# plot_npv_sediment_comparison()
#
# Side-by-side bar chart comparing NPV with and without sediment-induced
# storage loss, across climate scenarios.
#
# Arguments
#   npv_df   data.frame   columns: scenario, sediment, npv_b_ntd
#                         sediment: "Without sediment" | "With sediment"
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_npv_sediment_comparison <- function(npv_df) {
  
  ggplot(npv_df,
         aes(x = scenario, y = npv_b_ntd, fill = sediment)) +
    geom_col(position = position_dodge(width = 0.65),
             width = 0.55, alpha = 0.85) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey20") +
    geom_text(
      aes(label = sprintf("%.1f", npv_b_ntd),
          vjust = if_else(npv_b_ntd >= 0, -0.4, 1.2)),
      position = position_dodge(width = 0.65),
      size = 3
    ) +
    scale_fill_manual(values = .SED_COLS, name = NULL) +
    labs(
      title    = "NPV With vs Without Sediment Impact",
      subtitle = "Sediment reduces effective storage and power output over 35 years",
      x        = "Climate scenario",
      y        = "NPV (NTD billion)",
      caption  = paste(
        "NPV computed at equity discount rate = 8%;",
        "LTV = 75%; Loan rate = 3%"
      )
    ) +
    theme_fengping()
}

plot_reservoir_operation <- function(sim_weir,
                                     date_start   = NULL,
                                     date_end     = NULL,
                                     weir_label   = "W1",
                                     e_flow_mode  = "recommended",
                                     show_storage = TRUE) {
  
  stopifnot(
    is.data.frame(sim_weir),
    all(c("date", "Q_in_cms", "Q_eflow_cms",
          "Q_power_cms", "Q_spill_cms", "S_end_m3") %in% names(sim_weir))
  )
  
  if (is.null(date_start)) date_start <- min(sim_weir$date)
  if (is.null(date_end))   date_end   <- min(sim_weir$date) + 89
  
  date_start <- as.Date(date_start)
  date_end   <- as.Date(date_end)
  
  df <- sim_weir |>
    filter(date >= date_start, date <= date_end)
  
  if (nrow(df) == 0)
    stop("plot_reservoir_operation: no data in selected date range")
  
  eflow_thresh <- if (e_flow_mode == "committed") {
    if (weir_label == "W1") 0.48 else 0.06
  } else if (e_flow_mode == "recommended") {
    if (weir_label == "W1") 5.0 else 1.0
  } else {
    NA_real_
  }
  
  df_flow <- df |>
    pivot_longer(
      cols      = c(Q_in_cms, Q_eflow_cms,
                    Q_power_cms, Q_spill_cms),
      names_to  = "component",
      values_to = "Q_cms"
    ) |>
    mutate(
      component = dplyr::recode(component,
                                Q_in_cms    = "Inflow",
                                Q_eflow_cms = "E-flow release",
                                Q_power_cms = "Power diversion",
                                Q_spill_cms = "Spillage"
      ),
      component = factor(component,
                         levels = c("Inflow", "Power diversion",
                                    "E-flow release", "Spillage"))
    )
  
  p_flow <- ggplot(df_flow,
                   aes(x = date, y = Q_cms, colour = component)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    scale_colour_manual(
      values = c(
        "Inflow"          = "#888780",
        "Power diversion" = "#1D4E89",
        "E-flow release"  = "#1D9E75",
        "Spillage"        = "#D85A30"
      ),
      name = NULL
    ) +
    labs(
      title    = sprintf("Reservoir Operation — %s (%s to %s)",
                         weir_label,
                         format(date_start, "%Y-%m-%d"),
                         format(date_end,   "%Y-%m-%d")),
      subtitle = sprintf("E-flow mode: %s | threshold = %.2f cms",
                         e_flow_mode,
                         if (!is.na(eflow_thresh)) eflow_thresh else 0),
      x = NULL,
      y = "Flow Q (cms)"
    ) +
    theme_fengping() +
    theme(legend.position = "bottom")
  
  if (!is.na(eflow_thresh)) {
    p_flow <- p_flow +
      geom_hline(yintercept = eflow_thresh,
                 linetype   = "dashed",
                 colour     = .EFLOW_COL,
                 linewidth  = 0.5) +
      annotate("text",
               x      = date_start + 1,
               y      = eflow_thresh * 1.08,
               label  = paste0("E-flow: ", eflow_thresh, " cms"),
               size   = 2.8,
               colour = .EFLOW_COL,
               hjust  = 0)
  }
  
  if (!show_storage) return(p_flow)
  
  s_max <- if (weir_label == "W1") 967400 else 237300
  
  p_stor <- ggplot(df, aes(x = date)) +
    geom_area(aes(y = S_end_m3 / 1e3),
              fill = "#AEC6E8", alpha = 0.6) +
    geom_line(aes(y = S_end_m3 / 1e3),
              colour = "#1D4E89", linewidth = 0.7) +
    geom_hline(yintercept = s_max / 1e3,
               linetype = "dashed",
               colour   = "grey40",
               linewidth = 0.4) +
    annotate("text",
             x      = date_start + 1,
             y      = s_max / 1e3 * 1.03,
             label  = "S_max",
             size   = 2.8,
             colour = "grey40",
             hjust  = 0) +
    labs(
      x       = NULL,
      y       = "Storage (thousand m³)",
      caption = sprintf("S_max = %s m³ | Pondage-type weir (調整池式)",
                        format(s_max, big.mark = ","))
    ) +
    theme_fengping() +
    theme(plot.caption = element_text(size = 8))
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    p_flow / p_stor + plot_layout(heights = c(3, 2))
  } else {
    message("Install patchwork: install.packages('patchwork')")
    p_flow
  }
}



# -----------------------------------------------------------------------------
# plot_energy_calendar_heatmap()
#
# Calendar heatmap of daily energy generation.
# x-axis      : month (Jan–Dec)
# y-axis      : year (most recent at top)
# fill        : daily energy (kWh), log scale optional
#
# Arguments
#   sim_weir     data.frame   output of simulate_one_weir()
#                             required columns: date, energy_kwh
#   weir_label   character    label for plot title (default "W1")
#   log_scale    logical      TRUE = log10 fill scale (default TRUE)
#   year_filter  integer or NULL
#                NULL            → all years
#                single integer  → e.g. 2010
#                vector          → e.g. c(2010, 2015, 2020)
#
# Usage
#   plot_energy_calendar_heatmap(sim_baseline$W1)
#   plot_energy_calendar_heatmap(sim_baseline$W1, year_filter = 2010)
#   plot_energy_calendar_heatmap(sim_baseline$W1, year_filter = 2010:2020)
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_energy_calendar_heatmap <- function(sim_weir,
                                         weir_label  = "W1",
                                         log_scale   = TRUE,
                                         year_filter = NULL) {
  
  stopifnot(
    is.data.frame(sim_weir),
    all(c("date", "energy_kwh") %in% names(sim_weir))
  )
  
  # Build base data
  df <- sim_weir |>
    mutate(
      year   = lubridate::year(date),
      doy    = lubridate::yday(date),
      energy = if_else(energy_kwh <= 0, 0.1, energy_kwh)
    )
  
  # Apply year filter
  if (!is.null(year_filter)) {
    df <- df |> filter(year %in% year_filter)
    if (nrow(df) == 0)
      stop(sprintf(
        "plot_energy_calendar_heatmap: no data for year(s): %s",
        paste(year_filter, collapse = ", ")
      ))
  }
  
  # Summary stats for subtitle
  annual_gwh <- df |>
    group_by(year) |>
    summarise(gwh = sum(energy_kwh, na.rm = TRUE) / 1e6,
              .groups = "drop")
  
  mean_gwh  <- mean(annual_gwh$gwh, na.rm = TRUE)
  total_gwh <- sum(annual_gwh$gwh, na.rm = TRUE)
  
  year_label <- if (is.null(year_filter)) {
    "All Years"
  } else if (length(year_filter) == 1) {
    as.character(year_filter)
  } else {
    paste0(min(year_filter), "–", max(year_filter))
  }
  
  fill_label <- if (log_scale) "Daily energy\n(kWh, log10)"
  else "Daily energy\n(kWh)"
  
  # Base plot
  p <- ggplot(
    df,
    aes(x    = doy,
        y    = factor(year,
                      levels = rev(sort(unique(year)))),
        fill = energy)
  ) +
    geom_tile(colour = NA) +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152,
                 182, 213, 244, 274, 305, 335),
      labels = c("Jan","Feb","Mar","Apr","May","Jun",
                 "Jul","Aug","Sep","Oct","Nov","Dec"),
      expand = c(0, 0)
    ) +
    labs(
      title    = sprintf(
        "Daily Energy Generation — %s (%s)", weir_label, year_label
      ),
      subtitle = sprintf(
        "Mean annual: %.1f GWh/yr | Shown total: %.0f GWh",
        mean_gwh, total_gwh
      ),
      x       = NULL,
      y       = "Year",
      fill    = fill_label,
      caption = paste(
        "Source: Lishan Station (01T230) daily flow;",
        "e-flow mode: recommended;",
        "grey = zero or below-detection generation"
      )
    ) +
    theme_fengping() +
    theme(
      axis.text.y       = element_text(size = 7),
      panel.grid        = element_blank(),
      legend.position   = "right",
      legend.key.width  = unit(0.4, "cm"),
      legend.key.height = unit(1.5, "cm")
    )
  
  # Colour scale
  if (log_scale) {
    p <- p +
      scale_fill_gradientn(
        colours   = c("#F0F4F8", "#AEC6E8",
                      "#1D9E75", "#1D4E89", "#0A1628"),
        trans     = "log10",
        na.value  = "grey90",
        labels    = scales::comma,
        guide     = guide_colourbar(
          barwidth  = 0.8,
          barheight = 8
        )
      )
  } else {
    p <- p +
      scale_fill_gradientn(
        colours  = c("#F0F4F8", "#AEC6E8",
                     "#1D9E75", "#1D4E89", "#0A1628"),
        na.value = "grey90",
        labels   = scales::comma
      )
  }
  
  p
}


# =============================================================================
# NEW: Plot selected reservoir operation window
# =============================================================================
#
# Purpose:
#   Visualise forecast-informed reservoir operation over a selected time window.
#
# Plot structure:
#   Panel 1: Inflow line + stacked operational flow areas
#            - Inflow = line
#            - E-flow, power diversion, spill = stacked area
#   Panel 2: Storage trajectory
#   Panel 3: Daily energy generation
#
# Input:
#   sim_df should be sim_forward24$W1 or sim_forward24$W2
# =============================================================================

plot_operation_window <- function(sim_df,
                                  start_date,
                                  end_date,
                                  title_text = "Reservoir operation window",
                                  show_storage = TRUE,
                                  show_energy = TRUE) {
  
  stopifnot(
    is.data.frame(sim_df),
    all(c(
      "date",
      "Q_in_cms",
      "Q_eflow_cms",
      "Q_power_cms",
      "Q_spill_cms",
      "S_end_m3",
      "energy_kwh"
    ) %in% names(sim_df))
  )
  
  window_df <- sim_df |>
    mutate(date = as.Date(date)) |>
    filter(
      date >= as.Date(start_date),
      date <= as.Date(end_date)
    )
  
  if (nrow(window_df) == 0) {
    stop("plot_operation_window: no data found in selected date range.")
  }
  
  # ---------------------------------------------------------------------------
  # Panel 1: flow operation
  # ---------------------------------------------------------------------------
  
  flow_area_df <- window_df |>
    select(
      date,
      `Ecological flow` = Q_eflow_cms,
      `Power diversion` = Q_power_cms,
      `Spill / overflow` = Q_spill_cms
    ) |>
    pivot_longer(
      cols = -date,
      names_to = "operation",
      values_to = "Q_cms"
    ) |>
    mutate(
      operation = factor(
        operation,
        levels = c(
          "Ecological flow",
          "Power diversion",
          "Spill / overflow"
        )
      )
    )
  
  p_flow <- ggplot() +
    geom_area(
      data = flow_area_df,
      aes(x = date, y = Q_cms, fill = operation),
      position = "stack",
      alpha = 0.80
    ) +
    geom_line(
      data = window_df,
      aes(x = date, y = Q_in_cms),
      linewidth = 1.0,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "Ecological flow" = "#9ACD32",   # olive bright green
        "Power diversion" = "#1F5AA6",   # royal blue
        "Spill / overflow" = "#F6B6C1"   # light peach pink
      )
    ) +
    labs(
      title = title_text,
      subtitle = "Line = inflow; stacked areas = ecological flow, power diversion, and spill",
      x = NULL,
      y = "Flow (cms)",
      fill = NULL
    ) +
    theme_minimal(base_size = 11)
  
  # ---------------------------------------------------------------------------
  # Panel 2: storage
  # ---------------------------------------------------------------------------
  
  p_storage <- ggplot(window_df, aes(x = date, y = S_end_m3 / 1e6)) +
    geom_line(linewidth = 0.9) +
    labs(
      x = NULL,
      y = expression("Storage (" * 10^6 * " m"^3 * ")")
    ) +
    theme_minimal(base_size = 11)
  
  # ---------------------------------------------------------------------------
  # Panel 3: energy generation
  # ---------------------------------------------------------------------------
  
  p_energy <- ggplot(window_df, aes(x = date, y = energy_kwh / 1000)) +
    geom_col(width = 0.8) +
    labs(
      x = NULL,
      y = "Energy (MWh/day)"
    ) +
    theme_minimal(base_size = 11)
  
  # ---------------------------------------------------------------------------
  # Combine panels
  # ---------------------------------------------------------------------------
  
  if (show_storage && show_energy) {
    p_flow / p_storage / p_energy +
      patchwork::plot_layout(heights = c(2.2, 1, 1))
    
  } else if (show_storage && !show_energy) {
    p_flow / p_storage +
      patchwork::plot_layout(heights = c(2.2, 1))
    
  } else if (!show_storage && show_energy) {
    p_flow / p_energy +
      patchwork::plot_layout(heights = c(2.2, 1))
    
  } else {
    p_flow
  }
}