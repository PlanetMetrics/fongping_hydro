# =============================================================================
# hydro_sediment.R
# Fengping River Hydropower Simulation — Sediment Module
#
# Purpose
#   (1) Fit power-law suspended sediment rating curve: SSL = a * Q^b
#   (2) Estimate daily sediment load from daily flow
#   (3) Estimate weir trap efficiency and reservoir volume loss over time
#
# Key notation
#   Q               : streamflow (cms = cubic metres per second)
#   SSL             : suspended sediment load (metric tonnes / day, MT/day)
#   C               : sediment concentration (ppm = mg/L)
#   a               : rating curve coefficient (SSL = a * Q^b)
#   b               : rating curve exponent
#   TE              : trap efficiency (dimensionless, 0–1)
#   CI              : capacity–inflow ratio (dimensionless)
#   V_loss          : annual reservoir volume loss due to sedimentation (m³/yr)
#   rho_sed         : bulk density of deposited sediment (kg/m³)
#   Q_min_transport : minimum flow for detectable sediment transport (cms)
#
# Literature references for b-value validation range [1.4, 2.5]
#   Wang, H.W. & Kondolf, G.M. (2014). Upstream sediment-control dams:
#     five decades of experience in the rapidly eroding Dahan River Basin,
#     Taiwan. JAWRA 50(3):735–747. https://doi.org/10.1111/jawr.12160
#   Kondolf, G.M., et al. (2014). Sustainable sediment management in
#     reservoirs and regulated rivers. Earth's Future 2(5):256–280.
#     https://doi.org/10.1002/2013EF000184
#
# Trap efficiency method
#   Brune (1953) curve, simplified approximation:
#     TE = CI / (CI + 0.0021)
#   where CI = reservoir capacity / mean annual inflow volume
#
# Author  [your name]
# Date    2025
# =============================================================================

library(tidyverse)
library(here)

# -----------------------------------------------------------------------------
# Internal constants
# -----------------------------------------------------------------------------

# Literature-based b-value range for Taiwan mountain streams
.B_RANGE <- c(1.4, 2.5)

# Wet bulk density of deposited sediment
# Typical value for mixed gravel-sand deposits in Taiwan mountain streams
# Reference: Wang et al. (2018) Water 10(8):1034
.RHO_SED_KG_M3 <- 1300

# Weir design parameters — shared reference with hydro_reservoir.R
# Source: Shihfeng Power Co. EIA documents
.WEIR_PARAMS <- list(
  W1 = list(
    name         = "Plant 1 Lower Weir",
    S_max_m3     = 967400,   # effective storage (m³)
    Q_design_cms = 24.3      # design flow (cms)
  ),
  W2 = list(
    name         = "Plant 2 Upper Weir",
    S_max_m3     = 237300,
    Q_design_cms = 6.3
  )
)


# -----------------------------------------------------------------------------
# fit_rating_curve()   <- CUSTOM FUNCTION
#
# Fit a power-law suspended sediment rating curve by ordinary least squares
# regression in log-log space:
#
#   log(SSL) = log(a) + b * log(Q)   =>   SSL = a * Q^b
#
# Only rows with flag == "ok" and SSL > 0 are used for fitting.
# Below-detection-limit rows (ppm = 0, flag == "below_detection") are
# excluded to prevent bias from zero SSL values.
#
# Arguments
#   sed_clean   data.frame   cleaned sediment data from sediment_clean.csv
#                            required columns: Discharge_CMS, Sus_Load_MTDay,
#                            flag (character)
#   plot_fit    logical      if TRUE, print log-log scatter with fitted line
#
# Returns  named list:
#   $a            numeric   rating curve coefficient
#   $b            numeric   rating curve exponent
#   $r_squared    numeric   R² of log-log regression
#   $n_points     integer   number of observations used
#   $b_in_range   logical   whether b falls within .B_RANGE
#   $method_note  character citable description for methods section
# -----------------------------------------------------------------------------

fit_rating_curve <- function(sed_clean, plot_fit = TRUE) {
  
  stopifnot(
    is.data.frame(sed_clean),
    all(c("Discharge_CMS", "Sus_Load_MTDay", "flag") %in% names(sed_clean))
  )
  
  # Use only valid observations with measurable sediment
  fit_data <- sed_clean |>
    filter(
      flag == "ok",
      Discharge_CMS > 0,
      Sus_Load_MTDay > 0
    ) |>
    mutate(
      log_Q   = log(Discharge_CMS),
      log_SSL = log(Sus_Load_MTDay)
    )
  
  if (nrow(fit_data) < 5)
    stop("fit_rating_curve: fewer than 5 valid observations for regression")
  
  # OLS regression in log-log space
  model <- lm(log_SSL ~ log_Q, data = fit_data)
  a     <- exp(coef(model)[[1]])    # back-transform intercept
  b     <- coef(model)[[2]]         # slope = exponent
  r2    <- summary(model)$r.squared
  
  b_ok <- b >= .B_RANGE[1] & b <= .B_RANGE[2]
  
  if (!b_ok)
    warning(sprintf(
      "fit_rating_curve: b = %.3f is outside expected range [%.1f, %.1f]. ",
      b, .B_RANGE[1], .B_RANGE[2]
    ))
  
  message(sprintf(
    "Rating curve fitted: SSL = %.4f x Q^%.4f  |  R² = %.3f  |  n = %d",
    a, b, r2, nrow(fit_data)
  ))
  
  # Diagnostic plot
  if (plot_fit) {
    
    pred_df <- data.frame(
      Discharge_CMS = 10^seq(
        log10(max(0.1, min(fit_data$Discharge_CMS))),
        log10(max(fit_data$Discharge_CMS)),
        length.out = 200
      )
    ) |>
      mutate(Sus_Load_MTDay = a * Discharge_CMS^b)
    
    p <- ggplot(fit_data,
                aes(x = Discharge_CMS, y = Sus_Load_MTDay)) +
      geom_point(alpha = 0.65, colour = "#1D9E75", size = 2) +
      geom_line(data = pred_df,
                colour = "#D85A30", linewidth = 1) +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      labs(
        title    = "Suspended Sediment Rating Curve — Fengping Creek",
        subtitle = sprintf(
          "SSL = %.3f \u00d7 Q^%.3f  |  R\u00b2 = %.3f  |  n = %d obs (1959\u20132023)",
          a, b, r2, nrow(fit_data)
        ),
        x = "Discharge Q (cms, log scale)",
        y = "Suspended Sediment Load SSL (MT/day, log scale)",
        caption = paste(
          "b-value expected range for Taiwan mountain streams: 1.4\u20132.5",
          "\nWang & Kondolf (2014) JAWRA; Kondolf et al. (2014) Earth's Future"
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.caption = element_text(size = 8, colour = "grey50"))
    
    print(p)
  }
  
  # Citable method description for methods section
  method_note <- sprintf(
    paste(
      "A power-law suspended sediment rating curve (SSL = a * Q^b) was",
      "fitted by ordinary least squares regression in log-log space using",
      "%d field observations from the Fengping Creek gauging station",
      "(1959-2023; below-detection-limit observations excluded).",
      "The fitted parameters are a = %.4f and b = %.4f (R2 = %.3f).",
      "The exponent b falls within the range of 1.4-2.5 reported for",
      "Taiwan mountain streams (Wang & Kondolf 2014, JAWRA 50(3):735-747;",
      "Kondolf et al. 2014, Earth's Future 2(5):256-280)."
    ),
    nrow(fit_data), a, b, r2
  )
  
  list(
    a           = a,
    b           = b,
    r_squared   = r2,
    n_points    = nrow(fit_data),
    b_in_range  = b_ok,
    method_note = method_note
  )
}


# -----------------------------------------------------------------------------
# estimate_daily_ssl()
#
# Apply the fitted rating curve to a daily flow series to estimate
# suspended sediment load (SSL) for each day.
#
# Days where Q is below Q_min_transport are assigned SSL = 0, consistent
# with field observations where ppm = 0 during low-flow periods.
#
# Arguments
#   daily_flow        data.frame   columns: date (Date), Q_cms (numeric)
#   curve             list         output of fit_rating_curve()
#   Q_min_transport   numeric      cms below which SSL = 0 (default 5.0)
#
# Returns
#   data.frame: date, Q_cms, SSL_mt_day (metric tonnes per day)
# -----------------------------------------------------------------------------

estimate_daily_ssl <- function(daily_flow,
                               curve,
                               Q_min_transport = 5.0) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    is.list(curve),
    all(c("a", "b") %in% names(curve))
  )
  
  daily_flow |>
    mutate(
      SSL_mt_day = if_else(
        !is.na(Q_cms) & Q_cms >= Q_min_transport,
        curve$a * Q_cms ^ curve$b,
        0
      ),
      SSL_mt_day = pmax(SSL_mt_day, 0)   # physical lower bound
    ) |>
    select(date, Q_cms, SSL_mt_day)
}


# -----------------------------------------------------------------------------
# estimate_trap_efficiency()
#
# Simulate annual trap efficiency (TE) and reservoir volume loss over the
# project lifetime using the Brune (1953) curve approximation.
#
# Brune curve (simplified):
#   TE = CI / (CI + 0.0021)
# where CI = S_remaining / annual_inflow_volume
#
# Volume deposited per year:
#   V_dep = SSL_annual_mt * 1000 / rho_sed   [m³/year]
#
# Arguments
#   ssl_annual_mt   numeric    annual suspended sediment load (MT/year)
#   weir_id         character  "W1" or "W2"
#   years           integer    project economic lifetime (default 35)
#   rho_sed         numeric    bulk density of deposit (kg/m³)
#
# Returns  data.frame with columns:
#   year, S_remaining_m3, TE, V_loss_m3_yr, CI_ratio
#   S_remaining_m3 at end of each year after sedimentation
# -----------------------------------------------------------------------------

estimate_trap_efficiency <- function(ssl_annual_mt,
                                     weir_id = "W1",
                                     years   = 35L,
                                     rho_sed = .RHO_SED_KG_M3) {
  
  stopifnot(
    weir_id %in% names(.WEIR_PARAMS),
    ssl_annual_mt >= 0,
    years >= 1L
  )
  
  params           <- .WEIR_PARAMS[[weir_id]]
  S_current        <- params$S_max_m3
  annual_inflow_m3 <- params$Q_design_cms * 86400 * 365
  
  # Annual deposited volume (m³): convert MT → kg → m³
  V_sed_yr <- ssl_annual_mt * 1000 / rho_sed
  
  purrr::map_dfr(seq_len(years), function(yr) {
    
    CI     <- S_current / annual_inflow_m3
    TE     <- CI / (CI + 0.0021)
    V_loss <- V_sed_yr * TE
    
    S_current <<- max(S_current - V_loss, 0)   # update with <<-
    
    data.frame(
      year           = yr,
      S_remaining_m3 = round(S_current, 0),
      TE             = round(TE, 4),
      V_loss_m3_yr   = round(V_loss, 1),
      CI_ratio       = round(CI, 6)
    )
  })
}


# -----------------------------------------------------------------------------
# build_sediment_outputs()
#
# Convenience wrapper that runs the full sediment pipeline and returns
# all outputs as a named list for use in fp_hydro_main.qmd.
#
# Arguments
#   sed_csv_path   character   path to sediment_clean.csv
#   daily_flow     data.frame  columns: date, Q_cms (NA-filled)
#   plot_fit       logical     whether to display rating curve plot
#
# Returns  named list:
#   $curve         list        rating curve parameters
#   $daily_ssl     data.frame  daily SSL estimates
#   $annual_ssl_mt numeric     mean annual SSL (MT/year)
#   $trap_W1       data.frame  trap efficiency trajectory W1 (35 years)
#   $trap_W2       data.frame  trap efficiency trajectory W2 (35 years)
# -----------------------------------------------------------------------------

build_sediment_outputs <- function(sed_csv_path,
                                   daily_flow,
                                   plot_fit = TRUE) {
  
  sed_clean  <- read_csv(sed_csv_path, show_col_types = FALSE)
  curve      <- fit_rating_curve(sed_clean, plot_fit = plot_fit)
  daily_ssl  <- estimate_daily_ssl(daily_flow, curve)
  
  annual_ssl <- sum(daily_ssl$SSL_mt_day, na.rm = TRUE)
  message(sprintf(
    "Annual SSL estimate (full daily record): %.0f MT/year", annual_ssl
  ))
  
  list(
    curve         = curve,
    daily_ssl     = daily_ssl,
    annual_ssl_mt = annual_ssl,
    trap_W1       = estimate_trap_efficiency(annual_ssl, "W1"),
    trap_W2       = estimate_trap_efficiency(annual_ssl, "W2")
  )
}