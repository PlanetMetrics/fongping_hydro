# =============================================================================
# hydro_climate.R
# Fengping River Hydropower Simulation — Climate & Flow Scenario Module
#
# Purpose
#   (1) Analyse historical flow variability and non-stationarity
#   (2) Generate Monte Carlo bootstrap ensembles by annual resampling
#   (3) Apply SSP-based climate perturbations for future scenarios
#   (4) Compute rolling statistics and flashiness indices
#
# This module replaces and consolidates hydro_samplar.R.
#
# Key notation
#   Q_cms          : streamflow (cubic metres per second)
#   RBI            : Richards-Baker Flashiness Index
#                    RBI = sum|q_i - q_{i-1}| / sum(q_i)
#                    Higher RBI = more flashy (typical Taiwan torrent: 0.4–0.8)
#   n_sim          : number of Monte Carlo annual resamples
#   delta_mu       : fractional annual mean-flow shift relative to ref_year
#   delta_cv       : fractional annual CV scaling relative to ref_year
#   proj_year      : target future year for SSP perturbation
#   P[xx]          : [xx]-th percentile of simulated distribution
#   MCFI           : Multi-demand Conflict Frequency Index (from hydroclimate)
#
# Bootstrap method
#   Whole-year block bootstrap: resample complete calendar years with
#   replacement from the observed record. This preserves the within-year
#   seasonal structure and typhoon event integrity, which is critical for
#   flashy-stream sediment and power generation calculations.
#   Reference period: 1958–2025 (all available daily data, n = 67 years)
#
# SSP perturbation model (applied to bootstrapped ensemble)
#   mu_new  = mu_obs * (1 + delta_mu * (proj_year - ref_year))
#   Q_pert  = mu_new + (Q_obs - mu_obs) * cv_mult
#   cv_mult = 1 + delta_cv * (proj_year - ref_year)
#   Minimum flow: 0.01 cms (physical lower bound)
#
# SSP parameter assumptions
#   Source: IPCC AR6 WGI Chapter 11 (Taiwan eastern region precipitation)
#   baseline : no perturbation
#   ssp245   : moderate warming; delta_mu = +0.2%/yr; delta_cv = +0.3%/yr
#   ssp585   : high-end warming; delta_mu = +0.5%/yr; delta_cv = +0.8%/yr
#
# !! RESEARCH NOTE !!
#   n_sim is set to 200 for assignment runs.
#   For publication-quality results, use n_sim = 1000.
#   Search for "!!! n_sim" in this file to find the parameter.
#
# Author  [your name]
# Date    2025
# =============================================================================

library(tidyverse)
library(here)
library(lubridate)
library(zoo)

# -----------------------------------------------------------------------------
# Internal constants
# -----------------------------------------------------------------------------

.REF_YEAR <- 2019L          # SSP perturbation reference year

.SSP_PARAMS <- list(
  baseline = list(delta_mu = 0.000, delta_cv = 0.000,
                  label = "Baseline (1958–2025)"),
  ssp245   = list(delta_mu = 0.002, delta_cv = 0.003,
                  label = "SSP2-4.5 (moderate warming)"),
  ssp585   = list(delta_mu = 0.005, delta_cv = 0.008,
                  label = "SSP5-8.5 (high-end warming)")
)

# !!! n_sim = 200 for assignment; change to 1000 for research publication !!!
.N_SIM_DEFAULT <- 200L


# =============================================================================
# BLOCK 1 — Historical variability analysis
# =============================================================================

# -----------------------------------------------------------------------------
# compute_rbi()
#
# Calculate the Richards-Baker Flashiness Index for a flow time series.
#
# RBI = sum|q_i - q_{i-1}| / sum(q_i)
#
# Values near 0 indicate stable baseflow-dominated rivers.
# Values > 0.4 indicate highly flashy streams (typical Taiwan torrents).
#
# Reference
#   Baker, D.B., Richards, R.P., Loftus, T.T., & Kramer, J.W. (2004).
#   A new flashiness index: characteristics and applications to midwestern
#   rivers and streams. JAWRA 40(2):503–522.
#
# Arguments
#   Q_vec   numeric vector   daily flow (cms); must be chronologically ordered
#
# Returns
#   numeric scalar   RBI value (dimensionless)
# -----------------------------------------------------------------------------

compute_rbi <- function(Q_vec) {
  
  stopifnot(is.numeric(Q_vec), length(Q_vec) >= 2L)
  
  Q_valid <- Q_vec[!is.na(Q_vec)]
  if (length(Q_valid) < 2L) return(NA_real_)
  
  sum(abs(diff(Q_valid))) / sum(Q_valid)
}


# -----------------------------------------------------------------------------
# compute_annual_stats()
#
# Compute year-by-year summary statistics from a daily flow data.frame.
# Used to assess long-term trends and non-stationarity.
#
# Arguments
#   daily_flow   data.frame   columns: date (Date), Q_cms (numeric)
#
# Returns
#   data.frame with columns:
#     year, mean_Q, median_Q, sd_Q, cv_Q, Q95, Q05, max_Q,
#     rbi, n_days, n_valid
#   Q95 = 95th percentile (high flow); Q05 = 5th percentile (low flow)
# -----------------------------------------------------------------------------

compute_annual_stats <- function(daily_flow) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow))
  )
  
  daily_flow |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarise(
      mean_Q   = mean(Q_cms,   na.rm = TRUE),
      median_Q = median(Q_cms, na.rm = TRUE),
      sd_Q     = sd(Q_cms,     na.rm = TRUE),
      cv_Q     = sd_Q / mean_Q,
      Q05      = quantile(Q_cms, 0.05, na.rm = TRUE),
      Q95      = quantile(Q_cms, 0.95, na.rm = TRUE),
      max_Q    = max(Q_cms,    na.rm = TRUE),
      rbi      = compute_rbi(Q_cms),
      n_days   = n(),
      n_valid  = sum(!is.na(Q_cms)),
      .groups  = "drop"
    )
}


# -----------------------------------------------------------------------------
# compute_rolling_stats()
#
# Compute rolling 30-year window statistics to detect non-stationarity.
# Each window centred on a given year produces one row.
#
# Arguments
#   annual_stats   data.frame   output of compute_annual_stats()
#   window         integer      rolling window width in years (default 30)
#
# Returns
#   data.frame: year_centre, roll_mean_Q, roll_sd_Q, roll_cv_Q, roll_rbi
# -----------------------------------------------------------------------------

compute_rolling_stats <- function(annual_stats, window = 30L) {
  
  stopifnot(
    is.data.frame(annual_stats),
    "mean_Q" %in% names(annual_stats),
    nrow(annual_stats) >= window
  )
  
  half <- floor(window / 2L)
  n    <- nrow(annual_stats)
  
  purrr::map_dfr(seq(half + 1L, n - half), function(i) {
    slice_idx <- (i - half):(i + half)
    w         <- annual_stats[slice_idx, ]
    data.frame(
      year_centre   = annual_stats$year[i],
      roll_mean_Q   = mean(w$mean_Q,  na.rm = TRUE),
      roll_sd_Q     = mean(w$sd_Q,    na.rm = TRUE),
      roll_cv_Q     = mean(w$cv_Q,    na.rm = TRUE),
      roll_rbi      = mean(w$rbi,     na.rm = TRUE),
      roll_Q05      = mean(w$Q05,     na.rm = TRUE),
      roll_Q95      = mean(w$Q95,     na.rm = TRUE),
      n_years       = nrow(w)
    )
  })
}


# =============================================================================
# BLOCK 2 — Monte Carlo annual bootstrap
# =============================================================================

# -----------------------------------------------------------------------------
# build_year_index()
#
# Extract the list of complete calendar years available in the daily record.
# Years with fewer than min_days valid observations are excluded.
#
# Arguments
#   daily_flow   data.frame   columns: date, Q_cms
#   min_days     integer      minimum valid days per year (default 180)
#
# Returns
#   integer vector of usable year indices
# -----------------------------------------------------------------------------

build_year_index <- function(daily_flow, min_days = 180L) {
  
  daily_flow |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarise(n_valid = sum(!is.na(Q_cms)), .groups = "drop") |>
    filter(n_valid >= min_days) |>
    pull(year)
}


# -----------------------------------------------------------------------------
# apply_ssp_perturbation()
#
# Shift the mean and inflate the variance of a daily flow vector according
# to the chosen SSP pathway and target projection year.
#
# Perturbation model
#   years_fwd = proj_year - .REF_YEAR
#   mu_new    = mu_obs * (1 + delta_mu * years_fwd)
#   cv_mult   = 1 + delta_cv * years_fwd
#   Q_pert    = mu_new + (Q_obs - mu_obs) * cv_mult
#
# Arguments
#   Q_vec       numeric vector   daily flows for one simulated year (cms)
#   ssp         character        "baseline" | "ssp245" | "ssp585"
#   proj_year   integer          target future year
#
# Returns
#   numeric vector (same length)  perturbed daily flows (cms)
# -----------------------------------------------------------------------------

apply_ssp_perturbation <- function(Q_vec,
                                   ssp       = "baseline",
                                   proj_year = 2050L) {
  
  stopifnot(
    ssp %in% names(.SSP_PARAMS),
    proj_year >= .REF_YEAR
  )
  
  params       <- .SSP_PARAMS[[ssp]]
  years_fwd    <- as.integer(proj_year) - .REF_YEAR
  mu_obs       <- mean(Q_vec, na.rm = TRUE)
  mu_new       <- mu_obs * (1 + params$delta_mu * years_fwd)
  cv_mult      <- 1      + params$delta_cv  * years_fwd
  
  Q_pert <- mu_new + (Q_vec - mu_obs) * cv_mult
  pmax(Q_pert, 0.01)   # physical lower bound: 0.01 cms
}


# -----------------------------------------------------------------------------
# run_annual_bootstrap()   <- CUSTOM FUNCTION
#
# Generate a Monte Carlo ensemble of synthetic annual flow series by
# whole-year block bootstrap with optional SSP perturbation.
#
# Method
#   1. Draw n_sim years with replacement from usable_years.
#   2. For each drawn year, extract all daily flows for that calendar year.
#   3. Apply SSP perturbation to the daily flow vector.
#   4. Store the result as one simulated year in the ensemble.
#
# This approach preserves within-year seasonal structure and typhoon
# event integrity, which is critical for flashy-stream sediment and
# power generation calculations.
#
# Arguments
#   daily_flow    data.frame   columns: date (Date), Q_cms (numeric)
#   ssp           character    "baseline" | "ssp245" | "ssp585"
#   proj_year     integer      target projection year for SSP perturbation
#   n_sim         integer      number of synthetic years to generate
#                              !!! use 1000 for research publication !!!
#   seed          integer      RNG seed for reproducibility
#   min_days      integer      minimum valid days to include a year
#
# Returns  named list:
#   $ensemble      list of n_sim data.frames, each with columns date, Q_cms
#   $usable_years  integer vector of years used in bootstrap pool
#   $ssp           character
#   $proj_year     integer
#   $n_sim         integer
# -----------------------------------------------------------------------------

run_annual_bootstrap <- function(daily_flow,
                                 ssp       = "baseline",
                                 proj_year = 2050L,
                                 n_sim     = .N_SIM_DEFAULT,
                                 seed      = 42L,
                                 min_days  = 180L) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    ssp %in% names(.SSP_PARAMS),
    n_sim >= 10L
  )
  
  set.seed(seed)
  
  usable_years <- build_year_index(daily_flow, min_days)
  
  if (length(usable_years) < 5L)
    stop("run_annual_bootstrap: fewer than 5 usable years in daily record")
  
  message(sprintf(
    "Bootstrap pool: %d usable years (%d–%d) | ssp = %s | proj_year = %d | n_sim = %d",
    length(usable_years),
    min(usable_years), max(usable_years),
    ssp, proj_year, n_sim
  ))
  
  # Pre-split daily flow by year for fast lookup
  daily_by_year <- daily_flow |>
    mutate(year = year(date)) |>
    filter(year %in% usable_years) |>
    group_by(year) |>
    group_split()
  
  names(daily_by_year) <- purrr::map_int(daily_by_year, ~.x$year[1])
  
  # Draw years with replacement
  drawn_years <- sample(usable_years, size = n_sim, replace = TRUE)
  
  # Build ensemble
  ensemble <- purrr::map(drawn_years, function(yr) {
    
    yr_data <- daily_by_year[[as.character(yr)]]
    Q_pert  <- apply_ssp_perturbation(yr_data$Q_cms, ssp, proj_year)
    
    data.frame(
      date      = yr_data$date,
      Q_cms     = Q_pert,
      year_src  = yr        # source year label for traceability
    )
  })
  
  list(
    ensemble     = ensemble,
    usable_years = usable_years,
    ssp          = ssp,
    proj_year    = proj_year,
    n_sim        = n_sim
  )
}


# =============================================================================
# BLOCK 3 — Scenario summary statistics
# =============================================================================

# -----------------------------------------------------------------------------
# summarise_ensemble()
#
# Compute annual summary statistics across all simulated years in an ensemble.
# Outputs percentile bands for plotting uncertainty envelopes.
#
# Arguments
#   bootstrap_out   list   output of run_annual_bootstrap()
#
# Returns  named list:
#   $annual_stats   data.frame [n_sim rows]
#                   sim_id, annual_mean_Q, annual_max_Q, annual_Q05,
#                   annual_Q95, rbi, total_days
#   $percentiles    data.frame
#                   P10, P25, P50, P75, P90 for mean_Q, max_Q, Q05, Q95
# -----------------------------------------------------------------------------

summarise_ensemble <- function(bootstrap_out) {
  
  stopifnot(is.list(bootstrap_out), "ensemble" %in% names(bootstrap_out))
  
  annual_stats <- purrr::imap_dfr(bootstrap_out$ensemble, function(yr_df, i) {
    
    data.frame(
      sim_id       = i,
      year_src     = yr_df$year_src[1],
      mean_Q       = mean(yr_df$Q_cms,   na.rm = TRUE),
      median_Q     = median(yr_df$Q_cms, na.rm = TRUE),
      max_Q        = max(yr_df$Q_cms,    na.rm = TRUE),
      Q05          = quantile(yr_df$Q_cms, 0.05, na.rm = TRUE),
      Q95          = quantile(yr_df$Q_cms, 0.95, na.rm = TRUE),
      rbi          = compute_rbi(yr_df$Q_cms),
      total_days   = nrow(yr_df)
    )
  })
  
  probs  <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  pnames <- paste0("P", c(10, 25, 50, 75, 90))
  
  percentiles <- purrr::map_dfr(
    c("mean_Q", "max_Q", "Q05", "Q95", "rbi"),
    function(var) {
      vals <- annual_stats[[var]]
      q    <- quantile(vals, probs, na.rm = TRUE)
      df   <- as.data.frame(t(q))
      names(df) <- pnames
      df$variable <- var
      df
    }
  ) |>
    select(variable, everything())
  
  list(
    annual_stats = annual_stats,
    percentiles  = percentiles,
    ssp          = bootstrap_out$ssp,
    proj_year    = bootstrap_out$proj_year
  )
}


# =============================================================================
# BLOCK 4 — Convenience wrapper
# =============================================================================

# -----------------------------------------------------------------------------
# build_climate_scenarios()
#
# Run all three scenarios (baseline, SSP2-4.5, SSP5-8.5) and return
# a named list for use in fp_hydro_main.qmd and hydro_reservoir.R.
#
# Arguments
#   daily_flow   data.frame   columns: date, Q_cms (NA-filled)
#   proj_year    integer      target future year for SSP perturbation
#   n_sim        integer      bootstrap resamples per scenario
#                             !!! use 1000 for research publication !!!
#   seed         integer      RNG seed
#
# Returns  named list with slots "baseline", "ssp245", "ssp585"
#   Each slot contains output of run_annual_bootstrap() +
#   summarise_ensemble() merged together.
# -----------------------------------------------------------------------------

build_climate_scenarios <- function(daily_flow,
                                    proj_year = 2050L,
                                    n_sim     = .N_SIM_DEFAULT,
                                    seed      = 42L) {
  
  scenarios <- c("baseline", "ssp245", "ssp585")
  
  message(sprintf(
    "\nBuilding climate scenarios: proj_year = %d | n_sim = %d",
    proj_year, n_sim
  ))
  
  # !!! For research publication, n_sim should be 1000 !!!
  if (n_sim < 1000L) {
    message(sprintf(
      "  NOTE: n_sim = %d. Use n_sim = 1000 for publication-quality results.",
      n_sim
    ))
  }
  
  purrr::set_names(scenarios) |>
    purrr::map(function(ssp) {
      
      boot <- run_annual_bootstrap(
        daily_flow = daily_flow,
        ssp        = ssp,
        proj_year  = proj_year,
        n_sim      = n_sim,
        seed       = seed
      )
      
      summ <- summarise_ensemble(boot)
      
      c(boot, list(summary = summ))
    })
}
