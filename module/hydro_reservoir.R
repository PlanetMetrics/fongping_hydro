# =============================================================================
# hydro_reservoir.R
# Fengping River Hydropower Simulation — Reservoir Water Balance Module
#
# Purpose
#   Simulate daily reservoir water balance for the two pondage-type weirs
#   on Fengping Creek (豐坪溪), tracking inflow, power diversion, ecological
#   flow release, spillage, and effective storage over time.
#
# Module inputs
#   daily_flow    data.frame   columns: date (Date), Q_cms (numeric)
#                              output of fill_daily_na() or raw lishan_daily.csv
#   plant_params  list         design parameters for W1 (lower) and W2 (upper)
#   e_flow_mode   character    "committed" | "recommended" | "custom"
#   e_flow_custom numeric      e-flow in cms, used only when mode = "custom"
#
# Module outputs — returned by run_reservoir_simulation()
#   data.frame with daily time series:
#     date, Q_in, Q_eflow, Q_power, Q_spill, S_start, S_end,
#     energy_kwh, eflow_satisfied, power_satisfied
#
# Key notation
#   Q_in      : daily inflow (cms)
#   Q_eflow   : ecological flow release (cms) — mandatory minimum
#   Q_avail   : flow available after e-flow (cms) = Q_in - Q_eflow
#   Q_power   : flow diverted for power generation (cms)
#   Q_spill   : overflow when storage is full (cms)
#   S         : reservoir storage (m³)
#   S_max     : maximum effective storage (m³)
#   dt        : time step = 86400 seconds (1 day)
#
# Ecological flow modes
#   "committed"   : use values committed in EIA documents
#                   W1 = 0.48 cms, W2 = 0.06 cms
#   "recommended" : expert-recommended minimum to avoid dry-out risk
#                   W1 = 5.0  cms, W2 = 1.0  cms
#   "custom"      : user-defined value applied to both weirs
#
# Note on ppm time-pricing
#   A time-of-use (TOU) pricing multiplier column is reserved in the output
#   but set to 1.0 throughout. Future extension: pass a tou_schedule
#   data.frame (date, hour, multiplier) to apply peak/off-peak pricing.
#
# Author  [your name]
# Date    2025
# =============================================================================

library(tidyverse)
library(zoo)      # for na.approx (linear interpolation)
library(here)

# -----------------------------------------------------------------------------
# Plant design parameters
# Source: Shihfeng Power Co. EIA documents and investor presentation (2025)
# -----------------------------------------------------------------------------
.PLANT_PARAMS <- list(
  
  W1 = list(                          # Lower weir — Plant 1 (Francis turbine)
    name             = "Plant 1 (Lower Weir)",
    S_max_m3         = 967400,        # effective storage (m³)
    Q_design_cms     = 24.3,          # design flow (cms)
    capacity_kw      = 18100,         # installed capacity (kW)
    efficiency       = 0.88,          # turbine-generator efficiency
    e_flow_commit    = 0.48,          # committed e-flow (cms), EIA document
    e_flow_recommend = 5.0            # expert-recommended minimum (cms)
  ),
  
  W2 = list(                          # Upper weir — Plant 2 (Pelton turbine)
    name             = "Plant 2 (Upper Weir)",
    S_max_m3         = 237300,        # effective storage (m³)
    Q_design_cms     = 6.3,           # design flow (cms)
    capacity_kw      = 19000,         # installed capacity (kW)
    efficiency       = 0.88,
    e_flow_commit    = 0.06,          # committed e-flow (cms), EIA document
    e_flow_recommend = 1.0            # expert-recommended minimum (cms)
  )
)

.SECONDS_PER_DAY <- 86400L
.J_PER_KWH      <- 3.6e6
.RHO_WATER      <- 1000              # kg/m³
.GRAVITY        <- 9.81              # m/s²

# FIT rate for small hydro > 20,000 kW (114年度, 經濟部能源署)
# Source: 中華民國114年度再生能源電能躉購費率 (2025)
# Note: TOU pricing multiplier reserved — see module header
.FIT_NTD_PER_KWH <- 2.8599


# -----------------------------------------------------------------------------
# fill_daily_na()
#
# Fill missing daily flow values using linear interpolation for short gaps
# and a fallback reference (decadal data or period mean) for long gaps.
#
# This function is kept separate from the reservoir simulation so that
# different data sources (daily JSON, decadal CSV) can be pre-processed
# with the same interface before entering the model.
#
# Arguments
#   daily_flow      data.frame   columns: date (Date), Q_cms (numeric)
#   method          character    "linear" | "decadal" | "zero"
#   max_gap_linear  integer      gaps <= this length use linear interpolation
#                                gaps > this length use fallback method
#   decadal_ref     data.frame   optional; columns: tendays (1-36), Q_cms
#                                used as fallback when method = "decadal"
#
# Returns
#   data.frame with same columns; NA values filled; a logical column
#   `was_filled` marks imputed rows for transparency
#
# Note
#   method = "zero" is provided for sensitivity testing only and is
#   physically unrealistic for Fengping Creek (flashy stream with baseflow).
# -----------------------------------------------------------------------------

fill_daily_na <- function(daily_flow,
                          method         = "linear",
                          max_gap_linear = 7L,
                          decadal_ref    = NULL) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    method %in% c("linear", "decadal", "zero")
  )
  
  df <- daily_flow |>
    arrange(date) |>
    mutate(was_filled = is.na(Q_cms))
  
  if (method == "zero") {
    df$Q_cms[is.na(df$Q_cms)] <- 0
    return(df)
  }
  
  # Identify gap lengths for each NA run
  na_runs <- rle(is.na(df$Q_cms))
  gap_lengths <- rep(na_runs$lengths, na_runs$lengths)
  short_gap   <- is.na(df$Q_cms) & gap_lengths <= max_gap_linear
  
  # Short gaps: linear interpolation
  if (any(short_gap)) {
    df$Q_cms <- zoo::na.approx(df$Q_cms, na.rm = FALSE)
  }
  
  # Long gaps: fallback strategy
  long_gap <- is.na(df$Q_cms)
  
  if (any(long_gap)) {
    if (method == "decadal" && !is.null(decadal_ref)) {
      # Map each date to its ten-day period (1-36) and use decadal mean
      df <- df |>
        mutate(
          month          = lubridate::month(date),
          day            = lubridate::day(date),
          period_in_month = case_when(
            day <= 10 ~ 1L,
            day <= 20 ~ 2L,
            TRUE      ~ 3L
          ),
          tendays = (month - 1L) * 3L + period_in_month
        ) |>
        left_join(
          decadal_ref |> rename(Q_decadal = Q_cms),
          by = "tendays"
        ) |>
        mutate(
          Q_cms = if_else(is.na(Q_cms) & !is.na(Q_decadal),
                          Q_decadal, Q_cms)
        ) |>
        select(-month, -day, -period_in_month, -tendays, -Q_decadal)
      
    } else {
      # Fallback: use 30-day rolling median (forward + backward)
      roll_med <- zoo::rollapply(df$Q_cms, width = 30, FUN = median,
                                 na.rm = TRUE, fill = NA, align = "center")
      df$Q_cms[is.na(df$Q_cms)] <- roll_med[is.na(df$Q_cms)]
    }
  }
  
  # Last resort: remaining NAs get overall median
  overall_med <- median(df$Q_cms, na.rm = TRUE)
  df$Q_cms[is.na(df$Q_cms)] <- overall_med
  
  df
}


# -----------------------------------------------------------------------------
# get_eflow()
#
# Return the ecological flow threshold (cms) for a given weir and mode.
# Called internally by run_reservoir_simulation().
#
# Arguments
#   weir_id       character   "W1" or "W2"
#   e_flow_mode   character   "committed" | "recommended" | "custom"
#   e_flow_custom numeric     cms value used only when mode = "custom"
#
# Returns
#   numeric scalar  e-flow threshold (cms)
# -----------------------------------------------------------------------------

get_eflow <- function(weir_id, e_flow_mode, e_flow_custom = NA_real_) {
  
  params <- .PLANT_PARAMS[[weir_id]]
  
  switch(e_flow_mode,
         committed    = params$e_flow_commit,
         recommended  = params$e_flow_recommend,
         custom       = {
           stopifnot(!is.na(e_flow_custom), e_flow_custom >= 0)
           e_flow_custom
         },
         stop("get_eflow: e_flow_mode must be 'committed', 'recommended', or 'custom'")
  )
}


# -----------------------------------------------------------------------------
# simulate_one_weir()   <- CUSTOM FUNCTION 1
#
# Simulate daily water balance for a single pondage weir over the full
# daily flow time series using a discrete ODE (Euler forward integration).
#
# Water balance equation (per day):
#   S(t+1) = S(t) + [Q_in - Q_eflow - Q_power - Q_spill] × dt
#
# Dispatch rules (priority order):
#   1. Ecological flow  : always released first (hard constraint)
#   2. Power generation : min(Q_avail, Q_design) routed to turbine
#   3. Spillage         : any excess when S >= S_max
#
# Arguments
#   daily_flow    data.frame   columns: date, Q_cms (already NA-filled)
#   weir_id       character    "W1" | "W2"
#   e_flow_mode   character    "committed" | "recommended" | "custom"
#   e_flow_custom numeric      cms, used only when mode = "custom"
#   S_init        numeric      initial storage (m³); default = 0.5 × S_max
#   net_head_m    numeric      effective hydraulic head (m); if NA, estimated
#                              from installed capacity and design flow
#
# Returns
#   data.frame with columns:
#     date, Q_in_cms, Q_eflow_cms, Q_power_cms, Q_spill_cms,
#     S_start_m3, S_end_m3, energy_kwh,
#     eflow_satisfied (logical), power_satisfied (logical),
#     tou_multiplier (numeric, reserved = 1.0)
# -----------------------------------------------------------------------------

simulate_one_weir <- function(daily_flow,
                              weir_id       = "W1",
                              e_flow_mode   = "recommended",
                              e_flow_custom = NA_real_,
                              S_init        = NULL,
                              net_head_m    = NA_real_) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    weir_id %in% names(.PLANT_PARAMS)
  )
  
  params   <- .PLANT_PARAMS[[weir_id]]
  S_max    <- params$S_max_m3
  Q_design <- params$Q_design_cms
  eta      <- params$efficiency
  cap_kw   <- params$capacity_kw
  
  # Estimate net head if not supplied
  if (is.na(net_head_m) || net_head_m <= 0) {
    net_head_m <- (cap_kw * 1000) /
      (eta * .RHO_WATER * .GRAVITY * Q_design)
    message(sprintf(
      "simulate_one_weir [%s]: net_head estimated as %.1f m",
      weir_id, net_head_m
    ))
  }
  
  e_flow <- get_eflow(weir_id, e_flow_mode, e_flow_custom)
  S      <- if (is.null(S_init)) 0.5 * S_max else S_init
  n      <- nrow(daily_flow)
  
  # Pre-allocate output vectors
  Q_eflow_v  <- numeric(n)
  Q_power_v  <- numeric(n)
  Q_spill_v  <- numeric(n)
  S_start_v  <- numeric(n)
  S_end_v    <- numeric(n)
  energy_v   <- numeric(n)
  eflow_ok_v <- logical(n)
  power_ok_v <- logical(n)
  
  for (i in seq_len(n)) {
    
    Q_in    <- daily_flow$Q_cms[i]
    S_start <- S
    
    # --- Priority 1: ecological flow ---
    Q_ef   <- min(e_flow, Q_in)
    eflow_ok_v[i] <- Q_in >= e_flow
    Q_avail <- max(Q_in - Q_ef, 0)
    
    # --- Priority 2: power generation (limited by design flow) ---
    Q_pw   <- min(Q_avail, Q_design)
    power_ok_v[i] <- (Q_avail >= Q_design)
    
    # Volume routed to turbine (m³) and residual
    V_power  <- Q_pw * .SECONDS_PER_DAY
    Q_resid  <- Q_avail - Q_pw
    
    # --- Storage update ---
    S_new <- S + Q_resid * .SECONDS_PER_DAY
    
    # --- Priority 3: spillage ---
    Q_sp <- 0
    if (S_new > S_max) {
      spill_vol <- S_new - S_max
      Q_sp      <- spill_vol / .SECONDS_PER_DAY
      S_new     <- S_max
    }
    
    # --- Power output (kWh) ---
    # P (kW) = η × ρ × g × Q × H / 1000
    # E (kWh) = P × hours_per_day
    P_kw      <- eta * .RHO_WATER * .GRAVITY * Q_pw * net_head_m / 1000
    P_kw      <- min(P_kw, cap_kw)       # cap at installed capacity
    E_kwh     <- P_kw * 24               # 24 hours per day
    
    # --- Store results ---
    Q_eflow_v[i] <- Q_ef
    Q_power_v[i] <- Q_pw
    Q_spill_v[i] <- Q_sp
    S_start_v[i] <- S_start
    S_end_v[i]   <- S_new
    energy_v[i]  <- E_kwh
    S            <- S_new
  }
  
  data.frame(
    date           = daily_flow$date,
    Q_in_cms       = daily_flow$Q_cms,
    Q_eflow_cms    = Q_eflow_v,
    Q_power_cms    = Q_power_v,
    Q_spill_cms    = Q_spill_v,
    S_start_m3     = S_start_v,
    S_end_m3       = S_end_v,
    energy_kwh     = energy_v,
    eflow_satisfied = eflow_ok_v,
    power_satisfied = power_ok_v,
    tou_multiplier = 1.0    # reserved for future time-of-use pricing
  )
}


# -----------------------------------------------------------------------------
# run_reservoir_simulation()
#
# Convenience wrapper: fill NA values, run both weirs, return combined results.
#
# Arguments
#   daily_flow      data.frame   raw lishan_daily.csv (date, Q_cms)
#   e_flow_mode     character    "committed" | "recommended" | "custom"
#   e_flow_custom   numeric      cms, used only when mode = "custom"
#   na_method       character    passed to fill_daily_na()
#   max_gap_linear  integer      passed to fill_daily_na()
#   decadal_ref     data.frame   optional, passed to fill_daily_na()
#
# Returns
#   named list:
#     $W1          data.frame  daily results for lower weir
#     $W2          data.frame  daily results for upper weir
#     $combined    data.frame  W1 + W2 energy summed by date
#     $e_flow_mode character
#     $filled_flow data.frame  NA-filled inflow used as input
# -----------------------------------------------------------------------------

run_reservoir_simulation <- function(daily_flow,
                                     e_flow_mode    = "recommended",
                                     e_flow_custom  = NA_real_,
                                     na_method      = "linear",
                                     max_gap_linear = 7L,
                                     decadal_ref    = NULL) {
  
  filled <- fill_daily_na(daily_flow, na_method, max_gap_linear, decadal_ref)
  
  w1 <- simulate_one_weir(filled, "W1", e_flow_mode, e_flow_custom)
  w2 <- simulate_one_weir(filled, "W2", e_flow_mode, e_flow_custom)
  
  combined <- data.frame(
    date              = w1$date,
    Q_in_cms          = w1$Q_in_cms,
    energy_kwh_W1     = w1$energy_kwh,
    energy_kwh_W2     = w2$energy_kwh,
    energy_kwh_total  = w1$energy_kwh + w2$energy_kwh,
    eflow_ok_W1       = w1$eflow_satisfied,
    eflow_ok_W2       = w2$eflow_satisfied,
    power_ok_W1       = w1$power_satisfied,
    power_ok_W2       = w2$power_satisfied
  )
  
  list(
    W1          = w1,
    W2          = w2,
    combined    = combined,
    e_flow_mode = e_flow_mode,
    filled_flow = filled
  )
}





# -----------------------------------------------------------------------------
# scale_flow_to_watershed()
#
# Scale observed streamflow from the Lishan gauge (proxy for W1 catchment,
# A = 242 km²) to the upper weir W2 catchment (A = 31 km²) using the
# area-discharge power-law relationship Q = K * A^n.
#
# Method
#   Q_W2 = Q_W1 * (A_W2 / A_W1)^n_month
#
# The exponent n varies by month, reflecting seasonal differences in
# catchment response to precipitation.
#
# Source
#   Table 4.1-5, Fengping Creek Hydropower EIA Report (1999).
#   Relationship Q = K * A^n between monthly mean daily flow and
#   catchment area across gauging stations in the Xiuguluan River basin.
#
# Arguments
#   daily_flow    data.frame   columns: date (Date), Q_cms (numeric)
#                              representing W1 catchment (Lishan proxy, 242 km²)
#   A_from_km2    numeric      source catchment area in km² (default 242, W1)
#   A_to_km2      numeric      target catchment area in km² (default 31, W2)
#
# Returns
#   data.frame with same structure: date, Q_cms (scaled to target catchment)
# -----------------------------------------------------------------------------

# Monthly n exponents from EIA Table 4.1-5
# Q = K * A^n where A = catchment area (km²), Q = monthly mean daily flow (cms)
# Correlation coefficients range from 0.919 to 0.989 across months
.EIA_N_BY_MONTH <- c(
  0.770,   # January
  0.782,   # February
  0.772,   # March
  0.765,   # April
  0.701,   # May
  0.825,   # June
  0.782,   # July
  0.980,   # August  <- highest n; typhoon season
  0.873,   # September
  0.829,   # October
  0.830,   # November
  0.795    # December
)

scale_flow_to_watershed <- function(daily_flow,
                                    A_from_km2 = 242,
                                    A_to_km2   = 31) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    A_from_km2 > 0,
    A_to_km2   > 0
  )
  
  daily_flow |>
    mutate(
      month      = lubridate::month(date),
      n_exponent = .EIA_N_BY_MONTH[month],
      Q_cms      = Q_cms * (A_to_km2 / A_from_km2) ^ n_exponent
    ) |>
    select(date, Q_cms)
}









# =============================================================================
# NEW: 24-hour forecast-informed reservoir operation
# =============================================================================
#
# This section adds a forecast-informed operation rule while preserving the
# original immediate-dispatch reservoir model.
#
# Original model:
#   simulate_one_weir()
#   run_reservoir_simulation()
#
# New model:
#   simulate_one_weir_forward24()
#   run_reservoir_simulation_forward24()
#
# The new model uses a 24-hour look-ahead rule:
#   - ecological flow first
#   - turbine flow capped by design flow
#   - storage can support low-flow generation
#   - if tomorrow is wet, the reservoir leaves storage space to reduce spill
#   - water above S_max becomes overflow/spill
# =============================================================================


.estimate_net_head <- function(weir_id) {
  
  params <- .PLANT_PARAMS[[weir_id]]
  
  (params$capacity_kw * 1000) /
    (params$efficiency * .RHO_WATER * .GRAVITY * params$Q_design_cms)
}


simulate_one_weir_forward24 <- function(daily_flow,
                                        weir_id       = "W1",
                                        e_flow_mode   = "recommended",
                                        e_flow_custom = NA_real_,
                                        S_init        = NULL,
                                        net_head_m    = NA_real_) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    weir_id %in% names(.PLANT_PARAMS)
  )
  
  df <- daily_flow |>
    arrange(date) |>
    mutate(Q_cms = pmax(Q_cms, 0))
  
  params   <- .PLANT_PARAMS[[weir_id]]
  S_max    <- params$S_max_m3
  Q_design <- params$Q_design_cms
  eta      <- params$efficiency
  cap_kw   <- params$capacity_kw
  
  if (is.na(net_head_m) || net_head_m <= 0) {
    net_head_m <- .estimate_net_head(weir_id)
  }
  
  e_flow <- get_eflow(weir_id, e_flow_mode, e_flow_custom)
  S      <- if (is.null(S_init)) 0.5 * S_max else min(max(S_init, 0), S_max)
  n      <- nrow(df)
  
  Q_eflow_v   <- numeric(n)
  Q_power_v   <- numeric(n)
  Q_spill_v   <- numeric(n)
  Q_overmax_v <- numeric(n)
  S_start_v   <- numeric(n)
  S_target_v  <- numeric(n)
  S_end_v     <- numeric(n)
  energy_v    <- numeric(n)
  eflow_ok_v  <- logical(n)
  power_ok_v  <- logical(n)
  
  V_design_day <- Q_design * .SECONDS_PER_DAY
  
  for (i in seq_len(n)) {
    
    Q_in    <- df$Q_cms[i]
    Q_next  <- if (i < n) df$Q_cms[i + 1L] else Q_in
    S_start <- S
    
    # Priority 1: ecological flow
    Q_ef <- min(e_flow, Q_in)
    eflow_ok_v[i] <- Q_in >= e_flow
    
    V_in_today <- max(Q_in - Q_ef, 0) * .SECONDS_PER_DAY
    V_available_today <- S_start + V_in_today
    
    # 24-hour forecast
    Q_next_ef <- min(e_flow, Q_next)
    V_in_next <- max(Q_next - Q_next_ef, 0) * .SECONDS_PER_DAY
    
    storage_needed_for_tomorrow <- max(0, V_design_day - V_in_next)
    empty_space_needed_tomorrow <- max(0, V_in_next - V_design_day)
    
    # Target storage:
    # keep water if tomorrow is dry;
    # leave empty space if tomorrow is wet.
    S_target <- min(storage_needed_for_tomorrow, S_max)
    S_target <- min(S_target, max(0, S_max - empty_space_needed_tomorrow))
    
    # Power release today
    V_power <- min(
      V_design_day,
      max(0, V_available_today - S_target)
    )
    
    Q_pw <- V_power / .SECONDS_PER_DAY
    power_ok_v[i] <- Q_pw >= (Q_design - 1e-6)
    
    S_new <- V_available_today - V_power
    
    # Overflow / spill
    Q_sp <- 0
    if (S_new > S_max) {
      spill_vol <- S_new - S_max
      Q_sp  <- spill_vol / .SECONDS_PER_DAY
      S_new <- S_max
    }
    
    # Diagnostic: water available above turbine design intake
    Q_overmax <- max(0, (V_available_today / .SECONDS_PER_DAY) - Q_design)
    
    # Energy generation
    P_kw  <- eta * .RHO_WATER * .GRAVITY * Q_pw * net_head_m / 1000
    P_kw  <- min(P_kw, cap_kw)
    E_kwh <- P_kw * 24
    
    Q_eflow_v[i]   <- Q_ef
    Q_power_v[i]   <- Q_pw
    Q_spill_v[i]   <- Q_sp
    Q_overmax_v[i] <- Q_overmax
    S_start_v[i]   <- S_start
    S_target_v[i]  <- S_target
    S_end_v[i]     <- S_new
    energy_v[i]    <- E_kwh
    
    S <- S_new
  }
  
  data.frame(
    date              = df$date,
    Q_in_cms          = df$Q_cms,
    Q_eflow_cms       = Q_eflow_v,
    Q_power_cms       = Q_power_v,
    Q_spill_cms       = Q_spill_v,
    Q_over_design_cms = Q_overmax_v,
    S_start_m3        = S_start_v,
    S_target_m3       = S_target_v,
    S_end_m3          = S_end_v,
    energy_kwh        = energy_v,
    eflow_satisfied   = eflow_ok_v,
    power_satisfied   = power_ok_v,
    operation_mode    = "forward24",
    tou_multiplier    = 1.0
  )
}


run_reservoir_simulation_forward24 <- function(daily_flow,
                                               e_flow_mode    = "recommended",
                                               e_flow_custom  = NA_real_,
                                               na_method      = "linear",
                                               max_gap_linear = 7L,
                                               decadal_ref    = NULL,
                                               scale_w2       = TRUE) {
  
  filled_w1 <- fill_daily_na(daily_flow, na_method, max_gap_linear, decadal_ref)
  
  filled_w2 <- if (scale_w2) {
    scale_flow_to_watershed(filled_w1 |> select(date, Q_cms))
  } else {
    filled_w1 |> select(date, Q_cms)
  }
  
  w1 <- simulate_one_weir_forward24(
    daily_flow    = filled_w1,
    weir_id       = "W1",
    e_flow_mode   = e_flow_mode,
    e_flow_custom = e_flow_custom
  )
  
  w2 <- simulate_one_weir_forward24(
    daily_flow    = filled_w2,
    weir_id       = "W2",
    e_flow_mode   = e_flow_mode,
    e_flow_custom = e_flow_custom
  )
  
  combined <- data.frame(
    date                 = w1$date,
    Q_in_cms_W1          = w1$Q_in_cms,
    Q_in_cms_W2          = w2$Q_in_cms,
    Q_power_cms_W1       = w1$Q_power_cms,
    Q_power_cms_W2       = w2$Q_power_cms,
    Q_spill_cms_W1       = w1$Q_spill_cms,
    Q_spill_cms_W2       = w2$Q_spill_cms,
    energy_kwh_W1        = w1$energy_kwh,
    energy_kwh_W2        = w2$energy_kwh,
    energy_kwh_total     = w1$energy_kwh + w2$energy_kwh,
    eflow_ok_W1          = w1$eflow_satisfied,
    eflow_ok_W2          = w2$eflow_satisfied,
    power_ok_W1          = w1$power_satisfied,
    power_ok_W2          = w2$power_satisfied
  )
  
  list(
    W1          = w1,
    W2          = w2,
    combined    = combined,
    e_flow_mode = e_flow_mode,
    filled_flow = filled_w1,
    scale_w2    = scale_w2
  )
}


summarise_annual_generation <- function(sim_out,
                                        fit_ntd_per_kwh = .FIT_NTD_PER_KWH) {
  
  sim_out$combined |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(
      energy_gwh_W1      = sum(energy_kwh_W1, na.rm = TRUE) / 1e6,
      energy_gwh_W2      = sum(energy_kwh_W2, na.rm = TRUE) / 1e6,
      energy_gwh_total   = sum(energy_kwh_total, na.rm = TRUE) / 1e6,
      power_water_m3_W1  = sum(Q_power_cms_W1, na.rm = TRUE) * .SECONDS_PER_DAY,
      power_water_m3_W2  = sum(Q_power_cms_W2, na.rm = TRUE) * .SECONDS_PER_DAY,
      spill_m3_W1        = sum(Q_spill_cms_W1, na.rm = TRUE) * .SECONDS_PER_DAY,
      spill_m3_W2        = sum(Q_spill_cms_W2, na.rm = TRUE) * .SECONDS_PER_DAY,
      revenue_ntd        = sum(energy_kwh_total, na.rm = TRUE) * fit_ntd_per_kwh,
      eflow_ok_pct_W1    = mean(eflow_ok_W1, na.rm = TRUE) * 100,
      eflow_ok_pct_W2    = mean(eflow_ok_W2, na.rm = TRUE) * 100,
      .groups = "drop"
    )
}