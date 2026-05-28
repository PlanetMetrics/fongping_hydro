# =============================================================================
# hydro_finance.R
# Fengping River Hydropower Simulation — Finance Module
#
# Purpose
#   Convert physical power output to financial performance indicators:
#   annual revenue, NPV, IRR, payback period, and interest rate sensitivity.
#   Supports comparison between with-sediment and without-sediment scenarios.
#
# Core functions (original, preserved from team member)
#   energy_production()   daily kWh from flow volume and hydraulic head
#   money_value()         annualised revenue from daily kWh and FIT rate
#   npv_profits()         NPV of a uniform annual profit stream
#
# Extended functions (added for research module)
#   build_cash_flows()    year-by-year cash flows including CAPEX
#   compute_npv()         NPV from arbitrary cash flow vector
#   compute_irr()         IRR by bisection search
#   compute_payback()     simple payback period
#   run_sensitivity()     NPV across a grid of interest rates and LTV ratios
#   build_finance_outputs() convenience wrapper
#
# Key notation
#   Q_power_cms    : flow through turbines (cms)
#   H              : effective hydraulic head (m)
#   eta            : turbine-generator efficiency (dimensionless)
#   E_kwh          : energy generated (kWh)
#   FIT            : Feed-In Tariff (NTD / kWh) — regulated purchase price
#   CAPEX          : capital expenditure (NTD) — total construction cost
#   OPEX           : annual operating and maintenance cost (NTD / year)
#   LTV            : loan-to-value ratio (fraction, e.g. 0.75 = 75% debt)
#   r_loan         : annual loan interest rate (fraction)
#   r_equity       : equity discount rate (fraction)
#   WACC           : weighted average cost of capital
#   NPV            : net present value (NTD)
#   IRR            : internal rate of return (fraction)
#   T              : project economic lifetime (years)
#
# Financial assumptions
#   FIT rate       : NTD 2.8599 / kWh  (114年度, 裝置容量 > 20,000 kW)
#                    Source: 經濟部能源署, 中華民國114年度再生能源電能躉購費率
#   CAPEX          : NTD 9.7 billion (法說會 2025/9; up from original NTD 6.4B)
#   OPEX           : NTD 97 million / year (1% of CAPEX rule of thumb)
#   Project life   : 35 years (法說會 2025/9)
#   Typical LTV    : 70–80% for utility-scale renewable energy in Taiwan
#   Loan rate range: 1%–5% (floating, current Taiwan syndicated loan range)
#   Equity rate    : 8% (assumed required return for private utility investor)
#
# Note on time-of-use (TOU) pricing
#   A TOU multiplier column is reserved in outputs but set to 1.0 throughout.
#   Future extension: pass a tou_schedule data.frame to apply peak pricing.
#
# Author  [your name + team member name]
# Date    2025
# =============================================================================

library(tidyverse)
library(here)

# -----------------------------------------------------------------------------
# Financial constants
# Update when official cost estimates or FIT rates change.
# -----------------------------------------------------------------------------

.FIN <- list(
  fit_ntd_per_kwh  = 2.8599,      # FIT rate NTD/kWh (114年度, > 20 MW class)
  capex_ntd        = 9.7e9,       # Total CAPEX NTD (法說會 2025)
  opex_ntd_yr      = 9.7e7,       # Annual O&M = 1% of CAPEX
  project_life_yr  = 35L,         # Economic lifetime (years)
  ltv_default      = 0.75,        # Default loan-to-value ratio
  r_loan_default   = 0.03,        # Default loan interest rate
  r_equity         = 0.08,        # Required equity return
  kwh_per_gwh      = 1e6          # Unit conversion
)


# =============================================================================
# ORIGINAL TEAM FUNCTIONS (preserved unchanged)
# =============================================================================

# -----------------------------------------------------------------------------
# energy_production()
#
# Calculate daily energy generation from daily flow volumes and hydraulic head.
#
# Physics: E = rho * g * H * V * eta   (joules)
#          V = volume in m³
#          Convert J → kWh by dividing by 3.6e6
#
# Arguments
#   power_flows      : numeric vector  daily flow volume (m³/day)
#   height           : numeric         hydraulic head (m)
#   rho              : numeric         water density (kg/m³), default 1000
#   g                : numeric         gravity (m/s²), default 9.8
#   Keff             : numeric         turbine efficiency (0–1), default 0.8
#
# Returns
#   numeric vector   daily energy generated (kWh/day)
# -----------------------------------------------------------------------------

energy_production <- function(power_flows, height, rho = 1000, g = 9.8,
                              Keff = 0.8) {
  
  seconds_per_day <- 86400
  j_per_kwh       <- 3.6e6
  daily_kwh       <- numeric(length(power_flows))
  
  for (i in seq_along(power_flows)) {
    flow          <- power_flows[i] / seconds_per_day   # m³/day → m³/s
    watts         <- rho * height * flow * g * Keff
    daily_kwh[i]  <- watts * seconds_per_day / j_per_kwh
  }
  return(daily_kwh)
}


# -----------------------------------------------------------------------------
# money_value()
#
# Annualise revenue from a daily energy vector and a per-kWh price.
#
# Arguments
#   daily_kwh      : numeric vector  daily energy (kWh)
#   price_per_kwh  : numeric         electricity sale price (NTD/kWh)
#   days_per_year  : numeric         scaling factor, default 365
#
# Returns
#   numeric   yearly revenue (NTD/year)
# -----------------------------------------------------------------------------

money_value <- function(daily_kwh, price_per_kwh, days_per_year = 365) {
  
  total <- 0
  for (i in seq_along(daily_kwh)) {
    total <- total + daily_kwh[i] * price_per_kwh
  }
  n_days        <- length(daily_kwh)
  yearly_profit <- total * (days_per_year / n_days)
  return(yearly_profit)
}


# -----------------------------------------------------------------------------
# npv_profits()
#
# Calculate NPV of a uniform annual profit stream over a fixed horizon.
# Year-0 CAPEX is NOT included here; use compute_npv() for full cash flows.
#
# Arguments
#   yearly_profit  : numeric   annual profit (NTD/year)
#   discount_rate  : numeric   annual discount rate (e.g. 0.05 = 5%)
#   years          : integer   valuation horizon, default 100
#
# Returns
#   numeric   net present value of the profit stream (NTD)
# -----------------------------------------------------------------------------

npv_profits <- function(yearly_profit, discount_rate, years = 100) {
  
  npv <- 0
  for (t in seq_len(years)) {
    npv <- npv + yearly_profit / (1 + discount_rate)^t
  }
  return(npv)
}


# =============================================================================
# EXTENDED FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# build_cash_flows()
#
# Construct year-by-year net cash flow vector including CAPEX at year 0.
#
# Convention
#   Year 0 : equity portion of CAPEX  (negative)
#   Year 1…T : annual_revenue − OPEX − debt_service  (positive if viable)
#
# Debt service is calculated as a level annuity (equal annual payment):
#   annuity = loan_amount * r / (1 - (1+r)^-T)
#
# Arguments
#   annual_revenue_ntd   numeric   mean annual revenue (NTD)
#   ltv                  numeric   loan-to-value ratio (0–1)
#   r_loan               numeric   annual loan interest rate (fraction)
#   capex_ntd            numeric   total CAPEX (NTD)
#   opex_ntd_yr          numeric   annual O&M cost (NTD)
#   project_life_yr      integer   project lifetime (years)
#
# Returns
#   numeric vector of length (project_life_yr + 1)
# -----------------------------------------------------------------------------

build_cash_flows <- function(annual_revenue_ntd,
                             ltv             = .FIN$ltv_default,
                             r_loan          = .FIN$r_loan_default,
                             capex_ntd       = .FIN$capex_ntd,
                             opex_ntd_yr     = .FIN$opex_ntd_yr,
                             project_life_yr = .FIN$project_life_yr) {
  
  stopifnot(
    ltv >= 0, ltv <= 1,
    r_loan >= 0,
    capex_ntd > 0,
    project_life_yr >= 1L
  )
  
  equity_ntd <- capex_ntd * (1 - ltv)
  loan_ntd   <- capex_ntd * ltv
  
  # Level debt service annuity
  debt_service <- if (r_loan == 0 || loan_ntd == 0) {
    loan_ntd / project_life_yr
  } else {
    loan_ntd * r_loan / (1 - (1 + r_loan)^(-project_life_yr))
  }
  
  net_annual <- annual_revenue_ntd - opex_ntd_yr - debt_service
  c(-equity_ntd, rep(net_annual, project_life_yr))
}


# -----------------------------------------------------------------------------
# compute_npv()
#
# Calculate NPV from an arbitrary cash flow vector using equity discount rate.
#
# Formula
#   NPV = sum_{t=0}^{T}  CF_t / (1 + r_equity)^t
#
# Arguments
#   cash_flows      numeric vector   first element = year-0 equity outflow
#   discount_rate   numeric          equity discount rate (default 8%)
#
# Returns
#   numeric scalar   NPV in NTD
# -----------------------------------------------------------------------------

compute_npv <- function(cash_flows,
                        discount_rate = .FIN$r_equity) {
  
  stopifnot(is.numeric(cash_flows), length(cash_flows) >= 2L)
  
  t  <- seq_along(cash_flows) - 1L
  pv <- cash_flows / (1 + discount_rate)^t
  sum(pv)
}


# -----------------------------------------------------------------------------
# compute_irr()
#
# Estimate IRR by bisection search (root of NPV = 0).
#
# Arguments
#   cash_flows   numeric vector
#   tol          numeric   convergence tolerance (default 1e-6)
#   max_iter     integer   maximum iterations
#
# Returns
#   numeric scalar   IRR as a fraction; NA if no root found
# -----------------------------------------------------------------------------

compute_irr <- function(cash_flows, tol = 1e-6, max_iter = 1000L) {
  
  stopifnot(is.numeric(cash_flows), length(cash_flows) >= 2L)
  
  npv_at_r <- function(r) compute_npv(cash_flows, r)
  
  lo <- -0.5; hi <- 5.0
  if (sign(npv_at_r(lo)) == sign(npv_at_r(hi))) return(NA_real_)
  
  for (i in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    if (abs(hi - lo) < tol) break
    if (sign(npv_at_r(mid)) == sign(npv_at_r(lo))) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}


# -----------------------------------------------------------------------------
# compute_payback()
#
# Simple (undiscounted) payback period: years until cumulative cash flows
# first turn positive.
#
# Arguments
#   cash_flows   numeric vector
#
# Returns
#   numeric scalar   payback years; Inf if never recovered
# -----------------------------------------------------------------------------

compute_payback <- function(cash_flows) {
  
  cumulative <- cumsum(cash_flows)
  which_pos  <- which(cumulative > 0)
  if (length(which_pos) == 0L) return(Inf)
  which_pos[1L] - 1L
}


# -----------------------------------------------------------------------------
# run_sensitivity()
#
# Compute NPV across a grid of loan interest rates and LTV ratios.
# Returns a data.frame suitable for plotting a feasibility frontier.
#
# The feasibility frontier is the NPV = 0 contour on the (r_loan, LTV) plane:
#   above the line → financially viable
#   below the line → NPV negative
#
# Arguments
#   annual_revenue_ntd   numeric      mean annual revenue (NTD)
#   r_loan_seq           numeric vec  interest rate grid (default 1%–5%)
#   ltv_seq              numeric vec  LTV grid (default 50%–90%)
#   scenario_label       character    label for output column
#
# Returns
#   data.frame: scenario, r_loan, ltv, npv_ntd, irr_pct, payback_yr, viable
# -----------------------------------------------------------------------------

run_sensitivity <- function(annual_revenue_ntd,
                            r_loan_seq    = seq(0.01, 0.05, by = 0.005),
                            ltv_seq       = seq(0.50, 0.90, by = 0.05),
                            scenario_label = "baseline") {
  
  grid <- expand.grid(r_loan = r_loan_seq, ltv = ltv_seq)
  
  purrr::pmap_dfr(grid, function(r_loan, ltv) {
    
    cf      <- build_cash_flows(annual_revenue_ntd,
                                ltv    = ltv,
                                r_loan = r_loan)
    npv_val <- compute_npv(cf)
    irr_val <- compute_irr(cf)
    pb_val  <- compute_payback(cf)
    
    data.frame(
      scenario    = scenario_label,
      r_loan_pct  = r_loan * 100,
      ltv_pct     = ltv * 100,
      npv_ntd     = npv_val,
      npv_b_ntd   = npv_val / 1e9,   # NTD billions for plotting
      irr_pct     = if (!is.na(irr_val)) irr_val * 100 else NA_real_,
      payback_yr  = pb_val,
      viable      = npv_val > 0
    )
  })
}


# -----------------------------------------------------------------------------
# build_finance_outputs()
#
# Convenience wrapper: run full finance pipeline for one scenario and return
# a named list for use in fp_hydro_main.qmd.
#
# Arguments
#   annual_gwh        numeric   annual energy (GWh)
#   scenario_label    character label for sensitivity table
#   fit_ntd_per_kwh   numeric   FIT rate override (default from .FIN)
#   ltv               numeric   LTV override
#   r_loan            numeric   loan rate override
#
# Returns  named list:
#   $annual_revenue_ntd   $cash_flows   $npv_ntd
#   $irr_pct   $payback_yr   $sensitivity
# -----------------------------------------------------------------------------

build_finance_outputs <- function(annual_gwh,
                                  scenario_label  = "baseline",
                                  fit_ntd_per_kwh = .FIN$fit_ntd_per_kwh,
                                  ltv             = .FIN$ltv_default,
                                  r_loan          = .FIN$r_loan_default) {
  
  annual_rev <- annual_gwh * .FIN$kwh_per_gwh * fit_ntd_per_kwh
  cf         <- build_cash_flows(annual_rev, ltv, r_loan)
  npv_val    <- compute_npv(cf)
  irr_val    <- compute_irr(cf)
  pb_val     <- compute_payback(cf)
  sens       <- run_sensitivity(annual_rev,
                                scenario_label = scenario_label)
  
  message(sprintf(
    "[%s] Revenue: NTD %.2fB/yr | NPV: NTD %.2fB | IRR: %.1f%% | Payback: %.0f yr",
    scenario_label,
    annual_rev / 1e9,
    npv_val / 1e9,
    if (!is.na(irr_val)) irr_val * 100 else NA,
    pb_val
  ))
  
  list(
    annual_revenue_ntd = annual_rev,
    cash_flows         = cf,
    npv_ntd            = npv_val,
    irr_pct            = if (!is.na(irr_val)) irr_val * 100 else NA_real_,
    payback_yr         = pb_val,
    sensitivity        = sens
  )
}


# =============================================================================
# NEW: 20-year NPV trajectory and financial line plot
# =============================================================================
#
# Purpose:
#   Evaluate 20-year financial performance using annual hydropower revenue.
#
# Inputs:
#   annual_generation_df:
#     output from summarise_annual_generation()
#
# Key outputs:
#   - annual revenue
#   - annual net cash flow
#   - discounted cash flow
#   - cumulative NPV trajectory
#   - NPV financial line plot
# =============================================================================


build_cash_flows_20yr <- function(annual_revenue_ntd,
                                  ltv = .FIN$ltv_default,
                                  r_loan = .FIN$r_loan_default,
                                  capex_ntd = .FIN$capex_ntd,
                                  opex_ntd_yr = .FIN$opex_ntd_yr,
                                  project_life_yr = 20L) {
  
  stopifnot(
    is.numeric(annual_revenue_ntd),
    length(annual_revenue_ntd) >= 1L,
    ltv >= 0,
    ltv <= 1,
    r_loan >= 0,
    capex_ntd > 0,
    project_life_yr == 20L
  )
  
  # Use exactly 20 years.
  annual_revenue_ntd <- rep(annual_revenue_ntd, length.out = project_life_yr)
  
  equity_ntd <- capex_ntd * (1 - ltv)
  loan_ntd   <- capex_ntd * ltv
  
  # Level debt service annuity over 20 years
  debt_service_ntd <- if (r_loan == 0 || loan_ntd == 0) {
    loan_ntd / project_life_yr
  } else {
    loan_ntd * r_loan / (1 - (1 + r_loan)^(-project_life_yr))
  }
  
  annual_net_cash_flow_ntd <- annual_revenue_ntd - opex_ntd_yr - debt_service_ntd
  
  tibble::tibble(
    year = 0:project_life_yr,
    revenue_ntd = c(0, annual_revenue_ntd),
    opex_ntd = c(0, rep(opex_ntd_yr, project_life_yr)),
    debt_service_ntd = c(0, rep(debt_service_ntd, project_life_yr)),
    cash_flow_ntd = c(-equity_ntd, annual_net_cash_flow_ntd)
  )
}


compute_npv_trajectory_20yr <- function(cash_flow_df,
                                        discount_rate = .FIN$r_equity) {
  
  stopifnot(
    is.data.frame(cash_flow_df),
    all(c("year", "cash_flow_ntd") %in% names(cash_flow_df)),
    discount_rate >= 0
  )
  
  cash_flow_df |>
    mutate(
      discount_factor = 1 / (1 + discount_rate) ^ year,
      discounted_cash_flow_ntd = cash_flow_ntd * discount_factor,
      cumulative_npv_ntd = cumsum(discounted_cash_flow_ntd),
      cumulative_npv_b_ntd = cumulative_npv_ntd / 1e9
    )
}


run_npv_20yr_from_generation <- function(annual_generation_df,
                                         ltv = .FIN$ltv_default,
                                         r_loan = .FIN$r_loan_default,
                                         discount_rate = .FIN$r_equity,
                                         capex_ntd = .FIN$capex_ntd,
                                         opex_ntd_yr = .FIN$opex_ntd_yr,
                                         project_life_yr = 20L,
                                         revenue_col = "revenue_ntd") {
  
  stopifnot(
    is.data.frame(annual_generation_df),
    revenue_col %in% names(annual_generation_df)
  )
  
  annual_revenue <- annual_generation_df[[revenue_col]]
  
  # If the historical record has more than 20 years, use the mean revenue
  # as the representative annual revenue stream.
  representative_revenue <- mean(annual_revenue, na.rm = TRUE)
  
  cash_flow_df <- build_cash_flows_20yr(
    annual_revenue_ntd = representative_revenue,
    ltv = ltv,
    r_loan = r_loan,
    capex_ntd = capex_ntd,
    opex_ntd_yr = opex_ntd_yr,
    project_life_yr = project_life_yr
  )
  
  compute_npv_trajectory_20yr(
    cash_flow_df = cash_flow_df,
    discount_rate = discount_rate
  )
}


plot_npv_financial_line_20yr <- function(npv_df,
                                         title_text = "20-year NPV financial trajectory") {
  
  stopifnot(
    is.data.frame(npv_df),
    all(c("year", "cumulative_npv_b_ntd") %in% names(npv_df))
  )
  
  ggplot(npv_df, aes(x = year, y = cumulative_npv_b_ntd)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 1.1, color = "#1F5AA6") +
    geom_point(size = 2, color = "#1F5AA6") +
    labs(
      title = title_text,
      subtitle = "Year 0 includes equity CAPEX; years 1–20 show cumulative discounted cash flow",
      x = "Project year",
      y = "Cumulative NPV (NTD billion)"
    ) +
    theme_minimal(base_size = 11)
}




plot_sensitivity_frontier <- function(frontier_df,
                                      title_text = "Feasibility frontier") {
  
  ggplot(frontier_df, aes(x = r_loan_pct, y = ltv_pct, z = npv_b_ntd)) +
    geom_tile(aes(fill = npv_b_ntd)) +
    geom_contour(
      breaks = 0,
      color = "black",
      linewidth = 0.9
    ) +
    facet_wrap(~ scenario) +
    scale_fill_viridis_c(option = "C") +
    labs(
      title = title_text,
      subtitle = "Black contour = NPV = 0",
      x = "Loan interest rate (%)",
      y = "Loan-to-value ratio (%)",
      fill = "NPV\n(NTD billion)"
    ) +
    theme_minimal(base_size = 11)
}


plot_npv_comparison_20yr <- function(npv_instant,
                                     npv_forward24,
                                     title_text = "20-year NPV comparison") {
  
  plot_df <- dplyr::bind_rows(
    npv_instant |>
      dplyr::mutate(operation = "Immediate dispatch"),
    npv_forward24 |>
      dplyr::mutate(operation = "24-hour forecast-informed")
  )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = year,
      y = cumulative_npv_b_ntd,
      color = operation
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title = title_text,
      subtitle = "Year 0 includes equity CAPEX; years 1–20 show cumulative discounted cash flow",
      x = "Project year",
      y = "Cumulative NPV (NTD billion)",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}