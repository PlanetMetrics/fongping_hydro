## energy & financial gains module 

# 1) energy_production.R
# ------------------------
# Contract 
# Inputs:
#     power_flows      : numeric   daily flow through turbines (m^3/day)
#     reservoir_height : numeric   hydraulic head (m)
#     efficiency       : numeric   turbine efficiency, default 0.8
# Parameters
#' @param height reservoir height (m)
#' @param flow rate (m/s)
#' @param rho density of water = 1000 kg/m3
#' @param g acceleration due to gravity = 9.81 m/s^2
#' @param Keff 0.8 reservoir efficiency (0-1) 
#' @return Energy (E)
# Outputs:
#   numeric vector (E) daily energy generated (kWh/day)


# Physics: E = rho * g * h * V * eta   (joules)
# V = volume in m^3
# Convert J -> kWh by dividing by 3.6e6


energy_production <- function(power_flows, height, rho=1000, g=9.8,
                           Keff=0.8) {
  # here's where I calculate energy from daily flow volumes
  seconds_per_day = 86400
  j_per_kwh = 3.6e6
  
  daily_kwh = numeric(length(power_flows))
  
  for (i in seq_along(power_flows)) {
    # convert daily volume (m^3/day) to flow rate (m^3/s)
    flow = power_flows[i] / seconds_per_day
    # power in watts, same formula as power_gen_orig
    watts = rho * height * flow * g * Keff
    # energy = power * time, converted to kWh
    daily_kwh[i] = watts * seconds_per_day / j_per_kwh
  }
  return(daily_kwh)
}

# 2) money_value.R
# -------------
# Contract
#   Inputs:
#     daily_kwh     : numeric vector  daily energy generated (kWh)
#     price_per_kwh : numeric         electricity sale price ($/kWh)
#   Parameters:
#' @param days_per_year : numeric, default 365 — calendar days per year,
#                     used to scale the input period to yearly profit
#   Outputs:
#     numeric  yearly profit ($/year)

money_value = function(daily_kwh, price_per_kwh, days_per_year = 365) {
  total = 0
  for (i in seq_along(daily_kwh)) {
    total = total + daily_kwh[i] * price_per_kwh
  }
  
  # scale to yearly profit
  n_days = length(daily_kwh)
  yearly_profit = total * (days_per_year / n_days)
  
  return(yearly_profit)
}


# 3) npv_profits.R
# -------------
# Contract
#   Inputs:
#     yearly_profit : numeric  profit produced each year ($)
#     discount_rate : numeric  annual discount rate (e.g., 0.05 for 5%)
#   Parameters:
#' @param years : integer, 100 — valuation horizon (years)
#   Outputs:
#     numeric  net present value of the profit stream ($)

npv_profits <- function(yearly_profit, discount_rate, years = 100) {
  npv <- 0
  for (t in seq_len(years)) {
    # Discount year-t profit back to year 0
    npv <- npv + yearly_profit / (1 + discount_rate) ^ t
  }
  return(npv)
}