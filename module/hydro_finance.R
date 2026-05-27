## energy & financial gains module 

# 1) energy_production.R
# ------------------------
# Contract 
#   Inputs:
#     power_flows      : numeric   daily flow through turbines (m^3/day)
#     reservoir_height : numeric   hydraulic head (m)
#     efficiency       : numeric   turbine efficiency, default 0.8
#   Outputs:
#     numeric vector (E) daily energy generated (kWh/day)
#
# Physics: E = rho * g * h * V * eta   (joules)
#   rho = 1000 kg/m^3, g = 9.81 m/s^2, V in m^3
#   Convert J -> kWh by dividing by 3.6e6

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
#     daily_kwh     : numeric   daily energy generated (kWh)
#     price_per_kwh : numeric   electricity sale price ($/kWh)
#   Outputs:
#     numeric  total price generated over the input period

money_value = function(daily_kwh, price_per_kwh) {
  total = 0
  for (i in seq_along(daily_kwh)) {
    total = total + daily_kwh[i] * price_per_kwh
  }
  
  # yearly profit
  n_days = length(daily_kwh)
  yearly_profit = total * (365 / n_days)
  
  return(yearly_profit)
}


# 3) npv_profits.R
# -------------
# Contract 
#   Inputs:
#     yearly_profit : numeric  profit produced in year 1 ($)
#     discount_rate : numeric  annual discount rate (i.e. 0.05)
#     years         : integer  valuation horizon, years (100)
#     growth_rate   : numeric  annual profit growth (default 0 — flat)
#   Outputs:
#     numeric  net present value of the profit stream ($)

npv_profits <- function(yearly_profit,
                        discount_rate,
                        years       = 100,
                        growth_rate = 0) {
  npv <- 0
  for (t in seq_len(years)) {
    # Year-t profit grows at growth_rate, discounted back to year 0
    year_profit <- yearly_profit * (1 + growth_rate) ^ (t - 1)
    npv         <- npv + year_profit / (1 + discount_rate) ^ t
  }
  npv
}


