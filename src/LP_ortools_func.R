
VIRTUALENV_NAME = 'LP'
Sys.setenv(RETICULATE_PYTHON = paste0('~/.virtualenvs/', VIRTUALENV_NAME, '/bin/python'))


library(reticulate)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(readr)


# Sys.setenv(RETICULATE_PYTHON='~/.conda/envs/optimization/bin/python')
# Sys.setenv(RETICULATE_PYTHON='/opt/python/3.6.5/bin/python')

# cannot use conda environmet with shinyapps.io
# use virtual envs instead
# if('LP' %in% virtualenv_list()){
#   use_virtualenv('LP', required = T)
# } else {
#   virtualenv_create('LP')
#   virtualenv_install('LP', c('ortools', 'pandas', 'seaborn'))
#   use_virtualenv('LP', required = T)}



# conda_list()
# use_condaenv('optimization', required = T)

# use_condaenv('test', required = T)
# virtualenv_list()
# use_virtualenv('LP', required = T)
reticulate::py_config()

setwd('~/Documents/LP/src')

# access the python main module via the 'py' object
# py$x accesses the 'x' object created in python

# r.x would access the 'x' object created within R from Python


# flag for using outside energy to meet load obligations
use_outside_energy = F
outside_energy_cost = 10000

# use diverse profile
diverse_profile = 'Yes'

# peak load: annual peak in mw
peak_load = 1000

# restrict gas to a portion of total load, 0-1 or None
# e.g. 5 -> 5% limit on gas generation 
# and  100 -> no limit on gas generation
restrict_gas = 100

# divide by 100 since input is in percentages
# restrict_gas = restrict_gas/100


# battery parameters
min_charge_level = 0.1
init_ch_level = 0.5
batt_hours = 4
batt_eff = 0.85


# amount of CO2 per MW and MWh
# CO2 values from CSU study
gas_co2_ton_per_mwh = (411 + 854)/2000
# assumed 30 year life
gas_co2_ton_per_mw = (5981+1000+35566+8210+10165+1425)/(6*18)/30

wind_co2_ton_per_mwh = 0.2/2000
# assumed 30 year life
wind_co2_ton_per_mw = (754+10-241)/30

solar_co2_ton_per_mwh = 2.1/2000
# assumed 30 year life
solar_co2_ton_per_mw = (1202+250-46)/30

# battery C02 given in lbs
# assumed 15 year life
batt_co2_ton_per_mw = (1940400-83481+4903)/2000/15


# cost of C02 in $/ton
#co2_cost = 60
co2_cost = 1000

# calculate costs in $/MWh and $/MW 
# these costs include the cost of carbon

cc_gas_mwh = co2_cost * gas_co2_ton_per_mwh
cc_gas_mw = co2_cost * gas_co2_ton_per_mw

cc_wind_mwh = co2_cost * wind_co2_ton_per_mwh
cc_wind_mw = co2_cost * wind_co2_ton_per_mw

cc_solar_mwh = co2_cost * solar_co2_ton_per_mwh
cc_solar_mw = co2_cost * solar_co2_ton_per_mw

cc_batt_mw = co2_cost * batt_co2_ton_per_mw


# capacity cost in $/kw-mo
gas_cap_cost = 10.94 # CAPEX + FOM + Firm gas in $/kw-mo
gas_mw_cost = cc_gas_mw + gas_cap_cost*12*1000 # converted to $/MW-yr

heat_rate = 8403 # btu/kwh
vom = 6.73 # $/mwh
gas_fuel_cost = 4.37 # fuel cost in $/mmbtu
gas_trans_cost = 1.06 # transportation cost in $/mmbtu
gas_mwh_cost = cc_gas_mwh + (gas_fuel_cost+gas_trans_cost)*heat_rate/1000 + vom

batt_cost = 8.25
batt_mw_cost = cc_batt_mw + batt_cost*12*1000  # converted to $/MW-yr

wind_kw_mo = 4.85 # includes PSCo schedule 3VER 
wind_mw_cost = cc_wind_mw + wind_kw_mo*12*1000  # converted to $/MW-yr
wind_mwh_cost = cc_wind_mwh + 26.25 # $/mwh


solar_kw_mo = 0.19 # PSCo schedule 3VER
solar_mw_cost = cc_solar_mw + solar_kw_mo*12*1000  # converted to $/MW-yr
solar_mwh_cost = cc_solar_mwh + 33.68 # $/mwh


source_python('LP_ortools_func.py')
# def run_lp(diverse_profile
#            ,peak_load
#            ,restrict_gas
#            ,min_charge_level
#            ,init_ch_level
#            ,batt_hours
#            ,batt_eff
#            ,use_outside_energy
#            ,outside_energy_cost
#            ,gas_mw_cost
#            ,gas_mwh_cost
#            ,batt_mw_cost
#            ,wind_mw_cost
#            ,wind_mwh_cost
#            ,solar_mw_cost
#            ,solar_mwh_cost
# ):


# read demand and generation profiles
# note r. to get inputs from r
# and pandas is not imported as pd since this causes 
# issues with converting data
# see: https://github.com/rstudio/reticulate/issues/319
# however this is not a problem with py_to_r so
# we can use import pandas as pd

results = run_lp(diverse_profile
                 ,peak_load
                 ,restrict_gas
                 ,min_charge_level
                 ,init_ch_level
                 ,batt_hours
                 ,batt_eff
                 ,use_outside_energy
                 ,outside_energy_cost
                 ,gas_mw_cost
                 ,gas_mwh_cost
                 ,batt_mw_cost
                 ,wind_mw_cost
                 ,wind_mwh_cost
                 ,solar_mw_cost
                 ,solar_mwh_cost
)

results$obj_val

results$cap_mw

final_df = py_to_r(results$final_df) %>%
  as_tibble(rownames = 'time') 

final_df %>% glimpse()

# rename columns
final_df = final_df %>%
  mutate(time = ymd_hms(time)) %>%
  rename('load_2030' = '2030_load',
         'load_net' = 'net_load')

final_df %>% glimpse()
for(c in colnames(final_df)[c(-1,-2)]){
  final_df[,c] = round(final_df[,c], 2)
}
final_df %>% glimpse()

class(unlist(final_df[,'time'])) == 'numeric'
# function to plot hourly
plot_hourly = function(start_day = NULL, num_days = 7){
  # start_day = '2030-01-01'
  # num_days = 5
  
  df_days = date(final_df$time) %>% 
    unique() %>%
    ymd() 
  
  if(any(ymd(start_day) %in% df_days)){
    # t0 = start time, t1 = end time
    t0 = which(date(final_df$time) == ymd(start_day, tz='UTC')) %>% min()
  } else {
    
    start_day = sample(df_days[1:(length(df_days)-num_days)], 1)
    t0 = which(final_df$time == ymd(start_day))[1]
  }
  
  t1 = t0 + 24*num_days
  
  
  # source columns to plot
  # cols = c('time', 'crsp', 'lap', 'solar', 'wind', 'batt_discharge', 'gas', 'outside_energy')
  cols = c('time', 'solar', 'wind', 'batt_discharge', 'gas', 'outside_energy')
  
  col_idx = colnames(final_df) %in% cols
  
  # for stacked area chart
  pdata = final_df[t0:t1,col_idx] %>%
    pivot_longer(-time, names_to = 'Source', values_to = 'MW')
  
  # reorder levels for area plot
  pdata$Source = factor(pdata$Source, 
                        levels = c('solar', 'wind', 'batt_discharge', 'gas', 'outside_energy') %>%
                          rev())
  
  # for load lines
  df_load_2030 = final_df[t0:t1,c('time', 'load_2030', 'load_and_charge')] %>%
    pivot_longer(-time, names_to = 'Load', values_to = 'MW')
  
  
  # display.brewer.all(colorblindFriendly = TRUE)
  # The palette with black:
  cbp2 <- c("#000000", # black
            "#E69F00", # orange
            "#56B4E9", # light blue
            "#009E73", # green
            "#F0E442", # light yellow
            "#0072B2", # blue
            "#D55E00", # dark orange
            "#CC79A7") # purple
  
  
  n_source = pdata$Source %>% unique() %>% length()
  # idxs = c(2,7,8,4,5,3,6)
  idxs = c(2,7,8,4,5) # remove blues for hydro
  # cbp2 = cbp2[c(6,3,5,4,8,7,2)]
  cbp2 = cbp2[idxs[(length(idxs)-n_source+1):length(idxs)]]
  
  pdata %>%
    ggplot(aes(x = time)) +
    geom_area(aes(y = MW, fill = Source)) +
    geom_line(data = df_load_2030,
              aes(y = MW, linetype = Load), size = 0.75) +
    xlab('') +
    theme_minimal(base_size = 16) +
    scale_fill_manual(values = cbp2) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5))
  # scale_fill_viridis(discrete = TRUE)
  # scale_fill_brewer(palette = "Paired")
  # scale_fill_brewer(palette = "Set2")
  
}

plot_hourly('2030-12-30')

plot_hourly('2030-05-01')

plot_hourly('2030-01-01')

plot_hourly()

# test export capacity limit to place a value on excess energy
max_mw_export_cap = 150

excess_for_sale = ifelse(final_df$load_net > max_mw_export_cap,
                         max_mw_export_cap,
                         final_df$load_net)


summary(excess_for_sale)

10*sum(excess_for_sale)/1000000

