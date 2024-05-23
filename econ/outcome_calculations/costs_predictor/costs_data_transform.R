## LOADING PACKAGES
#setwd("~/Desktop/research asst/Global Code")
library(readr)
library(dplyr)
library(tidyverse)
library(data.table)
library(WDI)
library(ggplot2)

age_colors <- c('all' = '#d91818', 'adults' = '#e2790e', 'elderly' = '#eacb2c', 'children' = '#62adc1')

costs <- read_csv('econ/costs_predictor/national_costs.csv')[,1:14] %>% filter(!is.na(country))
costs_dt <- data.table(costs)
costs_dt <- costs_dt[, c('country','iso3c','currency','currency_iso3c','study_year','currency_year','study_pop',
                         'outcome','cost','study')][study_pop%in%c('adults','elderly','children'),]

## putting all costs into USD$2022

## back into local currency
# WDIsearch('exchange')
if(!exists('lcu_rates')){lcu_rates <- data.table(WDI(indicator='PA.NUS.FCRF', start=min(costs_dt$currency_year), end=2022))}
lcu_rates[, currency_year := year]
costs_dt <- costs_dt[lcu_rates, on = c('iso3c', 'currency_year'), lcu_rate := i.PA.NUS.FCRF]
costs_dt <- costs_dt[lcu_rates[currency_year==2022], on = c('iso3c'), lcu_rate_22 := i.PA.NUS.FCRF]
for(i in 1:nrow(costs_dt)){
  if(costs_dt$iso3c[i] %in% c('AUT','BEL','DEU','ESP','FIN','ITA','NLD')){
    costs_dt$lcu_rate[i] <- lcu_rates[iso3c=='EMU' & currency_year == costs_dt$currency_year[i]]$PA.NUS.FCRF
    costs_dt$lcu_rate_22[i] <- lcu_rates[iso3c=='EMU' & currency_year == 2022]$PA.NUS.FCRF
  }
  if(costs_dt$iso3c[i] %in% c('TWN') & costs_dt$currency_year[i] == 2010){
    costs_dt$lcu_rate[i] <- 31.5
    costs_dt$lcu_rate_22[i] <- 29.8
  }
}

# ggplot(costs_dt) + 
#   geom_point(aes(x=lcu_rate, y=lcu_rate_22, col=iso3c)) +
#   geom_line(aes(x=lcu_rate, y=lcu_rate), lty=2) +
#   theme_bw() + scale_y_log10() + scale_x_log10()

## inflate to $2022
if(!exists('deflate_rates')){deflate_rates <- data.table(WDI(indicator='NY.GDP.DEFL.ZS', start=min(costs$currency_year), end=2022))}
deflate_rates[, currency_year := year][, currency_iso3c := iso3c]
costs_dt <- costs_dt[deflate_rates, on = c('currency_iso3c', 'currency_year'), curr_yr_defl_rate := i.NY.GDP.DEFL.ZS]
costs_dt <- costs_dt[deflate_rates[year==2022,], on = c('currency_iso3c'), defl_rate_22 := i.NY.GDP.DEFL.ZS]
costs_dt[, inflator := defl_rate_22/curr_yr_defl_rate]
costs_dt[currency_iso3c == 'USA', cost_lcu := cost*lcu_rate]
costs_dt[currency_iso3c == iso3c, cost_lcu := cost]
costs_dt[, cost_lcu_2022 := cost_lcu*inflator][, cost_usd_2022 := cost_lcu_2022/lcu_rate]

## adding GDP per capita
# using 2022 GDP per capita as assuming that changes in costs are 
# absorbed into inflation...? check
if(!exists('gdp_data')){gdp_data <- WDI(indicator='NY.GDP.PCAP.KD', start=2022, end=2022)}
costs_gdp <- merge(costs_dt, gdp_data, by = c('iso3c'))
setnames(costs_gdp, 'NY.GDP.PCAP.KD', 'gdpcap')
costs_gdp[,c('country.y','iso2c','year', 'curr_yr_defl_rate','defl_rate_22') := NULL]

ggplot(costs_gdp) +
  geom_point(aes(x=gdpcap, y=cost_usd_2022, col=as.factor(study))) +
  theme_bw() + labs(col='Study') +
  facet_grid(outcome~study_pop, scales='free')  + 
  ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14))

## weighting the data for each country/study_pop

costs_weight <- copy(costs_gdp)

for(iso3c_i in unique(costs_gdp$iso3c)){
  for(study_pop_i in unique(costs_gdp$study_pop)){
    for(outcome_i in unique(costs_gdp$outcome)){
      dt <- costs_gdp[iso3c == iso3c_i & 
                           study_pop == study_pop_i &
                           outcome == outcome_i]
      if(nrow(dt)>0){
        weight_i <- 1/sqrt(nrow(dt))
        costs_weight[iso3c == iso3c_i &
                       study_pop == study_pop_i &
                       outcome == outcome_i,
                     weight := weight_i]
      }
    }
  }
}

## linear (with log) models
## joint fit

costs_jt_fit <- copy(costs_weight)

lmodel <- lm(log(cost_usd_2022) ~ 
               log(log(gdpcap))*relevel(factor(outcome),ref='outpatient')*relevel(factor(study_pop),ref='children'), 
             data = costs_jt_fit, weights = costs_jt_fit$weight) 

costs_jt_fit$study_pop <- factor(costs_jt_fit$study_pop, 
                                  levels = c('children','adults','elderly'))

costs_jt_fit[, pred_cost_jt := exp(unname(predict(lmodel, costs_jt_fit)))]

max_gdpcap <- 70000
# predicting_joint_costs <- data.table(
#   outcome = rep(c('outpatient','hospital'), each=max_gdpcap*4),
#   study_pop = rep(c('adults','all','children','elderly'), each=max_gdpcap),
#   gdpcap = rep(1:max_gdpcap, 8)
# )
predicting_joint_costs <- data.table(
  outcome = rep(c('outpatient','hospital'), each=max_gdpcap*3),
  study_pop = rep(c('adults','children','elderly'), each=max_gdpcap),
  gdpcap = rep(1:max_gdpcap, 6)
)
predicting_joint_costs[, pred_cost := exp(unname(predict(lmodel, predicting_joint_costs)))]
write_csv(predicting_joint_costs, file = 'econ/outcome_calculations/data/cost_predictor.csv')

ggplot(costs_jt_fit) +
  geom_point(aes(x=gdpcap, y=cost_usd_2022, col=as.factor(study)), size=2) +
  geom_line(data = predicting_joint_costs, aes(x=gdpcap, y=pred_cost)) +
  # geom_text(aes(x=gdpcap, y=costUSD2022, label=country.x), hjust = 0, nudge_x = 1000, check_overlap = T) +
  xlim(c(0,NA)) + labs(col='Study') + 
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14)) +
  facet_grid(outcome~study_pop, scales='free')

ggplot(costs_jt_fit) +
  geom_point(aes(x=gdpcap, y=cost_usd_2022, col=as.factor(iso3c)), size=2) +
  geom_line(data = predicting_joint_costs, aes(x=gdpcap, y=pred_cost)) +
  # geom_text(aes(x=gdpcap, y=costUSD2022, label=country.x), hjust = 0, nudge_x = 1000, check_overlap = T) +
  xlim(c(0,NA)) + labs(col='Study') + 
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14),
        legend.position='none') +
  scale_color_viridis(discrete=T) +
  facet_grid(outcome~study_pop, scales='free')

ggplot(costs_jt_fit) +
  geom_point(aes(x=(gdpcap), y=(cost_usd_2022), col=as.factor(study)), size=2) +
  geom_line(data = predicting_joint_costs, aes(x=(gdpcap), y=(pred_cost))) +
  # geom_text(aes(x=gdpcap, y=costUSD2022, label=country.x), hjust = 0, nudge_x = 1000, check_overlap = T) +
  labs(col='Study') + 
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14)) +
  facet_grid(outcome~study_pop, scales='free') +
  scale_x_log10(limits=c(1000,NA)) 

ggplot(predicting_joint_costs) +
  geom_line(aes(x=gdpcap, y=pred_cost, col=study_pop), lwd=1) +
  # geom_point(data = costs_jt_fit, aes(x=gdpcap, y=costUSD2022, col=study_pop)) +
  xlim(c(0,NA)) + labs(col='Study population') + 
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14)) +
  scale_colour_manual(values = age_colors) +
  facet_grid(outcome~., scales='free')










