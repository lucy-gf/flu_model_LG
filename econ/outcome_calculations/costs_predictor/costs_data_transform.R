## LOADING PACKAGES
#setwd("~/Desktop/research asst/Global Code")
library(readr)
library(dplyr)
library(tidyverse)
library(data.table)
library(WDI)
library(ggplot2)
library(viridis)

age_colors <- c('all' = '#d91818', 'adults' = '#e2790e', 'elderly' = '#eacb2c', 'children' = '#62adc1')

costs <- read_csv('econ/outcome_calculations/costs_predictor/national_costs.csv',show_col_types=F)[,1:14] %>% filter(!is.na(country))
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
gdp_data <- WDI(indicator='NY.GDP.PCAP.KD', start=2022, end=2022)
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
  # scale_color_viridis(discrete=T) +
  facet_grid(outcome~study_pop, scales='free')

ggplot(costs_jt_fit) +
  geom_point(aes(x=(gdpcap), y=(cost_usd_2022)), size=2, shape=1) +
  geom_line(data = predicting_joint_costs, aes(x=(gdpcap), y=(pred_cost))) +
  # geom_text(aes(x=gdpcap, y=costUSD2022, label=country.x), hjust = 0, nudge_x = 1000, check_overlap = T) +
  labs(col='Study') + 
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14),
        legend.position='none') +
  facet_grid(outcome~study_pop, scales='free') +
  scale_x_log10(limits=c(1000,70000),breaks=c(1000,3000,10000,30000,70000))
ggsave(("econ/outcome_calculations/costs_predictor/regression_model.png"),
       width=30,height=16,units="cm")

ggplot(predicting_joint_costs) +
  geom_line(aes(x=gdpcap, y=pred_cost, col=study_pop), lwd=1) +
  # geom_point(data = costs_jt_fit, aes(x=gdpcap, y=costUSD2022, col=study_pop)) +
  xlim(c(0,NA)) + labs(col='Study population') + 
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14)) +
  scale_colour_manual(values = age_colors) +
  facet_grid(outcome~., scales='free')

### 100 samples with uncertainty from lm

# which GDP per capita values do i need to fit for?
countries <- unique(read_csv('econ/outcome_calculations/data/national_ifrs.csv', 
                             show_col_types=F)$country_code)

gdp_data <- data.table(WDI(indicator='NY.GDP.PCAP.KD', start=2022, end=2022))
gdp_data_filt <- gdp_data[iso3c %in% unique(pred_costs$iso3c)]
setnames(gdp_data_filt, 'NY.GDP.PCAP.KD','gdpcap')
gdp_data_filt[iso3c == 'AFG', gdpcap := 355.78]
gdp_data_filt[iso3c == 'BTN', gdpcap := 3560.20]
gdp_data_filt[iso3c == 'ERI', gdpcap := 643.82]
gdp_data_filt[iso3c == 'PRK', gdpcap := 590]
gdp_data_filt[iso3c == 'LBN', gdpcap := 4136.10]
gdp_data_filt[iso3c == 'NCL', gdpcap := 35745.20]
gdp_data_filt[iso3c == 'SSD', gdpcap := 550.86]
gdp_data_filt[iso3c == 'SYR', gdpcap := 420.62]
gdp_data_filt[iso3c == 'LBN', gdpcap := 4136.10]
gdp_data_filt[iso3c == 'TON', gdpcap := 4426]
gdp_data_filt[iso3c == 'VEN', gdpcap := 15975.73]
gdp_data_filt[, c('iso2c','year') := NULL]
gdp_data_filt <- rbind(gdp_data_filt,
                       data.table(
                         country = c('French Guiana', 'Taiwan'),
                         iso3c = c('GUF','TWN'),
                         gdpcap = c(15600, 32679)
                       ))

pred_costs <- data.table(
  iso3c = rep(countries, each = 100*6),
  outcome = rep(c('outpatient','hospital'), each = 100*3),
  study_pop = rep(c('adults','children','elderly'), each=100),
  simulation_index = 1:100
)
pred_costs <- pred_costs[gdp_data_filt[,c('iso3c','gdpcap')], on='iso3c']
pred_costs <- cbind(pred_costs, (unname(predict(lmodel, pred_costs, interval = 'confidence'))))
setnames(pred_costs, 'V1','median')
setnames(pred_costs, 'V2','l95')
setnames(pred_costs, 'V3','u95')

# then need to fit beta/normal and sample from it
# then save as df and load into global_econ.R

f.normal <- function(mean, sd, x) {
  p <- pnorm(x, mean, sd)
  # return both
  return(c(p))
}
delta <- function(fit, actual) sum((fit-actual)^2)
objective <- function(theta, x, prob, ...) {
  ab <- (theta) 
  fit <- f.normal(ab[1], ab[2], x=as.numeric(x),...)
  # fit <- f.beta(ab[1], ab[2], x=as.numeric(x),...)
  return (delta(fit, prob))
}
fcn_fitting <- function(rates,
                        probs){
  x <- c(unlist(unname(rates)))
  sol <- suppressWarnings(optim(f=objective,p=c(1,1),
                                # method="BFGS",
                                x=x,
                                prob=c(probs),
                                control = list(reltol = 1e-15)
  ))
  parms <- (sol$par)       
  return(parms)
}

pb <- txtProgressBar(min = 0, max = 186*6, style = 3, width = 50, char = "=")

for(i in 1:(nrow(pred_costs)/100)){
  parms <- fcn_fitting(pred_costs[(100*i - 99),c('median','l95','u95')], c(0.5, 0.025, 0.975))
  pred_costs[((100*i - 99):(100*i)),"mean"] <- parms[1]
  pred_costs[((100*i - 99):(100*i)),"sd"] <- parms[2]
  pred_costs[((100*i - 99):(100*i)),"med_fit"] <- qnorm(p=c(0.5), mean=parms[1], sd=parms[2])
  pred_costs[((100*i - 99):(100*i)),"l95_fit"] <- qnorm(p=c(0.025), mean=parms[1], sd=parms[2])
  pred_costs[((100*i - 99):(100*i)),"u95_fit"] <- qnorm(p=c(0.975), mean=parms[1], sd=parms[2])
  setTxtProgressBar(pb, i)
}

pred_costs[,c('med_fit','l95_fit','u95_fit') := NULL]
pred_costs[, med_cost := exp(median)]

pb <- txtProgressBar(min = 0, max = 186*6*100, style = 3, width = 50, char = "=")

for(i in 1:(nrow(pred_costs))){
  sample_val <- exp(rnorm(1, mean = unlist(pred_costs[i,"mean"]), sd = unlist(pred_costs[i,"sd"])))
  pred_costs[i, sample_cost := sample_val]
  setTxtProgressBar(pb, i)
}

ggplot(pred_costs) +
  geom_point(aes(x=(gdpcap), y=(sample_cost), col=gdpcap)) +
  geom_line(data = predicting_joint_costs, aes(x=(gdpcap), y=(pred_cost))) +
  theme_bw() + ylab('Cost of care, $2022') + 
  xlab('GDP per capita, $2022') + 
  theme(text=element_text(size=14)) +
  facet_grid(outcome~study_pop, scales='free') #+
  scale_x_log10() 

save(pred_costs, file = 'econ/outcome_calculations/data/predicted_costs')




