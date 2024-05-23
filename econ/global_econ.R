#### Economic analysis - calculating outcomes 
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_vaccine_programs.R")
source("BS/BS_colors.R")

library(readr)
library(dplyr)
library(data.table)
library(odin)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(viridis)
library(patchwork)
library(wpp2022)
library(WDI)

## SENSITIVITY ANALYSES
outp_include <- F # are we including outpatient/non-hospitalisation visits?
WTP_GDP_ratio <- 1 # what proportion of GDP per capita is the willingness_to_pay threshold?

## PARAMETERS
cost_discount_rate_val <- 0.03
DALY_discount_rate_val <- 0.03
threshold <- 30000 ## NEED TO DECIDE 

## LOADING DATA (BASE CASE)
print('loading data')

scenario_name <- 'base'
if(exists('econ_cases_agg')){rm(econ_cases_agg)};if(exists('econ_cases')){rm(econ_cases)}
for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0('data/vacc_output_', scenario_name,'/econ_inp_', c_code, '.Rdata'))){
    print(paste0(c_code, ' exists'))
    load(paste0('data/vacc_output_', scenario_name,'/econ_inp_', c_code, '.Rdata'))
    econ_cases_agg1 <- melt(econ_inp, 
         id.vars=c('country_code','simulation_index','year','scenario'))
    econ_cases_agg1[,age_grp:=as.numeric(substr(variable, 3,3))][,variable:=NULL]
    econ_cases_agg1 <- econ_cases_agg1[, lapply(.SD,sum), 
                                     by=c('country_code','simulation_index','year','scenario','age_grp')]
    setnames(econ_cases_agg1, 'value','infections')
    econ_cases_agg1[, ct := substr(scenario,4,4)]
    econ_cases_agg1[, vt := substr(scenario,9,9)]
    if(exists('econ_cases_agg')){
      econ_cases_agg <- rbind(econ_cases_agg, econ_cases_agg1)
    }else{
      econ_cases_agg <- copy(econ_cases_agg1)
    }
  }
}

rm(econ_cases_agg1);rm(econ_inp)

print('IFRs')
## load ifrs
national_ifrs <- data.table(read_csv('econ/outcome_calculations/data/national_ifrs.csv', 
                                     show_col_types=F))

econ_cases_agg <- econ_cases_agg[national_ifrs[,c('country_code','simulation_index','age_grp','ifr')], on=c('country_code','simulation_index','age_grp')]
econ_cases_agg[,deaths := ifr*infections]

print('IHRs')
## load ihrs
global_ihrs <- data.table(read_csv('econ/outcome_calculations/data/global_ihrs.csv', 
                                     show_col_types=F))

econ_cases_agg <- econ_cases_agg[global_ihrs, on=c('simulation_index','age_grp')]
econ_cases_agg[,hospitalisations := ihr*infections]

## FOR NOW ADDING PROXY RISK OF OUTPATIENT (50x HOSP)
if(outp_include == T){
  econ_cases_agg[, outpatients := 50*hospitalisations]
}

## proportion symptomatic and fever
## from carrat (DOI: 10.1093/aje/kwm375)
print('symptomatic, fever')
symp_probs <- data.table(
  outcome = c('symptoms','fever'),
  med = c(66.9, 34.9)/100,
  l95 = c(58.3, 26.7)/100,
  u95 = c(74.5, 44.2)/100
  )
## determining distribution
f.gamma <- function(shape, rate, x) {
  p <- pgamma(x, shape, rate)
  # return both
  return(c(p))
}
delta <- function(fit, actual) sum((fit-actual)^2)
objective <- function(theta, x, prob, ...) {
  ab <- (theta) 
  fit <- f.gamma(ab[1], ab[2], x=as.numeric(x),...)
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

for(i in 1:nrow(symp_probs)){
  parms <- fcn_fitting(symp_probs[i,2:4], c(0.5, 0.025, 0.975))
  symp_probs[i,"shape"] <- parms[1]
  symp_probs[i,"rate"] <- parms[2]
  symp_probs[i,"med_fit"] <- qgamma(p=c(0.5), shape=parms[1], rate=parms[2])
  symp_probs[i,"l95_fit"] <- qgamma(p=c(0.025), shape=parms[1], rate=parms[2])
  symp_probs[i,"u95_fit"] <- qgamma(p=c(0.975), shape=parms[1], rate=parms[2])
}

symp_samples <- data.table(
  simulation_index = 1:100,
  symp_prob = rgamma(100, shape = unlist(symp_probs[outcome=='symptoms','shape']), rate = unlist(symp_probs[outcome=='symptoms','rate'])),
  fever_prob = rgamma(100, shape = unlist(symp_probs[outcome=='fever','shape']), rate = unlist(symp_probs[outcome=='fever','rate'])))

econ_cases_agg <- econ_cases_agg[symp_samples, on=c('simulation_index')]
econ_cases_agg[, symptomatics := symp_prob*infections][, fevers := fever_prob*infections][, milds := symptomatics - fevers]

## YLL per death
## currently no discounting, using 'lx' method

print('YLLs')

yll_df <- national_ifrs[simulation_index==1, c('country_code','age_grp','cluster_name')]
source('econ/outcome_calculations/calc_ylls.R')
pb <- txtProgressBar(min = 0, max = 186, style = 3, width = 50, char = "=")

for(i in 1:length(unique(yll_df$country_code))){
  iso3c_i <- unique(yll_df$country_code)[i]
  yll_df[country_code == iso3c_i, 'yll'] <- yll(LT = UNLT[ISO3_code == iso3c_i & MidPeriod == 2022.5], 
                                                r = 0, # discount rate
                                                smr = 1, 
                                                weight_method = "lx", # weight method to average LE by age group: "lx" "lxqx" "equal" "pop_ifr"
                                                model_ages = c(5,20,65))$d_LEx
  setTxtProgressBar(pb, i)
}

econ_cases_agg <- econ_cases_agg[yll_df[,c('country_code','age_grp','yll')], on=c('country_code','age_grp')]
econ_cases_agg[,YLLs := yll*deaths]

## GDP per capita in 2022
print('GDP per capita')

gdp_data <- data.table(WDI(indicator='NY.GDP.PCAP.KD', start=2022, end=2022))
setnames(gdp_data, 'iso3c','country_code'); setnames(gdp_data, 'NY.GDP.PCAP.KD','GDPpc')
gdp_data <- gdp_data[country_code %in% unique(yll_df$country_code),]
econ_cases_agg <- econ_cases_agg[gdp_data[,c('country_code','GDPpc')], on=c('country_code')]

### adding missing 2022 GPD per capita
econ_cases_agg[country_code == 'AFG', GDPpc := 355.78]
econ_cases_agg[country_code == 'VEN', GDPpc := 15975.73]

econ_cases_agg[,YLL_cost := YLLs*GDPpc]

## removing probabilities
econ_cases_agg[, c('ifr','ihr','symp_prob','fever_prob','yll'):=NULL]

## YLDs
print('YLDs')

flu_duration <- 4/365 ## IS THIS CORRECT?
## weights from Dawa et al. - is this ok?
DALY_weights <- data.table(
  outcome = c('mild','fever','hospitalisation'),
  med = c(0.006, 0.051, 0.133),
  l95 = c(0.002, 0.032, 0.088),
  u95 = c(0.012, 0.074, 0.190)
)
for(i in 1:nrow(DALY_weights)){
  parms <- fcn_fitting(DALY_weights[i,2:4], c(0.5, 0.025, 0.975))
  DALY_weights[i,"shape"] <- parms[1]
  DALY_weights[i,"rate"] <- parms[2]
  DALY_weights[i,"med_fit"] <- qgamma(p=c(0.5), shape=parms[1], rate=parms[2])
  DALY_weights[i,"l95_fit"] <- qgamma(p=c(0.025), shape=parms[1], rate=parms[2])
  DALY_weights[i,"u95_fit"] <- qgamma(p=c(0.975), shape=parms[1], rate=parms[2])
}
DALY_weight_samples <-  data.table(
  simulation_index = 1:100,
  mild_DALY = rgamma(100, shape = unlist(DALY_weights[outcome=='mild','shape']), rate = unlist(DALY_weights[outcome=='mild','rate'])),
  fever_DALY = rgamma(100, shape = unlist(DALY_weights[outcome=='fever','shape']), rate = unlist(DALY_weights[outcome=='fever','rate'])),
  hosp_DALY = rgamma(100, shape = unlist(DALY_weights[outcome=='hospitalisation','shape']), rate = unlist(DALY_weights[outcome=='hospitalisation','rate']))
)

econ_cases_agg <- econ_cases_agg[DALY_weight_samples, on=c('simulation_index')]
econ_cases_agg[, mild_DALYs := flu_duration*mild_DALY*milds][, fever_DALYs := flu_duration*fever_DALY*fevers][, hosp_DALYs := flu_duration*hosp_DALY*hospitalisations]
econ_cases_agg[, c('mild_DALY','fever_DALY','hosp_DALY'):=NULL]

## healthcare costs from regression
print('healthcare costs')

cost_predic <- data.table(read_csv('econ/outcome_calculations/data/cost_predictor.csv', show_col_types=F))
cost_predic_c <- dcast(cost_predic, study_pop + gdpcap ~ outcome, value.var = 'pred_cost')
cost_predic_c[study_pop == 'adults', age_grp := 3]
cost_predic_c[study_pop == 'children', age_grp := 1]
cost_predic_c[study_pop == 'elderly', age_grp := 4]
cost_predic_c <- rbind(cost_predic_c, cost_predic_c[study_pop == 'adults',][,age_grp := 2])

econ_cases_agg[, gdpcap := floor(GDPpc)]

cost_predic_c_filter <- cost_predic_c[gdpcap %in% unique(econ_cases_agg$gdpcap)]
cost_predic_c_filter[, study_pop := NULL]
setnames(cost_predic_c_filter, 'hospital', 'hosp_cost')
setnames(cost_predic_c_filter, 'outpatient', 'outp_cost')

econ_cases_agg <- econ_cases_agg[cost_predic_c_filter, on=c('age_grp','gdpcap')]
econ_cases_agg[, total_hosp_cost := hosp_cost*hospitalisations]
if(outp_include == T){
  econ_cases_agg[, total_outp_cost := outp_cost*outpatients]
}
econ_cases_agg[, c('gdpcap','hosp_cost','outp_cost'):=NULL]

## adding YLLs and YLDs to make DALYs 
econ_cases_agg[, total_DALYs := mild_DALYs + fever_DALYs + hosp_DALYs + YLLs]

## discounting
print('discounting')

econ_cases_agg[, discount_year := year - 2025]
econ_cases_agg[, cost_discount_rate := (1 + cost_discount_rate_val)^(-discount_year)]
econ_cases_agg[, DALY_discount_rate := (1 + DALY_discount_rate_val)^(-discount_year)]
if(outp_include == T){
  econ_cases_agg[, discounted_costs := (total_hosp_cost + total_outp_cost)*cost_discount_rate] 
}else{
  econ_cases_agg[, discounted_costs := (total_hosp_cost)*cost_discount_rate] 
}
econ_cases_agg[, discounted_DALYs := total_DALYs*DALY_discount_rate]

econ_cases_agg[, c('discount_year','cost_discount_rate','DALY_discount_rate'):=NULL]

## calculating DALY*WTP
econ_cases_agg[, cost_of_DALYs := WTP_GDP_ratio*discounted_DALYs*GDPpc]

########################################################################

## loading vaccine doses
scenario_name <- 'base'

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_base/vacc_doses_", c_code, "_",
                                                    scenario_name, ".csv"), show_col_types=F) %>% 
                          mutate(cluster_code = c_code))
}

vacc_doses <- data.table(vacc_doses)

vacc_doses <- vacc_doses[,!c('country','cluster_code')]
vacc_doses_m <- melt(vacc_doses, id.vars = c('country_code','vacc_program','year'))
vacc_doses_m[, c('country_code','variable'):=NULL]
vacc_doses_m <- vacc_doses_m[, lapply(.SD, sum), by=c('vacc_program','year')]
vacc_doses_m[, ct := ceiling(vacc_program/5)]
vacc_doses_m[, vt := (vacc_program %% 5)]
vacc_doses_m[vt==0, vt := 5]
setnames(vacc_doses_m, 'value','doses')

## implicitly discounting vaccine prices at the cost-discount-rate
vacc_doses_m[, discount_year := year - 2025]
vacc_doses_m[, discount_rate := (1 + cost_discount_rate_val)^(-discount_year)]
vacc_doses_m[, discounted_doses := doses*discount_rate]

## total
vacc_doses_sum <- vacc_doses_m[, c('year','ct','vt','discounted_doses')][, lapply(.SD, sum), by=c('year','ct','vt')]
vacc_doses_sum <- rbind(vacc_doses_sum, data.table(year=2025:2054,ct='v',vt='',discounted_doses=0))

########################################################################

## TOTAL OUTCOMES, COSTS, AND DOSES IN EACH SCENARIO

total_out <- econ_cases_agg[,c('simulation_index','year','ct','vt','infections','deaths','hospitalisations',
                               'YLLs','YLL_cost','discounted_costs','discounted_DALYs','cost_of_DALYs')][, lapply(.SD, sum), 
                                                                                         by = c('simulation_index','year','ct','vt')]
## adding doses
total_out <- total_out[vacc_doses_sum, on=c('year','ct','vt')]

agg_out <- total_out[, lapply(.SD, sum), by=c('simulation_index','ct','vt')][, year:=NULL]

base_scenario <- agg_out[vt=='',c('simulation_index','discounted_costs','cost_of_DALYs')]
setnames(base_scenario, 'discounted_costs', 'base_costs'); setnames(base_scenario, 'cost_of_DALYs', 'base_DALYs')
agg_out <- agg_out[base_scenario, on='simulation_index']
agg_out[, incremental_costs := base_costs - discounted_costs][, incremental_DALYs := base_DALYs - cost_of_DALYs]
agg_out <- agg_out[!vt=='',]

agg_out[, threshold_price := (incremental_DALYs + incremental_costs)/discounted_doses]
agg_out[ct==1, ct_n:='0-5']
agg_out[ct==2, ct_n:='0-12']
agg_out[ct==3, ct_n:='0-18']
agg_out[ct==4, ct_n:='65+']
agg_out[ct==5, ct_n:='0-18, 65+']
agg_out$ct_n <- factor(agg_out$ct_n, levels = unique(agg_out$ct_n))

write_csv(agg_out, file = 'econ/agg_out.csv')

ggplot(agg_out) + 
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt))) +
  ylab('Threshold cost') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  ylim(c(0,NA)) +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Coverage level')

ggplot(agg_out) + 
  # geom_line(aes(x=incremental_DALYs/1000000, y=incremental_costs/1000000, group=simulation_index), alpha=0.5) + 
  geom_point(aes(x=incremental_DALYs/1000000, y=incremental_costs/1000000, col=as.factor(vt))) + 
  facet_grid(ct~., scales='fixed',
                        labeller = labeller(vt = supp.labs,
                                            ct = supp.labs.age)) +
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(color = 'Vaccine type') + 
  theme(text=element_text(size=14)) +
  theme_bw() + xlim(c(0,NA)) + ylim(c(0,NA)) +
  xlab('Incremental DALYs (millions)') + 
  ylab('Incremental costs (millions)')



########################################################################

## AGE-SPECIFIC OUTCOMES PLOT

as_out <- econ_cases_agg[,c('simulation_index','age_grp','ct','vt','infections','deaths','hospitalisations',
                               'fevers','milds')][, lapply(.SD, sum), by = c('simulation_index','ct','vt','age_grp')]

as_out_summ <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- as_out[, lapply(.SD, get(meas)), by = c('ct','vt','age_grp')]
  dt[, measure := meas]
  as_out_summ <- rbind(as_out_summ, dt)
}
as_out_summ[,simulation_index:=NULL]
as_out_summ_m <- melt(as_out_summ, id.vars = c('vt','ct','age_grp','measure'))
as_out_summ_c <- dcast(as_out_summ_m, vt + ct + age_grp + variable ~ measure, value.var = 'value')[order(age_grp)]
as_out_summ_c[age_grp == 1, age_grp_n := '0-4']
as_out_summ_c[age_grp == 2, age_grp_n := '5-19']
as_out_summ_c[age_grp == 3, age_grp_n := '20-64']
as_out_summ_c[age_grp == 4, age_grp_n := '65+']
as_out_summ_c[variable == 'milds', variable := 'mild infections']

ggplot(as_out_summ_c[ct==5]) + 
  geom_col(aes(x=age_grp_n, y=median/(30*1000000), fill=vt, group=vt),
           position='dodge', col=1) +
  geom_errorbar(aes(x = age_grp_n, ymin = eti95L/(30*1000000), ymax = eti95U/(30*1000000), group=vt),
                position=position_dodge(0.9)) +
  facet_grid(variable~., scales='free',
             labeller = labeller(vt = supp.labs,
                                        ct = supp.labs.age)) + 
  ylab('Mean annual incidence, in millions') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Age group') 





