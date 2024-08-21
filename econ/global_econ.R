#### Economic analysis - calculating outcomes
source("BS/BS_vaccine_programs.R")
source("BS/BS_colors.R")

library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
library(odin)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
suppressPackageStartupMessages(library(tidyverse))
library(viridis)
library(patchwork)
library(wpp2022)
library(WDI)

#### SETTING SCENARIO ####

scenario_name <- 'BASE50_DEPTH'
outp_include <- F # including outpatient/non-hospitalisation visits T/F
disease_modification <- F
mod_val <- 0.5
WTP_choice <- c('lancet','gdp')[1]
WTP_GDP_ratio <- 0.5 # proportion of GDP per capita for the willingness_to_pay threshold
discount_SA <- F

econ_folder_name <- paste0(ifelse(disease_modification, '_disease_mod',''),
                           ifelse(outp_include, '_outpatient',''),
                           ifelse((WTP_choice=='gdp'), '_gdp_',''),
                           ifelse((WTP_choice=='gdp'), WTP_GDP_ratio,''),
                           ifelse(discount_SA, '_discount0', ''))
print(paste0('Folder = ', scenario_name,econ_folder_name))

if((disease_modification+discount_SA+outp_include+(WTP_choice=='gdp'))>1){
  print('BTW there is more than one econ sensitivity analysis on')
}

## PARAMETERS
cost_discount_rate_val <- c(0.03, 0.03)[1 + discount_SA]
DALY_discount_rate_val <- c(0.03, 0)[1 + discount_SA]
flu_duration <- 4/365
wastage <- 0.1

print(paste0(scenario_name, ', outpatient = ', outp_include,
             ', disease_mod = ', disease_modification, ', wtp_thresh = ', WTP_choice,
             ifelse((WTP_choice=='gdp'), paste0(', GDP% = ', 100*WTP_GDP_ratio), ''),
             ', discount rates = ', 100*cost_discount_rate_val, '% & ', 100*DALY_discount_rate_val, '%'))

## making directory for saving outputs if doesn't exist
if(!file.exists(paste0('econ/', scenario_name,econ_folder_name))){
  dir.create(paste0('econ/', scenario_name,econ_folder_name))
  dir.create(paste0('econ/', scenario_name,econ_folder_name,'_plots'))
}

## LOADING DATA (BASE CASE)
print('loading data')

if(exists('econ_cases_agg')){rm(econ_cases_agg)};if(exists('econ_cases')){rm(econ_cases)}
for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0('data/vacc_output_', scenario_name,'/no-sync.nosync/econ_inp_', c_code, '.Rdata'))){
    print(paste0(c_code, ' exists'))
    load(paste0('data/vacc_output_', scenario_name,'/no-sync.nosync/econ_inp_', c_code, '.Rdata'))
    if(exists('econ_inp_new')){
      econ_inp <- copy(econ_inp_new)
      rm(econ_inp_new)
    }
    econ_cases_agg1 <- melt(econ_inp,
         id.vars=c('country_code','simulation_index','year','scenario'))
    if(disease_modification==T){
      econ_cases_agg1[,vaccinated:=(substr(variable,2,2)=='V')]
      for(scenario_i in unique(econ_cases_agg1$scenario)){
        print(scenario_i)
        print(100*sum(econ_cases_agg1[scenario==scenario_i & vaccinated==T,]$value)/sum(econ_cases_agg1[scenario==scenario_i]$value))
      }
    }
    econ_cases_agg1[,age_grp:=as.numeric(substr(variable, 3,3))][,variable:=NULL]
    if(disease_modification==T){
      econ_cases_agg1 <- econ_cases_agg1[, lapply(.SD,sum),
                                         by=c('country_code','simulation_index','year','scenario','vaccinated','age_grp')]
    }else{
      econ_cases_agg1 <- econ_cases_agg1[, lapply(.SD,sum),
                                         by=c('country_code','simulation_index','year','scenario','age_grp')]
    }
    setnames(econ_cases_agg1, 'value','infections')
    econ_cases_agg1[, ct := substr(scenario,4,4)]
    econ_cases_agg1[, vt := substr(scenario,9,9)]
    n_row <- nrow(econ_cases_agg1[infections<0])
    # if(n_row>0){print(n_row)}
    econ_cases_agg1[infections<0, infections:=0]
    if(exists('econ_cases_agg')){
      econ_cases_agg <- rbind(econ_cases_agg, econ_cases_agg1)
    }else{
      econ_cases_agg <- copy(econ_cases_agg1)
    }
  }
}

rm(econ_cases_agg1);rm(econ_inp)
econ_aside <- copy(econ_cases_agg)

#### IFRS ####

print('IFRs')
## load ifrs
national_ifrs <- data.table(read_csv('econ/outcome_calculations/data/national_ifrs.csv',
                                     show_col_types=F))

econ_cases_agg <- econ_cases_agg[national_ifrs[,c('country_code','simulation_index','age_grp','ifr')], on=c('country_code','simulation_index','age_grp')]
econ_cases_agg[,deaths := ifr*infections]
if(disease_modification==T){
  econ_cases_agg[vaccinated==T,deaths := mod_val*deaths]
}

#### IHRS ####

print('IHRs')
## load ihrs
global_ihrs <- data.table(read_csv('econ/outcome_calculations/data/global_ihrs.csv',
                                     show_col_types=F))
outpatient_ratios <- data.table(read_csv("econ/outcome_calculations/data/outpatient_ratios.csv", show_col_types=F))

econ_cases_agg <- econ_cases_agg[global_ihrs, on=c('simulation_index','age_grp')]
econ_cases_agg[,hospitalisations := ihr*infections]
if(disease_modification==T){
  econ_cases_agg[vaccinated==T,hospitalisations := mod_val*hospitalisations]
  econ_cases_agg[, vaccinated:=NULL]
  econ_cases_agg <- econ_cases_agg[, lapply(.SD,sum), by = c('country_code','simulation_index','year','scenario','age_grp','ct','vt')]
}

if(outp_include == T){
  print('outpatients')
  econ_cases_agg <- econ_cases_agg[outpatient_ratios, on=c('country_code','simulation_index'),
                                   ratio:=ratio]
  econ_cases_agg[, outpatients := ratio*hospitalisations]
  econ_cases_agg[, ratio:=NULL]
}

#### SYMPTOMATIC + FEVER ####

## from carrat (DOI: 10.1093/aje/kwm375)
print('symptomatic, fever')
# symp_probs <- data.table(
#   outcome = c('symptoms','fever'),
#   med = c(66.9, 34.9)/100,
#   l95 = c(58.3, 26.7)/100,
#   u95 = c(74.5, 44.2)/100
#   )
# ## determining distribution
# f.gamma <- function(shape, rate, x) {
#   p <- pgamma(x, shape, rate)
#   # return both
#   return(c(p))
# }
# delta <- function(fit, actual) sum((fit-actual)^2)
# objective <- function(theta, x, prob, ...) {
#   ab <- (theta)
#   fit <- f.gamma(ab[1], ab[2], x=as.numeric(x),...)
#   # fit <- f.beta(ab[1], ab[2], x=as.numeric(x),...)
#   return (delta(fit, prob))
# }
# fcn_fitting <- function(rates,
#                         probs){
#
#   x <- c(unlist(unname(rates)))
#   sol <- suppressWarnings(optim(f=objective,p=c(1,1),
#                                 # method="BFGS",
#                                 x=x,
#                                 prob=c(probs),
#                                 control = list(reltol = 1e-15)
#   ))
#   parms <- (sol$par)
#   return(parms)
# }
#
# for(i in 1:nrow(symp_probs)){
#   parms <- fcn_fitting(symp_probs[i,2:4], c(0.5, 0.025, 0.975))
#   symp_probs[i,"shape"] <- parms[1]
#   symp_probs[i,"rate"] <- parms[2]
#   symp_probs[i,"med_fit"] <- qgamma(p=c(0.5), shape=parms[1], rate=parms[2])
#   symp_probs[i,"l95_fit"] <- qgamma(p=c(0.025), shape=parms[1], rate=parms[2])
#   symp_probs[i,"u95_fit"] <- qgamma(p=c(0.975), shape=parms[1], rate=parms[2])
# }
#
# symp_samples <- data.table(
#   simulation_index = 1:100,
#   symp_prob = rgamma(100, shape = unlist(symp_probs[outcome=='symptoms','shape']), rate = unlist(symp_probs[outcome=='symptoms','rate'])),
#   fever_prob = rgamma(100, shape = unlist(symp_probs[outcome=='fever','shape']), rate = unlist(symp_probs[outcome=='fever','rate'])))
# write_csv(symp_samples, file='econ/outcome_calculations/data/symp_samples.csv')

symp_samples <- data.table(read_csv('econ/outcome_calculations/data/symp_samples.csv', show_col_types=F))

econ_cases_agg <- econ_cases_agg[symp_samples, on=c('simulation_index')]
econ_cases_agg[, symptomatics := symp_prob*infections][, fevers := fever_prob*infections][, non_fevers := symptomatics - fevers]


#### YLLS ####

## discounting at same rate as DALYs, using 'lxqx' method

print('YLLs')

yll_df <- national_ifrs[simulation_index==1, c('country_code','age_grp','cluster_name')]
source('econ/outcome_calculations/calc_ylls.R')
pb <- txtProgressBar(min = 0, max = 186, style = 3, width = 50, char = "=")

for(i in 1:length(unique(yll_df$country_code))){
  iso3c_i <- unique(yll_df$country_code)[i]
  yll_df[country_code == iso3c_i, 'yll'] <- yll(LT = UNLT[ISO3_code == iso3c_i & MidPeriod == 2022.5],
                                                r = DALY_discount_rate_val,
                                                smr = 1,
                                                weight_method = "lxqx", # weight method to average LE by age group: "lx" "lxqx" "equal" "pop_ifr"
                                                model_ages = c(5,20,65))$d_LEx
  setTxtProgressBar(pb, i)
}

econ_cases_agg <- econ_cases_agg[yll_df, on=c('country_code','age_grp'), yll := yll]
econ_cases_agg[,YLLs := yll*deaths]

## removing probabilities etc.
econ_cases_agg[, c('ifr','ihr','symp_prob','fever_prob','yll'):=NULL]

#### YLDS ####
print('YLDs')

## weights from GBD
# DALY_weights <- data.table(
#   outcome = c('non_fever','fever','hospitalisation'),
#   med = c(0.006, 0.051, 0.133),
#   l95 = c(0.002, 0.032, 0.088),
#   u95 = c(0.012, 0.074, 0.190)
# )
# for(i in 1:nrow(DALY_weights)){
#   parms <- fcn_fitting(DALY_weights[i,2:4], c(0.5, 0.025, 0.975))
#   DALY_weights[i,"shape"] <- parms[1]
#   DALY_weights[i,"rate"] <- parms[2]
#   DALY_weights[i,"med_fit"] <- qgamma(p=c(0.5), shape=parms[1], rate=parms[2])
#   DALY_weights[i,"l95_fit"] <- qgamma(p=c(0.025), shape=parms[1], rate=parms[2])
#   DALY_weights[i,"u95_fit"] <- qgamma(p=c(0.975), shape=parms[1], rate=parms[2])
# }
# DALY_weight_samples <-  data.table(
#   simulation_index = 1:100,
#   non_fever_DALY = rgamma(100, shape = unlist(DALY_weights[outcome=='non_fever','shape']), rate = unlist(DALY_weights[outcome=='non_fever','rate'])),
#   fever_DALY = rgamma(100, shape = unlist(DALY_weights[outcome=='fever','shape']), rate = unlist(DALY_weights[outcome=='fever','rate'])),
#   hosp_DALY = rgamma(100, shape = unlist(DALY_weights[outcome=='hospitalisation','shape']), rate = unlist(DALY_weights[outcome=='hospitalisation','rate']))
# )
# write_csv(DALY_weight_samples, file='econ/outcome_calculations/data/DALY_weight_samples.csv')

DALY_weight_samples <- data.table(read_csv('econ/outcome_calculations/data/DALY_weight_samples.csv', show_col_types=F))

econ_cases_agg <- econ_cases_agg[DALY_weight_samples, on=c('simulation_index')]
econ_cases_agg[, non_fever_DALYs := flu_duration*non_fever_DALY*non_fevers][, fever_DALYs := flu_duration*fever_DALY*fevers][, hosp_DALYs := flu_duration*hosp_DALY*hospitalisations]
econ_cases_agg[, c('non_fever_DALY','fever_DALY','hosp_DALY'):=NULL]

## adding YLLs and YLDs to make DALYs
econ_cases_agg[, total_DALYs := non_fever_DALYs + fever_DALYs + hosp_DALYs + YLLs]

#### HEALTHCARE COSTS ####
print('healthcare costs')

load('econ/outcome_calculations/data/predicted_costs')
cost_predic_c <- dcast(pred_costs[,c('iso3c','outcome','study_pop',
                                     'simulation_index','gdpcap','sample_cost')], iso3c + study_pop + gdpcap + simulation_index ~ outcome, value.var = 'sample_cost')
cost_predic_c[study_pop == 'adults', age_grp := 3]
cost_predic_c[study_pop == 'children', age_grp := 1]
cost_predic_c[study_pop == 'elderly', age_grp := 4]
cost_predic_c <- rbind(cost_predic_c, cost_predic_c[study_pop == 'adults',][,age_grp := 2])
cost_predic_c[, study_pop := NULL]
setnames(cost_predic_c, 'hospital', 'hosp_cost')
setnames(cost_predic_c, 'outpatient', 'outp_cost')
setnames(cost_predic_c, 'iso3c', 'country_code')

econ_cases_agg <- econ_cases_agg[cost_predic_c, on=c('age_grp','country_code','simulation_index')]
econ_cases_agg[, total_hosp_cost := hosp_cost*hospitalisations]
if(outp_include == T){
  econ_cases_agg[, total_outp_cost := outp_cost*outpatients]
}
econ_cases_agg[, c('hosp_cost','outp_cost'):=NULL]

#### DISCOUNTING ####
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

## adding WTP
wtp_thresh <- data.table(read_csv('econ/outcome_calculations/data/WTP_thresholds.csv', show_col_type=F))
setnames(wtp_thresh, 'iso3c','country_code')
if(WTP_choice == 'lancet'){
  econ_cases_agg <- econ_cases_agg[wtp_thresh, on=c('country_code'), cet := cet]
  econ_cases_agg[, cost_of_DALYs := discounted_DALYs*cet]
  econ_cases_agg[, cet := NULL]
}
if(WTP_choice == 'gdp'){
  econ_cases_agg[, cost_of_DALYs := WTP_GDP_ratio*discounted_DALYs*gdpcap]
}


#### VACCINE DOSES ####
print('vaccine doses')

if(grepl('0',scenario_name)){
  vacc_scenario_name <- 'BASE50'
  if(scenario_name == 'LOW20'){
    vacc_scenario_name <- scenario_name
  }
}else{
  vacc_scenario_name <- 'base'
  if(scenario_name == 'low_cov'){
    vacc_scenario_name <- scenario_name
  }
}; print(paste0('Vaccination doses scenario: ', vacc_scenario_name))

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_", vacc_scenario_name, "/vacc_doses_", c_code, "_",
                                                  vacc_scenario_name, ".csv"), show_col_types=F) %>%
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

if(scenario_name %like% '(?i)breadth'){
  vacc_doses_m <- rbind(vacc_doses_m[vacc_doses_m$vt==1,],
                        vacc_doses_m[vacc_doses_m$vt==1,],
                        vacc_doses_m[vacc_doses_m$vt==1,],
                        vacc_doses_m[vacc_doses_m$vt==1,],
                        vacc_doses_m[vacc_doses_m$vt==1,])
  vacc_doses_m$vt <- rep(1:5, each = nrow(vacc_doses_m)/5)
}

## implicitly discounting vaccine prices at the cost discount rate
vacc_doses_m[, discount_year := year - 2025]
vacc_doses_m[, discount_rate := (1 + cost_discount_rate_val)^(-discount_year)]
vacc_doses_m[, discounted_doses := doses*discount_rate]

## total
vacc_doses_sum <- vacc_doses_m[, c('year','ct','vt','doses','discounted_doses')][, lapply(.SD, sum), by=c('year','ct','vt')]
vacc_doses_sum <- rbind(vacc_doses_sum, data.table(year=2025:2054,ct='v',vt='',doses=0,discounted_doses=0))

## national
vacc_doses_nat <- melt(vacc_doses, id.vars = c('country_code','vacc_program','year'))
vacc_doses_nat[, c('variable'):=NULL]
vacc_doses_nat <- vacc_doses_nat[, lapply(.SD, sum), by=c('country_code','vacc_program','year')]
vacc_doses_nat[, ct := ceiling(vacc_program/5)]
vacc_doses_nat[, vt := (vacc_program %% 5)]
vacc_doses_nat[vt==0, vt := 5]
if(scenario_name %like% '(?i)breadth'){
  vacc_doses_nat <- rbind(vacc_doses_nat[vacc_doses_nat$vt==1,],
                          vacc_doses_nat[vacc_doses_nat$vt==1,],
                          vacc_doses_nat[vacc_doses_nat$vt==1,],
                          vacc_doses_nat[vacc_doses_nat$vt==1,],
                          vacc_doses_nat[vacc_doses_nat$vt==1,])
  vacc_doses_nat$vt <- rep(1:5, each = nrow(vacc_doses_nat)/5)
}
vacc_doses_nat[,vt:=as.character(vt)][,ct:=as.character(ct)]
setnames(vacc_doses_nat, 'value','doses')
vt0 <- vacc_doses_nat[ct==1&vt==1,]
vt0[, vt:=''][, ct:='v'][,doses:=0][,vacc_program:=0]
vacc_doses_nat <- rbind(vacc_doses_nat, vt0)
vacc_doses_nat[, discount_year := year - 2025]
vacc_doses_nat[, discount_rate := (1 + cost_discount_rate_val)^(-discount_year)]
vacc_doses_nat[, discounted_doses := doses*discount_rate]

#### GLOBAL OUTCOMES, COSTS, AND DOSES ####
print('global outcomes')

total_out <- econ_cases_agg[,c('simulation_index','year','ct','vt','infections','deaths','hospitalisations',
                               'YLLs','total_DALYs','discounted_costs','discounted_DALYs','cost_of_DALYs')][, lapply(.SD, sum),
                                                                                         by = c('simulation_index','year','ct','vt')]
## adding doses
total_out <- total_out[vacc_doses_sum, on=c('year','ct','vt')]

agg_out <- total_out[, lapply(.SD, sum), by=c('simulation_index','ct','vt')][, year:=NULL]

base_scenario <- agg_out[vt=='',c('simulation_index','discounted_costs','discounted_DALYs','cost_of_DALYs')]
setnames(base_scenario, 'discounted_costs', 'base_costs'); setnames(base_scenario, 'cost_of_DALYs', 'base_DALYs_cost')
setnames(base_scenario, 'discounted_DALYs', 'base_DALYs')
agg_out <- agg_out[base_scenario, on='simulation_index']
agg_out[, incremental_costs := base_costs - discounted_costs][, incremental_DALYs := base_DALYs_cost - cost_of_DALYs]
agg_out <- agg_out[!vt=='',]

# agg_out[, threshold_price := (incremental_DALYs + incremental_costs)/discounted_doses] # removing as delivery not costed

agg_out[ct==1, ct_n:='0-4']
agg_out[ct==2, ct_n:='0-10']
agg_out[ct==3, ct_n:='0-17']
agg_out[ct==4, ct_n:='65+']
agg_out[ct==5, ct_n:='0-17, 65+']

save(agg_out, file = paste0('econ/',scenario_name,econ_folder_name,'/agg_out'))

#### WHO REGION OUTCOMES, COSTS, AND DOSES ####
print('region-specific outcomes')

total_who <- econ_cases_agg[,c('country_code','simulation_index','year','ct','vt','infections','deaths','hospitalisations',
                               'YLLs','discounted_costs','discounted_DALYs','cost_of_DALYs')][, lapply(.SD, sum),
                                                                                              by = c('country_code','simulation_index','year','ct','vt')]

## adding doses
total_who <- total_who[vacc_doses_nat, on=c('country_code','year','ct','vt')]

## adding WHO regions
who_regions <- data.table(read_csv('econ/outcome_calculations/data/WHO_regions.csv', show_col_types = F))
who_regions <- who_regions[country_code %in% total_who$country_code,]
total_who <- total_who[who_regions[,c('country_code','WHOREGION')], on=c('country_code')]

total_who_sum <- copy(total_who)
total_who_sum[,c('country_code','year','vacc_program'):=NULL]
total_who_sum <- total_who_sum[, lapply(.SD, sum, na.rm=T), by=c('WHOREGION','simulation_index','ct','vt')]

base_scenario <- total_who_sum[vt=='',c('WHOREGION','simulation_index','discounted_costs','discounted_DALYs','cost_of_DALYs')]
setnames(base_scenario, 'discounted_costs', 'base_costs'); setnames(base_scenario, 'cost_of_DALYs', 'base_DALYs_cost')
setnames(base_scenario, 'discounted_DALYs', 'base_DALYs')
total_who_sum <- total_who_sum[base_scenario, on=c('simulation_index','WHOREGION')]
total_who_sum[, incremental_costs := base_costs - discounted_costs][, incremental_DALYs := base_DALYs_cost - cost_of_DALYs]
total_who_sum <- total_who_sum[!vt=='',]

# total_who_sum[, threshold_price := (incremental_DALYs + incremental_costs)/discounted_doses] # removing as delivery not costed

total_who_sum[ct==1, ct_n:='0-4']
total_who_sum[ct==2, ct_n:='0-10']
total_who_sum[ct==3, ct_n:='0-17']
total_who_sum[ct==4, ct_n:='65+']
total_who_sum[ct==5, ct_n:='0-17, 65+']

save(total_who_sum, file = paste0('econ/',scenario_name,econ_folder_name,'/total_who_sum'))

#### NATIONAL OUTCOMES, COSTS, AND DOSES IN EACH SCENARIO ####
print('national outcomes')

total_nat_sum <- copy(total_who)
total_nat_sum[,c('year','vacc_program'):=NULL]
total_nat_sum <- total_nat_sum[, lapply(.SD, sum, na.rm=T), by=c('country_code','WHOREGION','ct','vt','simulation_index')]

base_scenario <- total_nat_sum[vt=='',c('country_code','simulation_index','discounted_costs','discounted_DALYs','cost_of_DALYs')]
setnames(base_scenario, 'discounted_costs', 'base_costs'); setnames(base_scenario, 'cost_of_DALYs', 'base_DALYs_cost')
setnames(base_scenario, 'discounted_DALYs', 'base_DALYs')
total_nat_sum <- total_nat_sum[base_scenario, on=c('simulation_index','country_code')]
total_nat_sum[, incremental_costs := base_costs - discounted_costs][, incremental_DALYs := base_DALYs_cost - cost_of_DALYs]
total_nat_sum <- total_nat_sum[!vt=='',]

# costs of doses
delivery_cost_samples <- data.table(read_csv('econ/outcome_calculations/data/delivery_cost_samples.csv', show_col_types=F))
setnames(delivery_cost_samples, 'iso3c','country_code')
total_nat_sum <- total_nat_sum[delivery_cost_samples, on=c('country_code','simulation_index')]

# threshold price, with wastage and delivery cost taken into account
total_nat_sum[, threshold_price := (incremental_DALYs + incremental_costs - discounted_doses*delivery_cost)/(discounted_doses/(1-wastage))]

total_nat_sum[ct==1, ct_n:='0-4']
total_nat_sum[ct==2, ct_n:='0-10']
total_nat_sum[ct==3, ct_n:='0-17']
total_nat_sum[ct==4, ct_n:='65+']
total_nat_sum[ct==5, ct_n:='0-17, 65+']

total_nat_sum[,c('discount_year','discount_rate'):=NULL]

save(total_nat_sum, file = paste0('econ/',scenario_name,econ_folder_name,'/total_nat_sum'))

#### AGE-SPECIFIC OUTCOMES ####
print('age-specific outcomes')

as_out <- econ_cases_agg[,c('simulation_index','age_grp','ct','vt','infections','deaths','hospitalisations',
                               'fevers','non_fevers','total_DALYs')][, lapply(.SD, sum), by = c('simulation_index','ct','vt','age_grp')]

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

save(as_out_summ_c, file=paste0('econ/',scenario_name,econ_folder_name,'/as_out_summ_c'))

as_nv <- as_out[ct=='v']
as_nv[, c('ct','vt'):=NULL]
setnames(as_nv, 'infections', 'nv_infections')
setnames(as_nv, 'deaths', 'nv_deaths')
setnames(as_nv, 'hospitalisations', 'nv_hospitalisations')
setnames(as_nv, 'fevers', 'nv_fevers')
setnames(as_nv, 'non_fevers', 'nv_non_fevers')
setnames(as_nv, 'total_DALYs', 'nv_total_DALYs')
as_out_av <- as_out[as_nv, on=c('simulation_index','age_grp')]
as_out_av[, infections_av := nv_infections - infections][, deaths_av := nv_deaths - deaths][, hospitalisations_av := nv_hospitalisations - hospitalisations]
as_out_av[, fevers_av := nv_fevers - fevers][, non_fevers_av := nv_non_fevers - non_fevers][, total_DALYs_av := nv_total_DALYs - total_DALYs]

out_av <- copy(as_out_av)
out_av[, age_grp := NULL]
out_av <- out_av[, lapply(.SD, sum), by = c('ct','vt','simulation_index')]
save(out_av, file=paste0('econ/',scenario_name,econ_folder_name,'/as_out_averted_100'))

as_out_av[, c('infections','deaths','hospitalisations','fevers','non_fevers','total_DALYs',
              'nv_infections','nv_deaths','nv_hospitalisations','nv_fevers','nv_non_fevers','nv_total_DALYs') := NULL]
as_out_av <- as_out_av[!ct=='v']
as_out_summ <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- as_out_av[, lapply(.SD, get(meas)), by = c('ct','vt','age_grp')]
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

save(as_out_summ_c, file=paste0('econ/',scenario_name,econ_folder_name,'/as_out_averted'))

as_out <- econ_cases_agg[,c('country_code','simulation_index','age_grp','ct','vt','infections','deaths','hospitalisations',
                  'fevers','non_fevers','total_DALYs')][, lapply(.SD, sum), by = c('country_code','simulation_index','ct','vt','age_grp')]
who_regions <- data.table(read_csv('econ/outcome_calculations/data/WHO_regions.csv', show_col_types = F))
who_regions <- who_regions[country_code %in% as_out$country_code,]
as_out <- as_out[who_regions[,c('country_code','WHOREGION')], on='country_code']
as_out[, country_code:=NULL]
as_out <- as_out[, lapply(.SD, sum), by = c('WHOREGION','simulation_index','ct','vt','age_grp')]
as_nv <- as_out[ct=='v']
as_nv[, c('ct','vt'):=NULL]
setnames(as_nv, 'infections', 'nv_infections')
setnames(as_nv, 'deaths', 'nv_deaths')
setnames(as_nv, 'hospitalisations', 'nv_hospitalisations')
setnames(as_nv, 'fevers', 'nv_fevers')
setnames(as_nv, 'non_fevers', 'nv_non_fevers')
setnames(as_nv, 'total_DALYs', 'nv_total_DALYs')
as_out_av <- as_out[as_nv, on=c('WHOREGION','simulation_index','age_grp')]
as_out_av[, infections_av := nv_infections - infections][, deaths_av := nv_deaths - deaths][, hospitalisations_av := nv_hospitalisations - hospitalisations]
as_out_av[, fevers_av := nv_fevers - fevers][, non_fevers_av := nv_non_fevers - non_fevers][, total_DALYs_av := nv_total_DALYs - total_DALYs]
as_out_av[, c('infections','deaths','hospitalisations','fevers','non_fevers','total_DALYs',
              'nv_infections','nv_deaths','nv_hospitalisations','nv_fevers','nv_non_fevers','nv_total_DALYs') := NULL]
as_out_av <- as_out_av[!ct=='v']
as_out_summ <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- as_out_av[, lapply(.SD, get(meas)), by = c('WHOREGION','ct','vt','age_grp')]
  dt[, measure := meas]
  as_out_summ <- rbind(as_out_summ, dt)
}
as_out_summ[,simulation_index:=NULL]
as_out_summ_m <- melt(as_out_summ, id.vars = c('WHOREGION','vt','ct','age_grp','measure'))
as_out_summ_c <- dcast(as_out_summ_m, WHOREGION + vt + ct + age_grp + variable ~ measure, value.var = 'value')[order(age_grp)]
as_out_summ_c[age_grp == 1, age_grp_n := '0-4']
as_out_summ_c[age_grp == 2, age_grp_n := '5-19']
as_out_summ_c[age_grp == 3, age_grp_n := '20-64']
as_out_summ_c[age_grp == 4, age_grp_n := '65+']

save(as_out_summ_c, file=paste0('econ/',scenario_name,econ_folder_name,'/as_out_averted_WHOREGION'))



