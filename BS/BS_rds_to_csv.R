## PLOTTING PROJECTION OUTPUTS
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_vaccine_programs.R")

library(readr)
library(dplyr)
library(data.table)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(bayestestR)

# ci50L <- function(x){return(unname(unlist(ci(x, ci=0.50)[2])))}
# ci50U <- function(x){return(unname(unlist(ci(x, ci=0.50)[3])))}
# ci95L <- function(x){return(unname(unlist(ci(x, ci=0.95)[2])))}
# ci95U <- function(x){return(unname(unlist(ci(x, ci=0.95)[3])))}
eti50L <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[25] + 0.75*sort(x)[26])
    }else{stop('length(x) != 100')}
}
eti50U <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[76] + 0.75*sort(x)[75])
  }else{stop('length(x) != 100')}
}
eti95L <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[3] + 0.475*sort(x)[4])
  }else{stop('length(x) != 100')}
}
eti95U <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[98] + 0.475*sort(x)[97])
  }else{stop('length(x) != 100')}
}

scenario_name <- c('none', 'base', 'low_cov', 'rel_inf', 'depth', 'breadth', 'direct')[6]

start_time <- Sys.time()

c_number <- 1
c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
                        "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
                        "Southern America")[c_number]
c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]
print(paste0(c_code, ', ', scenario_name))

# load no vacc cases
itz_cases_no_vacc <- data.table(readRDS(paste0("data/vacc_output_base/vacc_", c_code, '_none_ct_1.rds'))[[1]])

itz_cases_no_vacc_y <- copy(itz_cases_no_vacc)
itz_cases_no_vacc_y <- itz_cases_no_vacc_y[, year := year(week)]
itz_cases_no_vacc_y[, c("country","week") := NULL]  
no_vacc <- itz_cases_no_vacc_y[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'simulation_index', 'year')]
no_vacc[, scenario := 'no_vacc']
no_vacc_m <- no_vacc[, lapply(.SD, median, na.rm=T), by = c('country_code', 'year', 'scenario')]
no_vacc_m[, measure := 'median']
for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
  no_vacc_m2 <- no_vacc[, lapply(.SD, get(meas)), by = c('country_code', 'year', 'scenario')]
  no_vacc_m2[, measure := meas]
  no_vacc_m <- rbind(no_vacc_m, no_vacc_m2)
}
no_vacc_m[, simulation_index := NULL]

nat_ann <- copy(no_vacc_m)

itz_cases_no_vacc_w <- itz_cases_no_vacc[, c('country','country_code') := NULL]
global_weekly_nv <- itz_cases_no_vacc_w[, lapply(.SD, sum, na.rm=T), by = c('simulation_index', 'week')]
global_weekly_nv[, scenario := 'no_vacc']

global_weekly <- copy(global_weekly_nv)

print('no_vacc done')
print(Sys.time() - start_time)

for(ct in 1:5){
  list <-  readRDS(paste0("data/vacc_output_", scenario_name, "/no-sync.nosync/vacc_", c_code, '_', scenario_name, '_ct_', ct,'.rds'))
  for(vt in 1:5){
    itz_cases <- data.table(list[[vt]])
    
    itz_cases_y <- copy(itz_cases)
    itz_cases_y <- itz_cases_y[, year := year(week)]
    itz_cases_y[, c("country","week") := NULL]  
    dt <- itz_cases_y[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'simulation_index', 'year')]
    dt[, scenario := paste0('ct_',ct,'_vt_',vt)]
    dt_m <- dt[, lapply(.SD, median, na.rm=T), by = c('country_code', 'year', 'scenario')]
    dt_m[, measure := 'median']
    for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
      dt_m2 <- dt[, lapply(.SD, get(meas)), by = c('country_code', 'year', 'scenario')]
      dt_m2[, measure := meas]
      dt_m <- rbind(dt_m, dt_m2)
    }
    dt_m[, simulation_index := NULL]
    nat_ann <- rbind(nat_ann, dt_m)
    
    itz_cases_w <- itz_cases[, c('country','country_code') := NULL]
    global_weekly_dt <- itz_cases_w[, lapply(.SD, sum, na.rm=T), by = c('simulation_index', 'week')]
    global_weekly_dt[, scenario := paste0('ct_',ct,'_vt_',vt)]
    
    global_weekly <- rbind(global_weekly, global_weekly_dt)
    
    print(paste0(ct, ' ', vt, ' done'))
    print(Sys.time() - start_time)
  }
}

nat_ann_c <- copy(nat_ann)
global_weekly_c <- copy(global_weekly)

save(nat_ann_c, file = paste0('data/vacc_output_',scenario_name,'/nat_ann_', c_code, '.Rdata'))
save(global_weekly_c, file = paste0('data/vacc_output_',scenario_name,'/global_weekly_', c_code, '.Rdata'))

print(Sys.time() - start_time)














