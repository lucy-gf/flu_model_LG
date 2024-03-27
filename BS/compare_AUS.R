
## PLOTTING PROJECTION OUTPUTS
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_vaccine_programs.R")
source("BS/BS_colors.R")

library(readr)
library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(data.table)
library(odin)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(bayestestR)
library(viridis)
library(patchwork)
library(countrycode)

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

c_number <- 6
c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
            "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
            "Southern America")[c_number]
c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]

# load no vacc cases
itz_cases_no_vacc <- data.table(readRDS(paste0("data/vacc_output_scaling/vacc_", c_code, '_none_ct_1.rds'))[[1]])

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

for(ct in 1:5){
  list <-  readRDS(paste0("data/vacc_output_scaling/vacc_", c_code, '_', scenario_name, '_ct_', ct,'.rds'))
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




load(paste0('data/vacc_output/nat_ann_', c_code, '.Rdata'))
load(paste0('data/vacc_output/global_weekly_', c_code, '.Rdata'))

global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]
global_weekly_cum <- copy(global_weekly_m)
global_weekly_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario)]
global_weekly_cum[, value:=NULL]
global_weekly_nv <- global_weekly_cum[scenario=='no_vacc',]
global_weekly_cum[global_weekly_nv, on = c("simulation_index","week"), nv_cum := i.cum.sum]
global_weekly_cum[,averted := nv_cum-cum.sum]
global_weekly_cum[,nv_cum:=NULL]
global_weekly_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- global_weekly_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario')]
  dt[, measure := meas]
  global_weekly_meas <- rbind(global_weekly_meas, dt)
}
global_weekly_meas[,simulation_index:=NULL]
global_cumulative <- copy(global_weekly_meas)
global_cumulative[,averted:=NULL]
global_cumulative_wide <- dcast.data.table(global_cumulative, 
                                           week+scenario~measure,
                                           value.var = "cum.sum")
global_cumulative_wide[, ct := substr(scenario, 4,4)]
global_cumulative_wide[, vt := substr(scenario, 9,9)]

ggplot(global_cumulative_wide[!scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
  scale_fill_manual(values = vt_colors) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Cumulative cases, millions') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14))

orig_m <- data.table(melt(global_weekly_c, id.vars = c("simulation_index", "week", "scenario")))
orig_m[, variable:=NULL]
orig_m <- orig_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]
orig_cum <- copy(orig_m)
orig_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario)]
orig_cum[, value:=NULL]
orig_nv <- orig_cum[scenario=='no_vacc',]
orig_cum[orig_nv, on = c("simulation_index","week"), nv_cum := i.cum.sum]
orig_cum[,averted := nv_cum-cum.sum]
orig_cum[,nv_cum:=NULL]
orig_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- orig_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario')]
  dt[, measure := meas]
  orig_meas <- rbind(orig_meas, dt)
}
orig_meas[,simulation_index:=NULL]
orig_cumulative <- copy(orig_meas)
orig_cumulative[,averted:=NULL]
orig_cumulative_wide <- dcast.data.table(orig_cumulative, 
                                           week+scenario~measure,
                                           value.var = "cum.sum")
orig_cumulative_wide[, ct := substr(scenario, 4,4)]
orig_cumulative_wide[, vt := substr(scenario, 9,9)]

ggplot(orig_cumulative_wide[!scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
  scale_fill_manual(values = vt_colors) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Cumulative cases, millions') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14))


orig_df <- orig_cumulative_wide[orig_cumulative_wide$week ==
                       max(orig_cumulative_wide$week)]
new_df <- global_cumulative_wide[orig_cumulative_wide$week ==
                                 max(global_cumulative_wide$week)]

orig_df[,3:7] - new_df[,3:7]

## JUST NZL

## NEW (SCALED) DATA
none <- data.table(readRDS(paste0("data/vacc_output_scaling/vacc_", c_code, '_none_ct_1.rds'))[[1]])
cases <- none[none$country_code == 'NZL',]
cases[,scenario:='no_vacc']
for(ct in 1:5){
  list <-  readRDS(paste0("data/vacc_output_scaling/vacc_", c_code, '_', scenario_name, '_ct_', ct,'.rds'))
  for(vt in 1:5){
    itz_cases <- data.table(list[[vt]])
    itz_cases <- itz_cases[itz_cases$country_code == 'NZL',]
    itz_cases[,scenario:=paste0('ct_',ct,'_vt_',vt)]
    cases <- rbind(cases, itz_cases)
  }
}

## ORIGINAL DATA
none <- data.table(readRDS(paste0("data/vacc_output/vacc_", c_code, '_none_ct_1.rds'))[[1]])
orig_cases <- none[none$country_code == 'NZL',]
orig_cases[,scenario:='no_vacc']
for(ct in 1:5){
  list <-  readRDS(paste0("data/vacc_output/vacc_", c_code, '_', scenario_name, '_ct_', ct,'.rds'))
  for(vt in 1:5){
    itz_cases <- data.table(list[[vt]])
    itz_cases <- itz_cases[itz_cases$country_code == 'NZL',]
    itz_cases[,scenario:=paste0('ct_',ct,'_vt_',vt)]
    orig_cases <- rbind(orig_cases, itz_cases)
  }
}

cases[,c('country','country_code'):=NULL]
orig_cases[,c('country','country_code'):=NULL]
cases_m <- data.table(melt(cases, id.vars = c("simulation_index", "week", "scenario")))
cases_m[, variable:=NULL]
cases_m <- cases_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]
cases_m[, version:='old']
cases_m2 <- data.table(melt(orig_cases, id.vars = c("simulation_index", "week", "scenario")))
cases_m2[, variable:=NULL]
cases_m2 <- cases_m2[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]
cases_m2[, version:='new']
cases_m <- rbind(cases_m, cases_m2)
cases_m[,ct:=substr(scenario,4,4)]
cases_m[,vt:=substr(scenario,9,9)]

ggplot(cases_m[year(cases_m$week)==2035 & 
                 cases_m$simulation_index==1 &
                 # cases_m$version=='old'&
                 !cases_m$scenario == 'no_vacc',]) + 
  geom_line(aes(x=week, y=value, group=interaction(simulation_index, version), col=version), 
            alpha=0.5) +
  facet_grid(ct~vt, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  theme_minimal()

cases_m[, cum.sum := cumsum(value), by=list(simulation_index, scenario, version)]

ggplot(cases_m[!cases_m$scenario == 'no_vacc',]) + 
  geom_line(aes(x=week, y=cum.sum, group=interaction(simulation_index, version), col=version), 
            alpha=0.2) +
  facet_grid(ct~vt, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  theme_minimal()

cases_cum <- copy(cases_m)
cases_cum[, c('value','ct','vt'):=NULL]
cases_cum_meas <- cases_cum[, lapply(.SD, median, na.rm=T), by = c('version', 'week', 'scenario')]
cases_cum_meas[, measure := 'median']
for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- cases_cum[, lapply(.SD, get(meas)), by = c('version', 'week', 'scenario')]
  dt[, measure := meas]
  cases_cum_meas <- rbind(cases_cum_meas, dt)
}
cases_cum_meas[, simulation_index := NULL]

cases_meas_wide <- dcast.data.table(cases_cum_meas, 
                 week+scenario+version~measure,
                 value.var = "cum.sum")
cases_meas_wide[,ct:=substr(scenario,4,4)]
cases_meas_wide[,vt:=substr(scenario,9,9)]

ggplot(cases_meas_wide[!cases_meas_wide$scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, 
                  fill = version, group=version), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                  fill = version, group = version), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000, group = version), lwd=0.6) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Cumulative influenza infections, millions') +
  theme(text=element_text(size=14)) + ggtitle('New Zealand')

ggplot(cases_meas_wide[!cases_meas_wide$scenario=='no_vacc'&
                         year(cases_meas_wide$week)==2054 &
                         cases_meas_wide$ct==1 & cases_meas_wide$vt==1,]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000,
                  fill = version, group=version), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                  fill = version, group = version), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000, group = version), lwd=0.6) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Cumulative influenza infections, millions') +
  theme(text=element_text(size=14)) + ggtitle('New Zealand')

cases_meas_wide[cases_meas_wide$week == max(cases_meas_wide$week)]














