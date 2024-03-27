
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

cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/exemplar_output/exemplar_", cntr, 
                                             "_2010_2019.rds"))[[1]]))
}

cases_dt[,country:=NULL]
cases_dt_m <- data.table(melt(cases_dt, id.vars = c("country_code", "simulation_index", "week")))
cases_dt_m[grepl('A', variable), strain := "A"]
cases_dt_m[grepl('B', variable), strain := "B"]
cases_tot <- copy(cases_dt_m)
cases_tot[,variable:=NULL]
cases_tot <- cases_tot[, lapply(.SD, sum, na.rm=T), by=c('country_code','simulation_index','week','strain')]

ggplot(cases_tot) +
  geom_line(aes(x=week, y=value/1000000, group=simulation_index, col=strain), alpha=0.2) +
  theme_minimal() + ylab('Cases, millions') + 
  scale_color_manual(values = strain_colors1) +
  facet_grid(country_code~strain, scales='free')

ggplot(cases_tot[cases_tot$country_code=='CAN'&
                   cases_tot$strain=='A'&
                        year(cases_tot$week)%in%2016:2018]) +
  geom_line(aes(x=week, y=value, group=simulation_index), alpha=0.1) +
  theme_minimal()
  xlab('Year') + ylab('Influenza infections, millions') +
  theme(text=element_text(size=14)) + 
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) 


cases_cis <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- cases_tot[, lapply(.SD, get(meas)), by = c('country_code', 'week', 'strain')]
  dt[, measure := meas]
  cases_cis <- rbind(cases_cis, dt)
}
cases_cis[,simulation_index:=NULL]
cases_cis_wide <- dcast.data.table(cases_cis, 
                                   country_code+week+strain~measure,
                                   value.var = "value")
  
ggplot(cases_cis_wide) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, 
                  fill = strain), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                  fill = strain), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000, 
                col = strain), lwd=0.6) +
  scale_color_manual(values = strain_colors1) +
  scale_fill_manual(values = strain_colors1) +
  theme_bw() + facet_grid(country_code~strain, scales='free') +
  xlab('Year') + ylab('Influenza infections, millions') +
  theme(text=element_text(size=14)) + labs(fill = 'Strain', col='Strain') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) 

cases_cum <- copy(cases_tot)
cases_cum[, cum.sum := cumsum(value), by=list(country_code, simulation_index, strain)]
cases_cum[, value:=NULL]
cases_cum_cis <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- cases_cum[, lapply(.SD, get(meas)), by = c('country_code', 'week', 'strain')]
  dt[, measure := meas]
  cases_cum_cis <- rbind(cases_cum_cis, dt)
}
cases_cum_cis[, simulation_index:=NULL]

cases_cum_cis_wide <- dcast.data.table(cases_cum_cis, 
                                   country_code+week+strain~measure,
                                   value.var = "cum.sum")

ggplot(cases_cum_cis_wide) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, 
                  fill = strain), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                  fill = strain), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000, 
                col = strain), lwd=0.6) +
  scale_color_manual(values = strain_colors1) +
  scale_fill_manual(values = strain_colors1) +
  theme_bw() + facet_grid(country_code~strain, scales='free') +
  xlab('Year') + ylab('Cumulative influenza infections, millions') +
  theme(text=element_text(size=14)) + labs(fill = 'Strain', col='Strain') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) 

cases_fin <- copy(cases_tot)
cases_fin[, strain:=NULL]
cases_fin <- cases_fin[, lapply(.SD, sum, na.rm=T), by=c('country_code','simulation_index','week')]
cases_fin[, cum.sum := cumsum(value), by=list(country_code, simulation_index)]
cases_fin[, value:=NULL]
cases_fin_cis <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- cases_fin[, lapply(.SD, get(meas)), by = c('country_code', 'week')]
  dt[, measure := meas]
  cases_fin_cis <- rbind(cases_fin_cis, dt)
}
cases_fin_cis[, simulation_index:=NULL]

cases_fin_cis_wide <- dcast.data.table(cases_fin_cis, 
                                       country_code+week~measure,
                                       value.var = "cum.sum")

ggplot(cases_fin_cis_wide) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000), 
              fill=2,alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000), 
              fill=2,alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000), col=2,lwd=0.6) +
  theme_bw() + facet_grid(country_code~., scales='free') +
  xlab('Year') + ylab('Cumulative influenza infections, millions') +
  theme(text=element_text(size=14)) + 
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) 

















