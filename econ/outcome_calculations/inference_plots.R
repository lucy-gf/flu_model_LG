
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

n_simulations <- 100

### first: exemplar countries, and brazil
ifr_method <- c('exemplar','whole_itz','brazil')[1]
## load data
cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                        "_2010_2019_",n_simulations,".rds"))){
    print(cntr)
  cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/exemplar_", cntr, 
                                             "_2010_2019_",n_simulations,".rds"))[[1]]))
  }
}
cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/brazil_BRA_2010_2019_",n_simulations,".rds"))[[1]]))
## strain-specific cases
cases_dt[,country:=NULL]
cases_dt_m <- data.table(melt(cases_dt, id.vars = c("country_code", "simulation_index", "week")))
cases_dt_m[grepl('A', variable), strain := "A"]
cases_dt_m[grepl('B', variable), strain := "B"]
cases_tot <- copy(cases_dt_m)
cases_tot[,variable:=NULL]
cases_tot <- cases_tot[, lapply(.SD, sum, na.rm=T), by=c('country_code','simulation_index','week','strain')]
## with CIs
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
## cumulative strain-specific cases
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
## cumulative aggregated cases
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
## total number of cases in 2010-2015
cases_end <- cases_fin[cases_fin$week==max(cases_fin$week[year(cases_fin$week)<2016]),]
cases_end[,week:=NULL]
## attack rates
ages <- data.table(
  country_code = rep(c('ARG','AUS','CAN','CHN','GBR','GHA','TUR','BRA'), each=4),
  age_grp = rep(1:4,8), 
  age_pop = c(pop_age(wpp_age('Argentina', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('Australia', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('Canada', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('China', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('United Kingdom', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('Ghana', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('Turkey', 2015), age.limits=c(0,5,20,65))$population,
              pop_age(wpp_age('Brazil', 2015), age.limits=c(0,5,20,65))$population)) 
ages_agg <- ages[, lapply(.SD, sum), by=c('country_code')]
cases_as <- cases_end[ages_agg, on=c('country_code')]

### PLOTS
ggplot(cases_tot) +
  geom_line(aes(x=week, y=value/1000000, group=simulation_index, col=strain), alpha=0.2) +
  theme_minimal() + ylab('Cases, millions') + 
  scale_color_manual(values = strain_colors1) +
  facet_grid(country_code~strain, scales='free')

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
ggsave(paste0("econ/outcome_calculations/plots/exemplar_cases.png"),
       width=26,height=24,units="cm")

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
ggsave(paste0("econ/outcome_calculations/plots/exemplar_cumcases.png"),
       width=26,height=24,units="cm")

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
ggsave(paste0("econ/outcome_calculations/plots/exemplar_cumcasesagg.png"),
       width=26,height=24,units="cm")

ggplot(cases_end) +
  geom_histogram(aes(x=cum.sum/1000000, fill=country_code), bins=30) +
  scale_fill_manual(values = exemplar_colors) +
  theme_minimal() + xlab('Cumulative influenza cases, millions') +
  facet_grid(.~country_code, scales='free') +
  labs(fill='Country') 

ggplot(cases_as) +
  # geom_histogram(aes(x=cum.sum/(age_pop*10), fill=country_code), bins=30) +
  geom_density(aes(x=cum.sum/(age_pop*5), fill=country_code)) +
  scale_fill_manual(values = exemplar_colors) +
  theme_bw() + xlab('Annual attack rate') +
  # facet_grid(.~country_code, scales='free') +
  labs(fill='Country') #+ xlim(c(0.1,0.45))
ggsave(paste0("econ/outcome_calculations/plots/exemplar_attack_rates.png"),
       width=26,height=24,units="cm")

### okay much too high attack rate in brazil

## looking at whole itzs
ifr_method <- c('exemplar','whole_itz','brazil')[2]
## load data
cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                        "_2010_2019_",n_simulations,".rds"))){
    print(cntr)
    cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                                                          "_2010_2019_",n_simulations,".rds"))[[1]]))
  }
}

unique(cases_dt$country_code)
cases_dt[,country:=NULL]
cases_box <- melt(cases_dt, id.vars=c('country_code','simulation_index','week'))
cases_box[, year:=year(week)][,week:=NULL]
cases_box <- cases_box[,variable := NULL][,lapply(.SD,sum), by=c(c('country_code','simulation_index','year'))]
cases_box15 <- cases_box[year<2016,]
cases_box15 <- cases_box15[, lapply(.SD,mean), by=c('country_code','simulation_index')]
cases_box15[clusters, on=c('country_code'), itz := i.cluster_code]
pop_proj_WPP_data <- data.table(read_csv('data_for_BS/pop_hist_WPP_data.csv'))[Year==2015]
for(name_i in unique(pop_proj_WPP_data$name)){
  for(j in 1:nrow(clusters)){
    if(name_i %in% clusters[j, c('country', 'country_altern', 'country_altern_2')]){
      pop_proj_WPP_data[name == name_i, country_code := clusters[j, codes]]
    }
  }
}
pop_proj_WPP_data$tot_pop <- rowSums(pop_proj_WPP_data[,4:24])*1000
cases_box15[pop_proj_WPP_data, on=c('country_code'), tot_pop := i.tot_pop]
ggplot(cases_box15) + 
  geom_boxplot(aes(x=country_code, y=value/tot_pop, fill=itz)) +
  theme_minimal() + ylim(c(0,max(cases_box15$ar)))












