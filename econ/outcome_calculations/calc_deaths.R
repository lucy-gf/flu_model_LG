
## HOW MANY DEATHS IN BASE CASE UNDER OUR IFRS?
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
library(wpp2022)

## loading base cases
scenario_name <- 'base'
rm(econ_cases)
for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0('data/vacc_output_', scenario_name,'/econ_inp_', c_code, '.Rdata'))){
    print(paste0(c_code, ' exists'))
    load(paste0('data/vacc_output_', scenario_name,'/econ_inp_', c_code, '.Rdata'))
    if(exists('econ_cases')){
      econ_cases <- rbind(econ_cases, econ_inp)
    }else{
      econ_cases <- copy(econ_inp)
    }
  }
}

econ_cases_c <- melt(econ_cases, 
                     id.vars=c('country_code','simulation_index','year','scenario'))
econ_cases_c[,age_grp:=as.numeric(substr(variable, 3,3))][,variable:=NULL]
econ_cases_agg <- econ_cases_c[, lapply(.SD,sum), 
                               by=c('country_code','simulation_index','year','scenario','age_grp')]
setnames(econ_cases_agg, 'value','infections')

## load ifrs
national_ifrs <- data.table(read_csv('econ/outcome_calculations/data/national_ifrs.csv', 
                                     show_col_types=F))

econ_cases_agg <- econ_cases_agg[national_ifrs, on=c('country_code','simulation_index','age_grp')]

econ_cases_agg[,deaths := ifr*infections]
econ_cases_agg[, ct := substr(scenario,4,4)]
econ_cases_agg[, vt := substr(scenario,9,9)]

# tester <- econ_cases_agg[scenario=='ct_5_vt_1']
# tester[, c('scenario','country_code','attach_code','ct','vt','ifr'):=NULL]
# tester <- tester[,lapply(.SD,sum), by=c('cluster_name','year','simulation_index')]
# tester_mean <- tester[,lapply(.SD,mean), by=c('cluster_name','year')]
# tester_mean[, c('simulation_index','age_grp'):=NULL]
# tester_mean_sums <- tester_mean[,c('cluster_name','year','infections','deaths')][, lapply(.SD,sum), by=c('cluster_name','year')]
# tester_mean_sums_w <- dcast(tester_mean_sums, year ~ cluster_name, value.var = c('deaths','infections'))

## global deaths
g_d <- econ_cases_agg[,c('simulation_index','year','infections','deaths','ct','vt')][, lapply(.SD, sum, na.rm=T),
                                                              by = c('simulation_index','year','ct','vt')]

ggplot(g_d) +
  geom_boxplot(aes(x=year,y=deaths/100000, group=year, fill=vt)) + 
  facet_grid(ct~vt) + ylab('Deaths, 100,000s') +
  scale_fill_manual(values = vt_colors,
                   labels = c('Current','Improved (minimal)','Improved (efficacy)',
                              'Improved (breadth)','Universal')) +
  theme_minimal() 

ggplot(g_d[!vt%in%2:5]) +
  geom_boxplot(aes(x=year,y=deaths/100000, group=year, fill=vt)) + 
  facet_grid(ct~.) + ylab('Deaths, 100,000s') +
  scale_fill_manual(values = vt_colors,
                    labels = c('Current','Improved (minimal)','Improved (efficacy)',
                               'Improved (breadth)','Universal')) +
  theme_bw()

g_d_ann_av <- data.table()
for(meas in c('mean','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- g_d[, lapply(.SD, get(meas)), by = c('year', 'ct','vt')]
  dt[, measure := meas]
  g_d_ann_av <- rbind(g_d_ann_av, dt)
}

g_d_ann_av_plot <- dcast(g_d_ann_av[,c('year','ct','vt','deaths','measure')], year+ct+vt~measure, value.var = 'deaths')

ggplot(g_d_ann_av_plot) +
  geom_line(aes(x=year, y=mean/1000000, group=vt, col=vt), lwd=1) +
  geom_ribbon(aes(x=year, group=vt, ymin=eti50L/1000000, ymax=eti50U/1000000, fill=vt), alpha=0.4) +
  geom_ribbon(aes(x=year, group=vt, ymin=eti95L/1000000, ymax=eti95U/1000000, fill=vt), alpha=0.2) +
  scale_color_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + xlab('') + ylim(c(0,NA)) +
  facet_grid(ct~vt, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Mean annual deaths (millions)') +
  labs(fill = 'Vaccine type', col='Vaccine type') + 
  theme(text=element_text(size=14))
ggsave(paste0("econ/outcome_calculations/plots/base_mean_deaths_annual.png"),
       width=30,height=22,units="cm")

### LIMITING TO THE FIRST TEN YEARS FOR MEANS TO LIMIT INFLUENCE OF POPULATION GROWTH ###
g_d_mean <- g_d[year<2035][,lapply(.SD, mean), by=c('simulation_index','ct','vt')][,year:=NULL]
ggplot(g_d_mean) +
  geom_boxplot(aes(x=vt,y=deaths/100000, group=vt, fill=vt)) + 
  facet_grid(ct~.) + ylab('Deaths, 100,000s') +
  scale_fill_manual(values = vt_colors,
                    labels = c('Current','Improved (minimal)','Improved (efficacy)',
                               'Improved (breadth)','Universal')) +
  theme_minimal() 

g_d_m <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- g_d_mean[, lapply(.SD, get(meas)), by = c('ct','vt')]
  dt[, measure := meas]
  g_d_m <- rbind(g_d_m, dt)
}
g_d_m[, simulation_index:=NULL]
g_d_melt <- melt(g_d_m, id.vars=c('ct','vt','measure'))
g_d_cast <- dcast(g_d_melt, ct + vt + variable ~ measure, value.var=c('value'))
g_d_cast[ct=='v',ct:='None']

g_d_cast$ct <- factor(g_d_cast$ct, levels=c('None',1:5))

ggplot(g_d_cast[variable=='deaths',]) +
  geom_bar(aes(x=vt, y=median/1000000, group=vt, fill=vt), position='dodge',stat='identity') +
  geom_errorbar(aes(x=vt, group=vt, ymin=eti95L/1000000, ymax=eti95U/1000000)) +
  # geom_ribbon(aes(x=vt, ymin=eti50L/100000, ymax=eti50U/100000, fill = as.factor(vt)), alpha=0.4) +
  scale_fill_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + xlab('') +
  facet_grid(.~ct, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  ylab('Annual deaths (millions)') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14))
ggsave(paste0("econ/outcome_calculations/plots/base_mean_deaths.png"),
       width=30,height=22,units="cm")

g_d_cast[ct=='None', ]
g_d_cast[ct==5&vt==1, ]




# g_d_mean[,year:=NULL]
# g_d_nv <- g_d_mean[ct=='v',]
# g_d_mean[g_d_nv, on = c("simulation_index"), nv_deaths := i.deaths]
# g_d_mean[g_d_nv, on = c("simulation_index"), nv_infections := i.infections]
# g_d_mean[,averted_d := nv_deaths-deaths]
# g_d_mean[,averted_i := nv_infections-infections]
# g_d_m <- data.table()
# for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
#   dt <- g_d_mean[, lapply(.SD, get(meas)), by = c('ct','vt')]
#   dt[, measure := meas]
#   g_d_m <- rbind(g_d_m, dt)
# }
# g_d_m[,simulation_index:=NULL]
# g_d_melt <- melt(g_d_m, id.vars=c('ct','vt','measure'))
# g_d_cast <- dcast(g_d_melt, ct + vt + variable ~ measure, value.var=c('value'))
# 
# ggplot(g_d_cast[!ct=='v' & variable=='averted_i',]) +
#   geom_bar(aes(x=vt, y=median/5000000, group=vt, fill=vt), position='dodge',stat='identity') +
#   geom_errorbar(aes(x=vt, group=vt, ymin=eti95L/5000000, ymax=eti95U/5000000)) +
#   # geom_ribbon(aes(x=vt, ymin=eti50L/100000, ymax=eti50U/100000, fill = as.factor(vt)), alpha=0.4) +
#   scale_fill_manual(values = vt_colors, 
#                     labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
#   theme_bw() + xlab('') +
#   facet_grid(.~ct, scales='fixed',
#              labeller = labeller(vt = supp.labs,
#                                  ct = c(supp.labs.age, 'No \nvaccination'))) +
#   ylab('Annual deaths averted') +
#   labs(fill = 'Vaccine type') + 
#   theme(text=element_text(size=14))
























































