
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
library(wpp2022)

## CI functions 

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

## loading base cases
scenario_name <- 'base'
rm(econ_cases)
for(c_code in c('AUS','ARG','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0('data/vacc_output_', scenario_name,'/no-sync.nosync/econ_inp_', c_code, '.Rdata'))){
    print(paste0(c_code, ' exists'))
    load(paste0('data/vacc_output_', scenario_name,'/no-sync.nosync/econ_inp_', c_code, '.Rdata'))
    if(exists('econ_cases')){
      econ_cases <- rbind(econ_cases, econ_inp)
    }else{
      econ_cases <- copy(econ_inp)
    }
  }
}

# save(econ_cases, file='IHME_exemplar_IFRs/data/econ_cases')
# load('IHME_exemplar_IFRs/data/econ_cases')

econ_cases_c <- melt(econ_cases, 
                     id.vars=c('country_code','simulation_index','year','scenario'))
econ_cases_c[,age_grp:=as.numeric(substr(variable, 3,3))][,variable:=NULL]
econ_cases_agg <- econ_cases_c[, lapply(.SD,sum), 
                               by=c('country_code','simulation_index','year','scenario','age_grp')]
setnames(econ_cases_agg, 'value','infections')
setnames(econ_cases_agg, 'country_code','iso3c')

## load ifrs
load('IHME_exemplar_IFRs/data/nat_ifr')

econ_cases_agg[nat_ifr, on=c('iso3c','simulation_index','age_grp'),
               ifr := i.scaled_ifr]
econ_cases_agg[nat_ifr, on=c('iso3c','simulation_index','age_grp'),
               itz := i.itz]
econ_cases_agg[,deaths := ifr*infections]
econ_cases_agg[, ct := substr(scenario,4,4)]
econ_cases_agg[, vt := substr(scenario,9,9)]

econ_cases_agg2 <- copy(econ_cases_agg)[,c('scenario','iso3c'):=NULL]
econ_cases_agg2 <- econ_cases_agg2[,lapply(.SD,sum,na.rm=T),
                                  by = c('simulation_index','year','vt','ct','itz')]
econ_cases_agg2 <- econ_cases_agg2[!is.na(itz)]

econ_cases_agg2 <- econ_cases_agg2[, lapply(.SD,mean), by = c('simulation_index','vt','ct','itz')]

ggplot(econ_cases_agg2[year>2025&ct=='v',]) +
  geom_boxplot(aes(x=itz, y=deaths/100000, group=itz, fill=itz)) +
  # facet_grid(ct~vt) +
  xlab('ITZ') + labs(fill='ITZ') +
  ylab('Deaths, 100,000s') +
  scale_fill_manual(values = exemplar_colors) +
  theme_minimal() 

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

g_d_mean <- g_d[,lapply(.SD, mean), by=c('simulation_index','ct','vt')][,year:=NULL]
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

g_d_cast[ct=='v', ]
g_d_cast[ct==5&vt==1, ]

g_d_mean[,year:=NULL]
g_d_nv <- g_d_mean[ct=='v',]
g_d_mean[g_d_nv, on = c("simulation_index"), nv_deaths := i.deaths]
g_d_mean[g_d_nv, on = c("simulation_index"), nv_infections := i.infections]
g_d_mean[,averted_d := nv_deaths-deaths]
g_d_mean[,averted_i := nv_infections-infections]
g_d_m <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- g_d_mean[, lapply(.SD, get(meas)), by = c('ct','vt')]
  dt[, measure := meas]
  g_d_m <- rbind(g_d_m, dt)
}
g_d_m[,simulation_index:=NULL]
g_d_melt <- melt(g_d_m, id.vars=c('ct','vt','measure'))
g_d_cast <- dcast(g_d_melt, ct + vt + variable ~ measure, value.var=c('value'))

ggplot(g_d_cast[!ct=='v' & variable=='averted_i',]) +
  geom_bar(aes(x=vt, y=median/30000000, group=vt, fill=vt), position='dodge',stat='identity') +
  geom_errorbar(aes(x=vt, group=vt, ymin=eti95L/30000000, ymax=eti95U/30000000)) +
  # geom_ribbon(aes(x=vt, ymin=eti50L/100000, ymax=eti50U/100000, fill = as.factor(vt)), alpha=0.4) +
  scale_fill_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + xlab('') +
  facet_grid(.~ct, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = c(supp.labs.age, 'No \nvaccination'))) +
  ylab('Annual deaths averted') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14))

#### OLD STUFF 

## loading IFR estimates
data_source <- c('GBD','CDC')[2]
load(paste0('IHME_exemplar_IFRs/data/ifrs_', data_source))
ifrs <- copy(get(paste0('ifrs_', data_source)))

no_vacc_cases <- no_vacc_cases[ifrs, on = c('country_code','simulation_index','age_grp')]
no_vacc_cases[,deaths_annual:=ifr*cases_annual]

agg_cases_deaths <- no_vacc_cases[,c('country_code','simulation_index','cases_annual','deaths_annual')]
agg_cases_deaths <- agg_cases_deaths[, lapply(.SD, sum), by=c('country_code','simulation_index')]

ggplot(agg_cases_deaths) +
  geom_boxplot(aes(x=country_code, y=deaths_annual/1000000, group=country_code, fill=country_code)) +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  labs(fill = 'ITZ') + xlab('') + ylab('Annual deaths, millions') +
  theme(text=element_text(size=14)) + ggtitle(paste0(data_source,' death projections'))
ggsave(paste0("IHME_exemplar_IFRs/plots/", data_source, "_deaths_dist.png"),
       width=26,height=24,units="cm")
ggplot(agg_cases_deaths) +
  geom_boxplot(aes(x=country_code, y=cases_annual/1000000, group=country_code, fill=country_code)) +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  labs(fill = 'ITZ') + xlab('') + ylab('Annual cases, millions') +
  theme(text=element_text(size=14))
ggsave(paste0("IHME_exemplar_IFRs/plots/", data_source, "_cases_dist.png"),
       width=26,height=24,units="cm")
ggplot(agg_cases_deaths) +
  geom_bar(aes(x=simulation_index, y=deaths_annual/1000000, fill=country_code),
           position='fill',stat='identity') +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  labs(fill = 'ITZ') + xlab('') + ylab('Proportion of deaths') +
  theme(text=element_text(size=14))

# pop_proj_WPP_data <- data.table(read_csv('data_for_BS/pop_proj_WPP_data.csv'))[Year==2040]
# for(name_i in unique(pop_proj_WPP_data$name)){
#   for(j in 1:nrow(clusters)){
#     if(name_i %in% clusters[j, c('country', 'country_altern', 'country_altern_2')]){
#       pop_proj_WPP_data[name == name_i, iso3c := clusters[j, codes]]
#     }
#   }
# }
# pops <- pop_proj_WPP_data[!is.na(iso3c),][,c('name','iso3c')]
# for(i in 1:nrow(pops)){
#   pops[i, pop := 1000000*sum(pop_proj_WPP_data[i,4:24])]
# }
# pops <- pops[clusters[, c('iso3c','cluster_name')], on='iso3c']
# itz_pops <- pops[,c('cluster_name','pop')][,lapply(.SD, sum), by='cluster_name']
# ggplot(pops) +
#   geom_bar(aes(x=cluster_name, y=pop/1000000000, group=cluster_name, fill=cluster_name),
#            position='stack',stat='identity') +
#   theme_bw() + scale_fill_manual(values=cluster_colors2) +
#   labs(fill = 'ITZ') + xlab('') + ylab('2040 population, billions') +
#   theme(text=element_text(size=14))

global_nv_cases <- no_vacc_cases[,c('simulation_index','cases_annual','deaths_annual')][
  ,lapply(.SD,sum),by=c('simulation_index')
]

ggplot(global_nv_cases) +
  geom_histogram(aes(x=deaths_annual/100000), bins=30) + theme_minimal() +
  xlab('Annual deaths, 100,000s')
quantile(global_nv_cases$deaths_annual, probs=c(0.025,0.5,0.975))


## what if we use non-age-specific IFRs?
data_source <- c('GBD','CDC')[1]
agg_ifrs <- get(paste0('agg_',data_source))
load('IHME_exemplar_IFRs/data/no_vacc_cases')
no_vacc_cases <- no_vacc_cases[agg_ifrs[,c('country_code','simulation_index','GBD_ifr')],
                               on = c('country_code','simulation_index')]
no_vacc_cases[, deaths_annual := cases_annual*GBD_ifr]
over_ages <- no_vacc_cases[,lapply(.SD,sum), by=c('country_code','simulation_index')][,c('age_grp','cases5yr'):=NULL]
global_agg <- over_ages[,c('simulation_index','cases_annual','deaths_annual')][, lapply(.SD,sum), by=c('simulation_index')]
quantile(global_agg$deaths_annual, probs=c(0.025,0.5,0.975))

data_source <- c('GBD','CDC')[2]
agg_ifrs <- get(paste0('agg_',data_source))
load('IHME_exemplar_IFRs/data/no_vacc_cases')
no_vacc_cases[, age_grp2 := ifelse((age_grp%in%1:3), '<65', '65+')]
no_vacc_cases <- no_vacc_cases[agg_ifrs[,c('country_code','simulation_index','CDC_ifr','age_grp2')],
                               on = c('country_code','simulation_index','age_grp2')]
no_vacc_cases[, deaths_annual := cases_annual*CDC_ifr]
over_ages <- no_vacc_cases[,c('country_code','simulation_index','cases_annual','deaths_annual')][,lapply(.SD,sum), by=c('country_code','simulation_index')][,c('age_grp','cases5yr'):=NULL]
global_agg2 <- over_ages[,c('simulation_index','cases_annual','deaths_annual')][, lapply(.SD,sum), by=c('simulation_index')]
quantile(global_agg2$deaths_annual, probs=c(0.025,0.5,0.975))


## what if we just multiply the CDC/GBD mortality rates by 2040ish pop?
# using data from maps.R

clusters <- data.table(read_csv('data_for_BS/new_clustering.csv'))
clusters[, iso3c:=codes]
# pop_proj_WPP_data <- data.table(read_csv('data_for_BS/pop_proj_WPP_data.csv'))[Year==2040]
pop_proj_WPP_data <- data.table(read_csv('data_for_BS/pop_hist_WPP_data.csv'))[Year==2015]
for(name_i in unique(pop_proj_WPP_data$name)){
  for(j in 1:nrow(clusters)){
    if(name_i %in% clusters[j, c('country', 'country_altern', 'country_altern_2')]){
      pop_proj_WPP_data[name == name_i, iso3c := clusters[j, codes]]
    }
  }
}
pops <- pop_proj_WPP_data[,c('name','iso3c')]
for(i in 1:nrow(pops)){
  if(!is.na(pops[i,'iso3c'])){
    pops[i, pop := 1000*sum(pop_proj_WPP_data[i,4:24])]
    pops[i, popu65 := 1000*sum(pop_proj_WPP_data[i,4:16])]
    pops[i, pop6574 := 1000*sum(pop_proj_WPP_data[i,17:18])]
    pops[i, popo75 := 1000*sum(pop_proj_WPP_data[i,19:24])]
  }
}
pops <- pops[!is.na(iso3c),]
pops <- pops[clusters[, c('iso3c','cluster_name')], on='iso3c']
pops <- pops[ihme_data, on = 'iso3c', ihme_mort_rate := i.ihme_mort_rate/100000]
pops <- pops[cdc_data, on = 'iso3c', cdc_mort_u65_med := i.cdc_mort_u65_med/100000]
pops <- pops[cdc_data, on = 'iso3c', cdc_mort_6574_med := i.cdc_mort_6574_med/100000]
pops <- pops[cdc_data, on = 'iso3c', cdc_mort_o75_med := i.cdc_mort_o75_med/100000]
pops[, gbd_deaths := pop*ihme_mort_rate]
pops[, cdc_deaths := popu65*cdc_mort_u65_med + pop6574*cdc_mort_6574_med + popo75*cdc_mort_o75_med]

itz_indiv_deaths <- melt(pops[, c('cluster_name','cdc_deaths','gbd_deaths')][, lapply(.SD, sum, na.rm=T), by=c('cluster_name')], id.vars='cluster_name')
ggplot(itz_indiv_deaths) +
  geom_bar(aes(x=cluster_name, y=value/100000, group=variable, fill=cluster_name, lty=variable),
           col=1, position='dodge',stat='identity') +
  theme_bw() + scale_fill_manual(values = cluster_colors2) +
  ylab('Deaths, 100,000s') + labs(fill='ITZ',lty='Data source') +
  xlab('') + theme(text=element_text(size=14),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

global_deaths <- itz_indiv_deaths[,c('variable','value')][, lapply(.SD, sum), by = 'variable']
global_deaths[,value_huntho := value/100000]

## how do these mort rates compare to the exemplars?
pops[cluster_name=='Asia-Europe', ihme_ex := pops[iso3c=='TUR',ihme_mort_rate]]
pops[cluster_name=='Africa', ihme_ex := pops[iso3c=='GHA',ihme_mort_rate]]
pops[cluster_name=='Europe', ihme_ex := pops[iso3c=='GBR',ihme_mort_rate]]
pops[cluster_name=='Southern America', ihme_ex := pops[iso3c=='ARG',ihme_mort_rate]]
pops[cluster_name=='Eastern and Southern Asia', ihme_ex := pops[iso3c=='CHN',ihme_mort_rate]]
pops[cluster_name=='Northern America', ihme_ex := pops[iso3c=='CAN',ihme_mort_rate]]
pops[cluster_name=='Oceania-Melanesia-Polynesia', ihme_ex := pops[iso3c=='AUS',ihme_mort_rate]]

pops[cluster_name=='Asia-Europe', cdc_u65_ex := pops[iso3c=='TUR',cdc_mort_u65_med]]
pops[cluster_name=='Africa', cdc_u65_ex := pops[iso3c=='GHA',cdc_mort_u65_med]]
pops[cluster_name=='Europe', cdc_u65_ex := pops[iso3c=='GBR',cdc_mort_u65_med]]
pops[cluster_name=='Southern America', cdc_u65_ex := pops[iso3c=='ARG',cdc_mort_u65_med]]
pops[cluster_name=='Eastern and Southern Asia', cdc_u65_ex := pops[iso3c=='CHN',cdc_mort_u65_med]]
pops[cluster_name=='Northern America', cdc_u65_ex := pops[iso3c=='CAN',cdc_mort_u65_med]]
pops[cluster_name=='Oceania-Melanesia-Polynesia', cdc_u65_ex := pops[iso3c=='AUS',cdc_mort_u65_med]]

pops[cluster_name=='Asia-Europe', cdc_6574_ex := pops[iso3c=='TUR',cdc_mort_6574_med]]
pops[cluster_name=='Africa', cdc_6574_ex := pops[iso3c=='GHA',cdc_mort_6574_med]]
pops[cluster_name=='Europe', cdc_6574_ex := pops[iso3c=='GBR',cdc_mort_6574_med]]
pops[cluster_name=='Southern America', cdc_6574_ex := pops[iso3c=='ARG',cdc_mort_6574_med]]
pops[cluster_name=='Eastern and Southern Asia', cdc_6574_ex := pops[iso3c=='CHN',cdc_mort_6574_med]]
pops[cluster_name=='Northern America', cdc_6574_ex := pops[iso3c=='CAN',cdc_mort_6574_med]]
pops[cluster_name=='Oceania-Melanesia-Polynesia', cdc_6574_ex := pops[iso3c=='AUS',cdc_mort_6574_med]]

pops[cluster_name=='Asia-Europe', cdc_o75_ex := pops[iso3c=='TUR',cdc_mort_o75_med]]
pops[cluster_name=='Africa', cdc_o75_ex := pops[iso3c=='GHA',cdc_mort_o75_med]]
pops[cluster_name=='Europe', cdc_o75_ex := pops[iso3c=='GBR',cdc_mort_o75_med]]
pops[cluster_name=='Southern America', cdc_o75_ex := pops[iso3c=='ARG',cdc_mort_o75_med]]
pops[cluster_name=='Eastern and Southern Asia', cdc_o75_ex := pops[iso3c=='CHN',cdc_mort_o75_med]]
pops[cluster_name=='Northern America', cdc_o75_ex := pops[iso3c=='CAN',cdc_mort_o75_med]]
pops[cluster_name=='Oceania-Melanesia-Polynesia', cdc_o75_ex := pops[iso3c=='AUS',cdc_mort_o75_med]]

ggplot(pops) + 
  geom_point(aes(x=ihme_mort_rate, y=ihme_ex, col=cluster_name)) + 
  geom_line(aes(x=ihme_mort_rate, y=ihme_mort_rate), lty=2) +
  theme_bw() + scale_color_manual(values=cluster_colors2)

pops[, gbd_ex_deaths := pop*ihme_ex]
pops[, cdc_ex_deaths := popu65*cdc_u65_ex + pop6574*cdc_6574_ex + popo75*cdc_o75_ex]

ggplot(pops) + 
  geom_point(aes(x=cdc_deaths, y=cdc_ex_deaths, col=cluster_name)) + 
  geom_line(aes(x=cdc_deaths, y=cdc_deaths), lty=2) +
  theme_bw() + scale_color_manual(values=cluster_colors2) +
  scale_y_log10() + scale_x_log10()

itz_ex_deaths <- melt(pops[, c('cluster_name','cdc_ex_deaths','gbd_ex_deaths')][, lapply(.SD, sum, na.rm=T), by=c('cluster_name')], id.vars='cluster_name')
ggplot(itz_ex_deaths) +
  geom_bar(aes(x=cluster_name, y=value/100000, group=variable, fill=cluster_name, lty=variable),
           col=1, position='dodge',stat='identity') +
  theme_bw() + scale_fill_manual(values = cluster_colors2) +
  ylab('Deaths, 100,000s') + labs(fill='ITZ',lty='Data source') +
  xlab('') + theme(text=element_text(size=14),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

global_ex_deaths <- itz_ex_deaths[,c('variable','value')][, lapply(.SD, sum), by = 'variable']
global_ex_deaths[,value_huntho := value/100000]

cbind(itz_indiv_deaths,itz_ex_deaths[,2:3])
## over-estimating in GBR and ARG, underestimating in TUR (and CAN a bit)

## what are annual cases/attack rate in exemplars in 2010:2019 under current vaccinations
## and 2025:2054 under no vaccinations?

## 2010:2019: (from exemplar_plots.R)
source('IHME_exemplar_IFRs/exemplar_plots.R')
mean_inference_cases <- cases_annual[, lapply(.SD,mean), by=c('country_code')]

nv_cases <- data.table()
for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  nv_cases <- rbind(nv_cases,
                    (data.table(readRDS(paste0("data/vacc_output_base/vacc_", c_code, '_none_ct_1.rds'))[[1]])[country_code == c_code,]))
  print(c_code)
}
nv_cases[, year:=year(week)]
nv_cases[,c('country','week'):=NULL]
nv_cases <- nv_cases[, lapply(.SD, sum), by=c('country_code', 'simulation_index','year')]
nv_cases_l <- melt(nv_cases, id.vars=c('country_code','simulation_index','year'))
nv_cases_l[, variable:=NULL]
no_vacc_cases <- nv_cases_l[, lapply(.SD, sum), by=c('country_code','simulation_index','year')]
mean_nv_cases <- no_vacc_cases[, lapply(.SD,mean), by=c('country_code')]

mean_nv_cases[mean_inference_cases, on=c('country_code'), mean_inference_cases := i.value]

ggplot(mean_nv_cases) +
  geom_point(aes(x=mean_inference_cases, y=value, col=country_code)) +
  geom_line(aes(x=mean_inference_cases, y=mean_inference_cases),lty=2) +
  theme_bw() + scale_y_log10() + scale_x_log10()

mean_nv_cases[, mult := value/mean_inference_cases]

to_show <- mean_nv_cases[,c('country_code','value','mean_inference_cases','mult')]
setnames(to_show, 'mean_inference_cases','mean_inference_cases_2010_2019')
setnames(to_show, 'value','mean_predic_cases_2025_2054')
setnames(to_show, 'mult','multiplier')


## OK what about 70% 0-18 & 65+ cases?

cov_cases <- data.table()
for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  load(paste0("data/vacc_output_base/global_weekly_", c_code, '.Rdata'))
  cov_cases <- rbind(cov_cases,cbind(global_weekly_c[scenario=='ct_5_vt_1',],c_code))
  print(c_code)
}

cov_cases[, year:=year(week)]
cov_cases[, c('scenario','week'):=NULL]
cov_cases <- cov_cases[,lapply(.SD,sum), by=c('c_code','year','simulation_index')]
cov_cases_m <- melt(cov_cases, id.vars = c('c_code','year','simulation_index'))
cov_cases_m[, variable:=NULL]
cov_cases_m <- cov_cases_m[,lapply(.SD,sum), by=c('c_code','year','simulation_index')]
cov_cases_m <- cov_cases_m[year>2025,]

ggplot(cov_cases_m) +
  geom_boxplot(aes(x=c_code, y=value/1000000, group=c_code, fill=c_code)) +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  labs(fill = 'ITZ') + xlab('') + ylab('Annual cases, millions') +
  theme(text=element_text(size=14))

cov_cases_med <- cov_cases_m[, lapply(.SD,median), by=c('c_code')]
sum(cov_cases_med$value)/1000000000 # median global annual cases, billions

cov_cases_m_as <- melt(cov_cases, id.vars = c('c_code','year','simulation_index'))
cov_cases_m_as[, age_grp := as.numeric(substr(variable, 3,3))]
cov_cases_m_as[, variable:=NULL]
cov_cases_m_as <- cov_cases_m_as[,lapply(.SD,sum), by=c('c_code','year','simulation_index','age_grp')]
cov_cases_m_as <- cov_cases_m_as[year>2025,]
cov_cases_m_as <- cov_cases_m_as[,lapply(.SD,mean),by=c('c_code','simulation_index','age_grp')]
setnames(cov_cases_m_as, 'c_code','country_code')
setnames(cov_cases_m_as, 'value','cases')
cov_cases_m_as[,year:=NULL]

## adding ifrs:
load('IHME_exemplar_IFRs/data/ifrs_CDC')
load('IHME_exemplar_IFRs/data/ifrs_GBD')

cov_cases_m_as <- cov_cases_m_as[ifrs_CDC, on=c('country_code','simulation_index','age_grp'),
                                 cdc_ifr := i.ifr]
cov_cases_m_as[, cdc_deaths := cases*cdc_ifr]

ggplot(cov_cases_m_as) +
  geom_boxplot(aes(x=country_code, y=cases/1000000, group=country_code, fill=country_code)) +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  facet_grid(age_grp~., scales='free', labeller=labeller(age_grp=supp.labs.agegrps)) + 
  labs(fill = 'ITZ') + xlab('') + ylab('Annual cases, millions') +
  theme(text=element_text(size=14))

ggplot(cov_cases_m_as[,lapply(.SD,sum),by=c('country_code','simulation_index')]) +
  geom_boxplot(aes(x=country_code, y=cases/1000000, group=country_code, fill=country_code)) +
  geom_point(data=mean_inference_cases, aes(x=country_code,y=value, col=country_code)) +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  scale_color_manual(values=exemplar_colors) +
  facet_grid(age_grp~., scales='free', labeller=labeller(age_grp=supp.labs.agegrps)) + 
  labs(fill = 'ITZ') + xlab('') + ylab('Annual cases, millions') +
  theme(text=element_text(size=14)) 

ggplot(cov_cases_m_as) +
  geom_boxplot(aes(x=country_code, y=cdc_deaths, group=country_code, fill=country_code)) +
  theme_bw() + scale_fill_manual(values=exemplar_colors) +
  facet_grid(age_grp~., scales='free', labeller=labeller(age_grp=supp.labs.agegrps)) + 
  labs(fill = 'ITZ') + xlab('') + ylab('Annual deaths') +
  theme(text=element_text(size=14))

ggplot(cov_cases_m_as) +
  geom_boxplot(aes(x=country_code, y=cdc_deaths, group=interaction(age_grp, country_code), 
                   fill=as.factor(age_grp))) +
  theme_bw() + scale_fill_manual(values=age_colors1,labels=supp.labs.agegrps) +
  # facet_grid(age_grp~., scales='free', labeller=labeller(age_grp=supp.labs.agegrps)) + 
  labs(fill = 'Age group') + xlab('ITZ') + ylab('Annual deaths') +
  theme(text=element_text(size=14))

## global deaths
global_as_d <- copy(cov_cases_m_as)
global_as_d[,c('country_code','age_grp','cdc_ifr'):=NULL]
global_as_d <- global_as_d[,lapply(.SD,sum), by=c('simulation_index')]

quantile(global_as_d$cdc_deaths, c(0.025,0.5,0.975))

ggplot(melt(global_as_d, id.vars='simulation_index')) +
  geom_histogram(aes(x=value/1000000), bins=50) +
  theme_bw() + 
  facet_wrap(variable~., scales='free', nrow=2)

# per person?
pop_proj_WPP_data40 <- data.table(read_csv('data_for_BS/pop_proj_WPP_data.csv'))[Year==2040]
pop_proj_WPP_data15 <- data.table(read_csv('data_for_BS/pop_hist_WPP_data.csv'))[Year==2015]
pop40 <- sum(pop_proj_WPP_data40[,4:24])*1000
pop15 <- sum(pop_proj_WPP_data15[,4:24])*1000

cov_cases_m_as <- melt(cov_cases, id.vars = c('c_code','year','simulation_index'))
cov_cases_m_as[, age_grp := as.numeric(substr(variable, 3,3))]
cov_cases_m_as[, variable:=NULL]
cov_cases_m_as <- cov_cases_m_as[,lapply(.SD,sum), by=c('c_code','year','simulation_index','age_grp')]
cov_cases_m_as <- cov_cases_m_as[year>2025,]
setnames(cov_cases_m_as, 'c_code','country_code')
setnames(cov_cases_m_as, 'value','cases')
cov_cases_m_as <- cov_cases_m_as[ifrs_CDC, on=c('country_code','simulation_index','age_grp'),
                                 cdc_ifr := i.ifr]
cov_cases_m_as[, cdc_deaths := cases*cdc_ifr]

global_as_d <- copy(cov_cases_m_as)
global_as_d[,c('country_code','age_grp','cdc_ifr'):=NULL]
global_as_d <- global_as_d[,lapply(.SD,sum), by=c('simulation_index','year')]

for(year_i in unique(global_as_d$year)){
  global_as_d[year==year_i, global_pop:= sum(data.table(read_csv('data_for_BS/pop_proj_WPP_data.csv'))[Year==year_i,4:24])*1000]
}

global_as_d[,annual_d_pp := cdc_deaths/global_pop]
c(300000, 650000)*pop40/pop15

ggplot(global_as_d) +
  geom_histogram(aes(x=annual_d_pp),bins=50) +
  geom_vline(xintercept=300000/pop15, lty=2) +
  geom_vline(xintercept=650000/pop15, lty=2) +
  theme_bw()

pop40/pop15




# ## is my population too big..?
# source('BS/BS_data_fcns.R')
# 
# global_pop <- data.table(country=clusters$country)
# for(i in 1:nrow(clusters)){
#   global_pop[i, iso3c:=clusters[i,'codes']]
#   demog <- data.table(fcn_weekly_demog(country = clusters[i,2:4],
#                                        pop_coverage = c(0,0,0,0),
#                                        weeks_vaccinating = 12,
#                                        first_year_all = T,
#                                        imm_duration = 1, # in years 
#                                        coverage_pattern = 1,
#                                        hemisphere = 'NH',
#                                        start_year = 2025, years = 30))
#   global_pop[i, pop54 := sum(demog[week=="2054-12-28",value])]
#   print(paste0(unname(unlist(clusters[i,'codes'])), i))
# }
# 
# wpp_proj <- data.table(read_csv('data_for_BS/pop_proj_WPP_data.csv'))[Year==2054,]
# 
# for(i in 1:nrow(clusters)){
#   global_pop[i, real_pop54:= sum(wpp_proj[name %in% clusters[i,2:4],4:24])*1000]
# }
# 
# ggplot(global_pop) +
#   geom_point(aes(x=real_pop54/1000000, y=pop54/1000000)) +
#   geom_line(aes(x=real_pop54/1000000, y=real_pop54/1000000), lty=2) +
#   theme_bw() + scale_x_log10() + scale_y_log10()
# 
# sum(global_pop$pop54)/1000000000
# sum(global_pop$real_pop54)/1000000000
# 
# clusters <- data.table(clusters)
# clusters[,iso3c:=codes]
# global_pop <- global_pop[clusters, on=c('iso3c'), itz := i.cluster_name]
# itz_pop <- global_pop[,c('itz','pop54')][,lapply(.SD,sum),by=c('itz')]




















