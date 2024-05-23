
## CALCULATING HOSPITALISATION RATES
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_colors.R")

library(readr)
library(dplyr)
library(socialmixr)
library(data.table)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(wpp2022)
library(qs)

clusters <- read_csv('data/new_clustering.csv', show_col_types=F)
pop_hist_WPP_data <- data.table(read_csv('econ/outcome_calculations/data/pop_hist_WPP_data.csv', show_col_types=F))
pops_dt <- data.table(country_code = rep(clusters$codes, each=4),
                      age_grp = rep(1:4, 186))
for(i in 1:nrow(pops_dt)){
  country_names <- unlist(unname(clusters[clusters$codes==pops_dt$country_code[i], 
                                          c('country','country_altern','country_altern_2')]))
  if(pops_dt$country_code[i]=='TUR'){country_names <- c(country_names, 'Turkiye')}
  if(pops_dt$age_grp[i] == 1){
    pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:4])))]
  }
  if(pops_dt$age_grp[i] == 2){
    pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 5:7])))]
  }
  if(pops_dt$age_grp[i] == 3){
    pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 8:16])))]
  }
  if(pops_dt$age_grp[i] == 4){
    pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 17:24])))]
  }
  pops_dt[i, tot_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))]
}

n_simulations <- 100

### Exemplars & Brazil
ifr_method <- c('exemplar','whole_itz','brazil')[1]
cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                        "_2010_2019_",n_simulations,".rds"))){
    print(cntr)
    cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                                                          "_2010_2019_",n_simulations,".rds"))[[1]]))
  }
}
ifr_method <- c('exemplar','whole_itz','brazil')[3]
cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_BRA_2010_2019_",
                                                      n_simulations,".rds"))[[1]]))

## turning into national mean age-specific annual number of infections
cases_m <- melt(cases_dt, id.vars=c('country','country_code','simulation_index','week'))
cases_m[, age_grp := as.numeric(substr(variable, 3,3))][, year := year(week)][, c('country','week','variable'):=NULL]
cases_m <- cases_m[year<2016]
cases_m <- cases_m[,lapply(.SD,sum), by=c('country_code','simulation_index','year','age_grp')]
cases_m <- cases_m[,lapply(.SD,mean), by=c('country_code','simulation_index','age_grp')][, year:=NULL]


### Whole ITZs
ifr_method <- c('exemplar','whole_itz','brazil')[2]
cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                        "_2010_2019_",n_simulations,".rds"))){
    print(cntr)
    cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                                                          "_2010_2019_",n_simulations,".rds"))[[1]]))
  }
}

## turning into national mean age-specific annual number of infections
cases_dt15 <- cases_dt[year(week)<2016,]
cases_m2 <- melt(cases_dt15, id.vars=c('country','country_code','simulation_index','week'))
cases_m2[, age_grp := as.numeric(substr(variable, 3,3))][, year := year(week)][, c('country','week','variable'):=NULL]
cases_m2 <- cases_m2[,lapply(.SD,sum), by=c('country_code','simulation_index','year','age_grp')]
cases_m2 <- cases_m2[,lapply(.SD,mean), by=c('country_code','simulation_index','age_grp')][, year:=NULL]

global_cases_ex <- rbind(cases_m[country_code %in% c('AUS','ARG','BRA',
                                                  'CAN','GBR')], cases_m2)

global_cases_ex <- global_cases_ex[pops_dt[country_code %in% unique(global_cases_ex$country_code)], on=c('country_code','age_grp')]

global_cases_ex[, attack := value/age_pop]

ggplot(global_cases_ex) + 
  geom_histogram(aes(x=attack, fill=country_code), bins=100) +
  xlim(c(0,NA)) + theme(legend.position='none') +
  facet_grid(age_grp~.)

global_cases_ex[, attach_code := country_code]
global_cases <- data.table(clusters[,c('codes','cluster_name')])
setnames(global_cases, 'codes','country_code')
global_cases[cluster_name %in% c('Africa', 'Asia-Europe', 'Eastern and Southern Asia'), attach_code := country_code]
global_cases[cluster_name %in% c('Oceania-Melanesia-Polynesia'), attach_code := 'AUS']
global_cases[cluster_name %in% c('Southern America'), attach_code := 'BRA']
global_cases[cluster_name %in% c('Northern America'), attach_code := 'CAN']
global_cases[cluster_name %in% c('Europe'), attach_code := 'GBR']
global_cases[country_code %in% c('ARG'), attach_code := 'ARG']
global_cases <- global_cases[, lapply(.SD, rep, 400)]
global_cases[, age_grp := rep(rep(1:4, each=186), 100)][, simulation_index := rep(1:100, each=4*186)]
global_cases <- global_cases[global_cases_ex[,!'country_code'], on=c('attach_code', 'simulation_index', 'age_grp')]
setnames(global_cases, 'age_pop', 'attached_age_pop')
setnames(global_cases, 'tot_pop', 'attached_tot_pop')

global_cases <- global_cases[pops_dt, on=c('country_code','age_grp')]

global_cases[, infections := attack*age_pop]
global_agg <- global_cases[,c('simulation_index','age_grp','age_pop','infections')][, lapply(.SD, sum), by=c('simulation_index','age_grp')]
global_agg[, attack := infections/age_pop]

ggplot(global_agg) + 
  geom_boxplot(aes(x=age_grp, y=attack, group=age_grp, fill=as.factor(age_grp))) +
  ylim(c(0,0.5)) + theme_bw() + theme(text=element_text(size=14)) +
  ylab('Attack rate') + xlab('Age group') + labs(fill='Age group')

### hospitalisation rates per 100,000 
### from doi: 10.7189/jogh.13.04003

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

hosp_rates <- data.table(
  age_grp = c('<5','65+','all'),
  median = c(224.0,96.8,40.5)/100000,
  l95 = c(118.8, 57.0, 24.3)/100000,
  u95 = c(420.0, 164.3, 67.4)/100000
)

for(i in 1:nrow(hosp_rates)){
  parms <- fcn_fitting(mort_rates = hosp_rates[i,c(2:4)], 
                       mort_probs = c(0.5, 0.025, 0.975))
  hosp_rates[i,"shape"] <- parms[1]
  hosp_rates[i,"rate"] <- parms[2]
  hosp_rates[i,"med_fit"] <- qgamma(p=c(0.5), shape=parms[1], rate=parms[2])
  hosp_rates[i,"l95_fit"] <- qgamma(p=c(0.025), shape=parms[1], rate=parms[2])
  hosp_rates[i,"u95_fit"] <- qgamma(p=c(0.975), shape=parms[1], rate=parms[2])
}

## randomly sampling

hosp_samples <- data.table(
  quantile = runif(100, 0, 1)
)
hosp_samples[, '<5' := qgamma(hosp_samples$quantile, shape=unlist(hosp_rates[age_grp=='<5', 'shape']), rate=unlist(hosp_rates[age_grp=='<5', 'rate']))]
hosp_samples[, '65+' := qgamma(hosp_samples$quantile, shape=unlist(hosp_rates[age_grp=='65+', 'shape']), rate=unlist(hosp_rates[age_grp=='65+', 'rate']))]
hosp_samples[, 'all' := qgamma(hosp_samples$quantile, shape=unlist(hosp_rates[age_grp=='all', 'shape']), rate=unlist(hosp_rates[age_grp=='all', 'rate']))]

## <5:

under5 <- global_agg[age_grp==1]
under5[, hosp := hosp_samples[,'<5']][, ihr := hosp/attack]

## 65+:

over65 <- global_agg[age_grp==4]
over65[, hosp := hosp_samples[,'65+']][, ihr := hosp/attack]

## all pop

all_pop <- global_agg[, lapply(.SD, sum), by=c('simulation_index')]
all_pop[, c('age_grp', 'attack') := NULL]
setnames(all_pop, 'age_pop','tot_pop')
all_pop[, attack := infections/tot_pop]
all_pop[, hosp := hosp_samples[,'all']][, ihr := hosp/attack]

## 5-64
rest_of_pop <- data.table(simulation_index = 1:100)
rest_of_pop[, pop_u5 := under5$age_pop[1]]
rest_of_pop[, pop_o65 := over65$age_pop[1]]
rest_of_pop[, pop_all := all_pop$tot_pop[1]]
rest_of_pop[, pop_rest := pop_all - pop_u5 - pop_o65]
rest_of_pop[, hosp_u5 := under5$hosp]
rest_of_pop[, hosp_o65 := over65$hosp]
rest_of_pop[, hosp_all := all_pop$hosp]
rest_of_pop[, infec_rest := global_agg[age_grp %in% 2:3][, lapply(.SD,sum), by=c('simulation_index')]$infections]

rest_of_pop[, ihr := (pop_all*hosp_all - (pop_u5*hosp_u5) - (pop_o65*hosp_o65))/infec_rest]
rest_of_pop[, hosp := (pop_all*hosp_all - (pop_u5*hosp_u5) - (pop_o65*hosp_o65))/pop_rest]

quantile(100000*under5$ihr, c(0.5,0.025,0.975))
quantile(100000*over65$ihr, c(0.5,0.025,0.975))
quantile(100000*all_pop$ihr, c(0.5,0.025,0.975))
quantile(100000*rest_of_pop$ihr, c(0.5,0.025,0.975))

ggplot(rest_of_pop) + 
  geom_histogram(aes(x=100000*ihr), bins=30) + 
  theme_bw() + geom_vline(xintercept=0, lty=2, col='black') +
  scale_x_continuous(breaks=10*(-12:15)) +
  geom_vline(xintercept=100000*median(rest_of_pop$ihr), lty=2, col='red') +
  xlab('IHR in ages 5-64, per 100,000 infections') + 
  theme(text=element_text(size=14))

ggplot(under5) + 
  geom_histogram(aes(x=100000*ihr), bins=30) + 
  theme_bw() + #geom_vline(xintercept=0, lty=2, col='black') +
  # scale_x_continuous(breaks=10*(-12:15)) +
  geom_vline(xintercept=100000*median(under5$ihr), lty=2, col='red') +
  xlab('IHR in ages <5, per 100,000 infections') + 
  theme(text=element_text(size=14))

ggplot(over65) + 
  geom_histogram(aes(x=100000*ihr), bins=30) + 
  theme_bw() + #geom_vline(xintercept=0, lty=2, col='black') +
  # scale_x_continuous(breaks=10*(-12:15)) +
  geom_vline(xintercept=100000*median(over65$ihr), lty=2, col='red') +
  xlab('IHR in ages 65+, per 100,000 infections') + 
  theme(text=element_text(size=14))

## re-weighting 5-64 into 5-19, 20-64

## cromer data
## from http://dx.doi.org/10.1016/j.jinf.2013.11.013
cromer <- data.table(
  age_grp = c('[5,15)','[15,45)','[45,65)'),
  hosp_rate = c(14, 12, 27)/100000
)

# data from 2005 so reweight by 2005 UK population
GBR_vec <- pop_hist_WPP_data[name == 'United Kingdom' & Year == 2005, ]

cromer_weighted <- data.table(
  age_grp = c('[5,20)','[20,65)'),
  hosp_rate = c((14*sum(GBR_vec[,5:6]) + 12*sum(GBR_vec[,7]))/sum(GBR_vec[,5:7]), 
                (12*sum(GBR_vec[,8:12]) + 27*sum(GBR_vec[,13:16]))/sum(GBR_vec[,8:16]))/100000
)

global_age_dist <- global_agg[simulation_index==1, c('age_grp','age_pop')]

fcn_ihr_weight <- function(ihr, hosp_rates, age_dist){
  
  alpha <- hosp_rates[1]/hosp_rates[2] # ihr1 <- alpha*ihr2
  
  ihr2 <- unname(ihr*sum(age_dist[1:2])/(alpha*age_dist[1] + age_dist[2]))
  
  return(c(alpha*ihr2, ihr2))
}

fcn_ihr_weight(ihr, hosp_rates, age_dist)

## making into x400 data.table
global_ihrs <- rbind(under5[,c('simulation_index','age_grp','ihr')],
                     data.table(simulation_index = rep(1:100, 2),
                          age_grp = c(rep(2,100),rep(3,100)),
                          ihr = NA),
                     over65[,c('simulation_index','age_grp','ihr')])

for(i in 1:100){
  index_ihrs <- fcn_ihr_weight(rest_of_pop$ihr[i], unlist(cromer_weighted$hosp_rate), unname(unlist(global_age_dist[age_grp %in% 2:3,'age_pop'])))
  global_ihrs[age_grp == 2 & simulation_index == i, ihr := index_ihrs[1]]
  global_ihrs[age_grp == 3 & simulation_index == i, ihr := index_ihrs[2]]
}

write_csv(global_ihrs, file='econ/outcome_calculations/data/global_ihrs.csv')

ggplot(global_ihrs) + 
  geom_density(aes(x=ihr*100000, fill=age_grp)) + 
  theme_bw() + facet_grid(age_grp~., scales='free', labeller = labeller(age_grp = supp.labs.agegrps)) + 
  xlab('IHR, per 100,000 infections') + scale_x_log10(limits=c(10,2000),
                                                      breaks=c(10,25,50,100,200,300,
                                                               600,1000,1500,2000)) +
  theme(text=element_text(size=14),
        legend.position='none')  

global_ihrs[, lapply(.SD, eti95L), by=c('age_grp')]$ihr*100000
global_ihrs[, lapply(.SD, median), by=c('age_grp')]$ihr*100000
global_ihrs[, lapply(.SD, eti95U), by=c('age_grp')]$ihr*100000




