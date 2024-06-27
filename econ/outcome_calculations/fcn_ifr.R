
## FUNCTION TO CALCULATE AGE-SPECIFIC IFRS
## FROM MORTALITY + INCIDENCE DATA (2010-2015)
#setwd("~/Desktop/research asst/Global Code")

library(readr)
library(dplyr)
library(data.table)
library(countrycode)
library(readxl)
library(tidyverse)
library(wpp2022)
library(fitdistrplus)

f.beta <- function(alpha, beta, x) {
  p <- pbeta(x, alpha, beta)
  # return both
  return(c(p))
}
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
               control = list(reltol = 1e-20)
  ))
  parms <- (sol$par)       

  return(parms)
}

clusters <- data.table(read_csv('econ/outcome_calculations/data/new_clustering.csv', show_col_types=F))[, iso3c := codes]
data(mx1dt)
mx1dt$iso3c <- suppressWarnings(countrycode(mx1dt$country_code, origin = 'un', destination = 'iso3c'))
pop_hist_WPP_data <- data.table(read_csv('econ/outcome_calculations/data/pop_hist_WPP_data.csv', show_col_types=F))

age_transform <- function(vec){
  return(c(rep((vec[1:(length(vec)-1)]/5), each=5), vec[length(vec)]))
}

## fcn_ifr - make all data into age-specific ifrs

fcn_ifr <- function(
  country_code, ## iso3c code for country, used to obtain age-specific population sizes and  UNLT
  mortality_rates, ## mortality rates (median, 95% CI) - will be a data.table with 3 age categories
  mortality_probs = c(0.5, 0.025, 0.975), ## corresponding quantile probabilities
  incidence, ## age-specific influenza incidence rates in 2010-2015, 100 simulations
  incidence_ages = c(0,5,20,65) ## lower bound of each age group for the incidence counts
){
  
  ## make mortality_rates into a data.table if it isn't already
  mortality_rates <- data.table(mortality_rates)
  ## check that they are reasonably a probability, not rate per 100,000
  if(sum(mortality_rates[,!'age_grp']>1)>0){
    stop('Mortality rates must be individual-level, not per 100,000')
  }
  ## check mortality rates are in correct order (not very thorough)
  if(((mortality_probs[1] - mortality_probs[2] > 0) &
     (mortality_probs[3] - mortality_probs[1] > 0))){
    if(((sum((mortality_rates[,2] - mortality_rates[,3])<0)>0 |
         sum((mortality_rates[,4] - mortality_rates[,2])<0)>0))){
      stop('Mortality rates must be individual-level, not per 100,000')
    }
  }
  
  ## make incidence a data.table if not already
  incidence <- data.table(incidence)
  
  ## 1: make mortality rates interval into continuous distributions
  for(i in 1:nrow(mortality_rates)){
    parms <- fcn_fitting(mort_rates = mortality_rates[i,c(2:4)], 
                      mort_probs = mortality_probs[c(1:3)])
    mortality_rates[i,"shape"] <- parms[1]
    mortality_rates[i,"rate"] <- parms[2]
  }
  
  ## 2: get 100 mortality rate samples for each age group
  mortality_samples <- data.table(age_grp = rep(c('<65','65-75','75+'), each=100),
                                  simulation_index = rep(1:100,3))
  for(age_grp_i in unique(mortality_samples$age_grp)){
    # mortality_samples[age_grp==age_grp_i,
    #               value:= c(rbeta(n = 100, shape1 = as.numeric(mortality_rates[age_grp==age_grp_i, 'shape']),
    #                               shape2 = as.numeric(mortality_rates[age_grp==age_grp_i, 'rate'])))]
    mortality_samples[age_grp==age_grp_i,
                      value:= c(rgamma(n = 100, shape = as.numeric(mortality_rates[age_grp==age_grp_i, 'shape']),
                                       rate = as.numeric(mortality_rates[age_grp==age_grp_i, 'rate'])))]
  }
  
  # check quantiles of samples 
  # ignoring lower quantile as very hard to fit
  if(100*abs(median(mortality_samples[age_grp=='<65']$value) - mortality_rates[age_grp=='<65','med'])/mortality_rates[age_grp=='<65','med'] > 15){
    print(paste0('Median of sample more than 15% out (<65): ', 
                 round(100*abs(median(mortality_samples[age_grp=='<65']$value) - mortality_rates[age_grp=='<65','med'])/mortality_rates[age_grp=='<65','med'], digits=2), '%'))
  }
  if(100*abs(median(mortality_samples[age_grp=='65-75']$value) - mortality_rates[age_grp=='65-75','med'])/mortality_rates[age_grp=='65-75','med'] > 15){
    print(paste0('Median of sample more than 15% out (65-75): ',
                 round(100*abs(median(mortality_samples[age_grp=='65-75']$value) - mortality_rates[age_grp=='65-75','med'])/mortality_rates[age_grp=='65-75','med'], digits=2), '%'))
  }
  if(100*abs(median(mortality_samples[age_grp=='75+']$value) - mortality_rates[age_grp=='75+','med'])/mortality_rates[age_grp=='75+','med'] > 15){
    print(paste0('Median of sample more than 15% out (75+): ',
                 round(100*abs(median(mortality_samples[age_grp=='75+']$value) - mortality_rates[age_grp=='75+','med'])/mortality_rates[age_grp=='75+','med'], digits=2), '%'))
  }
  if(100*abs(quantile(mortality_samples[age_grp=='<65']$value, probs=0.975) - mortality_rates[age_grp=='<65','u95'])/mortality_rates[age_grp=='<65','u95'] > 25){
    print(paste0('U95 of sample more than 25% out (<65): ',
                 round(100*abs(quantile(mortality_samples[age_grp=='<65']$value, probs=0.975) - mortality_rates[age_grp=='<65','u95'])/mortality_rates[age_grp=='<65','u95'], digits=2), '%'))
  }
  if(100*abs(quantile(mortality_samples[age_grp=='65-75']$value, probs=0.975) - mortality_rates[age_grp=='65-75','u95'])/mortality_rates[age_grp=='65-75','u95'] > 25){
    print(paste0('U95 of sample more than 25% out (65-75): ',
                 round(100*abs(quantile(mortality_samples[age_grp=='65-75']$value, probs=0.975) - mortality_rates[age_grp=='65-75','u95'])/mortality_rates[age_grp=='65-75','u95'], digits=2), '%'))
  }
  if(100*abs(quantile(mortality_samples[age_grp=='75+']$value, probs=0.975) - mortality_rates[age_grp=='75+','u95'])/mortality_rates[age_grp=='75+','u95'] > 25){
    print(paste0('U95 of sample more than 25% out (75+): ',
                 round(100*abs(quantile(mortality_samples[age_grp=='75+']$value, probs=0.975) - mortality_rates[age_grp=='75+','u95'])/mortality_rates[age_grp=='75+','u95'], digits=2), '%'))
  }
  
  ## 3: obtaining 65+ mortality rate by weighting by 2015 age-specific population
  country_names <- unlist(unname(clusters[iso3c==country_code, c('country','country_altern','country_altern_2')]))
  if(country_code=='TUR'){country_names <- c(country_names, 'Turkiye')}
  age_pop <- data.table(country = country_names[1],
                        lower.age.limit = seq(0,100,by=5),
                        year = 2015,
                        population = 1000*unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))
  
  if(!nrow(age_pop) == 21){stop('age_pop wrong number of rows!')}
  
  weighted_vals <- (mortality_samples[age_grp=='65-75']$value * sum(age_pop[lower.age.limit %in% c(65,70), 'population']) + 
                      mortality_samples[age_grp=='75+']$value * sum(age_pop[lower.age.limit > 71, 'population']))/sum(age_pop[lower.age.limit>64, 'population'])
  
  mortality_samples <- rbind(mortality_samples, data.table(age_grp = rep('65+', 100),
                                                           simulation_index = 1:100,
                                                           value = weighted_vals))
  setnames(mortality_samples, 'value', 'mortality_rate')
  # ggplot(mortality_samples) +
  #   geom_boxplot(aes(x=age_grp, y=value)) + scale_y_log10()
  
  ## 4: obtaining age-specific all-cause mortality distribution (<65s)
  unique(mx1dt[(grepl('Turk', mx1dt$name))]$name)
  UNLT <- data.table(mx1dt[name %in% country_names & year == 2015, 
                           c('iso3c','age','mxB')])
  setnames(UNLT, 'mxB', 'mx')
  if(!nrow(UNLT) == 101){stop('UNLT wrong number of rows!')}
  
  UNLT$population <- c(age_transform(age_pop$population))
  
  UNLT[age < incidence_ages[2], age_grp := 1]
  UNLT[age %in% incidence_ages[2]:(incidence_ages[3] - 1), age_grp := 2]
  UNLT[age %in% incidence_ages[3]:(incidence_ages[4] - 1), age_grp := 3]
  UNLT[age >= incidence_ages[4], age_grp := 4]
  
  age_grp_pop <- data.table(age_grp = c(1:4), pop = c(sum(age_pop[lower.age.limit < incidence_ages[2]]$population),
                                                      sum(age_pop[lower.age.limit %in% incidence_ages[2]:(incidence_ages[3] - 1)]$population),
                                                      sum(age_pop[lower.age.limit %in% incidence_ages[3]:(incidence_ages[4] - 1)]$population),
                                                      sum(age_pop[lower.age.limit >= incidence_ages[4]]$population)))
  if(!sum(age_grp_pop$pop) == sum(age_pop$population)){stop("Age-specific populations don't match")}
  
  UNLT <- UNLT[age_grp_pop, on=c('age_grp'), pop_agg := i.pop][,prop_of_age_grp := population/pop_agg]
  for(i in unique(UNLT$age_grp)){
    UNLT[age_grp==i, weighted_mx := sum(mx*prop_of_age_grp)]
  }
  
  # ggplot(UNLT) + 
  #   geom_line(aes(x=age, y=mx), lwd=1) +
  #   geom_line(aes(x=age, y=(weighted_mx)), lwd=1, lty=2) +
  #   scale_y_log10() 
  
  all_cause_weighted <- UNLT[age %in% incidence_ages, c('age_grp','weighted_mx')]
  
  ## 5: making over-65s ifrs
  
  over65 <- mortality_samples[age_grp=='65+',]
  over65 <- over65[incidence[age_grp==4,], on=c('simulation_index'),
                   infections := i.value]
  over65[, population := age_grp_pop[age_grp==4,'pop']]
  over65[, incid_rate := infections/population][, ifr := mortality_rate/incid_rate]
    
  ## 6: making under-65s ifrs
  
  u65incidence <- incidence[age_grp %in% 1:3,]
  u65incidence <- u65incidence[all_cause_weighted[age_grp %in% 1:3], on = c('age_grp')]
  u65incidence[, unscaled_deaths := value*weighted_mx]
  
  under65scaling <- u65incidence[, c('simulation_index','unscaled_deaths')]
  under65scaling <- under65scaling[, lapply(.SD, sum), by=c('simulation_index')]
  under65scaling[, u65pop := sum(age_grp_pop[age_grp %in% 1:3, 'pop'])]
  under65scaling[, unscaled_mort_rate := unscaled_deaths/u65pop]
  under65scaling <- under65scaling[mortality_samples[age_grp=='<65'], on = c('simulation_index')]
  under65scaling[, scaled_deaths := mortality_rate*u65pop]
  under65scaling[, scaling_factor := mortality_rate/unscaled_mort_rate]
  
  under65ifr <- incidence[age_grp %in% 1:3, c('simulation_index','age_grp','value')]
  under65ifr <- under65ifr[all_cause_weighted[age_grp %in% 1:3], on=c('age_grp'), 
                           weighted_mx := i.weighted_mx]
  under65ifr <- under65ifr[under65scaling, on = c('simulation_index'),
                           scaling_factor := i.scaling_factor]
  under65ifr <- under65ifr[age_grp_pop[age_grp %in% 1:3], on = c('age_grp')]
  under65ifr[, ifr:=weighted_mx*scaling_factor][, deaths:=ifr*value][,mortality_rate:=deaths/pop]
  
  ## 7: merging ifrs
  over65[,age_grp:=4]
  ifrs <- rbind(under65ifr[, c('age_grp','simulation_index','ifr')],
                over65[, c('age_grp','simulation_index','ifr')])
  # ggplot(ifrs) +
  #   geom_boxplot(aes(x=age_grp, y=ifr*100000, fill=age_grp)) +
  #   scale_y_log10()
  
  # check mortality rate estimates are right:
  if(abs(sum(quantile(unlist(unname((data.frame(under65ifr) %>% group_by(simulation_index) %>% summarise(sum(deaths)))[,2]/sum(age_grp_pop[age_grp %in% 1:3]$pop)))) -
             quantile(mortality_samples[age_grp=='<65']$mortality_rate))) > 0.00001){stop("Mortality rates don't line up (<65)")}
  if(abs(sum(quantile(unlist(unname((data.frame(over65) %>% group_by(simulation_index) %>% summarise(sum(infections*ifr)))[,2]/sum(age_grp_pop[age_grp == 4]$pop)))) -
    quantile(mortality_samples[age_grp=='65+']$mortality_rate))) > 0.00001){stop("Mortality rates don't line up (65+)")}
  
  return(ifrs)
}










