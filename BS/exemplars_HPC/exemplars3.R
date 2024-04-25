
### ADDING INITIAL R0 TO SAMPLED EPIDEMICS ###
library(fluEvidenceSynthesis)
library(dplyr)
library(data.table)
library(ggplot2)
source("BS/BS_colors.R")
source("BS/BS_data_fcns.R")
source("BS/BS_vaccine_programs.R")

length_simulation <- 10
n_simulations <- 1000

matches <- read.csv("data_for_BS/matches.csv") # saved from vaccine_inputs.R

sampled_years_df <- data.frame(index = rep(1:n_simulations,each=10),
                               years = rep(2010:2019, n_simulations),
                               N_A_match = rep(matches$match[matches$hemisphere == 'NH' &
                                                               matches$strain_match == 'INF_A'], n_simulations),
                               N_B_match = rep(matches$match[matches$hemisphere == 'NH' &
                                                               matches$strain_match == 'INF_B'], n_simulations),
                               S_A_match = rep(matches$match[matches$hemisphere == 'SH' &
                                                               matches$strain_match == 'INF_A'], n_simulations),
                               S_B_match = rep(matches$match[matches$hemisphere == 'SH' &
                                                               matches$strain_match == 'INF_B'], n_simulations))

cntr_num <- 3

  country_code_index <- c("ARG", "AUS", "CAN", "CHN", "GBR", "GHA", "TUR")[cntr_num]
  hemisphere_input <- c('SH', 'SH', 'NH', 'NH', 'NH', 'NH', 'NH')[cntr_num]

  sampled_epidemics <- data.frame() # will be around 45,000 epidemics
  print(country_code_index)
  for(index in (1:n_simulations)){
    epids_to_add <- fcn_sample_epidemics(country_code = country_code_index,
                                         sampled_years = sampled_years_df[sampled_years_df$index == index,]) %>%
      mutate(simulation_index = index, pushback = NA, init_ageing_date = NA, init_nye = NA, pushback_T = F)
    code <- epids_to_add$exemplar_code[1]
    country_name <- countrycode(code, origin = 'iso3c', destination = 'country.name')
    country_altern <- clusters$country_altern[clusters$codes == code]
    contact_matrix <- as.matrix(unname(list_contact_matr[[country_name]]))
    yr_res_pop <- unlist(lapply(pop_age(wpp_age(country_name, 2015))$population/5, function(x) rep(x,5)))
    group_pop <- pop_age(wpp_age(country_name, 2015), age.limits=c(0,5,20,65))$population # size of each age group
    contact_matrix_small <- t(t(contact_matrix)/group_pop)
    for(i in 1:nrow(epids_to_add)){
      original_date <- as.Date(paste0(as.numeric(epids_to_add$day[i]), '-', epids_to_add$month[i], '-',
                                      epids_to_add$year[i]), '%d-%m-%Y')
      pb_input <- fcn_identify_start(sus_input = unlist(epids_to_add$sus[i]),
                                     trans_input = unlist(epids_to_add$trans[i]),
                                     infected = unlist(epids_to_add$infected[i]),
                                     date = original_date,
                                     country = country_name,
                                     strain_input = epids_to_add$strain[i],
                                     yr_res_pop = yr_res_pop,
                                     contact_matrix_small = contact_matrix_small,
                                     ageing_dates = c('01-04','01-10'),
                                     hemisphere = hemisphere_input,
                                     simulation_year = unlist(epids_to_add$simulation_cal_year[i]))
      if(length(pb_input) == 1){
        epids_to_add$pushback[i] <- pb_input
      }else{
        epids_to_add$pushback_T[i] <- F
        if(pb_input[2] == 'ageing_date'){
          epids_to_add$init_ageing_date[i] <- pb_input[1]
        }else{
          epids_to_add$init_nye[i] <- pb_input[1]
        }}
    }
    sampled_epidemics <- rbind(sampled_epidemics,epids_to_add)
    if(index %% 10 == 0){
      print(paste0('# simulation = ', index))
      write_csv(sampled_epidemics, file = paste0("data_for_BS/original_epidemics_", length_simulation,
                                                 "_", n_simulations,'_',country_code_index,".csv"))
    }
  }


start_year <- 2010
years <- 10


  country_code_index <- c("ARG", "AUS", "CAN", "CHN", "GBR", "GHA", "TUR")[cntr_num]
  hemisphere_input <- c('SH', 'SH', 'NH', 'NH', 'NH', 'NH', 'NH')[cntr_num]
  hemisphere_vacc <- hemisphere_input

  sampled_epidemics <- read_csv(paste0("data_for_BS/original_epidemics_", length_simulation,
                                             "_", n_simulations,'_',country_code_index,".csv"))

  ITZ <- clusters %>% filter(codes == country_code_index)
  n <- nrow(ITZ)
  
  monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                          to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                          by = 1)) == 'Monday')
  n_weeks <- length(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                             to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                             by = 7))
  
  scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_none')]][1]
  scenarios[[1]]$init_vaccinated <- unlist(unname(cov_data[cov_data$country==ITZ$country &
                                               cov_data$year==start_year,3:6]))
  
  cases_rds <- list()
  
  for(vacc_prog_num in 1:length(scenarios)){
    
    cases_30_100 <- data.frame(country = rep(rep(unique(ITZ$country), each=n_weeks), n_simulations),
                               country_code = rep(rep(unique(ITZ$codes), each=n_weeks), n_simulations),
                               simulation_index = rep(1:n_simulations, each = n*n_weeks),
                               week = rep(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                                                   to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                                                   by = 7), length(unique(ITZ$country))*n_simulations),
                               IU1A = 0, IU2A = 0, IU3A = 0, IU4A = 0,
                               IV1A = 0, IV2A = 0, IV3A = 0, IV4A = 0,
                               IU1B = 0, IU2B = 0, IU3B = 0, IU4B = 0,
                               IV1B = 0, IV2B = 0, IV3B = 0, IV4B = 0)
    start_time <- Sys.time()
    for(sims in 1:n_simulations){
      for(country_code_index in unique(cases_30_100$country_code)){
        country_time <- Sys.time()
        cases_30_100[cases_30_100$simulation_index == sims &
                       cases_30_100$country_code == country_code_index, 5:20] <-
          fcn_simulate_epidemics_exemplar(country_code_input = country_code_index,
                                 hemisphere = hemisphere_vacc,
                                 sampled_epidemics_input =
                                   sampled_epidemics[sampled_epidemics$simulation_index == sims,],
                                 vacc_prog = scenarios[[vacc_prog_num]]) %>% select(!week)
        #print(paste0(countrycode(country_code_index, origin='iso3c', destination = 'country.name')))
        #print(Sys.time() - country_time)
      }
      if(sims %% 50 == 0){
        cases_rds[[vacc_prog_num]] <- cases_30_100
        if(vacc_prog_num == length(scenarios)){
          names(cases_rds) <- names(scenarios)
        }
        saveRDS(cases_rds, file = paste0("data/exemplar_output/exemplar_", country_code_index,
                                         "_2010_2019.rds"))
      }
      if(sims %% 10 == 0){
        start_time <- Sys.time()
        print(paste0('Simulations complete: ', sims,
                     ', time taken: ', round(Sys.time() - start_time, digits = 2))) 
      }
    }
}




