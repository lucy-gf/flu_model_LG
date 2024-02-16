
## SIMULATING X-YEAR TIME SERIES FOR EACH CLUSTER/COUNTRY IN INF A+B ##

## RUNNING WITH VARIOUS VACCINATION SCENARIOS ##
#setwd("~/Desktop/research asst/Global Code")
source("BS/BS_data_fcns.R")
source("BS/BS_vaccine_programs.R")
length_simulation <- 30
n_simulations <- 100

start_year <- 2025 
years <- 30

## MAKE SURE THAT THE ORDER OF LOOPS ENSURES THE LEAST COMPUTATIONAL TIME

# loops through each ITZ and vaccination program:
for(c_number in c(4)){ #1:7){
    c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
                "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
                "Southern America")[c_number]
    cluster_code <- c('GHA','TUR','CHN','GBR','CAN','AUS','ARG')[c_number]
    hemisphere_vacc <- c('NH','NH','NH','NH','NH','SH','SH')[c_number]

    sampled_epidemics <- rbind(read_csv(paste0("data_for_BS/sampled_epidemics_", length_simulation,
                                               "_", n_simulations,"_",cluster_code, ".csv")))

    ITZ <- clusters %>% filter(cluster_name == c_name)
    n <- nrow(ITZ)

    monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                            to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                            by = 1)) == 'Monday')
    n_weeks <- length(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                               to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                               by = 7))
    
    cases_rds <- list()
    
    for(vacc_prog_num in 1:length(vaccine_programs)){
      
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

    for(sims in 1:n_simulations){
      start_time <- Sys.time()
      for(country_code_index in unique(cases_30_100$country_code)){
        country_time <- Sys.time()
        cases_30_100[cases_30_100$simulation_index == sims &
                       cases_30_100$country_code == country_code_index, 5:20] <-
          fcn_simulate_epidemics(country_code_input = country_code_index,
                                 hemisphere = hemisphere_vacc,
                                 sampled_epidemics_input =
                                   sampled_epidemics[sampled_epidemics$simulation_index == sims,],
                                 vacc_prog = vaccine_programs[[vacc_prog_num]]) %>% select(!week)
        print(paste0(countrycode(country_code_index, origin='iso3c', destination = 'country.name')))
        print(Sys.time() - country_time)
      }
      if(sims %% 50 == 0){
        cases_rds[[vacc_prog_num]] <- cases_30_100
        saveRDS(cases_rds, file = paste0("data/vacc_output/vacc_", cluster_code, "_waning_changed_high_cov.rds"))
        # write_csv(cases_30_100, file = paste0("data/vacc_output/vacc_", cluster_code, "_vp_", vacc_prog_num, ".csv"))
      }
      print(paste0('Vaccination program: ', vacc_prog_num,
                   ', Simulations complete: ', sims,
                   ', time taken: ', round(Sys.time() - start_time, digits = 2)))
    }
  }
}


# ## ONLY RUNNING ONE COUNTRY
# c_number <- 5
# n_simulations <- 100
#   c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
#               "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
#               "Southern America")[c_number]
#   cluster_code <- c('GHA','TUR','CHN','GBR','CAN','AUS','ARG')[c_number]
#   hemisphere_vacc <- c('NH','NH','NH','NH','NH','SH','SH')[c_number]
# 
#   sampled_epidemics <- rbind(read_csv(paste0("data_for_BS/sampled_epidemics_30_100_",cluster_code, ".csv")))
# 
#   ITZ <- clusters %>% filter(cluster_name == c_name,
#                              codes=='CAN')
#   n <- nrow(ITZ)
# 
#   monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
#                                           to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
#                                           by = 1)) == 'Monday')
#   n_weeks <- length(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
#                              to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
#                              by = 7))
# 
#   for(vacc_prog_num in 1:length(vaccine_programs)){
#   cases_30_100 <- data.frame(country = rep(rep(unique(ITZ$country), each=n_weeks), n_simulations),
#                              country_code = rep(rep(unique(ITZ$codes), each=n_weeks), n_simulations),
#                              simulation_index = rep(1:n_simulations, each = n*n_weeks),
#                              week = rep(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
#                                                  to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
#                                                  by = 7), length(unique(ITZ$country))*n_simulations),
#                              IU1A = 0, IU2A = 0, IU3A = 0, IU4A = 0,
#                              IV1A = 0, IV2A = 0, IV3A = 0, IV4A = 0,
#                              IU1B = 0, IU2B = 0, IU3B = 0, IU4B = 0,
#                              IV1B = 0, IV2B = 0, IV3B = 0, IV4B = 0)
# 
#   for(sims in 1:n_simulations){
#     start_time <- Sys.time()
#     for(country_code_index in unique(cases_30_100$country_code)){
#       country_time <- Sys.time()
#       cases_30_100[cases_30_100$simulation_index == sims &
#                    cases_30_100$country_code == country_code_index, 5:20] <-
#         fcn_simulate_epidemics(country_code_input = country_code_index,
#                                hemisphere = hemisphere_vacc,
#                                sampled_epidemics_input =
#                                sampled_epidemics[sampled_epidemics$simulation_index == sims,],
#                                vacc_prog = vaccine_programs[[vacc_prog_num]]) %>% select(!week)
#       print(paste0(countrycode(country_code_index, origin='iso3c', destination = 'country.name')))
#       print(Sys.time() - country_time)
#     }
#     write_csv(cases_30_100, file = paste0("data/vacc_output/vacc_", cluster_code, "_vp_", vacc_prog_num, ".csv"))
#     print(paste0('Vaccination program: ', vacc_prog_num,
#                  ', Simulations complete: ', sims,
#                  ', time taken: ', round(Sys.time() - start_time, digits = 2)))
#   }
# }



# # ## CHECK:
# cases_30_100 %>%
#   filter(simulation_index == sims,
#          year(week)>2030) %>%
#   ggplot() +
#   geom_line(aes(x=week, y=IV1A)) +
#   theme_minimal() +
#   geom_vline(xintercept = as.Date(paste0("01-10-", 2026:2054), '%d-%m-%Y'),
#              lty=3, col=2) +
#   geom_vline(xintercept = as.Date(paste0("01-04-", 2026:2054), '%d-%m-%Y'),
#              lty=3, col=3)






