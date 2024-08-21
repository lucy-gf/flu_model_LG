
### RUNNING ORIGINAL (2010-2019) EPIDEMICS ###
# setwd("~/Desktop/research asst/Global Code")
library(fluEvidenceSynthesis)
library(dplyr)
library(data.table)
library(ggplot2)
source("BS/BS_colors.R")
source("BS/BS_vaccine_programs.R")
source("BS/BS_data_fcns.R")

start_year <- 2010
years <- 10
length_simulation <- years
n_simulations <- 100

ifr_method <- c('exemplar','whole_itz','brazil')[3]

vec <- case_when(ifr_method=='exemplar'~list(c(1:7)), 
                 ifr_method=='whole_itz'~list(c(4,6,7)),
                 ifr_method=='brazil'~list(c(1)))[[1]]

for(cntr_num in vec){

country_code_index <- c("ARG", "AUS", "CAN", "CHN", "GBR", "GHA", "TUR")[cntr_num]
cluster_code <- country_code_index
hemisphere_input <- c('SH', 'SH', 'NH', 'NH', 'NH', 'NH', 'NH')[cntr_num]
hemisphere_vacc <- hemisphere_input
if(ifr_method=='brazil'){
  country_code_index <- "ARG"
  hemisphere_input <- 'SH'
  hemisphere_vacc <- hemisphere_input
}
if(ifr_method=='whole_itz' & !country_code_index %in% c('CHN','GHA','TUR')){
  stop("ITZ not suitable for 'whole_itz' method")
}

sampled_epidemics <- read_csv(paste0("data_for_BS/original_epidemics_", length_simulation,
                                     "_100_",country_code_index,"_wr0.csv"))

ITZ <- NULL

if(ifr_method == 'exemplar'){
  ITZ <- clusters %>% filter(codes == country_code_index)
}
if(ifr_method == 'whole_itz'){
  ITZ <- clusters %>% filter(cluster_code == country_code_index)
}
if(ifr_method == 'brazil'){
  ITZ <- clusters %>% filter(codes == 'BRA')
}

n <- nrow(ITZ)

monday_start <- which(weekdays(seq.Date(from = as.Date(paste0("01-01-", start_year), '%d-%m-%Y'),
                                        to = as.Date(paste0("07-01-", start_year), '%d-%m-%Y'),
                                        by = 1)) == 'Monday')
n_weeks <- length(seq.Date(from = as.Date(paste0(monday_start, "-01-", start_year), '%d-%m-%Y'),
                           to = as.Date(paste0("31-12-", (start_year + years - 1)), '%d-%m-%Y'),
                           by = 7))

scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_none')]][1]

print(paste0(ifr_method, ', ', country_code_index))

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
        fcn_simulate_original_epidemics(country_code_input = country_code_index,
                                        hemisphere = hemisphere_vacc,
                                        sampled_epidemics_input =
                                          sampled_epidemics[sampled_epidemics$simulation_index == sims,],
                                        vacc_prog = scenarios[[vacc_prog_num]]) %>% select(!week)
      #print(Sys.time() - country_time)
    }
    if(sims %% 10 == 0){
      cases_rds[[vacc_prog_num]] <- cases_30_100
      if(vacc_prog_num == length(scenarios)){
        names(cases_rds) <- names(scenarios)
      }
      saveRDS(cases_rds, file = paste0("data/original_epids_output/",
                                       ifr_method, "_", cluster_code,
                                       "_2010_2019_",n_simulations,".rds"))
    }
    if(sims %% 10 == 0){
      print(paste0('Simulations complete: ', sims,
                   ', time taken: ', round(Sys.time() - start_time, digits = 2)))
      start_time <- Sys.time()
    }
  }
}
}











