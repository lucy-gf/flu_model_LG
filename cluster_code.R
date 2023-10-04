## LOADING PACKAGES (which are manually installed on the cluster)

library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(data.table)
library(odin)
library(parallel)

## ** GOING TO LOAD IN ALL OF THE NECESSARY DATASETS ETC. BY HAND **

# data needed: 

df_cntr_table <- read.csv("data_for_cluster/df_cntr_table.csv") # saved from epid_identif_cont_matrs.R
list_contact_matr <- readRDS("data_for_cluster/list_contact_matr.RDS") # saved from output
df_epid_threshold <- read.csv("data_for_cluster/df_epid_threshold.csv") # saved from output
matches <- read.csv("data_for_cluster/matches.csv") # saved from vaccine_inputs.R
cov_data <- read.csv("data_for_cluster/cov_data.csv") # saved from vaccine_inputs.R

age_limits <- c(0,5,20,65); age_group_names <- paste0(age_limits,"-", c(age_limits[2:length(age_limits)],99))
infection_delays <- c( 0.8, 1.8 ) # 0.8 and 1.8 day

# also used "scp -r fcns lshlg12@hpclogin.lshtm.ac.uk:" to upload the fcns folder
# will then probably need to source these:
source("fcns/inference_function.R")
source("fcns/2_1b_model_epidemic_yearcross.R")

# MCMC parameters

post_size <- 1000 #5000 #3000 
thinning_steps <- 1 #100
burn_in <- 1 #100000
seed_to_use <- 55 #99

strain <- 'INF_A'

# ** FUNCTION TO BE USED WITH LAPPLY/MCLAPPLY **

cluster_mcmc <- function(k){ # k %in% 1:7
  
  sel_cntr <- df_cntr_table$country[k] # setting the selected country
  
  contact_matrix <- as.matrix(unname(list_contact_matr[[sel_cntr]]))
  
  yr_res_pop <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))
  group_pop <- pop_age(wpp_age(sel_cntr, 2015), age.limits=c(0,5,20,65))$population # size of each age group
  contact_matrix_small <- t(t(contact_matrix)/group_pop) ### CHECK THIS
  
  # lets first subset data for one country and NONSENT data only
  
  data_fitting <- df_epid_threshold %>% 
    filter(grepl("NONSENT",metasource) & country %in% sel_cntr & STRAIN %in% "INF_A") %>%
    select(!c(over_peak,flu_included,over_inclusion,positivity,flu_peak,seq))

  epidemics_to_fit <- list()
  for (epid_index_val in unique(data_fitting$epid_index[data_fitting$epidem_inclusion>0])) {
    xx <- data_fitting %>% filter(epid_index %in% epid_index_val)
    epidemics_to_fit[[epid_index_val]] <- list(
      start=min(xx$ISO_WEEKSTARTDATE),
      end = max(xx$ISO_WEEKSTARTDATE),
      initial_params = c(-8, 10, 0.4, 4, 0, 0),# c(-7, 10, 0.5, 3, 0, 0), # c(-8, 10, 0.4, 4, 0, 0),
      data_points = xx$value,
      type = "A"
    )
  }

  dates_to_run <- list()
  for(epid_index_val in 1:length(epidemics_to_fit)){
    dates_to_run[[epid_index_val]] <- c(
      as.Date(epidemics_to_fit[[epid_index_val]]$start),
      as.Date(epidemics_to_fit[[epid_index_val]]$end) + 7 # need end date to be the date AFTER the last week we want to model
    )
  }
  # dates_to_run <- as.Date(dates_to_run,origin="1970-01-01")

  # vacc calendar (no fitting in the vaccination)
  vaccine_calendar_list <- list()
  for(epid_index_val in 1:length(epidemics_to_fit)){
    vaccine_calendar_list[[epid_index_val]] <- as_vaccination_calendar(
      efficacy = rep(0,length(age_group_names)),
      dates = dates_to_run[[epid_index_val]],
      coverage = matrix(
        0,
        nrow = 3, #length(dates_to_run[[epid_index_val]]),
        ncol = length(age_group_names)
      ),
      no_age_groups = length(age_group_names),
      no_risk_groups=1
    )
  }
  
  initial_parameters <- list()
  for(epid_index_val in 1:length(epidemics_to_fit)){
    initial_parameters[[epid_index_val]] <- epidemics_to_fit[[epid_index_val]]$initial_params
    names(initial_parameters[[epid_index_val]]) <- c("reporting", "transmissibility","susceptibility","initial_infected")
  }
  
  ##### Run the fit #####

  post_size <- post_size #3000 
  thinning_steps <- thinning_steps #100
  burn_in <- burn_in #100000
  seed_to_use <- seed_to_use #99
  
  #for (epidemic_to_run in 1:length(epidemics_to_fit)){
  set.seed(seed_to_use)
  global_index <<- 1
  output <- custom_inference( # this calls "fcns/inference_function.R"
    input_demography = yr_res_pop,
    vaccine_calendar_list = vaccine_calendar_list, # vaccine_calendar_list[[epidemic_to_run]],
    input_polymod = contact_matrix_small, #contacts_matrixformat, 
    ili = NULL,
    mon_pop = NULL,
    # n_pos = epidemics_to_fit[[epidemic_to_run]]$data_points,
    epidemics_to_fit = epidemics_to_fit, #epidemics_to_fit[[epidemic_to_run]], # passing in all the epidemics instead of n_pos
    n_samples = NULL,
    initial = unlist(initial_parameters), #initial_parameters[[epidemic_to_run]],
    mapping = NULL,
    nbatch = post_size,
    nburn = burn_in,
    blen = thinning_steps,
    sel_cntr = sel_cntr
  )
  
  return(output) # output of the function 'cluster_mcmc'
  
}

results <- mclapply(1:7, cluster_mcmc)#, mc.cores = 4) # will update to include both strains/sources in future
# mclapply took 50mins for 7x500 steps
# 11hrs for 7x10,000 steps

names(results) <- df_cntr_table$country

saveRDS(results, file = "cluster_data_29-09.rds")

# this is binomial likelihood with no vaccinations and
# updated priors, over 200,000 steps






