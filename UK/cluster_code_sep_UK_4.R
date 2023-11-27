
### EXAMPLE OF FITTING FOR 2 EPIDEMICS IN 1 COUNTRY

## LOADING PACKAGES

library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(data.table)
library(odin)
library(parallel)

## ** LOADING IN DATA **

df_cntr_table <- read.csv("data_for_cluster/df_cntr_table.csv") # saved from epid_identif_cont_matrs.R
list_contact_matr <- readRDS("data_for_cluster/list_contact_matr.RDS") # saved from output
df_epid_threshold <- read.csv("data_for_cluster/df_epid_threshold.csv") # saved from output
matches <- read.csv("data_for_cluster/matches.csv") # saved from vaccine_inputs.R
cov_data <- read.csv("data_for_cluster/cov_data.csv") # saved from vaccine_inputs.R

age_limits <- c(0,5,20,65); age_group_names <- paste0(age_limits,"-", c(age_limits[2:length(age_limits)],99))
infection_delays <- c( 0.8, 1.8 ) # 0.8 and 1.8 day

# ** SOURCE FILES **
source("fcns/inference_function.R")
source("fcns/2_1b_model_epidemic_yearcross.R")

# MCMC parameters

post_size <- 10000 #5000 #3000 
thinning_steps <- 50 #50
burn_in <- 500000
seed_to_use <- 55 #99

strain <- c('INF_A', 'INF_B')[2]

# ** FUNCTION TO BE USED WITH LAPPLY/MCLAPPLY **

cluster_mcmc <- function(k, epidemic_to_run){ # k %in% 1:7
  
  sel_cntr <- df_cntr_table$country[k] # setting the selected country
  
  contact_matrix <- as.matrix(unname(list_contact_matr[[sel_cntr]]))
  
  yr_res_pop <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))
  group_pop <- pop_age(wpp_age(sel_cntr, 2015), age.limits=c(0,5,20,65))$population # size of each age group
  contact_matrix_small <- t(t(contact_matrix)/group_pop) 
  
  # lets first subset data for one country and NONSENT data only
  
  data_fitting <- df_epid_threshold %>% 
    filter(grepl("NONSENT",metasource) & country %in% sel_cntr & STRAIN %in% strain) %>%
    select(!c(over_peak,flu_included,over_inclusion,flu_peak,seq))

  epidemics_to_fit <- list()
  for (epid_index_val in unique(data_fitting$epid_index[data_fitting$epidem_inclusion>0])) {
    xx <- data_fitting %>% filter(epid_index %in% epid_index_val)
    epidemics_to_fit[[epid_index_val]] <- list(
      start=min(xx$ISO_WEEKSTARTDATE),
      end = max(xx$ISO_WEEKSTARTDATE),
      initial_params = c(-8, 10, 0.4, 4, 0, 0),
      weeks = xx$ISO_WEEKSTARTDATE,
      data_points = xx$value,
      type = strain
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

  post_size <- post_size 
  thinning_steps <- thinning_steps 
  burn_in <- burn_in
  seed_to_use <- seed_to_use 

  set.seed(seed_to_use)
  global_index <<- 1
  output <- custom_inference_sep(
    input_demography = yr_res_pop,
    vaccine_calendar = vaccine_calendar_list[[epidemic_to_run]],
    input_polymod = contact_matrix_small, #contacts_matrixformat,
    ili = NULL,
    mon_pop = NULL,
    # n_pos = epidemics_to_fit[[epidemic_to_run]]$data_points,
    epidemics_to_fit = epidemics_to_fit[[epidemic_to_run]], # passing in all the epidemics instead of n_pos      n_samples = NULL,
    initial = initial_parameters[[epidemic_to_run]],
    mapping = NULL,
    nbatch = post_size,
    nburn = burn_in,
    blen = thinning_steps,
    sel_cntr = sel_cntr,
    hemisphere_input = df_cntr_table[df_cntr_table$country == sel_cntr,]$hemisphere
  )
  
  return(output) # each epidemic's batch & llikelihoods
}

a <- 7 # CHANGE THIS TO CHANGE EXEMPLAR COUNTRY
epidemic_to_run <- 1 # CHANGE THIS TO RUN A DIFFERENT EPIDEMIC

## run first epidemic
results <- cluster_mcmc(a, epidemic_to_run)

sel_cntr <- df_cntr_table$country[a]

saveRDS(results, file = paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
                               epidemic_to_run, "_",
                               post_size, "_", thinning_steps, "_", burn_in, ".rds"))

## run second epidemic
epidemic_to_run <- epidemic_to_run + 1

results2 <- cluster_mcmc(a, epidemic_to_run)

saveRDS(results2, file = paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
                               epidemic_to_run, "_",
                               post_size, "_", thinning_steps, "_", burn_in, ".rds"))

## run third epidemic
epidemic_to_run <- epidemic_to_run + 1

results3 <- cluster_mcmc(a, epidemic_to_run)

saveRDS(results3, file = paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
                                epidemic_to_run, "_",
                                post_size, "_", thinning_steps, "_", burn_in, ".rds"))

# ## run fourth epidemic
# epidemic_to_run <- epidemic_to_run + 1
# 
# results4 <- cluster_mcmc(a, epidemic_to_run)
# 
# saveRDS(results4, file = paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
#                                 epidemic_to_run, "_",
#                                 post_size, "_", thinning_steps, "_", burn_in, ".rds"))
# 
# ## run fifth epidemic
# epidemic_to_run <- epidemic_to_run + 1
# 
# results5 <- cluster_mcmc(a, epidemic_to_run)
# 
# saveRDS(results5, file = paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
#                                 epidemic_to_run, "_",
#                                 post_size, "_", thinning_steps, "_", burn_in, ".rds"))






