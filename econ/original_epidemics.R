
### ADDING INITIAL R0 TO SAMPLED EPIDEMICS ###
library(fluEvidenceSynthesis)
library(dplyr)
library(data.table)
library(ggplot2)
source("BS/BS_colors.R")
source("BS/BS_data_fcns.R")
source("BS/BS_vaccine_programs.R")

length_simulation <- 10
n_simulations <- 100

matches <- read.csv("data_for_BS/matches.csv") # saved from vaccine_inputs.R
post_samples_merged_simulation <- readRDS("data_for_BS/post_samples_merged_simulation.rds")

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

cntr_num <- 1

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

  ### ADDING R0 ###
  
  for(c_num in 1:7){
    
    cluster_code <- c('ARG','AUS','CAN','CHN','GHA','TUR','GBR')[c_num]
    original_epidemics <- read_csv(paste0("data_for_BS/original_epidemics_10_100_",cluster_code, ".csv"))
    sel_cntr <- c('Argentina','Australia','Canada','China','Ghana','Turkey','United Kingdom')[c_num]
    
    contact_matrix <- as.matrix(unname(list_contact_matr[[c_num]]))
    yr_res_pop <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))
    group_ages <- stratify_by_age(yr_res_pop, limits = age.group.limits)
    contact_matrix_small <- t(t(contact_matrix)/group_ages) 
    
    original_epidemics$r0 <- NA
    
    for(k_row in 1:nrow(original_epidemics)){
      original_epidemics$r0[k_row] <- fluEvidenceSynthesis::as_R0(
        transmission_rate = original_epidemics$trans[k_row],
        contact_matrix = contact_matrix_small, 
        # age_groups = group_ages*c((0.2*1 + 0.8*sampled_epidemics$sus[k_row]),
        #                           rep(sampled_epidemics$sus[k_row],3))
        age_groups = group_ages
      )
    }
    write_csv(original_epidemics, file = paste0("data_for_BS/original_epidemics_10_100_",cluster_code,"_wr0.csv"))
  }
  
  original_epidemics <- data.frame()
  for(cluster_code in c('ARG','AUS','CAN','CHN','GHA','TUR','GBR')){
    original_epidemics <- rbind(original_epidemics,
                                read_csv(paste0("data_for_BS/original_epidemics_10_100_",cluster_code, "_wr0.csv")))
  }
  
  original_epidemics %>% ggplot() +
    geom_point(aes(x=trans, y=r0, col=exemplar_code)) +
    scale_color_manual(values=exemplar_colors) +
    theme_minimal() + xlim(c(0,NA)) + ylim(c(0,NA)) +
    # xlab('Susceptibility') +
    ylab('R0') +
    theme_minimal() +
    # facet_grid(exemplar_code~strain, scales='fixed') +
    labs(col='Country') +
    theme(text=element_text(size=14))
  
  original_epidemics %>% ggplot() +
    geom_histogram(aes(x=r0, fill=exemplar_code), bins=100) +
    scale_fill_manual(values=exemplar_colors) +
    theme_bw() +
    # facet_grid(exemplar_code~., scales='fixed') +
    xlab('R0') + ylab('Density') +
    theme(text=element_text(size=14))
  




