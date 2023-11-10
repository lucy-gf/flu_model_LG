
# READING IN AND VIEWING RESULTS FROM THE COMMAND LINE

library(data.table)
library(bayestestR)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(odin)
library(patchwork)
library(viridis)

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

seed_to_use <- 55 
strain <- c('INF_A','INF_B')[1]
post_size <- 10000 #5000 #3000 
thinning_steps <- 50 #100
burn_in <- 200000 #100000
sel_cntr <- 'Turkey'

# DATAFRAME   

post_samples_merged <- data.table()

for(i in 1:df_cntr_table[,paste0(strain,'_NONSENT')][df_cntr_table$country==sel_cntr]){
  epid <- i
  if(file.exists(paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
                        epid, "_",
                        post_size, "_", thinning_steps, "_", burn_in, ".rds"))){
    results_RDS <- readRDS(paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
                                  epid, "_",
                                  post_size, "_", thinning_steps, "_", burn_in, ".rds"))
    post_samples <- data.table(results_RDS$batch)
    # colnames(post_samples) <- rep(c("reporting","trans","sus","infected","blank","blank"),length(epidemics_to_fit))
    colnames(post_samples) <- paste0(c("reporting","trans","sus","infected","blank","blank"))
    # remove the blank ones
    cols_to_delete <- which(names(post_samples) %like% "blank")
    post_samples[, (cols_to_delete) := NULL ]
    # make into real values
    post_samples$reporting <- exp(post_samples$reporting)
    post_samples$trans <- post_samples$trans/100
    post_samples$infected <- 10^(post_samples$infected)
    # add the likelihoods
    post_samples[, ll := results_RDS$llikelihoods]
    post_samples$timestep <- 1:nrow(post_samples)
    post_samples_m <- melt.data.table(post_samples, id.vars = "timestep")
    post_samples_m[, c("variable", "epidemic_id") := tstrsplit(variable, "_", fixed = TRUE)]
    post_samples_m$epidemic_id #factor(post_samples_m$epidemic_id, levels=unique(post_samples_m$epidemic_id))
    post_samples_m <- post_samples_m %>% mutate(epidemic_id = epid)
    post_samples_merged <- rbind(post_samples_merged,post_samples_m)
  }
}

#write_csv(post_samples_merged, file = paste0('MCMC_output/','merged_',post_size,'_',thinning_steps,'_',burn_in,'.csv'))

demography <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))
contact_matrix_input <- t(t(as.matrix(unname(list_contact_matr[[sel_cntr]])))/pop_age(wpp_age(sel_cntr, 2015), age.limits=c(0,5,20,65))$population)

## PLOTS

print('traces and histograms')

supp.labs <- c("Initial infected", "Reporting rate", "Susceptibility", "Transmissibility")
names(supp.labs) <- c('infected','reporting','sus','trans')

p1 <- ggplot(post_samples_merged[!post_samples_merged$variable == 'll'],
       aes(x = timestep, y = value, col = epidemic_id)) +
  facet_grid(variable~epidemic_id, scales = "free_y") +
  geom_line() + scale_x_continuous(breaks=seq(0,nrow(post_samples),
                                              nrow(post_samples)/4)) +
  labs(title = paste0("Parameter traces, ", sel_cntr, ', ', strain)) + theme_bw() +
  scale_color_viridis(breaks=seq(1,max(post_samples_merged$epidemic_id),1)) +
  labs(color = "Epidemic") + ylab('Value') + xlab('Timestep') 
p1
ggsave(paste0("output/plots/inference/traces/", sel_cntr, '_', strain, '_',
              post_size, "_", thinning_steps, "_", burn_in, "60.png"), 
       width=36,height=24,units="cm")

p2 <- ggplot(post_samples_merged[post_samples_merged$variable == 'll'],
       aes(x = timestep, y = value, col = epidemic_id)) +
  facet_grid(variable~epidemic_id, scales = "free_y") +
  geom_line() + scale_x_continuous(breaks=seq(0,nrow(post_samples),
                                              nrow(post_samples)/4)) +
  labs(title = paste0("Log-likelihood")) + theme_bw() +
  scale_color_viridis(breaks=seq(1,max(post_samples_merged$epidemic_id),1)) +
  labs(color = "Epidemic") + ylab('Log-likelihood') + xlab('Timestep') + ylim(NA, 0)

p3 <- ggplot(post_samples_merged[!post_samples_merged$variable=='ll']) +
  facet_wrap(epidemic_id~variable, scales="free", ncol=8) +
  geom_histogram(aes(x=value,fill=variable), bins=500) +
  labs(title = paste0("Epidemic posterior parameters, ", sel_cntr, ', ', strain)) +
  theme_bw()
p3
ggsave(paste0("output/plots/inference/hist/", sel_cntr, '_', strain, '_',
              post_size, "_", thinning_steps, "_", burn_in, "60.png"), 
       width=36,height=36,units="cm")


p4 <- post_samples_merged %>% filter(!variable %in% c('ll')) %>%
  ggplot() +
  facet_wrap(variable~., scales="free", ncol=1,
             labeller = labeller(variable = supp.labs)) +
  geom_density(aes(x=value,fill=epidemic_id, group=interaction(variable, epidemic_id)),
               alpha=0.2) +
  labs(title = paste0("Epidemic posterior parameters")) +
  theme_bw() + scale_fill_viridis(breaks=seq(1,max(post_samples_merged$epidemic_id),1)) +
  labs(fill = "Epidemic") + ylab('Density') + xlab('') 


## PAIRWISE VARIABLE ANALYSIS

# post_samples_merged %>% pivot_wider(names_from = variable, values_from = value) %>% 
#   ggplot(aes(x = sus, y = trans, group = epidemic_id, col = timestep)) +
#   geom_point(pch=4,size=0.5) + facet_wrap(.~epidemic_id, scales='free') +
#   labs(title = paste0("Parameter traces")) +
#   scale_color_viridis(option='A', direction=-1) + theme_minimal()


## EPIDEMICS USING THESE PARAMETERS:

print('trajectory')

weekly_fit_function <- function(e, sel_cntr){
  fit_pars <- c(unname(unlist(post_samples_merged %>% filter(timestep==max(post_samples_merged$timestep),
                                                             #country==sel_cntr,
                                                             epidemic_id==e,
                                                             !variable=='ll') %>% 
                                select(value))),0,0)
  
  fit_pars[1] <- log(fit_pars[1])
  fit_pars[2] <- 100*fit_pars[2]
  fit_pars[4] <- log10(fit_pars[4])
  
  #fit_pars <- c(-7, 10, 0.5, 3, 0, 0)#c(-8, 10, 0.4, 4, 0, 0)
  
  data_fitting <- df_epid_threshold %>% 
    filter(grepl("NONSENT",metasource) & country %in% sel_cntr & STRAIN %in% strain) %>%
    select(!c(over_peak,flu_included,over_inclusion,positivity,flu_peak,seq))
  epidemics_to_fit <- list()
  for (epid_index_val in unique(data_fitting$epid_index[data_fitting$epidem_inclusion>0])) {
    xx <- data_fitting %>% filter(epid_index %in% epid_index_val)
    epidemics_to_fit[[epid_index_val]] <- list(
      start=min(xx$ISO_WEEKSTARTDATE),
      end = max(xx$ISO_WEEKSTARTDATE),
      initial_params = c(-8, 10, 0.4, 4, 0, 0), # c(-8, 10, 0.4, 4, 0, 0), 
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
  vaccine_calendar <- vaccine_calendar_list[[e]]
  
  fit_pars[2] <- fit_pars[2]/100 
  age.group.limits <- c(5,20,65)
  hemisphere_input <- df_cntr_table[df_cntr_table$country == sel_cntr,]$hemisphere
  
  year_index = as.numeric(year(vaccine_calendar$dates[1]))
  # if in NH, has vaccinations and epidemic starts in the first half of the year,
  # replace with previous year's vaccination program/match
  if(hemisphere_input == 'NH' & 
     sel_cntr %in% c('Canada', 'United Kingdom') & # ignore countries with no vacc
     month(vaccine_calendar$dates[1]) %in% 1:7){
    year_index <- year_index - 1
  }
  coverage = as.numeric(c(unname(cov_data %>% filter(country==sel_cntr,
                                                     year==year_index) %>% select(!country:year))))
  efficacy = unlist(unname(matches %>% filter(year == year_index,
                                              strain_match == strain,
                                              hemisphere == hemisphere_input
  ) %>% select(!hemisphere:match)))
  
  prop_vacc_start <- list(
    prop_vaccine_compartments = rep(coverage, 3), #proportion of individuals vaccinated
    prop_R_vaccinated = rep(efficacy, 3), #proportion of vaccinated individuals in Rv compartment (as opposed to Sv)
    prop_R = rep(0,12) # proportion of individuals in R (unvaccinated)
    # set to 0 as immunity is not assumed to be related to the previous season
  )
  
  fit <- incidence_function_fit(
    demography_input = demography,
    parameters = fit_pars,
    calendar_input = vaccine_calendar,
    contact_ids_sample = NA, 
    contacts = contact_matrix_input,
    waning_rate = 0,
    vaccination_ratio_input = prop_vacc_start,
    begin_date = vaccine_calendar$dates[1], 
    end_date = vaccine_calendar$dates[length(vaccine_calendar$dates)],  
    year_to_run = year(vaccine_calendar$dates[1]), 
    efficacy_now =rep(0,12), 
    efficacy_next=rep(0,12),
    efficacy_next2 =rep(0,12), 
    previous_summary = NA, 
    age_groups_model = age.group.limits
  )
  fit <- data.table(fit)
  weekly_cases <- fit[, sum(.SD, na.rm=TRUE), by="time"]
  weekly_cases$V1 <- weekly_cases$V1*exp(fit_pars[1])
  return(weekly_cases)
}

data_fitting <- df_epid_threshold %>% 
  filter(grepl("NONSENT",metasource) & country %in% sel_cntr & STRAIN %in% strain) %>%
  select(!c(over_peak,flu_included,over_inclusion,positivity,flu_peak,seq))
epidemics_to_fit <- list()
for (epid_index_val in unique(data_fitting$epid_index[data_fitting$epidem_inclusion>0])) {
  xx <- data_fitting %>% filter(epid_index %in% epid_index_val)
  epidemics_to_fit[[epid_index_val]] <- list(
    start=min(xx$ISO_WEEKSTARTDATE),
    end = max(xx$ISO_WEEKSTARTDATE),
    initial_params = c(-7, 9, 0.25, 3, 0, 0), # c(-8, 10, 0.4, 4, 0, 0), 
    data_points = xx$value,
    type = "A"
  )
}

index <- c(); week <- c(); value <- c(); data <- c()
for(e in unique(post_samples_merged$epidemic_id)){
  weekly_fit_data <- weekly_fit_function(e, sel_cntr)
  index <- c(index, rep(e, length(weekly_fit_data$time)))
  week <- c(week, weekly_fit_data$time)
  value <- c(value, weekly_fit_data$V1)
}
fitting <- data.frame(country = rep(sel_cntr,length(index)),
                      index = index,
                      week = week,
                      value = value,
                      data = NA) 
for(i in 1:nrow(fitting)){
  val <- unname(unlist(df_epid_threshold %>% 
                         filter(country==sel_cntr, STRAIN==strain,
                                ISO_WEEKSTARTDATE==as.Date(fitting$week[i]),
                                metasource=='NONSENT') %>% 
                         select(value)))
  if(length(val) > 0){
    fitting$data[i] <- as.numeric(val)
  }
}

# line plots with *last* parameters of fit, not median etc. (just an example)
p5 <- fitting %>% 
  ggplot(aes(x=as.Date(week), y=value, group=index)) +
  geom_line(lwd=1) +
  theme_bw() + facet_wrap(.~index, scales='free') +
  geom_point(aes(x=as.Date(week), y=data, group=index), size=1.5, col=2) +
  geom_line(aes(x=as.Date(week), y=data, group=index), col=2, lwd=0.6) +
  ggtitle(paste0(sel_cntr, ', ', strain)) + 
  theme(axis.text.x = element_text(angle = 90))
p5
ggsave(paste0("output/plots/inference/trajectory/", sel_cntr, '_', strain, '_',
              post_size, "_", thinning_steps, "_", burn_in, "60.png"), 
       width=36,height=24,units="cm")


## Credible intervals
# i.e. simulate 1000 of the epidemics and geom_ribbon the x% HDI of epidemics
# at each time point 

print('CIs')

weekly_fit_function_CI <- function(pars, e, sel_cntr){
  
  pars[1] <- log(pars[1])
  pars[2] <- 100*pars[2]
  pars[4] <- log10(pars[4])
  
  data_fitting <- df_epid_threshold %>% 
    filter(grepl("NONSENT",metasource) & country %in% sel_cntr & STRAIN %in% strain) %>%
    select(!c(over_peak,flu_included,over_inclusion,positivity,flu_peak,seq))
  epidemics_to_fit <- list()
  for (epid_index_val in unique(data_fitting$epid_index[data_fitting$epidem_inclusion>0])) {
    xx <- data_fitting %>% filter(epid_index %in% epid_index_val)
    epidemics_to_fit[[epid_index_val]] <- list(
      start = min(xx$ISO_WEEKSTARTDATE),
      end = max(xx$ISO_WEEKSTARTDATE),
      initial_params = c(-8, 10, 0.4, 4, 0, 0), # c(-8, 10, 0.4, 4, 0, 0), 
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
  vaccine_calendar <- vaccine_calendar_list[[e]]
  
  pars[2] <- pars[2]/100 
  age.group.limits <- c(5,20,65)
  hemisphere_input <- df_cntr_table[df_cntr_table$country == sel_cntr,]$hemisphere
  
  year_index = as.numeric(year(vaccine_calendar$dates[1]))
  # if in NH, has vaccinations and epidemic starts in the first half of the year,
  # replace with previous year's vaccination program/match
  if(hemisphere_input == 'NH' & 
     sel_cntr %in% c('Canada', 'United Kingdom') & # ignore countries with no vacc
     month(vaccine_calendar$dates[1]) %in% 1:7){
    year_index <- year_index - 1
  }
  coverage = as.numeric(c(unname(cov_data %>% filter(country==sel_cntr,
                                                     year==year_index) %>% select(!country:year))))
  efficacy = unlist(unname(matches %>% filter(year == year_index,
                                              strain_match == strain,
                                              hemisphere == hemisphere_input
  ) %>% select(!hemisphere:match)))
  
  prop_vacc_start <- list(
    prop_vaccine_compartments = rep(coverage, 3), #proportion of individuals vaccinated
    prop_R_vaccinated = rep(efficacy, 3), #proportion of vaccinated individuals in Rv compartment (as opposed to Sv)
    prop_R = rep(0,12) # proportion of individuals in R (unvaccinated)
    # set to 0 as immunity is not assumed to be related to the previous season
  )
  
  fit <- incidence_function_fit(
    demography_input = demography,
    parameters = pars,
    calendar_input = vaccine_calendar,
    contact_ids_sample = NA, 
    contacts = contact_matrix_input,
    waning_rate = 0,
    vaccination_ratio_input = prop_vacc_start,
    begin_date = vaccine_calendar$dates[1], 
    end_date = vaccine_calendar$dates[length(vaccine_calendar$dates)],  
    year_to_run = year(vaccine_calendar$dates[1]), 
    efficacy_now =rep(0,12), 
    efficacy_next=rep(0,12),
    efficacy_next2 =rep(0,12), 
    previous_summary = NA, 
    age_groups_model = age.group.limits
  )
  fit <- data.table(fit)
  weekly_cases <- fit[, sum(.SD, na.rm=TRUE), by="time"]
  weekly_cases$V1 <- weekly_cases$V1*exp(pars[1])
  return(weekly_cases)
}

# this takes a few minutes as simulates ~10,000 epidemics
CI95_epis <- data.frame()
track <- data.frame()
for(epid_index in unique(fitting$index)){
  set.seed(seed_to_use)
  sample_epis <- data.frame(timestep = sort(sample(unique(post_samples_merged$timestep), 
                                                   size=1000, replace=F)),
                            reporting = NA,
                            trans = NA,
                            sus = NA,
                            infected = NA)
  sample_epis <- cbind(sample_epis, matrix(NA, nrow=nrow(sample_epis), 
                                           ncol=nrow(fitting[fitting$index==epid_index,])))
  for(i in 1:nrow(sample_epis)){
    sample_epis[i,2:5] <- unname(unlist(post_samples_merged %>% 
                                          filter(timestep == sample_epis$timestep[i],
                                                 epidemic_id == epid_index,
                                                 !variable=='ll') %>% select(value)))
    predic_epi <- weekly_fit_function_CI(unname(unlist(sample_epis[i,2:5])), 
                                         epid_index, sel_cntr)$V1
    sample_epis[i,6:ncol(sample_epis)] <- predic_epi
  }
  CIdf <- data.frame(index = epid_index,
                     lower = rep(NA, nrow(fitting[fitting$index==epid_index,])),
                     upper = rep(NA, nrow(fitting[fitting$index==epid_index,])))
  for(i in 1:nrow(CIdf)){
    CIdf$lower[i] <- unname(unlist(hdi(sample_epis[,i+5], ci = 0.95)[2]))
    CIdf$upper[i] <- unname(unlist(hdi(sample_epis[,i+5], ci = 0.95)[3]))
  }
  CI95_epis <- rbind(CI95_epis, CIdf)
  track_input <- sample_epis %>% select(!c(reporting:infected)) %>% pivot_longer(!timestep) %>% 
    rename(epid = timestep, week = name, cases = value)
  track <- rbind(track, cbind(epid_index, track_input))
} 

names(CI95_epis) <- c('index','lower','upper')
fitting$lower <- CI95_epis$lower
fitting$upper <- CI95_epis$upper

p6 <- fitting %>% 
  ggplot() +
  geom_ribbon(aes(x = as.Date(week), ymin = lower, ymax = upper),
              fill='cornflowerblue', alpha = 0.75) +
  theme_bw() + facet_wrap(.~index, scales='free') + 
  geom_point(aes(x=as.Date(week), y=data, group=index), size=1.5, col=2) +
  geom_line(aes(x=as.Date(week), y=data, group=index), col=2, lwd=0.6) +
  ggtitle(paste0(sel_cntr, ', ', strain)) + 
  theme(axis.text.x = element_text(angle = 90))
p6
ggsave(paste0("output/plots/inference/95CI/", sel_cntr, '_', strain, '_',
              post_size, "_", thinning_steps, "_", burn_in, "60.png"), 
       width=36,height=24,units="cm")

## Plotting 1000 epidemics instead

track <- track %>% rename(index = epid_index) 
track <- track %>% 
  arrange(index, as.numeric(week)) %>% 
  mutate(week_date = rep(fitting$week, each=1000))
track_joint <- track %>% mutate(data = rep(fitting$data, each=1000))

p7 <- ggplot(track_joint, aes(x = as.Date(week_date), y = cases, group=epid)) +
  geom_line(alpha = 0.2) +
  theme_bw() + facet_wrap(.~index, scales='free') + ylim(c(0,NA)) +
  ggtitle(paste0(sel_cntr, ', ', strain)) +
  geom_point(aes(x=as.Date(week_date), y=data, group=epid, col=index), size=1.5) +
  geom_line(aes(x=as.Date(week_date), y=data, group=epid, col=index), lwd=0.6) +
  scale_color_viridis(breaks=seq(1,max(post_samples_merged$epidemic_id),1)) +
  labs(color = "Epidemic") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Date') + ylab('Influenza cases') +
  labs(title = paste0("Epidemic fitting"))
p7
ggsave(paste0("output/plots/inference/95CI/", sel_cntr, '_', strain, '_',
              post_size, "_", thinning_steps, "_", burn_in, "OVERLAY60.png"), 
       width=36,height=24,units="cm")



p4 + (p2 + p7 + plot_layout(ncol = 1, heights = c(1, 2))) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ',
                  title = paste0(sel_cntr, ', ', strain),
                  theme = theme(plot.title = element_text(size = 16))) + 
                    plot_layout(ncol = 2, guides='collect') 
ggsave(paste0("output/plots/inference/patchwork/", sel_cntr, '_', strain, '_',
              post_size, "_", thinning_steps, "_", burn_in, ".png"), 
       width=50,height=40,units="cm")



