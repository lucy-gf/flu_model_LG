
# READING IN AND VIEWING RESULTS FROM THE COMMAND LINE

strain <- 'INF_A'
post_size <- 10000 #5000 #3000 
thinning_steps <- 30 #100
burn_in <- 200000 #100000
sel_cntr <- 'Turkey'

# DATAFRAME   

post_samples_merged <- data.table()

for(i in 1:10){
  epid <- i
  results_RDS <- readRDS(paste0("command_line_runs/", strain, '_sep/', sel_cntr,"_epid",
                                epid, "_",
                                post_size, "_", thinning_steps, "_", burn_in, ".rds"))
  post_samples <- data.table(results_RDS$batch)
  # colnames(post_samples) <- rep(c("reporting","trans","sus","infected","blank","blank"),length(epidemics_to_fit))
  colnames(post_samples) <- paste0(c("reporting","trans","sus","infected","blank","blank"))
  # remove the blank ones
  cols_to_delete <- which(names(post_samples) %like% "blank")
  post_samples[, (cols_to_delete) := NULL ]
  # add the likelihoods
  post_samples[, ll := results_RDS$llikelihoods]
  post_samples$timestep <- 1:nrow(post_samples)
  post_samples_m <- melt.data.table(post_samples, id.vars = "timestep")
  post_samples_m[, c("variable", "epidemic_id") := tstrsplit(variable, "_", fixed = TRUE)]
  post_samples_m$epidemic_id #factor(post_samples_m$epidemic_id, levels=unique(post_samples_m$epidemic_id))
  post_samples_m <- post_samples_m %>% mutate(epidemic_id = epid)
  post_samples_merged <- rbind(post_samples_merged,post_samples_m)
}

#write_csv(post_samples_merged, file = paste0('MCMC_output/','merged_',post_size,'_',thinning_steps,'_',burn_in,'.csv'))

saved_vals <- post_samples_merged %>% filter(timestep==max(post_samples_merged$timestep)) %>% 
  select(epidemic_id, variable, value)

write_csv(saved_vals, file='saved_vals.csv')


## PLOTS

ggplot(post_samples_merged,
       aes(x = timestep, y = value, col = epidemic_id)) +
  facet_grid(variable~epidemic_id, scales = "free_y") +
  geom_line() + scale_x_continuous(breaks=seq(0,nrow(post_samples),
                                              nrow(post_samples)/4)) +
  labs(title = paste0("Parameter traces")) + theme_bw() 

ggplot(post_samples_merged[post_samples_merged$variable == 'll'],
       aes(x = timestep, y = value, col = epidemic_id)) +
  facet_grid(variable~epidemic_id, scales = "free_y") +
  geom_line() + scale_x_continuous(breaks=seq(0,nrow(post_samples),
                                              nrow(post_samples)/4)) +
  labs(title = paste0("Parameter traces")) + theme_bw() 

ggplot(post_samples_merged[!post_samples_merged$variable=='ll']) +
  facet_grid(epidemic_id~variable, scales="free" ) +
  geom_histogram(aes(x=value,fill=variable), bins=500) +
  labs(title = paste0("Epidemic posterior parameters")) +
  theme_bw()


## PAIRWISE VARIABLE ANALYSIS

post_samples_merged %>% pivot_wider(names_from = variable, values_from = value) %>% 
  ggplot(aes(x = trans, y = sus, group = epidemic_id, col = timestep)) +
  geom_point(pch=4,size=0.5) + facet_wrap(.~epidemic_id) +
  labs(title = paste0("Parameter traces")) + theme_bw() +
  scale_color_viridis(option='A', direction=-1) + theme_minimal()

post_samples_merged %>% pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(prod = sus*trans) %>% 
  ggplot(aes(x=timestep, y=prod, group=epidemic_id, col=epidemic_id)) +
  geom_point(size=0.7) + theme_minimal() + facet_wrap(.~epidemic_id) +
  scale_y_continuous(limits=c(0,5))
  

## EPIDEMICS USING THESE PARAMETERS:

weekly_fit_function <- function(e, sel_cntr){
  fit_pars <- c(unname(unlist(post_samples_merged %>% filter(timestep==max(post_samples_merged$timestep),
                                                             #country==sel_cntr,
                                                             epidemic_id==e,
                                                             !variable=='ll') %>% 
                                select(value))),0,0)
  
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
  
  year_index = as.numeric(year(vaccine_calendar$dates[1]))
  # if in NH, has vaccinations and epidemic starts in the first half of the year,
  # replace with previous year's vaccination program/match
  if(hemisphere_input == 'NH' & 
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
    demography_input = unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5))),
    parameters = fit_pars,
    calendar_input = vaccine_calendar,
    contact_ids_sample = NA, 
    contacts = t(t(as.matrix(unname(list_contact_matr[[sel_cntr]])))/pop_age(wpp_age(sel_cntr, 2015), age.limits=c(0,5,20,65))$population),
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
for(e in 1:length(epidemics_to_fit)){
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
    fitting$data[i] <- val
  }
}

fitting %>% 
  ggplot(aes(x=week, y=value, group=index)) +
  geom_line(lwd=1) +
  theme_minimal() + facet_wrap(.~index, scales='free') + 
  geom_point(aes(x=week, y=data, group=index), size=1.5, col=2) +
  geom_line(aes(x=week, y=data, group=index), col=2, lwd=0.6) +
  ggtitle(sel_cntr) + 
  theme(axis.text.x = element_text(angle = 90))



















