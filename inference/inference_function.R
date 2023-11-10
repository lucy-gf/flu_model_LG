### FUNCTION: custom_inference ###
## for running all epidemics in one country at once (joint log-likelihood)
custom_inference <- function(
    input_demography, 
    vaccine_calendar_list, 
    input_polymod, 
    ili = NULL, 
    mon_pop = NULL, 
    # n_pos, 
    epidemics_to_fit, # passing in all epidemics instead of n_pos 
    n_samples, 
    initial,
    mapping,
    nbatch, 
    nburn, 
    blen,
    sel_cntr,
    hemisphere_input
) {

  age.group.limits <- c(5,20,65)  # 4 age groups used in the model
  age.groups <- stratify_by_age(input_demography, age.group.limits)

  # Define combined log likelihood which calls llikelihood for each epidemic 
  # and generates the combined likelihood.
  
  #names(initial_parameters) <- c("reporting", "transmissibility","susceptibility","initial_infected")
  
  out_rate_cent <<- 0
  
  llikelihood_all <- function(pars){
    
    # browser()
    out_num <- 1000
    if(global_index %% out_num == 0){
      print(paste0("Step = ",global_index,', Acceptance = ', 
                   round((out_rate - out_rate_cent)/(out_num/100), digits = 1),
                   '%'))
      out_rate_cent <<- out_rate
      # old_time <<- new_time
      # new_time <<- Sys.time()
      # print(new_time - old_time)
    }
    global_index <<- global_index+1
    #if(global_index==10){browser()}
    total_ll <- 0
    all_pars <- pars
    for (e in 1:length(epidemics_to_fit)){
      pars <- all_pars[(((e-1))*6+1):(e*6)] # select parameter set for epidemic e (currently no shared pars)
      vaccine_calendar <- vaccine_calendar_list[[e]] # select vaccine calendar for epidemic e
      n_pos <- data.frame(time = as.Date(epidemics_to_fit[[e]]$weeks),
                          data = epidemics_to_fit[[e]]$data_points)
      # if("NaN" %in% pars) {browser}
      # if(global_index==11652){browser()}
      ll_epidemic <- llikelihood(pars,n_pos,vaccine_calendar,sel_cntr,hemisphere_input) # generate log-likelihood for epidemic e
      total_ll <- total_ll + ll_epidemic # then add it to the total log-likelihood
      if(total_ll=="-Inf"){break}
    }
    #browser()
    return(total_ll)
  }
  
  
  # Define the actual log likelihood function
  llikelihood <- function(pars, n_pos, vaccine_calendar, sel_cntr,hemisphere_input) {
    
    # contacts <- fluEvidenceSynthesis::contact_matrix(
    #   as.matrix(input_polymod),
    #   input_demography, 
    #   age.group.limits
    # ) # Now using the original Prem et al. matrix, divided by population
    
    # Population size initially infected by age
    initial.infected <- rep(10^pars[4], length(age.groups))
    pars[2] <- pars[2]/100 # transmissibility
    
    # WITH EXISTING VACCINATIONS
    year_index = as.numeric(year(vaccine_calendar$dates[1]))
    # if in NH, has vaccinations and epidemic starts in the first half of the year,
    # replace with previous year's vaccination program/match
    if(hemisphere_input == 'NH' & 
       month(vaccine_calendar$dates[1]) < 7){
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
    
    # Run simulation
    # Note that to reduce complexity 
    # we are using the same susceptibility parameter for multiple age groups
    # browser()
    odes <- incidence_function_fit(
      demography_input = input_demography,
      parameters = pars,
      calendar_input = vaccine_calendar,
      contact_ids_sample = NA, 
      contacts = input_polymod,
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
    # browser()
    odes <- data.table(odes)
    
    # break if it goes too wrong
    if(sum(is.na(odes)) > 0){return(as.numeric("-Inf"))}
    
    weekly_cases <- odes[, sum(.SD, na.rm=TRUE), by="time"]
    weekly_cases <- left_join(weekly_cases, n_pos, by="time")
    
    total_ll <- 0
    for(i in 1:nrow(weekly_cases)){
      if(is.na(weekly_cases$data[i])){weekly_ll <- 0} else{
        if (round(as.numeric(weekly_cases$V1[i])) < weekly_cases$data[i]) {
          #print('weekly < n_pos')
          return(-Inf) #(-1e+10)
        } else
        {
          weekly_ll <-  dbinom(
            x = weekly_cases$data[i],
            size = round(as.numeric(weekly_cases[i,"V1"])),
            prob = exp(pars[1]),
            log = T
          )
          if(is.nan(weekly_ll)) {return(as.numeric("-Inf"))}
        }
      }
      total_ll <- total_ll + weekly_ll  
    }
    #browser()
    return(total_ll)
  }
  
  ### FUNCTION: llprior ### 
  llprior <- function(pars) {

    # 1 is reporting prob, 2/10 is transmission, 3 is susceptibility, 10^4 is initial infections
    r0_gamma_pars <- c(11.082591, 9.248767)
    names(r0_gamma_pars) = c('shape', 'rate')
    sus_beta_pars <- c(50.19925, 32.55043)
    names(sus_beta_pars) <- c('shape1', 'shape2')
    
    for (e in 1:(length(pars)/6)){
      if(
        exp(pars[(e-1)*6+1]) < 0 || # exp(rep) is a probability
        exp(pars[(e-1)*6+1]) > 1 ||
        pars[(e-1)*6+2] < 0 || # transmissibility must be geq 0
        pars[(e-1)*6+3] < 0 || # susceptibility is a probability
        pars[(e-1)*6+3] > 1  ||
        pars[(e-1)*6+4] < -1 || # minimum of 1 infected individual in each age group
        pars[(e-1)*6+4] > log10(min(age.groups)) # demography-specific upper bound for init_inf
      ) {print('Parameter out of prior bounds'); return(-Inf)}
    }

    lprob <- 0

    for(e in 1:(length(pars)/6)){
      R0 <- fluEvidenceSynthesis::as_R0(
        transmission_rate = pars[(e-1)*6+2]/100,
        contact_matrix = input_polymod, # contact_matrix = contacts_matrixformat,
        age_groups = stratify_by_age(input_demography, limits = age.group.limits)
      )
      # prior on R0
      lprob <- lprob + dgamma(R0 - 1, shape = r0_gamma_pars[1], rate = r0_gamma_pars[2], log = T)
      # prior on susceptibility
      lprob <- lprob + unname(dbeta(pars[(e-1)*6+3], shape1 = sus_beta_pars[1], shape2 = sus_beta_pars[2], log = T))
    }
    
    return(lprob)
    #return(0)
  }
  ### END OF FUNCTION: llprior ###
  
  ### should pre-define a dataframe of 'output' for each country with x rows 
  ### otherwise will take a really long time to 'extend' a dataframe
  ### instead of taking a shorter amount of time to fill in an existing one

  # Store the contact ids used during inference
  contact.ids <- list()
  
  out_rate <<- 0
  # Run adaptive.mcmc
  mcmc.result <- adaptive.mcmc(
    lprior = llprior, 
    llikelihood = llikelihood_all, #llikelihood, 
    nburn = nburn, 
    initial = initial,
    nbatch = nbatch, 
    blen = blen,
    acceptfun = function() {
      out_rate <<- out_rate + 1
    }
  )

  #mcmc.result$contact.ids <- t(data.frame(contact.ids))
  mcmc.result
}
### END OF FUNCTION: custom_inference ##


### FUNCTION: custom_inference_sep ###
## for running each epidemic individually (fits quicker)
custom_inference_sep <- function(
    input_demography, 
    vaccine_calendar, 
    input_polymod, 
    ili = NULL, 
    mon_pop = NULL, 
    # n_pos, 
    epidemics_to_fit, # passing in all epidemics but will choose one 
    n_samples, 
    initial,
    mapping,
    nbatch, 
    nburn, 
    blen,
    sel_cntr,
    hemisphere_input
) {
  
  global_index <<- 1
  age.group.limits <- c(5,20,65)  # 4 age groups used in the model
  age.groups <- stratify_by_age(input_demography, age.group.limits)
  
  #names(initial_parameters) <- c("reporting", "transmissibility","susceptibility","initial_infected")
  
  out_rate_cent <<- 0
  
  llikelihood_print <- function(pars){
    
    out_num <- 1000
    if(global_index %% out_num == 0){
      print(paste0("Step = ",global_index,', Acceptance = ', 
                   round((out_rate - out_rate_cent)/(out_num/100), digits = 1),
                   '% (Epidemic ', year(vaccine_calendar$dates[1]), ')'))
      out_rate_cent <<- out_rate
      # old_time <<- new_time
      # new_time <<- Sys.time()
      # print(new_time - old_time)
    }
    global_index <<- global_index+1
  
    ll_epidemic <- 0
    n_pos <- data.frame(time = as.Date(epidemics_to_fit$weeks),
                        data = epidemics_to_fit$data_points)

    ll_epidemic <- llikelihood(pars,n_pos,vaccine_calendar,sel_cntr,hemisphere_input) 
    #browser()
    return(ll_epidemic)
  }
  
  
  # Define the actual log likelihood function
  llikelihood <- function(pars, n_pos, vaccine_calendar, sel_cntr, hemisphere_input) {
    
    # Population size initially infected by age
    initial.infected <- rep(10^pars[4], length(age.groups))
    pars[2] <- pars[2]/100 # transmissibility
    
    # WITH EXISTING VACCINATIONS
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
    
    # Run simulation
    # Note that to reduce complexity 
    # we are using the same susceptibility parameter for multiple age groups
    # browser()
    odes <- incidence_function_fit(
      demography_input = input_demography,
      parameters = pars,
      calendar_input = vaccine_calendar,
      contact_ids_sample = NA, 
      contacts = input_polymod,
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
    # browser()
    odes <- data.table(odes)
    
    # break if it goes too wrong
    if(sum(is.na(odes)) > 0){return(as.numeric("-Inf"))}
    
    weekly_cases <- odes[, sum(.SD, na.rm=TRUE), by="time"]
    weekly_cases <- left_join(weekly_cases, n_pos, by="time")
    
    total_ll <- 0
    for(i in 1:nrow(weekly_cases)){
      if(is.na(weekly_cases$data[i])){weekly_ll <- 0} else{
      if (round(as.numeric(weekly_cases$V1[i])) < weekly_cases$data[i]) {
        #print('weekly < n_pos')
        return(-Inf) #(-1e+10)
      } else
      {
        weekly_ll <-  dbinom(
          x = weekly_cases$data[i],
          size = round(as.numeric(weekly_cases[i,"V1"])),
          prob = exp(pars[1]),
          log = T
        )
        if(is.nan(weekly_ll)) {return(as.numeric("-Inf"))}
      }
      }
      total_ll <- total_ll + weekly_ll  
    }
    #browser()
    return(total_ll)
  }
  
  ### FUNCTION: llprior ### 
  llprior <- function(pars) {
    
    # 1 is reporting prob, 2/10 is transmission, 3 is susceptibility, 10^4 is initial infections
    r0_gamma_pars <- c(11.082591, 9.248767)
    names(r0_gamma_pars) = c('shape', 'rate')
    sus_beta_pars <- c(50.19925, 32.55043)
    names(sus_beta_pars) <- c('shape1', 'shape2')
    
    if(
      exp(pars[1]) < 0 || # exp(rep) is a probability
      exp(pars[1]) > 1 ||
      pars[2] < 0 || # transmissibility must be geq 0
      pars[3] < 0 || # susceptibility is a probability
      pars[3] > 1  ||
      pars[4] < -1 || # minimum of 1 infected individual in each age group
      pars[4] > log10(min(age.groups)) # demography-specific upper bound for init_inf
    ) {print('Parameter out of prior bounds'); return(-Inf)}
    
    lprob <- 0
    
    R0 <- fluEvidenceSynthesis::as_R0(
        transmission_rate = pars[2]/100,
        contact_matrix = input_polymod, # contact_matrix = contacts_matrixformat,
        age_groups = stratify_by_age(input_demography, limits = age.group.limits)
    )
    # prior on R0
    lprob <- lprob + dgamma(R0 - 1, shape = r0_gamma_pars[1], rate = r0_gamma_pars[2], log = T)
    # prior on susceptibility
    lprob <- lprob + unname(dbeta(pars[3], shape1 = sus_beta_pars[1], shape2 = sus_beta_pars[2], log = T))
    
    return(lprob)
    #return(0)
  }
  ### END OF FUNCTION: llprior ###
  
  # Store the contact ids used during inference
  contact.ids <- list()
  
  out_rate <<- 0
  # Run adaptive.mcmc
  mcmc.result <- adaptive.mcmc(
    lprior = llprior, 
    llikelihood = llikelihood_print, 
    nburn = nburn, 
    initial = initial,
    nbatch = nbatch, 
    blen = blen,
    acceptfun = function() {
      out_rate <<- out_rate + 1
    }
  )
  
  #mcmc.result$contact.ids <- t(data.frame(contact.ids))
  mcmc.result
}
### END OF FUNCTION: custom_inference_sep ##



### FUNCTION: incidence_function_fit ###
incidence_function_fit <- function(
    demography_input,
    parameters,
    calendar_input,
    contact_ids_sample,
    contacts,
    waning_rate,
    vaccination_ratio_input,
    begin_date, 
    end_date,  
    year_to_run, 
    efficacy_now, 
    efficacy_next,
    efficacy_next2, 
    previous_summary, 
    age_groups_model
){
  
  risk_ratios_input <- matrix(c(rep(0,8)), ncol = 4 , byrow = T) # Not using risk groups so setting this here for now.
  age_group_sizes <- stratify_by_age(demography_input, age_groups_model)
  population_stratified <- stratify_by_risk(age_group_sizes, risk_ratios_input)
 
  # define model timings
  interval = 7
  t <- as.numeric(seq(begin_date, end_date, interval))
  # define age group inputs
  no_groups <- length(population_stratified)
  
  initial_infected <- rep(10^parameters[4],length(age_group_sizes)) #6)
  initial_infected <- stratify_by_risk(initial_infected, risk_ratios_input)
  susceptibility <- c((0.2*1 + 0.8*parameters[3]), rep(parameters[3],3)) # 1 in [0,1), parameters[3] in [1,\infty)
  transmissibility = parameters[2]
  
  initial_vaccinated_prop = unlist(vaccination_ratio_input[[1]])
  initial_Rv_prop = unlist(vaccination_ratio_input[[2]])
  initial_R_prop = unlist(vaccination_ratio_input[[3]])
  #Assume that all R become susceptible again at the start of each posterior
  initial_R_prop <- rep(0,no_groups) ### these are the same now but may not be later, check which is correct
  
  # Contacts matrix only covers one set of age groups, here we "repeat" it to also cover 
  # risk groups
  # new_cij <- matrix(rep(0,18*18), nrow = 18)
  # for (k in 1:3) {
  #   for (l in 1:3) {
  #     lk <- (k - 1)*6 + 1
  #     ll <- (l - 1)*6 + 1
  #     new_cij[lk:(lk + 6 - 1), ll:(ll + 6 - 1)] <- contacts_matrixformat
  #   }
  # }
  new_cij <- matrix(rep(0,12*12), nrow = 12)
  for (k in 1:3) {
    for (l in 1:3) {
      lk <- (k - 1)*4 + 1
      ll <- (l - 1)*4 + 1
      new_cij[lk:(lk + 4 - 1), ll:(ll + 4 - 1)] <- contacts
    }
  }
  
  # specify the model

  mod <- gen_seeiir_ag_vacc_waning$new(
    no_groups = no_groups,
    cij = new_cij,
    trans = transmissibility,
    pop = population_stratified,
    I0 = initial_infected,
    V0 = initial_vaccinated_prop,
    R0 = initial_R_prop,
    RV0 = initial_Rv_prop,
    susc = rep(susceptibility,3),
    alpha = calendar_input$efficacy[1:no_groups],
    omega = waning_rate,
    dates = calendar_input$dates,
    calendar = matrix(calendar_input$calendar, ncol = 4*3),
    gamma1 = 2/infection_delays[1],
    gamma2 = 2/infection_delays[2], 
    num_vac_start = rep(0,no_groups) # don't need to be tracked
  )
  
  # run the model
  y <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  # calculate the cumulative values
  y <- mod$transform_variables(y)$cumI
  # Returning the differences in cumulative infections from one timestep to the other
  y <- data.table(y[2:(nrow(y)), ] - y[1:(nrow(y) - 1), ])
  # add the dates
  y$time <- seq(begin_date,length.out = nrow(y), by = interval)
  y<-data.table(y)
  # reformat the dates

  return(y)
}
### END OF FUNCTION: incidence_function_fit ###








