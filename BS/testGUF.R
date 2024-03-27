
scenario_name <- c('none', 'base', 'low_cov', 'rel_inf')[2]
c_number <- 5
c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
            "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
            "Southern America")[c_number]
c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]

annual_country <- data.frame()
weekly_agg <- data.frame()

# load no vacc cases
itz_cases_no_vacc <- readRDS(paste0("data/vacc_output/vacc_", c_code, '_none_ct_1.rds'))[[1]]
large_list <- readRDS(paste0("data/vacc_output/vacc_", c_code, '_', scenario_name, '_ct_4.rds'))

country_code_index <- 'GTM'

guf_inc <- itz_cases_no_vacc %>% filter(country_code==country_code_index) %>% mutate(scenario = 'no_vacc')
for(i in 1:5){
  tester <- data.frame(large_list[i]) %>% 
      rename(country=paste0('vt_',i,'_ct_4.country'),country_code=paste0('vt_',i,'_ct_4.country_code'),
             simulation_index=paste0('vt_',i,'_ct_4.simulation_index'),week=paste0('vt_',i,'_ct_4.week'),
             IU1A=paste0('vt_',i,'_ct_4.IU1A'),IU2A=paste0('vt_',i,'_ct_4.IU2A'),
             IU3A=paste0('vt_',i,'_ct_4.IU3A'),IU4A=paste0('vt_',i,'_ct_4.IU4A'),
             IV1A=paste0('vt_',i,'_ct_4.IV1A'),IV2A=paste0('vt_',i,'_ct_4.IV2A'),
             IV3A=paste0('vt_',i,'_ct_4.IV3A'),IV4A=paste0('vt_',i,'_ct_4.IV4A'),
             IU1B=paste0('vt_',i,'_ct_4.IU1B'),IU2B=paste0('vt_',i,'_ct_4.IU2B'),
             IU3B=paste0('vt_',i,'_ct_4.IU3B'),IU4B=paste0('vt_',i,'_ct_4.IU4B'),
             IV1B=paste0('vt_',i,'_ct_4.IV1B'),IV2B=paste0('vt_',i,'_ct_4.IV2B'),
             IV3B=paste0('vt_',i,'_ct_4.IV3B'),IV4B=paste0('vt_',i,'_ct_4.IV4B')) %>% 
      filter(country_code==country_code_index) %>% 
      mutate(scenario = names(large_list)[i]) 
    guf_inc <- rbind(guf_inc, tester)
}

guf_inc <- guf_inc %>% mutate(vt = substr(scenario, 4,4), ct = substr(scenario, 9,9))
guf_inc$tot <- rowSums(guf_inc[,5:20])

guf_inc %>% filter(simulation_index==1, year(week) > 2050) %>% ggplot() +
  geom_line(aes(x=week, y=tot, group=vt, color=vt)) +
  theme_minimal() +
  scale_color_manual(values=vt_colors)

guf_inc <- guf_inc %>% group_by(simulation_index, scenario) %>% 
  mutate(cum=cumsum(tot)) %>% ungroup()

sims <- 1

onedrive <- guf_inc %>% filter(simulation_index == sims) %>% ggplot() +
  geom_line(aes(x=week, y=cum, group=vt, color=vt), lwd=0.7) +
  theme_minimal() +
  scale_color_manual(values=vt_colors)

length_simulation <- 30
n_simulations <- 100
start_year <- 2025 
years <- 30
sampled_epidemics <- rbind(read_csv(paste0("data_for_BS/sampled_epidemics_", length_simulation,
                                           "_", n_simulations,"_",c_code, ".csv")))

cases_30_100 <- data.frame()
k <- 4 # age-specific targeting strategy, in 1:5
targeting <- paste0('ct_', k)
scenario_name <- c('none', 'base', 'low_cov', 'rel_inf')[2]
if(scenario_name == 'none'){
  scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_', scenario_name)]]
}else{
  scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_', scenario_name)]][1:5 + (k-1)*5]
}
for(vacc_prog_num in 1:5){
  cases_30_100 <- rbind(cases_30_100,
                        fcn_simulate_epidemics(country_code_input = country_code_index,
                                               hemisphere = hemisphere,
                                               sampled_epidemics_input =
                                                 sampled_epidemics[sampled_epidemics$simulation_index == sims,],
                                               vacc_prog = scenarios[[vacc_prog_num]]) %>% mutate(vt = vacc_prog_num))
  print(vacc_prog_num)
}
scenario_name <- c('none', 'base', 'low_cov', 'rel_inf')[1]
if(scenario_name == 'none'){
  scenarios <- vaccine_programs_merged[[paste0('vaccine_programs_', scenario_name)]]
}
vacc_prog_num <- 1
cases_30_100 <- rbind(cases_30_100,
                      fcn_simulate_epidemics(country_code_input = country_code_index,
                                             hemisphere = hemisphere,
                                             sampled_epidemics_input =
                                               sampled_epidemics[sampled_epidemics$simulation_index == sims,],
                                             vacc_prog = scenarios[[vacc_prog_num]]) %>% mutate(vt = 0))

cases_30_100$tot <- rowSums(cases_30_100[,2:17])
cases_30_100 <- cases_30_100 %>% group_by(vt) %>% 
  mutate(cum = cumsum(tot)) %>% ungroup()

mycom <- cases_30_100 %>% ggplot() +
  geom_line(aes(x=week, y=cum, group=as.factor(vt), color=as.factor(vt)), lwd=0.7) +
  theme_minimal() +
  scale_color_manual(values=vt_colors)
  
onedrive + mycom

similar <- guf_inc %>% filter(simulation_index == sims) %>% 
  arrange(vt, week) %>% 
  select(!c(country, country_code, simulation_index, scenario, ct))

similar %>% select(week, vt, tot, cum) %>% 
  mutate(tot_diff = cases_30_100$tot - tot,cum_diff = cases_30_100$cum - cum) %>% 
  filter(vt==1) %>% 
  ggplot() +
  geom_line(aes(x=week, y=tot, group=vt, col=2)) +
  geom_line(data = cases_30_100[cases_30_100$vt==1,],
            aes(x=week, y=tot, group=as.factor(vt)), col=3, lty=2) +
  theme_minimal()

similar %>% select(week, vt, tot, cum) %>% mutate(
  totcom = cases_30_100$tot,
  cumcom = cases_30_100$cum,
  tot_diff = cases_30_100$tot - tot,
  cum_diff = cases_30_100$cum - cum) %>% 
  #filter(!vt=='v') %>% 
  ggplot() +
  geom_line(aes(x=week, y=tot_diff, group=vt, col=vt)) +
  theme_minimal()

similar %>% select(week, vt, tot, cum) %>% 
  mutate(tot_diff = cases_30_100$tot - tot,cum_diff = cases_30_100$cum - cum) %>% 
  filter(vt=='v') %>% 
  ggplot() +
  geom_line(aes(x=week, y=tot, group=vt), lwd=1, col=2) +
  geom_line(data = cases_30_100[cases_30_100$vt==0,],
            aes(x=week, y=tot, group=as.factor(vt)), col=3, lty=2, lwd=0.5) +
  theme_minimal()








