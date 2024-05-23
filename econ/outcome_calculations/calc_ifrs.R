
## CALCULATING ALL NECESSARY IFRS
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_colors.R")
source("econ/outcome_calculations/fcn_ifr.R")

n_simulations <- 100

### Exemplars & Brazil
ifr_method <- c('exemplar','whole_itz','brazil')[1]
## load data
print('loading data')
cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                        "_2010_2019_",n_simulations,".rds"))){
    print(cntr)
    cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                                                          "_2010_2019_",n_simulations,".rds"))[[1]]))
  }
}
ifr_method <- c('exemplar','whole_itz','brazil')[3]
cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_BRA_2010_2019_",
                                                      n_simulations,".rds"))[[1]]))

## turning into national mean age-specific annual number of infections
cases_m <- melt(cases_dt, id.vars=c('country','country_code','simulation_index','week'))
cases_m[, age_grp:=as.numeric(substr(variable,3,3))][, year := year(week)][, c('country','week','variable'):=NULL]
cases_m <- cases_m[year<2016]
cases_m <- cases_m[,lapply(.SD,sum), by=c('country_code','simulation_index','age_grp','year')]
cases_m <- cases_m[,lapply(.SD,mean), by=c('country_code','simulation_index','age_grp')][, year:=NULL]

exemplar_ifrs <- cases_m[,c('country_code','simulation_index','age_grp')]

## mortality rates
cdc_data_l <- data.table(read_csv('econ/outcome_calculations/data/extracted_data/ALL_Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv', show_col_types=F))

print('fitting exemplars')
for(country_code_index in unique(exemplar_ifrs$country_code)){
  morts <-  data.table(
    age_grp = c('<65','65-75','75+'),
    med = cdc_data_l[iso3c == country_code_index & meas=='median',]$value,
    l95 = cdc_data_l[iso3c == country_code_index & meas=='L95',]$value,
    u95 = cdc_data_l[iso3c == country_code_index & meas=='U95',]$value)
  print(country_code_index)
  ifrs_i <- fcn_ifr(
    country_code = country_code_index,
    mortality_rates = morts,
    mortality_probs = c(0.5, 0.025, 0.975), 
    incidence = cases_m[country_code == country_code_index], 
    incidence_ages = c(0,5,20,65))
  
  exemplar_ifrs[country_code == country_code_index, ifr := ifrs_i$ifr]
}

# ggplot(exemplar_ifrs) + 
#   geom_boxplot(aes(x=country_code, y=100000*ifr, fill=country_code)) + 
#   facet_grid(age_grp~., scales='free') + 
#   scale_y_log10() +
#   theme_bw() + ylab('Deaths per 100,000 infections') + 
#   xlab('Country')


### Whole ITZs
ifr_method <- c('exemplar','whole_itz','brazil')[2]
## load data
print('loading whole itzs')
cases_dt <- data.table()
for(cntr in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  if(file.exists(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                        "_2010_2019_",n_simulations,".rds"))){
    print(cntr)
    cases_dt <- rbind(cases_dt, data.table(readRDS(paste0("data/original_epids_output/", ifr_method, "_", cntr, 
                                                          "_2010_2019_",n_simulations,".rds"))[[1]]))
  }
}

## turning into national mean age-specific annual number of infections
cases_dt15 <- cases_dt[year(week)<2016,]
cases_m <- melt(cases_dt15, id.vars=c('country','country_code','simulation_index','week'))
cases_m[, age_grp:=as.numeric(substr(variable,3,3))][, year := year(week)][, c('country','week','variable'):=NULL]
cases_m <- cases_m[,lapply(.SD,sum), by=c('country_code','simulation_index','age_grp','year')]
cases_m <- cases_m[,lapply(.SD,mean), by=c('country_code','simulation_index','age_grp')][, year:=NULL]

## mortality rates
cdc_data_l <- data.table(read_csv('econ/outcome_calculations/data/extracted_data/ALL_Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv', show_col_types=F))

## ifrs
itz_ifrs <- cases_m[,c('country_code','simulation_index','age_grp')]
print('fitting whole itzs')
for(country_code_index in unique(itz_ifrs$country_code)){
  morts <-  data.table(
    age_grp = c('<65','65-75','75+'),
    med = cdc_data_l[iso3c == country_code_index & meas=='median',]$value,
    l95 = cdc_data_l[iso3c == country_code_index & meas=='L95',]$value,
    u95 = cdc_data_l[iso3c == country_code_index & meas=='U95',]$value)
  
  if(country_code_index == 'MAC'){
    morts <- data.table(
      age_grp = c('<65','65-75','75+'),
      med = cdc_data_l[iso3c == 'CHN' & meas=='median',]$value,
      l95 = cdc_data_l[iso3c == 'CHN' & meas=='L95',]$value,
      u95 = cdc_data_l[iso3c == 'CHN' & meas=='U95',]$value)
  }
  if(country_code_index == 'PSE'){
    morts <- data.table(
      age_grp = c('<65','65-75','75+'),
      med = cdc_data_l[iso3c == 'JOR' & meas=='median',]$value,
      l95 = cdc_data_l[iso3c == 'JOR' & meas=='L95',]$value,
      u95 = cdc_data_l[iso3c == 'JOR' & meas=='U95',]$value)
  }
  
  print(country_code_index)
  ifrs_i <- fcn_ifr(
    country_code = country_code_index,
    mortality_rates = morts,
    mortality_probs = c(0.5, 0.025, 0.975), 
    incidence = cases_m[country_code == country_code_index], 
    incidence_ages = c(0,5,20,65))
  
  itz_ifrs[country_code == country_code_index, ifr := ifrs_i$ifr]
}

itz_ifrs[, iso3c := country_code]
itz_ifrs <- itz_ifrs[clusters, on=c('iso3c'), itz := i.cluster_name]

# ggplot(itz_ifrs) + 
#   geom_boxplot(aes(x=country_code, y=100000*ifr, fill=itz)) + 
#   facet_grid(age_grp~itz, scales='free') +
#   scale_y_log10() +
#   theme_bw() + ylab('Deaths per 100,000 infections') + 
#   xlab('Country') 

## maps

national_ifrs_vals <- rbind(exemplar_ifrs[!country_code %in% c('CHN','GHA','TUR')], itz_ifrs[,1:4])
national_ifrs_vals[, attach_code := country_code]
national_ifrs <- clusters[,c('codes','cluster_name')]
setnames(national_ifrs, 'codes','country_code')
national_ifrs[cluster_name %in% c('Africa', 'Asia-Europe', 'Eastern and Southern Asia'), attach_code := country_code]
national_ifrs[cluster_name %in% c('Oceania-Melanesia-Polynesia'), attach_code := 'AUS']
national_ifrs[cluster_name %in% c('Southern America'), attach_code := 'BRA']
national_ifrs[cluster_name %in% c('Northern America'), attach_code := 'CAN']
national_ifrs[cluster_name %in% c('Europe'), attach_code := 'GBR']
national_ifrs[country_code %in% c('ARG'), attach_code := 'ARG']
national_ifrs <- national_ifrs[, lapply(.SD, rep, 400)]
national_ifrs[, age_grp := rep(rep(1:4, each=186), 100)][, simulation_index := rep(1:100, each=4*186)]
national_ifrs <- national_ifrs[national_ifrs_vals[,!'country_code'], on=c('attach_code', 'simulation_index','age_grp')]

write_csv(national_ifrs, file='econ/outcome_calculations/data/national_ifrs.csv')
national_ifrs <- data.table(read_csv('econ/outcome_calculations/data/national_ifrs.csv',
                                show_col_types=F))

national_med_ifrs <- national_ifrs[,c('country_code','age_grp','ifr')][, lapply(.SD, median), by=c('country_code','age_grp')]
national_med_ifrs[, age_grp1 := paste0('age', age_grp)]
national_med_ifrs_w <- dcast(national_med_ifrs[,c('country_code','age_grp1','ifr')], country_code ~ age_grp1, value.var = 'ifr')

require(maps)
world_map <- data.table(map_data("world") %>% rename(country = region) %>%
                          mutate(country_code = countrycode(country, origin = 'country.name', destination = 'iso3c')))

world_map <- world_map[national_med_ifrs_w, on=c('country_code')]
world_map[,c('order','country','subregion'):=NULL]
world_map_l <- melt(world_map, id.vars = c('long','lat','group','country_code'))

supp.labs.agegrps <- c('0-4','5-19','20-64','65+')
names(supp.labs.agegrps) <- c('age1','age2','age3','age4')

ggplot(world_map_l, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = 100000*value),
               col = 'black', lwd=0.4) +
  theme_bw() +
  facet_wrap(variable~., nrow=2,
             labeller = labeller(variable = supp.labs.agegrps)) +
  # scale_fill_viridis(option='A',direction=-1, breaks=300*0:4) +
  scale_fill_gradientn(trans = 'log',
                       colors = viridis(15, direction=-1, option='A'),
                       breaks = c(0.1,1,10,100,1000, 2000)) +
  xlab('') + ylab('') +
  labs(fill='IFR \nper 100,000') +
  theme(text=element_text(size=14))
ggsave(paste0("econ/outcome_calculations/maps/age_spec_ifrs.png"),
       width=50,height=28,units="cm")

ggplot(world_map_l[variable=='age4'], aes(long, lat, group = group)) +
  geom_polygon(aes(fill = 100000*value),
               col = 'black', lwd=0.4) +
  theme_bw() +
  facet_wrap(variable~., nrow=2,
             labeller = labeller(variable = supp.labs.agegrps)) +
  # scale_fill_viridis(option='A',direction=-1, breaks=300*0:4) +
  scale_fill_gradientn(colors = viridis(15, direction=-1, option='A'),
                       breaks = 2*(seq(0,600,by=100))) +
  xlab('') + ylab('') +
  labs(fill='IFR \nper 100,000') +
  theme(text=element_text(size=14))
ggsave(paste0("econ/outcome_calculations/maps/age_spec_ifrs_o65.png"),
       width=25,height=14,units="cm")

cdc_data <- data.table(read_csv('econ/outcome_calculations/data/extracted_data/Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv'))
cdc_data[, country_code := iso3c]
world_map <- world_map[cdc_data[,c('country_code','cdc_mort_u65_med','cdc_mort_6574_med','cdc_mort_o75_med')],
                       on=c('country_code')]

world_map_l2 <- melt(world_map, id.vars = c('long','lat','group','country_code'))

supp.labs.agegrps <- c('<65','65-75','75+')
names(supp.labs.agegrps) <- c('cdc_mort_u65_med','cdc_mort_6574_med','cdc_mort_o75_med')

ggplot(world_map_l2[grepl('cdc',variable)], aes(long, lat, group = group)) +
  geom_polygon(aes(fill = value),
               col = 'black', lwd=0.4) +
  theme_bw() +
  facet_wrap(variable~., nrow=2,
             labeller = labeller(variable = supp.labs.agegrps)) +
  # scale_fill_viridis(option='A',direction=-1, breaks=300*0:4) +
  # scale_fill_gradientn(trans = 'log',
  #                      colors = viridis(15, direction=-1, option='A'),
  #                      breaks = c(0.1,1,10,100,1000, 2000)) +
  scale_fill_gradientn(trans = 'log',
                       colors = viridis(15, direction=-1, option='A'),
                       breaks = c(0.1,1,10,100,1000, 2000)) +
  xlab('') + ylab('') +
  labs(fill='Influenza-associated \nmortality rate \nper 100,000') +
  theme(text=element_text(size=14))
ggsave(paste0("econ/outcome_calculations/maps/age_spec_cdc_mort.png"),
       width=50,height=28,units="cm")


## age-aggregated national ifrs

national_ifrs <- data.table(read_csv('econ/outcome_calculations/data/national_ifrs.csv',
                                     show_col_types=F))
pop_hist_WPP_data <- data.table(pop_hist_WPP_data)
agg_ifrs <- national_ifrs[,!'attach_code'][, lapply(.SD, median), by=c('country_code','cluster_name','age_grp')]

for(i in 1:nrow(agg_ifrs)){
  country_names <- unlist(unname(clusters[clusters$codes==agg_ifrs$country_code[i], 
                                          c('country','country_altern','country_altern_2')]))
  agg_ifrs[i, name := country_names[1]]
  if(agg_ifrs$country_code[i]=='TUR'){country_names <- c(country_names, 'Turkiye')}
  if(agg_ifrs$age_grp[i] == 1){
    agg_ifrs[i, pop_prop := sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:4])))/
               sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))]
  }
  if(agg_ifrs$age_grp[i] == 2){
    agg_ifrs[i, pop_prop := sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 5:7])))/
               sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))]
  }
  if(agg_ifrs$age_grp[i] == 3){
    agg_ifrs[i, pop_prop := sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 8:16])))/
               sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))]
  }
  if(agg_ifrs$age_grp[i] == 4){
    agg_ifrs[i, pop_prop := sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 17:24])))/
               sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))]
  }
  if(i %% 100 == 0){print(paste0(round(100*i/nrow(agg_ifrs),2), '%'))}
}

agg_ifrs[, weighted_ifr := ifr*pop_prop]
agg_ifrs_sum <- agg_ifrs[,c('name','country_code','cluster_name','weighted_ifr')][, lapply(.SD, sum), 
                                                                           by=c('name','country_code','cluster_name')]
setnames(agg_ifrs_sum, 'weighted_ifr','aggr_ifr')

ggplot(agg_ifrs_sum) + 
  geom_histogram(aes(x=aggr_ifr*100000, fill=cluster_name), bins=100) +
  # facet_wrap(cluster_name~.) +
  theme_bw() + xlab('Age-aggregated IFR, per 100,000 infections') +
  labs(fill='ITZ') + 
  scale_fill_manual(values = cluster_colors2)

agg_ifrs_sum[, ifr_per_100000 := aggr_ifr*100000]
write_csv(agg_ifrs_sum, file='econ/outcome_calculations/data/aggregated_ifrs.csv')

standard_ifrs <- national_ifrs[,!'attach_code'][, lapply(.SD, median), by=c('country_code','cluster_name','age_grp')]

library(popEpi)

standard_ifrs[age_grp == 1, pop_prop := sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:4])))/
               sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2015, 4:24])))]

# https://seer.cancer.gov/stdpopulations/world.who.html
standard_ifrs[age_grp == 1, pop_prop := (8.860)/100.035]
standard_ifrs[age_grp == 2, pop_prop := (8.690 + 8.600 + 8.470)/100.035]
standard_ifrs[age_grp == 3, pop_prop := (8.220 + 7.930 + 7.610 + 7.150 + 6.590 + 6.040 + 
                                           5.370 + 4.550 + 3.720)/100.035]
standard_ifrs[age_grp == 4, pop_prop := (2.960 + 2.210 + 1.520 + 0.910 + 0.440 + 0.150 + 
                                           0.040 + 0.005)/100.035]

standard_ifrs[, weighted_ifr := ifr*pop_prop]
standard_ifrs_sum <- standard_ifrs[,c('country_code','cluster_name','weighted_ifr')][, lapply(.SD, sum), 
                                                                                  by=c('country_code','cluster_name')]
setnames(standard_ifrs_sum, 'weighted_ifr','aggr_ifr')

ggplot(standard_ifrs_sum) + 
  geom_histogram(aes(x=aggr_ifr*100000, fill=cluster_name), bins=100) +
  # facet_wrap(cluster_name~.) +
  theme_bw() + xlab('Age-aggregated IFR, per 100,000 infections') +
  labs(fill='ITZ') + theme(text = element_text(size=14)) +
  scale_fill_manual(values = cluster_colors2)


ggplot(agg_ifrs_sum) + 
  geom_boxplot(aes(x=cluster_name, y=100000*aggr_ifr, group=cluster_name, fill=cluster_name)) + 
  theme_bw() + scale_y_log10(breaks=c(0,5,10,20,30,50,100)) + labs(fill='ITZ') +
  xlab('ITZ') + ylab('IFR, per 100,000 infections') 

ggplot(national_ifrs) + 
  geom_boxplot(aes(x=age_grp, y=100000*ifr, group=age_grp, fill=age_grp)) + 
  theme_bw() + scale_y_log10()


quantile(national_ifrs[country_code=='BRA'&age_grp==4]$ifr*100000)





