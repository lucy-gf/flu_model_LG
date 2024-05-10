
## CALCULATING ALL NECESSARY IFRS
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_colors.R")
source("IHME_exemplar_IFRs/fcn_ifr.R")

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

## mortality rates
cdc_data_l <- data.table(read_csv('IHME_exemplar_IFRs/data/extracted_data/ALL_Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv', show_col_types=F))

exemplar_ifrs <- cases_m[,c('country_code','simulation_index','age_grp')]

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

ggplot(exemplar_ifrs) + 
  geom_boxplot(aes(x=country_code, y=100000*ifr, fill=country_code)) + 
  facet_grid(age_grp~., scales='free') + 
  scale_y_log10() +
  theme_bw() + ylab('Deaths per 100,000 infections') + 
  xlab('Country')


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
cdc_data_l <- data.table(read_csv('IHME_exemplar_IFRs/extracted_data/ALL_Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv', show_col_types=F))

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

write_csv(itz_ifrs, file='IHME_exemplar_IFRs/data/itz_ifrs.csv')
itz_ifrs <- data.table(read_csv('IHME_exemplar_IFRs/data/itz_ifrs.csv', 
                                show_col_types=F))

ggplot(itz_ifrs) + 
  geom_boxplot(aes(x=country_code, y=100000*ifr, fill=itz)) + 
  facet_grid(age_grp~itz, scales='free') +
  scale_y_log10() +
  theme_bw() + ylab('Deaths per 100,000 infections') + 
  xlab('Country') 

## maps

national_ifrs <- rbind(exemplar_ifrs, itz_ifrs[,1:4])
national_med_ifrs <- national_ifrs[, lapply(.SD, median), by=c('country_code','age_grp')]
national_med_ifrs[, age_grp1 := paste0('age', age_grp)]
national_med_ifrs_w <- dcast(national_med_ifrs[,c('country_code','age_grp1','ifr')], country_code ~ age_grp1, value.var = 'ifr')

require(maps)
world_map <- data.table(map_data("world") %>% rename(country = region) %>% 
                          mutate(country_code = countrycode(country, origin = 'country.name', destination = 'iso3c')))

all_nations <- clusters[,c('codes','cluster_name')]
setnames(all_nations, 'codes','country_code')
national_med_ifrs_w[, attach_code := country_code]
all_nations[cluster_name %in% c('Africa', 'Asia-Europe', 'Eastern and Southern Asia'), attach_code := country_code]
all_nations[cluster_name %in% c('Oceania-Melanesia-Polynesia'), attach_code := 'AUS']
all_nations[cluster_name %in% c('Southern America'), attach_code := 'BRA']
all_nations[cluster_name %in% c('Northern America'), attach_code := 'CAN']
all_nations[cluster_name %in% c('Europe'), attach_code := 'GBR']
all_nations[country_code %in% c('ARG'), attach_code := 'ARG']
all_nations <- all_nations[national_med_ifrs_w[,!'country_code'], on='attach_code']
  
world_map <- world_map[all_nations, on=c('country_code')]

world_map[,c('order','country','subregion'):=NULL]

world_map_l <- melt(world_map, id.vars = c('long','lat','group','country_code','cluster_name','attach_code'))

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
ggsave(paste0("IHME_exemplar_IFRs/maps/age_spec_ifrs.png"),
       width=50,height=28,units="cm")


cdc_data <- data.table(read_csv('IHME_exemplar_IFRs/extracted_data/Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv'))
cdc_data[, country_code := iso3c]
world_map <- world_map[cdc_data[,c('country_code','cdc_mort_u65_med','cdc_mort_6574_med','cdc_mort_o75_med')], 
                       on=c('country_code')]

world_map_l2 <- melt(world_map, id.vars = c('long','lat','group','country_code','cluster_name','attach_code'))

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
  labs(fill='IFR \nper 100,000') +
  theme(text=element_text(size=14)) 
ggsave(paste0("IHME_exemplar_IFRs/maps/age_spec_cdc_mort.png"),
       width=50,height=28,units="cm")






















