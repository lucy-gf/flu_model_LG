
## making maps of national mortality rates etc.
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_colors.R")

library(readr)
library(dplyr)
library(socialmixr)
library(fluEvidenceSynthesis)
library(data.table)
library(odin)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(bayestestR)
library(viridis)
library(patchwork)
library(countrycode)
library(wpp2022)

## IHME estimates 

ihme_data <- read_xlsx('econ/outcome_calculations/extracted_data/Lancet_Resp_Med_2018_Influenza_LRTI_burden_table1.xlsx')
ihme_data <- data.table(ihme_data)

## adding iso3c codes
ihme_data[, iso3c := countrycode(Region, origin = 'country.name', destination='iso3c')]

# which iso3c conversions didn't work?
# ihme_data[is.na(iso3c),Region] # all income levels, regions, countries within UK, small islands etc.
ihme_data <- ihme_data[!is.na(iso3c),]

## converting mortality rate text into numbers
ihme_data[, ihme_mort_text := gsub(" ","",substr(`Deaths per 100 000 (95% UI)`,1,4))]
ihme_data[, nchar := nchar(ihme_mort_text)]
ihme_data[nchar==3, ihme_mort_rate := as.numeric(paste0(substr(ihme_mort_text,1,1),'.',substr(ihme_mort_text,3,3)))]
ihme_data[nchar==4, ihme_mort_rate := as.numeric(paste0(substr(ihme_mort_text,1,2),'.',substr(ihme_mort_text,4,4)))]

## taking only relevant data
map_data <- ihme_data[,c('iso3c','ihme_mort_rate')]

require(maps)

## making map dataset 
world_map <- data.table(map_data("world") %>% rename(country = region) %>% 
  mutate(iso3c = countrycode(country, origin = 'country.name', destination = 'iso3c')))

## adding ihme mortality rates, and tagging exemplar countries
world_map <- world_map[map_data, on=c('iso3c')]
world_map[, exemplar := 0]
world_map[iso3c %in% c('ARG','AUS','CAN','CHN','GBR','GHA','TUR'),
          exemplar := 1]
world_map[,c('order','country','subregion'):=NULL]

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = ihme_mort_rate), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  scale_fill_viridis(option='A',direction=-1,
                                    breaks = c(0, 4, 8, 12)) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='Influenza \nmortality \nper 100,000') +
  theme(text=element_text(size=14)) +
  ggtitle('GBD mortality estimates (2017)')
ggsave(paste0("econ/outcome_calculations/maps/IHME_mort.png"),
       width=30,height=16,units="cm")

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = as.factor(exemplar)), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  # scale_fill_viridis(option='A',direction=-1,
  #                    breaks = c(0, 4, 8, 12)) +
  xlab('Longitude') + ylab('Latitude') + 
  scale_fill_manual(values = c('0'='white','1'='maroon')) +
  theme(text=element_text(size=14),
        legend.position='none') +
  ggtitle('Exemplar countries')
ggsave(paste0("econ/outcome_calculations/maps/exemplars.png"),
       width=26,height=16,units="cm")

### CDC mortality rates

cdc_data <- read_xlsx('econ/outcome_calculations/extracted_data/Lancet_CDC_2018_Influenza_mortality_Table_S6.xlsx')
cdc_data <- data.table(cdc_data)

## cleaning data
library(stringr)
cdc_data[, country := word(cdc_data$string, 1)] 
cdc_data[, cdc_mort_u65_med := as.numeric(gsub('·', '.', word(cdc_data$string, 2)))] 
cdc_data[, cdc_mort_6574_med := as.numeric(gsub('·', '.', word(cdc_data$string, 6)))] 
cdc_data[, cdc_mort_o75_med := as.numeric(gsub('·', '.', word(cdc_data$string, 10)))] 
cdc_data[, iso3c := countrycode(gsub('_',' ',cdc_data$country), origin='country.name',destination='iso3c')]
write_csv(cdc_data, 'econ/outcome_calculations/extracted_data/Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv')
cdc_data <- data.table(read_csv('econ/outcome_calculations/extracted_data/Clean_Lancet_CDC_2018_Influenza_mortality_Table_S6.csv'))
cdc_data[, continent := countrycode(cdc_data$iso3c, origin='iso3c',destination='continent')]

## adding to map data
world_map <- world_map[cdc_data[,c('iso3c','cdc_mort_u65_med','cdc_mort_6574_med','cdc_mort_o75_med')], 
                       on=c('iso3c')]

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cdc_mort_u65_med), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  scale_fill_viridis(option='A',direction=-1) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='Influenza \nmortality \nper 100,000') +
  theme(text=element_text(size=14))  +
  ggtitle('<65 CDC mortality estimates (2017)')
ggsave(paste0("econ/outcome_calculations/maps/CDC_mort_u65.png"),
       width=30,height=16,units="cm")

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cdc_mort_6574_med), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  scale_fill_viridis(option='A',direction=-1) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='Influenza \nmortality \nper 100,000') +
  theme(text=element_text(size=14)) +
  ggtitle('65-74 CDC mortality estimates (2017)')
ggsave(paste0("econ/outcome_calculations/maps/CDC_mort_6574.png"),
       width=30,height=16,units="cm")

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cdc_mort_o75_med), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  scale_fill_viridis(option='A',direction=-1) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='Influenza \nmortality \nper 100,000') +
  theme(text=element_text(size=14)) +
  ggtitle('75+ CDC mortality estimates (2017)')
ggsave(paste0("econ/outcome_calculations/maps/CDC_mort_o75.png"),
       width=30,height=16,units="cm")

## ITZs 

clusters <- data.table(read_csv('data_for_BS/new_clustering.csv'))
clusters[, iso3c := codes]
world_map <- world_map[clusters[,c('iso3c','cluster_name')], 
                       on=c('iso3c')]

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cluster_name), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  scale_fill_manual(values=cluster_colors2) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='ITZ') +
  theme(text=element_text(size=14)) +
  ggtitle('Influenza Transmission Zones')
ggsave(paste0("econ/outcome_calculations/maps/ITZs.png"),
       width=34,height=16,units="cm")

ggplot(melt(world_map[cluster_name=='Southern America'][,!c('ihme_mort_rate','exemplar')], id.vars = c('long','lat','group','iso3c','cluster_name')), 
       aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (value)), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  facet_grid(.~variable, scales='free') +
  scale_fill_viridis(option='A',direction=-1) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='Influenza \nmortality \nper 100,000') +
  theme(text=element_text(size=14)) +
  ggtitle('CDC mortality estimates')


## no_vacc incidence/attack rates

# nv_cases <- data.table()
# for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
#   nv_cases <- rbind(nv_cases,
#                     cbind(data.table(readRDS(paste0("data/vacc_output_base/vacc_", c_code, '_none_ct_1.rds'))[[1]]), c_code))
#   print(c_code)
# }
# 
# annual_inc <- copy(nv_cases)
# annual_inc <- annual_inc[,c('week','c_code',
#                             'simulation_index'):=NULL]
# annual_inc <- annual_inc[,lapply(.SD, sum), by=c('country','country_code')]
# annual_inc_l <- melt(annual_inc, id.vars = c('country','country_code'))
# annual_inc_l[, variable:=NULL]
# annual_inc_l <- annual_inc_l[, lapply(.SD,sum), by=c('country','country_code')]
# annual_inc_l[, annual_cases := value/(30*100)]
# pop_proj_WPP_data <- data.table(read_csv('data_for_BS/pop_proj_WPP_data.csv'))
# for(name_i in unique(pop_proj_WPP_data$name)){
#   for(j in 1:nrow(clusters)){
#     if(name_i %in% clusters[j, c('country', 'country_altern', 'country_altern_2')]){
#       pop_proj_WPP_data[name == name_i, iso3c := clusters[j, codes]]
#     }
#   }
# }
# 
# for(i in 1:nrow(annual_inc_l)){
#   annual_inc_l[i, pop := 1000*sum(pop_proj_WPP_data[pop_proj_WPP_data$Year == 2040 &
#                                                  pop_proj_WPP_data$iso3c == unlist(annual_inc_l[i,country_code]),
#                                                4:24])]
# }
# annual_inc_l[, attack_r := annual_cases/pop]
# save(annual_inc_l, file="econ/outcome_calculations/data/annual_inc_l")

load("econ/outcome_calculations/data/annual_inc_l")

## adding to map data
annual_inc_l[, iso3c := country_code]
world_map <- world_map[annual_inc_l[,c('iso3c','attack_r')], 
                       on=c('iso3c')]

ggplot(world_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = attack_r), 
               col = 'black', lwd=0.4) +
  theme_bw() +
  scale_fill_viridis(option='A',direction=-1) +
  xlab('Longitude') + ylab('Latitude') + 
  labs(fill='Mean annual \nattack rate') +
  theme(text=element_text(size=14)) +
  ggtitle('Annual attack rate under no vaccinations')
ggsave(paste0("econ/outcome_calculations/maps/no_vacc_attack_rates.png"),
       width=30,height=16,units="cm")








