
## PLOTTING PROJECTION OUTPUTS
#setwd("~/Desktop/research asst/Global Code")

source("BS/BS_vaccine_programs.R")

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

vacc_doses <- data.frame()
for(c_code in c('GBR')){#c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses/vacc_doses_", c_code, ".csv")) %>% 
                        mutate(cluster_code = c_code))
}

vacc_doses <- vacc_doses %>% filter(country_code == 'GBR')

clusters <- read_csv("data/new_clustering.csv")

supp.labs <- names(vacc_type_list)
names(supp.labs) <- c(1:length(vacc_type_list))
supp.labs2 <- c("Asia-\nEurope", "Africa", "Europe", "Southern\nAmerica",           
                 "Oceania-\nMelanesia-\nPolynesia", "Eastern and\nSouthern Asia",  
                 "Northern\nAmerica")
names(supp.labs2) <-  c('TUR','GHA','GBR','ARG','AUS','CHN','CAN')

pop2025 <- data.frame(cluster_name = unique(clusters$cluster_name),
                      cluster_code = c('TUR','GHA','GBR','ARG','AUS','CHN','CAN'),
                      tot_pop = NA)
for(i in 1:7){
  pop2025$tot_pop[i] <- 1000*sum(pop_proj_WPP_data %>% filter(Year == 2025) %>% 
                                  filter(name %in% unname(unlist(c(clusters[clusters$cluster_name == pop2025$cluster_name[i],
                                                                            c('country', 'country_altern', 'country_altern_2')])))) %>% 
                                  select(!c(name, Type, Year)))
}


## GRAPHS

vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  group_by(cluster_code, year, name, vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(wasted = grepl('w', name), age_grp = substr(name, 2, 2))

vacc_doses_g %>% 
  ggplot(aes(fill=wasted, alpha=age_grp, y=value/1000000, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = 'Set1', labels = c('Unvaccinated','Vaccinated \n(Wasted)')) +
  xlab('Year') + ylab('Vaccine doses, millions') +
  ggtitle("Annual age-specific vaccine doses given in UK, millions") +
  scale_alpha_discrete(range = c(0.4, 1),
                       labels = c('0-5','5-20','20-65','65+')) +
  facet_grid(.~vacc_program, scales='free_y',
             labeller = labeller(vacc_program = supp.labs,
                                 cluster_code = supp.labs2)) + 
  theme_bw() + labs(fill = "Vaccine recipient",
                    alpha = 'Age group') +
  xlab("Year") + 
  theme(text=element_text(size=14))
# ggsave(paste0("output/plots/BS_plots/doses_given_by_ITZ.png"),
#        width=36,height=24,units="cm")

# vacc_doses_plots <- left_join(vacc_doses_g, pop2025) %>% mutate(dose_pp = value/tot_pop)
# 
# vacc_doses_plots %>% 
#   ggplot(aes(fill=name, y=dose_pp, x=year)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_viridis(discrete = T) +
#   xlab('Year') + ylab('Vaccine doses per individual') +
#   ggtitle("Annual age-specific vaccine doses given per individual (based on 2025 population)") +
#   facet_grid(cluster_code~vacc_program, scales='free_y') + 
#   theme_bw() +
#   xlab("")


## vaccinations over the period

wd_ex <- data.frame()
for(i in 1:length(vaccine_programs)){
  wd_ex <- rbind(wd_ex, (fcn_weekly_demog(country = 'United Kingdom',
                                         hemisphere = 'NH',
                                         pop_coverage = coverage_vec,
                                         imm_duration = vacc_type_list[[vaccine_programs[[i]]$vacc_type]]$imm_duration, # in years 
                                         coverage_pattern = vacc_type_list[[vaccine_programs[[i]]$vacc_type]]$coverage_pattern,
                                         weeks_vaccinating = vaccine_programs[[i]]$weeks_vaccinating) %>% mutate(vp=i)))
}

wd_ex %>% filter(V==T, year < 2035) %>% 
  ggplot() +
  geom_line(aes(x=week, y=value/total_as, group=vp, col=as.factor(vp)), lwd=0.8) +
  theme_bw() + ylab('Vaccinated proportion of population') + xlab('Year') +
  ggtitle(paste0('Population projections')) + 
  facet_grid(age_grp~., scales='free') + 
  theme(text=element_text(size=14)) + 
  labs(color = "Vaccine type") +
  scale_color_manual(labels = c("Current", "Improved \nminimal",
                                "Improved \nefficacy", "Improved \nbreadth",
                                "Universal"), values = 2:6)

# wd_ex %>% filter(V==T, year %in% 2026:2035, age_grp=='0-5') %>% 
#   ggplot() +
#   geom_line(aes(x=week, y=value, group=vp, col=as.factor(vp)), lwd=0.8) +
#   # geom_line(aes(x=week, y=total_as, group=vp, col=as.factor(vp)), lwd=0.8) +
#   theme_bw() + ylab('Vaccinated proportion of population') + xlab('Year') +
#   ggtitle(paste0('Population projections')) + 
#   facet_grid(age_grp~., scales='free') + 
#   theme(text=element_text(size=14)) + 
#   labs(color = "Vaccine type") +
#   scale_color_manual(labels = c("Current", "Improved \nminimal",
#                                 "Improved \nefficacy", "Improved \nbreadth",
#                                 "Universal"), values = 2:6)












