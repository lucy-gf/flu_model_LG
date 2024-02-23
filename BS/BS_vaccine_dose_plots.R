
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

supp.labs <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs) <- c(1:length(vacc_type_list))
supp.labs2 <- c("Asia-\nEurope", "Africa", "Europe", "Southern\nAmerica",           
                 "Oceania-\nMelanesia-\nPolynesia", "Eastern and\nSouthern Asia",  
                 "Northern\nAmerica")
names(supp.labs2) <-  c('TUR','GHA','GBR','ARG','AUS','CHN','CAN')

## GRAPHS

vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  group_by(cluster_code, year, name, vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(wasted = grepl('w', name), age_grp = substr(name, 2, 2))

vacc_doses_g %>% filter(!age_grp == 3) %>% 
  ggplot(aes(fill=wasted, alpha=age_grp, y=value/1000000, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = 'Set1', labels = c('Unvaccinated','Vaccinated \n(Ineffective)')) +
  xlab('Year') + ylab('Vaccine doses, millions') +
  ggtitle("Annual age-specific vaccine doses given in UK, millions") +
  scale_alpha_discrete(range = c(0.4, 1),
                       labels = c('0-5','5-20','65+')) +
  facet_grid(.~vacc_program, scales='free_y',
             labeller = labeller(vacc_program = supp.labs,
                                 cluster_code = supp.labs2)) + 
  theme_bw() + labs(fill = "Vaccine recipient",
                    alpha = 'Age group') +
  xlab("Year") + 
  theme(text=element_text(size=14))
# ggsave(paste0("output/plots/BS_plots/doses_given_by_ITZ.png"),
#        width=36,height=24,units="cm")

vacc_doses_g %>% filter(!age_grp == 3) %>% select(!c(wasted, age_grp)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(prop1 = w1/(w1+v1), prop2 = w2/(w2+v2), prop4 = w4/(w4+v4), 
         tot_prop = (w1 + w2 + w4)/(v1 + v2 + v4 + w1 + w2 + w4)) %>%
  pivot_longer(!c(cluster_code, year, vacc_program)) %>% 
  filter(grepl('prop', name)) %>% mutate(age_grp = case_when(
    name == 'prop1' ~ '0-5',  name == 'prop2' ~ '5-20',
    name == 'prop4' ~ '65+',  name == 'tot_prop' ~ 'Total'
  )) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = value, group = age_grp, col = age_grp), lwd=1) +
  scale_color_brewer(palette = 'Set2') +
  xlab('Year') + ylab('Proportion of vaccine doses ineffective') +
  ggtitle("Proportion of annual age-specific vaccine doses ineffective") +
  facet_grid(.~vacc_program, scales='free_y',
             labeller = labeller(vacc_program = supp.labs,
                                 cluster_code = supp.labs2)) + 
  theme_bw() + xlab("Year") + 
  theme(text=element_text(size=14))
# ggsave(paste0("output/plots/BS_plots/doses_given_by_ITZ.png"),
#        width=36,height=24,units="cm")




## vaccinations over the period

wd_ex <- data.frame()
for(i in 1:length(vaccine_programs_base)){
  wd_ex <- rbind(wd_ex, (fcn_weekly_demog(country = 'Canada',
                                         hemisphere = 'NH',
                                         pop_coverage = vaccine_programs_base[[i]]$pop_coverage,
                                         imm_duration = vacc_type_list[[vaccine_programs_base[[i]]$vacc_type]]$imm_duration, # in years 
                                         coverage_pattern = vacc_type_list[[vaccine_programs_base[[i]]$vacc_type]]$coverage_pattern,
                                         weeks_vaccinating = vaccine_programs_base[[i]]$weeks_vaccinating) %>% mutate(vp = names(vaccine_programs_base)[i])))
}

wd_ex %>% filter(V==T, year < 2035) %>% 
  mutate(vacc_type = substr(vp, 4, 4), targeting = substr(vp, 9, 9)) %>% 
  ggplot() +
  geom_line(aes(x=week, y=value/total_as, group=vacc_type, col=as.factor(vacc_type)), lwd=0.8) +
  theme_bw() + ylab('Vaccinated proportion of population') + xlab('Year') +
  ggtitle(paste0('Population projections')) + 
  facet_grid(age_grp~targeting, scales='free') + 
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












