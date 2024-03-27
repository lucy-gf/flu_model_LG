
## PLOTTING PROJECTION OUTPUTS
#setwd("~/Desktop/research asst/Global Code")
#source("BS/BS_data_fcns.R")
source("BS/BS_vaccine_programs.R")
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

scenario_name <- c('base', 'low_cov', 'rel_inf')[1]

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses/vacc_doses_", c_code, "_",
                                                  scenario_name, ".csv")) %>% 
                        mutate(cluster_code = c_code))
}

clusters <- read_csv("data/new_clustering.csv")

## GRAPHS

vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  group_by(cluster_code, year, name, vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(wasted = grepl('w', name), age_grp = substr(name, 2, 2)) %>%
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))

vacc_doses_g %>% filter(!age_grp == 3, age_cov==5) %>% 
  ggplot(aes(fill=wasted, y=value/1000000, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = 'Set1', labels = c('Unvaccinated','Vaccinated \n(Ineffective)')) +
  xlab('Year') + ylab('Vaccine doses, millions') +
  #ggtitle("Annual age-specific vaccine doses given in UK, millions") +
  # scale_alpha_discrete(range = c(0.4, 1),
  #                      labels = c('0-5','5-20','65+')) +
  facet_grid(.~vacc_type, scales='free_y',
             labeller = labeller(vacc_type = supp.labs)) + 
  theme_bw() + labs(fill = "Vaccine recipient",
                    alpha = 'Age group') +
  xlab("Year") + theme(text=element_text(size=14))
ggsave(paste0("output/plots/BS_plots/vacc_doses/doses_given_by_VS.png"),
       width=30,height=7,units="cm")

vacc_doses_g %>% filter(!age_grp == 3, year==2054) %>% select(!c(wasted, age_grp)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(prop1 = w1/(w1+v1), prop2 = w2/(w2+v2), prop4 = w4/(w4+v4),
         tot_prop = (w1 + w2 + w4)/(v1 + v2 + v4 + w1 + w2 + w4)) %>%
  pivot_longer(!c(cluster_code, year, vacc_program, vacc_type, age_cov)) %>%
  filter(grepl('prop', name)) %>% mutate(age_grp = case_when(
    name == 'prop1' ~ '0-5',  name == 'prop2' ~ '5-20',
    name == 'prop4' ~ '65+',  name == 'tot_prop' ~ 'Total'
  )) %>% filter(age_cov == 5) %>% 
  ggplot() +
  geom_bar(aes(fill=age_grp, x=as.factor(vacc_type), y=value), position="dodge", stat="identity") +
  scale_fill_manual(values = age_colors) +
  # facet_grid(age_cov~vacc_type, labeller = labeller(vacc_type = supp.labs,
  #                                                   age_cov = supp.labs.cov)) +
  xlab('') + ylab('Proportion of vaccine doses ineffective') +
  theme_bw() + labs(fill='Age group') + 
  scale_y_continuous(limits = c(0,0.6), breaks=seq(0,0.6,0.1)) + 
  theme(text=element_text(size=14)) +
  scale_x_discrete(labels=c("Current", "Improved\n(minimal)",
                            "Improved\n(efficacy)", "Improved\n(breadth)",
                            "Universal")) +
  theme(axis.text.x = element_text(color=1))
ggsave(paste0("output/plots/BS_plots/vacc_doses/doses_ineffective_by_age.png"),
       width=22,height=10,units="cm")


vacc_doses_g %>% filter(!age_grp == 3) %>% 
  ggplot(aes(fill=age_grp, y=value/1000000, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  xlab('Year') + ylab('Annual vaccine doses, millions') +
  scale_fill_manual(values = age_colors1, labels = c('0-5','5-20', '65+')) +
  facet_grid(age_cov~vacc_type, scales='fixed',
             labeller = labeller(vacc_type = supp.labs,
                                 age_cov = supp.labs.cov)) + 
  theme_bw() + labs(fill = "Age group") +
  xlab("Year") + 
  theme(text=element_text(size=14))
ggsave(paste0("output/plots/BS_plots/vacc_doses/doses_given_by_age.png"),
       width=26,height=24,units="cm")



## vaccinated population over the 30-year period

wd_ex <- data.frame()
for(i in 1:length(vaccine_programs_base)){
  wd_ex <- rbind(wd_ex, (fcn_weekly_demog(country = 'Canada',
                                         hemisphere = 'NH',
                                         pop_coverage = vaccine_programs_base[[i]]$pop_coverage,
                                         imm_duration = vacc_type_list[[vaccine_programs_base[[i]]$vacc_type]]$imm_duration, # in years 
                                         coverage_pattern = vacc_type_list[[vaccine_programs_base[[i]]$vacc_type]]$coverage_pattern,
                                         weeks_vaccinating = vaccine_programs_base[[i]]$weeks_vaccinating) %>% mutate(vp = names(vaccine_programs_base)[i])))
}



wd_ex %>% filter(V==T, year<2030) %>% 
  mutate(vacc_type = substr(vp, 4, 4), targeting = substr(vp, 9, 9)) %>% 
  filter(targeting=='1', age_grp=='0-5') %>% 
  ggplot() +
  geom_line(aes(x=week, y=value/total_as, group=vacc_type, col=as.factor(vacc_type)), lwd=0.9) +
  theme_bw() + ylab('Vaccinated proportion') + xlab('Year') +
  #ggtitle(paste0('Vaccinated proportion of 0-5 age group')) + 
  #facet_grid(targeting~age_grp, scales='free') + 
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.1)) +
  geom_hline(yintercept=0.7, lty=2, alpha=0.5) + 
  theme(text=element_text(size=14)) + 
  labs(color = "Vaccine type") +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
  scale_color_manual(labels = c("Current", "Improved\n(minimal)",
                                "Improved\n(efficacy)", "Improved\n(breadth)",
                                "Universal"), values = vt_colors)
ggsave(paste0("output/plots/BS_plots/vacc_doses/pop_immunity_base.png"),
       width=28,height=12,units="cm")

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












