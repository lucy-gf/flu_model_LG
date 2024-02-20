
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

c_number <- 4
c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
            "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
            "Southern America")[c_number]
c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]
clusters <- read_csv("data/new_clustering.csv")
ITZ <- clusters %>% filter(cluster_name == c_name)

supp.labs <- names(vacc_type_list)
names(supp.labs) <- c(1:5)
supp.labs1 <- names(vacc_type_list)
names(supp.labs1) <- c(0:4)
supp.labs2 <- c('Overall', 'VE varying', 'Waning varying')
names(supp.labs2) <- c(1:3)

## dataframe of cases and vaccination status- and age-specific populations
n_sims <- 100
cases_list <- readRDS(paste0("data/vacc_output/vacc_", c_code, "_waning_changed_high_cov.rds")) 

cases <- data.frame()
for(i in 1:length(vaccine_programs)){
  cases_list[[i]] <- cases_list[[i]] %>% filter(country_code == 'GBR') %>% mutate(vp = i,
                                                                                  vacc_type = ((i-1) %% 5),
                                                                                  scen = (vp - vacc_type + 4)/5)
  cases <- rbind(cases, cases_list[[i]])
}

cases_plot <- cases #%>%
#   group_by(vp, simulation_index, country) %>% mutate(inout = sum(IU1A)) %>% 
#   ungroup() %>% filter(inout > 0)

## GRAPHS

cases_plot$tot <- rowSums(cases_plot[,5:20])
cases_plot <- cases_plot %>% 
  group_by(simulation_index, vp, country) %>% mutate(cum_tot = cumsum(tot)) %>% 
  ungroup() %>% 
  group_by(week, vp, country) %>% mutate(med_cum = median(cum_tot),
                                L50 = unname(unlist(ci(cum_tot, ci=0.50)[2])),
                                U50 = unname(unlist(ci(cum_tot, ci=0.50)[3])),
                                L95 = unname(unlist(ci(cum_tot, ci=0.95)[2])),
                                U95 = unname(unlist(ci(cum_tot, ci=0.95)[3]))) %>% 
  ungroup()


cases_plot %>% filter(simulation_index==1) %>% ggplot() +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000), fill = 2, alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000), fill = 2, alpha=0.4) +
  geom_line(aes(x=week, y=med_cum/1000000), lwd=0.6) +
  facet_grid(.~vp, labeller = labeller(vp = supp.labs)) +
  theme_bw() + 
  scale_fill_brewer(labels = c("Current", "Improved \nminimal",
                               "Improved \nefficacy", "Improved \nbreadth",
                               "Universal")) +
  xlab('Year') + ylab('Cumulative influenza infections (millions)') +
  ggtitle('Cumulative influenza infections') +
  labs(fill='Vaccine type') + 
  theme(text=element_text(size=14)) +
  ylim(c(0,1000))

cases_plot %>% filter(simulation_index==1) %>% ggplot() +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, fill = as.factor(vacc_type)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, fill = as.factor(vacc_type)), alpha=0.4) +
  geom_line(aes(x=week, y=med_cum/1000000), lwd=0.6) +
  facet_grid(scen~vacc_type, labeller = labeller(scen = supp.labs2, 
                                                 vacc_type = supp.labs1)) +
  theme_bw() + 
  scale_fill_brewer(labels = c("Current", "Improved \nminimal",
                               "Improved \nefficacy", "Improved \nbreadth",
                               "Universal"), palette = 'Set1') +
  xlab('Year') + ylab('Cumulative influenza infections (millions)') +
  ggtitle('Cumulative influenza infections') +
  labs(fill='Vaccine type') + 
  theme(text=element_text(size=14)) +
  ylim(c(0,800))

cases_plot %>% filter(simulation_index==1) %>% ggplot() +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, fill = as.factor(vp)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, fill = as.factor(vp)), alpha=0.4) +
  geom_line(aes(x=week, y=med_cum/1000000), lwd=0.6) +
  facet_grid(.~vp, labeller = labeller(vp = supp.labs)) +
  theme_bw() + 
  scale_fill_brewer(labels = c("Current", "Improved \nminimal",
                                             "Improved \nefficacy", "Improved \nbreadth",
                                             "Universal"), palette = 'Set1') +
  xlab('Year') + ylab('Cumulative influenza infections (millions)') +
  ggtitle('Cumulative influenza infections') +
  labs(fill='Vaccine type') +
  theme(text=element_text(size=14))












