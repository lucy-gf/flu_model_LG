
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

k <- 1 # age-specific targeting strategy, in 1:5
targeting <- paste0('ct_', k)
scenario_name <- c('none', 'base', 'low_cov', 'rel_inf')[1]
n_sims <- 100
supp.labs <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs) <- c(1:5)
supp.labs1 <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs1) <- c(0:4)
supp.labs.strain <- c('Total','Influenza A','Influenza B')
names(supp.labs.strain) <- c('tot','totA','totB')
supp.labs.ITZ <- c("Africa", "Asia-Europe", "Eastern and\nSouthern Asia",
                   "Europe", "Northern\nAmerica", "Oceania-\nMelanesia-\nPolynesia",
                   "Southern\nAmerica")
names(supp.labs.ITZ) <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")
# add labels for age-specific targetings too

cases_plot <- data.frame()
for(c_number in c(1:3, 5:7)){
  c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
              "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
              "Southern America")[c_number]
  c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
  hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]
  clusters <- read_csv("data/new_clustering.csv")
  ITZ <- clusters %>% filter(cluster_name == c_name)
  
  itz_cases <- readRDS(paste0("data/vacc_output/vacc_", c_code, '_', scenario_name, '_', 
                                  targeting, ".rds"))[[1]]
  
  itz_cases$tot <- rowSums(itz_cases[,5:20])
  itz_cases$totA <- rowSums(itz_cases[,5:12])
  itz_cases$totB <- rowSums(itz_cases[,13:20])
  
  itz_cases_agg <- itz_cases %>% select(country, country_code, simulation_index, 
                                        week, tot, totA, totB) %>% 
    pivot_longer(!c(country, country_code, simulation_index, week)) %>% 
    rename(strain = name, total = value) %>% 
    group_by(simulation_index, country, strain) %>% 
    mutate(cum_tot = cumsum(total)) %>% ungroup() %>% 
    group_by(simulation_index, week, strain) %>% 
    mutate(cum_tot_agg = sum(cum_tot)) %>% ungroup() %>% 
    filter(country_code == c_code) %>% 
    select(country_code, simulation_index, strain,
           week, cum_tot_agg) %>% 
    rename(cum_tot = cum_tot_agg) %>%  
    group_by(week, strain) %>% mutate(med_cum = median(cum_tot),
                              L50 = unname(unlist(ci(cum_tot, ci=0.50)[2])),
                              U50 = unname(unlist(ci(cum_tot, ci=0.50)[3])),
                              L95 = unname(unlist(ci(cum_tot, ci=0.95)[2])),
                              U95 = unname(unlist(ci(cum_tot, ci=0.95)[3]))) %>% 
    ungroup()
  
  ## dataframe of cases and vaccination status- and age-specific populations
  cases_plot <- rbind(cases_plot, itz_cases_agg)
}

cases_plot_global <- cases_plot %>% group_by(simulation_index, strain, week) %>% 
  mutate(cum_tot_agg = sum(cum_tot)) %>% ungroup() %>% 
  filter(country_code == 'ARG') %>% 
  select(simulation_index, week, strain, cum_tot_agg) %>% rename(cum_tot = cum_tot_agg) %>% 
  group_by(week, strain) %>% mutate(med_cum = median(cum_tot),
                            L50 = unname(unlist(ci(cum_tot, ci=0.50)[2])),
                            U50 = unname(unlist(ci(cum_tot, ci=0.50)[3])),
                            L95 = unname(unlist(ci(cum_tot, ci=0.95)[2])),
                            U95 = unname(unlist(ci(cum_tot, ci=0.95)[3]))) %>% 
  ungroup() %>% filter(simulation_index == 1)

cases_plot %>% filter(simulation_index == 1) %>% 
  ggplot() +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, fill = strain), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, fill = strain), alpha=0.4) +
  geom_line(aes(x=week, y=med_cum/1000000), lwd=0.6) +
  scale_fill_brewer(labels = c('Total', 'INF_A', 'INF_B'), palette='Set2') +
  theme_bw() + facet_grid(country_code~strain, scales='free',
                          labeller = labeller(strain = supp.labs.strain,
                                              country_code = supp.labs.ITZ)) +
  xlab('Year') + ylab('Cumulative cases, millions') +
  labs(fill = 'Strain') + 
  theme(text=element_text(size=14))

cases_plot %>% 
  ggplot() +
  geom_line(aes(x=week, y=cum_tot/1000000, group=simulation_index, col=strain), lwd=0.3, alpha=0.2) +
  theme_bw() + facet_grid(country_code~strain, scales='free',
                          labeller = labeller(strain = supp.labs.strain,
                                              country_code = supp.labs.ITZ)) +
  scale_color_brewer(labels = c('Total', 'INF_A', 'INF_B'), palette='Set2') +
  labs(color = 'Strain') + 
  xlab('Year') + ylab('Cumulative cases, millions') +
  theme(text=element_text(size=14))

cases_plot_global %>% 
  ggplot() +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, fill = strain), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, fill = strain), alpha=0.4) +
  geom_line(aes(x=week, y=med_cum/1000000), lwd=0.6) +
  scale_fill_brewer(labels = c('Total', 'INF_A', 'INF_B'), palette='Set2') +
  theme_bw() + facet_grid(strain~., scales='free',
                          labeller = labeller(strain = supp.labs.strain,
                                              country_code = supp.labs.ITZ)) +
  xlab('Year') + ylab('Cumulative cases, millions') +
  labs(fill = 'Strain') + 
  theme(text=element_text(size=14))



















