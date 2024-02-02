
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

c_number <- 6
c_name <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
            "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
            "Southern America")[c_number]
c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
hemisphere <- c("NH", "NH", "NH", "NH", "NH", "SH", "SH")[c_number]
clusters <- read_csv("data/new_clustering.csv")
ITZ <- clusters %>% filter(cluster_name == c_name)

supp.labs <- names(vacc_type_list)
names(supp.labs) <- c(1:5)

## dataframe of cases and vaccination status- and age-specific populations
n_sims <- 100
cases <- c()
for(vp_index in 1:length(vaccine_programs)){
  cases <- rbind(cases, 
                        read_csv(paste0("data/vacc_output/vacc_", c_code, "_vp_",
                                        vp_index, ".csv")) %>% 
                          arrange(simulation_index) %>% mutate(vp = vp_index))
  
}

cases_plot <- cases %>% 
  group_by(vp, simulation_index, country) %>% mutate(inout = sum(IU1A)) %>% 
  ungroup() %>% filter(inout > 0)

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

# cases_plot %>% filter(country=='Australia') %>% ggplot() + 
#   geom_line(aes(x=week, y=cum_tot, group=simulation_index, col=as.factor(vp)), 
#             lwd=0.6, alpha=0.3) +
#   geom_line(aes(x=week, y=med_cum, group=simulation_index), lwd=1) +
#   facet_grid(.~vp, labeller = labeller(vp = supp.labs)) +
#   theme_bw() + theme(legend.position = 'none') + 
#   xlab('Year') + ylab('Cumulative influenza infections')

cases_plot %>% ggplot() + 
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, fill=as.factor(vp)),
              alpha=0.4) +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, fill=as.factor(vp)),
              alpha=0.2) +
  geom_line(aes(x=week, y=med_cum/1000000, group=simulation_index), lwd=0.6) +
  facet_grid(.~vp, labeller = labeller(vp = supp.labs)) +
  theme_bw() + theme(legend.position = 'none') + 
  xlab('Year') + ylab('Cumulative influenza infections (millions)') +
  ggtitle("Cumulative influenza infections (millions)") +
  theme(text=element_text(size=14)) + scale_fill_brewer()
# ggsave(paste0("output/plots/BS_plots/cum_cases.png"), 
#        width=32,height=16,units="cm")

hist(unname(unlist((cases_plot %>% filter(vp==3) %>% select(med_cum)) -
  (cases_plot %>% filter(vp==5) %>% select(med_cum)))))

graph_data <- cases_plot %>% select(simulation_index, week, vp, cum_tot) %>% 
  pivot_wider(names_from = vp, values_from = cum_tot) %>% 
  mutate(diff = `3` - `5`) 

graph_data %>% ggplot() +
  geom_line(aes(x=week, y=diff, group=simulation_index, col=simulation_index)) +
  theme_minimal() 



## further plots

cases_plot %>% filter(vp==1, year(week)>2030,
                      simulation_index == 1) %>% ggplot() + 
  geom_line(aes(x=week, y=(IU3A+IU3B+IV3A+IV3B), group=simulation_index, col=simulation_index),
            lwd=0.6) +
  facet_grid(.~vp) +
  theme_minimal() +
  geom_vline(xintercept = as.Date(paste0("01-10-", 2026:2054), '%d-%m-%Y'),
                              lty=3, col=2) +
  geom_vline(xintercept = as.Date(paste0("01-04-", 2026:2054), '%d-%m-%Y'),
                              lty=3, col=3)

cases_plot %>% filter(simulation_index == 1) %>% 
  select(week, vp, U1, V1, U2, V2, U3, V3, U4, V4) %>% 
  pivot_longer(!c(week, vp)) %>%
  mutate(U = grepl('U', name),
         V = grepl('V', name),
         age_grp = case_when(grepl('1', name) ~ '0-5',
                             grepl('2', name) ~ '5-20',
                             grepl('3', name) ~ '20-65',
                             grepl('4', name) ~ '65+')) %>% 
  group_by(week, vp, age_grp) %>% 
  mutate(total_as = sum(value)) %>% ungroup() %>% 
  filter(V==T) %>% 
  ggplot() + 
  geom_line(aes(x=week, y=value/1000000, group=V, col=vp), 
            lwd=0.8) +
  geom_line(aes(x=week, y=total_as/1000000), lwd=0.8) +
  theme_bw() + ylab('Vaccinated population, millions') + xlab('Year') +
  ggtitle(paste0('Vaccinated projections by vaccine type')) + 
  facet_grid(age_grp~vp, scales = 'free',
             labeller = labeller(vp = supp.labs)) + 
  scale_color_viridis(option='D') + 
  theme(legend.position = 'none') +
  theme(text=element_text(size=14))






## MERGING CASES

# age-aggregated
cases_30_100 <- cases_30_100 %>% mutate(tot_A = V1A + V2A + V3A + V4A,
                                        tot_B = V1B + V2B + V3B + V4B,
                                        tot = tot_A + tot_B)
sel_cntr <- countrycode(c_code, origin='iso3c', destination='country.name')

p1 <- cases_30_100 %>% filter(country==sel_cntr,
                              simulation_index == 1) %>% ggplot() +
  geom_line(aes(x=week, y=tot_A)) +
  theme_minimal() + ylab('Influenza A cases') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
  xlab('') + labs(fill = "Strain", color='Strain') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        title=element_text(size=16),
        text=element_text(size=14)) +
  ggtitle(paste0('Example influenza cases, ', sel_cntr))
p2 <- cases_30_100 %>% filter(country==sel_cntr,
                              simulation_index == 1) %>% ggplot() +
  geom_line(aes(x=week, y=tot_B)) +
  theme_minimal() + ylab('Influenza B cases') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
  xlab('') + labs(fill = "Strain", color='Strain') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        title=element_text(size=16),
        text=element_text(size=14)) 
p3 <- cases_30_100 %>% filter(country==sel_cntr,
                              simulation_index == 1) %>% ggplot() +
  geom_line(aes(x=week, y=tot)) +
  theme_minimal() + ylab('Total influenza cases') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
  xlab('') + labs(fill = "Strain", color='Strain') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        title=element_text(size=16),
        text=element_text(size=14)) 
p1 + p2 + p3 +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ',
                  theme = theme(plot.title = element_text(size = 16))) + 
  plot_layout(nrow = 3) 
ggsave(paste0("output/plots/ITZ_no_vacc/ITZ_", c_number, "/exemplar_time_series.png"), 
       width=32,height=24,units="cm")


#cumulative
cases_30_100 <- cases_30_100 %>% 
  group_by(country, simulation_index) %>% arrange(week) %>% 
  mutate(CV1A = cumsum(V1A),CV2A = cumsum(V2A),
         CV3A = cumsum(V3A),CV4A = cumsum(V4A),
         CV1B = cumsum(V1B),CV2B = cumsum(V2B),
         CV3B = cumsum(V3B),CV4B = cumsum(V4B),
         CV1 = CV1A + CV1B, CV2 = CV2A + CV2B,
         CV3 = CV3A + CV3B, CV4 = CV4A + CV4B,
         cum_tot_A = cumsum(tot_A), cum_tot_B = cumsum(tot_B),
         cum_tot = cumsum(tot)) %>% ungroup() %>% 
  arrange(simulation_index, country, week)

# country-specific plots:

## CHANGE TO BE MEDIAN INSTEAD OF MEAN!

for(sel_cntr in ITZ$country){
  cases_30_100_means <- cases_30_100 %>% filter(country==sel_cntr) %>%
    select(!c(country, country_code, V1A:CV4B)) %>%
    pivot_longer(!c(week, simulation_index)) %>% group_by(week, name) %>%
    mutate(L50 = unname(unlist(ci(value, method = "ETI", ci = 0.5)[2])),
           U50 = unname(unlist(ci(value, method = "ETI", ci = 0.5)[3])),
           L95 = unname(unlist(ci(value, method = "ETI", ci = 0.95)[2])),
           U95 = unname(unlist(ci(value, method = "ETI", ci = 0.95)[3])),
           mean = mean(value)) %>%
    ungroup() %>% filter(simulation_index == 1)

  country_strain <- cases_30_100_means %>% filter(name %in% c('cum_tot_A', 'cum_tot_B', 'cum_tot')) %>% ggplot() +
    geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, group=name, fill=name), alpha=0.3) +
    geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, group=name, fill=name), alpha=0.2) +
    geom_line(aes(x=week, y=mean/1000000, group=name, col=name), lwd=1) +
    theme_minimal() + ylab('Cumulative influenza cases (millions)') +
    scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
    xlab('') + labs(fill = "Strain", color='Strain') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          title=element_text(size=16),
          text=element_text(size=14)) +
    ggtitle(paste0('Cumulative influenza cases, ', sel_cntr, ', mean projection, 50% and 95% CI')) +
    scale_color_manual(values = c("black", "darkorange", "mediumorchid"),
                       labels = c('Total', 'Influenza A', 'Influenza B')) +
    scale_fill_manual(values = c("black", "darkorange", "mediumorchid"),
                      labels = c('Total', 'Influenza A', 'Influenza B'))

  country_age <- cases_30_100_means %>% filter(name %in% c('cum_tot', 'CV1', 'CV2', 'CV3', 'CV4')) %>%
    ggplot() +
    geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, group=name, fill=name), alpha=0.3) +
    geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, group=name, fill=name), alpha=0.2) +
    geom_line(aes(x=week, y=mean/1000000, group=name, col=name), lwd=1) +
    theme_minimal() + ylab('Cumulative influenza cases (millions)') +
    scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
    xlab('') + labs(fill = "Age group", color='Age group') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          title=element_text(size=16),
          text=element_text(size=14)) +
    #ggtitle(paste0('Cumulative influenza cases, ', c_name, ', mean projection, 50% and 95% CI')) +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                                  "#F0E442"),
                       labels = c('Total', '0-5', '5-20', '20-65', '65+')) +
    scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                      labels = c('Total', '0-5', '5-20', '20-65', '65+'))

  country_strain + country_age +
    plot_annotation(tag_levels = 'a', tag_prefix = '(',
                    tag_suffix = ')  ',
                    theme = theme(plot.title = element_text(size = 16))) +
    plot_layout(nrow = 2)
  ggsave(paste0("output/plots/ITZ_no_vacc/ITZ_", c_number, "/",
         sel_cntr, ".png"),
         width=32,height=24,units="cm")

  # if(sel_cntr == ITZ$country[1]){
  #   countries_strain <- country_strain
  #   countries_age <- country_age
  # }else{
  #   countries_strain <- countries_strain + country_strain
  #   countries_age <- countries_age + country_age
  # }
  print(paste0(sel_cntr))
}

# countries_strain + plot_layout(guides='collect')
# ggsave(paste0("output/plots/ITZ_no_vacc/ITZ_", c_number, "/all_strain.png"),
#        width=32,height=24,units="cm")
#
# countries_age  + plot_layout(guides='collect')
# ggsave(paste0("output/plots/ITZ_no_vacc/ITZ_", c_number, "/all_age.png"),
#        width=32,height=24,units="cm")

### ITZ totals ###

print("onto ITZs")

ITZ_cases <- cases_30_100 %>% 
  group_by(simulation_index, week) %>% 
  summarise(sum(cum_tot_A), sum(cum_tot_B), sum(cum_tot),
            sum(CV1), sum(CV2), sum(CV3), sum(CV4)) %>% 
  rename(SCTA = `sum(cum_tot_A)`, SCTB = `sum(cum_tot_B)`, SCT = `sum(cum_tot)`, 
         SCV1 = `sum(CV1)`, SCV2 = `sum(CV2)`, SCV3 = `sum(CV3)`, SCV4 = `sum(CV4)`)

# CIs

ITZ_cases1 <- ITZ_cases %>% 
  pivot_longer(!c(week, simulation_index)) %>% group_by(week, name) %>% 
  mutate(L50 = unname(unlist(ci(value, method = "ETI", ci = 0.5)[2])),
         U50 = unname(unlist(ci(value, method = "ETI", ci = 0.5)[3])),
         L95 = unname(unlist(ci(value, method = "ETI", ci = 0.95)[2])),
         U95 = unname(unlist(ci(value, method = "ETI", ci = 0.95)[3])),
         mean = mean(value)) %>% 
  ungroup() %>% filter(simulation_index == 1)

# ITZ-aggregated plots

ITZ_strain <- ITZ_cases1 %>% filter(name %in% c('SCT', 'SCTA', 'SCTB')) %>% ggplot() +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, group=name, fill=name), alpha=0.3) +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, group=name, fill=name), alpha=0.2) +
  geom_line(aes(x=week, y=mean/1000000, group=name, col=name), lwd=1) +
  theme_minimal() + ylab('Cumulative influenza cases (millions)') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
  xlab('') + labs(fill = "Strain", color='Strain') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        title=element_text(size=16),
        text=element_text(size=14)) +
  ggtitle(paste0('Cumulative influenza cases, ', c_name, ', mean projection, 50% and 95% CI')) +
  scale_color_manual(values = c("black", "darkorange", "mediumorchid"),
                     labels = c('Total', 'Influenza A', 'Influenza B')) +
  scale_fill_manual(values = c("black", "darkorange", "mediumorchid"),
                     labels = c('Total', 'Influenza A', 'Influenza B'))

ITZ_age <- ITZ_cases1 %>% filter(name %in% c('SCT', 'SCV1', 'SCV2', 'SCV3', 'SCV4')) %>% 
  ggplot() +
  geom_ribbon(aes(x=week, ymin=L50/1000000, ymax=U50/1000000, group=name, fill=name), alpha=0.3) +
  geom_ribbon(aes(x=week, ymin=L95/1000000, ymax=U95/1000000, group=name, fill=name), alpha=0.2) +
  geom_line(aes(x=week, y=mean/1000000, group=name, col=name), lwd=1) +
  theme_minimal() + ylab('Cumulative influenza cases (millions)') +
  scale_x_date(date_labels="%Y",breaks="1 year",expand=expansion(0.01,0)) +
  xlab('') + labs(fill = "Age group", color='Age group') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        title=element_text(size=16),
        text=element_text(size=14)) +
  #ggtitle(paste0('Cumulative influenza cases, ', c_name, ', mean projection, 50% and 95% CI')) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                "#F0E442"),
                     labels = c('Total', '0-5', '5-20', '20-65', '65+')) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                    labels = c('Total', '0-5', '5-20', '20-65', '65+'))

ITZ_strain + ITZ_age + 
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                                       tag_suffix = ')  ',
                                       theme = theme(plot.title = element_text(size = 16))) + 
  plot_layout(nrow = 2) 
ggsave(paste0("output/plots/ITZ_no_vacc/ITZ_", c_number, "/ITZ_totals.png"), 
       width=32,height=24,units="cm")










