
## PLOTTING PROJECTION OUTPUTS
#setwd("~/Desktop/research asst/Global Code")

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

clusters <- data.table(read_csv("data/new_clustering.csv"))

eti50L <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[25] + 0.75*sort(x)[26])
  }else{stop('length(x) != 100')}
}
eti50U <- function(x){
  if(length(x) == 100){
    return(0.25*sort(x)[76] + 0.75*sort(x)[75])
  }else{stop('length(x) != 100')}
}
eti95L <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[3] + 0.475*sort(x)[4])
  }else{stop('length(x) != 100')}
}
eti95U <- function(x){
  if(length(x) == 100){
    return(0.525*sort(x)[98] + 0.475*sort(x)[97])
  }else{stop('length(x) != 100')}
}

## VACC DOSES
scenario_name <- c('base', 'low_cov', 'rel_inf')[1]
vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses/vacc_doses_", c_code, "_",
                                                  scenario_name, ".csv")) %>% 
                        mutate(cluster_code = c_code))
}

vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  group_by(cluster_code, year, name, vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(wasted = grepl('w', name), age_grp = substr(name, 2, 2)) %>%
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))

# vacc_doses_g %>% filter(!age_grp == 3, age_cov==5) %>% 
#   ggplot(aes(fill=wasted, y=value/1000000, x=year)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_brewer(palette = 'Set1', labels = c('Unvaccinated','Vaccinated \n(Ineffective)')) +
#   xlab('Year') + ylab('Vaccine doses, millions') +
#   #ggtitle("Annual age-specific vaccine doses given in UK, millions") +
#   # scale_alpha_discrete(range = c(0.4, 1),
#   #                      labels = c('0-5','5-20','65+')) +
#   facet_grid(.~vacc_type, scales='free_y',
#              labeller = labeller(vacc_type = supp.labs)) + 
#   theme_bw() + labs(fill = "Vaccine recipient",
#                     alpha = 'Age group') +
#   xlab("Year") + theme(text=element_text(size=14))
# ggsave(paste0("output/plots/BS_plots/vacc_doses/doses_given_by_VS.png"),
#        width=30,height=7,units="cm")
# 
# vacc_doses_g %>% filter(!age_grp == 3, year==2054) %>% select(!c(wasted, age_grp)) %>%
#   pivot_wider(names_from = name, values_from = value) %>%
#   mutate(prop1 = w1/(w1+v1), prop2 = w2/(w2+v2), prop4 = w4/(w4+v4),
#          tot_prop = (w1 + w2 + w4)/(v1 + v2 + v4 + w1 + w2 + w4)) %>%
#   pivot_longer(!c(cluster_code, year, vacc_program, vacc_type, age_cov)) %>%
#   filter(grepl('prop', name)) %>% mutate(age_grp = case_when(
#     name == 'prop1' ~ '0-5',  name == 'prop2' ~ '5-20',
#     name == 'prop4' ~ '65+',  name == 'tot_prop' ~ 'Total'
#   )) %>% filter(age_cov == 5) %>% 
#   ggplot() +
#   geom_bar(aes(fill=age_grp, x=as.factor(vacc_type), y=value), position="dodge", stat="identity") +
#   scale_fill_manual(values = age_colors) +
#   # facet_grid(age_cov~vacc_type, labeller = labeller(vacc_type = supp.labs,
#   #                                                   age_cov = supp.labs.cov)) +
#   xlab('') + ylab('Proportion of vaccine doses ineffective') +
#   theme_bw() + labs(fill='Age group') + 
#   scale_y_continuous(limits = c(0,0.6), breaks=seq(0,0.6,0.1)) + 
#   theme(text=element_text(size=14)) +
#   scale_x_discrete(labels=c("Current", "Improved\n(minimal)",
#                             "Improved\n(efficacy)", "Improved\n(breadth)",
#                             "Universal")) +
#   theme(axis.text.x = element_text(color=1))
# ggsave(paste0("output/plots/BS_plots/vacc_doses/doses_ineffective_by_age.png"),
#        width=22,height=10,units="cm")

doses_plot <- vacc_doses_g %>% filter(!age_grp == 3) %>% 
  ggplot(aes(fill=age_grp, y=value/1000000, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  xlab('Year') + ylab('Annual vaccine doses, millions') +
  scale_fill_manual(values = age_colors1, labels = c('0-5','5-20', '65+')) +
  facet_grid(age_cov~vacc_type, scales='fixed',
             labeller = labeller(vacc_type = supp.labs,
                                 age_cov = supp.labs.age)) + 
  theme_bw() + labs(fill = "Age group") +
  xlab("Year") + 
  theme(text=element_text(size=14)); doses_plot
# ggsave(paste0("output/plots/BS_plots/vacc_doses/doses_given_by_age.png"),
#        width=26,height=24,units="cm")


### NATIONAL-LEVEL PLOTS

load('data/vacc_output/nat_ann.Rdata') # nrow 186*30*26*5 = 725400
load('data/vacc_output/global_weekly.Rdata') # nrow 1565*100*26 = 4069000

nat_ann_m <- data.table(melt(nat_ann, id.vars = c("country_code", "year", "scenario", "measure")))
nat_ann_m[grepl('U', variable), vacc_status := "U"]
nat_ann_m[grepl('V', variable), vacc_status := "V"]
nat_ann_m[grepl('A', variable), strain := "A"]
nat_ann_m[grepl('B', variable), strain := "B"]
nat_ann_m[grepl('1', variable), age_grp := "0-5"]
nat_ann_m[grepl('2', variable), age_grp := "5-20"]
nat_ann_m[grepl('3', variable), age_grp := "20-65"]
nat_ann_m[grepl('4', variable), age_grp := "65+"]
nat_ann_m[grepl('no_vacc', scenario), ct := "0"]
nat_ann_m[grepl('no_vacc', scenario), vt := "0"]
nat_ann_m[grepl('ct_1', scenario), ct := "1"]
nat_ann_m[grepl('ct_2', scenario), ct := "2"]
nat_ann_m[grepl('ct_3', scenario), ct := "3"]
nat_ann_m[grepl('ct_4', scenario), ct := "4"]
nat_ann_m[grepl('ct_5', scenario), ct := "5"]
nat_ann_m[grepl('vt_1', scenario), vt := "1"]
nat_ann_m[grepl('vt_2', scenario), vt := "2"]
nat_ann_m[grepl('vt_3', scenario), vt := "3"]
nat_ann_m[grepl('vt_4', scenario), vt := "4"]
nat_ann_m[grepl('vt_5', scenario), vt := "5"]

nat_ann_novs <- nat_ann_m[measure == 'median',]
nat_ann_novs[, c('ct','vt','variable','vacc_status','measure','strain','age_grp','year') := NULL]
nat_ann_novs <- nat_ann_novs[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'scenario')]
nat_ann_averted <- copy(nat_ann_novs)
nat_ann_nv <- nat_ann_novs[scenario == 'no_vacc',]
nat_ann_averted[nat_ann_nv, on = c("country_code"), nv := i.value]
nat_ann_averted[, averted := nv - value]
nat_ann_max <- nat_ann_averted[scenario == 'ct_5_vt_5',]
nat_ann_averted[nat_ann_max, on = c("country_code"), max_av := i.averted]
nat_ann_averted[, prop_of_max := averted/max_av]
nat_ann_averted[, ct := substr(scenario, 4,4)]
nat_ann_averted[, vt := substr(scenario, 9,9)]

ggplot(nat_ann_averted[!scenario=='no_vacc',]) +
  geom_bar(aes(x=country_code, y=prop_of_max, group=vt, fill=vt), 
           position="dodge", stat="identity") +
  facet_grid(ct~., scales='free', labeller = labeller(ct = supp.labs.age)) +
  theme_bw() + ylab('Relative averted influenza cases') + xlab('Country') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type')

clusters[,country_code := codes]
nat_ann_averted[clusters, on = c("country_code"), itz := i.cluster_name]

for(c_name in unique(clusters$cluster_name)){
  ggplot(nat_ann_averted[!scenario=='no_vacc' & itz==c_name,]) +
    geom_bar(aes(x=country_code, y=prop_of_max, group=vt, fill=vt), 
             position="dodge", stat="identity") +
    facet_grid(ct~., scales='fixed', labeller = labeller(ct = supp.labs.age)) +
    theme_bw() + ylab('Relative averted influenza cases') + xlab('Country') +
    scale_fill_manual(values=vt_colors, labels = supp.labs) +
    labs(fill = 'Vaccine type') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0('output/plots/BS_plots/cumulative/rel_averted/rel_av_',c_name,'.png'),
         width=38,height=24,units="cm")
}

### GLOBAL PLOTS

global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]

global_weekly_cum <- copy(global_weekly_m)
global_weekly_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario)]
global_weekly_cum[, value:=NULL]

global_weekly_nv <- global_weekly_cum[scenario=='no_vacc',]
global_weekly_cum[global_weekly_nv, on = c("simulation_index","week"), nv_cum := i.cum.sum]
global_weekly_cum[,averted := nv_cum-cum.sum]
global_weekly_cum[,nv_cum:=NULL]

global_weekly_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- global_weekly_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario')]
  dt[, measure := meas]
  global_weekly_meas <- rbind(global_weekly_meas, dt)
}
global_weekly_meas[,simulation_index:=NULL]

global_cumulative <- copy(global_weekly_meas)
global_cumulative[,averted:=NULL]
global_cumulative_wide <- dcast.data.table(global_cumulative, 
                                            week+scenario~measure,
                                            value.var = "cum.sum")
global_cumulative_wide[, ct := substr(scenario, 4,4)]
global_cumulative_wide[, vt := substr(scenario, 9,9)]

global_averted <- copy(global_weekly_meas)
global_averted[,cum.sum:=NULL]
global_averted_wide <- dcast.data.table(global_averted, 
                                           week+scenario~measure,
                                           value.var = "averted")
global_averted_wide[, ct := substr(scenario, 4,4)]
global_averted_wide[, vt := substr(scenario, 9,9)]


ggplot(global_cumulative_wide[!scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
  scale_fill_manual(values = vt_colors) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Cumulative cases, millions') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14))

ggsave(paste0("output/plots/BS_plots/cumulative/cumulative_cases.png"),
       width=26,height=24,units="cm")

averted_plot <- ggplot(global_averted_wide[!scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
  scale_fill_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Averted cases, millions') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14)); averted_plot

ggsave(paste0("output/plots/BS_plots/cumulative/averted_cases.png"),
  width=26,height=24,units="cm")

ggplot(global_cumulative_wide[scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
  scale_fill_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + 
  xlab('Year') + ylab('Cumulative cases, millions') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14))

ggsave(paste0("output/plots/BS_plots/cumulative/no_vacc_global.png"),
       width=13,height=12,units="cm")

averted_single_facet <- ggplot(global_averted_wide[!scenario=='no_vacc',]) +
  geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, group=vt, fill = as.factor(vt)), alpha=0.2) +
  geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, group=vt, fill = as.factor(vt)), alpha=0.4) +
  geom_line(aes(x=week, y=median/1000000, group=vt, col=as.factor(vt)), lwd=0.6) +
  scale_fill_manual(values = vt_colors,
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_color_manual(values = vt_colors,
                     labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + facet_grid(.~ct, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Averted cases, millions') +
  labs(fill = 'Vaccine type', col = 'Vaccine type') + 
  theme(text=element_text(size=14)); averted_single_facet

ggsave(paste0("output/plots/BS_plots/cumulative/averted_cases2.png"),
       width=26,height=24,units="cm")

averted_single_facet + doses_plot + plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')

ggsave(paste0("output/plots/BS_plots/cumulative/averted_and_doses.png"),
       width=52,height=24,units="cm")

print_df <- global_cumulative_wide[week==max(global_cumulative_wide$week)] 
write_csv(print_df, file="data/vacc_output/cases2054.csv")

global_annual <- copy(global_weekly_m)
global_annual[,year:=as.numeric(substr(week,1,4))]
global_annual[,week:=NULL]
global_annual <- global_annual[, lapply(.SD, sum, na.rm=T), by = c('simulation_index','scenario','year')]

global_annual_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- global_annual[, lapply(.SD, get(meas)), by = c('year', 'scenario')]
  dt[, measure := meas]
  global_annual_meas <- rbind(global_annual_meas, dt)
}
global_annual_meas[,simulation_index:=NULL]
global_annual_wide <- dcast.data.table(global_annual_meas, 
                                           year+scenario~measure,
                                           value.var = "value")
global_annual_wide[, ct := substr(scenario, 4,4)]
global_annual_wide[, vt := substr(scenario, 9,9)]

ggplot(global_annual_wide[!scenario=='no_vacc',]) +
  geom_bar(aes(x=year, y=median/1000000, fill = vt), position="stack", stat="identity") +
  scale_color_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values = vt_colors, 
                    labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
  theme_bw() + facet_grid(ct~vt, scales='fixed',
                          labeller = labeller(vt = supp.labs,
                                              ct = supp.labs.age)) +
  xlab('Year') + ylab('Annual cases, millions') +
  labs(fill = 'Vaccine type') + 
  theme(text=element_text(size=14)) +
  geom_errorbar(aes(x = year, ymin = eti95L/1000000, ymax = eti95U/1000000), 
                width = 0.2, alpha=0.6) 

ggsave(paste0("output/plots/BS_plots/cumulative/annual_cases.png"),
       width=30,height=24,units="cm")


## CUMULATIVE/AVERTED PLOTS FOR EACH ITZ

for(c_number in 1:7){
  c_code <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")[c_number]
  if(file.exists(paste0('data/vacc_output/global_weekly_', c_code, '.Rdata'))){
    load(paste0('data/vacc_output/global_weekly_', c_code, '.Rdata'))
  }  
  
  global_weekly_m <- data.table(melt(global_weekly_c, id.vars = c("simulation_index", "week", "scenario")))
  global_weekly_m[, variable:=NULL]
  global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]
  
  global_weekly_cum <- copy(global_weekly_m)
  global_weekly_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario)]
  global_weekly_cum[, value:=NULL]
    
  global_weekly_nv <- global_weekly_cum[scenario=='no_vacc',]
  global_weekly_cum[global_weekly_nv, on = c("simulation_index","week"), nv_cum := i.cum.sum]
  global_weekly_cum[,averted := nv_cum-cum.sum]
  global_weekly_cum[,nv_cum:=NULL]
  
  global_weekly_meas <- data.table()
  for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
    dt <- global_weekly_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario')]
    dt[, measure := meas]
    global_weekly_meas <- rbind(global_weekly_meas, dt)
  }
  global_weekly_meas[,simulation_index:=NULL]
  
  global_cumulative <- copy(global_weekly_meas)
  global_cumulative[,averted:=NULL]
  global_cumulative_wide <- dcast.data.table(global_cumulative, 
                                             week+scenario~measure,
                                             value.var = "cum.sum")
  global_cumulative_wide[, ct := substr(scenario, 4,4)]
  global_cumulative_wide[, vt := substr(scenario, 9,9)]
  
  global_averted <- copy(global_weekly_meas)
  global_averted[,cum.sum:=NULL]
  global_averted_wide <- dcast.data.table(global_averted, 
                                          week+scenario~measure,
                                          value.var = "averted")
  global_averted_wide[, ct := substr(scenario, 4,4)]
  global_averted_wide[, vt := substr(scenario, 9,9)]
  
  ggplot(global_cumulative_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
    scale_fill_manual(values = vt_colors) +
    theme_bw() + facet_grid(ct~vt, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Cumulative cases, millions') +
    # ggtitle(paste0(countrycode(c_code, origin='iso3c',destination='country.name'))) +
    labs(fill = 'Vaccine type') + 
    theme(text=element_text(size=14))
  
  ggsave(paste0("output/plots/BS_plots/cumulative/ITZ_weekly/cumulative_cases/",c_code,".png"),
         width=26,height=24,units="cm")
  
  averted_plot <- ggplot(global_averted_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
    scale_fill_manual(values = vt_colors, 
                      labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
    theme_bw() + facet_grid(ct~vt, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Averted cases, millions') +
    labs(fill = 'Vaccine type') + 
    theme(text=element_text(size=14)); averted_plot
  
  ggsave(paste0("output/plots/BS_plots/cumulative/ITZ_weekly/averted_cases/",c_code,".png"),
         width=26,height=24,units="cm")
  
  averted_single_facet <- ggplot(global_averted_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, group=vt, fill = as.factor(vt)), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, group=vt, fill = as.factor(vt)), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000, group=vt, col=as.factor(vt)), lwd=0.6) +
    scale_fill_manual(values = vt_colors,
                      labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
    scale_color_manual(values = vt_colors,
                       labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
    theme_bw() + facet_grid(.~ct, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Averted cases, millions') +
    labs(fill = 'Vaccine type', col = 'Vaccine type') + 
    theme(text=element_text(size=14)); averted_single_facet
  
  ggsave(paste0("output/plots/BS_plots/cumulative/ITZ_weekly/averted_cases_single/",c_code,".png"),
         width=26,height=24,units="cm")
}

## CUMULATIVE/AVERTED PLOTS FOR EACH STRAIN

load('data/vacc_output/global_weekly.Rdata')  

for(strain_index in c('A','B')){
  global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
  global_weekly_m <- global_weekly_m[grepl(strain_index, variable),]
  global_weekly_m[, variable:=NULL]
  global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario")]
  
  global_weekly_cum <- copy(global_weekly_m)
  global_weekly_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario)]
  global_weekly_cum[, value:=NULL]
  
  global_weekly_nv <- global_weekly_cum[scenario=='no_vacc',]
  global_weekly_cum[global_weekly_nv, on = c("simulation_index","week"), nv_cum := i.cum.sum]
  global_weekly_cum[,averted := nv_cum-cum.sum]
  global_weekly_cum[,nv_cum:=NULL]
  
  global_weekly_meas <- data.table()
  for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
    dt <- global_weekly_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario')]
    dt[, measure := meas]
    global_weekly_meas <- rbind(global_weekly_meas, dt)
  }
  global_weekly_meas[,simulation_index:=NULL]
  
  global_cumulative <- copy(global_weekly_meas)
  global_cumulative[,averted:=NULL]
  global_cumulative_wide <- dcast.data.table(global_cumulative, 
                                             week+scenario~measure,
                                             value.var = "cum.sum")
  global_cumulative_wide[, ct := substr(scenario, 4,4)]
  global_cumulative_wide[, vt := substr(scenario, 9,9)]
  
  global_averted <- copy(global_weekly_meas)
  global_averted[,cum.sum:=NULL]
  global_averted_wide <- dcast.data.table(global_averted, 
                                          week+scenario~measure,
                                          value.var = "averted")
  global_averted_wide[, ct := substr(scenario, 4,4)]
  global_averted_wide[, vt := substr(scenario, 9,9)]
  
  ggplot(global_cumulative_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
    scale_fill_manual(values = vt_colors) +
    theme_bw() + facet_grid(ct~vt, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Cumulative cases, millions') +
    # ggtitle(paste0(countrycode(c_code, origin='iso3c',destination='country.name'))) +
    labs(fill = 'Vaccine type') + 
    theme(text=element_text(size=14))
  
  ggsave(paste0("output/plots/BS_plots/cumulative/by_strain/cumulative_cases/",strain_index,".png"),
         width=26,height=24,units="cm")
  
  averted_plot <- ggplot(global_averted_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill = as.factor(vt)), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill = as.factor(vt)), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000), lwd=0.6) +
    scale_fill_manual(values = vt_colors, 
                      labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
    theme_bw() + facet_grid(ct~vt, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Averted cases, millions') +
    labs(fill = 'Vaccine type') + 
    theme(text=element_text(size=14)); averted_plot
  
  ggsave(paste0("output/plots/BS_plots/cumulative/by_strain/averted_cases/",strain_index,".png"),
         width=26,height=24,units="cm")
  
  averted_single_facet <- ggplot(global_averted_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, group=vt, fill = as.factor(vt)), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, group=vt, fill = as.factor(vt)), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000, group=vt, col=as.factor(vt)), lwd=0.6) +
    scale_fill_manual(values = vt_colors,
                      labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
    scale_color_manual(values = vt_colors,
                       labels = c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')) +
    theme_bw() + facet_grid(.~ct, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Averted cases, millions') +
    labs(fill = 'Vaccine type', col = 'Vaccine type') + 
    theme(text=element_text(size=14)); averted_single_facet
  
  ggsave(paste0("output/plots/BS_plots/cumulative/by_strain/averted_cases_single/",strain_index,".png"),
         width=26,height=24,units="cm")
}

## AGE-SPECIFIC CUMULATIVE/AVERTED PLOTS 

load('data/vacc_output/global_weekly.Rdata')  

  global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
  global_weekly_m[grepl('1', variable), age_grp := "0-5"]
  global_weekly_m[grepl('2', variable), age_grp := "5-20"]
  global_weekly_m[grepl('3', variable), age_grp := "20-65"]
  global_weekly_m[grepl('4', variable), age_grp := "65+"]
  global_weekly_m[, variable:=NULL]
  global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "week", "scenario", "age_grp")]
  
  global_weekly_cum <- copy(global_weekly_m)
  global_weekly_cum[, cum.sum := cumsum(value), by=list(simulation_index, scenario, age_grp)]
  global_weekly_cum[, value:=NULL]
  
  global_weekly_nv <- global_weekly_cum[scenario=='no_vacc',]
  global_weekly_cum[global_weekly_nv, on = c("simulation_index","week",'age_grp'), nv_cum := i.cum.sum]
  global_weekly_cum[,averted := nv_cum-cum.sum]
  global_weekly_cum[,nv_cum:=NULL]
  
  global_weekly_meas <- data.table()
  for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
    dt <- global_weekly_cum[, lapply(.SD, get(meas)), by = c('week', 'scenario', 'age_grp')]
    dt[, measure := meas]
    global_weekly_meas <- rbind(global_weekly_meas, dt)
  }
  global_weekly_meas[,simulation_index:=NULL]
  
  global_cumulative <- copy(global_weekly_meas)
  global_cumulative[,averted:=NULL]
  global_cumulative_wide <- dcast.data.table(global_cumulative, 
                                             week+scenario+age_grp~measure,
                                             value.var = "cum.sum")
  global_cumulative_wide[, ct := substr(scenario, 4,4)]
  global_cumulative_wide[, vt := substr(scenario, 9,9)]
  
  global_averted <- copy(global_weekly_meas)
  global_averted[,cum.sum:=NULL]
  global_averted_wide <- dcast.data.table(global_averted, 
                                          week+scenario+age_grp~measure,
                                          value.var = "averted")
  global_averted_wide[, ct := substr(scenario, 4,4)]
  global_averted_wide[, vt := substr(scenario, 9,9)]
  
  ggplot(global_cumulative_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, 
                    group = age_grp, fill = age_grp), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                    group = age_grp, fill = age_grp), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000, 
                  group = age_grp, col = age_grp), lwd=0.6) +
    scale_fill_manual(values = age_colors) +
    scale_color_manual(values = age_colors) +
    theme_bw() + facet_grid(ct~vt, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Cumulative cases, millions') +
    # ggtitle(paste0(countrycode(c_code, origin='iso3c',destination='country.name'))) +
    labs(fill = 'Age group', color = 'Age group') + 
    theme(text=element_text(size=14))
  
  ggsave(paste0("output/plots/BS_plots/cumulative/by_age/cumulative_cases.png"),
         width=26,height=24,units="cm")
  
  ggplot(global_averted_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, 
                    group = age_grp, fill = age_grp), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                    group = age_grp, fill = age_grp), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000, 
                  group = age_grp, col = age_grp), lwd=0.6) +
    scale_fill_manual(values = age_colors) +
    scale_color_manual(values = age_colors) +
    theme_bw() + facet_grid(ct~vt, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Cumulative cases, millions') +
    # ggtitle(paste0(countrycode(c_code, origin='iso3c',destination='country.name'))) +
    labs(fill = 'Age group', color = 'Age group') + 
    theme(text=element_text(size=14))
  
  ggsave(paste0("output/plots/BS_plots/cumulative/by_age/averted_cases.png"),
         width=26,height=24,units="cm")
  
  ggplot(global_averted_wide[!scenario=='no_vacc',]) +
    geom_ribbon(aes(x=week, ymin=eti95L/1000000, ymax=eti95U/1000000, 
                    group = ct, fill = ct), alpha=0.2) +
    geom_ribbon(aes(x=week, ymin=eti50L/1000000, ymax=eti50U/1000000, 
                    group = ct, fill = ct), alpha=0.4) +
    geom_line(aes(x=week, y=median/1000000, 
                  group = ct, col = ct), lwd=0.6) +
    scale_fill_manual(values = vt_colors) +
    scale_color_manual(values = vt_colors) +
    theme_bw() + facet_grid(vt~age_grp, scales='fixed',
                            labeller = labeller(vt = supp.labs,
                                                ct = supp.labs.age)) +
    xlab('Year') + ylab('Cumulative cases, millions') +
    # ggtitle(paste0(countrycode(c_code, origin='iso3c',destination='country.name'))) +
    labs(fill = 'Coverage plan', color = 'Coverage plan') + 
    theme(text=element_text(size=14))
  
  ggsave(paste0("output/plots/BS_plots/cumulative/by_age/agexvacc.png"),
         width=26,height=24,units="cm")


## NNV
  
nnv <- nat_ann[measure=='median',]
nnv[, c('measure','year'):=NULL]
nnv <- nnv[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'scenario')]
nnv_m <- data.table(melt(nnv, id.vars = c("country_code", "scenario")))
nnv_m[, variable:=NULL]
nnv_m <- nnv_m[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'scenario')]
nnv_m[, ct := substr(scenario, 4,4)]
nnv_m[, vt := substr(scenario, 9,9)]

vd <- data.table(vacc_doses)
vd[, c('country','cluster_code'):=NULL]
vd[,vt := (vacc_program %% 5)]
vd[vt == 0, vt := 5]
vd[,ct := ceiling(vacc_program/5)]
vd[,c('year','vacc_program'):=NULL] 
vd <- vd[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'vt', 'ct')]
vd_m <- data.table(melt(vd, id.vars = c('country_code','vt','ct')))
vd_m[, vt := as.character(vt)]
vd_m[, ct := as.character(ct)]
vd_v <- vd_m[grepl('v', variable)]
vd_m[,variable:=NULL]
vd_v[,variable:=NULL]
vd_all <- vd_m[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'vt', 'ct')]
vd_nonwasted <- vd_v[, lapply(.SD, sum, na.rm=T), by = c('country_code', 'vt', 'ct')]

nnv_m[vd_all, on = c("country_code", 'ct','vt'), all_vd := i.value]
nnv_m[vd_nonwasted, on = c("country_code", 'ct','vt'), nonwasted_vd := i.value]
nnv_m[is.na(all_vd), all_vd:=0]
nnv_m[is.na(nonwasted_vd), nonwasted_vd:=0]

novacc_nnv <- nnv_m[scenario=='no_vacc',]
nnv_m[novacc_nnv, on = c("country_code"), novacc_cases := i.value]
nnv_m[, averted := novacc_cases - value]
nnv_fin <- nnv_m[!scenario=='no_vacc',]
nnv_fin[, all_nnv := all_vd/averted]
nnv_fin[, nw_nnv := nonwasted_vd/averted]

nnv_fin[clusters, on = c("country_code"), itz := i.cluster_name]
nnv_fin[, country := countrycode(country_code, origin='iso3c',destination='country.name')]

for(c_name in unique(nnv_fin$itz)){
  ggplot(nnv_fin[itz==c_name,]) +
    geom_bar(aes(x=country, y=all_nnv, group=vt, fill=vt), 
             position="dodge", stat="identity") +
    facet_grid(ct~., scales='free', labeller = labeller(ct = supp.labs.age)) +
    theme_bw() + ylab('Number needed to vaccinate') + xlab('Country') +
    scale_fill_manual(values=vt_colors, labels = supp.labs) +
    labs(fill = 'Vaccine type') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0('output/plots/BS_plots/NNV/all/',c_name,'.png'),
         width=38,height=24,units="cm")
  
  ggplot(nnv_fin[itz==c_name,]) +
    geom_bar(aes(x=country, y=nw_nnv, group=vt, fill=vt), 
             position="dodge", stat="identity") +
    facet_grid(ct~., scales='free', labeller = labeller(ct = supp.labs.age)) +
    theme_bw() + ylab('Number needed to vaccinate') + xlab('Country') +
    scale_fill_manual(values=vt_colors, labels = supp.labs) +
    labs(fill = 'Vaccine type') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0('output/plots/BS_plots/NNV/nonwasted/',c_name,'.png'),
         width=38,height=24,units="cm")
}


ggplot(nnv_fin) +
  geom_boxplot(aes(x=vt, y=all_nnv, fill=vt, group=vt)) +
  facet_grid(ct~itz, scales='free', labeller = labeller(ct = supp.labs.age,
                                                        itz = supp.labs.ITZ2)) +
  theme_bw() + ylab('Number needed to vaccinate') + xlab('') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type') +
  theme(text=element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(paste0('output/plots/BS_plots/NNV/nnv_boxplot.png'),
       width=30,height=30,units="cm")


## WTF IS WRONG WITH GERMANY DATA

ggplot(nat_ann_m[country_code %in% c('DEU','FRA') &
                   scenario %in% c('no_vacc','ct_5_vt_5'),]) +
  geom_bar(aes(x=year, y=value/1000000, fill=age_grp, group=year),
           position='stack',stat='identity') +
  theme_minimal() + ylab('') + xlab('Vaccine type') +
  scale_fill_manual(values=age_colors) +
  facet_grid(scenario~country_code) +
  labs(fill = 'Age') + ylim(c(0,NA)) +
  theme(text=element_text(size=14))

itz_cases_no_vacc <- data.table(readRDS(paste0("data/vacc_GBR_none_ct_1.rds"))[[1]])

deu <- itz_cases_no_vacc[country %in% c('Germany', 'France'),]
deu[,country_code:=NULL]
deu_m <- data.table(melt(deu, id.vars=c('country','simulation_index','week')))
deu_m[,variable:=NULL]
deu_m2 <- deu_m[, lapply(.SD, sum, na.rm=T), by=c('country', 'simulation_index','week')]

deu_cum <- copy(deu_m2)
deu_cum[, cum.sum := cumsum(value), by=list(country, simulation_index)]

ggplot(deu[country=='Germany'&simulation_index<2&year(week)<2030,]) +
  geom_line(aes(x=week, y=IU3A, col=simulation_index, group=simulation_index)) +
  theme_minimal() + ylab('') + xlab('Week') +
  facet_grid(country~.) +
  ylim(c(0,NA)) +
  theme(text=element_text(size=14))













