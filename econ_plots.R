#### Economic analysis - plots
setwd("~/Desktop/research asst/Global Code")

source("BS/BS_colors.R")

library(readr)
library(dplyr)
library(data.table)
library(odin)
library(parallel)
library(countrycode)
library(ggplot2)
library(readxl)
library(tidyverse)
library(viridis)
library(patchwork)
library(wpp2022)
library(WDI)
options(scipen=1000000)

scenario_name <- 'base'
econ_folder_name <- '_discount0'
print(paste0(scenario_name,econ_folder_name))

########################################################
### GLOBAL THRESHOLDS ###

print('global')
load(paste0('econ/',scenario_name,econ_folder_name,'/agg_out'))

agg_out$ct_n <- factor(agg_out$ct_n, levels = unique(agg_out$ct_n))

ggplot(agg_out) + 
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt))) +
  ylab('Threshold cost ($)') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  ylim(c(0,NA)) +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/global_thresholds.png"),
       width=25,height=14,units="cm")

bill <- 1000000000

ggplot(agg_out) + 
  # geom_line(aes(x=incremental_DALYs/1000000, y=incremental_costs/1000000, group=simulation_index), alpha=0.5) + 
  geom_point(aes(x=incremental_DALYs/bill, y=incremental_costs/bill, col=as.factor(vt)), shape=1) + 
  facet_grid(ct~., scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                                'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(color = 'Vaccine type') + 
  theme(text=element_text(size=14)) +
  theme_bw() + xlim(c(0,NA)) + ylim(c(0,NA)) +
  xlab('Averted discounted DALYs (billions)') + 
  ylab('Averted discounted costs (billions)')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/global_averted.png"),
       width=25,height=14,units="cm")

### looking at order of expensiveness
agg_agg <- agg_out[,lapply(.SD,median), by=c('ct','vt','ct_n')]
agg_agg$ct_n <- factor(agg_agg$ct_n, levels = unique(agg_agg$ct_n))
ggplot(agg_agg) + 
  geom_bar(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt)),
               position='dodge',stat='identity',col='black') +
  ylab('Median threshold cost ($)') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  ylim(c(0,NA)) +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')

ggplot(agg_agg) + 
  geom_point(aes(x=incremental_costs, y=discounted_DALYs, 
               col=as.factor(vt))) +
  geom_line(aes(x=incremental_costs, y=discounted_DALYs, 
                col=as.factor(vt))) +
  # ylab('Median threshold cost ($)') + 
  scale_color_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(color = 'Vaccine type') + 
  # facet_grid(vt~.) +
  ylim(c(0,NA)) +
  theme(text=element_text(size=14)) +
  theme_bw() 
# so no 'domination'


########################################################
### WHO REGION THRESHOLDS ###
print('WHO regions')

load(paste0('econ/',scenario_name,econ_folder_name,'/total_who_sum'))

total_who_sum$ct_n <- factor(total_who_sum$ct_n, levels = unique(total_who_sum$ct_n))

ggplot(total_who_sum) + 
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt))) +
  ylab('Threshold cost ($)') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  ylim(c(0,NA)) +
  facet_wrap(WHOREGION~., scales='free', nrow=3, labeller = labeller(WHOREGION=who_region_labs)) + 
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/who_region_thresholds.png"),
       width=25,height=20,units="cm")

ggplot(total_who_sum) + 
  # geom_line(aes(x=incremental_DALYs/1000000, y=incremental_costs/1000000, group=simulation_index), alpha=0.5) + 
  geom_point(aes(x=incremental_DALYs/bill, y=incremental_costs/bill, col=as.factor(vt)),shape=1) + 
  facet_grid(ct~., scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                                'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(color = 'Vaccine type') + 
  theme(text=element_text(size=14)) +
  facet_wrap(WHOREGION~., scales='free', nrow=3, labeller = labeller(WHOREGION=who_region_labs)) + 
  theme_bw() + xlim(c(0,NA)) + ylim(c(0,NA)) +
  xlab('Averted DALYs (billions)') + 
  ylab('Averted costs (billions)')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/who_region_averted.png"),
       width=25,height=20,units="cm")

########################################################
### NATIONAL THRESHOLDS ###
print('national')

load(paste0('econ/',scenario_name,econ_folder_name,'/total_nat_sum'))

total_nat_sum$ct_n <- factor(total_nat_sum$ct_n, levels = unique(total_nat_sum$ct_n))
total_nat_sum[vt==1, vt_n:='Current']
total_nat_sum[vt==2, vt_n:='Improved (minimal)']
total_nat_sum[vt==3, vt_n:='Improved (efficacy)']
total_nat_sum[vt==4, vt_n:='Improved (breadth)']
total_nat_sum[vt==5, vt_n:='Universal']
total_nat_sum$vt_n <- factor(total_nat_sum$vt_n, levels=unique(total_nat_sum$vt_n))
  
WHO_regions <- data.table(read_csv('econ/outcome_calculations/data/WHO_regions.csv', show_col_types=F))
total_nat_sum <- total_nat_sum[WHO_regions, on='country_code', WHOREGION := WHOREGION]

## adding gdp!
gdp_data <- data.table(WDI(indicator='NY.GDP.PCAP.KD', start=2022, end=2022))
gdp_data_filt <- gdp_data[iso3c %in% unique(total_nat_sum$country_code)]
setnames(gdp_data_filt, 'NY.GDP.PCAP.KD','gdpcap')
gdp_data_filt[iso3c == 'AFG', gdpcap := 355.78]
gdp_data_filt[iso3c == 'BTN', gdpcap := 3560.20]
gdp_data_filt[iso3c == 'ERI', gdpcap := 643.82]
gdp_data_filt[iso3c == 'PRK', gdpcap := 590]
gdp_data_filt[iso3c == 'LBN', gdpcap := 4136.10]
gdp_data_filt[iso3c == 'NCL', gdpcap := 35745.20]
gdp_data_filt[iso3c == 'SSD', gdpcap := 550.86]
gdp_data_filt[iso3c == 'SYR', gdpcap := 420.62]
gdp_data_filt[iso3c == 'LBN', gdpcap := 4136.10]
gdp_data_filt[iso3c == 'TON', gdpcap := 4426]
gdp_data_filt[iso3c == 'VEN', gdpcap := 15975.73]
gdp_data_filt[, c('iso2c','year') := NULL]
gdp_data_filt <- rbind(gdp_data_filt,
                       data.table(
                         country = c('French Guiana', 'Taiwan'),
                         iso3c = c('GUF','TWN'),
                         gdpcap = c(15600, 32679)
                       ))
setnames(gdp_data_filt, 'iso3c','country_code')
total_nat_sum <- total_nat_sum[gdp_data_filt, on='country_code']
total_nat_sum <- total_nat_sum[order(gdpcap)]

ggplot(total_nat_sum[ct==2]) + 
  geom_boxplot(aes(x=as.factor(gdpcap), y=threshold_price, 
                   fill=WHOREGION), outlier.shape = NA) +
  ylab('Threshold cost ($)') + 
  facet_grid(vt_n~., scales='free', labeller = labeller(WHOREGION=who_region_labs)) +
  guides(fill="none") + scale_y_log10() +
  theme_bw() + xlab('Vaccine type') +
  theme(text=element_text(size=14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank()) 
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/boxplot_all_countries.png"),
       width=40,height=25,units="cm")

# median box-plot
total_nat_sum_limited <- total_nat_sum[,c('country_code','ct','vt','ct_n','vt_n','WHOREGION','threshold_price')]
nat_meds <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- total_nat_sum_limited[, lapply(.SD, get(meas)), by = c('country_code','ct','vt','ct_n','vt_n','WHOREGION')]
  dt[, measure := meas]
  nat_meds <- rbind(nat_meds, dt)
}

total_nat_sum_nnvd <- total_nat_sum[,c('country_code','ct','vt','ct_n','vt_n','discounted_doses','incremental_DALYs')]
total_nat_sum_nnvd[, nnvd := incremental_DALYs/discounted_doses]
total_nat_sum_nnvd[, c('discounted_doses','incremental_DALYs'):=NULL]

ggplot(nat_meds[measure=='median']) + 
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt))) +
  # geom_jitter(aes(x=ct_n, y=threshold_price), alpha = 0.3) +
  ylab('Threshold cost ($)') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  ylim(c(0,NA)) + 
  facet_wrap(WHOREGION~., scales='free', nrow=3, labeller = labeller(WHOREGION=who_region_labs)) + 
  theme_bw() + xlab('Vaccination targets') +
  theme(text=element_text(size=14), legend.position='none') 
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/boxplot_medians_WHOfacet.png"),
       width=35,height=25,units="cm")

### which age-targeting strategy is preferable in diff countries?
pref_strat_dt <- nat_meds[measure=='median']
pref_strat <- data.table(ct=rep(1:5, each=5), vt=rep(1:5,5))
for(i in 1:nrow(pref_strat)){
  pref_strat[i, pref := 0]
  dt <- pref_strat_dt[vt==pref_strat[i,vt]]
  for(c_code in unique(dt$country_code)){
    dt_c <- dt[country_code==c_code][order(-threshold_price)]
    if(pref_strat[i,ct]==as.numeric(dt_c[1,ct])){
      pref_strat[i, pref := pref + 1]
    }
  }
}
pref_strat <- pref_strat[order(vt)]
pref_strat[vt==1, vt_n:='Current']
pref_strat[vt==2, vt_n:='Improved (minimal)']
pref_strat[vt==3, vt_n:='Improved (efficacy)']
pref_strat[vt==4, vt_n:='Improved (breadth)']
pref_strat[vt==5, vt_n:='Universal']
pref_strat$vt_n <- factor(pref_strat$vt_n, levels=unique(pref_strat$vt_n))
ggplot(pref_strat) + 
  geom_col(aes(x=vt_n, y=pref, group=ct, fill=as.factor(ct)),
           position='dodge', col=1) +
  scale_fill_manual(values=age_targ_colors,
                    labels = c('0-4','0-10',
                               '0-17','65+','0-17, 65+')) +
  theme_bw() + xlab('') +
  theme(text=element_text(size=14)) +
  labs(fill='Age-targeting\nstrategy') +
  ylab('Preferred strategy in n countries')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/preferred_strategy.png"),
       width=25,height=15,units="cm")

## saving range of 95% CIs
nat_meds_less <- nat_meds[,c('country_code','ct','vt','WHOREGION','threshold_price','measure')]
save_dt <- nat_meds_less[measure=='median'][, .(min = min(threshold_price)), by=c('ct','vt','WHOREGION')]
save_vec_min <- save_dt$min
save_dt <- nat_meds_less[measure=='median'][, .(max = max(threshold_price)), by=c('ct','vt','WHOREGION')]
save_vec_max <- save_dt$max
nat_meds_w <- dcast(nat_meds_less, country_code + ct + vt + WHOREGION ~ measure, value.var = 'threshold_price')
nat_meds_w <- nat_meds_w[median %in% c(save_vec_min,save_vec_max)]
if(!nrow(nat_meds_w)==(length(save_vec_min) + length(save_vec_max))){print('nat_meds_w too long!')}
nat_meds_w[median %in% c(save_vec_min), type := 'min']
nat_meds_w[median %in% c(save_vec_max), type := 'max']
nat_meds_w <- nat_meds_w[order(WHOREGION,-type,ct,vt)]
nat_meds_w[vt==1, vt_n:='Current']
nat_meds_w[vt==2, vt_n:='Improved (minimal)']
nat_meds_w[vt==3, vt_n:='Improved (efficacy)']
nat_meds_w[vt==4, vt_n:='Improved (breadth)']
nat_meds_w[vt==5, vt_n:='Universal']
nat_meds_w[ct==1, ct_n:='0-4']
nat_meds_w[ct==2, ct_n:='0-10']
nat_meds_w[ct==3, ct_n:='0-17']
nat_meds_w[ct==4, ct_n:='65+']
nat_meds_w[ct==5, ct_n:='0-17, 65+']
nat_meds_w[, c('country_code','ct','vt'):=NULL]
nat_meds_w[, median:=round(median,2)]
nat_meds_w[, '95% CI':=paste0('(',round(eti95L,2),', ',round(eti95U,2),')')]
nat_meds_w[, c('eti50L','eti50U','eti95L','eti95U'):=NULL]
nat_meds_ww <- dcast(nat_meds_w, WHOREGION + ct_n + vt_n ~ type, value.var=c('median','95% CI'))
setnames(nat_meds_ww, 'WHOREGION','WHO Region')
nat_meds_ww <- nat_meds_ww[,c('WHO Region','ct_n','vt_n','median_min','95% CI_min','median_max','95% CI_max')]
write_csv(nat_meds_ww, paste0("econ/",scenario_name,econ_folder_name,"_plots/min_max_WHO.csv"))

nat_meds <- nat_meds[gdp_data_filt, on='country_code']

ggplot(nat_meds[measure=='median']) +
  geom_point(aes(x=gdpcap, y=threshold_price, 
                 col=as.factor(vt)),alpha=1,shape=1,size=1) +
  geom_line(aes(x=gdpcap, y=threshold_price, 
                 col=as.factor(vt))) +
  ylab('Median threshold cost ($)') + 
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                                'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(color = 'Vaccine type') +
  facet_grid(ct~., scales='fixed', labeller = labeller(vt = supp.labs,
                                                        ct = supp.labs.age)) +
  scale_x_log10(limits=c(262,110000),breaks=c(300,1000,3000,10000,30000,100000)) + scale_y_log10() +
  theme_bw() + xlab('GDP per capita') + theme(text=element_text(size=14)) 
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/gdp_vs_threshold.png"),
       width=35,height=25,units="cm")

## adding income-level (world bank)
print('income groups')

income <- data.table(read_csv('econ/world-bank-income-groups.csv',show_col_types=F))
setnames(income, 'Code','country_code')
setnames(income, "World Bank's income classification",'income_class')
income_filt <- income[country_code %in% nat_meds$country_code & Year==2022][, c('country_code','income_class')]
income_filt <- rbind(income_filt, data.table(country_code = c('GUF','VEN','XKX'),
                                             income_class = c('Upper-middle-income countries',
                                                              'Upper-middle-income countries',
                                                              'Upper-middle-income countries')))

## adding income groups
nat_meds <- nat_meds[income_filt, on=c('country_code')]
nat_meds$income_class <- factor(nat_meds$income_class,levels=c('Low-income countries','Lower-middle-income countries',
                                                                   'Upper-middle-income countries','High-income countries'))


incomefacet <- ggplot(nat_meds[measure=='median']) + 
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt), col=as.factor(vt))) +
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt)), outlier.shape=NULL, alpha=0) +
  ylab('Threshold cost ($)') + 
  scale_color_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type', color = 'Vaccine type') + 
  ylim(c(0,NA)) + 
  facet_wrap(income_class~., scales='free', nrow=2, labeller = labeller(WHOREGION=who_region_labs)) + 
  theme_bw() + xlab('Vaccination targets') + theme(text=element_text(size=12)); incomefacet
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/boxplot_medians_incomefacet.png"),
       width=40,height=22,units="cm")

nat_meds <- nat_meds[gdp_data_filt,on='country_code']
ggplot(nat_meds[measure=='median' & ct==2]) + 
  geom_point(aes(x=gdpcap, y=threshold_price, 
                   color=income_class)) +
  ylab('Threshold cost ($)') + xlab('GDP per capita, $') +
  scale_color_manual(values = income_colors) +
  labs(col = 'Income group') + 
  scale_y_log10(breaks=c(0.1,1,10,100,1000,10000)) + scale_x_log10(limits=c(200,110000), breaks=c(300,1000,10000,100000)) +
  facet_grid(.~vt_n, scales='free') +
  theme_bw() + theme(text=element_text(size=14)) 
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/scatter_gdp_cost.png"),
       width=42,height=15,units="cm")

nat_meds_wide <- nat_meds[,c('country_code','gdpcap','threshold_price','WHOREGION','ct','vt','ct_n','vt_n','measure')]
nat_meds_wide <- dcast(nat_meds_wide, country_code + gdpcap + WHOREGION + ct + vt + ct_n + vt_n ~ measure, value.var = 'threshold_price')
scatterCIs <- ggplot(nat_meds_wide[ct==2]) + 
  geom_point(aes(x=gdpcap, y=median, 
                 color=WHOREGION)) +
  geom_segment(aes(x=gdpcap, xend=gdpcap, y=eti95L, yend=eti95U, 
                 color=WHOREGION),alpha=0.5) +
  ylab('Threshold cost ($)') + xlab('GDP per capita ($)') +
  labs(col = 'WHO Region') + 
  scale_color_manual(values = WHO_colors, labels = who_region_labs2) +
  scale_y_log10(breaks=c(0.1,1,10,100,1000,10000),labels=c(0.1,1,10,100,1000,10000)) + 
  scale_x_log10(limits=c(200,110000), breaks=c(300,1000,10000,100000)) +
  facet_grid(.~vt_n, scales='free') +
  theme_bw() + theme(text=element_text(size=12)); scatterCIs
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/scatter_gdp_cost_CIs.png"),
       width=42,height=15,units="cm")

wtp_thresh <- data.table(read_csv('econ/outcome_calculations/data/WTP_thresholds.csv', show_col_type=F))
itzs <- data.table(read_csv('econ/outcome_calculations/data/new_clustering.csv', show_col_type=F))
setnames(wtp_thresh, 'iso3c','country_code')
setnames(itzs, 'codes','country_code')
nat_meds_wide <- nat_meds_wide[wtp_thresh, on='country_code']
nat_meds_wide <- nat_meds_wide[itzs, on='country_code',itz:=cluster_name]
ggplot(nat_meds_wide[ct==2 & vt==5]) + 
  geom_point(aes(x=gdpcap, y=median, 
                 color=itz)) +
  geom_segment(aes(x=gdpcap, xend=gdpcap, y=eti95L, yend=eti95U,
                   color=itz),alpha=0.5) +
  ylab('Threshold cost ($)') + xlab('GDP per capita ($)') +
  labs(col = 'ITZ') + 
  scale_color_manual(values = cluster_colors2) +
  scale_y_log10(breaks=c(0.1,1,10,100,1000,10000),labels=c(0.1,1,10,100,1000,10000)) +
  # scale_x_log10(limits=c(50,100000), breaks=c(30,100,1000,10000,100000)) +
  scale_x_log10() +
  facet_wrap(.~itz, scales='free', nrow=4) +
  theme_bw() + theme(text=element_text(size=12))
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/scatter_wtp_cost_CIs.png"),
       width=42,height=15,units="cm")

# ggplot(nat_meds_wide[ct==2]) + 
#   geom_point(aes(x=gdpcap, y=median, 
#                  color=gdp_prop)) +
#   geom_segment(aes(x=gdpcap, xend=gdpcap, y=eti95L, yend=eti95U, 
#                    color=gdp_prop),alpha=0.5) +
#   ylab('Threshold cost ($)') + xlab('GDP per capita ($)') +
#   labs(col = 'WTP as prop\nof GDP') +
#   scale_color_viridis() +
#   # scale_color_manual(values = cluster_colors2) +
#   scale_y_log10(breaks=c(0.1,1,10,100,1000,10000),labels=c(0.1,1,10,100,1000,10000)) +
#   scale_x_log10(limits=c(200,110000), breaks=c(300,1000,10000,100000)) +
#   facet_grid(.~vt_n, scales='free') +
#   theme_bw() + theme(text=element_text(size=12))
# ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/scatter_prop_cost_CIs.png"),
#        width=42,height=15,units="cm")

library(patchwork)
scatterCIs + incomefacet + plot_layout(nrow=2, guides='collect', heights = c(2, 2)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/Fig4.png"),
       width=35,height=30,units="cm")

## saving range of 95% CIs
nat_meds_less <- nat_meds[,c('country_code','ct','vt','income_class','threshold_price','measure')]
save_dt <- nat_meds_less[measure=='median'][, .(min = min(threshold_price)), by=c('ct','vt','income_class')]
save_vec_min <- save_dt$min
save_dt <- nat_meds_less[measure=='median'][, .(max = max(threshold_price)), by=c('ct','vt','income_class')]
save_vec_max <- save_dt$max
nat_meds_w <- dcast(nat_meds_less, country_code + ct + vt + income_class ~ measure, value.var = 'threshold_price')
nat_meds_w <- nat_meds_w[median %in% c(save_vec_min,save_vec_max)]
if(!nrow(nat_meds_w)==(length(save_vec_min) + length(save_vec_max))){print('nat_meds_w too long!')}
nat_meds_w[median %in% c(save_vec_min), type := 'min']
nat_meds_w[median %in% c(save_vec_max), type := 'max']
nat_meds_w <- nat_meds_w[order(income_class,-type,ct,vt)]
nat_meds_w[vt==1, vt_n:='Current']
nat_meds_w[vt==2, vt_n:='Improved (minimal)']
nat_meds_w[vt==3, vt_n:='Improved (efficacy)']
nat_meds_w[vt==4, vt_n:='Improved (breadth)']
nat_meds_w[vt==5, vt_n:='Universal']
nat_meds_w[ct==1, ct_n:='0-4']
nat_meds_w[ct==2, ct_n:='0-10']
nat_meds_w[ct==3, ct_n:='0-17']
nat_meds_w[ct==4, ct_n:='65+']
nat_meds_w[ct==5, ct_n:='0-17, 65+']
nat_meds_w[, c('country_code','ct','vt'):=NULL]
nat_meds_w[, median:=round(median,2)]
nat_meds_w[, '95% CI':=paste0('(',round(eti95L,2),', ',round(eti95U,2),')')]
nat_meds_w[, c('eti50L','eti50U','eti95L','eti95U'):=NULL]
nat_meds_ww <- dcast(nat_meds_w, income_class + ct_n + vt_n ~ type, value.var=c('median','95% CI'))
setnames(nat_meds_ww, 'income_class','Income group')
nat_meds_ww <- nat_meds_ww[,c('Income group','ct_n','vt_n','median_min','95% CI_min','median_max','95% CI_max')]
write_csv(nat_meds_ww, paste0("econ/",scenario_name,econ_folder_name,"_plots/min_max_income.csv"))


total_nat_sum <- total_nat_sum[income_filt,on='country_code']
income_agg <- copy(total_nat_sum)
income_agg[,c('country_code','WHOREGION','country'):=NULL]
income_agg <- income_agg[,lapply(.SD,sum),by=c('income_class','ct','ct_n','vt','vt_n','simulation_index')]
income_agg[,threshold_price:=(incremental_costs+incremental_DALYs)/discounted_doses]
income_agg$income_class <- factor(income_agg$income_class,levels=c('Low-income countries','Lower-middle-income countries',
                                                                   'Upper-middle-income countries','High-income countries'))
income_agg$ct_n <- factor(income_agg$ct_n,levels=unique(income_agg$ct_n))
ggplot(income_agg) + 
  geom_boxplot(aes(x=ct_n, y=threshold_price, 
                   fill=as.factor(vt))) +
  ylab('Threshold cost ($)') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  ylim(c(0,NA)) + 
  facet_wrap(income_class~., scales='free', nrow=4) + 
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/income_boxplot.png"),
       width=20,height=28,units="cm")

### HEATMAPS ###
nat_map <- nat_meds[measure=='median']
nat_map <- nat_map[income_filt, on='country_code']

classes <- unique(nat_map$income_class)
classes_n <- c('Low-income','Lower-middle-income',
                  'Upper-middle-income','High-income')

plot1 <- ggplot(nat_map[income_class==classes[1]]) + 
  ggtitle(paste(classes_n[1])) +
  geom_tile(aes(x=country_code, y=vt_n, fill=threshold_price)) +
  coord_flip() +
  facet_grid(ct~., labeller = labeller(vt = supp.labs,
                                       ct = supp.labs.age)) +
  scale_fill_viridis(option='A', trans='log', 
                     limits=c(0.01,5000), breaks=c(1,10,100,1000,5000)) +
  theme_bw() + xlab('') + labs(fill='Cost ($)') +
  ylab('') +
  theme(text = element_text(size=14),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot2 <- ggplot(nat_map[income_class==classes[2]]) + 
  ggtitle(paste(classes_n[2])) +
  geom_tile(aes(x=country_code, y=vt_n, fill=threshold_price)) +
  coord_flip() +
  facet_grid(ct~., labeller = labeller(vt = supp.labs,
                                       ct = supp.labs.age)) +
  scale_fill_viridis(option='A', trans='log', 
                     limits=c(0.01,5000), breaks=c(1,10,100,1000,5000)) +
  theme_bw() + xlab('') + labs(fill='Cost ($)') +
  ylab('') +
  theme(text = element_text(size=14),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot3 <- ggplot(nat_map[income_class==classes[3]]) + 
  ggtitle(paste(classes_n[3])) +
  geom_tile(aes(x=country_code, y=vt_n, fill=threshold_price)) +
  coord_flip() +
  facet_grid(ct~., labeller = labeller(vt = supp.labs,
                                       ct = supp.labs.age)) +
  scale_fill_viridis(option='A', trans='log', 
                     limits=c(0.01,5000), breaks=c(1,10,100,1000,5000)) +
  theme_bw() + xlab('') + labs(fill='Threshold cost ($)') +
  ylab('') +
  theme(text = element_text(size=14),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot4 <- ggplot(nat_map[income_class==classes[4]]) + 
  ggtitle(paste(classes_n[4])) + 
  geom_tile(aes(x=country_code, y=vt_n, fill=threshold_price)) +
  coord_flip() +
  facet_grid(ct~., labeller = labeller(vt = supp.labs,
                                       ct = supp.labs.age)) +
  scale_fill_viridis(option='A', trans='log', 
                     limits=c(0.01,5000), breaks=c(1,10,100,1000,5000)) +
  theme_bw() + xlab('') + labs(fill='Threshold cost ($)') +
  ylab('') +
  theme(text = element_text(size=14),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

library(patchwork)

plot1 + plot3 + plot2 + plot4 + plot_layout(guides='collect',nrow=1)

ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/heatmap_income_x_ct.png"),
       width=33,height=30,units="cm")

# ########################################################
# ### CHANGING WTP thresholds ###
# 
# load(paste0('econ/',scenario_name,'/wtp_thresh10'))
# 
# ggplot(wtp_thresh10[ct==5&vt==5&simulation_index==1]) + 
#   geom_line(aes(x=wtp_prop, y=threshold_price, group=country_code, color=country_code))
# 
# wtp_thresh10 <- wtp_thresh10[income_filt, on='country_code']
# itzs <- data.table(read_csv('econ/outcome_calculations/data/ITZs.csv', show_col_types = F))
# setnames(itzs, 'codes','country_code')
# setnames(itzs, 'cluster_name','itz')
# wtp_thresh10 <- wtp_thresh10[itzs[,c('country_code','itz')], on=c('country_code')]
# 
# wtp_thresh10_med <- wtp_thresh10[,lapply(.SD,median),by=c('country_code','ct','vt','wtp_prop','income_class','itz')]
# wtp_thresh10_med[, simulation_index:=NULL]
# classes <- unique(wtp_thresh10_med$income_class)
# 
# rm(hmps)
# for(i in 1:length(classes)){
#   plot_i <- ggplot(wtp_thresh10_med[income_class==classes[i]]) + 
#     ggtitle(paste(classes[i])) + 
#     geom_tile(aes(x=country_code, y=wtp_prop, fill=threshold_price)) +
#     coord_flip() +
#     facet_grid(ct~vt, labeller = labeller(vt = supp.labs,
#                                           ct = supp.labs.age)) +
#     scale_fill_viridis(option='A', trans='log', 
#                        limits=c(0.01,20000), breaks=c(1,10,100,1000,5000,20000)) +
#     theme_bw() + xlab('') + labs(fill='Threshold price') +
#     scale_y_continuous(breaks=1:5/5) + 
#     ylab('WTP threshold as proportion of GDP per capita') +
#     theme(text = element_text(size=14),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
#   
#   if(i == 1){
#     hmps <- plot_i
#   }else{
#     hmps <- hmps + plot_i
#   }
#   
#   plot_i 
#   ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/wtp_prop_", i,".png"),
#          width=30,height=26,units="cm")
#   
#   print(classes[i])
# }
# 
# hmps + plot_layout(guides='collect')
# ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/wtp_prop_all.png"),
#        width=45,height=40,units="cm")
# 
# ## probability of cost-effectiveness
# wtp_thresh10_prob <- copy(wtp_thresh10)
# ## is this right? guess for now
# wtp_thresh10_prob[grepl('Low-in',income_class), est_cost := 5]
# wtp_thresh10_prob[grepl('Lower-mid',income_class), est_cost := 10]
# wtp_thresh10_prob[grepl('Upper-middle',income_class), est_cost := 20]
# wtp_thresh10_prob[grepl('High',income_class), est_cost := 30]
# wtp_thresh10_prob[, over := (threshold_price>=est_cost)]
# wtp_thresh10_prob <- wtp_thresh10_prob[,lapply(.SD,sum),by=c('country_code','ct','vt','wtp_prop','income_class','itz','est_cost')]
# wtp_thresh10_prob[,c('simulation_index','threshold_price'):=NULL]
# wtp_thresh10_prob[, prob := over/100]
# 
# rm(hmps_prob)
# for(i in 1:length(classes)){
#   plot_i <- ggplot(wtp_thresh10_prob[income_class==classes[i]]) + 
#     ggtitle(paste0(classes[i])) + 
#     geom_tile(aes(x=country_code, y=wtp_prop, fill=prob)) +
#     coord_flip() +
#     facet_grid(ct~vt, labeller = labeller(vt = supp.labs,
#                                           ct = supp.labs.age)) +
#     scale_fill_viridis(option='A', 
#                        limits=c(0,1)) +
#     theme_bw() + xlab('') + labs(fill='Probability of \ncost-effectiveness') +
#     ylab('WTP threshold as proportion of GDP per capita') +
#     theme(text = element_text(size=14),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
#   
#   if(i == 1){
#     hmps_prob <- plot_i
#   }else{
#     hmps_prob <- hmps_prob + plot_i
#   }
#   
#   # plot_i 
#   # ggsave(paste0("econ/",scenario_name,"_plots/wtp_prob_", i,".png"),
#   #        width=30,height=26,units="cm")
#   
#   print(classes[i])
# }
# 
# wtp_thresh10_prob$income_class <- factor(wtp_thresh10_prob$income_class, levels=unique(wtp_thresh10_prob$income_class))
# ggplot(wtp_thresh10_prob[wtp_prop==1]) + 
#   # ggtitle(paste0(classes[i])) + 
#   geom_tile(aes(x=country_code, y=vt, fill=prob)) +
#   coord_flip() +
#   facet_grid(income_class~ct, labeller = labeller(vt = supp.labs,
#                                         ct = supp.labs.age)) +
#   scale_fill_viridis(option='mako',
#                      limits=c(0,1)) +
#   theme_bw() + xlab('') + labs(fill='Probability of \ncost-effectiveness') +
#   ylab('WTP threshold as proportion of GDP per capita') +
#   theme(text = element_text(size=14),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# wtp_heat <- wtp_thresh10_prob[wtp_prop==1]
# for(income_class_i in unique(wtp_heat$income_class)){
#   k <- nrow(wtp_heat[income_class==income_class_i,])
#   wtp_heat[income_class==income_class_i,num:=k/25]
# }
# wtp_heat_sort <- wtp_heat[ct==1&vt==1][order(income_class)]
# index <- 0
# id_num <- 0
# for(i in 1:nrow(wtp_heat_sort)){
#   if(index == wtp_heat_sort[i,]$num){
#     id_num <- id_num + 1
#   }else{
#     id_num <- 1
#   }
#   wtp_heat_sort[i, id := id_num]
#   index <- wtp_heat_sort[i,]$num
# }
# wtp_heat <- wtp_heat[wtp_heat_sort,on='country_code',id:=id]
# wtp_heat[, id_prop:=id/num]
# 
# ggplot(wtp_heat) + 
#   # ggtitle(paste0(classes[i])) + 
#   geom_tile(aes(x=(id), y=vt, fill=prob)) +
#   coord_flip() +
#   facet_grid(income_class~ct, labeller = labeller(vt = supp.labs,
#                                                   ct = supp.labs.age)) +
#   scale_fill_viridis(option='mako',
#                      limits=c(0,1)) +
#   theme_bw() + xlab('') + labs(fill='Probability of \ncost-effectiveness') +
#   ylab('WTP threshold as proportion of GDP per capita') +
#   theme(text = element_text(size=14))
# 
# 
# wtp_prob <- data.table(ct = rep(1:5, each=5),
#                        vt = 1:5,
#                        wtp_prop = rep(1:10/10, each=25))
# for(i in 1:nrow(wtp_prob)){
#   df <- wtp_thresh10_prob[ct==wtp_prob[i,]$ct &
#                             vt==wtp_prob[i,]$vt &
#                             wtp_prop==wtp_prob[i,]$wtp_prop]
#   sum <- as.numeric(nrow(df[prob>0.8]))
#   wtp_prob[i, countries := sum]
# }
# 
# ggplot(wtp_prob) +
#   geom_line(aes(x=wtp_prop, y=countries, group=as.factor(vt), col=as.factor(vt))) + 
#   scale_color_manual(values=vt_colors,
#                       labels = c('Current','Improved (minimal)',
#                                   'Improved (efficacy)','Improved (breadth)','Universal')) +
#   theme_bw() + theme(text=element_text(size=14)) +
#   scale_x_continuous(breaks=1:10/10) +
#   scale_y_continuous(limits=c(0,186)) +
#   xlab('WTP threshold as proportion of GDP per capita') +
#   labs(col='Vaccine type') +
#   ylab('Countries where vaccine is > 80% 
#        probability of being cost-effective') +
#   facet_grid(ct~., scales='free',
#              labeller = labeller(vt = supp.labs,
#                                  ct = supp.labs.age)) 
# 
# clusters <- read_csv('data/new_clustering.csv', show_col_types=F)
# pop_hist_WPP_data <- data.table(read_csv('econ/outcome_calculations/data/pop_hist_WPP_data.csv', show_col_types=F))
# pops_dt <- data.table(country_code = rep(clusters$codes, each=4),
#                       age_grp = rep(1:4, 186))
# for(i in 1:nrow(pops_dt)){
#   country_names <- unlist(unname(clusters[clusters$codes==pops_dt$country_code[i], 
#                                           c('country','country_altern','country_altern_2')]))
#   if(pops_dt$country_code[i]=='TUR'){country_names <- c(country_names, 'Turkiye')}
#   if(pops_dt$age_grp[i] == 1){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 4:4])))]
#   }
#   if(pops_dt$age_grp[i] == 2){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 5:7])))]
#   }
#   if(pops_dt$age_grp[i] == 3){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 8:16])))]
#   }
#   if(pops_dt$age_grp[i] == 4){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 17:24])))]
#   }
#   pops_dt[i, tot_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 4:24])))]
# }
# 
# wtp_thresh10_prob <- wtp_thresh10_prob[pops_dt, on = 'country_code', tot_pop := tot_pop]
# global_pop <- sum(pops_dt$tot_pop)/4
# 
# wtp_pop <- data.table(ct = rep(1:5, each=5),
#                        vt = 1:5,
#                        wtp_prop = rep(1:10/10, each=25))
# for(i in 1:nrow(wtp_pop)){
#   df <- wtp_thresh10_prob[ct==wtp_pop[i,]$ct &
#                             vt==wtp_pop[i,]$vt &
#                             wtp_prop==wtp_pop[i,]$wtp_prop]
#   sum <- sum(df[prob>0.8]$tot_pop)
#   wtp_pop[i, population := sum/global_pop]
# }
# 
# ggplot(wtp_pop) +
#   geom_line(aes(x=wtp_prop, y=100*population, group=as.factor(vt), col=as.factor(vt))) + 
#   scale_color_manual(values=vt_colors,
#                      labels = c('Current','Improved (minimal)',
#                                 'Improved (efficacy)','Improved (breadth)','Universal')) +
#   theme_bw() + theme(text=element_text(size=14)) +
#   scale_x_continuous(breaks=1:5/5) +
#   scale_y_continuous(limits=c(0,NA)) +
#   xlab('WTP threshold as proportion of GDP per capita') +
#   labs(col='Vaccine type') +
#   ylab('% of global population where vaccine is
#        more than 80% probability cost-effective') +
#   facet_grid(.~ct, scales='free',
#              labeller = labeller(vt = supp.labs,
#                                  ct = supp.labs.age)) 
# 
# 
# hmps_prob + plot_layout(guides='collect')
# ggsave(paste0("econ/",scenario_name,"_plots/wtp_prob_all.png"),
#        width=45,height=40,units="cm")
# 
# inmb <- copy(total_nat_sum)
# inmb <- inmb[income_filt, on=c('country_code')]
# inmb[grepl('Low-in',income_class), est_cost := 5]
# inmb[grepl('Lower-mid',income_class), est_cost := 10]
# inmb[grepl('Upper-middle',income_class), est_cost := 20]
# inmb[grepl('High',income_class), est_cost := 30]
# inmb[,vacc_costs := discounted_doses*est_cost]
# inmb[, inmb := - vacc_costs + incremental_costs + incremental_DALYs]
# 
# inmb_med <- inmb[,c('country_code','ct','vt','ct_n','income_class','inmb')]
# inmb_med <- inmb_med[,lapply(.SD,median),by=c('country_code','ct','vt','ct_n','income_class')]
# 
# inmb_med$income_class <- factor(inmb_med$income_class, levels=unique(inmb_med$income_class))
# ggplot(inmb_med) + 
#   geom_boxplot(aes(x=ct_n, y=inmb/1000000000, 
#                    fill=as.factor(vt))) +
#   ylab('Incremental net monetary benefit (billions)') + 
#   scale_fill_manual(values=vt_colors,
#                     labels = c('Current','Improved (minimal)',
#                                'Improved (efficacy)','Improved (breadth)','Universal')) + 
#   labs(fill = 'Vaccine type') + 
#   facet_grid(income_class~., scales='free') + 
#   theme(text=element_text(size=14)) +
#   theme_bw() + xlab('Vaccination targets')
# 
# inmb_agg <- inmb[,c('income_class','simulation_index','ct','vt','ct_n','vacc_costs','incremental_costs','incremental_DALYs')]
# inmb_agg <- inmb_agg[,lapply(.SD,sum), by=c('income_class','simulation_index','ct','vt','ct_n')]
# inmb_agg[, inmb := - vacc_costs + incremental_costs + incremental_DALYs]
# 
# inmb_agg$ct_n <- factor(inmb_agg$ct_n, levels=unique(inmb_agg$ct_n))
# inmb_agg$income_class <- factor(inmb_agg$income_class, levels=c("Low-income countries", "Lower-middle-income countries", "Upper-middle-income countries", "High-income countries"   ))
# ggplot(inmb_agg) + 
#   geom_boxplot(aes(x=ct_n, y=inmb/1000000000, 
#                    fill=as.factor(vt))) +
#   ylab('Incremental net monetary benefit (billions)') + 
#   scale_fill_manual(values=vt_colors,
#                     labels = c('Current','Improved (minimal)',
#                                'Improved (efficacy)','Improved (breadth)','Universal')) + 
#   labs(fill = 'Vaccine type') + 
#   geom_hline(yintercept=0,lty=2,alpha=0.5) +
#   facet_grid(income_class~., scales='free') + 
#   theme(text=element_text(size=14)) +
#   theme_bw() + xlab('Vaccination targets')
# ggsave(paste0("econ/",scenario_name,"_plots/INMB.png"),
#        width=25,height=20,units="cm")
# 
# # ggplot(inmb_agg) + 
# #   geom_point(aes(x = inmb/1000000000,
# #                  y=(vacc_costs - incremental_costs)/1000000000, 
# #                    col=as.factor(vt))) +
# #   ylab('Incremental costs (billions)') + 
# #   scale_color_manual(values=vt_colors,
# #                     labels = c('Current','Improved (minimal)',
# #                                'Improved (efficacy)','Improved (breadth)','Universal')) + 
# #   labs(col = 'Vaccine type') + 
# #   facet_grid(ct_n~income_class, scales='free') + 
# #   theme(text=element_text(size=14)) +
# #   theme_bw() + xlab('Incremental DALYs (billions)')
# # ggsave(paste0("econ/",scenario_name,"_plots/INMB_point.png"),
# #        width=25,height=20,units="cm")
# 


########################################################
### AGE-SPECIFIC OUTCOMES ###
print('age-specific outcomes')

load(paste0('econ/',scenario_name,econ_folder_name,'/as_out_summ_c'))

as_out_summ_c[variable == 'non_fevers', variable := 'non-fever cases']
as_out_summ_c[variable == 'total_DALYs', variable := 'total DALYs']
as_out_summ_c$age_grp_n <- factor(as_out_summ_c$age_grp_n, levels = unique(as_out_summ_c$age_grp_n))
age_name_vec <- c('0-4','0-10','0-17','65+','0-17, 65+')
for(i in 1:5){
# using ct=2 as it is likely to be the recommendation of the paper?
ggplot(as_out_summ_c[ct==i]) +
  geom_col(aes(x=age_grp_n, y=median/(30*1000000), fill=vt, group=vt),
           position='dodge', col=1) +
  geom_errorbar(aes(x = age_grp_n, ymin = eti95L/(30*1000000), ymax = eti95U/(30*1000000), group=vt),
                position=position_dodge(0.9)) +
  facet_wrap(variable~., scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age),
             nrow=3) +
  ggtitle(paste0('Outcome incidence, under ', age_name_vec[i],' age-targeting strategy')) + 
  ylab('Mean annual incidence (millions)') +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(fill = 'Vaccine type') +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Age group')

ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/age_spec_outcomes_",i,".png"),
       width=25,height=14,units="cm")
}

as_outcomes <- ggplot(as_out_summ_c[!grepl('fever',variable) & !ct=='v']) +
  geom_col(aes(x=age_grp_n, y=median/(30*1000000), fill=vt, group=vt),
           position='dodge', col=1) +
  geom_errorbar(aes(x = age_grp_n, ymin = eti95L/(30*1000000), ymax = eti95U/(30*1000000), group=vt),
                position=position_dodge(0.9)) +
  facet_grid(variable~ct, scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Mean annual incidence (millions)') +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(fill = 'Vaccine type') + theme_bw() + 
  theme(text=element_text(size=14)) +
  xlab('Age group'); as_outcomes
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/age_spec_outcomes.png"),
       width=25,height=14,units="cm")

load(paste0('econ/',scenario_name,econ_folder_name,'/as_out_averted'))
as_out_summ_c$age_grp_n <- factor(as_out_summ_c$age_grp_n, levels = unique(as_out_summ_c$age_grp_n))
as_out_summ_c[, variable := gsub('_av','',variable)]
as_out_summ_c$variable <- factor(as_out_summ_c$variable, levels = unique(as_out_summ_c$variable))
av_outcomes <- ggplot(as_out_summ_c[!grepl('fever',variable)]) +
  geom_col(aes(x=age_grp_n, y=median/(30*1000000), fill=vt, group=vt),
           position='dodge', col=1) +
  geom_errorbar(aes(x = age_grp_n, ymin = eti95L/(30*1000000), ymax = eti95U/(30*1000000), group=vt),
                position=position_dodge(0.9)) +
  facet_grid(variable~ct, scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Mean averted annual outcomes (millions)') +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(fill = 'Vaccine type') + theme_bw() + 
  theme(text=element_text(size=14)) +
  xlab('Age group'); av_outcomes
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/age_spec_averted.png"),
       width=32,height=18,units="cm")

load(paste0('econ/',scenario_name,econ_folder_name,'/as_out_averted_WHOREGION'))
as_out_summ_c$age_grp_n <- factor(as_out_summ_c$age_grp_n, levels = unique(as_out_summ_c$age_grp_n))
as_out_summ_c[, variable := gsub('_av','',variable)]
as_out_summ_c$variable <- factor(as_out_summ_c$variable, levels = unique(as_out_summ_c$variable))
for(region in unique(as_out_summ_c$WHOREGION)){
  ggplot(as_out_summ_c[!grepl('fever',variable) & WHOREGION==region]) +
    geom_col(aes(x=age_grp_n, y=median/(30*1000000), fill=vt, group=vt),
             position='dodge', col=1) +
    ggtitle(region) + 
    geom_errorbar(aes(x = age_grp_n, ymin = eti95L/(30*1000000), ymax = eti95U/(30*1000000), group=vt),
                  position=position_dodge(0.9)) +
    facet_grid(variable~ct, scales='free',
               labeller = labeller(vt = supp.labs,
                                   ct = supp.labs.age)) +
    ylab('Mean averted annual outcomes (millions)') +
    scale_fill_manual(values=vt_colors,
                      labels = c('Current','Improved (minimal)',
                                 'Improved (efficacy)','Improved (breadth)','Universal')) +
    labs(fill = 'Vaccine type') + theme_bw() + 
    theme(text=element_text(size=14)) +
    xlab('Age group')
  ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/age_spec_averted_",region,".png"),
         width=32,height=20,units="cm")
}


##### NNV ##### 

load(paste0('data/vacc_output_base/global_weekly.Rdata'))
global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "scenario")]
global_weekly_m[, week:=NULL]
global_weekly_nv <- global_weekly_m[scenario=='no_vacc',]
global_weekly_m[global_weekly_nv, on = c("simulation_index"), nv_inf := i.value]
global_weekly_m[,averted := nv_inf-value]
global_weekly_m[,nv_inf:=NULL]
global_weekly_m <- global_weekly_m[!scenario=='no_vacc']
global_weekly_m[, scen := 'base']
agg_out_short <- copy(global_weekly_m)
vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_base/vacc_doses_", c_code, "_",
                                                  "base.csv"),show_col_type=F) %>% 
                        mutate(cluster_code = c_code))
  
}
vacc_doses_g <- vacc_doses %>% dplyr::select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  dplyr::select(!c(name,year)) %>% 
  group_by(vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))
vacc_doses_g <- data.table(vacc_doses_g)
setnames(vacc_doses_g, 'vacc_type','vt'); setnames(vacc_doses_g, 'age_cov','ct')
setnames(vacc_doses_g, 'value','doses')
agg_out_short[, ct:=as.numeric(substr(scenario,4,4))]
agg_out_short[, vt:=as.numeric(substr(scenario,9,9))]
agg_out_short <- agg_out_short[vacc_doses_g[, c('ct','vt','doses')], on=c('ct','vt')]
agg_out_short[, nnv:=doses/averted]
agg_out_short[ct==1, ct_n:='0-4']
agg_out_short[ct==2, ct_n:='0-10']
agg_out_short[ct==3, ct_n:='0-17']
agg_out_short[ct==4, ct_n:='65+']
agg_out_short[ct==5, ct_n:='0-17, 65+']
agg_out_short$ct_n <- factor(agg_out_short$ct_n,levels=unique(agg_out_short$ct_n))

nnv_outcomes <- ggplot(agg_out_short) +
  # geom_boxplot(aes(x=ct_n, y=nnv, 
  #                  fill=as.factor(vt)), outlier.colour=NA) +
  # geom_jitter(aes(x=ct_n, y=nnv, group=as.factor(vt), color=as.factor(vt)), size=0.4, alpha=0.9) +
  geom_point(aes(x=vt, y=nnv, col=as.factor(vt)),
             position=position_jitterdodge(dodge.width=0.9,jitter.width=2.5), shape=4, alpha=0.7) +
  # geom_boxplot(aes(x=ct_n, y=nnv, fill=as.factor(vt)), alpha=0.1, outlier.colour=NA, position=position_dodge(width=0.9)) +
  # facet_grid(.~ct, scales='free', labeller = labeller(ct = supp.labs.age)) +
  theme_bw() + ylab('Number needed to vaccinate') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type') +
  # xlab('Age-targeting strategy') +
  # geom_hline(yintercept=1, lty=2, alpha=0.3) +
  # scale_y_continuous(breaks=0:8) +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10), limits=c(0.1,10)) +
  facet_grid(.~ct_n, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  xlab('') +
  # geom_hline(yintercept=1, lty=2, alpha=0.3) +
  # scale_y_continuous(limits=c(0,250)) +
  # scale_y_log10(limits=c(3,250)) +
  theme(text=element_text(size=14),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank());nnv_outcomes

nnv_outcomes_box <- ggplot(agg_out_short) +
  # geom_boxplot(aes(x=ct_n, y=nnv,
  #                  fill=as.factor(vt)), outlier.colour=NA) +
  # geom_jitter(aes(x=ct_n, y=nnv, group=as.factor(vt), color=as.factor(vt)), size=0.4, alpha=0.9) +
  geom_boxplot(aes(x=vt, y=nnv, fill=as.factor(vt))) +
  # geom_boxplot(aes(x=ct_n, y=nnv, fill=as.factor(vt)), alpha=0.1, outlier.colour=NA, position=position_dodge(width=0.9)) +
  # facet_grid(.~ct, scales='free', labeller = labeller(ct = supp.labs.age)) +
  theme_bw() + ylab('Number needed to vaccinate') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type') +
  xlab('Age-targeting strategy') +
  # geom_hline(yintercept=1, lty=2, alpha=0.3) +
  # scale_y_continuous(breaks=0:8) +
  facet_grid(.~ct_n, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10), limits=c(0.1,10)) +
  theme(text=element_text(size=14),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank());nnv_outcomes_box


# nnv_outcomes3 <- ggplot(agg_out_short) +
#   geom_density(aes(x=nnv, 
#                    fill=as.factor(vt)), alpha=0.5,) +
#   facet_grid(ct~., scales='free', labeller = labeller(ct = supp.labs.age)) +
#   theme_bw() + ylab('Density') +
#   scale_fill_manual(values=vt_colors, labels = supp.labs) +
#   labs(fill = 'Vaccine type') + 
#   xlab('Number needed to vaccinate') +
#   scale_x_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10), limits=c(0.1,10)) +
#   geom_vline(xintercept=1,lty=2,alpha=0.5) +
#   theme(text=element_text(size=14),
#         legend.position='none');nnv_outcomes3

# agg_out_meas <- data.table()
# for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
#   dt <- agg_out_short[,c('scenario','ct','vt','ct_n','nnv')][, lapply(.SD, get(meas)), by = c('scenario','ct','vt','ct_n')]
#   dt[, measure := meas]
#   agg_out_meas <- rbind(agg_out_meas, dt)
# }
# agg_out_meas[,simulation_index:=NULL]
# agg_out_meas_w <- dcast(agg_out_meas, scenario + ct + ct_n + vt ~ measure, value.var='nnv')
# agg_out_meas_w$ct_n <- factor(agg_out_meas_w$ct_n, levels=unique(agg_out_meas_w$ct_n))
# nnv_outcomes2 <- ggplot(agg_out_meas_w) +
#   geom_point(aes(x=ct_n, y=median,
#                    col=as.factor(vt)), size=2, position = position_dodge(0.2)) +
#   geom_errorbar(aes(x=ct_n,
#                    ymin=eti95L, ymax=eti95U, 
#                  col=as.factor(vt)), width=0.1, lwd=0.8, position = position_dodge(0.2)) +
#   # facet_grid(.~ct_n, scales='fixed') + 
#   theme_bw() + ylab('Number needed to vaccinate') +
#   scale_fill_manual(values=vt_colors, labels = supp.labs) +
#   scale_color_manual(values=vt_colors, labels = supp.labs) +
#   labs(col = 'Vaccine type') + 
#   xlab('Age-targeting strategy') +
#   coord_flip() +
#   # scale_y_continuous(breaks=0:8) +
#   scale_y_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10), limits=c(0.1,10)) +
#   theme(text=element_text(size=14),
#         legend.position='none'); nnv_outcomes2

### nnv DALYs
agg_out[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
agg_out[, nnvd := incremental_DALYs/discounted_doses]

nnvd_plot <- ggplot(agg_out) +
  # geom_point(aes(x=ct_n, y=nnvd, col=as.factor(vt)),
  # position=position_jitterdodge(dodge.width=0.9,jitter.width=1), shape=4) +
  geom_boxplot(aes(x=vt, y=nnvd, fill=as.factor(vt))) +
  theme_bw() + ylab('DALYs averted per dose given') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type') +
  facet_grid(.~ct_n, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  xlab('') +
  # geom_hline(yintercept=1, lty=2, alpha=0.3) +
  # scale_y_continuous(limits=c(0,250)) +
  # scale_y_log10(limits=c(3,250)) +
  theme(text=element_text(size=14),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank());nnvd_plot
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/dalys_per_dose.png"),
       width=20,height=15,units="cm")


nnv_outcomes + nnvd_plot + av_outcomes + plot_layout(nrow=3, guides='collect', heights = c(2, 2, 3)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/Fig3.png"),
       width=35,height=36,units="cm")

layout <- "
AAABBB
AAABBB
CCCCCC
CCCCCC
CCCCCC
"
nnv_outcomes + nnvd_plot +
    av_outcomes + plot_layout(design = layout, guides='collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/Fig3-1.png"),
       width=35,height=36,units="cm")
  

### nnv DALYs by WHO region
nnvd_WHO <- total_nat_sum[, c('WHOREGION','ct','ct_n','vt','simulation_index','incremental_DALYs','discounted_doses')]
nnvd_WHO <- nnvd_WHO[, lapply(.SD, sum), by=c('WHOREGION','ct','ct_n','vt','simulation_index')]
nnvd_WHO[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
nnvd_WHO[, nnvd := incremental_DALYs/discounted_doses]

nnvd_WHO_plot <- ggplot(nnvd_WHO) +
  # geom_point(aes(x=ct_n, y=nnvd, col=as.factor(vt)),
  # position=position_jitterdodge(dodge.width=0.9,jitter.width=1), shape=4) +
  geom_boxplot(aes(x=ct_n, y=nnvd, fill=as.factor(vt))) +
  theme_bw() + ylab('DALYs averted per dose given') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type') +
  facet_wrap(.~WHOREGION, scales='free',labeller = labeller(WHOREGION=who_region_labs),nrow=3) + 
  xlab('Age-targeting strategy') +
  theme(text=element_text(size=14));nnvd_WHO_plot
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/WHO_dalys_per_dose.png"),
       width=35,height=30,units="cm")


# ########################################################
# ### PROPORTION OF COUNTRIES/POPULATION ###
# print('proportion of countries/population')
# 
# price_prop <- nat_meds[measure=='median']
# clusters <- read_csv('data/new_clustering.csv', show_col_types=F)
# pop_hist_WPP_data <- data.table(read_csv('econ/outcome_calculations/data/pop_hist_WPP_data.csv', show_col_types=F))
# pops_dt <- data.table(country_code = rep(clusters$codes, each=4),
#                       age_grp = rep(1:4, 186))
# for(i in 1:nrow(pops_dt)){
#   country_names <- unlist(unname(clusters[clusters$codes==pops_dt$country_code[i],
#                                           c('country','country_altern','country_altern_2')]))
#   if(pops_dt$country_code[i]=='TUR'){country_names <- c(country_names, 'Turkiye')}
#   if(pops_dt$age_grp[i] == 1){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 4:4])))]
#   }
#   if(pops_dt$age_grp[i] == 2){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 5:7])))]
#   }
#   if(pops_dt$age_grp[i] == 3){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 8:16])))]
#   }
#   if(pops_dt$age_grp[i] == 4){
#     pops_dt[i, age_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 17:24])))]
#   }
#   pops_dt[i, tot_pop := 1000*sum(unlist(unname(pop_hist_WPP_data[name %in% country_names & Year == 2020, 4:24])))]
# }
# price_prop <- price_prop[pops_dt, on = 'country_code', tot_pop := tot_pop]
# 
# max_price <- max(price_prop$threshold_price)
# min_price <- min(price_prop$threshold_price)
# global_pop <- sum(pops_dt$tot_pop)/4
# prices <- data.table(price = rep((ceiling(min_price):ceiling(max_price)), 25),
#                      ct = rep(rep(1:5, each = (ceiling(max_price))),5),
#                      vt = rep(1:5, each = 5*(ceiling(max_price))))
# 
# pb <- txtProgressBar(min = 0, max = nrow(prices), style = 3, width = 50, char = "=")
# for(i in 1:nrow(prices)){
#   dt_i <- price_prop[ct==prices[i,ct] & vt==prices[i,vt] 
#                      & threshold_price>=prices[i,price],]
#   no_countries_i <- nrow(dt_i)
#   pop <- sum(dt_i$tot_pop)
#   prices[i, no_countries := no_countries_i]
#   prices[i, prop_pop := pop/global_pop]
#   setTxtProgressBar(pb, i)
# }
# 
# ggplot(prices) + 
#   geom_line(aes(x=price, y=no_countries, col=as.factor(vt)),lwd=0.8) +
#   facet_grid(ct~., labeller = labeller(vt = supp.labs,
#                                        ct = supp.labs.age)) +
#   theme_bw() +
#   scale_color_manual(values=vt_colors,
#                      labels = c('Current','Improved (minimal)',
#                                 'Improved (efficacy)','Improved (breadth)','Universal')) + 
#   labs(color = 'Vaccine type') + 
#   scale_y_continuous(breaks=c(0,50,100,150,186)) +
#   theme(text=element_text(size=14)) +
#   ylab('Number of countries cost-effective') +
#   xlab('Threshold cost ($)') + scale_x_log10(breaks=c(1,3,10,30,100,300,1000,3000))
# ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/countries_cet.png"),
#        width=25,height=20,units="cm")
# 
# ggplot(prices) + 
#   geom_line(aes(x=price, y=prop_pop, col=as.factor(vt)),lwd=0.8) +
#   facet_grid(ct~., labeller = labeller(vt = supp.labs,
#                                        ct = supp.labs.age)) +
#   theme_bw() +
#   scale_color_manual(values=vt_colors,
#                      labels = c('Current','Improved (minimal)',
#                                 'Improved (efficacy)','Improved (breadth)','Universal')) + 
#   labs(color = 'Vaccine type') + 
#   theme(text=element_text(size=14)) +
#   ylab('Proportion of global population cost-effective') +
#   xlab('Threshold cost ($)') + scale_x_log10(breaks=c(1,3,10,30,100,300,1000,3000))
# ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/population_cet.png"),
#        width=25,height=20,units="cm")









