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

#### SETTING SCENARIO ####

scenario_name <- 'BASE50'
econ_folder_name <- ''
print(paste0(scenario_name,econ_folder_name))

#### GLOBAL ####

print('global')
load(paste0('econ/',scenario_name,econ_folder_name,'/agg_out'))

agg_out$ct_n <- factor(agg_out$ct_n, levels = unique(agg_out$ct_n))

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


#### WHO REGION ####
print('WHO regions')

load(paste0('econ/',scenario_name,econ_folder_name,'/total_who_sum'))

total_who_sum$ct_n <- factor(total_who_sum$ct_n, levels = unique(total_who_sum$ct_n))

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

#### NATIONAL THRESHOLDS ####
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
                   fill=as.factor(vt), col=as.factor(vt))) +
  geom_boxplot(aes(x=ct_n, y=threshold_price, group=interaction(ct_n,vt)), fill = NA, outliers=F) +
  # geom_jitter(aes(x=ct_n, y=threshold_price), alpha = 0.3) +
  ylab('Threshold price ($)') +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_color_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(fill = 'Vaccine type', color = 'Vaccine type') +
  # ylim(c(0,NA)) +
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
  ylab('Threshold price ($)') +
  scale_color_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(fill = 'Vaccine type', color = 'Vaccine type') +
  # ylim(c(0,NA)) +
  facet_wrap(income_class~., scales='free', nrow=2, labeller = labeller(WHOREGION=who_region_labs)) +
  theme_bw() + xlab('Vaccination targets') + theme(text=element_text(size=12)); incomefacet
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/boxplot_medians_incomefacet.png"),
       width=40,height=22,units="cm")

nat_meds <- nat_meds[gdp_data_filt,on='country_code']

zero_val <- 0.01 # here, less than one cent = 0

nat_meds_wide <- nat_meds[,c('country_code','gdpcap','threshold_price','WHOREGION','ct','vt','ct_n','vt_n','measure')]
nat_meds_wide <- dcast(nat_meds_wide, country_code + gdpcap + WHOREGION + ct + vt + ct_n + vt_n ~ measure, value.var = 'threshold_price')
scatterCIs <- ggplot() +
  geom_hline(yintercept = zero_val, lty=2) + 
  geom_segment(data = nat_meds_wide[ct==2 & eti95L>zero_val],
               aes(x=gdpcap, xend=gdpcap, y=eti95L, yend=eti95U,
                 color=WHOREGION),alpha=0.4) +
  geom_segment(data = nat_meds_wide[ct==2 & median>zero_val & eti95L<=zero_val],
               aes(x=gdpcap, xend=gdpcap, y=median, yend=eti95U,
                   color=WHOREGION),alpha=0.4) +
  geom_segment(data = nat_meds_wide[ct==2 & median>zero_val & eti95L<=zero_val],
               aes(x=gdpcap, xend=gdpcap, y=zero_val, yend=eti95U,
                   color=WHOREGION),alpha=0.25) +
  geom_point(data = nat_meds_wide[ct==2 & median<=zero_val],
             aes(x=gdpcap, y=zero_val,
                 color=WHOREGION)) +
  geom_point(data = nat_meds_wide[ct==2 & median>zero_val],
             aes(x=gdpcap, y=median,
                 color=WHOREGION)) +
  ylab('Threshold price ($)') + xlab('GDP per capita ($)') +
  labs(col = 'WHO Region') +
  scale_color_manual(values = WHO_colors, labels = who_region_labs2) +
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),labels=c('< 0.01',0.1,1,10,100,1000,10000)) +
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

## adding number of countries not cost-effective
not_ce <- nat_meds_less[measure=='median']
not_ce[, c('country_code','measure'):=NULL]
not_ce <- not_ce[, paste0(round(100*nrow(.SD[threshold_price>0])/nrow(.SD)), '%'), by=c('income_class','ct','vt')]
nat_meds_w <- nat_meds_w[not_ce, on=c('ct','vt','income_class')]

nat_meds_w[, c('country_code','ct','vt'):=NULL]
nat_meds_w[, data:=paste0(signif(median,2),' (',signif(eti95L,2),', ',signif(eti95U,2),')')]
nat_meds_w[, c('eti50L','eti50U','eti95L','eti95U'):=NULL]
nat_meds_ww <- dcast(nat_meds_w, income_class + ct_n + vt_n + V1 ~ type, value.var=c('data'))
nat_meds_ww <- nat_meds_ww[,c('income_class','ct_n','vt_n','min','max', 'V1')]
nat_meds_ww <- nat_meds_ww[order(income_class,factor(ct_n, levels=c('0-4','0-10','0-17','65+','0-17, 65+')),
                                 factor(vt_n, levels=c('Current','Improved (minimal)',
                                                       'Improved (efficacy)','Improved (breadth)','Universal')))]
setnames(nat_meds_ww, 'income_class','Income group')
setnames(nat_meds_ww, 'ct_n','Age-targeting strategy')
setnames(nat_meds_ww, 'vt_n','Vaccine type')
setnames(nat_meds_ww, 'min','Minimum price ($)')
setnames(nat_meds_ww, 'max','Maximum price ($)')
setnames(nat_meds_ww, 'V1','Countries cost-effective')
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
                   fill=as.factor(vt), color=as.factor(vt))) +
  geom_boxplot(aes(x=ct_n, y=threshold_price, fill=as.factor(vt),
                   group=interaction(ct_n,vt)), outlier.shape=NULL, alpha=0) +
  ylab('Threshold price ($)') +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_color_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(fill = 'Vaccine type', color = 'Vaccine type') +
  ylim(c(0,NA)) +
  facet_wrap(income_class~., scales='free', nrow=4) +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/income_boxplot.png"),
       width=20,height=28,units="cm")


#### AGE-SPECIFIC OUTCOMES ####
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
  ylab('Annual incidence (millions)') +
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
  ylab('Annual incidence (millions)') +
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

## writing csv of total outcomes averted
save_csv <- as_out_summ_c[,!c('age_grp','age_grp_n')]
save_csv <- save_csv[, lapply(.SD, sum), by=c('vt','ct','variable')]
save_csv[,c('eti50L','eti50U'):=NULL]
save_csv <- save_csv[variable %in% c('infections','hospitalisations','deaths')]
save_csv[variable=='infections', scaling_factor:=1000000000]
save_csv[variable=='hospitalisations', scaling_factor:=1000000]
save_csv[variable=='deaths', scaling_factor:=1]
save_csv[, median := median/(30*scaling_factor)][, eti95L := eti95L/(30*scaling_factor)][, eti95U := eti95U/(30*scaling_factor)]
save_csv[, out := paste0(signif(median,3), ' (', signif(eti95L,3), ', ', signif(eti95U,3), ')')]
save_csv_w <- dcast(save_csv[,c('variable','ct','vt','out')], formula = ct + vt ~ variable, value.var = 'out')
save_csv_w[ct==1, ct_n:='0-4']
save_csv_w[ct==2, ct_n:='0-10']
save_csv_w[ct==3, ct_n:='0-17']
save_csv_w[ct==4, ct_n:='65+']
save_csv_w[ct==5, ct_n:='0-17, 65+']
save_csv_w[vt==1, vt_n:='Current']
save_csv_w[vt==2, vt_n:='Improved (minimal)']
save_csv_w[vt==3, vt_n:='Improved (efficacy)']
save_csv_w[vt==4, vt_n:='Improved (breadth)']
save_csv_w[vt==5, vt_n:='Universal']
save_csv_w[,c('ct','vt'):=NULL]
write_csv(save_csv_w[,c('ct_n','vt_n','infections','hospitalisations','deaths')], file=paste0("econ/",scenario_name,econ_folder_name,"_plots/outcomes_averted.csv"))


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

#### NNV #####
print('nnv')

wastage <- 0.1

total_nat_nnv <- copy(total_nat_sum)
no_vacc_infs <- data.table()
scen_name <- ifelse(scenario_name %like% 'disease|gdp|outpatient|discount',
                    (strsplit(scenario_name,'_')[[1]]), scenario_name)
for(c_code in c('ARG','AUS','CAN','CHN','GBR','GHA','TUR')){
  load(paste0('data/vacc_output_', scen_name,'/econ_inp_', c_code, '.Rdata'))
  econ_inp <- econ_inp[scenario=='no_vacc']
  econ_inp[, scenario := NULL]
  econ_inp_m <- melt(econ_inp, id.vars=c('country_code','simulation_index','year'))
  econ_inp_m[, c('variable','year') := NULL]
  econ_inp_m <- econ_inp_m[, lapply(.SD, sum), by = c('country_code','simulation_index')]
  no_vacc_infs <- rbind(no_vacc_infs, econ_inp_m)
  print(paste0(c_code,' done'))
}
setnames(no_vacc_infs, 'value', 'no_vacc_infections')

total_nat_nnv <- total_nat_nnv[no_vacc_infs, on = c('country_code','simulation_index')]
total_nat_nnv[, nnv := (doses/(1-wastage))/(no_vacc_infections - infections)]

total_nat_nnv[ct==1, ct_n:='0-4']
total_nat_nnv[ct==2, ct_n:='0-10']
total_nat_nnv[ct==3, ct_n:='0-17']
total_nat_nnv[ct==4, ct_n:='65+']
total_nat_nnv[ct==5, ct_n:='0-17, 65+']

total_global_nnv <- total_nat_nnv[,c('ct_n','vt','simulation_index','infections','no_vacc_infections','doses')]
total_global_nnv <- total_global_nnv[, lapply(.SD, sum), by=c('ct_n','vt','simulation_index')]
total_global_nnv[, nnv := (doses/(1-wastage))/(no_vacc_infections - infections)]

total_global_nnv$ct_n <- factor(total_global_nnv$ct_n, levels=unique(total_global_nnv$ct_n))

nnv_outcomes <- ggplot(total_global_nnv) +
  geom_point(aes(x=vt, y=nnv, col=as.factor(vt)),
             position=position_jitterdodge(dodge.width=0.9,jitter.width=2.5), shape=4, alpha=0.7) +
  theme_bw() + ylab('Number needed to vaccinate') +
  # scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type') +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10), limits=c(0.08,10)) +
  facet_grid(.~ct_n, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  xlab('') +
  theme(text=element_text(size=14),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank());nnv_outcomes

## infections under each ct x vt
cum_infs <- rbind(total_global_nnv,
                  data.table(ct_n = 'none', vt=0, simulation_index=1:100,
                             infections=total_global_nnv[ct_n=='0-4'&vt==1]$no_vacc_infections,
                             no_vacc_infections=total_global_nnv[ct_n=='0-4'&vt==1]$no_vacc_infections,
                             doses=0,nnv=0))
cum_infs[, c('no_vacc_infections','doses','nnv'):=NULL]
cum_meds <- data.table()
for(meas in c('median', 'eti95L', 'eti95U')){
  dt <- cum_infs[, lapply(.SD, get(meas)), by = c('ct_n','vt')]
  dt[, measure := meas]
  cum_meds <- rbind(cum_meds, dt)
}
cum_meds[vt==1, vt_n:='Current']
cum_meds[vt==2, vt_n:='Improved (minimal)']
cum_meds[vt==3, vt_n:='Improved (efficacy)']
cum_meds[vt==4, vt_n:='Improved (breadth)']
cum_meds[vt==5, vt_n:='Universal']
cum_meds_csv <- cum_meds[, c('ct_n','vt_n','measure','infections')]
cum_meds_csv <- dcast(cum_meds_csv, ct_n + factor(vt_n, levels=unique(vt_n)) ~ measure, value.var = 'infections')
cum_meds_csv[, column := paste0(signif(median/(30*1000000),3),
                                  ' (',signif(eti95L/(30*1000000),3),
                                  ', ',signif(eti95U/(30*1000000),3),')')]
cum_meds_csv <- cum_meds_csv[, c('ct_n','vt_n','column')]
setnames(cum_meds_csv, 'ct_n', 'Age-targeting strategy')
setnames(cum_meds_csv, 'vt_n', 'Vaccine type')
setnames(cum_meds_csv, 'column', 'Annual infections (millions)')
write_csv(cum_meds_csv, file=paste0("econ/",scenario_name,econ_folder_name,"_plots/cumulative_infections.csv"))

## proportion of infections averted by each ct x vt
total_global_nnv[, averted := no_vacc_infections - infections]
total_global_nnv[, prop_averted := (no_vacc_infections - infections)/no_vacc_infections]
prop_av <- data.table()
for(meas in c('median', 'eti95L', 'eti95U')){
  dt <- total_global_nnv[, lapply(.SD, get(meas)), by = c('ct_n','vt')]
  dt[, measure := meas]
  prop_av <- rbind(prop_av, dt)
}
prop_av[vt==1, vt_n:='Current']
prop_av[vt==2, vt_n:='Improved (minimal)']
prop_av[vt==3, vt_n:='Improved (efficacy)']
prop_av[vt==4, vt_n:='Improved (breadth)']
prop_av[vt==5, vt_n:='Universal']
prop_av_csv <- prop_av[, c('ct_n','vt_n','measure','averted','prop_averted')]
prop_av_csv <- dcast(prop_av_csv, ct_n + factor(vt_n, levels=unique(vt_n)) ~ measure, value.var = c('averted','prop_averted'))
prop_av_csv[, av_column := paste0(signif(averted_median/(30*1000000),3),
                               ' (',signif(averted_eti95L/(30*1000000),3),
                               ', ',signif(averted_eti95U/(30*1000000),3),')')]
prop_av_csv[, prop_column := paste0(round(prop_averted_median*100),'%')]
prop_av_csv <- prop_av_csv[, c('ct_n','vt_n','av_column','prop_column')]
setnames(prop_av_csv, 'ct_n', 'Age-targeting strategy')
setnames(prop_av_csv, 'vt_n', 'Vaccine type')
setnames(prop_av_csv, 'av_column', 'Annual infections averted (millions)')
setnames(prop_av_csv, 'prop_column', 'Median proportion of infections averted')
write_csv(prop_av_csv, file=paste0("econ/",scenario_name,econ_folder_name,"_plots/infections_averted.csv"))

## write doses.csv
doses_csv <- total_global_nnv[simulation_index==1,][, c('ct_n','vt','doses')]
doses_csv[vt==1, vt_n:='Current']
doses_csv[vt==2, vt_n:='Improved (minimal)']
doses_csv[vt==3, vt_n:='Improved (efficacy)']
doses_csv[vt==4, vt_n:='Improved (breadth)']
doses_csv[vt==5, vt_n:='Universal']
doses_csv[,vt:=NULL]
doses_csv[, doses := signif(doses/(30*1000000),3)]
write_csv(doses_csv, file=paste0("econ/",scenario_name,econ_folder_name,"_plots/doses.csv"))


### nnv by WHO region
nnv_who <- total_nat_nnv[, c('country_code','WHOREGION','ct','vt','simulation_index','nnv')]
nnv_who_med <- nnv_who[, lapply(.SD, median), by=c('country_code','WHOREGION','ct','vt')]

ggplot(nnv_who_med) +
  geom_boxplot(aes(x=vt, y=nnv, fill=as.factor(vt), col=as.factor(vt))) +
  geom_boxplot(aes(x=vt, y=nnv, group=as.factor(vt),
                   fill = as.factor(vt)), outliers = F, alpha=0) +
  theme_bw() + ylab('Number needed to vaccinate') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type', col = 'Vaccine type') +
  facet_grid(WHOREGION~ct, scales='fixed',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age,
                                 WHOREGION=who_region_labs)) +
  xlab('') +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10)) +
  theme(text=element_text(size=14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/NNV_by_WHO.png"),
       width=30,height=34,units="cm")

### nnv DALYs
agg_out[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
agg_out[, nnvd := (base_DALYs - discounted_DALYs)/discounted_doses]

nnvd_plot <- ggplot(agg_out) +
  # geom_point(aes(x=ct_n, y=nnvd, col=as.factor(vt)),
  # position=position_jitterdodge(dodge.width=0.9,jitter.width=1), shape=4) +
  geom_boxplot(aes(x=vt, y=nnvd, fill=as.factor(vt), col=as.factor(vt))) +
  geom_boxplot(aes(x=vt, y=nnvd, group=as.factor(vt)), fill=NA, outliers = F) +
  theme_bw() + ylab('DALYs averted per dose given') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  # labs(fill = 'Vaccine type') +
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
nnvd_WHO <- total_nat_sum[, c('country_code','WHOREGION','ct','ct_n','vt','simulation_index','base_DALYs','discounted_DALYs','discounted_doses')]
nnvd_WHO[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
nnvd_WHO[, nnvd := (base_DALYs - discounted_DALYs)/discounted_doses]
nnvd_WHO_med <- nnvd_WHO[, c('country_code','WHOREGION','ct','ct_n','vt','nnvd')][, lapply(.SD, median), by=c('country_code','WHOREGION','ct','ct_n','vt')]

nnvd_WHO_plot <- ggplot(nnvd_WHO_med) +
  # geom_point(aes(x=ct_n, y=nnvd, col=as.factor(vt)),
  # position=position_jitterdodge(dodge.width=0.9,jitter.width=1), shape=4) +
  geom_boxplot(aes(x=ct_n, y=nnvd, fill=as.factor(vt), col=as.factor(vt))) +
  geom_boxplot(aes(x=ct_n, y=nnvd, group=interaction(ct_n, as.factor(vt)), fill=as.factor(vt)), outlier.shape=NULL, alpha=0) +
  theme_bw() + ylab('DALYs averted per dose given') +
  scale_fill_manual(values=vt_colors, labels = supp.labs) +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  labs(fill = 'Vaccine type', color = 'Vaccine type') +
  facet_wrap(.~WHOREGION, scales='fixed',labeller = labeller(WHOREGION=who_region_labs),nrow=3) +
  xlab('Age-targeting strategy') + 
  theme(text=element_text(size=14));nnvd_WHO_plot
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/WHO_dalys_per_dose.png"),
       width=35,height=30,units="cm")

### deaths prevented per dose

vacc_doses_m <- total_nat_sum[simulation_index==1][,c('ct','vt','doses')]
vacc_doses_m <- vacc_doses_m[, lapply(.SD, sum), by=c('ct','vt')]
vacc_doses_m[, ct := as.numeric(ct)][, vt := as.numeric(vt)]

load(paste0("econ/",scenario_name,econ_folder_name,"/as_out_averted_100"))
out_av <- out_av[!ct=='v']
out_av[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
out_av <- out_av[vacc_doses_m, on=c('ct','vt')]
deaths_dt <- out_av[, c('simulation_index','ct','vt','infections_av','deaths_av','hospitalisations_av','doses')]
deaths_dt_l <- melt(deaths_dt, id.vars = c('simulation_index','ct','vt','doses'))

ggplot(deaths_dt_l[variable %in% c('deaths_av','hospitalisations_av')]) +
  geom_point(aes(x=vt, y=doses/value, col=as.factor(vt)),
             position=position_jitterdodge(dodge.width=0.9,jitter.width=2.5), shape=4, alpha=0.7) +
  # geom_boxplot(aes(x=vt, y=doses/value, fill=as.factor(vt))) +
  theme_bw() + ylab('Number needed to vaccinate') +
  scale_color_manual(values=vt_colors, labels = supp.labs) +
  # scale_fill_manual(values=vt_colors, labels = supp.labs) +
  labs(col = 'Vaccine type', fill = 'Vaccine type') +
  scale_y_log10() +
  # scale_y_log10(breaks=c(0.1,0.3,1,3,10), labels=c(0.1,0.3,1,3,10), limits=c(0.1,10)) +
  facet_grid(variable~ct, scales='free',
             labeller = labeller(variable = var_labs,
                                 ct = supp.labs.age)) +
  xlab('') +
  theme(text=element_text(size=14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(paste0("econ/",scenario_name,econ_folder_name,"_plots/NNV_deaths_hosp.png"),
       width=26,height=15,units="cm")

# ## number of outcomes under no vaccinations
# load(paste0("econ/",scenario_name,econ_folder_name,"/as_out_averted_100"))
# nov <- out_av[ct=='v']
# paste0(signif(median(nov$infections)/(30*1000000000), 3), ' (', signif(eti95L(nov$infections)/(30*1000000000), 3),
#               ', ', signif(eti95U(nov$infections)/(30*1000000000), 3), ')')
# paste0(signif(median(nov$hospitalisations)/(30*1000000), 3), ' (', signif(eti95L(nov$hospitalisations)/(30*1000000), 3),
#        ', ', signif(eti95U(nov$hospitalisations)/(30*1000000), 3), ')')
# paste0(signif(median(nov$deaths)/(30), 3), ' (', signif(eti95L(nov$deaths)/(30), 3),
#        ', ', signif(eti95U(nov$deaths)/(30), 3), ')')






