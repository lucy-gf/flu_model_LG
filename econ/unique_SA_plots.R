#### Economic analysis - sensitivity analysis plots
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

#### BASE VS BREADTH VS DEPTH ####

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_BASE50/vacc_doses_", c_code, "_",
                                                  "BASE50.csv"),show_col_type=F) %>% 
                        mutate(cluster_code = c_code))
  
}
vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  select(!c(name,year)) %>% 
  group_by(vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))

wastage <- 0.1

###
vacc_print <- data.table(vacc_doses_g)
vacc_print[vacc_type==1, vt:='Current']
vacc_print[vacc_type==2, vt:='Improved (minimal)']
vacc_print[vacc_type==3, vt:='Improved (efficacy)']
vacc_print[vacc_type==4, vt:='Improved (breadth)']
vacc_print[vacc_type==5, vt:='Universal']
vacc_print[age_cov==1, ct:='0-4']
vacc_print[age_cov==2, ct:='0-10']
vacc_print[age_cov==3, ct:='0-17']
vacc_print[age_cov==4, ct:='65+']
vacc_print[age_cov==5, ct:='0-17, 65+']
vacc_print[,doses := signif(value/((1-wastage)*30*1000000), 3)]
write_csv(vacc_print,'data/vacc_doses_BASE50/vacc_print.csv')
###

vacc_doses_g <- rbind(cbind(vacc_doses_g,scenario='depth'),
                      cbind(vacc_doses_g,scenario='base'),
                      cbind(rbind(vacc_doses_g[vacc_doses_g$vacc_type==1,],
                                  vacc_doses_g[vacc_doses_g$vacc_type==1,],
                                  vacc_doses_g[vacc_doses_g$vacc_type==1,],
                                  vacc_doses_g[vacc_doses_g$vacc_type==1,],
                                  vacc_doses_g[vacc_doses_g$vacc_type==1,]), scenario='breadth'))
vacc_doses_g <- data.table(vacc_doses_g)
vacc_doses_g[scenario=='breadth',]$vacc_type <- rep(1:5, each=5)
setnames(vacc_doses_g, 'vacc_type','vt'); setnames(vacc_doses_g, 'age_cov','ct'); setnames(vacc_doses_g, 'value','doses');
vacc_doses_g[, vacc_program := NULL]

load(paste0('econ/BASE50/as_out_averted_100'))
comparison <- cbind(out_av[,c('ct','vt','simulation_index','infections_av')], scenario='base')
load(paste0('econ/BASE50_DEPTH/as_out_averted_100'))
comparison <- rbind(comparison, cbind(out_av[,c('ct','vt','simulation_index','infections_av')], scenario='depth'))
load(paste0('econ/BASE50_BREADTH/as_out_averted_100'))
comparison <- rbind(comparison, cbind(out_av[,c('ct','vt','simulation_index','infections_av')], scenario='breadth'))
comparison <- comparison[!ct=='v']
comparison[, ct := as.numeric(ct)][, vt := as.numeric(vt)]

comparison <- comparison[vacc_doses_g, on = c('ct','vt','scenario')]
comparison[, nnv := (doses/(1-wastage))/(infections_av)]

comparison[ct==1, ct_n:='0-4']
comparison[ct==2, ct_n:='0-10']
comparison[ct==3, ct_n:='0-17']
comparison[ct==4, ct_n:='65+']
comparison[ct==5, ct_n:='0-17, 65+']

comparison_cis <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- comparison[, lapply(.SD, get(meas)), by = c('scenario','ct','vt','ct_n')]
  dt[, measure := meas]
  comparison_cis <- rbind(comparison_cis, dt)
}
comparison_cis[,simulation_index:=NULL]
comparison_cis_w <- dcast(comparison_cis[,c('scenario','ct','vt','ct_n','nnv','measure')], scenario+ct+vt+ct_n~measure, value.var='nnv')
comparison_cis_w[vt==1, vt_n:='Current']
comparison_cis_w[vt==2, vt_n:='Improved (minimal)']
comparison_cis_w[vt==3, vt_n:='Improved (efficacy)']
comparison_cis_w[vt==4, vt_n:='Improved (breadth)']
comparison_cis_w[vt==5, vt_n:='Universal']
comparison_cis_w$vt_n <- factor(comparison_cis_w$vt_n, levels=unique(comparison_cis_w$vt_n))
comparison_cis_w$ct_n <- factor(comparison_cis_w$ct_n, levels=unique(comparison_cis_w$ct_n))
comparison_cis_w[scenario=='base', scenario:='Base']
comparison_cis_w[scenario=='breadth', scenario:='Breadth']
comparison_cis_w[scenario=='depth', scenario:='Depth']

ggplot(comparison_cis_w) + 
  geom_line(aes(x=vt_n, y=median, group=scenario,
                 col=scenario),lwd=0.8) +
  geom_ribbon(aes(x=vt, ymin=eti95L, ymax=eti95U, fill = scenario), alpha=0.2) +
  geom_ribbon(aes(x=vt, ymin=eti50L, ymax=eti50U, fill = scenario), alpha=0.4) +
  ylab('Number needed to vaccinate') + 
  labs(color = 'Scenario') + 
  labs(fill = 'Scenario') + 
  facet_grid(.~ct_n, scales='fixed') +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10,30), labels=c(0.1,0.3,1,3,10,30)) +
  theme_bw() + xlab('') +
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette='Dark2') +
  scale_fill_brewer(palette='Dark2')
ggsave(paste0("econ/BASE50_DEPTH_plots/depth_breadth.png"),
       width=30,height=14,units="cm")



#### BASE VS REL INF ####

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_BASE50/vacc_doses_", c_code, "_",
                                                  "BASE50.csv"),show_col_type=F) %>% 
                        mutate(cluster_code = c_code))
  
}
vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  select(!c(name,year)) %>% 
  group_by(vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))

vacc_doses_g <- rbind(cbind(vacc_doses_g,scenario='base'),
                      cbind(vacc_doses_g,scenario='rel_inf'))
vacc_doses_g <- data.table(vacc_doses_g)
setnames(vacc_doses_g, 'vacc_type','vt'); setnames(vacc_doses_g, 'age_cov','ct'); setnames(vacc_doses_g, 'value','doses');
vacc_doses_g[, vacc_program := NULL]

load(paste0('econ/BASE50/as_out_averted_100'))
comparison <- cbind(out_av[,c('ct','vt','simulation_index','infections_av')], scenario='base')
load(paste0('econ/BASE50_REL_INF/as_out_averted_100'))
comparison <- rbind(comparison, cbind(out_av[,c('ct','vt','simulation_index','infections_av')], scenario='rel_inf'))
comparison <- comparison[!ct=='v']
comparison[, ct := as.numeric(ct)][, vt := as.numeric(vt)]

wastage <- 0.1

comparison <- comparison[vacc_doses_g, on = c('ct','vt','scenario')]
comparison[, nnv := (doses/(1-wastage))/(infections_av)]

comparison[ct==1, ct_n:='0-4']
comparison[ct==2, ct_n:='0-10']
comparison[ct==3, ct_n:='0-17']
comparison[ct==4, ct_n:='65+']
comparison[ct==5, ct_n:='0-17, 65+']

comparison_cis <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- comparison[, lapply(.SD, get(meas)), by = c('scenario','ct','vt','ct_n')]
  dt[, measure := meas]
  comparison_cis <- rbind(comparison_cis, dt)
}
comparison_cis[,simulation_index:=NULL]
comparison_cis_w <- dcast(comparison_cis[,c('scenario','ct','vt','ct_n','nnv','measure')], scenario+ct+vt+ct_n~measure, value.var='nnv')
comparison_cis_w[vt==1, vt_n:='Current']
comparison_cis_w[vt==2, vt_n:='Improved (minimal)']
comparison_cis_w[vt==3, vt_n:='Improved (efficacy)']
comparison_cis_w[vt==4, vt_n:='Improved (breadth)']
comparison_cis_w[vt==5, vt_n:='Universal']
comparison_cis_w$vt_n <- factor(comparison_cis_w$vt_n, levels=unique(comparison_cis_w$vt_n))
comparison_cis_w$ct_n <- factor(comparison_cis_w$ct_n, levels=unique(comparison_cis_w$ct_n))
comparison_cis_w[scenario=='base', scenario:='Base']
comparison_cis_w[scenario=='rel_inf', scenario:='Reduced relative \ninfectiousness']

ggplot(comparison_cis_w) + 
  geom_line(aes(x=vt_n, y=median, group=scenario,
                col=scenario),lwd=0.8) +
  geom_ribbon(aes(x=vt, ymin=eti95L, ymax=eti95U, fill = scenario), alpha=0.2) +
  geom_ribbon(aes(x=vt, ymin=eti50L, ymax=eti50U, fill = scenario), alpha=0.4) +
  ylab('Number needed to vaccinate') + 
  labs(color = 'Scenario') + 
  labs(fill = 'Scenario') + 
  facet_grid(.~ct_n, scales='fixed') +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10,30), labels=c(0.1,0.3,1,3,10,30)) +
  theme_bw() + xlab('') +
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette='Dark2') +
  scale_fill_brewer(palette='Dark2')
ggsave(paste0("econ/BASE50_REL_INF_plots/rel_inf_nnv.png"),
       width=30,height=14,units="cm")


#### 20% VS 50% VS 70% ####
load(paste0('econ/BASE50/as_out_averted_100'))
cov_comp <- cbind(out_av, coverage=0.5)
load(paste0('econ/LOW20/as_out_averted_100'))
cov_comp <- rbind(cov_comp, cbind(out_av, coverage=0.2))
load(paste0('econ/base/as_out_averted_100'))
cov_comp <- rbind(cov_comp, cbind(out_av, coverage=0.7))
cov_comp <- cov_comp[!ct=='v']

cov_comp_s <- cov_comp[, c('ct','vt','coverage','infections_av','hospitalisations_av','deaths_av','total_DALYs_av')]
cov_comp_l <- melt(cov_comp_s, id.vars = c('ct','vt','coverage'))
cov_comp_median <- cov_comp_l[, lapply(.SD, median), by=c('vt','ct','coverage','variable')]
setnames(cov_comp_median, 'value','median')
l95 <- cov_comp_l[, lapply(.SD, eti95L), by=c('vt','ct','coverage','variable')]
u95 <- cov_comp_l[, lapply(.SD, eti95U), by=c('vt','ct','coverage','variable')]
cov_comp_median[, l95 := l95$value][, u95 := u95$value]
cov_comp_median[variable == 'infections_av', variable := 'Infections']
cov_comp_median[variable == 'hospitalisations_av', variable := 'Hospitalisations']
cov_comp_median[variable == 'deaths_av', variable := 'Deaths']
cov_comp_median[variable == 'total_DALYs_av', variable := 'Total DALYs']

ggplot(cov_comp_median) +
  geom_ribbon(aes(x=100*coverage, ymin=l95/(30*1000000), ymax=u95/(30*1000000), 
                fill=vt), alpha=0.2) +
  geom_line(aes(x=100*coverage, y=median/(30*1000000), col=vt, group=vt), lwd=0.8) +
  geom_point(aes(x=100*coverage, y=median/(30*1000000), col=vt, group=vt)) +
  facet_grid(variable~ct, scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Global annual incidence averted (millions)') +
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                                'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(col = 'Vaccine type', fill = 'Vaccine type', shape = 'Coverage level') + theme_bw() +
  theme(text=element_text(size=14)) +
  xlab('') + ylim(c(0,NA)) + scale_x_continuous(breaks=c(20,50,70), limits=c(10,80))

ggplot(cov_comp_median) +
  geom_bar(aes(x=as.factor(100*coverage), y=median/(30*1000000), fill=vt, group=vt),
               position='dodge', stat='identity', col='black') +
  geom_errorbar(aes(x=as.factor(100*coverage), ymin=l95/(30*1000000), ymax=u95/(30*1000000), group=vt), 
                position=position_dodge(0.9)) +
  facet_grid(variable~ct, scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Global annual incidence averted (millions)') +
  # scale_color_manual(values=rep('black',5),
  #                    labels = c('Current','Improved (minimal)',
  #                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(col = 'Vaccine type', fill = 'Vaccine type', shape = 'Coverage level') + theme_bw() +
  theme(text=element_text(size=16)) +
  xlab('Vaccination coverage (%)') + ylim(c(0,NA))
ggsave(paste0("econ/LOW20_plots/coverage_vs_averted.png"),
       width=35,height=20,units="cm")

## saving data for main text table 2
coverage_table <- cov_comp_median[ct==2]
coverage_table[vt==1, vt_n:='Current']
coverage_table[vt==2, vt_n:='Improved (minimal)']
coverage_table[vt==3, vt_n:='Improved (efficacy)']
coverage_table[vt==4, vt_n:='Improved (breadth)']
coverage_table[vt==5, vt_n:='Universal']
coverage_table[, c('ct','vt'):=NULL]
coverage_table <- coverage_table[!grepl('DALY',variable)]
coverage_table[, median := median/30][, l95 := l95/30][, u95 := u95/30]
coverage_table[variable=='Deaths', data:=paste0(signif(median,3),
                                                ' (',signif(l95,3),
                                                ', ',signif(u95,3),')')]
coverage_table[variable=='Hospitalisations', data:=paste0(signif(median/1000000,3),
                                                ' (',signif(l95/1000000,3),
                                                ', ',signif(u95/1000000,3),')')]
coverage_table[variable=='Hospitalisations',variable:='Hospitalisations (millions)']
coverage_table[variable=='Infections', data:=paste0(signif(median/1000000000,3),
                                                          ' (',signif(l95/1000000000,3),
                                                          ', ',signif(u95/1000000000,3),')')]
coverage_table[variable=='Infections',variable:='Infections (billions)']
coverage_table[,coverage := paste0(100*coverage, '%')]
coverage_table[, c('median','l95','u95'):=NULL]
coverage_table_w <- dcast(coverage_table, 
                          factor(variable, levels=unique(variable)) + coverage ~ factor(vt_n, levels=unique(vt_n)), 
                          value.var = 'data')
write_csv(coverage_table_w, file='econ/LOW20_plots/coverage_table.csv')

#### COVERAGE VS NNV ####
load('econ/LOW20/agg_out')
low20 <- agg_out[simulation_index==1,][, c('ct','vt','doses')]
load('econ/BASE50/agg_out')
base50 <- agg_out[simulation_index==1][, c('ct','vt','doses')]
load('econ/base/agg_out')
high70 <- agg_out[simulation_index==1][, c('ct','vt','doses')]
doses <- data.table(rbind(cbind(low20, coverage=0.2),
                          cbind(base50, coverage=0.5),
                          cbind(high70, coverage=0.7)))
  
cov_nnv <- cov_comp[, c('ct','vt','coverage','infections_av','hospitalisations_av','deaths_av','total_DALYs_av')]
cov_nnv <- melt(cov_nnv, id.vars=c('ct','vt','coverage'))
cov_nnv <- cov_nnv[doses, on=c('ct','vt','coverage')]
cov_nnv[variable == 'infections_av', variable := 'Infections']
cov_nnv[variable == 'hospitalisations_av', variable := 'Hospitalisations']
cov_nnv[variable == 'deaths_av', variable := 'Deaths']
cov_nnv[variable == 'total_DALYs_av', variable := 'Total DALYs']

cov_nnv_plot <- ggplot(cov_nnv[ct==2]) +
  # geom_point(aes(x=as.factor(100*coverage), y=doses/value, col=vt, group=vt),
  #            position=position_jitterdodge(dodge.width=0.9,jitter.width=0.2), shape=4, alpha=0.5) +
  geom_boxplot(aes(x=as.factor(100*coverage), y=doses/value, fill=vt, col=vt,
                   group=interaction(as.factor(vt),as.factor(coverage))),
               width=1.1) +
  geom_boxplot(aes(x=as.factor(100*coverage), y=doses/value, fill=vt,
                   group=interaction(as.factor(vt),as.factor(coverage))),
               alpha=0, outliers=F, width=1.1) +
  # geom_errorbar(aes(x=as.factor(100*coverage), ymin=doses/l95, ymax=doses/u95, 
  #                   col = vt, group=vt), 
  #               position=position_dodge(0.9)) +
  facet_grid(variable~coverage, scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Number needed to vaccinate') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                                'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(col = 'Vaccine type', fill = 'Vaccine type', shape = 'Coverage level') + theme_bw() +
  theme(text=element_text(size=14),
        strip.background.x= element_blank(),
        strip.text.x = element_blank()) +
  xlab('Vaccination coverage (%)'); cov_nnv_plot
ggsave(paste0("econ/LOW20_plots/coverage_vs_nnv.png"),
       width=25,height=25,units="cm")

ggplot(cov_nnv[ct==2]) +
  # geom_point(aes(x=as.factor(100*coverage), y=doses/value, col=vt, group=vt),
  #            position=position_jitterdodge(dodge.width=0.9,jitter.width=0.2), shape=4, alpha=0.5) +
  geom_boxplot(aes(x=as.factor(100*coverage), y=doses/value, fill=vt, col=vt,
                   group=interaction(as.factor(vt),as.factor(coverage))),
               width=1) +
  geom_boxplot(aes(x=as.factor(100*coverage), y=doses/value, fill=vt,
                   group=interaction(as.factor(vt),as.factor(coverage))),
               alpha=0, outliers=F, width=1) +
  # geom_errorbar(aes(x=as.factor(100*coverage), ymin=doses/l95, ymax=doses/u95, 
  #                   col = vt, group=vt), 
  #               position=position_dodge(0.9)) +
  facet_grid(variable~vt, scales='free',
             labeller = labeller(vt = supp.labs,
                                 ct = supp.labs.age)) +
  ylab('Number needed to vaccinate') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_manual(values=vt_colors,
                     labels = c('Current','Improved (minimal)',
                                'Improved (efficacy)','Improved (breadth)','Universal')) +
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) +
  labs(col = 'Vaccine type', fill = 'Vaccine type', shape = 'Coverage level') + theme_bw() +
  theme(text=element_text(size=14),
        strip.background.x= element_blank(),
        strip.text.x = element_blank()) +
  xlab('Vaccination coverage (%)') + scale_y_log10()

cov_nnv_plot + scale_y_log10()
ggsave(paste0("econ/LOW20_plots/coverage_vs_nnv_log.png"),
       width=35,height=20,units="cm")

cov_nnv[, per := value/doses]
cov_nnv_med <- cov_nnv[, 1000*median(per), by=c('ct','vt','coverage','variable')]
cov_nnv_med <- cov_nnv_med[order(ct,vt,variable,coverage)]
write_csv(cov_nnv_med, file='econ/LOW20_plots/outcomes_av_per_dose.csv')

#### INCIDENCE PLOTS #####
incidence_ex <- data.table()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  data <- data.table(readRDS(paste0('data/vacc_output_base/vacc_', c_code,'_none_ct_1.rds'))[[1]])
  data <- data[country_code==c_code]
  data[, c('country'):=NULL]
  data <- melt(data, id.vars = c('country_code','simulation_index','week'))
  data[, strain := paste0('Influenza ', substr(variable, 4,4))]
  data[, variable:=NULL]
  data <- data[, lapply(.SD, sum), by=c('country_code','simulation_index','week','strain')]
  incidence_ex <- rbind(incidence_ex, data)
  print(c_code)
}

ggplot(incidence_ex[simulation_index<=10]) + 
  geom_line(aes(x=week,y=value/1000000,group=simulation_index,col=strain),
            lwd=0.6, alpha=0.3) +
  facet_grid(country_code~strain, scales='free', 
             labeller = labeller(country_code=supp.labs.country)) +
  theme_bw() + xlab('') + ylab('Infections (millions)')  +
  theme(text=element_text(size=14),
        legend.position='none') + 
  scale_color_manual(values=c('#7d66ac', '#e483a4'))

ggsave(paste0("SM_plots/incidence_ex.png"),
       width=30,height=25,units="cm")




