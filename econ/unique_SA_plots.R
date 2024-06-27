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
load(paste0('data/vacc_output_breadth/global_weekly.Rdata'))
global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "scenario")]
global_weekly_m[, week:=NULL]
global_weekly_nv <- global_weekly_m[scenario=='no_vacc',]
global_weekly_m[global_weekly_nv, on = c("simulation_index"), nv_inf := i.value]
global_weekly_m[,averted := nv_inf-value]
global_weekly_m[,nv_inf:=NULL]
global_weekly_m <- global_weekly_m[!scenario=='no_vacc']
global_weekly_m[, scen := 'breadth']
agg_out_short <- rbind(agg_out_short, global_weekly_m)
load(paste0('data/vacc_output_depth/global_weekly.Rdata'))
global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "scenario")]
global_weekly_m[, week:=NULL]
global_weekly_nv <- global_weekly_m[scenario=='no_vacc',]
global_weekly_m[global_weekly_nv, on = c("simulation_index"), nv_inf := i.value]
global_weekly_m[,averted := nv_inf-value]
global_weekly_m[,nv_inf:=NULL]
global_weekly_m <- global_weekly_m[!scenario=='no_vacc']
global_weekly_m[, scen := 'depth']
agg_out_short <- rbind(agg_out_short, global_weekly_m)
agg_out_short[, ct:=as.numeric(substr(scenario,4,4))]
agg_out_short[, vt:=as.numeric(substr(scenario,9,9))]
agg_out_short[, scenario:=NULL]

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_base/vacc_doses_", c_code, "_",
                                                  "base.csv"),show_col_type=F) %>% 
                          mutate(cluster_code = c_code))
  
}
vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  select(!c(name,year)) %>% 
  group_by(vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))

vacc_doses_g <- rbind(cbind(vacc_doses_g,scen='depth'),
                      cbind(vacc_doses_g,scen='base'),
                        cbind(rbind(vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,],
                        vacc_doses_g[vacc_doses_g$vacc_type==1,]), scen='breadth'))
vacc_doses_g <- data.table(vacc_doses_g)
vacc_doses_g[scen=='breadth',]$vacc_type <- rep(1:5, each=5)
setnames(vacc_doses_g, 'vacc_type','vt'); setnames(vacc_doses_g, 'age_cov','ct')
agg_out_short[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
agg_out_short <- agg_out_short[vacc_doses_g, on=c('ct','vt','scen'), doses:=value]
agg_out_short[,nnv := doses/averted]
agg_out_short[ct==1, ct_n:='0-4']
agg_out_short[ct==2, ct_n:='0-10']
agg_out_short[ct==3, ct_n:='0-17']
agg_out_short[ct==4, ct_n:='65+']
agg_out_short[ct==5, ct_n:='0-17, 65+']
agg_out_short$ct_n <- factor(agg_out_short$ct_n, levels=unique(agg_out_short$ct_n))
agg_out_short[,vacc_program := NULL]

ggplot(agg_out_short) + 
  geom_boxplot(aes(x=as.factor(ct_n), y=nnv, 
                   fill=as.factor(vt))) +
  ylab('Number needed to vaccinate') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  facet_grid(scen~., scales='fixed') + 
  scale_y_log10(breaks=c(0.3,1,3,10)) +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')

agg_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- agg_out_short[, lapply(.SD, get(meas)), by = c('scen','ct','vt','ct_n')]
  dt[, measure := meas]
  agg_meas <- rbind(agg_meas, dt)
}
agg_meas[,simulation_index:=NULL]
agg_meas_w <- dcast(agg_meas[,c('scen','ct','vt','ct_n','nnv','measure')], scen+ct+vt+ct_n~measure, value.var='nnv')
agg_meas_w[vt==1, vt_n:='Current']
agg_meas_w[vt==2, vt_n:='Improved (minimal)']
agg_meas_w[vt==3, vt_n:='Improved (efficacy)']
agg_meas_w[vt==4, vt_n:='Improved (breadth)']
agg_meas_w[vt==5, vt_n:='Universal']
agg_meas_w$vt_n <- factor(agg_meas_w$vt_n, levels=unique(agg_meas_w$vt_n))
agg_meas_w[scen=='base', scen:='Base']
agg_meas_w[scen=='breadth', scen:='Breadth']
agg_meas_w[scen=='depth', scen:='Depth']

ggplot(agg_meas_w) + 
  geom_line(aes(x=vt_n, y=median, group=scen,
                 col=scen),lwd=0.8) +
  geom_ribbon(aes(x=vt, ymin=eti95L, ymax=eti95U, fill = scen), alpha=0.2) +
  geom_ribbon(aes(x=vt, ymin=eti50L, ymax=eti50U, fill = scen), alpha=0.4) +
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
ggsave(paste0("econ/depth_plots/depth_breadth.png"),
       width=30,height=14,units="cm")



############

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
load(paste0('data/vacc_output_rel_inf/global_weekly.Rdata'))
global_weekly_m <- data.table(melt(global_weekly, id.vars = c("simulation_index", "week", "scenario")))
global_weekly_m[, variable:=NULL]
global_weekly_m <- global_weekly_m[, lapply(.SD, sum, na.rm=T), c("simulation_index", "scenario")]
global_weekly_m[, week:=NULL]
global_weekly_nv <- global_weekly_m[scenario=='no_vacc',]
global_weekly_m[global_weekly_nv, on = c("simulation_index"), nv_inf := i.value]
global_weekly_m[,averted := nv_inf-value]
global_weekly_m[,nv_inf:=NULL]
global_weekly_m <- global_weekly_m[!scenario=='no_vacc']
global_weekly_m[, scen := 'rel_inf']
agg_out_short <- rbind(agg_out_short, global_weekly_m)
agg_out_short[, ct:=as.numeric(substr(scenario,4,4))]
agg_out_short[, vt:=as.numeric(substr(scenario,9,9))]
agg_out_short[, scenario:=NULL]

vacc_doses <- data.frame()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  vacc_doses <- rbind(vacc_doses, read_csv(paste0("data/vacc_doses_base/vacc_doses_", c_code, "_",
                                                  "base.csv"),show_col_type=F) %>% 
                        mutate(cluster_code = c_code))
  
}
vacc_doses_g <- vacc_doses %>% select(!c(country)) %>% pivot_longer(!c(country_code, vacc_program, year, cluster_code)) %>% 
  select(!c(name,year)) %>% 
  group_by(vacc_program) %>% summarise(value = sum(value)) %>% arrange(vacc_program) %>% 
  mutate(vacc_type = (vacc_program %% 5), age_cov = ceiling(vacc_program/5)) %>% 
  mutate(vacc_type = case_when(vacc_type == 0 ~ 5, .default = vacc_type))

vacc_doses_g <- rbind(cbind(vacc_doses_g,scen='base'),
                      cbind(vacc_doses_g,scen='rel_inf'))
vacc_doses_g <- data.table(vacc_doses_g)
setnames(vacc_doses_g, 'vacc_type','vt'); setnames(vacc_doses_g, 'age_cov','ct')
agg_out_short[, ct := as.numeric(ct)][, vt := as.numeric(vt)]
agg_out_short <- agg_out_short[vacc_doses_g, on=c('ct','vt','scen'), doses:=value]
agg_out_short[,nnv := doses/averted]
agg_out_short[ct==1, ct_n:='0-4']
agg_out_short[ct==2, ct_n:='0-10']
agg_out_short[ct==3, ct_n:='0-17']
agg_out_short[ct==4, ct_n:='65+']
agg_out_short[ct==5, ct_n:='0-17, 65+']
agg_out_short$ct_n <- factor(agg_out_short$ct_n, levels=unique(agg_out_short$ct_n))
agg_out_short[,vacc_program := NULL]

ggplot(agg_out_short) + 
  geom_boxplot(aes(x=as.factor(ct_n), y=nnv, 
                   fill=as.factor(vt))) +
  ylab('Number needed to vaccinate') + 
  scale_fill_manual(values=vt_colors,
                    labels = c('Current','Improved (minimal)',
                               'Improved (efficacy)','Improved (breadth)','Universal')) + 
  labs(fill = 'Vaccine type') + 
  facet_grid(scen~., scales='fixed') + 
  scale_y_log10(breaks=c(0.3,1,3,10)) +
  theme(text=element_text(size=14)) +
  theme_bw() + xlab('Vaccination targets')

agg_meas <- data.table()
for(meas in c('median','eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt <- agg_out_short[, lapply(.SD, get(meas)), by = c('scen','ct','vt','ct_n')]
  dt[, measure := meas]
  agg_meas <- rbind(agg_meas, dt)
}
agg_meas[,simulation_index:=NULL]
agg_meas_w <- dcast(agg_meas[,c('scen','ct','vt','ct_n','nnv','measure')], scen+ct+vt+ct_n~measure, value.var='nnv')
agg_meas_w[vt==1, vt_n:='Current']
agg_meas_w[vt==2, vt_n:='Improved (minimal)']
agg_meas_w[vt==3, vt_n:='Improved (efficacy)']
agg_meas_w[vt==4, vt_n:='Improved (breadth)']
agg_meas_w[vt==5, vt_n:='Universal']
agg_meas_w$vt_n <- factor(agg_meas_w$vt_n, levels=unique(agg_meas_w$vt_n))
agg_meas_w[scen=='base', scen:='Base']
agg_meas_w[scen=='rel_inf', scen:='Reduced\nrelative\ninfectiousness']

ggplot(agg_meas_w) + 
  geom_line(aes(x=vt_n, y=median, group=scen,
                col=scen),lwd=0.8) +
  geom_ribbon(aes(x=vt, ymin=eti95L, ymax=eti95U, fill = scen), alpha=0.2) +
  geom_ribbon(aes(x=vt, ymin=eti50L, ymax=eti50U, fill = scen), alpha=0.4) +
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
ggsave(paste0("econ/rel_inf_plots/rel_inf_nnv.png"),
       width=30,height=14,units="cm")

####### INCIDENCE PLOTS ######
incidence_ex <- data.table()
for(c_code in c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")){
  data <- data.table(readRDS(paste0('data/vacc_output_base/vacc_', c_code,'_none_ct_1.rds'))[[1]])
  data <- data[simulation_index==1 & country_code==c_code]
  data[, c('country', 'simulation_index'):=NULL]
  data <- melt(data, id.vars = c('country_code','week'))
  data[, strain := paste0('Influenza ', substr(variable, 4,4))]
  data[, variable:=NULL]
  data <- data[, lapply(.SD, sum), by=c('country_code','week','strain')]
  incidence_ex <- rbind(incidence_ex, data)
  print(c_code)
}

ggplot(incidence_ex) + 
  geom_line(aes(x=week,y=value/1000000,group=strain,col=strain),lwd=0.6) +
  facet_grid(country_code~strain, scales='free') +
  theme_bw() + xlab('') + ylab('Infections (millions)')  +
  theme(text=element_text(size=14)) +
  scale_color_manual(values=c('#7d66ac', '#e483a4'))

ggsave(paste0("SM_plots/incidence_ex.png"),
       width=30,height=25,units="cm")













