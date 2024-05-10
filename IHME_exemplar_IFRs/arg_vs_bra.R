
## ARGENTINA VS ALTERNATIVES [BRAZIL] FOR IFRs

#setwd("~/Desktop/research asst/Global Code")

# load libraries, functions
# this file loads libraries and some custom-made functions
source("fcns/fcns.R")

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

# LOAD DATA
# FLUNET
flunet_data <- read_csv("data/VIW_FNT.csv") 
# POPULATION 
cntr_pop_2020 = pop[-(1:29),] %>% select(name,`2020`) %>% rename(pop_size=`2020`) %>% 
  mutate(pop_size=round(pop_size/1e3,2)) %>% rename(COUNTRY_AREA_TERRITORY=name)

flunet_data_UK_summed <- left_join(
  flunet_data %>% 
    mutate(COUNTRY_AREA_TERRITORY=ifelse(grepl("United Kingdom",COUNTRY_AREA_TERRITORY),
                                         "United Kingdom",COUNTRY_AREA_TERRITORY),
           COUNTRY_CODE=ifelse(grepl("United Kingdom",COUNTRY_AREA_TERRITORY),"UK",COUNTRY_CODE)) %>%
    filter(COUNTRY_AREA_TERRITORY %in% gsub(", England","",unique(flu_ITZ_clusters$country_altern_name)) & 
             ISO_YEAR>=2008 & ISO_WEEKSTARTDATE<as.Date("2020-04-01")) %>%
    group_by(COUNTRY_AREA_TERRITORY,COUNTRY_CODE,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE) %>%
    summarise(INF_A=sum(INF_A),INF_B=sum(INF_B),SPEC_PROCESSED_NB=sum(SPEC_PROCESSED_NB),
              ISO_WEEKSTARTDATE=unique(ISO_WEEKSTARTDATE)) %>% ungroup(),
  flu_ITZ_clusters %>% select(country,country_altern_name) %>% 
    rename(COUNTRY_AREA_TERRITORY=country_altern_name) %>% unique() ) %>% 
  mutate(country=ifelse(COUNTRY_AREA_TERRITORY %in% "United Kingdom","United Kingdom",country))

cntr_pop_2020 %>% filter(COUNTRY_AREA_TERRITORY %in% c('Argentina','Brazil'))

sa_flu <- flunet_data_UK_summed %>% filter(COUNTRY_CODE %in% c('ARG','BRA'),
                                           ISO_YEAR %in% 2010:2019) %>% 
  pivot_longer(cols = c('INF_A','INF_B')) %>% rename(STRAIN=name) %>% 
  mutate(pop=ifelse(country=='Argentina', 45.20, 212.56))

sa_flu %>% ggplot() + 
  geom_line(aes(x=ISO_WEEKSTARTDATE, y=value, group=country, color=country), lwd=1) +
  facet_wrap(STRAIN~country, scales='free') + xlab('Week') +
  theme_bw() + ylab('Reported cases')

sa_flu %>% ggplot() + 
  geom_line(aes(x=ISO_WEEKSTARTDATE, y=value, group=country, color=country), lwd=1) +
  facet_grid(STRAIN~., scales='free') + xlab('Week') +
  theme_bw() + ylab('Reported cases')

sa_flu %>% ggplot() + 
  geom_line(aes(x=ISO_WEEKSTARTDATE, y=value/(pop*1000000), group=country, color=country), lwd=1) +
  facet_grid(STRAIN~., scales='free') + xlab('Week') +
  theme_bw() + ylab('Reported cases per person')


## i guess the question is: when you apply Argentina's flu parameters/timing to Brazil but use 
## their demography & vaccination coverage, do you get a similar shape?

# assumed coverage = c(0.70,0.1,0.15,0.73)
    # confirm/validate/source this on Tuesday
# runs in 'IHME_exemplar_IFRs/brazil_runs'

braz <- data.table(readRDS("data/original_epids_output/brazil_BRA_2010_2019_100.rds")[[1]])
braz[, c('country','country_code'):=NULL]
braz_l <- melt(braz, id.vars = c('simulation_index','week'))[,strain := substr(variable,4,4)][,variable:=NULL]
braz_l <- braz_l[,lapply(.SD,sum),by=c('simulation_index','week','strain')]

ggplot(braz_l) + 
  geom_line(aes(week, value/1000000, group=simulation_index, col=strain)) +
  facet_grid(strain~.) + scale_x_continuous(breaks=2010:2020) +
  scale_color_manual(values = strain_colors1) +
  ylab('Infections, millions') + xlab('Year') + 
  theme_minimal()
  
braz_med <- braz_l[, lapply(.SD, median, na.rm=T), by = c('week','strain')]
braz_med[, measure := 'median']
for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt_m2 <- braz_l[, lapply(.SD, get(meas)), by = c('week','strain')]
  dt_m2[, measure := meas]
  braz_med <- rbind(braz_med, dt_m2)
}
braz_med[, simulation_index := NULL]
braz_med <- dcast(braz_med, week + strain ~ measure, value.vars=value)

ggplot(braz_med) + 
  geom_line(aes(week, median/1000000)) +
  geom_ribbon(aes(week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill=strain), alpha=0.3) +
  geom_ribbon(aes(week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill=strain), alpha=0.5) +
  facet_grid(strain~.) + scale_x_continuous(breaks=2010:2020) +
  scale_color_manual(values = strain_colors1) +
  scale_fill_manual(values = strain_colors1) +
  ylab('Infections, millions') + xlab('Year') + 
  theme_minimal()

braz_cum15 <- braz_l[year(week) < 2016,]
braz_cum15 <- braz_cum15[,strain:=NULL][, lapply(.SD,sum), by=c('simulation_index','week')]
braz_cum15 <- braz_cum15[,cum.sum := cumsum(value), by=c('simulation_index')]

ggplot(braz_cum15) + 
  geom_line(aes(week, cum.sum/1000000, group=simulation_index)) +
  # geom_ribbon(aes(week, ymin=eti95L/1000000, ymax=eti95U/1000000, fill=strain), alpha=0.3) +
  # geom_ribbon(aes(week, ymin=eti50L/1000000, ymax=eti50U/1000000, fill=strain), alpha=0.5) +
  # scale_x_continuous(breaks=2010:2016) +
  scale_color_manual(values = strain_colors1) +
  scale_fill_manual(values = strain_colors1) +
  ylab('Infections, millions') + xlab('Year') + 
  theme_minimal()

braz_cum15med <- braz_cum15[, lapply(.SD, median, na.rm=T), by = c('week')]
braz_cum15med[, measure := 'median']
for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt_m2 <- braz_cum15[, lapply(.SD, get(meas)), by = c('week')]
  dt_m2[, measure := meas]
  braz_cum15med <- rbind(braz_cum15med, dt_m2)
}
braz_cum15med[, c('simulation_index','value') := NULL]
braz_cum15med_w <- dcast(braz_cum15med, week ~ measure, value.var='cum.sum')

ggplot(braz_cum15med_w) + 
  geom_line(aes(week, median/1000000)) +
  geom_ribbon(aes(week, ymin=eti95L/1000000, ymax=eti95U/1000000), alpha=0.3) +
  geom_ribbon(aes(week, ymin=eti50L/1000000, ymax=eti50U/1000000), alpha=0.5) +
  # facet_grid(strain~.) + 
  # scale_x_continuous(breaks=2010:2015) +
  # scale_color_manual(values = strain_colors1) +
  # scale_fill_manual(values = strain_colors1) +
  ylab('Infections, millions') + xlab('Year') + 
  theme_minimal()

arg <- data.table(readRDS("data/original_epids_output/exemplar_ARG_2010_2019_100.rds")[[1]])
arg[, c('country','country_code'):=NULL]
arg_l <- melt(arg, id.vars = c('simulation_index','week'))[,strain := substr(variable,4,4)][,variable:=NULL]
arg_l <- arg_l[,lapply(.SD,sum),by=c('simulation_index','week','strain')]
arg_med <- arg_l[, lapply(.SD, median, na.rm=T), by = c('week','strain')]
arg_med[, measure := 'median']
for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt_m2 <- arg_l[, lapply(.SD, get(meas)), by = c('week','strain')]
  dt_m2[, measure := meas]
  arg_med <- rbind(arg_med, dt_m2)
}
arg_med[, simulation_index := NULL]
arg_med <- dcast(arg_med, week + strain ~ measure, value.vars=value)
arg_cum15 <- arg_l[year(week) < 2016,]
arg_cum15 <- arg_cum15[,strain:=NULL][, lapply(.SD,sum), by=c('simulation_index','week')]
arg_cum15 <- arg_cum15[,cum.sum := cumsum(value), by=c('simulation_index')]
arg_cum15med <- arg_cum15[, lapply(.SD, median, na.rm=T), by = c('week')]
arg_cum15med[, measure := 'median']
for(meas in c('eti50L', 'eti50U', 'eti95L', 'eti95U')){
  dt_m2 <- arg_cum15[, lapply(.SD, get(meas)), by = c('week')]
  dt_m2[, measure := meas]
  arg_cum15med <- rbind(arg_cum15med, dt_m2)
}
arg_cum15med[, c('simulation_index','value') := NULL]
arg_cum15med_w <- dcast(arg_cum15med, week ~ measure, value.var='cum.sum')

cum15med <- rbind(cbind(braz_cum15med_w,'BRA'), cbind(arg_cum15med_w,'ARG'))

ggplot(cum15med) + 
  geom_line(aes(week, median/1000000, group=V2)) +
  geom_ribbon(aes(week, ymin=eti95L/1000000, ymax=eti95U/1000000, group=V2, fill=V2), alpha=0.3) +
  geom_ribbon(aes(week, ymin=eti50L/1000000, ymax=eti50U/1000000, group=V2, fill=V2), alpha=0.5) +
  # facet_grid(strain~.) +
  # scale_x_continuous(breaks=2010:2015) +
  # scale_color_manual(values = strain_colors1) +
  # scale_fill_manual(values = strain_colors1) +
  ylab('Infections, millions') + xlab('Year') + 
  theme_minimal()

cum15med[, pop := ifelse(V2=='BRA', 212.56*1000000, 45.20*1000000)]
ggplot(cum15med) + 
  geom_line(aes(week, median/pop, group=V2)) +
  geom_ribbon(aes(week, ymin=eti95L/pop, ymax=eti95U/pop, group=V2, fill=V2), alpha=0.3) +
  geom_ribbon(aes(week, ymin=eti50L/pop, ymax=eti50U/pop, group=V2, fill=V2), alpha=0.5) +
  # facet_grid(strain~.) +
  # scale_x_continuous(breaks=2010:2015) +
  # scale_color_manual(values = strain_colors1) +
  # scale_fill_manual(values = strain_colors1) +
  ylab('Infections, per person (ish)') + xlab('Year') + 
  theme_minimal()

## what does the cumulative flu data look like?
sa_flu_dt <- data.table(sa_flu)
sa_flu_cum <- sa_flu_dt[,cum.sum := cumsum(value), by=c('country','STRAIN')]

ggplot(sa_flu_cum) + 
  geom_line(aes(ISO_WEEKSTARTDATE, cum.sum/1000000, group=country, col=country),lwd=1) +
  facet_grid(STRAIN~.) +
  theme(text=element_text(size=14)) +
  ylab('Infections, millions') + xlab('Year') + 
  theme_bw()

sa_flu_cum_agg <- sa_flu_dt[,c('country','ISO_WEEKSTARTDATE','value')][,lapply(.SD,sum), by=c('country','ISO_WEEKSTARTDATE')][,cum.sum := cumsum(value), by=c('country')]

ggplot(sa_flu_cum_agg[year(ISO_WEEKSTARTDATE)<2016]) + 
  geom_line(aes(ISO_WEEKSTARTDATE, cum.sum/1000000, group=country, col=country),lwd=1) +
  ylab('Infections, millions') + xlab('Year') + 
  theme(text=element_text(size=14)) +
  theme_bw()

## comparing the two
comb_data <- rbind(cbind(braz_med,'BRA'), cbind(arg_med,'ARG'))
setnames(comb_data, 'V2','COUNTRY_CODE')
setnames(sa_flu_cum, 'ISO_WEEKSTARTDATE','week')
comb_data[, STRAIN := paste0('INF_', strain)]
comb_data <- comb_data[sa_flu_cum, on=c('week','COUNTRY_CODE','STRAIN'), flunet := i.value]
comb_data_l <- melt(comb_data[,c('week','STRAIN','median','COUNTRY_CODE','flunet')], id.vars=c('week','STRAIN','COUNTRY_CODE'))

ggplot(comb_data_l[COUNTRY_CODE=='BRA']) + 
  geom_line(aes(week, value, col=variable), lwd=1) +
  facet_grid(variable~STRAIN, scales='free') + 
  scale_x_continuous(breaks=2010:2020) +
  ylab('Infections') + xlab('Year') + 
  theme_minimal() + theme(text=element_text(size=14))

ggplot(comb_data_l[variable=='flunet']) + 
  geom_line(aes(week, value, col=COUNTRY_CODE), lwd=1) +
  facet_grid(COUNTRY_CODE~STRAIN, scales='free') + 
  # scale_x_continuous(breaks=2010:2020) +
  ylab('Infections') + xlab('Year') + 
  theme_bw() + theme(text=element_text(size=14),
                     legend.position='none')

ggplot(comb_data_l[,value2 := value][variable=='flunet',value2:=5000*value]) + 
  geom_line(aes(week, value2, group=variable, col=variable), lwd=1) +
  facet_grid(COUNTRY_CODE~STRAIN, scales='free') + 
  scale_x_continuous(breaks=2010:2020) +
  theme(text=element_text(size=14)) +
  ylab('Infections') + xlab('Year') + 
  theme_minimal()

cov_data %>% filter(country%in%c('Argentina','Brazil')) %>% pivot_longer(!c(country, year)) %>%
  ggplot() + ylab('Vaccination coverage (%)') + xlab('Year') +
  geom_line(aes(x=year, y=value*100, group=country, col=country), lwd=1) +
  facet_wrap(name~.) + theme_bw() + ylim(c(0,100)) + 
  scale_x_continuous(breaks=2010:2019)

cov_data %>% filter(country %in% c('Argentina','Brazil')) %>% 
  pivot_longer(!c(country, year)) %>%
  group_by(country) %>% summarise(mean(value))


(100000*0.000802/(cases_as[country_code=='BRA']$cum.sum/(5*203000000)))
cdc_data_l[iso3c=='BRA']























