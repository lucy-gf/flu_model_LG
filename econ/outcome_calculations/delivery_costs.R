
## VACCINE DELIVERY COST EXTRACTION AND PREDICTION
setwd("~/Desktop/research asst/Global Code")

options(scipen=100000)
source('BS/BS_colors.R')
library(countrycode)
library(ggplot2)
library(data.table)
library(readr)
library(WDI)
library(patchwork)
library(readxl)

## data from:
## Producing Standardized Country-Level Immunization Delivery Unit Cost Estimates
## https://doi.org/10.1007/s40273-020-00930-6

# economic or financial cost per dose
costing_type <- c('economic','financial')[1]

deliv_raw <- data.table(read_csv(paste0('econ/outcome_calculations/data/extracted_data/delivery_costs_', costing_type,'.csv'),
                        show_col_types=F))
deliv_raw[, iso3c := suppressWarnings(countrycode(Country, origin='country.name', destination = 'iso3c'))][grepl('Micronesia', Country), iso3c := 'FSM']
deliv_raw[, gdppc2018 := as.numeric(gsub('[$,]', '', unname(unlist(deliv_raw[,4]))))]
deliv_raw[, cost_per_dose := as.numeric(gsub('[$ ]', '', substr(unname(unlist(deliv_raw[,8])),1,6)))]
deliv_raw[, cost_per_dose_l95 := as.numeric(gsub('[$ ()]', '', unname(unlist(strsplit(unname(unlist(deliv_raw[,8])), split = "[ –]")))[3*(1:nrow(deliv_raw)) - 1]))]
deliv_raw[, cost_per_dose_u95 := as.numeric(gsub('[$ ()]', '', unname(unlist(strsplit(unname(unlist(deliv_raw[,8])), split = "[ –]")))[3*(1:nrow(deliv_raw))]))]

setnames(deliv_raw, 'World Bank Income Levela [10]', 'income_level')

deliv_data <- deliv_raw[, c('iso3c','income_level','gdppc2018','cost_per_dose',
                            'cost_per_dose_l95', 'cost_per_dose_u95')]

ggplot(deliv_data) +
  geom_point(aes(x=gdppc2018, y=cost_per_dose, col=income_level)) +
  geom_segment(aes(x=gdppc2018, xend=gdppc2018,
                   y=cost_per_dose_l95, yend=cost_per_dose_u95, col=income_level),
               alpha=0.5) + theme_bw() + labs(col='Income group') +
  scale_x_log10('GDP per capita (USD 2018)') + scale_y_log10('Delivery cost (USD 2018)')
ggsave(paste0('econ/outcome_calculations/plots/portnoy_data.png'),
       width=20,height=15,units="cm")

deliv_data[, c('income_level','gdppc2018') := NULL]

## filtering out countries not in analysis
ex_data <- data.table(read_csv('econ/outcome_calculations/data/new_clustering.csv', show_col_types = F))
setnames(ex_data, 'codes','iso3c')
deliv_data <- deliv_data[iso3c %in% ex_data$iso3c]
n_portnoy_countries <- nrow(deliv_data) #129

## adding in countries without estimates
deliv_data <- rbind(deliv_data,
                    data.table(
                      iso3c = setdiff(ex_data$iso3c, deliv_data$iso3c),
                      cost_per_dose = NA, cost_per_dose_l95 = NA, cost_per_dose_u95 = NA
                      ))

## inflating delivery cost estimates to $2022
inflation_rates <- data.table(WDI(indicator='NY.GDP.DEFL.ZS', start=2010, end=2022))
USD_rates <- inflation_rates[iso3c=='USA',]
USD_inflation_2010 <- as.numeric(USD_rates[year==2022]$NY.GDP.DEFL.ZS/ USD_rates[year==2010]$NY.GDP.DEFL.ZS)
USD_inflation_2017_to_2021 <- as.numeric(USD_rates[year==2021]$NY.GDP.DEFL.ZS/ USD_rates[year==2017]$NY.GDP.DEFL.ZS)
USD_inflation_2020_to_2021 <- as.numeric(USD_rates[year==2021]$NY.GDP.DEFL.ZS/ USD_rates[year==2020]$NY.GDP.DEFL.ZS)
USD_inflation_2018 <- as.numeric(USD_rates[year==2022]$NY.GDP.DEFL.ZS/ USD_rates[year==2018]$NY.GDP.DEFL.ZS)
USD_inflation_2021 <- as.numeric(USD_rates[year==2022]$NY.GDP.DEFL.ZS/ USD_rates[year==2021]$NY.GDP.DEFL.ZS)
deliv_data[, cost_per_dose_2022 := cost_per_dose*USD_inflation_2018]
deliv_data[, cost_per_dose_l95 := cost_per_dose_l95*USD_inflation_2018]
deliv_data[, cost_per_dose_u95 := cost_per_dose_u95*USD_inflation_2018]
setnames(deliv_data, 'cost_per_dose', 'cost_per_dose_portnoy')

######################################################
#### EXTRAPOLATING VALUES TO HICS WITH REGRESSION ####
######################################################

# exchange rates to USD:
lcu_rates <- data.table(WDI(indicator='PA.NUS.FCRF', start=2010, end=2022))
gbp_to_usd_2018 <- 1/as.numeric(lcu_rates[iso3c=='GBR' & year==2018]$PA.NUS.FCRF)
eur_to_usd_2010 <- 1/as.numeric(lcu_rates[iso3c=='EMU' & year==2010]$PA.NUS.FCRF)
eur_to_usd_2017 <- 1/as.numeric(lcu_rates[iso3c=='EMU' & year==2017]$PA.NUS.FCRF)

## adding in some HIC data

deliv_data[iso3c=='USA', cost_per_dose_2022 := USD_inflation_2021*25]
## https://www.healthsystemtracker.org/chart-collection/where-do-americans-get-vaccines-and-how-much-does-it-cost-
## to-administer-them/#Where%20adults%20receive%20flu%20vaccine,%20by%20race%20and%20ethnicity,%202018

deliv_data[iso3c=='GBR', cost_per_dose_2022 := USD_inflation_2018*gbp_to_usd_2018*14.05]
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6930088/

deliv_data[iso3c=='ESP', cost_per_dose_2022 := USD_inflation_2010*eur_to_usd_2010*(2 + 4.88 + 7)/3]
## https://www.gacetasanitaria.org/es-analisis-coste-efectividad-vacunacion-antineumococica-espana-articulo-S0213911111000963
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4994732/

#############
#############

## healthcare expenditure per capita
health_data <- data.table(WDI(indicator='SH.XPD.CHEX.PC.CD', start=2010, end=2022))
health_data <- health_data[iso3c %in% deliv_data$iso3c]
setnames(health_data, 'SH.XPD.CHEX.PC.CD', 'health_expen')
health_data_2021 <- health_data[year==2021][,c('iso3c','health_expen')]
setdiff(ex_data$iso3c, health_data_2021$iso3c)

## FILLING IN MISSING DATA
health_data_2021 <- rbind(health_data_2021, data.table(iso3c=c('GUF','TWN'),
                                             health_expen=c(NA, USD_inflation_2017_to_2021*1595)))

missing <- health_data_2021[is.na(health_expen)]$iso3c

health_data_2021[iso3c=='CUB',health_expen := as.numeric(health_data[iso3c=='CUB' & year==2020]$health_expen)*USD_inflation_2020_to_2021]
health_data_2021[iso3c=='XKX',health_expen := eur_to_usd_2017*158*USD_inflation_2017_to_2021]
# https://msh.rks-gov.net/Documents/DownloadDocument?fileName=Raporti%20p%C3%ABr%20NHA%20%E2%80%93%20ENG52070994.0793.pdf

# matching further missing countries to countries with similar GDP per capita
load('econ/outcome_calculations/data/predicted_costs')
gdp_data <- pred_costs[simulation_index==1 & outcome=='outpatient' & study_pop=='adults'][, c('iso3c','gdpcap')]
missing <- health_data_2021[is.na(health_expen)]$iso3c
gdp_data <- gdp_data[order(gdpcap)]
missing_closest <- data.table(country = missing, closest = NULL)
for(i in 1:nrow(missing_closest)){
  row <- which(gdp_data$iso3c %in% c(missing_closest[i, 'country']))
  closest_in <- ifelse(abs(gdp_data[row - 1]$gdpcap - gdp_data[row]$gdpcap) <= abs(gdp_data[row + 1]$gdpcap - gdp_data[row]$gdpcap),
    gdp_data[row - 1]$iso3c, gdp_data[row + 1]$iso3c)
  missing_closest[i, closest := closest_in]
  health_data_2021[iso3c==missing_closest[i,'country'],
                   health_expen := health_data_2021[iso3c==missing_closest[i,'closest']]$health_expen]
}

health_data_2022 <- copy(health_data_2021)
# inflating to USD 2022
health_data_2022[, health_expen := USD_inflation_2021*health_expen]

####################
#### REGRESSION ####
####################

deliv_data <- deliv_data[health_data_2022, on='iso3c']

regression_plot <- ggplot(deliv_data) +
  geom_point(aes(x=health_expen, y=cost_per_dose_2022, color=is.na(cost_per_dose_portnoy))) +
  scale_x_log10('Healthcare expenditure per capita (USD 2022)') + scale_y_log10('Delivery cost (USD 2022)') +
  theme_bw() + theme(legend.position='none'); regression_plot
ggsave(paste0('econ/outcome_calculations/plots/healthcare_pre_regression.png'),
       width=20,height=15,units="cm")

## adding weights
# 29 studies informed the portnoy study
deliv_data[!is.na(cost_per_dose_portnoy), weight := 29/n_portnoy_countries]
deliv_data[is.na(cost_per_dose_portnoy), weight := 1]

## training the model
lmodel <- lm((cost_per_dose_2022) ~ ((health_expen)), data = deliv_data,
             weights = deliv_data$weight)

deliv_data[, pred_cost := (predict(lmodel, deliv_data))]
deliv_data[, pred_cost_conf_l95 := (predict(lmodel, deliv_data, interval='confidence')[,2])]
deliv_data[, pred_cost_conf_u95 := (predict(lmodel, deliv_data, interval='confidence')[,3])]
deliv_data[, pred_cost_pred_l95 := suppressWarnings(predict(lmodel, deliv_data, interval='prediction')[,2])]
deliv_data[, pred_cost_pred_u95 := suppressWarnings(predict(lmodel, deliv_data, interval='prediction')[,3])]

interval <- c('confidence','prediction')[2]

deliv <- ggplot(deliv_data) +
  geom_point(aes(x=health_expen, y=cost_per_dose_2022, alpha=weight)) +
  geom_line(aes(x=health_expen, y=pred_cost), col=1, lty=2) +
  geom_ribbon(aes(x=health_expen, ymin=get(paste0('pred_cost_',substr(interval,1,4), '_l95')),
                  ymax=get(paste0('pred_cost_',substr(interval,1,4), '_u95'))), fill=2, lty=2, alpha=0.3) +
  # scale_x_log10() + scale_y_log10() +
  theme_bw() + ylab('Delivery cost') +
  xlab('Healthcare expenditure per capita ($ 2022)') +
  theme(panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"),
        legend.position='none'); deliv

ggsave(paste0('econ/outcome_calculations/plots/deliv_prediction_',interval,'.png'),
       width=20,height=15,units="cm")


###########################################
#### SAMPLING VALUES FOR ALL COUNTRIES ####
###########################################

delivery_samples <- data.table(iso3c = rep(deliv_data$iso3c, each=100),
                               simulation_index = 1:100,
                               delivery_cost = 0)

z <- qnorm(0.975, 0, 1)

for(c_code in unique(delivery_samples$iso3c)){
  if(!is.na(deliv_data[iso3c==c_code,'cost_per_dose_portnoy'])){
    mean_i <- log(unlist(unname(deliv_data[iso3c==c_code,'cost_per_dose_2022'])))
    sd_i <- (log(unlist(unname(deliv_data[iso3c==c_code,'cost_per_dose_u95']))) -
             log(unlist(unname(deliv_data[iso3c==c_code,'cost_per_dose_l95']))))/(2*z)
    samples_vec <- rlnorm(100, meanlog = mean_i, sdlog = sd_i)
    delivery_samples[iso3c==c_code, 'delivery_cost'] <- samples_vec
  }else{
    mean_i <- (unlist(unname(deliv_data[iso3c==c_code,'pred_cost'])))
    sd_i <- ((unlist(unname(deliv_data[iso3c==c_code,get(paste0('pred_cost_',substr(interval,1,4), '_u95'))]))) -
               (unlist(unname(deliv_data[iso3c==c_code,get(paste0('pred_cost_',substr(interval,1,4), '_l95'))]))))/(2*z)
    samples_vec <- rnorm(100, mean = mean_i, sd = sd_i)
    delivery_samples[iso3c==c_code, 'delivery_cost'] <- samples_vec
  }
  if(sum(samples_vec <= 0) > 0){print('Less than 0 cost')}
}

write_csv(delivery_samples, file='econ/outcome_calculations/data/delivery_cost_samples.csv')

delivery_samples <- delivery_samples[deliv_data, on='iso3c']

ggplot(delivery_samples) +
  geom_jitter(aes(x=health_expen, y=delivery_cost, group=iso3c), alpha = 0.2) + theme_bw() +
  geom_point(aes(x=health_expen, y=cost_per_dose_2022), col='red', shape=4) +
  geom_segment(data = delivery_samples[simulation_index==1],
               aes(x=health_expen, xend=health_expen,
                   y=cost_per_dose_l95, yend=cost_per_dose_u95), col='red', alpha=0.4) +
  scale_x_log10('Healthcare expenditure per capita (USD 2022)') +
  ylab('Delivery cost (USD 2022)')
ggsave(paste0('econ/outcome_calculations/plots/delivery_cost_samples.png'),
       width=20,height=15,units="cm")


portnoy_plot <- ggplot(delivery_samples[!is.na(cost_per_dose_portnoy)]) +
  theme_bw() +
  geom_point(aes(x=health_expen, y=cost_per_dose_2022), col='black') +
  geom_segment(data = delivery_samples[!is.na(cost_per_dose_portnoy) & simulation_index==1],
               aes(x=health_expen, xend=health_expen,
                   y=cost_per_dose_l95, yend=cost_per_dose_u95), col='black', alpha=0.4) +
  scale_x_log10('Healthcare expenditure per capita (USD 2022)', limits=c(15,1520), breaks=c(30,100,300,1000)) +
  scale_y_log10() + 
  ylab('Delivery cost (USD 2022)');portnoy_plot
ggsave(paste0('econ/outcome_calculations/plots/portnoy_data_health_exp.png'),
       width=20,height=15,units="cm")

ggplot() +
  theme_bw() +
  geom_point(data = delivery_samples[!is.na(cost_per_dose_portnoy)], aes(x=health_expen, y=cost_per_dose_2022), col='black') +
  geom_segment(data = delivery_samples[!is.na(cost_per_dose_portnoy) & simulation_index==1],
               aes(x=health_expen, xend=health_expen,
                   y=cost_per_dose_l95, yend=cost_per_dose_u95), col='black', alpha=0.4) +
  geom_point(data = delivery_samples[iso3c %in% c('USA','GBR','ESP')], aes(x=health_expen, y=cost_per_dose_2022), col='red') +
  scale_x_log10('Healthcare expenditure per capita (USD 2022)', breaks=c(10,30,100,300,1000,3000,10000)) +
  scale_y_log10(breaks=c(0.1,0.3,1,3,10,30)) + 
  ylab('Delivery cost (USD 2022)')
ggsave(paste0('econ/outcome_calculations/plots/all_delivery_data.png'),
       width=20,height=15,units="cm")



