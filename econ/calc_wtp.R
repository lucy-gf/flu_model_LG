
## CALCULATING CURRET WTP THRESHOLD PER DALY
#setwd("~/Desktop/research asst/Global Code")

library(data.table)
library(readr)
library(ggplot2)

raw_data <- read_csv("econ/outcome_calculations/data/RAW_WTP_thresholds.csv", show_col_types=F)

names <- c()
cet <- c()
gdp_prop <- c()

for(i in 1:4){
  # print(i)
  transf_char <- gsub(' ', ',', raw_data[5*(i-1) + 1,])
  data_i <- strsplit(transf_char,',')[[1]]
  data_i <- gsub('_', ' ', data_i)
  names <- c(names, data_i)
  # print(length(data_i))
  
  transf_char <- gsub(' ', '_', raw_data[5*(i-1) + 2,])
  transf_char <- gsub(',', '', transf_char)
  data_i <- strsplit(transf_char,'_')[[1]]
  data_i <- as.numeric(substr(data_i, 2, 1000))
  cet <- c(cet, data_i)
  # print(length(data_i))
  
  transf_char <- gsub(' ', '_', raw_data[5*(i-1) + 3,])
  data_i <- strsplit(transf_char,'_')[[1]]
  data_i <- data_i[1:(length(data_i)/4)*4]
  data_i <- as.numeric(data_i)
  gdp_prop <- c(gdp_prop, data_i)
  # print(length(data_i))
}

data <- data.table(
  country = names,
  cet = cet,
  gdp_prop = gdp_prop
)

## adding iso3c
itzs <- data.table(read_csv('econ/outcome_calculations/data/new_clustering.csv',show_col_types=F))
for(i in 1:nrow(data)){
  name <- data[i,]$country
  iso3c_i <- itzs[country==name |
                    country_altern==name |
                    country_altern_2==name, codes]
  if(length(iso3c_i)>0){
    data[i, iso3c := iso3c_i]  
  }
}

data[grepl('Cote', country), iso3c := 'CIV']
data[grepl('Rep Congo', country), iso3c := 'COG']
data[grepl('Dem Rep Congo', country), iso3c := 'COD']
data[grepl('Gambia', country), iso3c := 'GMB']
data[grepl('Bissau', country), iso3c := 'GNB']
data[grepl('Turkiye', country), iso3c := 'TUR']
data[grepl('Slovak', country), iso3c := 'SVK']
data[grepl('Lao', country), iso3c := 'LAO']
data[grepl('Bahamas', country), iso3c := 'BHS']
data <- data[!is.na(iso3c)]
data[, country:=NULL]

data <- rbind(data, data.table(cet = NA, gdp_prop = NA, iso3c = setdiff(itzs$codes, data$iso3c)))

## adding GDP 2022
library(WDI)
gdp_data <- data.table(WDI(indicator='NY.GDP.PCAP.KD', start=2022, end=2022))
gdp_data_filt <- gdp_data[iso3c %in% unique(data$iso3c)]
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
                         gdpcap = c(15600, 32679)))
data <- data[gdp_data_filt, on=c('iso3c'), gdpcap := gdpcap]

data[, new_cet := gdpcap*gdp_prop]

ggplot(data) + 
  geom_point(aes(gdpcap,cet)) +
  geom_line(aes(gdpcap,gdpcap),lty=2)

## filling in 17 missing values
ggplot(data) + 
  geom_point(aes(gdpcap,gdp_prop)) +
  scale_x_log10()

lmodel <- lm(gdp_prop ~ log(gdpcap), data=data)

plot_dt <- data.table(gdpcap = 200:110000, gdp_prop = NA)
plot_dt[, gdp_prop := predict(lmodel, plot_dt)]
ggplot(data) + 
  geom_point(aes(gdpcap,gdp_prop)) +
  geom_line(data = plot_dt, aes(gdpcap,gdp_prop), lty=2) +
  scale_x_log10(breaks = c(300,1000,3000,10000,30000,100000))

data_na <- data[is.na(gdp_prop)]
data_na[, gdp_prop := predict(lmodel, data_na)]
data <- data[!is.na(cet)]
data <- rbind(data, data_na)
data[, new_cet := gdpcap*gdp_prop]

ggplot(data) + 
  geom_point(aes(gdpcap,new_cet, col = is.na(cet))) +
  geom_line(aes(gdpcap,gdpcap),lty=2)

data_save <- data[,c('iso3c','gdpcap','gdp_prop','new_cet')]
setnames(data_save, 'new_cet','cet')
write_csv(data_save, file = "econ/outcome_calculations/data/WTP_thresholds.csv")

inc <- read_csv('econ/world-bank-income-groups.csv',show_col_types=F)
inc <- data.table(inc)
inc <- inc[Year==2022 & Code %in% data$iso3c]
setnames(inc, 'Code','iso3c')
setnames(inc, "World Bank's income classification", "income")
data <- data[inc, on=c('iso3c'), income:=income]

ggplot(data[!is.na(income)]) + 
  geom_point(aes(gdpcap,new_cet, col = income)) +
  geom_line(aes(gdpcap,gdpcap),lty=2) + theme_minimal() +
  scale_y_log10() + scale_x_log10(breaks = c(300,1000,3000,10000,30000,100000))

ggplot(data[!is.na(income)]) + 
  geom_point(aes(gdpcap,gdp_prop, col = income)) +
  theme_minimal() +
  scale_x_log10(breaks = c(300,1000,3000,10000,30000,100000))



