
## COMPARING WTP THRESHOLDS BETWEEN PICHON-RIVIERE AND 50% OF GDP PER CAPITA

library(data.table)
library(readr)
library(ggplot2)
options(scipen=1000000)
setwd("~/Desktop/research asst/Global Code")
source("BS/BS_colors.R")

wtp_thresh <- data.table(read_csv('econ/outcome_calculations/data/WTP_thresholds.csv', show_col_type=F))
setnames(wtp_thresh, 'iso3c','country_code')
wtp_thresh[, gdp_thresh := 0.5*gdpcap]
income <- data.table(read_csv('econ/world-bank-income-groups.csv',show_col_types=F))
setnames(income, 'Code','country_code')
setnames(income, "World Bank's income classification",'income_class')
income <- income[country_code %in% wtp_thresh$country_code & Year==2022,]
wtp_thresh <- wtp_thresh[income, on='country_code']
wtp_thresh$income_class <- factor(wtp_thresh$income_class,levels=c('Low-income countries','Lower-middle-income countries',
                                                               'Upper-middle-income countries','High-income countries'))

ggplot(wtp_thresh) + 
  geom_point(aes(x=gdp_thresh, y=cet, col=income_class)) +
  scale_y_log10() + scale_x_log10() +
  theme_bw() + geom_line(aes(x=cet,y=cet),lty=2) +
  scale_color_manual(values = income_colors) +
  xlab('50% of GDP per capita') + labs(col='World Bank income group') +
  ylab('Willingness-to-pay threshold from Pichon-Riviere et al.')
ggsave(paste0("econ/wtp_comparison.png"),
       width=22,height=16,units="cm")


