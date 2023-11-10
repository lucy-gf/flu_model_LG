
## EFFICACY DATA  

# match/mismatch efficacy taken from UK paper
# age-specific, lower for >65s

match_efficacy <- c(rep(0.7,3),0.46) 
mismatch_efficacy <- c(rep(0.42,3),0.28)

# need annual Inf_A and Inf_B matches for each hemisphere

matches <- data.frame(hemisphere = rep(c('NH','SH'), each=20),
                      year = rep(2010:2019, 2),
                      strain_match = rep(c('INF_A', 'INF_B'), each=10),
                      match = rep(NA, 40))

# T for match, F for mismatch
matches[matches$hemisphere == 'NH' & matches$strain_match == 'INF_A',]$match <-
  c(T, F, F, T, F, T, F, F, T, T)
matches[matches$hemisphere == 'NH' & matches$strain_match == 'INF_B',]$match <-
  c(T, T, T, T, F, T, T, F, T, T)
matches[matches$hemisphere == 'SH' & matches$strain_match == 'INF_A',]$match <-
  c(T,T,T,T,T,F,T,F,T,F)
matches[matches$hemisphere == 'SH' & matches$strain_match == 'INF_B',]$match <-
  c(T,T,T,T,T,T,F,T,T,T)

matches <- matches %>% mutate('0-5' = match*match_efficacy[1] + (1-match)*mismatch_efficacy[1], 
                              '5-20' = match*match_efficacy[2] + (1-match)*mismatch_efficacy[2], 
                              '20-65' = match*match_efficacy[3] + (1-match)*mismatch_efficacy[3], 
                              '65+' = match*match_efficacy[4] + (1-match)*mismatch_efficacy[4]) 

#write.csv(matches, "data_for_cluster/matches.csv", row.names=F)

## COVERAGE DATA-FRAME FOR EACH COUNTRY + YEAR + AGE GROUP

cov_data <- data.frame(country = rep(df_cntr_table$country, each=10),
                       year = rep(2010:2019, length(df_cntr_table$country)))

cov_data <- cov_data %>% mutate('0-5' = NA, '5-20' = NA, '20-65' = NA, '65+' = NA) 
cov_data[cov_data$country %in% c('China', 'Ghana', 'Turkey'), 3:6] <- 0 # should Turkey be =0?
cov_data[cov_data$country %in% c('Argentina'), 3:6] <- rep(c(0.72,0.5,0.5,0.55), each = length(2010:2019))
cov_data[cov_data$country %in% c('Australia'), 3:6] <- rep(c(0.32,0.3,0.23,0.72), each = length(2010:2019))
cov_data[cov_data$country %in% c('Canada'), 3:6] <- rep(c(0.32,0.23,0.24,0.64), each = length(2010:2019))
cov_data[cov_data$country %in% c('United Kingdom') & cov_data$year %in% 2010:2012, 3:6] <- rep(c(0,0,0.08,0.73), each = length(2010:2012))
cov_data[cov_data$country %in% c('United Kingdom') & cov_data$year %in% 2013:2019, 3:6] <- rep(c(0.35,0.48,0.15,0.73), each = length(2013:2019))

level_order <- c('0-5', '5-20', '20-65', '65+') 

# write.csv(cov_data, "data_for_cluster/cov_data.csv", row.names=F)

# cov_data %>% pivot_longer(!country:year) %>% rename(Coverage = value,
#                                                      Age = name) %>%
#   ggplot(aes(x=factor(Age, levels = level_order), y=Coverage, group=interaction(country,year),
#              col=country, fill=country)) +
#   theme_minimal() + scale_y_continuous(limits=c(0,1)) +
#   #geom_bar(position="dodge", stat="identity") +
#   geom_line(lwd=1) + geom_point(size=2) + xlab('Age')









