
# source functions, call libraries
source("fcns/fcns.R")
# load data
source("fcns/load_flunet.R")

# load filtered data
# this csv file is saved version of `df_posit_sel_cntrs` dataframe in the flunet_data.R file
n_sel=1
df_posit_sel_cntrs <- read_csv(c("output/df_positivity_counts_sel_cntrs.csv",
                                 "output/df_positivity_counts_sel_cntrs_ALTERN.csv")[n_sel])
df_posit_sel_cntrs <- df_posit_sel_cntrs %>% 
  mutate(ISO_WEEKSTARTDATE = as.Date(ISO_WEEKSTARTDATE, format = "%d/%m/%Y"))

countries <- c("Argentina", "Australia", "Canada", "China", "Ghana", "Turkey", 'United Kingdom', 'UK')

# PLOT #1
# IDENTIFY THE EPIDEMICS
# set thresholds (length of season, lower and upper percentiles to include as an epidemic)
length_lim_val=8
low_thresh_val=data.frame(country = countries[1:7], val = rep(0.5,7))
up_thresh_val=data.frame(country = countries[1:7], val = rep(0.65,7))
end_thresh_val=low_thresh_val
up_thresh_val$val[up_thresh_val$country=='Argentina'] <- 0.80
up_thresh_val$val[up_thresh_val$country=='Ghana'] <- 0.80
up_thresh_val$val[up_thresh_val$country=='Turkey'] <- 0.80
up_thresh_val$val[up_thresh_val$country=='United Kingdom'] <- 0.80
end_thresh_val$val[end_thresh_val$country%in% c('Australia')] <- 0.64
end_thresh_val$val[end_thresh_val$country%in% c('China')] <- 0.6

sel_variable_val=c("value","positivity")[1] # value means counts
df_epid_threshold <- fcn_identify_seasons_three_thresh(
  df_input = df_posit_sel_cntrs %>% 
                filter(ISO_WEEKSTARTDATE>=as.Date("2010-01-01") & ISO_WEEKSTARTDATE<=as.Date("2020-04-01")) %>%
                group_by(country,ISO_YEAR,ISO_WEEK,ISO_WEEKSTARTDATE,cluster_name,STRAIN,metasource) %>%
                summarise(SPEC_PROCESSED_NB=sum(SPEC_PROCESSED_NB),value=sum(value)) %>% ungroup() %>% 
                filter(metasource=='NONSENT'),
  sel_variable=sel_variable_val, source_varname="metasource",
  low_thresh_df=low_thresh_val,up_thresh_df = up_thresh_val,length_lim=length_lim_val,
  end_thresh_df=end_thresh_val, rolling_yrs=NA, cushioning=0)
# save
#write_csv(df_epid_threshold,file = "output/df_epid_threshold.csv")
# write_csv(df_epid_threshold,file = "data_for_cluster/df_epid_threshold.csv")
# df_epid_threshold <- read.csv("output/df_epid_threshold.csv")

df_epid_threshold <- df_epid_threshold %>% 
  mutate(ISO_WEEKSTARTDATE = as.Date(ISO_WEEKSTARTDATE, format = "%Y-%m-%d"))

# PLOTTING epidemics
y_log_flag <- F
df_epid_lims <- fcn_find_bloc_lims(df_epid_threshold,log_flag=y_log_flag)
#write_csv(df_epid_lims,file = "data_for_cluster/df_epid_lims.csv")

df_epid_threshold %>% filter(metasource=='NONSENT') %>% 
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  #filter(country %in% countries) %>% 
  ggplot() + facet_grid(country~STRAIN,scales="free_y") +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=metasource)) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_included,group=metasource), col=1, lty=2) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_end,group=metasource), col=3, lty=4) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_peak,group=metasource), col=2, lty=3) +
  xlab("") + ylab("Positive influenza tests") + labs(colour="epidemic",fill="") +
  scale_color_manual(values=c("grey32","red","grey46","darkred")) +
  scale_fill_manual(values=c("red","darkred")) +
  scale_x_date(date_labels="%Y\n%m",breaks="6 month",expand=expansion(0.01,0))  +
  theme_bw() + standard_theme + theme(legend.position="none",strip.text=element_text(size=13)) +
  guides(color=guide_legend(nrow=2))
ggsave("SM_plots/epid_values_unidentif.png",
       width=36,height=24,units="cm")

strain_labs <- c('Influenza A','Influenza B')
names(strain_labs) <- c('INF_A','INF_B')
country_colours <- c("TUR"="#FDC328FF", "ARG" = '#65156EFF', "AUS" = '#A92395FF', 
                     "GBR" = '#E56B5DFF', "CHN" = '#F89441FF', "CAN" = '#CC4678FF', 
                     "GHA" = '#FCFFA4FF') 
strain_colors <- c('INF_A' = '#7d66ac', 'INF_B' = '#e483a4')

values <- df_epid_threshold %>% filter(metasource=='NONSENT') %>% 
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  #filter(country %in% countries) %>% 
  ggplot() + facet_grid(country~STRAIN,scales="free_y", 
                        labeller = labeller(STRAIN = strain_labs)) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_included,group=metasource), col=1, lty=2) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_end,group=metasource), col=3, lty=4) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_peak,group=metasource), col=2, lty=3) +
  xlab("") + ylab("Positive influenza tests") + labs(colour="epidemic",fill="") +
  scale_x_date(date_labels="%Y",breaks="2 year",expand=expansion(0.01,0))  +
  geom_rect(data=df_epid_lims[df_epid_lims$metasource=='NONSENT',], aes(xmin=start_date,xmax=end_date,
                                   fill=STRAIN,ymin=min_val,ymax=max_val),alpha=1/3) +
  scale_fill_manual(values=strain_colors, labels = c('Influenza A',
                                                     'Influenza B')) +
  # scale_color_manual(values=strain_colors) +
  # scale_color_manual(values=c("grey32","red","grey46","darkred")) +
  # scale_fill_manual(values=c("red","darkred")) +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=metasource)) +
  # geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=metasource)) +
  theme_bw() + 
  labs(fill='Strain') + theme(text=element_text(size=14)) +
  guides(color=guide_legend(nrow=2)); values
ggsave("SM_plots/epid_values.png",
       width=36,height=24,units="cm")

for(i in c(1,2,3,4,5,6,8)){
df_epid_threshold %>% filter(metasource=='NONSENT') %>% 
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  filter(country %in% countries[i]) %>% 
  ggplot() + facet_grid(STRAIN~.,scales="free_y") +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=metasource,colour=factor(source_epid))) +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_included,group=metasource), col=1, lty=2) +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_end,group=metasource), col=3, lty=4) +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_peak,group=metasource), col=2, lty=3) +
  xlab("") + ylab("Positive influenza tests") + labs(colour="epidemic",fill="") +
  scale_color_manual(values=c("grey32","red","grey46","darkred")) +
  scale_fill_manual(values=c("red","darkred")) +
  scale_x_date(date_labels="%Y\n%m",breaks="6 month",expand=expansion(0.01,0)) +
  geom_rect(data=df_epid_lims[df_epid_lims$metasource=='NONSENT'&
                                df_epid_lims$country==countries[i],], aes(xmin=as.Date(start_date),xmax=as.Date(end_date),
                                                                        fill=metasource,ymin=min_val,ymax=max_val),alpha=1/4) +
  theme_bw() + standard_theme + theme(legend.position="none",strip.text=element_text(size=13)) +
  guides(color=guide_legend(nrow=2)) + ggtitle(countries[i])
  ggsave(paste0("epid_identif/", countries[i], ".png"), 
         width=36,height=20,units="cm")

df_epid_threshold %>% filter(metasource=='NONSENT') %>% 
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  filter(country %in% countries[i], STRAIN %in% c('INF_A', 'INF_B')[1:2],
         epid_index > 0) %>% 
  ggplot() + facet_wrap(STRAIN~epid_index,scales="free") +
  #geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_included,group=metasource), col=1, lty=2) +
  #geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_peak,group=metasource), col=2, lty=3) +
  #geom_line(aes(x=ISO_WEEKSTARTDATE,y=flu_end,group=metasource), col=3, lty=4) +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=metasource,colour=factor(source_epid))) +
  xlab("") + ylab("Positive influenza tests") + labs(colour="epidemic",fill="") +
  scale_color_manual(values=c("grey32","red","grey46","darkred")) +
  scale_fill_manual(values=c("red","darkred")) + ylim(c(0,NA)) +
  theme_bw() + standard_theme + theme(legend.position="none",strip.text=element_text(size=13)) +
  guides(color=guide_legend(nrow=2)) + ggtitle(countries[i])
ggsave(paste0("epid_identif/", countries[i], "_sep.png"), 
       width=36,height=24,units="cm")
}

# background rates, avoiding the beginning and end of the inference period
# to exclude any 'half-epidemics'
bg_rates <- df_epid_threshold %>% filter(epidem_inclusion==0,
                                              ISO_YEAR %in% 2011:2019) %>% 
  group_by(country, STRAIN) %>% summarise(mean(value))

positivity <- df_epid_threshold %>% filter(metasource=='NONSENT') %>% 
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  ggplot() + facet_grid(country~STRAIN,scales="free_y") +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value/SPEC_PROCESSED_NB,group=metasource,colour=factor(source_epid))) +
  xlab("") + ylab("Test positivity") + labs(colour="epidemic",fill="") +
  scale_color_manual(values=c("grey32","red","grey46","darkred")) +
  scale_fill_manual(values=c("red","darkred")) +
  scale_x_date(date_labels="%Y\n%m",breaks="6 month",expand=expansion(0.01,0))  +
  geom_rect(data=df_epid_lims[df_epid_lims$metasource=='NONSENT',], aes(xmin=start_date,xmax=end_date,
                                                                        fill=metasource,ymin=min_val,ymax=1),alpha=1/4) +
  theme_bw() + standard_theme + theme(legend.position="none",strip.text=element_text(size=13)) +
  guides(color=guide_legend(nrow=2)) 
positivity
ggsave("SM_plots/epid_positivity.png", 
       width=36,height=24,units="cm")

positivity + values + plot_layout(nrow=2) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')  ')
ggsave("SM_plots/epid_both.png", 
       width=40,height=50,units="cm")

# both strains on one graph
df_epid_threshold %>% filter(metasource=='NONSENT') %>% 
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  ggplot() + facet_grid(country~.,scales="free_y") + #facet_grid(country~STRAIN,scales="free_y") +
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=STRAIN,colour=STRAIN)) +
  xlab("") + ylab("# positive tests") + labs(colour="epidemic",fill="") +
  scale_color_manual(values=c("grey32","red","grey46","darkred")) +
  scale_fill_manual(values=c("red","darkred")) +
  # scale_x_date(date_labels="%Y\n%m",breaks="6 month",expand=expansion(0.01,0))  +
  geom_rect(data=df_epid_lims[df_epid_lims$metasource=='NONSENT',], aes(xmin=start_date,xmax=end_date,
                                  fill=metasource,ymin=min_val,ymax=max_val),alpha=1/4) +
  theme_bw() + standard_theme + theme(legend.position="top",strip.text=element_text(size=13)) +
  guides(color=guide_legend(nrow=2))

# SAVE
# ggsave(paste0("output/plots/epid_identif/",ifelse(n_sel>1,"ALT_CNTRS/",""),
#               "low_thresh",low_thresh_val, "_up_thresh",up_thresh_val,
#               "_length_lim",length_lim_val,ifelse(grepl("pos",sel_variable_val),"_byposit","_bycount"),
#               ifelse(y_log_flag,"_ylog",""),".png"),
#        width=36,height=24,units="cm")


### ### ### ### ### ### ### ### ### ### ### ### 
# list of epidemics
#
# there are no epidemics identified based on sentinel data only, 
# so lets use the season limits from the nonsentinel (this includes nonsentinel and notdefined in FluNet)
# because this has higher counts, more years etc

# lets first subset data for 1 country, 1 strain, 1 datatype ---> fitting 
data_fitting <- df_epid_threshold %>% 
  filter(country %in% "Canada" & STRAIN %in% "INF_A" & metasource %in% "NONSENT") %>%
  select(!c(flu_peak,over_peak,flu_included,over_inclusion,seq))

if (all(data_fitting$epidem_inclusion==(data_fitting$epid_index>0))) {
  data_fitting <- data_fitting %>% select(!epidem_inclusion)   }

# df_epid_lims %>% filter(country %in% "Canada" & STRAIN %in% "INF_A" & metasource %in% "NONSENT")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Get contact matrices and aggregate to our age groups:
# [0-5), [5-20), [20-65), [65-]
library(socialmixr)
# FluEvidenceSynthesis
library(fluEvidenceSynthesis)

# load CONTACT MATRICES from [Prem 2021]
# https://github.com/kieshaprem/synthetic-contact-matrices/tree/master/output/syntheticmatrices
load("data/contact_all.rdata") 
# matrices are `contact_all$GHA` ...
df_cntr_table = data.frame(country=unique(df_posit_sel_cntrs$country),
                          COUNTRY_CODE=unique(df_posit_sel_cntrs$COUNTRY_CODE)) %>%
                mutate(country_sub=case_when(
                  COUNTRY_CODE %in% "UK" ~ "GBR",
                  COUNTRY_CODE %in% "AUS" ~ "NZL", TRUE ~ COUNTRY_CODE))
# extract matrices from `contact_all`
cm_list <- lapply(df_cntr_table$country_sub, function(x) contact_all[[x]])
names(cm_list) <- df_cntr_table$country

# adding number of epidemics to fit to each country:
x <- 4
for(strain in c('A','B')){
  for(source in c('NONSENT','SENTINEL')){
    df_cntr_table[,x] <- rep.int(NA,nrow(df_cntr_table))
    colnames(df_cntr_table)[x] <- paste0('INF_',strain,'_',source) 
    x <- x + 1
  }
}
for(i in 1:nrow(df_cntr_table)){
  for(strain_loop in c('A','B')){
    for(source in c('NONSENT','SENTINEL')){
      suppressWarnings(df_cntr_table[i,paste0('INF_',strain_loop,'_',source)] <- 
                         max(df_epid_lims$index[df_epid_lims$country %in% c(df_cntr_table$country[i],df_cntr_table$COUNTRY_CODE[i]) &
                                                  df_epid_lims$STRAIN == paste0('INF_',strain_loop) &
                                                  df_epid_lims$metasource == source]))
      if(df_cntr_table[i,paste0('INF_',strain_loop,'_',source)]==-Inf){
        df_cntr_table[i,paste0('INF_',strain_loop,'_',source)] <- 0
      }
    }
  }
}

# adding hemispheres
df_cntr_table <- df_cntr_table %>% 
  mutate(hemisphere = c('SH', 'SH', 'NH', 'NH', 'NH', 'NH', 'NH'))
write.csv(df_cntr_table, "data_for_cluster/df_cntr_table.csv", row.names=F)

# data.frame of number of epidemics starting in each year, in each country
# year_epis <- data.frame(country = rep(df_cntr_table$country, each = 2*length(2010:2019)),
#                         code = rep(df_cntr_table$COUNTRY_CODE, each = 2*length(2010:2019)),
#                         year = rep(2010:2019, 2*7),
#                         strain = rep(c('INF_A', 'INF_B'), each = length(2010:2019)),
#                         epidemic = NA,
#                         hemisphere = NA)
# for(i in 1:nrow(year_epis)){
#   year_epis$hemisphere[i] <- unname(unlist(df_cntr_table %>% filter(country == year_epis$country[i]) %>% 
#                                              select(hemisphere)))
#   year_epis$epidemic[i] <- sum(as.numeric(unlist(df_epid_lims %>% filter(country %in% c(year_epis$country[i], year_epis$code[i]),
#                           STRAIN == year_epis$strain[i],
#                           year(start_date) == year_epis$year[i]) %>% ungroup %>% 
#                             select(index))) > 0, na.rm=T)
# }
# year_epis <- year_epis %>% mutate(binary = 1 - (epidemic==0))
# ggplot(year_epis, aes(x = year, y = country, fill = factor(binary))) + 
#   geom_tile() + theme_minimal() + scale_x_continuous(breaks=2010:2019) +
#   facet_grid(.~strain)



#write_csv(df_cntr_table, file='data_for_cluster/df_cntr_table.csv')

# aggregate into our age groups
age_limits <- c(0,5,20,65); age_group_names <- paste0(age_limits,"-", c(age_limits[2:length(age_limits)],99))
# contact_matrix(age.limits = age_limits)

# with socialmixr
orig_agegroups <- colnames(readRDS("data/UK_contact_matrix_sum.RDS"))
xx <- as.numeric(gsub("\\+","",unlist(lapply(orig_agegroups, function(x) strsplit(x, split ="-")[[1]][1]))))
cm_polymod_uk <- contact_matrix(polymod, countries="United Kingdom", age.limits=xx)$matrix
cm_polymod_uk_merged <- contact_matrix(polymod, countries="United Kingdom", age.limits = c(0, 5, 20,65))$matrix

# we need popul sizes by our age groups
for (sel_cntr in df_cntr_table$country) {
# age group sizes we need (pop_age from socialmixr)
model_age_groups <- data.frame(agegroup_name=age_group_names, duration=diff(c(age_limits,120)),
                               wpp_agegroup_low=c(1,2,5,14), wpp_agegroup_high=c(1,4,13,16),
                               popul=pop_age(wpp_age(sel_cntr, 2015), age.limit=age_limits)$population)
# age group population sizes corresponding to [Prem 2021] matrices
standard_age_groups <- fun_cntr_agestr(i_cntr = sel_cntr, i_year="2020",
                              age_low_vals = seq(0,75,5), age_high_vals = c(seq(4,74,5),120))

# modify contact matrix to correspond to our age groups
sel_cntr_code <- df_cntr_table$country_sub[df_cntr_table$country %in% sel_cntr]
C_m_merged_nonrecipr <- fun_create_red_C_m(C_m_full=contact_all[[sel_cntr_code]],
                                        model_agegroups=model_age_groups,
                                        orig_age_groups_duration=standard_age_groups$duration,
                                        orig_age_groups_sizes=standard_age_groups$values)
# fun_create_red_C_m(C_m_full=list(cm_polymod_uk,contact_all[[sel_cntr_code]])[[2]],
#                    model_agegroups=model_age_groups,
#                    orig_age_groups_duration=standard_age_groups$duration,
#                    orig_age_groups_sizes=standard_age_groups$values)

# make it reciprocal for the larger group
if (!exists("list_contact_matr")) { list_contact_matr <- list()}
list_contact_matr[[sel_cntr]] <- fun_recipr_contmatr(C_m_full = C_m_merged_nonrecipr,
                        age_group_sizes = model_age_groups$popul)
}

write_rds(list_contact_matr,"output/list_contact_matr.RDS")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# # Engl ILI dataset
# as.data.frame(ili$ili) %>% mutate(week=1:52) %>% pivot_longer(!week,names_to = "age") %>% 
#   ggplot() + geom_line(aes(x=week,y=value)) + facet_wrap(~age,scales="free_y") + 
#   theme_bw() + standard_theme
# 
# ### ###
# 
# cromer_samples_by_age <- qread("/home/lshmk17/Downloads/cromer_samples_by_age.qs")
# cromer_summ_stat <- cromer_samples_by_age %>% group_by(risk_group,subtype,outcome,age) %>% 
#   summarise(median=median(value),mean=mean(value),
#             ci95_l=quantile(value,probs=0.025),ci95_u=quantile(value,probs=0.975),
#             ci50_l=quantile(value,probs=0.25),ci50_u=quantile(value,probs=0.75))
# 
# cromer_summ_stat %>% filter(age<99) %>% # mutate(across(c(mean,ci95_l,ci95_u,ci50_l,ci50_u), ~ ifelse(. == 0, NA, .))) %>%
#   ggplot(aes(x=age)) + facet_grid(outcome~risk_group,scales="free_y") +
#   geom_line(aes(y=mean*100,color=subtype),size=1.2) + 
#   geom_ribbon(aes(ymin=ci50_l*100,ymax=ci50_u*100,fill=subtype),alpha=1/5) +
#   ylab("probability (%)") + # scale_y_log10() + 
#   theme_bw() + standard_theme # 
