
### expanding the ITZ clustering produced by Chen et al.
library(countrycode)
library(dplyr)
library(wpp2022)
library(ggtext)
library(ggplot2)

# CHEN ET AL.:

flu_ITZ_clusters <- read_csv("data/flu_ITZ_clusters.csv") %>% mutate(country_code = countrycode(country, 
                                                                                         origin = 'country.name',
                                                                                         destination = 'iso3c'))
length(unique(flu_ITZ_clusters$country)) # 109 countries

# Prem et al. 2021:
load("data/contact_all.rdata")
names(contact_all)
countrycode(names(contact_all), origin = 'iso3c', destination = 'country.name')

# countrycode(setdiff(flu_ITZ_clusters$country_code, names_df$codes), 
#             destination = 'country.name',
#             origin = 'iso3c') 

# data-frame of Prem et al.
names_df <- data.frame(codes = names(contact_all), country = countrycode(names(contact_all), 
                                                          destination = 'country.name',
                                                          origin = 'iso3c'), in_cluster = F, in_prem = T) 
names_df <- rbind(names_df, data.frame(codes = setdiff(flu_ITZ_clusters$country_code, names_df$codes), 
                                       country = countrycode(setdiff(flu_ITZ_clusters$country_code, names_df$codes), 
                                                             destination = 'country.name',
                                                             origin = 'iso3c') , 
                                       in_cluster = F, in_prem = F)) %>% arrange(codes)
  
#names_df <- names_df %>% add_row(codes = 'JPN', country = 'Japan', in_cluster = F) 

# Adding T if in Chen et al. clustering
for(i in 1:nrow(names_df)){
  if(sum(grepl(names_df$codes[i], flu_ITZ_clusters$country_code)) > 0){
    names_df$in_cluster[i] <- T
  }
}

sum(names_df$in_cluster) # 109
sum(names_df$in_prem) # 177

## alternative names
names_df <- names_df %>% mutate(country_altern = case_when(
  grepl('United States', country) ~ "United States of America",grepl('Palestinian', country) ~ "State of Palestine",
  grepl('Laos', country) ~ "Lao People's Democratic Republic",grepl('Iran', country) ~ "Iran (Islamic Republic of)",
  grepl('Moldova', country) ~ "Republic of Moldova",grepl('Macao', country) ~ "China, Macao SAR",
  grepl('Myanmar', country) ~ "Myanmar",grepl('Bolivia', country) ~ "Bolivia (Plurinational State of)",
  grepl('Venezuela', country) ~ "Venezuela (Bolivarian Republic of)",grepl('Brazzaville', country) ~ "Congo",
  grepl('Kinshasa', country) ~ "Democratic Republic of the Congo",grepl('Syria', country) ~ "Syrian Arab Republic",
  grepl('Vietnam', country) ~ "Viet Nam",grepl('Brunei', country) ~ "Brunei Darussalam",
  grepl('Turkey', country) ~ "Türkiye",grepl('South Korea', country) ~ "Republic of Korea",
  grepl('North Korea', country) ~ "Dem. People's Republic of Korea",grepl('Bosnia & Herzegovina', country) ~ "Bosnia and Herzegovina",
  grepl('Ivoir', country) ~ "Côte d'Ivoire",grepl('Cape Verde', country) ~ "Cabo Verde",
  grepl('Hong Kong SAR China', country) ~ "China, Hong Kong SAR",grepl('St. Lucia', country) ~ "Saint Lucia",
  grepl('Russia', country) ~ "Russian Federation",grepl('São Tomé & Príncipe', country) ~ "Sao Tome and Principe",
  grepl('Trinidad & Tobago', country) ~ "Trinidad and Tobago",grepl('Tanzania', country) ~ "United Republic of Tanzania",
  grepl('St. Vincent & Grenadines', country) ~ "Saint Vincent and the Grenadines",
  grepl('Eswatini', country) ~ "Swaziland",
  grepl('Timor-Leste', country) ~ "East Timor",
  T ~ country
)) %>% select(codes, country, country_altern, in_cluster)

# which countries aren't in  Chen et al. or Prem et al.?
pop_proj_WPP_data <- read_xlsx("data/WPP2022_POP_F02_1_POPULATION_5-YEAR_AGE_GROUPS_BOTH_SEXES.xlsx",
                               sheet = 2, skip = 16) %>% 
  select(!c(Index, Variant, Notes, `Location code`, `ISO3 Alpha-code`, `ISO2 Alpha-code`,
            `SDMX code**`, `Parent code`)) %>% 
  filter(Type == 'Country/Area') %>% rename(name = `Region, subregion, country or area *`)

WPP_tot_pop_df <- pop_proj_WPP_data %>% filter(Year == 2025) %>% 
  select(!c(Type, Year)) 
WPP_tot_pop_df[2:22] <- sapply(WPP_tot_pop_df[2:22],as.numeric)
WPP_tot_pop_df$tot_pop <- rowSums(WPP_tot_pop_df[,2:22])
WPP_tot_pop_df <- WPP_tot_pop_df %>% select(name, tot_pop)

not_in_CP <- WPP_tot_pop_df %>% filter(! name %in% unique(names_df$country) &
                            ! name %in% unique(names_df$country_altern),
                          tot_pop > 1000) %>% 
  arrange(desc(tot_pop))

names_df <- rbind(names_df, 
                  data.frame(codes = countrycode(not_in_CP$name, origin='country.name', destination='iso3c'),
                             country = not_in_CP$name,
                             country_altern = not_in_CP$name,
                             in_cluster = F)) %>% 
                    mutate(cluster = NA, capital = NA, latitude = NA, longitude = NA)
names_df[186, 1] <- 'XKX'; names_df[186, 3] <- 'Kosovo'

for(i in 1:nrow(names_df)){
  if(names_df$in_cluster[i] == T){
    names_df$cluster[i] <- unlist(flu_ITZ_clusters %>% 
      filter(country_code == names_df$codes[i],
             method == 'kmeans_cluster') %>% 
      select(cluster_name))
  }
}

# capital city package:
library(maps)
data(world.cities)
capitals <- world.cities %>% filter(capital == 1)

for(i in 1:nrow(names_df)){
  if(nrow(capitals %>% 
          filter(country.etc %in% c(names_df$country[i], names_df$country_altern[i])) %>% 
          select(name)) > 1){
    names_df$capital[i] <- unlist(capitals %>% 
                                    filter(country.etc %in% c(names_df$country[i], names_df$country_altern[i])) %>% 
                                    arrange(desc(pop)) %>% 
                                    select(name))[1]
  }
  if(nrow(capitals %>% 
          filter(country.etc %in% c(names_df$country[i], names_df$country_altern[i])) %>% 
          select(name)) == 1){
    names_df$capital[i] <- unlist(capitals %>% 
                                    filter(country.etc %in% c(names_df$country[i], names_df$country_altern[i])) %>% 
                                    select(name))
  }
  if(nrow(capitals %>%
          filter(country.etc %in% c(names_df$country[i], names_df$country_altern[i])) %>%
          select(name)) == 0){
    vec <- rep(0, nrow(capitals))
    for(j in 1:nrow(capitals)){
      if(grepl(capitals$country.etc[j], names_df$country[i])){
        vec[j] <- j
      }
    }
    if(sum(vec)>0){
      capital_name <- capitals$name[vec[vec > 0]]
      names_df$capital[i] <- capital_name
    }
  }
}

names_df %>% filter(is.na(capital))
names_df <- names_df %>% mutate(capital = case_when(
  grepl('Ivoir', country) ~ capitals$name[capitals$country.etc == 'Ivory Coast'],
  grepl('Czechia', country) ~ capitals$name[capitals$country.etc == 'Czech Republic'],
  grepl("North Korea", country) ~ capitals$name[capitals$country.etc == 'Korea North'],
  grepl("South Korea", country) ~ capitals$name[capitals$country.etc == 'Korea South'],
  grepl("United States", country) ~ capitals$name[capitals$country.etc == 'USA'],
  grepl("United Kingdom", country) ~ capitals$name[capitals$country.etc == 'UK'],
  grepl("St. Vincent & Grenadines", country) ~ capitals$name[capitals$country.etc == 'Saint Vincent and The Grenadines'],
  grepl("Fiji", country) ~ 'Suva',
  grepl("Montenegro", country) ~ 'Podgorica',
  grepl("Palestinian Territories", country) ~ 'Jerusalem',
  grepl("Kosovo", country) ~ 'Pristina',
  !is.na(capital) ~ capital
))
names_df %>% filter(is.na(capital)) # should be 0 rows

for(i in 1:nrow(names_df)){
  if(is.na(names_df$cluster[i])){names_df$cluster[i] <- 0}
  long <- unlist(capitals %>% filter(name == names_df$capital[i]) %>% 
    select(long))
  lat <- unlist(capitals %>% filter(name == names_df$capital[i]) %>% 
    select(lat))
  names_df$longitude[i] <- ifelse(length(long)==0, NA, as.numeric(long))
  names_df$latitude[i] <- ifelse(length(lat)==0, NA, as.numeric(lat))
}

names_df[names_df$capital=='Suva',]$longitude <- 178.42
names_df[names_df$capital=='Suva',]$latitude <- 18.14
names_df[names_df$capital=='Podgorica',]$longitude <- 19.26
names_df[names_df$capital=='Podgorica',]$latitude <- 42.43
names_df[names_df$capital=='Pristina',]$longitude <- 21.17
names_df[names_df$capital=='Pristina',]$latitude <- 42.67


names_df %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=latitude, 
                 color=as.factor(cluster),
                 shape=as.factor(cluster == 0)), size=3) +
  geom_text(aes(x=longitude, y=latitude, label=country),
            nudge_y = 2) +
  theme_bw() + coord_fixed(ratio = 1) +
  scale_color_brewer(palette = 'Set1') 
# ggsave("distances.png",
#        width=76,height=30,units="cm")

names_df %>% filter(cluster>0) %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=latitude, 
                 color=as.factor(cluster)), size=3) +
  theme_minimal() + coord_fixed(ratio = 1) +
  scale_color_brewer(palette = 'Set1')


## CLASSIFICATION

flu_ITZ_clusters <- flu_ITZ_clusters %>% mutate(capital = NA, longitude = NA, latitude = NA)
for(i in 1:nrow(flu_ITZ_clusters)){
  if(nrow(capitals %>% 
          filter(country.etc %in% c(flu_ITZ_clusters$country[i])) %>% 
          select(name)) > 1){
    flu_ITZ_clusters$capital[i] <- unlist(capitals %>% 
                                    filter(country.etc == flu_ITZ_clusters$country[i]) %>% 
                                    arrange(desc(pop)) %>% 
                                    select(name))[1]
  }
  if(nrow(capitals %>% 
          filter(country.etc %in% c(flu_ITZ_clusters$country[i])) %>% 
          select(name)) == 1){
    flu_ITZ_clusters$capital[i] <- unlist(capitals %>% 
                                    filter(country.etc %in% c(flu_ITZ_clusters$country[i])) %>% 
                                    select(name))
  }
  if(nrow(capitals %>%
          filter(country.etc %in% c(flu_ITZ_clusters$country[i])) %>%
          select(name)) == 0){
    vec <- rep(0, nrow(capitals))
    for(j in 1:nrow(capitals)){
      if(grepl(capitals$country.etc[j], flu_ITZ_clusters$country[i])){
        vec[j] <- j
      }
    }
    if(sum(vec)>0){
      capital_name <- capitals$name[vec[vec > 0]]
      flu_ITZ_clusters$capital[i] <- capital_name
    }
  }
}
flu_ITZ_clusters <- flu_ITZ_clusters %>% filter(method == 'kmeans_cluster') %>% mutate(capital = case_when(
  grepl('Ivoir', country) ~ capitals$name[capitals$country.etc == 'Ivory Coast'],
  grepl('Czechia', country) ~ capitals$name[capitals$country.etc == 'Czech Republic'],
  grepl("North Korea", country) ~ capitals$name[capitals$country.etc == 'Korea North'],
  grepl("Republic of Korea", country) ~ capitals$name[capitals$country.etc == 'Korea South'],
  grepl("United States", country) ~ capitals$name[capitals$country.etc == 'USA'],
  grepl("United Kingdom", country) ~ capitals$name[capitals$country.etc == 'UK'],
  grepl("Viet Nam", country) ~ capitals$name[capitals$country.etc == 'Vietnam'],
  grepl("St. Vincent & Grenadines", country) ~ capitals$name[capitals$country.etc == 'Saint Vincent and The Grenadines'],
  grepl("Fiji", country) ~ 'Suva',
  grepl("Montenegro", country) ~ 'Podgorica',
  grepl("Palestinian Territories", country) ~ 'Jerusalem',
  !is.na(capital) ~ capital
))

for(i in 1:nrow(flu_ITZ_clusters)){
  long <- unlist(capitals %>% filter(name == flu_ITZ_clusters$capital[i]) %>% 
                   select(long))
  lat <- unlist(capitals %>% filter(name == flu_ITZ_clusters$capital[i]) %>% 
                  select(lat))
  flu_ITZ_clusters$longitude[i] <- ifelse(length(long)==0, NA, as.numeric(long))
  flu_ITZ_clusters$latitude[i] <- ifelse(length(lat)==0, NA, as.numeric(lat))
}

flu_ITZ_clusters[flu_ITZ_clusters$capital=='Suva',]$longitude <- 178.42
flu_ITZ_clusters[flu_ITZ_clusters$capital=='Suva',]$latitude <- 18.14

cluster_means <- data.frame(cluster = unique(flu_ITZ_clusters$cluster_name),
                            long_mean = NA, lat_mean = NA)
for(i in 1:nrow(cluster_means)){
  xx <- flu_ITZ_clusters %>% filter(cluster_name == cluster_means$cluster[i],
                            !is.na(longitude))
  cluster_means$long_mean[i] <- mean(xx$longitude)
  cluster_means$lat_mean[i] <- mean(xx$latitude)
}

# must also link to Oceania-Melanesia-Polynesia by wrap-around
cluster_means <- rbind(cluster_means, c('Oceania-Melanesia-Polynesia', - 360 + 157.78500, 
                                        -28.7900000))
cluster_means$long_mean <- as.numeric(cluster_means$long_mean)
cluster_means$lat_mean <- as.numeric(cluster_means$lat_mean)

names_df %>% filter(cluster>0) %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=latitude,
                 color=as.factor(cluster)), size=3) +
  scale_color_manual(values = cluster_colors2) + 
  geom_point(data = cluster_means[abs(cluster_means$long_mean)<180,], 
             aes(x=long_mean, y=lat_mean, color=as.factor(cluster)), size=10, shape=4) +
  theme_bw() + coord_fixed(ratio = 1) 

## function to assign cluster which minimises euclidean distance
min_dist_cluster <- function(long, lat){
  if(is.na(long)){return(0)}
  minimiser <- cluster_means %>% 
    mutate(dist = sqrt((long_mean-long)^2 + 
                         (lat_mean-lat)^2)) %>% arrange(dist)
  return(minimiser$cluster[1])
}

names_df <- names_df %>% mutate(new_cluster = cluster)
for(i in 1:nrow(names_df)){
  if(names_df$cluster[i] == 0){
    names_df$new_cluster[i] <- min_dist_cluster(long = names_df$longitude[i], 
                                                lat = names_df$latitude[i])
  }
}

scatter_plot <- names_df %>% filter(new_cluster>0) %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=latitude,
                 color=as.factor(new_cluster),
                 shape=as.factor(cluster == 0)), size=3) +
  geom_point(data = cluster_means[abs(cluster_means$long_mean)<180,], 
             aes(x=long_mean, y=lat_mean, color=as.factor(cluster)), 
             size=8, shape=4) +
  theme_bw() + coord_fixed(ratio = 1) + labs(shape='Source', color='ITZ') +
  scale_color_manual(values=cluster_colors2) + xlab('Longitude') + ylab('Latitude') +
  scale_shape_manual(values=c(19,17), labels = c('Chen et al.', 'Addition')) +
  theme(axis.title.x = element_markdown(size = 14),
        axis.title.y = element_markdown(size = 14),
        text=element_text(size=14)); scatter_plot
ggsave("SM_plots/reclustered.png",
       width=38,height=15,units="cm")


# Validation on pre-clustered countries - would any be re-clustered?

names_df <- names_df %>% mutate(validation = NA)
for(i in 1:nrow(names_df)){
  names_df$validation[i] <- min_dist_cluster(long = names_df$longitude[i], 
                                              lat = names_df$latitude[i])
}

names_df %>% filter(!new_cluster == validation) %>% select(country, new_cluster, validation)

names_df %>% filter(!new_cluster == validation) %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=latitude,
                 color=as.factor(new_cluster)), size=3) +
  geom_point(data = cluster_means[abs(cluster_means$long_mean)<180,], 
             aes(x=long_mean, y=lat_mean, color=as.factor(cluster)), size=10, shape=4) +
  geom_text(aes(x=longitude, y=latitude, label=country),
            nudge_y = -4) +
  theme_bw() + coord_fixed(ratio = 1) +
  scale_color_brewer(palette = 'Set1') 


## printing the table for SM

print_df <- names_df %>% mutate(source = case_when(
  in_cluster == T ~ 'Chen et al.',
  in_cluster == F ~ 'Addition'
)) %>% filter(!new_cluster == 0) %>% select(country, new_cluster, source)

# library(googlesheets4)
# print_df %>% write_sheet(
#   ss = gs4_get(
#     "https://docs.google.com/spreadsheets/d/126D0R00iIu4GkODcyrS2bNl4A8ICcyqj7Rydg9IZADM/edit#gid=0"
#   ),
#   sheet = "Clusters"
# )


names_df_save <- names_df %>% select(codes, country, country_altern, country_altern_2, new_cluster, in_cluster) %>% 
  rename(cluster_name = new_cluster)

write_csv(names_df_save, file = paste0("data/new_clustering.csv"))
clusters <- read_csv("data/new_clustering.csv")

names_df_print <- data.table(names_df)[, c('codes','country','in_cluster','new_cluster')]
setnames(names_df_print, 'codes','ISO3C Code')
setnames(names_df_print, 'country','Name')
setnames(names_df_print, 'new_cluster','ITZ')
names_df_print[in_cluster==T, Source:='Chen et al.']
names_df_print[in_cluster==F, Source:='Addition']
names_df_print[,in_cluster:=NULL]

write_csv(names_df_print, file = paste0("data/itz_print.csv"))

## plotting as a map

require(maps)

world_map <- map_data("world") %>% rename(country = region) %>% 
  mutate(codes = countrycode(country, origin = 'country.name', destination = 'iso3c'))
cluster_map <- left_join(clusters, world_map, by = "codes")

cluster_map %>% filter(in_cluster == T) %>% 
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cluster_name), color = "white") +
  theme_bw() + scale_fill_manual(values=cluster_colors2) +
  xlab('Longitude') + ylab('Latitude') + labs(fill='ITZ') +
  theme(text=element_text(size=14))
ggsave("SM_plots/cluster_map_OG.png",
       width=50,height=20,units="cm")

full_map <- ggplot(cluster_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cluster_name), color = "white") +
  theme_bw() + scale_fill_manual(values=cluster_colors2) +
  xlab('Longitude') + ylab('Latitude') + labs(fill='ITZ') +
  theme(text=element_text(size=14), legend.position='none')
# ggsave("SM_plots/cluster_map.png",
#        width=50,height=20,units="cm")

library(patchwork) 
scatter_plot + full_map + plot_layout(guides='collect',nrow=2,heights=c(10,11)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')
ggsave("SM_plots/FigS1.png",
              width=30,height=25,units="cm")


cluster_colors2 <- c("Asia-Europe"="#21908CFF", "Southern America" = '#65156EFF', "Oceania-Melanesia-Polynesia" = '#CC4678FF', 
                    "Europe" = '#ffeba4', "Eastern and Southern Asia" = '#51C56AFF', "Northern America" = '#F89441FF', 
                    "Africa" = '#31688EFF')  
library(data.table)
cluster_map <- data.table(cluster_map)
cluster_map <- cluster_map[,exemplar := (codes %in% c('ARG','AUS','CAN','CHN','GBR','TUR','GHA'))]
centroids <- data.table(codes = c('ARG','AUS','CAN','CHN','GBR','TUR','GHA'),
                        long = c(-65.17328,135,-105,100,-2,35.667,-1.2),
                        lat = c(-35.38541,-26,58,35,53,39.167,8))
  
map <- ggplot(cluster_map, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = cluster_name), lwd=0.2, color='black') +
  geom_point(data = centroids, aes(x = long, y = lat, group=codes), fill = "white", 
             color='black', shape=21, size=2.5) + 
  theme_void() + scale_fill_manual(values=cluster_colors2) +
  xlab('Longitude') + ylab('Latitude') + 
  guides(color=NULL) +
  labs(fill='Influenza Transmission Zone') +
  theme(text=element_text(size=14)); map


design <- "
  1111
  1111
  1111
  2222
  2222
  2222
  2222
"

map + values +
  plot_layout(design = design) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')
ggsave(paste0("SM_plots/fig2.png"),
       width=30,height=26,units="cm")






