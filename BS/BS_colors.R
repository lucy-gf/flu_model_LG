
## COLORS PALETTES AND LABELS FOR OUTPUTS
#setwd("~/Desktop/research asst/Global Code")

strain_colors <- c('tot' = 'black', 'totA' = '#7d66ac', 'totB' = '#e483a4')
strain_colors1 <- c('A' = '#7d66ac', 'B' = '#e483a4')
vt_colors <- c('1' = '#d91818', '2' = '#e2790e', '3' = '#eacb2c', '4' = '#62adc1', '5' = '#324da0')
age_colors <- c('0-5' = '#FDE725FF', '5-20' = '#5DC863FF', '20-65' = '#21908CFF', '65+' = '#3B528BFF', 'Total' = '#440154FF')
age_colors1 <- c('1' = '#FDE725FF', '2' = '#5DC863FF', '3' = '#21908CFF', '4' = '#3B528BFF', '5' = '#440154FF')
cluster_colors2 <- c("Asia-Europe"="#21908CFF", "Southern America" = '#65156EFF', "Oceania-Melanesia-Polynesia" = '#CC4678FF', 
                     "Europe" = '#FCFFA4FF', "Eastern and Southern Asia" = '#51C56AFF', "Northern America" = '#F89441FF', 
                     "Africa" = '#31688EFF')  
exemplar_colors <- c("TUR"="#21908CFF", "ARG" = '#65156EFF', "AUS" = '#CC4678FF', 
                     "GBR" = '#FCFFA4FF', "CHN" = '#51C56AFF', "CAN" = '#F89441FF', 
                     "GHA" = '#31688EFF')  

supp.labs <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs) <- c(1:5)
supp.labs1 <- c('Current','Improved (minimal)','Improved (efficacy)','Improved (breadth)','Universal')
names(supp.labs1) <- c(0:4)

supp.labs.age <- c('0-5','0-11','0-18','65+','0-18, 65+')
names(supp.labs.age) <- c(1:5)
supp.labs.cov <- supp.labs.age

supp.labs.agegrps <- c('0-4','5-19','20-64','65+')
names(supp.labs.agegrps) <- c(1:4)

supp.labs.strain <- c('Total','Influenza A','Influenza B')
names(supp.labs.strain) <- c('tot','totA','totB')

supp.labs.ITZ <- c("Africa", "Asia-Europe", "Eastern and\nSouthern Asia",
                   "Europe", "Northern\nAmerica", "Oceania-\nMelanesia-\nPolynesia",
                   "Southern\nAmerica")
names(supp.labs.ITZ) <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")
supp.labs2 <- c("Asia-\nEurope", "Africa", "Europe", "Southern\nAmerica",           
                "Oceania-\nMelanesia-\nPolynesia", "Eastern and\nSouthern Asia",  
                "Northern\nAmerica")
names(supp.labs2) <-  c('TUR','GHA','GBR','ARG','AUS','CHN','CAN')

supp.labs.country <- c("Ghana", "Turkey", "China","United\nKingdom", "Canada", "Australia","Argentina")
names(supp.labs.country) <- c("GHA", "TUR", "CHN", "GBR", "CAN", "AUS", "ARG")

supp.labs.ITZ2 <- c("Africa", "Asia-Europe", "Eastern and\nSouthern Asia",
                   "Europe", "Northern\nAmerica", "Oceania-\nMelanesia-\nPolynesia",
                   "Southern\nAmerica")
names(supp.labs.ITZ2) <- c("Africa", "Asia-Europe", "Eastern and Southern Asia",
                           "Europe", "Northern America", "Oceania-Melanesia-Polynesia",
                           "Southern America")










