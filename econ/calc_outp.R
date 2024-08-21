
## CALCULATING OUTPATIENT/HOSPITALISATION RATIOS
#setwd("~/Desktop/research asst/Global Code")

library(data.table)
library(readr)
library(ggplot2)

rates <- data.table(read_csv("econ/outcome_calculations/data/lrti_outpat.csv", show_col_types=F))

## determining distribution
f.gamma <- function(shape, rate, x) {
  p <- pgamma(x, shape, rate)
  # return both
  return(c(p))
}
delta <- function(fit, actual) sum((fit-actual)^2)
objective <- function(theta, x, prob, ...) {
  ab <- (theta)
  fit <- f.gamma(ab[1], ab[2], x=as.numeric(x),...)
  # fit <- f.beta(ab[1], ab[2], x=as.numeric(x),...)
  return (delta(fit, prob))
}
fcn_fitting <- function(rates,
                        probs){

  x <- c(unlist(unname(rates)))
  sol <- suppressWarnings(optim(f=objective,p=c(1,1),
                                # method="BFGS",
                                x=x,
                                prob=c(probs),
                                control = list(reltol = 1e-15)
  ))
  parms <- (sol$par)
  return(parms)
}

rates[,3:5] <- rates[,3:5]/100000000

for(i in 1:nrow(rates)){
  parms <- fcn_fitting(rates[i,3:5], c(0.5, 0.025, 0.975))
  rates[i,"shape"] <- parms[1]
  rates[i,"rate"] <- parms[2]
  rates[i,"med_fit"] <- qgamma(p=c(0.5), shape=parms[1], rate=parms[2])
  rates[i,"l95_fit"] <- qgamma(p=c(0.025), shape=parms[1], rate=parms[2])
  rates[i,"u95_fit"] <- qgamma(p=c(0.975), shape=parms[1], rate=parms[2])
}

# ggplot(rates) + 
#   geom_point(aes(x=U95, y=u95_fit)) +
#   geom_line(aes(x=U95,y=U95),lty=2)

ratio_samples <- data.table(
  region = rep(unique(rates$region),each=100),
  simulation_index = 1:100
)

for(i in 1:nrow(ratio_samples)){
  shape <- rates[region==ratio_samples[i,region] &
                 outcome=='hospitalisation',]$shape 
  rate <- rates[region==ratio_samples[i,region] &
                   outcome=='hospitalisation',]$rate
  ratio_samples[i, hosp := 100000000*rgamma(1,shape,rate)]
  shape <- rates[region==ratio_samples[i,region] &
                   outcome=='outpatient',]$shape 
  rate <- rates[region==ratio_samples[i,region] &
                  outcome=='outpatient',]$rate
  ratio_samples[i, outp := 100000000*rgamma(1,shape,rate)]
}

# ggplot(ratio_samples) + 
#   geom_histogram(aes(x=hosp),bins=30) + 
#   facet_wrap(region~., scales='free')

ratio_samples[, ratio := outp/hosp]

ggplot(ratio_samples) +
  geom_density(aes(x=ratio)) +
  facet_wrap(region~., scales='free')

regions <- data.table(read_csv("econ/outcome_calculations/data/GBDsuperregions.csv", show_col_types=F))
itzs <- data.table(read_csv("econ/outcome_calculations/data/new_clustering.csv", show_col_types=F))
itzs[,c('cluster_name','in_cluster'):=NULL]

setnames(regions, 'val','region')
for(i in 1:nrow(itzs)){
  names <- c(unname(unlist(itzs[i,c('country', 'country_altern','country_altern_2')])))
  dt <- regions[country %in% names] 
  if(nrow(dt)==1){
    itzs[i,region := dt$region]
  }
}
setnames(itzs,'codes','country_code')

outpatient_ratios <- itzs[,c('country_code','region')]
outpatient_ratios <- outpatient_ratios[,lapply(.SD,rep,each=100)]
outpatient_ratios[, simulation_index := rep(1:100, 186)]
outpatient_ratios <- outpatient_ratios[ratio_samples, on=c('region','simulation_index'),
                                       ratio := ratio]
write_csv(outpatient_ratios, file = "econ/outcome_calculations/data/outpatient_ratios.csv")



