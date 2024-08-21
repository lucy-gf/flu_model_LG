
## CALCULATING YLLS PER DEATH
#setwd("~/Desktop/research asst/Global Code")

library(data.table) 

load("econ/outcome_calculations/data/WPP2022_Life_Table_Medium.csv")
UNLT  <- UNLT[MidPeriod == 2022.5 & Sex == "Total"]  # only need combined sexes for 2025

yll <- function(
    LT = NULL, # data table with UN life table (if using)
    qx = NULL, # vector of qx from a life table
    LT_ages = NULL, # vector of age group starting points corresponding to qx
    r = 0, # discount rate
    smr = 1, # standardised mortality rate adjustment for comorbidities
    weight_method = "lxqx", # weight method to average LE by age group: "lx" "lxqx" "equal" "pop_ifr"
    model_ages, # vector of age group lower bounds (may or may not include 0)
    func = NULL # function mapping age (1:100) to IFR, to weight to if weight_method=="pop_ifr"
){
  ## checking demographic inputs are reasonable
  if(is.null(LT) & is.null(qx)){stop('Must enter a qx vector or life table.')}
  if(is.null(LT) & !is.null(qx) & is.null(LT_ages)){stop('Must enter corresponding age vector.')}
  if(!is.null(qx) & !is.null(LT_ages) & (length(qx) != length(LT_ages))){
    stop('qx and LT_ages must be the same length')}
  
  ## building life table if using qx and not a life table
  if(!is.null(qx) & !is.null(LT_ages) & (length(qx)==length(LT_ages))){
    LT_out <- data.table(Age = LT_ages,
                         qx = qx)
    LT_out <- LT_out[, AgeGrpEnd := lead(Age)
    ][Age == 100, AgeGrpEnd := 101][,AgeGrpSpan := AgeGrpEnd - Age]
    
    # convert age group qx to estimated yearly qx
    LT_out[, qx := 1 - (1 - qx)^(1 / AgeGrpSpan)]
    # instantaneous death rate (Briggs et al)
    LT_out[, dx := -log((1 - qx))]
    
    # lx - number of people alive at the start of age x
    LT_out[Age == 0, lx := 100000] # starting radix
    for(a in 1:max(LT_out$Age)){
      LT_out[Age == a,
             lx := LT_out[Age == a - 1, lx] *        # alive at start of previous age group
               exp(-LT_out[Age == a - 1, dx] * smr)]  # not dead during previous age group
    }
    
    # Lx - assumed ax = 0.5 at all ages as no information on IMR etc.
    for (a in 0:(max(LT_out$Age) - 1)){
      LT_out[Age == a,
             Lx :=  0.5 * (lx + LT_out[Age == a + 1, lx])       
      ]
    }
    # Lx in final age group - USING BEST FIT VAL PREDICTOR FOR NOW...
    penul_qx <- LT_out[Age == sort(LT_out$Age)[nrow(LT_out)-1], qx]
    LT_out[Age == max(LT_out$Age), Lx := 0.6987219*lx/penul_qx]
  }
  
  ## matching format if input is life table, not qx
  if(is.null(qx) & !is.null(LT)){
    # finding the columns that we want
    age_cols <- c(which(grepl('age', names(LT), ignore.case	= T)))
    demog_cols <- c(which(grepl('qx', names(LT), ignore.case	= T)),
                    which(grepl('dx', names(LT), ignore.case	= T)),
                    which(grepl('lx', names(LT), ignore.case	= T)))
    cols <- c(age_cols, demog_cols)
    LT_out <- LT[, ..cols]
    
    ## renaming age column
    ## (must increase, but be the lowest)
    age_cols_new <- c(which(grepl('age', names(LT_out), ignore.case	= T)))
    age_col_maybe <- c()
    for(i in age_cols_new){
      # must be numeric
      if(sum(is.na(suppressWarnings(as.numeric(unlist(LT_out[, ..i])))))==0){
        vec_i <- as.numeric(unlist(LT_out[, ..i]))
        # must be increasing every time
        if(sum(vec_i[2:length(vec_i)] - vec_i[1:(length(vec_i) - 1)] > 0) == (length(vec_i) - 1)){
          age_col_maybe <- c(age_col_maybe, i)
        }
      }
    }
    if(length(age_col_maybe) > 1){
      age_col_maybe <- unname(which.min(LT_out[1,..age_col_maybe]))
    }
    if(length(age_col_maybe) == 0){stop("Please label a life table 'age' column")}
    setnames(LT_out, names(LT_out)[age_col_maybe], 'Age')
  }
  
  # discounted Tx
  for(a in 0:max(LT_out$Age)){
    LT_out[Age == a,
           dTx := sum(LT_out[Age >= a, Lx / (1 + r)^(Age - a)])]
  }
  
  # discounted LEx
  LT_out[, dLEx := dTx/lx]
  
  # Age groups for output 
  if(!0 %in% model_ages){model_ages <- c(0,model_ages)}
  model_ages_end <- c((model_ages - 1)[model_ages > 0], 100)
  AgeBands <- data.table(AgeBand = 1:length(model_ages),
                         start = model_ages,
                         end = model_ages_end)
  LT_out[, joinAge := Age]
  LT_out <- LT_out[AgeBands, on = .(joinAge >= start, joinAge <= end)][Age != 100] # drop age 100 -- WHY?
  
  # dLEx for age bands weighted by lx
  if (weight_method == "lx"){ # weighting based on people alive at that age in the life table
    
    out <- LT_out[,
                  .(d_LEx = sum(lx * dLEx)/sum(lx)),
                  by = .(AgeBand)][order(AgeBand)]
    
  } else if (weight_method == "lxqx"){# number of people alive and the risk of death at each age
    
    out <- LT_out[, 
                  .(d_LEx = sum(lx * qx * dLEx)/sum(lx * qx)),
                  by = .(AgeBand)][order(AgeBand)]
    
  } else if (weight_method == "equal") { # each age is weighted equally (i.e. what people usually use without thinking)
    
    out <- LT_out[, 
                  .(d_LEx = mean(dLEx)),
                  by=.(AgeBand)][order(AgeBand)]
    
  } else if (weight_method == "pop_ifr") { # actually population projections (UN World pop)
    
    LT_out[, ifr:=func(Age)]
    out <- LT_out[, 
                  .(d_LEx = sum(lx * ifr * dLEx)/sum(lx * ifr)),
                  by = .(AgeBand)][order(AgeBand)]
    
  } else {stop('Invalid argument for weight method')}
  
  out[, Age := paste0("[", AgeBands$start, ",", 
                      ifelse(AgeBands$end == max(AgeBands$end), '+',
                             AgeBands$end + 1), ")")]
  
  return(out) 
}

























