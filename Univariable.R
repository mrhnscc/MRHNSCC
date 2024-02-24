


# Loading libraries
library(tidyverse)
library(TwoSampleMR)
library(MVMR)
library(ggplot2)
library(simex)

# reading the summary stats

analysis_sets<- read_csv("analysis_sets.csv")

counter<- 1

#main results
UnivariableMRResults<- data.frame()
presso_results2<- list()

presso_results<- data.frame()
# Univariable MR heterogeneity

heterogeneity <- data.frame()

# Univariable MR Isq
isqs <- data.frame() # data frame for Isq stats.

# Univarable MR Intercepts

intercepts <- data.frame()

# I-squared function for is SIMEX required for MR-Egger 
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

calculate_isq <- function(dat ) {
  #calculate Isq wieghted and unweighted
  I2<-c()
  
  #Rename required columns
  BetaXG <- dat$beta.exposure[dat$mr_keep=="TRUE"]
  BetaYG <- dat$beta.outcome[dat$mr_keep=="TRUE"]
  seBetaXG <- dat$se.exposure[dat$mr_keep=="TRUE"]
  seBetaYG <- dat$se.outcome[dat$mr_keep=="TRUE"]
  BXG             = abs(BetaXG)         # gene--exposure estimates are positive  
  
  # Calculate F statistics
  # and I-squared statistics
  # to measure Instrument 
  # strength for MR-Egger
  
  F   = BXG^2/seBetaXG^2
  mF  = mean(F)
  Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
  Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
  
  #Save results
  I2<-data.frame(cbind(mF, Isq_unweighted, Isq_weighted))
  colnames(I2) <- c("mF", "Isq_unweighted", "Isq_weighted")
  return(I2)
  
  
  
}

exposures<- levels(factor(analysis_sets$exposure))
outcomes<- levels(factor(analysis_sets$outcome))


for ( i in seq_along(exposures)) {
  
  
  for ( j in seq_along(outcomes)) {
  
    
      analysis_set <- analysis_sets %>%filter(exposure== exposures[i] & outcome== outcomes[j])
      
      # performing univariable MR
      Univariable <- mr(analysis_set, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median" , "mr_weighted_mode" ) )
      
      # heterogeneity (supp table 2)
      heterogen <- mr_heterogeneity(analysis_set)
      
      heterogeneity <- rbind(heterogeneity, heterogen)
      
      
      # isq (supp table 1)
      isq <- calculate_isq(analysis_set)
      isq$outcome <- outcomes[j]
      isq$exposure <- exposures[i]
      isqs <- rbind(isqs, isq)
      
      #  intercepts (supp table 3)
      inter <- mr_pleiotropy_test(analysis_set)
      intercepts <- rbind(intercepts, inter)
      
  
      
      presso <- MRPRESSO:: mr_presso(BetaOutcome = "beta.outcome",
                                     BetaExposure= "beta.exposure",
                                     SdOutcome= "se.outcome",
                                     SdExposure = "se.exposure",
                                     OUTLIERtest = T,
                                     DISTORTIONtest = T,
                                     data= as.data.frame(analysis_set))
      
      

      presso[[1]]$outcome <- outcomes[j]
      presso[[1]]$Exposure <- exposures[i]
      
      

      
      
      
      
      
      # storing the results 
      UnivariableMRResults <- rbind( UnivariableMRResults, Univariable)
      
      
      presso_results <- rbind(presso_results, presso[[1]])
      
      presso_results2[[counter]] <- presso
      counter <- counter+1
     
    

  }
}

