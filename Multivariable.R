# this script performs the steps of MVMR 

library(tidyverse)
library(TwoSampleMR)
library(MVMR)
library(ggplot2)

source("Ridge.R")
library(tidyr)

library(MendelianRandomization)



outcomes <- c( "HPV+", "HPV-")



################################################################################
# STEP1 :: Estimating pairwise covariances between SNP associations
################################################################################


# not required as these are two different samples

################################################################################
# STEP2: Read  Data
################################################################################


analysis_sets<- read_csv("analysis_setsm.csv")



# defining empty data frames for results
ivw<- data.frame()


strength<- data.frame()


pleiotropy <- list()


qhet <- data.frame()



outcome_stats <- list()

ridge <- list()

mv_intercepts <- data.frame()

for ( i in seq_along(outcomes)) {
  csi_dpw <- analysis_sets%>%filter(Outcome== outcomes[i])
  
  #restrict to 108 SNPs for CSI and 60 SNPs for DPW 
  
  
  
  XGs_betas <- csi_dpw[,c(1, 4:5)]
  colnames(XGs_betas)[2:3] <- str_remove_all(colnames(XGs_betas)[2:3] , "xg_")
  
  XGs_se <- csi_dpw[,c(1,6:7)]
  colnames(XGs_se)[2:3] <- str_remove_all(colnames(XGs_se)[2:3] , "xgse_")
  
  #Remove NAs 
  XGs_betas <- na.omit(XGs_betas)
  XGs_se <- na.omit(XGs_se)
  
  
  YG <- csi_dpw[, 1:3]
  colnames(YG)[2:3] <- c("yg", "ygse")
  
  XGs_betas <- XGs_betas[order(XGs_betas$SNP),]
  XGs_se <- XGs_se[order(XGs_se$SNP),]
  YG <- YG[order(YG$SNP),]
  
  
  
  mvmr <- format_mvmr(XGs_betas[,c(2:3)], YG$yg, XGs_se[,c(2:3)], YG$ygse, XGs_betas$SNP)
  mvmr_res <- mvmr(mvmr, 0, 1)
  
  ################################################################################
  # step5: perform MVMR
  ################################################################################
  
  #MVMR-Egger 
  mr_mvivw <- mr_mvivw(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), 
                                  bxse = cbind(XGs_se$DPW, XGs_se$CSI), 
                                  by = YG$yg, YG$ygse))
  #save(mr_mvivw, file="mr_mvivw_dpwandcsi_final.Rdata")
  
  mr_mvegger <- mr_mvegger(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI),
                                      bxse = cbind(XGs_se$DPW, XGs_se$CSI), 
                                      by = YG$yg, YG$ygse), 
                           orientate = 1)
  #save(mr_mvegger, file="mr_mvegger_dpwandcsi_final.Rdata")
  
  mvmr_results_DPW <- c(exp(mr_mvivw$Estimate[1]),
                        exp(mr_mvivw$CILower[1]),
                        exp(mr_mvivw$CIUpper[1]),
                        mr_mvivw$Pvalue[1], 
                        exp(mr_mvegger$Estimate[1]),
                        exp(mr_mvegger$CILower.Est[1]),
                        exp(mr_mvegger$CIUpper.Est[1]), 
                        mr_mvegger$Pvalue.Est[1]
                        )
  
  mvmr_results_DPW <- c(format(mvmr_results_DPW, scientific=F))
  names(mvmr_results_DPW) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")
  
  mvmr_results_CSI <- c(exp(mr_mvivw$Estimate[2]*0.6940093),
                        exp(mr_mvivw$CILower[2]*0.6940093), 
                        exp(mr_mvivw$CIUpper[2]*0.6940093), 
                        mr_mvivw$Pvalue[2], exp(mr_mvegger$Estimate[2]*0.6940093),
                        exp(mr_mvegger$CILower.Est[2]*0.6940093), 
                        exp(mr_mvegger$CIUpper.Est[2]*0.6940093),
                        mr_mvegger$Pvalue.Est[2]
                        )
  
  mvmr_results_CSI <- c(format(mvmr_results_CSI, scientific=F))
  names(mvmr_results_CSI) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")
  
  mvmr_results <- data.frame(t(cbind(mvmr_results_DPW, mvmr_results_CSI) ) )
  #to get the outcome per standard deviation increase then multiply beta and se by 0.6940093
  
  mv_inter <- t(c(mr_mvegger$Intercept[1],
                  mr_mvegger$CILower.Int[1],
                  mr_mvegger$CIUpper.Int[1],
                  mr_mvegger$Pvalue.Int[1],
                  outcomes[i]
                  )
                )
  
  
  colnames(mv_inter)<- c("Intercept", "LB", "UB","Pvalue" , "Outcome")  
  mv_intercepts<- rbind(mv_intercepts, mv_inter)
  
  mvmr_results$exposure <- c("DPW", "CSI")
  mvmr_results$outcome <-outcomes[i]
  
  ivw<- rbind(ivw, mvmr_results)
  
  
  
  
  ################################################################################
  # STEP3: Check instrument strength
  ################################################################################
  streng <- data.frame(t(strength_mvmr(mvmr, gencov = 0)))
  
  
  streng$outcome <- outcomes[i]
  
  
  
  
  
  
  strength <- rbind(strength, streng)
  
  
  
  
  
  ################################################################################
  #Step4: Check Pleiotropy
  ################################################################################
  pleiotropy[[i]] <- pleiotropy_mvmr(mvmr, gencov=0)
  names(pleiotropy)[i] <- outcomes[i]
  
  ################################################################################
  #Step5: Ridge
  ################################################################################
  

  ridge[[i]] <- mvmr_ridge(F.data = mvmr)
  names(ridge)[i] <-outcomes[i]
  
  
  
  
}
