library(MVMR)
library(glmnet)
library(boot)
library(ggplot2)
library(ggsci)


mvmr_ridge <- function (F.data, sig_list= 0){
  # takes input as the formated fdata
  
  ## Estimate instrument strength
  sres<-strength_mvmr(r_input=F.data,gencov=sig_list)
  ## Take a look for heterogenity (pleio) 
  pres<-pleiotropy_mvmr(r_input=F.data,gencov=sig_list)
  ## perform normal MVMR (IVW)
  res<-ivw_mvmr(r_input=F.data, gencov = sig_list)
  ## Ridge Regression
  ## Create matrix for later
  X <- as.matrix(F.data[, str_detect(colnames(F.data),"^(betaX)")] )
  Y <- as.matrix(F.data[,"betaYG"])
  
  colnames(Y) <- "betaYG"
  ### Get SNP weights
  weight<-1/F.data$sebetaYG^2
  
  # Using cross validation glmnet
  set.seed(123)
  ## Set range for lambdas
  lambdas <- 10^seq(-1, -5, by = -.1)
  lambdas<-c(lambdas)
  
  ##### IVW
  ## Do the Cross validation
  ridge_cv <- cv.glmnet(X,Y, weights=weight, alpha = 0, intercept=FALSE, lambda = lambdas) # by default with intercept
  ## Get the min Lambda
  best_lambda <- ridge_cv$lambda.min
  ### grab the labda resutls to plot later
  ridge <- glmnet(X,Y, weights=weight, alpha = 0, intercept=FALSE, lambda = lambdas) # by default with intercept
  ## get nice tidied version of the data 
  tidy_ridge <- broom::tidy(ridge)
  ## Create bootstrap matrix 
  bootmat <- as.data.frame(cbind(X, Y, weight, as.character(F.data$SNP)))
  
  ncolx <- ncol(X)
  
  colnames(bootmat)[ ncol(bootmat)] <- "SNP"
  
  #############
  ### bootstrap the model function (this would need to be changed depending on N var in model)
  samplebeta_IVW <- function(x, d, actual_lambda) {
    ridge_cv <-  glmnet(x=as.matrix(x[d,1:ncolx]),y= as.matrix(x[d,ncolx+1]), weights=x[d,ncolx+2], alpha = 0, intercept=FALSE, lambda = actual_lambda) 
    
    
    unlist(ridge_cv$beta[1:ncolx, ])

  }
  
  CI_IVW<-c()
  
  for (l in lambdas) { 
    
    b = boot(bootmat, samplebeta_IVW, R=5000, actual_lambda=l) 
    
    ci <- data.frame()
    
    for (m in 1:ncolx) {
      myrow <- c(l,paste0("betaX", m), boot.ci(b, type="perc", index=m)$perc[1,4],boot.ci(b, type="perc", index=m)$perc[1,5])
      ci <- rbind(ci, myrow)
    }
    colnames(ci) <-""
    
    
    CI_IVW<-rbind(CI_IVW,ci)
  }
  
  CI_IVW<-as.data.frame(CI_IVW)
  names(CI_IVW)<-c("lambda","term","LB","UB")
  if (all(CI_IVW$UB <= CI_IVW$LB)) names(CI_IVW)[3:4] <- c("UB", "LB") # name depending on which is bigger. ##SMH FML
  ## Merge through the confidence intervals with the point estimates
  results<-merge(tidy_ridge,CI_IVW, by = c("lambda", "term"))
  
  ## Make sure things are numeric
  results[, c("LB")] <- as.numeric(results[, c("LB")])
  results[, c("UB")] <- as.numeric(results[, c("UB")])
  
  
  ### Calculate the CI for these to check which 
  
  #### Make a plot to  see where the 'best lambda; sits.
  plotIVW <- ggplot(results,aes(x=log(lambda), y = estimate, color = term)) +
    geom_point(position=position_dodge(width=.2)) +
    geom_errorbar(aes(ymin=LB, ymax=UB), position=position_dodge(width=.2)) +
    labs(x="Lambda",y="Beta Estimate",title = "")+   
    theme_bw() +
    scale_color_npg() +
    geom_vline(xintercept = log(best_lambda), linetype="dotted") +
    geom_hline(yintercept = 0) 
  ##################
  ##################
  ## LOO
  CI_LOO_Ridge<-data.frame()
  
  for(i in  1:nrow(F.data)) {
    print(i)
    if(i == 1){
      ## get the mean estimate 
      loo_ridge_cv <- glmnet(bootmat[,1:ncolx],bootmat[, ncolx+1], weights=bootmat$weight, 
                             alpha = 0, intercept=FALSE, lambda = best_lambda) # by default with intercept
      
      b = boot(bootmat, samplebeta_IVW, R=1000, actual_lambda=best_lambda) 
      
      
      ci <- data.frame()
      
      for (m in 1:ncolx) {
        
        myrow <- c("Overall", l,paste0("betaX", m),loo_ridge_cv$beta[m,], boot.ci(b, type="perc", index=m)$perc[1,4], boot.ci(b, type="perc", index=m)$perc[1,5])
        ci <- rbind(ci, myrow)
      }
      
      colnames(ci) <- ""
      CI_LOO_Ridge<-rbind(CI_LOO_Ridge,ci)
    }
    ## create temp dataframe
    bootmat_temp <- bootmat[-i,]
    
    ## get the mean estimate 
    loo_ridge_cv <- glmnet(bootmat_temp[,1:ncolx],bootmat_temp[, ncolx+1], weights=bootmat_temp$weight, 
                           alpha = 0, intercept=FALSE, lambda = best_lambda) # by default with intercept
    
    b = boot(bootmat_temp, samplebeta_IVW, R=1000, actual_lambda=best_lambda) 
    
    
    ci <- data.frame()
    for ( m in 1:ncolx) {
      
      myrow <- c(bootmat$SNP[i], l,paste0("betaX", m),loo_ridge_cv$beta[m,], boot.ci(b, type="perc", index=m)$perc[1,4], boot.ci(b, type="perc", index=1)$perc[1,5])
      ci <- rbind(ci, myrow)
    }
    
    colnames(ci) <- ""
    CI_LOO_Ridge<-rbind(CI_LOO_Ridge,ci)
    
  }
  
  colnames(CI_LOO_Ridge) <- c("Model","lambda","term","beta", "UB", "LB")
  
  CI_LOO_Ridge[, 4:6] <- apply(CI_LOO_Ridge[, 4:6], 2, as.numeric)
  ### Plot 
  plotridge <- ggplot(CI_LOO_Ridge, aes(beta, Model, color = term) ) +
    geom_point(position=position_dodge(width=.3)) +
    geom_pointrange(aes(xmin=LB, xmax=UB),position=position_dodge(width=.3)) +
    labs(x="Beta Estimate (95CI%)",y="SNP",title = "")+   
    theme_bw() +
    scale_color_npg() +
    geom_vline(xintercept = 0, linetype="dotted") 
  return(list(results=results, CI_LOO_Ridge= CI_LOO_Ridge, CI_IVW= CI_IVW, plotridge=plotridge, plotIVW=plotIVW, best_lambda= best_lambda))
}
















