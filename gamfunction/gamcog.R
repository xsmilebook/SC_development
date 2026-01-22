library(tidyr)
library(mgcv)
library(psych)
library(dplyr)
library(ecostats)

#### Fit the linear effects of a continuous variable with a smooth item as covariates ####
## Function to fit evaluate the linear effects of a continuous variable on dependent variables per edge.
gam.fit.cognition <- function(region, dataname, cognition_var, smooth_var, covariates, knots,corrmethod, set_fx = FALSE, stats_only = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  cognition<-gam.data[ ,cognition_var]
  outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
  gam.data[outlierindx1 ,cognition_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,cognition_var]) & !is.na(gam.data[ ,region]))
  
  gam.data <- gam.data[NonNANIndex,]
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ %s + s(%s, k = %s, fx = %s) + %s",region, cognition_var, smooth_var, knots, set_fx, covariates))
  modelformula.null<-as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s",region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  
  modelformula.int <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, cognition_var, covariates))
  gam.model.int <- gam(modelformula.int, method="REML", data = gam.data)
  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.t <- gam.results$p.table[2,3]
  gam.smooth.pvalue <- gam.results$p.table[2,4]
  
  #Full versus reduced model anova p-value
  #anova.cov.pvalue <- anova.gam(gam.model.null,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  if (stats_only){
    anova.cov.pvalue <- tryCatch(
      {
        # Avoid internal parallelism / serialization issues on clusters.
        # Observed failure mode: `object 'nbinom2' not found` when ncpus > 1.
        anovaPB(gam.model.null, gam.model, n.sim = 1000, test = 'Chisq', ncpus = 1)$`Pr(>Chi)`[2]
      },
      error = function(e) {
        warning("anovaPB failed (set anova.cov.pvalue=NA): ", conditionMessage(e))
        NA_real_
      }
    )
  }else{anova.cov.pvalue <- NA}
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  anova.cov.int.pvalue <- anova.gam(gam.model.null,gam.model.int,test='Chisq')$`Pr(>Chi)`[2]
  if(is.na(anova.cov.int.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.int.pvalue <- 1}
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  sse.model.int <- sum((gam.model.int$y - gam.model.int$fitted.values)^2)
  partialRsq.int <- (sse.nullmodel - sse.model.int)/sse.nullmodel
  ### effect direction
  if(gam.smooth.t < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1
    partialRsq.int <- partialRsq.int*-1}
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  varcorformula2 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", cognition_var, smooth_var, knots, set_fx, covariates))
  res1<-residuals(gam(varcorformula1,method="REML", data = gam.data))
  res2<-residuals(gam(varcorformula2,method="REML", data = gam.data))
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-as.numeric(PCorr_Test$r)
  corrp <- as.numeric(PCorr_Test$p)
  
  stats.results <- cbind(parcel, cognition_var, gam.smooth.t, gam.smooth.pvalue,anova.cov.pvalue, anova.cov.int.pvalue,partialRsq.int, partialRsq, correstimate, corrp)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
}




