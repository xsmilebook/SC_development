library(lme4)
library(mgcv)
library(gamm4)
library(gratia)
library(tidyverse)
library(dplyr)
library(pbkrtest)
set.seed(925)

pboot <- function(modelobj){
  numsims <- 1000
  
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  subID<-df$subID
  y <- df[,thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1|subID))
  m2 <- lmer(y ~ -1 + mat2 + (1|subID))
  refdist <- PBrefdist(m1, m2, nsim=numsims)#, cl=cl)
  pb <- PBmodcomp(m1, m2, ref = refdist)
  int_pval <- pb$test["PBtest","p.value"]
  
  return(int_pval)
}


#### FIT GAM SMOOTH ####
##Function to fit a GAMM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gamm.fit.smooth <- function(region, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = FALSE, resdata = FALSE){
  
  #Fit the gamm
  gam.data <- get(dataname)
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  #Extract gam input data
  df <- gamm.model$gam$model #extract the data used to build the gamm, i.e., a df of y + predictor values 
  np = 1000
  EPS <- 1e-07
  theseVars <- attr(gamm.model$gam$terms,"term.labels")
  varClasses <- attr(gamm.model$gam$terms,"dataClasses")
  thisResp <- as.character(gamm.model$gam$terms[[2]])
  #Create a prediction data frame, used to estimate (posterior) model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) {
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  
  #GAMM derivatives (gratia column names vary by version; normalize for stability)
  derv_raw <- derivatives(gamm.model$gam, term = sprintf("s(%s)", smooth_var), interval = "simultaneous", unconditional = FALSE, data = pred)
  derv_raw <- as.data.frame(derv_raw)
  age_col <- if ("age" %in% names(derv_raw)) "age" else if ("data" %in% names(derv_raw)) "data" else smooth_var
  deriv_col <- if (".derivative" %in% names(derv_raw)) ".derivative" else if ("derivative" %in% names(derv_raw)) "derivative" else NA_character_
  lower_col <- if (".lower_ci" %in% names(derv_raw)) ".lower_ci" else if ("lower" %in% names(derv_raw)) "lower" else NA_character_
  upper_col <- if (".upper_ci" %in% names(derv_raw)) ".upper_ci" else if ("upper" %in% names(derv_raw)) "upper" else NA_character_
  if (is.na(deriv_col) || is.na(lower_col) || is.na(upper_col)) {
    stop("Unsupported gratia::derivatives() output columns: ", paste(names(derv_raw), collapse = ", "))
  }
  derv <- data.frame(
    age = as.numeric(derv_raw[[age_col]]),
    derivative = as.numeric(derv_raw[[deriv_col]]),
    lower_ci = as.numeric(derv_raw[[lower_col]]),
    upper_ci = as.numeric(derv_raw[[upper_col]])
  )
  # Identify derivative significance window(s)
  derv <- derv %>% mutate(sig = !(0 > lower_ci & 0 < upper_ci))
  derv$sig_deriv <- derv$derivative * derv$sig
  derv$.derivative <- derv$derivative
  derv$.lower_ci <- derv$lower_ci
  derv$.upper_ci <- derv$upper_ci
  derv[[smooth_var]] <- derv$age
  
  derv2 <- derivatives(gamm.model$gam, order = 2, term = sprintf("s(%s)", smooth_var), type = "central", unconditional = FALSE, data = pred)
  derv2 <- as.data.frame(derv2)
  deriv2_col <- if (".derivative" %in% names(derv2)) ".derivative" else if ("derivative" %in% names(derv2)) "derivative" else NA_character_
  if (is.na(deriv2_col)) stop("Unsupported gratia::derivatives() output columns (order=2): ", paste(names(derv2), collapse = ", "))
  meanderv2 <- mean(as.numeric(derv2[[deriv2_col]]))

  #GAM statistics
  #F value for the smooth term and GAMM-based significance of the smooth term
  gamm.smooth.F <- gamm.results$s.table[3]
  gamm.smooth.pvalue <- gamm.results$s.table[4]
  gamm.edf <- gamm.results$s.table[1]
  
  #Send to bootstrap function
  if (mod_only==FALSE & resdata==FALSE){
    bootstrap_pvalue<-pboot(gamm.model) 
  }else{bootstrap_pvalue=1}
  
  #Calculate the magnitude and significance of the smooth term effect by comparing full and reduced models
  ##Compare a full model GAMM (with the smooth term) to a nested, reduced model (with covariates only)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, covariates)) #no smooth term
  gamm.nullmodel <- gamm4(nullmodel, random=~(1|subID), REML=FALSE, data = gam.data)
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.nullmodel$gam$y - gamm.nullmodel$gam$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  bootstrap.zvalue = qnorm(bootstrap_pvalue / 2, lower.tail=FALSE)
  df <- gamm.model$gam$model
  varcorformula1 <- as.formula(sprintf("%s ~ %s+%s", region, covariates, "(1|subID)"))
  varcorformula2 <- as.formula(sprintf("%s ~ %s+%s", smooth_var, covariates, "(1|subID)"))
  model.tmp <- lmer(varcorformula1, data = df)
  fit1 <- predict(model.tmp)
  res1<-residuals(lmer(varcorformula1, data = df))
  res2<-residuals(lmer(varcorformula2, data = df))
  PCorr_Test <- cor.test(res1, res2)
  correstimate<-PCorr_Test$estimate
  corrp <- PCorr_Test$p.value
  ### effect direction
  mean.derivative <- mean(derv$.derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    partialRsq <- partialRsq*-1
    bootstrap.zvalue<-bootstrap.zvalue*-1}
  
  #Derivative-based temporal characteristics
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    change.onset <- min(derv$age[derv$sig==T])} #find first age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ #if gam derivative is never significant
    change.onset <- NA} #assign NA
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$age[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA}  
  
  #Age of decrease onset and offset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$age[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0){
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
      decrease.offset <- max(decreasing.range)
    }
    if(length(decreasing.range) == 0){
      decrease.onset <- NA
      decrease.offset <-NA
    }}
  if(sum(derv$sig) == 0){
    decrease.onset <- NA
    decrease.offset <-NA}  
  
  
  #Age of increase onset and offset
  if(sum(derv$sig) > 0){ 
    increasing.range <- derv$age[derv$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
    if(length(increasing.range) > 0){
      increase.onset <- min(increasing.range)
      increase.offset <- max(increasing.range) #find oldest age with a significant positive derivative
      }
    if(length(increasing.range) == 0){
      increase.onset <- NA
      increase.offset <- NA
    }
  }
  if(sum(derv$sig) == 0){ 
    increase.onset <- NA
    increase.offset <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ 
    change.offset <- max(derv$age[derv$sig==T])} #find last age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ 
    change.offset <- NA}  
  
  full.results <- cbind(parcel, gamm.smooth.F, gamm.smooth.pvalue, partialRsq, bootstrap.zvalue,
                        bootstrap_pvalue,correstimate, corrp,gamm.edf, meanderv2, change.onset, peak.change, 
                        decrease.onset, decrease.offset, increase.onset, 
                        increase.offset, change.offset)
  stats.results <- cbind(parcel, gamm.smooth.F, gamm.smooth.pvalue, partialRsq, bootstrap.zvalue,
                         bootstrap_pvalue,correstimate, corrp,gamm.edf, meanderv2)
  if(mod_only == TRUE)
    return(gamm.model)
  if (resdata == TRUE)
    return(data.frame(SCres=res1, fit=fit1))
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}




