# Modified based on V. J. Sydnor et al. (2023)
library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
# NOTE: historical code used ecostats::anovaPB() (parametric bootstrap) to
# compare full vs reduced models. That path is fragile on the cluster and can
# spawn parallel backends unexpectedly. We instead use mgcv's built-in model
# comparison p-value (approximate) via anova.gam().

#### FIT GAM SMOOTH ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each edge and save out statistics and derivative-based characteristics.
gam.fit.smooth <- function(region, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  np = 1000
  theseVars <- attr(gam.model$terms,"term.labels")
  varClasses <- attr(gam.model$terms,"dataClasses")
  thisResp <- as.character(gam.model$terms[[2]])
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
  
  # GAM derivatives
  # NOTE: gratia output column names differ by version. Normalize to a common schema.
  derv_raw <- derivatives(gam.model, term = sprintf("s(%s)", smooth_var), interval = "simultaneous", unconditional = FALSE, data = pred)
  derv2_raw <- derivatives(gam.model, order = 2, term = sprintf("s(%s)", smooth_var), type = "central", unconditional = FALSE, data = pred)

  normalize_deriv <- function(df, smooth_var) {
    df <- as.data.frame(df)
    age_col <- if ("age" %in% names(df)) "age" else if ("data" %in% names(df)) "data" else smooth_var
    deriv_col <- if (".derivative" %in% names(df)) ".derivative" else if ("derivative" %in% names(df)) "derivative" else NA_character_
    se_col <- if (".se" %in% names(df)) ".se" else if ("se" %in% names(df)) "se" else NA_character_
    lower_col <- if (".lower_ci" %in% names(df)) ".lower_ci" else if ("lower" %in% names(df)) "lower" else NA_character_
    upper_col <- if (".upper_ci" %in% names(df)) ".upper_ci" else if ("upper" %in% names(df)) "upper" else NA_character_
    if (is.na(deriv_col) || is.na(se_col) || is.na(lower_col) || is.na(upper_col)) {
      stop("Unsupported gratia::derivatives() output columns: ", paste(names(df), collapse = ", "))
    }
    out <- data.frame(
      age = df[[age_col]],
      .derivative = df[[deriv_col]],
      .se = df[[se_col]],
      .lower_ci = df[[lower_col]],
      .upper_ci = df[[upper_col]]
    )
    out
  }

  derv <- normalize_deriv(derv_raw, smooth_var)
  derv2 <- as.data.frame(derv2_raw)
  derv2_col <- if (".derivative" %in% names(derv2)) ".derivative" else if ("derivative" %in% names(derv2)) "derivative" else names(derv2)[grepl("derivative", names(derv2))][1]
  if (is.na(derv2_col) || !nzchar(derv2_col)) stop("Unsupported gratia 2nd derivative output columns: ", paste(names(derv2), collapse = ", "))
  meanderv2 <- mean(derv2[[derv2_col]])
  
  #Identify derivative significance window(s
  derv <- derv %>%
    mutate(sig = !(0 > .lower_ci & 0 < .upper_ci))
  derv$sig_deriv <- derv$.derivative * derv$sig

  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  gam.edf <- gam.results$s.table[1]
  #Calculate the magnitude and significance of the smooth term effect by comparing full and reduced models
  ##Compare a full model GAM (with the smooth term) to a nested, reduced model (with covariates only)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, covariates)) #no smooth term
  gam.nullmodel <- gam(nullmodel, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ## Full versus reduced model p-value (mgcv approximate test; avoids bootstrap parallelism).
  if (mod_only == FALSE) {
    a <- anova(gam.nullmodel, gam.model, test = "Chisq")
    pcol <- intersect(c("Pr(>Chi)", "Pr(>F)"), colnames(a))
    if (length(pcol) == 0) stop("Unexpected anova() output columns: ", paste(colnames(a), collapse = ", "))
    anova.smooth.pvalue <- a[[pcol[[1]]]][2]
  } else {
    anova.smooth.pvalue <- NA
  }
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  anova.smooth.zvalue = qnorm(anova.smooth.pvalue / 2, lower.tail=FALSE)
  varcorformula1 <- as.formula(sprintf("%s ~ %s", region, covariates))
  varcorformula2 <- as.formula(sprintf("%s ~ %s", smooth_var, covariates))
  res1<-residuals(lm(varcorformula1, data = gam.data))
  res2<-residuals(lm(varcorformula2, data = gam.data))
  PCorr_Test <- cor.test(res1, res2)
  correstimate<-PCorr_Test$estimate
  corrp <- PCorr_Test$p.value
  ### effect direction
  mean.derivative <- mean(derv$.derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    partialRsq <- partialRsq*-1
    anova.smooth.zvalue<-anova.smooth.zvalue*-1}
  
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
  
  #Age of maximal increasing change
  if(sum(derv$sig) > 0){
    derv$sig_deriv_round<-round(derv$sig_deriv, 5)
    increasing.derv <- derv$sig_deriv_round[derv$sig_deriv > 0]
    if(length(increasing.derv) > 0){
      maxval.increase <- max(increasing.derv)
      window.peak.increase.change <- derv$age[derv$sig_deriv_round == maxval.increase]
      peak.increase.change <- mean(window.peak.increase.change)
    }else{peak.increase.change<-NA}}else{peak.increase.change<-NA}
    
  #Age of decrease onset and offset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$age[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0){
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
      decrease.offset <- max(decreasing.range)}
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
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.zvalue,
                        anova.smooth.pvalue,correstimate, corrp,gam.edf, meanderv2, change.onset, peak.change,peak.increase.change, 
                        decrease.onset, decrease.offset, increase.onset, 
                        increase.offset, change.offset)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq,anova.smooth.zvalue, anova.smooth.pvalue,correstimate, corrp,gam.edf, meanderv2)
  if(mod_only == TRUE)
    return(gam.model)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}
