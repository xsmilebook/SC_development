library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(ecostats)
set.seed(925)

#### PREDICT GAM SMOOTH FITTED VALUES FOR A SPECIFIED VALUE OF AN INTERACTING COVARIATE ####
## continuous interaction covariate
##Function to predict fitted values of a region for a given value of a covariate, using a varying coefficients smooth-by-linear covariate interaction
gam.smooth.predict.covariateinteraction <- function(region, dataname, smooth_var, int_var, int_var.predict.percentile, covariates, knots, set_fx = FALSE, increments, stats_only=TRUE){
  #Fit the gam
  gam.data <- get(dataname)
  int<-gam.data[ ,int_var]
  outlierindx1<-which(int<mean(int,na.rm = T)-3*sd(int,na.rm = T) | 
                      int>mean(int,na.rm = T)+3*sd(int,na.rm = T))
  gam.data[outlierindx1 ,int_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,int_var]) & !is.na(gam.data[ ,region]))
  gam.data <- gam.data[NonNANIndex,]
  int_var.predict<-quantile(gam.data[ ,int_var], c(int_var.predict.percentile))
  
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- region
  
  #Fit the gam
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+%5$s+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gam.model.null <- gam(modelformula.null, method = "REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  ##Full versus reduced model anova p-value
  #anova.int.pvalue <- anova(gam.model.null, gam.model,test='Chisq')$`Pr(>Chi)`[2]
  if (stats_only==TRUE){
    anova.int.pvalue <- anovaPB(gam.model.null,gam.model, n.sim = 1000,test='Chisq')$`Pr(>Chi)`[2]
  }else{anova.int.pvalue=NA}
  gam.int.pvalue <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = int_var),"p-value"][1]
  # interaction effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  Int.F <- gam.results$s.table[1,3]
  ### get predicted data
  ##############################
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on the modal number
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on the modal number
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  pred[,int_var] <- as.numeric(int_var.predict)
  
  #Generate fitted (predicted) values based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = gam.model, data = pred)
  predicted.smooth$fitted.centered <- scale(predicted.smooth$.fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth <- predicted.smooth %>% dplyr::select(all_of(smooth_var), .fitted, .se, .lower_ci, .upper_ci, fitted.centered)
  
  changed.range <- predicted.smooth$.fitted[which.max(predicted.smooth$age)]-predicted.smooth$.fitted[which.min(predicted.smooth$age)]
  changed.ratio <- predicted.smooth$.fitted[which.max(predicted.smooth$age)] / predicted.smooth$.fitted[which.min(predicted.smooth$age)]
  SCweight.age22 <- predicted.smooth$.fitted[which.max(predicted.smooth$age)]
  SCweight.age8 <- predicted.smooth$.fitted[which.min(predicted.smooth$age)]
  ############################################
  
  ### critical age in the case of different level of interaction variable
  ##############################################
  derv<-derivatives(gam.model, by=int_var, term = sprintf('s(%s)',smooth_var),
                    data=pred, interval = "simultaneous", unconditional = F)
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > .lower_ci & 0 < .upper_ci)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$.derivative*derv$sig
  
  # Age of increase onset and offset
  if (sum(derv$sig) > 0){
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
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$age[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA} 
  ############################################
  
  stats.reults<-cbind(parcel, int_var, anova.int.pvalue, gam.int.pvalue, IntpartialRsq,Int.F, changed.range, 
                      changed.ratio, SCweight.age22, SCweight.age8, increase.onset, increase.offset,
                      peak.change)

  full.results<-list(stats.reults, predicted.smooth)
  if(stats_only == TRUE)
    return(stats.reults)
  if(stats_only == FALSE)
    return(full.results)
}