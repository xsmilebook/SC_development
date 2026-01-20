#### DERIVATIVES ####
library(gratia)

##Function to compute smooth derivatives for a main GAM model and for individual draws from the simulated posterior distribution
gam.derivatives <- function(modobj,smooth_var, draws, increments, return_posterior_derivatives = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #extract the gam
  if (any(class(modobj)=="gam")) {
    gam.model <- modobj
  } else if (class(modobj$gam)=="gam") {
    gam.model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  region<-gam.model$terms[[2]]
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
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
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "character" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  pred2 <- pred #second prediction df
  pred2[,smooth_var] <- pred[,smooth_var] + EPS #finite differences
  
  # Estimate smooth derivatives
  # NOTE: gratia output column names differ by version. Normalize to a common schema.
  derivs_raw <- derivatives(gam.model, term = sprintf("s(%s)", smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL, data = pred)
  derivs_raw <- as.data.frame(derivs_raw)
  age_col <- if ("age" %in% names(derivs_raw)) "age" else if ("data" %in% names(derivs_raw)) "data" else smooth_var
  deriv_col <- if (".derivative" %in% names(derivs_raw)) ".derivative" else if ("derivative" %in% names(derivs_raw)) "derivative" else NA_character_
  se_col <- if (".se" %in% names(derivs_raw)) ".se" else if ("se" %in% names(derivs_raw)) "se" else NA_character_
  lower_col <- if (".lower_ci" %in% names(derivs_raw)) ".lower_ci" else if ("lower" %in% names(derivs_raw)) "lower" else NA_character_
  upper_col <- if (".upper_ci" %in% names(derivs_raw)) ".upper_ci" else if ("upper" %in% names(derivs_raw)) "upper" else NA_character_
  if (is.na(deriv_col) || is.na(se_col) || is.na(lower_col) || is.na(upper_col)) {
    stop("Unsupported gratia::derivatives() output columns: ", paste(names(derivs_raw), collapse = ", "))
  }

  derivs.fulldf <- data.frame(
    age = as.numeric(derivs_raw[[age_col]]),
    derivative = as.numeric(derivs_raw[[deriv_col]]),
    se = as.numeric(derivs_raw[[se_col]]),
    lower_ci = as.numeric(derivs_raw[[lower_col]]),
    upper = as.numeric(derivs_raw[[upper_col]])
  )
  derivs.fulldf$significant <- !(0 > derivs.fulldf$lower_ci & 0 < derivs.fulldf$upper)
  derivs.fulldf$significant.derivative <- derivs.fulldf$derivative * derivs.fulldf$significant
  colnames(derivs.fulldf)[1] <- sprintf("%s", smooth_var)
  
  #Estimate posterior smooth derivatives from simulated GAM posterior distribution
  if(return_posterior_derivatives == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
    X1 <- predict(gam.model, newdata = pred2, type = "lpmatrix") #get matrix of linear predictors for pred2
    Xp <- (X1 - X0) / EPS 
    posterior.derivs <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
    posterior.derivs <- as.data.frame(posterior.derivs)
    # posterior.derivs[1,1]<-sum(Xp[1,] * sims[1,]), posterior.derivs[1,2]<-sum(Xp[1,] * sims[2,])...
    # posterior.derivs is the change of posterior smooth term when there is an infinitesimal change of age
    # namely posterior derivatives
    colnames(posterior.derivs) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    posterior.derivs <- cbind(as.numeric(pred[,smooth_var]), posterior.derivs) #add smooth_var increments from pred df to first column
    colnames(posterior.derivs)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    posterior.derivs <- cbind(as.character(region), posterior.derivs) #add parcel label to first column
    colnames(posterior.derivs)[1] <- "label_ID" #label the column
    posterior.derivs.long <- posterior.derivs %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.derivative")
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.derivative)
  
  if(return_posterior_derivatives == FALSE)
    return(derivs.fulldf)
  if(return_posterior_derivatives == TRUE)
    return(posterior.derivs.long)
}
