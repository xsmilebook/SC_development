# Modified based on V. J. Sydnor et al. (2023)
# generate 1000 fitted values at 1000 equispaced values of smooth_var for GAM or GAMM models.
# non-interested variables were set as median or mode.
plotdata_generate <- function(modobj,smooth_var){
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  if (!inherits(model, "gam")) {
    stop("model is not a mgcv::gam object; got class: ", paste(class(model), collapse = ","))
  }
  
  np <- 1000 #number of predicted values
  df = model$model
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  # line plot with no interaction
  thisPred <- data.frame(init = rep(0,np))
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisClass=="character"){
      df[,thisVar] <- as.factor(df[,thisVar])
      thisClass <- "factor"
    }
    
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

  # Ensure newdata matches training data classes/levels (especially factors).
  for (v in theseVars) {
    if (!v %in% names(df) || !v %in% names(pred)) next
    if (is.factor(df[[v]])) {
      pred[[v]] <- factor(as.character(pred[[v]]), levels = levels(df[[v]]))
    } else if (is.ordered(df[[v]])) {
      pred[[v]] <- factor(as.character(pred[[v]]), levels = levels(df[[v]]), ordered = TRUE)
    } else if (is.numeric(df[[v]])) {
      pred[[v]] <- as.numeric(pred[[v]])
    }
  }
  
  #pred$Sex<-levels(Behavior$Sex)[2]
  #pred$handnessfactor<-levels(Behavior$handnessfactor)[3]
  # NOTE: In some environments/models, `predict(..., se.fit=TRUE)` can fail with:
  # `lm object does not have a proper 'qr' component` (rank zero / qr=FALSE).
  # For robustness we fall back to `se.fit=FALSE` and keep CI columns as NA.
  # Also force the correct mgcv method to avoid accidental dispatch to other
  # predict() methods from attached packages.
  p <- tryCatch(data.frame(mgcv::predict.gam(model, pred, se.fit = TRUE)), error = function(e) NULL)
  if (is.null(p) || !all(c("fit", "se.fit") %in% names(p))) {
    fit_only <- as.numeric(mgcv::predict.gam(model, pred, se.fit = FALSE))
    p <- data.frame(fit = fit_only, se.fit = rep(NA_real_, length(fit_only)))
  }
  pred <- cbind(pred,p)
  pred$selo <- pred$fit - 1.96*pred$se.fit
  pred$sehi <- pred$fit + 1.96*pred$se.fit
  pred$fit.C <- scale(pred$fit, center = T, scale = F) # zero-centerd fitted values
  pred$fit.Z <- scale(pred$fit, center = T, scale = T) # z-scored fitted values
  pred$fit.floor <- pred$fit-pred$fit[1] # fitted values minus the initial fit
  pred$fit.ratio <- pred$fit / pred$fit[1] # fitted values divided by the initial fit
  pred$selo.C <- pred$fit - 1.96*pred$se.fit-mean(pred$fit)
  pred$sehi.C <- pred$fit + 1.96*pred$se.fit-mean(pred$fit)
  pred[,thisResp] = 1
  return(pred)
}
