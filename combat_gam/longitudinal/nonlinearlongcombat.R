#Non-Linear LongCombat Function 
#Load Libraries as Dependencies (Also Depends on NeuroCombat.R script)
library(invgamma)
library(lme4)
library('pbkrtest')
library('mgcv')
library(gamm4)
library('matrixStats')
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(parallel)


nlongcombat=function(features,knot,ranef,X_nlin,X_lin,batch,use_random=TRUE,n_cores=1){
stand=c()
subid=ranef
maxi=matrix(0,nrow=nrow(features),ncol=ncol(features))
Fitted_site=matrix(0,nrow=nrow(features),ncol=ncol(features))
Residual=matrix(0,nrow=nrow(features),ncol=ncol(features))

fit_feature <- function(i) {
  y=features[,i]
  Xe=cbind(y,X_nlin,X_lin,batch,subid=ranef)
  Xe=as.data.frame(Xe)

  e=""
  for(j in 1:ncol(X_nlin)){
    e=paste(e,"+","s(",colnames(X_nlin)[j],",bs ='tp'",", k=",knot,", fx=TRUE)",sep="")
  }
  e=substring(e,2)

  #Do the same for linear coefficients 
  f=""
  for(z in 1:ncol(X_lin)){
    f=paste(f,"+",colnames(X_lin)[z],sep="")
  }
  f=substring(f,2)

if (use_random) {
  yur=gamm4(as.formula(paste('y','~',"as.factor(batch)",
                             "+",e,"+",f,"+1",sep='')),
            random=~(1|subid),data=data.frame(Xe))
  maxi_i=as.matrix(yur$gam$fitted.values) 
  mm=as.matrix(model.matrix(yur$gam))
  batch_cols=grepl("^as\\.factor\\(batch\\)", colnames(mm))
  fitted_site_i=mm[,batch_cols,drop=FALSE] %*% as.matrix(yur$gam$coefficients[batch_cols])
  stand_i=as.data.frame(VarCorr(yur$mer))[3,5]
} else {
  yur=gam(as.formula(paste('y','~',"as.factor(batch)",
                           "+",e,"+",f,"+1",sep='')),
          method="REML",data=data.frame(Xe))
  maxi_i=as.matrix(yur$fitted.values) 
  mm=as.matrix(model.matrix(yur))
  batch_cols=grepl("^as\\.factor\\(batch\\)", colnames(mm))
  fitted_site_i=mm[,batch_cols,drop=FALSE] %*% as.matrix(yur$coefficients[batch_cols])
  stand_i=sd(residuals(yur))
}
  if (is.na(stand_i) || stand_i == 0) {
    stand_i <- 1
  }
  residual_i=(y-(maxi_i-fitted_site_i))/stand_i
  residual_i[is.na(residual_i)] <- 0
  print(paste0('Fitting Feature',i,sep=' '))
  return(list(stand=stand_i, maxi=maxi_i, fitted_site=fitted_site_i, residual=residual_i))
}

if (n_cores > 1) {
  results <- mclapply(seq_len(ncol(features)), fit_feature, mc.cores=n_cores)
} else {
  results <- lapply(seq_len(ncol(features)), fit_feature)
}
stand <- vapply(results, function(x) x$stand, numeric(1))
maxi <- do.call(cbind, lapply(results, function(x) x$maxi))
Fitted_site <- do.call(cbind, lapply(results, function(x) x$fitted_site))
Residual <- do.call(cbind, lapply(results, function(x) x$residual))

data.harmonized=neuroCombat(t(Residual),as.factor(batch))
dat=t(data.harmonized$dat.combat)

for(i in 1:ncol(dat)){
  dat[,i]=dat[,i]*stand[i]
}
harmonized=dat+(maxi-Fitted_site)

return(list(harmonized=harmonized))
}


#---------------------------------
#Example Input and Function
#---------------------------------
#knot = 45
#ranef=simdata$subid
#X_nlin=as.matrix(simdata$age)
#colnames(X_nlin)='age'
#X_lin=cbind(simdata$time)
#colnames(X_lin)=c('time')
#batch=simdata$batch
#features=simdata[,6:25]
#colnames(features)=colnames(simdata)[6:25]


#dr=nlongcombat(features,knot,ranef,X_nlin,X_lin,batch)
  
  
  
  

  
  
