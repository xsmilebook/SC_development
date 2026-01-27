## This script is to generate derivative and posterior derivative values for gam models.
## Age will be sampled from 8 to 22, with a total of 1000 data points. 
## Covariables will be set as median or mode.
## The predicted derivatives will be generated.
## The data will be used to draw developmental trajectories and analyse the developmental alignment.
library(tidyverse)
library(R.matlab)
library(psych)
library(gratia)
library(mgcv)
library(parallel)

rm(list = ls())
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2
# set path
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ChineseCohort'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Figure_ChineseCohort_final'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/results_ChineseCohort'
}else{
  # in PC
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ChineseCohort'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ChineseCohort_final'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ChineseCohort'
}

#### load data
CVthr = 75
gamresultsum<-readRDS(paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_CV', CVthr,'_scale_TRUE.rds'))
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution,"_CV", CVthr,"_sumSCinvnode.sum.msmtcsd.combatAge.rds"))
SCrank <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    SCrank[x,y] = x^2 + y^2
  }
}
SCrank <- SCrank[lower.tri(SCrank, diag = T)]
SCrank.df <- data.frame(parcel=paste0("SC.", 1:elementnum, "_h"), SCrank=SCrank)
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))

## derivative
#######################################
Agevector <- SCdata.sum.merge$Age
# set cluster
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
clusterEvalQ(cl, {
  library(tidyverse)
  library(R.matlab)
  library(psych)
  library(gratia)
  library(mgcv)
  source(paste0(functionFolder, '/gamderivatives.R'))
})

derivative.sum<-parLapply(cl, 1:nrow(gamresultsum), function(x){
  SClabel.tmp <- gamresultsum$parcel[x]
  modobj<-gammodelsum[[x]]
  draws<-1
  increments<-1000
  derivdata<- gam.derivatives(modobj, "Age",Agevector, draws, increments, return_posterior_derivatives = FALSE)
  derivdata$label_ID<-gamresultsum$parcel[x]
  meanSC <- mean(modobj$model[,SClabel.tmp],na.rm=T)
  derivdata$meanSC<-meanSC
  return(derivdata)
})

derivative.df<-do.call(rbind, lapply(derivative.sum, function(x) data.frame(x)))
saveRDS(derivative.df, paste0(resultFolder, '/derivative.df_Yeo', elementnum, '_CV', CVthr,'.rds'))

## posterior derivative
############################
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
derivative.posterior.sum<-parLapply(cl, 1:nrow(gamresultsum), function(x){
  modobj<-gammodelsum[[x]]
  draws<-1000
  increments<-1000
  SClabel<-gamresultsum$parcel[x]
  derivdata<- gam.derivatives(modobj, "Age",Agevector,draws, increments, return_posterior_derivatives = TRUE)
  derivdata$SCrank<-SCrank.df$SCrank[SCrank.df$parcel==SClabel]
  meanSC <- mean(modobj$model[,SClabel],na.rm=T)
  derivdata$meanSC<-meanSC
  maxSC <- max(modobj$model[,SClabel],na.rm=T)
  derivdata$maxSC<-maxSC
  return(derivdata)
})
saveRDS(derivative.posterior.sum, paste0(resultFolder, '/derivative.posterior.df.Yeo', Yeoresolution, '_CV', CVthr,'.rds'))

