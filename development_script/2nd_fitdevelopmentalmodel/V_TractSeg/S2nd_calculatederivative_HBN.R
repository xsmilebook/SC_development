## This script is to generate derivative and posterior derivative values for gam models.
## Age will be sampled from 5 to 22, with a total of 1000 data points. 
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
# set resolution
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
# set path
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HBN'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Figure_HBN_final'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/results_HBN'
}else{
  # in PC
  demopath<-'D:/xuxiaoyu/open_dataset_information/HBN/demography'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HBN'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HBN_final'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HBN'
  source(paste0(functionFolder, '/colorbarvalue.R'))
}

#### load data
CVthr = 75
gamresultsum<-readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_CV', CVthr,'.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_CV', CVthr,'_scale_TRUE.rds'))
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_SA", ds.resolution,"_CV", CVthr,"_sumSCinvnode.TractSeg.combatage.rds"))
SCrank <- matrix(NA, ds.resolution, ds.resolution)
for (x in 1:ds.resolution){
  for (y in 1:ds.resolution){
    SCrank[x,y] = x^2 + y^2
  }
}
SCrank <- SCrank[lower.tri(SCrank, diag = T)]
SCrank.df <- data.frame(parcel=paste0("SC.", 1:elementnum, "_h"), SCrank=SCrank)
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))

## derivative
#######################################
agevector <- SCdata.sum.merge$age
# set cluster


derivative.sum<-mclapply(1:nrow(gamresultsum), function(x){
  SClabel.tmp <- gamresultsum$parcel[x]
  modobj<-gammodelsum[[x]]
  draws<-1
  increments<-1000
  derivdata<- gam.derivatives(modobj, "age",agevector, draws, increments, return_posterior_derivatives = FALSE)
  derivdata$label_ID<-gamresultsum$parcel[x]
  meanSC <- mean(modobj$model[,SClabel.tmp],na.rm=T)
  derivdata$meanSC<-meanSC
  return(derivdata)
}, mc.cores = 30)

derivative.df<-do.call(rbind, lapply(derivative.sum, function(x) data.frame(x)))
saveRDS(derivative.df, paste0(resultFolder, '/derivative.df', elementnum, '_TractSeg_CV', CVthr,'.rds'))

## posterior derivative
############################
derivative.posterior.sum<-mclapply(1:nrow(gamresultsum), function(x){
  modobj<-gammodelsum[[x]]
  draws<-1000
  increments<-1000
  SClabel<-gamresultsum$parcel[x]
  derivdata<- gam.derivatives(modobj, "age",agevector, draws, increments, return_posterior_derivatives = TRUE)
  derivdata$SCrank<-SCrank.df$SCrank[SCrank.df$parcel==SClabel]
  meanSC <- mean(modobj$model[,SClabel],na.rm=T)
  derivdata$meanSC<-meanSC
  maxSC <- max(modobj$model[,SClabel],na.rm=T)
  derivdata$maxSC<-maxSC
  return(derivdata)
}, mc.cores = 30)
saveRDS(derivative.posterior.sum, paste0(resultFolder, '/derivative.posterior.df.SA', ds.resolution, '_TractSeg_CV', CVthr,'.rds'))

