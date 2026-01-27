## ABCD dataset
## This script is to fit gam models for each connection.
## Statistical indexes and gam model files will be generated.  
rm(list=ls())
library(mgcv);
library(parallel)
library(tidyverse)
wdpath <- getwd()
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2
# set path
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
}else{
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
}

# input consistency threshold used to filter spurious streamlines
CVthr=75
SCdata.sum.merge<-readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
SCdata <- SCdata.sum.merge %>% drop_na(c(2:(1+elementnum)))
nrow(SCdata)
SCdata[,c("sex", "handness")] <- lapply(SCdata[,c("sex", "handness")], as.factor)
SCdata$age <- SCdata$age / 12
summary(SCdata[,c("sex", "handness", "mean_fd", "age")])
source(paste0(functionFolder, '/gammsmooth.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))

detectCores()
## calculate gam results
covariates<-"sex+mean_fd"
dataname<-"SCdata"
smooth_var<-"age"
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-grep("SC.",names(SCdata.sum.merge),value=T)[x]
    region<-SClabel
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
    gamresult<-as.data.frame(gamresult)
    return(gamresult)
  }, mc.cores = 50)
  
  gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  gamresultsum.df[,c(2:17)]<-lapply(gamresultsum.df[,c(2:17)], as.numeric)
  class(gamresultsum.df$gamm.smooth.pvalue)
  class(gamresultsum.df$partialRsq)
  # add fdr p values
  gamresultsum.df$pfdr <- p.adjust(gamresultsum.df$bootstrap_pvalue, method = "fdr")
  gamresultsum.df$sig<-(gamresultsum.df$pfdr<0.05)
  
  saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV',CVthr, '.rds'))
}else{
  gamresultsum.df <- readRDS(paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV',CVthr, '.rds'))
}

## calculate gam models
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- list()
  for (x in 1:elementnum){
    SClabel<-grep("SC.",names(SCdata.sum.merge),value=T)[x]
    region<-SClabel
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    resultsum[[x]] <- gamresult
  }
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV',CVthr,'.rds'))
}else{resultsum <- readRDS(paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV',CVthr,'.rds'))}

## plot data raw
#### generate function data
if (str_detect(wdpath, "cuizaixu_lab")){
  gammodelsum <- resultsum
  plotdatasum<-mclapply(1:elementnum, function(x){
    modobj<-gammodelsum[[x]]
    if (is.na(modobj[1])){
      plotdata<-NA
    }else{
      plotdata<- plotdata_generate(modobj, smooth_var="age")
      plotdata$SC_label <- names(plotdata)[14]
      plotdata[,14] <- NULL
    }
    return(plotdata)
  }, mc.cores = 2)
  plotdatasum.df <- do.call(rbind, lapply(plotdatasum, function(x) data.frame(x)))
  saveRDS(plotdatasum.df, paste0(interfileFolder, '/plotdatasum.df_Yeo', Yeoresolution, '_sumSCinvnode_siteall_CV',CVthr,'.rds'))
}else{plotdatasum.df <- readRDS(paste0(interfileFolder, '/plotdatasum.df_Yeo', Yeoresolution, '_sumSCinvnode_siteall_CV',CVthr,'.rds'))}

# To avoid the influence of averaged weight on derivative analyses, we divided SC strength of each edges by their
# weight at age of 8.9.
if (str_detect(wdpath, "cuizaixu_lab")){
  SCdata.diw <- SCdata
  for (i in 1:elementnum){
    SClabel <- grep("SC.",names(SCdata.sum.merge),value=T)[i]
    plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label==SClabel, ]
    SCdata.diw[ ,SClabel] <- SCdata[ ,SClabel] / plotdata.tmp$fit[1]
  }
  saveRDS(SCdata.diw, paste0(interfileFolder, "/SCdata.diw_Yeo", Yeoresolution, "CV", CVthr, ".rds"))
  SCdata.diw[,c("sex", "handness")] <- lapply(SCdata.diw[,c("sex", "handness")], as.factor)
  covariates<-"sex+mean_fd"
  dataname<-"SCdata.diw"
  smooth_var<-"age"
  resultsum <- list()
  for (x in 1:elementnum){
    SClabel<-grep("SC.",names(SCdata.sum.merge),value=T)[x]
    region<-SClabel
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    resultsum[[x]] <- gamresult
  }
  
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
  # gam results
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-grep("SC.",names(SCdata.sum.merge),value=T)[x]
    region<-SClabel
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
    gamresult<-as.data.frame(gamresult)
    return(gamresult)
  }, mc.cores = 50)
  
  gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  gamresultsum.df[,c(2:17)]<-lapply(gamresultsum.df[,c(2:17)], as.numeric)
  class(gamresultsum.df$gamm.smooth.pvalue)
  class(gamresultsum.df$partialRsq)
  summary(gamresultsum.df)
  saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
}
