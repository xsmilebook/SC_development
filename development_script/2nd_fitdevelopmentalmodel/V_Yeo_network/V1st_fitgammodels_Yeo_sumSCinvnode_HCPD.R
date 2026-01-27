## This script is to fit gam models for each edge.
## Statistical indexes and gam model files will be generated.  
rm(list=ls())
library(mgcv)
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
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
}else{
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
}
# set which consistency threshold used to filter spurious streamlines
CVthr=75
# load data
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution,"_CV", CVthr,"_sumSCinvnode.sum.msmtcsd.combatage.rds"))
nrow(SCdata.sum.merge)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
SCdata.sum.merge[,2:(elementnum+1)] <- lapply(SCdata.sum.merge[,2:(elementnum+1)], as.numeric)
# source function
source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))
detectCores()
## calculate gam results
covariates<-"sex+mean_fd"
dataname<-"SCdata.sum.merge"
smooth_var<-"age"
# set cluster
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
invisible(clusterEvalQ(cl, {
  library(mgcv)
  library(tidyverse)
  source(paste0(functionFolder, '/gamsmooth.R'))
  source(paste0(functionFolder, '/plotdata_generate.R'))
}))

resultsum <- parLapply(cl, 1:elementnum, function(x){
  SClabel<-grep("SC.", names(SCdata.sum.merge), value=T)[x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
  gamresult<-as.data.frame(gamresult)
  return(gamresult)
})

gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df[,c(2:18)]<-lapply(gamresultsum.df[,c(2:18)], as.numeric)
class(gamresultsum.df$gam.smooth.pvalue)
class(gamresultsum.df$partialRsq)
# add fdr p values
gamresultsum.df$pfdr<- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
gamresultsum.df$sig<-(gamresultsum.df$pfdr<0.05)
saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'.rds'))

## calculate gam models
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-names(SCdata.sum.merge)[1+x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 20)
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'.rds'))
}else{resultsum <- readRDS(paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'.rds'))}

## plot data raw
#### generate function data
if (str_detect(wdpath, "cuizaixu_lab")){
  gammodelsum <- resultsum
  plotdatasum<-mclapply(1:elementnum, function(x){
    modobj<-gammodelsum[[x]]
    plotdata<- plotdata_generate(modobj, dataname=NA,  "age")
    plotdata$SC_label <- names(plotdata)[14]
    plotdata[,14] <- NULL
    return(plotdata)
  }, mc.cores = 50)
  plotdatasum.df <- do.call(rbind, lapply(plotdatasum, function(x) data.frame(x)))
  saveRDS(plotdatasum.df, paste0(interfileFolder, '/plotdatasum.df_Yeo', Yeoresolution,'_sumSCinvnode_CV', CVthr,'.rds'))
}else{plotdatasum.df <- readRDS(paste0(interfileFolder, '/plotdatasum.df_Yeo', Yeoresolution,'_sumSCinvnode_CV', CVthr,'.rds'))}

# To avoid the influence of averaged weight on derivative analyses, we divided SC strength of each edges by their
# weight at age of 8.
SCdata.diw <- SCdata.sum.merge
for (i in 1:elementnum){
  SClabel <- names(SCdata.sum.merge)[i+1]
  plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label==SClabel, ]
  SCdata.diw[ ,SClabel] <- SCdata.sum.merge[ ,SClabel] / plotdata.tmp$fit[1]
}
saveRDS(SCdata.diw, paste0(interfileFolder, "/SCdata.diw_Yeo", Yeoresolution,"CV", CVthr, ".rds"))
if (str_detect(wdpath, "cuizaixu_lab")){
  covariates<-"sex+mean_fd"
  dataname<-"SCdata.diw"
  smooth_var<-"age"
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-names(SCdata.sum.merge)[1+x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 15)
  
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
}else{resultsum <- readRDS(paste0(interfileFolder, '/gammodel_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))}

# gam results
covariates<-"sex+mean_fd"
dataname<-"SCdata.diw"
smooth_var<-"age"
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
resultsum <- parLapply(cl, 1:elementnum, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
  gamresult<-as.data.frame(gamresult)
  return(gamresult)
})
gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df[,c(2:18)]<-lapply(gamresultsum.df[,c(2:18)], as.numeric)
summary(gamresultsum.df)
# add fdr p values
gamresultsum.df$pfdr<- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
gamresultsum.df$sig<-(gamresultsum.df$pfdr<0.05)
saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))


