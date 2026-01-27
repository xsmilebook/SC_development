## This script is to fit gam models for each edge.
## Statistical indexes and gam model files will be generated.  
rm(list=ls())
library(mgcv)
library(parallel)
library(tidyverse)
wdpath <- getwd()
# set resolution
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
# set path
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/HBN/processed/SC'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath/HBN'
  Volume_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/HBN/processed/schaefer400_nodevolume'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HBN'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Figure_HBN_final'
}else{
  # in PC
  demopath<-'D:/xuxiaoyu/open_dataset_information/HBN/demography'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HBN'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HBN_final'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
}

# set which consistency threshold used to filter spurious streamlines
CVthr=75
# load data
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_SA", ds.resolution,"_CV", CVthr,"_sumSCinvnode.TractSeg.combatage.rds"))
nrow(SCdata.sum.merge)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
# source function
source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))
detectCores()
## calculate gam results
covariates<-"sex+mean_fd"
dataname<-"SCdata.sum.merge"
smooth_var<-"age"

resultsum <- list()

for (x in 1:elementnum){
  SClabel<-grep("SC.", names(SCdata.sum.merge), value=T)[x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
  gamresult<-as.data.frame(gamresult)
  
  resultsum[[x]] <- gamresult
}

gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df[,c(2:18)]<-lapply(gamresultsum.df[,c(2:18)], as.numeric)
class(gamresultsum.df$gam.smooth.pvalue)
class(gamresultsum.df$partialRsq)
# add fdr p values
gamresultsum.df$pfdr<- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
gamresultsum.df$sig<-(gamresultsum.df$pfdr<0.05)
saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_TractSeg_CV', CVthr,'.rds'))

  
## calculate gam models
resultsum <- mclapply(1:elementnum, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
  return(gamresult)
}, mc.cores = 30)
saveRDS(resultsum, paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_TractSeg_CV', CVthr,'.rds'))

## plot data raw
#### generate function data
gammodelsum <- resultsum
# clusterExport(cl, varlist = ls(), envir = .GlobalEnv)

plotdatasum<-mclapply(1:elementnum, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, dataname=NA, "age")
  plotdata$SC_label <- names(plotdata)[14]
  plotdata[,14] <- NULL
  return(plotdata)
}, mc.cores = 30)
plotdatasum.df <- do.call(rbind, lapply(plotdatasum, function(x) data.frame(x)))
saveRDS(plotdatasum.df, paste0(interfileFolder, '/plotdatasum.df_SA', ds.resolution,'_sumSCinvnode_TractSeg_CV', CVthr,'.rds'))

# To avoid the influence of averaged weight on derivative analyses, we divided SC strength of each edges by their
# weight at age of 6.
SCdata.diw <- SCdata.sum.merge
for (i in 1:elementnum){
  SClabel <- names(SCdata.sum.merge)[i+1]
  plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label==SClabel, ]
  SCdata.diw[ ,SClabel] <- SCdata.sum.merge[ ,SClabel] / plotdata.tmp$fit[1]
}
saveRDS(SCdata.diw, paste0(interfileFolder, "/SCdata.diw_SA", ds.resolution,"_TractSeg_CV", CVthr, ".rds"))


covariates<-"sex+mean_fd"
dataname<-"SCdata.diw"
smooth_var<-"age"
# clusterExport(cl, varlist = ls(), envir = .GlobalEnv)

resultsum <- mclapply(1:elementnum, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
  return(gamresult)
}, mc.cores = 10)

saveRDS(resultsum, paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_TractSeg_CV', CVthr,'_scale_TRUE.rds'))

# gam results
covariates<-"sex+mean_fd"
dataname<-"SCdata.diw"
smooth_var<-"age"

resultsum <- mclapply(1:elementnum, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
  gamresult<-as.data.frame(gamresult)
  return(gamresult)
}, mc.cores = 30)
gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df[,c(2:18)]<-lapply(gamresultsum.df[,c(2:18)], as.numeric)
summary(gamresultsum.df)
saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_TractSeg_CV', CVthr,'_scale_TRUE.rds'))


