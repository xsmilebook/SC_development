## ABCD data
## This script is to generate a dataframe, in which each column is the strength for an edge.
## For schaefer 400 atlas, 70786 edges left after deleting edges connecting to limbic regions.
library(R.matlab)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(corrplot)
rm(list = ls())
# Set path and load data
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/qsiPrep/SC_matrix'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath'
  Volume_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/schaefer400_7_nodevolume'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD'
}else{
  # in PC
  demopath<-"/Users/xuxiaoyu_work/Cuilab/open_dataset_information/ABCD/info"
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final')
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  
  source(paste0(functionFolder, '/colorbarvalue.R'))
  
}

Behavior <- read.csv(paste0(demopath, '/DemodfScreenFinal.csv'))
Behavior$siteID <- gsub("site0", "site", Behavior$siteID)
tab <- table(Behavior$siteID)
tab <- tab[order(tab, decreasing = T)]

# organise behavior data
##################
summary(Behavior[,c("sex", "eventname", "handness")])
factorvar<-c("sex", "eventname", "handness")
for (i in 1:3){
  Behavior[, factorvar[i]] <-as.factor(Behavior[, factorvar[i]])
}
table(Behavior$eventname)
# 2YearFollowUpYArm1  baselineYear1Arm1 
# 3155               3949
Behavior$age <- Behavior$age / 12
write.csv(Behavior, paste0(demopath, "/DemodfScreenFinal_SCdev.csv"), row.names = F)

# import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
## output matrix is in Yeo 7 order, so reorder schaefer400 index to Yeo 7 order
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
orderSA_7<-order(schaefer400_index$finalrank.wholebrainrank)

#### import SC data
#### 376 regions, 377*376/2=70876 SCs
#################################################
colname <- character(length = 70876)
for (i in 1:70876){
  colname[i] <- paste0('SC.', as.character(i))
}

SCdata.sum <- mclapply(1:nrow(Behavior), function(i){
  scanID <- Behavior$scanID[i]
  siteID <- Behavior$siteID4[i]
  eventname <-strsplit(scanID, 'ses-')[[1]][2]
  SCname <- paste0(scanID, '_space-T1w_desc-preproc_msmtconnectome.mat')
  SC_file_path <- paste0(SC_path, '/', eventname, '/SIEMENS/', siteID, '/', SCname)
  if (file.exists(SC_file_path)){
    SCmat <- readMat(SC_file_path)
    SCmat <- SCmat$schaefer400.sift.invnodevol.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    SCmat <- SCmat[orderSA_7, orderSA_7]
    indexup <- upper.tri(SCmat)
    indexsave <- !indexup ###keep lower triangle and diagonal
    SCdat <- as.data.frame(c(SCmat[indexsave]))
    SCdat <- as.data.frame(t(SCdat), row.names = NULL)
    names(SCdat) <- colname
    row.names(SCdat) <- NULL
    SCdat$scanID[1] <- scanID
  }
  return(SCdat)
}, mc.cores = 50)
ncoldf <- lapply(SCdata.sum, function(x) ncol(x))
SCdata.df <- do.call(rbind, SCdata.sum)
saveRDS(SCdata.df, paste0(interfileFolder, '/SCdata.sum.msmtcsd.delLM.rds'))
SCdata.sum.merge <- merge(SCdata.df, Behavior, by="scanID")
## calculate CV
meanSC<-colMeans(SCdata.sum.merge[,2:70877])
sd.SC <- mclapply(1:70876, function(x) {
  sd.tmp<-sd(SCdata.sum.merge[,x+1])
  return(sd.tmp)
}, mc.cores = 40)
sd.SC<-as.numeric(sd.SC)
CV.SC<-sd.SC/meanSC
Perct.CV.SC <- quantile(CV.SC, probs=seq(0, 1, 0.25)) # extract 
# 0%         25%        50%        75%       100% 
# 0.4295531  1.5129073  1.9121141  2.5241666 71.1706369
## ID of edges over threshold
deleteindex.delLM <- which(CV.SC>Perct.CV.SC[4])
SCdata.sum.merge[,deleteindex.delLM+1] <- 0
meanSC[deleteindex.delLM] <-0
saveRDS(SCdata.sum.merge, paste0(interfileFolder, '/SCdata.sum.CV75.merge.SAorder.delLM.rds'))
saveRDS(deleteindex.delLM, paste0(interfileFolder, '/CV75_deleteindex.SAorder.delLM.rds'))
########################################################################

## plot
meanSC<-colMeans(SCdata.sum.merge.delLM[,2:70877])
meanSC25<-colMeans(SCdata.sum.merge.delLM_25[,2:70877])
#376
Matsize<-376
Matrix.376 <- matrix(NA, nrow=Matsize, ncol =Matsize)
indexup <- upper.tri(Matrix.376)
indexsave <- !indexup ###keep lower triangle and diagonal
index <- as.numeric(meanSC)
Matrix.376[indexsave] <- index
Matrix.376[indexup] <- t(Matrix.376)[indexup]
colnames(Matrix.376) <-seq(1, Matsize)
rownames(Matrix.376) <-seq(1, Matsize)
tiff( 
  filename = paste0(FigureFolder, '/SCmatrix/SClog376_CV75.tiff'),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(log(Matrix.376), col=rev(COL2(diverging = "RdBu", n=200)), axes = FALSE)
dev.off()


# sparcity
sparcity <- rep(0, nrow(SCdata.sum.merge))
sparcity.df<-mclapply(1:nrow(SCdata.sum.merge), function(i){
  SCmat.tmp <- SCdata.sum.merge[i,2:70877]
  nover0 <- length(which(SCmat.tmp>0))
  sparcity<-nover0/70876
  return(sparcity)
}, mc.cores=4)
sparcity.df <- unlist(sparcity.df)
summary(sparcity.df); sd(sparcity.df)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08687 0.55373 0.59869 0.58916 0.63423 0.71385 
# sd=0.0617236
