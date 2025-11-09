## This script is to generate a dataframe, in which each column is the strength for an edge.
## For schaefer 400 atlas, 70786 edges left after deleting edges connecting to limbic regions.
## The ID of edges which should be deleted according to consistency threshold will be saved out. 
library(R.matlab)
library(mgcv)
library(visreg)
library(ggplot2)
library(tidyverse)
library(parallel)
library(reshape)
library(corrplot)
rm(list = ls())
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/connectomeMatrix/defaultatlas'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath'
  Volume_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/schaefer400_nodevolume'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD'
}else{
  # in PC
  demopath<-'D:/xuxiaoyu/open_dataset_information/HCP/HCPD_behavior'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  source(paste0(functionFolder, '/colorbarvalue.R'))
}
Behavior <- read.csv(paste0(demopath, '/HCPD_demo_behav.csv')) # 622 subjects with complete dMRI & normal anat

factorvar<-c("sex", "race_ethnicity", "handnessfactor")
for (i in 1:3){
  Behavior[, factorvar[i]] <-as.factor(Behavior[, factorvar[i]])
}

#################################################
### exclude subjects with big head motion
Behavior <- Behavior[Behavior$age>8,] # 604 left
Behavior$subID <- paste0("sub-", Behavior$subID)
#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
## qsiprep output matrix is in Yeo 7 order, so reorder schaefer400 index to Yeo 7 order
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
orderSA_7<-order(schaefer400_index$finalrank.wholebrainrank)

#################################################
#### import SC data
#### 376 regions, 377*376/2=70876 SCs
#################################################
colname <- character(length = 70876)
for (i in 1:70876){
  colname[i] <- paste0('SC.', as.character(i))
}
SCdata.sum<- data.frame(t(rep(0,70876)))
names(SCdata.sum)<-colname
SCdata.sum$subID <- "NULL"
for (i in 1:nrow(Behavior)){
  subID <- Behavior$subID[i]
  SCname <- paste0(subID, '_space-T1w_desc-preproc_msmtconnectome.mat')
  volumefile <- paste0(Volume_path, '/', subID, '_Volume7.txt')
  if (file.exists(paste0(SC_path, '/', SCname))){
    SCmat <- readMat(paste0(SC_path, '/', SCname))
    ## Normalize the SC counts by volume geometric mean.
    SCmat <- SCmat$schaefer400.sift.radius2.count.connectivity[1:400, 1:400]
    nodevolume <- read_table(volumefile, col_names=F)
    nodevolume <- as.numeric(nodevolume$X2[1:400]) # Use volume rather than voxel size
    
    volumeSC <- matrix(NA, 400, 400)
    for (x in 1:400){
      for (y in 1:400){
        volumeSC[x,y] <- sqrt(nodevolume[x]*nodevolume[y])
      }
    }
    
    SCmat.invnode <- SCmat / volumeSC
    
    #Reorder the nodes
    SCmat.invnode <- SCmat.invnode[schaefer376_delLM, schaefer376_delLM]
    SCmat <- SCmat[orderSA_7, orderSA_7]
    indexup <- upper.tri(SCmat)
    indexsave <- !indexup ###keep lower triangle and diagonal
    SCdat <- as.data.frame(c(SCmat[indexsave]))
    SCdat <- as.data.frame(t(SCdat), row.names = NULL)
    names(SCdat) <- colname
    row.names(SCdat) <- NULL
    SCdat$subID[1] <- subID
    SCdata.sum<-rbind(SCdata.sum, SCdat)
  }
}
SCdata.sum<-SCdata.sum[-1,] # 1 hour required
saveRDS(SCdata.sum, paste0(interfileFolder, '/SCdata.sum.msmtcsd.delLMover8.rds'))#70876, S-A order, 604
SCdata.sum.merge <- merge(SCdata.sum, Behavior, by="subID")
mean(SCdata.sum.merge$mean_fd)+3*sd(SCdata.sum.merge$mean_fd) #1.44
SCdata.sum.merge <- SCdata.sum.merge[SCdata.sum.merge$mean_fd<1.44,] # excessive head motion 14 removed
saveRDS(SCdata.sum.merge, paste0(interfileFolder, '/SCdata.sum.msmtcsd.delLMover8.merge.rds')) #590

## calculate CV
meanSC<-colMeans(SCdata.sum.merge[,2:70877])
sd.SC <- mclapply(1:70876, function(x) {
  sd.tmp<-sd(SCdata.sum.merge[,x+1])
  return(sd.tmp)
}, mc.cores = 4)
sd.SC<-as.numeric(sd.SC)
CV.SC<-sd.SC/meanSC
Perct.CV.SC <- quantile(CV.SC, probs=seq(0, 1, 0.25))

## ID of edges over threshold
deleteindex.delLM <- which(CV.SC>Perct.CV.SC[4])
SCdata.sum.merge[,deleteindex.delLM+1] <- 0
meanSC[deleteindex.delLM] <-0
saveRDS(SCdata.sum.merge, paste0(interfileFolder, '/SCdata.sum.CV75.merge.SAorder.delLMover8.rds'))
saveRDS(deleteindex.delLM, paste0(interfileFolder, '/CV75_deleteindex.SAorder.delLMover8.rds'))
########################################################################

### plot matrix
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
#  Min.   1st Qu.  Median  Mean   3rd Qu.  Max. 
# 0.4228  0.6651  0.6810  0.6754  0.6929  0.7188 
# sd = 0.02742583
