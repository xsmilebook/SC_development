## This script is to generate a dataframe, in which each column is the strength for an edge in large-scale network.
## For schaefer 400 --> ds.resolution atlas, elementnum mean edges in large-scale network.
library(R.matlab)
library(ggplot2)
library(tidyverse)
library(parallel)
library(reshape)
library(openxlsx)
library(corrplot)
rm(list = ls())
wdpath <- getwd()
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/connectomeMatrix/defaultatlas'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/demopath'
  Volume_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/schaefer400_nodevolume'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_HCPD'
}else{
  # in PC
  demopath<-'/Users/xuxiaoyu_work/Cuilab/open_dataset_information/HCP/HCPD_behavior'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
  FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final'
  functionFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  source(paste0(functionFolder, '/colorbarvalue.R'))
}

Behavior <- read.csv(paste0(demopath, '/HCPD_demo_behav.csv'))
factorvar<-c("sex", "race_ethnicity", "handnessfactor")
for (i in 1:3){
  Behavior[, factorvar[i]] <-as.factor(Behavior[, factorvar[i]])
}
Behavior <- Behavior[Behavior$age>8,]
# load data
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
# schaefer400_index_SA: order the schaefer400 rows based on SA axis
schaefer400_index_SA <- schaefer400_index_SA[order(schaefer400_index_SA$finalrank.wholebrainrank),]
# schaefer400_index: order the schaefer400 rows based on schaefer400-7 index
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
# orderSA_7.SA: the order index of SA axis in schaefer400_index
orderSA_7.SA<-order(schaefer400_index$finalrank.wholebrainrank)
# filter index of P75th and of CV.
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.SAorder.delLMover8.rds'))

# delete limbic region and extract S-A order index
# limbicindex: limbic region index in schaefer400_index
# limbicindex.SA: limbic region index in schaefer400_index_SA
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
limbicindex.SA <- which(str_detect(schaefer400_index_SA$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index # schaefer376_delLM: index without limbic regions
orderSA_7<-order(schaefer400_index$finalrank.wholebrainrank) # the order index of SA axis in schaefer400_index without limbic regions
orderSA_7.delLM <- orderSA_7.SA[-limbicindex.SA] # the order index of SA axis in schaefer400_index_SA without limbic regions
# assign each region to ds.resolution*ds.resolution network.
schaefer400_index <- schaefer400_index[order(schaefer400_index$finalrank.wholebrainrank),]
schaefer400_index <- schaefer400_index %>%
  mutate(SAds.resolutionnode = ntile(finalrank.wholebrainrank, ds.resolution))
summary(schaefer400_index$SAds.resolutionnode)
SAds.resolutionnode <- schaefer400_index$SAds.resolutionnode

# SC 376*376 --> ds.resolution*ds.resolution
matrixds.resolution <- matrix(NA, ds.resolution, ds.resolution)
matrixds.resolution[lower.tri(matrixds.resolution, diag = T)] <- c(1:elementnum)
matrixds.resolution[upper.tri(matrixds.resolution)] <- t(matrixds.resolution)[upper.tri(matrixds.resolution)]
matrix_SCds.resolution <- matrix(NA, 376, 376)
for (x in 1:ds.resolution){
  for (y in 1:ds.resolution){
    xindex <- which(SAds.resolutionnode==x)
    yindex <- which(SAds.resolutionnode==y)
    matrix_SCds.resolution[xindex, yindex] <- matrixds.resolution[x,y]
  }
}
# an index telling how 376*376 map to 12*12
SAds.resolution <- matrix_SCds.resolution[lower.tri(matrix_SCds.resolution, diag = T)]
#################################################

#### import SC data
#### ds.resolution regions, (ds.resolution+1)*ds.resolution/2=elementnum SCs
#### extract a dataframe containing elementnum columns, each represents an edge.
#################################################
colname <- character(length = elementnum)
for (i in 1:elementnum){
  colname[i] <- paste0('SC.', as.character(i))
}

SCdata.sum<- data.frame(t(rep(0,elementnum)))
names(SCdata.sum)<-colname
SCdata.sum$subID <- "NULL"
SCdata.sum75 <- SCdata.sum
SCdata.sum75_noInvNode <- SCdata.sum


for (i in 1:nrow(Behavior)){
  subID <- Behavior$subID[i]
  SCname <- paste0(subID, '_space-T1w_desc-preproc_msmtconnectome.mat')
  volumefile <- paste0(Volume_path, '/', subID, '_Volume7.txt')
  # all the T1 parcellation for HCPD succeed.
  if (file.exists(paste0(SC_path, '/', SCname))){
    SCmat <- readMat(paste0(SC_path, '/', SCname)) 
    # load steamline counts matrix & fiber length matrix
    SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    SCmat_raw <- SCmat_raw[orderSA_7, orderSA_7] # 376*376 nodes sorted by S-A axis
    indexup <- upper.tri(SCmat_raw)
    indexsave <- !indexup
    SCmat_raw <- SCmat_raw[indexsave] # 1*70876 each element represents streamline counts
    SCmat_raw75 <- SCmat_raw
    SCmat_raw75[deleteindex75]<-0 # remove top 1/4 inconsistent connetions
    
    df <- data.frame(
      group = SAds.resolution,
      value75 = SCmat_raw75,
    )
    # compute the sum of streamline counts / length for each fraction, in total of elementnum.
    result <- df %>%
      group_by(group) %>%
      summarise(sum_value75 = sum(value75))
    sumSC.raw75 <- result$sum_value75[1:elementnum]

    ## node volume
    nodevolume <- read_table(volumefile, col_names=F)
    nodevolume <- as.numeric(nodevolume$X1[orderSA_7.delLM]) # sort as SA-axis without limbic
    df2 <- data.frame(
      group = SAds.resolutionnode,
      value = nodevolume
    )
    result2 <- df2 %>%
      group_by(group) %>%
      summarise(sum_value = sum(value))
    nodevolume_sum <- result2$sum_value # sum of nodes' volume for each node fraction (ds.resolution).
    
    ### ds.resolution*ds.resolution
    volumeSC <- matrix(NA, ds.resolution, ds.resolution)
    for (x in 1:ds.resolution){
      for (y in 1:ds.resolution){
        volumeSC[x,y] <- (nodevolume_sum[x]+nodevolume_sum[y])/2 # geometric mean
      }
    }
    volumeSC <- volumeSC[lower.tri(volumeSC, diag = T)] # the scale values of node volume for each edge.
    sumSC.invnode75 <- sumSC.raw75 / volumeSC
    
    SCdat75 <- as.data.frame(sumSC.invnode75)
    SCdat75 <- as.data.frame(t(SCdat75), row.names = NULL)
    names(SCdat75) <- colname
    row.names(SCdat75) <- NULL
    SCdat75$subID[1] <- subID
    SCdata.sum75<-rbind(SCdata.sum75, SCdat75)
    
    SCdat75_noinvnode <- as.data.frame(sumSC.raw75)
    SCdat75_noinvnode <- as.data.frame(t(SCdat75_noinvnode), row.names = NULL)
    names(SCdat75_noinvnode) <- colname
    row.names(SCdat75_noinvnode) <- NULL
    SCdat75_noinvnode$subID[1] <- subID
    SCdata.sum75_noInvNode <- rbind(SCdata.sum75_noInvNode, SCdat75_noinvnode)
    
  }
}
SCdata.sum75<-SCdata.sum75[-1,]
SCdata.sum75_noInvNode <- SCdata.sum75_noInvNode[-1,]
SCdata.sum75.merge <- merge(SCdata.sum75, Behavior, by="subID")

SCdata.sum75_noInvNode.merge <- merge(SCdata.sum75_noInvNode, Behavior, by="subID")
# exclude subjects with big head motion
mean(SCdata.sum75.merge$mean_fd)+3*sd(SCdata.sum75.merge$mean_fd) #1.438502
SCdata.sum75.merge <- SCdata.sum75.merge[SCdata.sum75.merge$mean_fd<1.44, ] #14 subjects were excluded
SCdata.sum75_noInvNode.merge <- SCdata.sum75_noInvNode.merge[SCdata.sum75_noInvNode.merge$mean_fd<1.44, ] #14 subjects were excluded

saveRDS(SCdata.sum75.merge, paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
saveRDS(SCdata.sum75_noInvNode.merge, paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV75_sumSC.sum.msmtcsd.merge.rds'))

### plot matrix
#ds.resolution
meanSC.75<-colMeans(SCdata.sum75.merge[,2:(elementnum+1)])
Matsize<-ds.resolution
Matrix.ds.resolution <- matrix(NA, nrow=Matsize, ncol =Matsize)
indexup <- upper.tri(Matrix.ds.resolution)
indexsave <- !indexup ###keep lower triangle and diagonal
index <- as.numeric(meanSC.75)
Matrix.ds.resolution[indexsave] <- index
Matrix.ds.resolution[indexup] <- t(Matrix.ds.resolution)[indexup]
colnames(Matrix.ds.resolution) <-seq(1, Matsize)
rownames(Matrix.ds.resolution) <-seq(1, Matsize)

tiff( 
  filename = paste0(FigureFolder, '/SCmatrix/SC', ds.resolution, '_CV75_sumSCinvnode.tiff'),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(Matrix.ds.resolution, col=rev(COL2(diverging = "RdBu", n=200)), axes = TRUE)
dev.off()

# sparcity
sparcity <- rep(0, nrow(SCdata.sum75.merge))
sparcity.df<-mclapply(1:nrow(SCdata.sum75.merge), function(i){
  SCmat.tmp <- SCdata.sum75.merge[i,2:(elementnum+1)]
  nover0 <- length(which(SCmat.tmp>0))
  sparcity<-nover0/elementnum
  return(sparcity)
}, mc.cores=4)
sparcity.df <- unlist(sparcity.df)
summary(sparcity.df)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1       1       1       1       1       1 


## output
#1. SCdata_SA'ds.resolution'_CV75_sumSCinvnode.sum.msmtcsd.merge.rds : elementnum variables + behavior variables + subID * obs
# elementnum edges connect ds.resolution fractions of cortical regions. Edge weights are the connectivity scales by the node volumes.


