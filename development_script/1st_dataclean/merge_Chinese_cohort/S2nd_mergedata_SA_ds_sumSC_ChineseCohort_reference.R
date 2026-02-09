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

SC_path_Cui <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/defaultatlas'
SC_path_SNU <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SNU_data/processed/SC'
SC_path_CCNP <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/SC'

demopath_Cui<-'/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/demography'
demopath_SNU<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath/SNU'
demopath_CCNP<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath/CCNP'
interfileFolder <- '/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort'
FigureFolder <- '/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/Figure_ChineseCohort_final'

Volume_path_Cui <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/schaefer400_nodevolume'
Volume_path_SNU <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SNU_data/processed/schaefer400_nodevolume'
Volume_path_CCNP <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/schaefer400_nodevolume'


# CuiBP
Behavior_Cui <- read.xlsx(paste0(demopath_Cui, '/basic_demo_merge_screen.xlsx'), detectDates=T) # 491 subjects with complete dMRI & normal anat
Behavior_Cui$Sex <- as.factor(dplyr::recode(as.character(Behavior_Cui$Sex), "M" = "Male", "F" = "Female"))
Behavior_Cui <- dplyr::rename(Behavior_Cui, subID = MRI_ID)
# SNU
Behavior_SNU <- read.csv(paste0(demopath_SNU, '/basic_demo_merge_screen.csv')) # 145 subjects with complete dMRI & normal anat
Behavior_SNU <- Behavior_SNU %>% dplyr::rename(all_of(c(subID="sub_ID", Sex="Gender_num")))
Behavior_SNU$Sex <- factor(Behavior_SNU$Sex, levels = c(1, 2), labels = c("Male", "Female"))
# CCNP
Behavior_CCNP <- read.csv(paste0(demopath_CCNP, '/basic_demo_devCCNPPEK.csv')) # 323 subjects with complete dMRI & normal anat
Behavior_CCNP$Sex <- as.factor(Behavior_CCNP$Sex)

# load data
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
# schaefer400_index_SA: order the schaefer400 rows based on SA axis
schaefer400_index_SA <- schaefer400_index_SA[order(schaefer400_index_SA$finalrank.wholebrainrank),]
# schaefer400_index: order the schaefer400 rows based on schaefer400-7 index
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
# orderSA_7.SA: the order index of SA axis in schaefer400_index
orderSA_7.SA<-order(schaefer400_index$finalrank.wholebrainrank)
# filter index of P75th and P25th of CV.
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.SAorder.delLM.rds'))
deleteindex25 <- readRDS(paste0(interfileFolder, '/CV25_deleteindex.SAorder.delLM.rds'))

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
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   6.000   6.457   9.000  12.000 
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
SCdata.sum75 <- SCdata.sum25 <- SCdata.sum
SCdata.sum75_noInvNode <- SCdata.sum25_noInvNode <- SCdata.sum

# length
colname2 <- character(length = elementnum)
for (i in 1:elementnum){
  colname2[i] <- paste0('length.', as.character(i))
}

SClength.sum<- data.frame(t(rep(0,elementnum)))
names(SClength.sum)<-colname2
SClength.sum$subID <- "NULL"
SClength.sum75 <- SClength.sum25 <- SClength.sum

# rbind demographic data
Behavior <- Behavior_Cui %>% select(subID, Age, Sex, Handedness, ICV, mean_fd)
Behavior$study <- "CuiBP"
Behavior2 <- Behavior_SNU %>% select(subID, Age, Sex, Handedness, ICV, mean_fd)
Behavior2$study <- "SNU"
Behavior3 <- Behavior_CCNP %>% select(scanID, ScanAge, Sex, handedness, ICV, mean_fd) %>% 
  dplyr::rename(all_of(c(subID="scanID", Handedness="handedness",  Age="ScanAge")))
Behavior3$study <- "CCNP"

Behavior <- rbind(Behavior, Behavior2, Behavior3)

for (i in 1:nrow(Behavior)){
  subID <- Behavior$subID[i]
  if (i <= nrow(Behavior_Cui)){
    SCname <- paste0(subID, '_dir-PA_space-T1w_desc-preproc_msmtconnectome.mat')
    SC_path <- SC_path_Cui
    Volume_path <- Volume_path_Cui
    
  }else if (i <= nrow(Behavior_Cui)+nrow(Behavior_SNU)){
    SCname <- paste0(subID, '_dir-AP_space-T1w_desc-preproc_dhollanderconnectome.mat')
    SC_path <- SC_path_SNU
    Volume_path <- Volume_path_SNU
    
  }else if (i <= nrow(Behavior_Cui)+nrow(Behavior_SNU)+nrow(Behavior_CCNP)){
    SCname <- paste0(subID, '_space-T1w_desc-preproc_dhollanderconnectome.mat')
    SC_path <- SC_path_CCNP
    Volume_path <- Volume_path_CCNP
    
  }

  volumefile <- paste0(Volume_path, '/', subID, '_Volume7.txt')
  # all the T1 parcellation for CuiBP succeed.
  if (file.exists(paste0(SC_path, '/', SCname))){
    SCmat <- readMat(paste0(SC_path, '/', SCname)) 
    # load steamline counts matrix & fiber length matrix
    SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    SCmat_raw <- SCmat_raw[orderSA_7, orderSA_7] # 376*376 nodes sorted by S-A axis
    length_raw <- SCmat$schaefer400.radius2.meanlength.connectivity[schaefer376_delLM, schaefer376_delLM]
    length_raw <- length_raw[orderSA_7, orderSA_7] # 376*376 nodes sorted by S-A axis
    totallength_raw <- length_raw * SCmat_raw
    indexup <- upper.tri(SCmat_raw)
    indexsave <- !indexup
    SCmat_raw <- SCmat_raw[indexsave] # 1*70876 each element represents streamline counts
    SCmat_raw75 <- SCmat_raw25 <- SCmat_raw
    SCmat_raw75[deleteindex75]<-0 # remove top 1/4 inconsistent connetions
    SCmat_raw25[deleteindex25]<-0 # remove top 3/4 inconsistent connetions
    totallength_raw <- totallength_raw[indexsave]
    totallength_raw75 <- totallength_raw25 <- totallength_raw
    totallength_raw75[deleteindex75]<-0
    totallength_raw25[deleteindex25]<-0
    df <- data.frame(
      group = SAds.resolution,
      value75 = SCmat_raw75,
      value25 = SCmat_raw25,
      length75 = totallength_raw75,
      length25 = totallength_raw25
    )
    # compute the sum of streamline counts / length for each fraction, in total of elementnum.
    result <- df %>%
      group_by(group) %>%
      summarise(sum_value75 = sum(value75), sum_value25 = sum(value25), sum_length75=sum(length75), 
                sum_length25=sum(length25))
    mean_length75 <- result$sum_length75 / result$sum_value75
    mean_length25 <- result$sum_length25 / result$sum_value25
    sumSC.raw75 <- result$sum_value75[1:elementnum]
    sumSC.raw25 <- result$sum_value25[1:elementnum]
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
        volumeSC[x,y] <- (nodevolume_sum[x]+nodevolume_sum[y])/2
      }
    }
    volumeSC <- volumeSC[lower.tri(volumeSC, diag = T)] # the scale values of node volume for each edge.
    sumSC.invnode75 <- sumSC.raw75 / volumeSC
    sumSC.invnode25 <- sumSC.raw25 / volumeSC
    
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
    
    SCdat25 <- as.data.frame(sumSC.invnode25)
    SCdat25 <- as.data.frame(t(SCdat25), row.names = NULL)
    names(SCdat25) <- colname
    row.names(SCdat25) <- NULL
    SCdat25$subID[1] <- subID
    SCdata.sum25<-rbind(SCdata.sum25, SCdat25)
    
    SCdat25_noinvnode <- as.data.frame(sumSC.raw25)
    SCdat25_noinvnode <- as.data.frame(t(SCdat25_noinvnode), row.names = NULL)
    names(SCdat25_noinvnode) <- colname
    row.names(SCdat25_noinvnode) <- NULL
    SCdat25_noinvnode$subID[1] <- subID
    SCdata.sum25_noInvNode <- rbind(SCdata.sum25_noInvNode, SCdat25_noinvnode)
    
    mean_length75 <- as.data.frame(t(mean_length75))
    names(mean_length75) <- colname2
    mean_length75$subID[1] <- subID
    SClength.sum75 <- rbind(SClength.sum75, mean_length75)
    
    mean_length25 <- as.data.frame(t(mean_length25))
    names(mean_length25) <- colname2
    mean_length25$subID[1] <- subID
    SClength.sum25 <- rbind(SClength.sum25, mean_length25)
  }
}
SCdata.sum75<-SCdata.sum75[-1,]
SCdata.sum25<-SCdata.sum25[-1,]
SCdata.sum75_noInvNode <- SCdata.sum75_noInvNode[-1,]
SCdata.sum25_noInvNode <- SCdata.sum25_noInvNode[-1,]
SClength.sum75 <- SClength.sum75[-1,]
SClength.sum25 <- SClength.sum25[-1,]
SCdata.sum75.merge <- merge(SCdata.sum75, Behavior, by="subID")
SCdata.sum25.merge <- merge(SCdata.sum25, Behavior, by="subID")
SCdata.sum75.merge <- merge(SCdata.sum75.merge, SClength.sum75, by="subID")
SCdata.sum25.merge <- merge(SCdata.sum25.merge, SClength.sum25, by="subID")

SCdata.sum75_noInvNode.merge <- merge(SCdata.sum75_noInvNode, Behavior, by="subID")
SCdata.sum25_noInvNode.merge <- merge(SCdata.sum25_noInvNode, Behavior, by="subID")

saveRDS(SCdata.sum75.merge, paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
saveRDS(SCdata.sum25.merge, paste0(interfileFolder, '/SCdata_SA', ds.resolution,'_CV25_sumSCinvnode.sum.msmtcsd.merge.rds'))
saveRDS(SCdata.sum75_noInvNode.merge, paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV75_sumSC.sum.msmtcsd.merge.rds'))
saveRDS(SCdata.sum25_noInvNode.merge, paste0(interfileFolder, '/SCdata_SA', ds.resolution,'_CV25_sumSC.sum.msmtcsd.merge.rds'))

### plot matrix
#ds.resolution
meanSC.75<-colMeans(SCdata.sum75.merge[,2:(elementnum+1)])
meanSC.25<-colMeans(SCdata.sum25.merge[,2:(elementnum+1)])
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
sparcity <- rep(0, nrow(SCdata.sum25.merge))
sparcity.df<-mclapply(1:nrow(SCdata.sum25.merge), function(i){
  SCmat.tmp <- SCdata.sum25.merge[i,2:(elementnum+1)]
  nover0 <- length(which(SCmat.tmp>0))
  sparcity<-nover0/elementnum
  return(sparcity)
}, mc.cores=4)
sparcity.df <- unlist(sparcity.df)
summary(sparcity.df)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9872       1       1       1       1       1 




