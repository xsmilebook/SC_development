## This script is to generate a dataframe, in which each column is the strength for an edge in large-scale network.
## For schaefer 400 --> ds.resolution atlas, elementnum mean edges in large-scale network.
library(R.matlab)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(reshape)
rm(list = ls())
# set the resolution of large-scale SC matrix
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/qsiPrep/SC_matrix'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/demopath'
  Volume_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/schaefer400_7_nodevolume'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  }else{
  # in PC
  demopath<-"/Users/xuxiaoyu_work/Cuilab/open_dataset_information/ABCD/info"
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution)
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
    
  source(paste0(functionFolder, '/colorbarvalue.R'))
    
  }

Behavior <- read.csv(paste0(demopath, '/DemodfScreenFinal.csv'))
Behavior$siteID <- gsub("site0", "site", Behavior$siteID)
summary(Behavior[,c("sex", "eventname", "handness")])
factorvar<-c("sex", "eventname", "handness")
for (i in 1:3){
  Behavior[, factorvar[i]] <-as.factor(Behavior[, factorvar[i]])
}
Behavior$age <- Behavior$age / 12
#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
schaefer400_index_SA <- schaefer400_index_SA[order(schaefer400_index_SA$finalrank.wholebrainrank),]
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
orderSA_7.SA<-order(schaefer400_index$finalrank.wholebrainrank)
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.SAorder.delLM.rds'))

# delete limbic region and extract S-A order index
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
limbicindex.SA <- which(str_detect(schaefer400_index_SA$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
orderSA_7<-order(schaefer400_index$finalrank.wholebrainrank) #1~376
orderSA_7.delLM <- orderSA_7.SA[-limbicindex.SA] #1~400, but without limbic regions
schaefer400_index <- schaefer400_index[order(schaefer400_index$finalrank.wholebrainrank),]
schaefer400_index <- schaefer400_index %>%
  mutate(SAds.resolutionnode = ntile(finalrank.wholebrainrank, ds.resolution))
summary(schaefer400_index$SAds.resolutionnode)
SAds.resolutionnode <- schaefer400_index$SAds.resolutionnode

# SC 376*376 --> ds.resolution*ds.resolution label
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
SAds.resolution <- matrix_SCds.resolution[lower.tri(matrix_SCds.resolution, diag = T)]
saveRDS(SAds.resolution, paste0(interfileFolder, '/SA',ds.resolution, '_index.rds'))

#### import SC data
#### ds.resolution regions, (ds.resolution+1)*ds.resolution/2=elementnum SCs
#################################################
colname <- character(length = elementnum)
for (i in 1:elementnum){
  colname[i] <- paste0('SC.', as.character(i))
}
if (str_detect(wdpath, "cuizaixu_lab")){
  SCdata.sum <- mclapply(1:nrow(Behavior), function(i){
    scanID <- Behavior$scanID[i]
    siteID <- Behavior$siteID4[i]
    eventname <-strsplit(scanID, 'ses-')[[1]][2]
    SCname <- paste0(scanID, '_space-T1w_desc-preproc_msmtconnectome.mat')
    SC_file_path <- paste0(SC_path, '/', eventname, '/SIEMENS/', siteID, '/', SCname)
    volumefile <- paste0(Volume_path, '/', scanID, '.txt')
    if (file.exists(SC_file_path)){
      SCmat <- readMat(SC_file_path)
      SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
      SCmat_raw <- SCmat_raw[orderSA_7, orderSA_7]
      indexup <- upper.tri(SCmat_raw)
      indexsave <- !indexup
      SCmat_raw <- SCmat_raw[indexsave]
      SCmat_raw75 <- SCmat_raw
      SCmat_raw75[deleteindex75]<-0
      df <- data.frame(
        group = SAds.resolution,
        value75 = SCmat_raw75
      )
      
      result <- df %>%
        group_by(group) %>%
        summarise(sum_value75 = sum(value75))
      sumSC.raw75 <- result$sum_value75[1:elementnum]
      ## node volume
      if (file.exists(volumefile)){
        nodevolume <- read_table(volumefile, col_names=F)
        
        if (nrow(nodevolume)==453){
          nodevolume <- as.numeric(nodevolume$X1[orderSA_7.delLM])
          df2 <- data.frame(
            group = SAds.resolutionnode,
            value = nodevolume
          )
          result2 <- df2 %>%
            group_by(group) %>%
            summarise(sum_value = sum(value))
          nodevolume_sum <- result2$sum_value
          
          ### ds.resolution*ds.resolution
          volumeSC <- matrix(NA, ds.resolution, ds.resolution)
          for (x in 1:ds.resolution){
            for (y in 1:ds.resolution){
              volumeSC[x,y] <- (nodevolume_sum[x]+nodevolume_sum[y])/2
            }
          }
          volumeSC <- volumeSC[lower.tri(volumeSC, diag = T)]
          sumSC.invnode75 <- sumSC.raw75 / volumeSC
          
        }else{
          sumSC.invnode75 <- rep(NA, elementnum)
        }}else{
          sumSC.invnode75 <- rep(NA, elementnum)
        }
      
      ###keep lower triangle and diagonal
      SCdat75 <- as.data.frame(sumSC.invnode75)
      SCdat75 <- as.data.frame(t(SCdat75), row.names = NULL)
      names(SCdat75) <- colname
      row.names(SCdat75) <- NULL
      SCdat75$scanID[1] <- scanID
      #SCdata.sum75<-rbind(SCdata.sum75, SCdat75)
      
      # not inverse node volume
      SCdat75_noInvNode <- as.data.frame(sumSC.raw75)
      SCdat75_noInvNode <- as.data.frame(t(SCdat75_noInvNode), row.names = NULL)
      names(SCdat75_noInvNode) <- colname
      row.names(SCdat75_noInvNode) <- NULL
      SCdat75_noInvNode$scanID[1] <- scanID
      
    }
    print(i)
    return(data.frame(SCdat75, SCdat75_noInvNode))
  },  mc.cores = 50)
  SCdata.sum75 <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c(1:(elementnum+1))])))
  SCdata.sum75_noInvNode <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum+2):(elementnum*2+2))])))
  names(SCdata.sum75_noInvNode) <- c(colname, "scanID")
    
  SCdata.sum75.merge <- merge(SCdata.sum75, Behavior, by="scanID")
  SCdata.sum75_noInvNode.merge <- merge(SCdata.sum75_noInvNode, Behavior, by="scanID")
  
  # convert variables' types
  SCdata.sum75.merge$subID <- as.factor(SCdata.sum75.merge$subID) ; SCdata.sum75.merge$siteID <- as.factor(SCdata.sum75.merge$siteID)
  SCdata.sum75_noInvNode.merge$subID <- as.factor(SCdata.sum75_noInvNode.merge$subID) ; SCdata.sum75_noInvNode.merge$siteID <- as.factor(SCdata.sum75_noInvNode.merge$siteID)
  
  saveRDS(SCdata.sum75.merge, paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum75_noInvNode.merge, paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV75_sumSC.sum.msmtcsd.merge.rds'))
}else{
  SCdata.sum75.merge <- readRDS(paste0(interfileFolder, '/SCdata_SA',ds.resolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
}

# plot

Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))
SCdata.sum75.merge <- SCdata.sum75.merge[!is.na(SCdata.sum75.merge$SC.28),]
# 75
age.time <- c(9, 11, 13)
meanSCdata.sepage75 <- list()
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum75.merge$age>=age.tmp & SCdata.sum75.merge$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum75.merge[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:(elementnum+1))])
  meanSCdata.sepage75[[i]] <- meanSCdata.tmp}

SCmin <- min(c(meanSCdata.sepage75[[1]], meanSCdata.sepage75[[2]], meanSCdata.sepage75[[3]]))
SCmax <- max(c(meanSCdata.sepage75[[1]], meanSCdata.sepage75[[2]], meanSCdata.sepage75[[3]]))
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum75.merge$age>=age.tmp & SCdata.sum75.merge$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum75.merge[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:(elementnum+1))])
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- meanSCdata.tmp
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, ds.resolution)
  rownames(Matrix.tmp) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  Fig<-ggplot(data =matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
    scale_fill_distiller(type="seq", palette = "RdBu",limit=c(SCmin-0.1, SCmax+0.1), na.value = "grey")+
    scale_color_distiller(type="seq", palette = "RdBu",limit=c(SCmin-0.1, SCmax+0.1),na.value = "grey")+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
    ggtitle(label = paste("CV75, age", age.tmp, "~", (age.tmp+1)))+labs(x=NULL, y=NULL)+
    scale_y_continuous(breaks=NULL, labels = NULL)+
    scale_x_continuous(breaks=NULL, labels = NULL)+
    theme(axis.line = element_blank(), 
          #axis.ticks=element_line(linewidth = 0),
          axis.text.x=element_text(size=ds.resolution, angle=45, hjust=1), 
          axis.text.y=element_text(size=ds.resolution, angle=315, hjust=1,vjust=1),
          axis.title =element_text(size=18),
          plot.title = element_text(size=18, hjust = 0.5),
          legend.title=element_text(size=18),
          legend.text=element_text(size=18), 
          panel.background=element_rect(fill=NA),
          panel.grid.major=element_line(linewidth = 0), 
          panel.grid.minor=element_line(linewidth = 1))
  Fig
  filename<-paste0(FigureFolder,"/CV75/Matrix",ds.resolution, "_sumSCinvnode_Age8_22/Age_", age.tmp, "_ds.resolutionnet_delLM_CV75.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}




