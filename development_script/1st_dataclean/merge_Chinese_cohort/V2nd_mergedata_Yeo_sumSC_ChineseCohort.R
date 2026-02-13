## This script is to generate a dataframe, in which each column is the strength for an edge in large-scale network.
## For schaefer 400 --> Yeo network atlas, elementnum mean edges in large-scale network.
library(R.matlab)
library(ggplot2)
library(tidyverse)
library(parallel)
library(reshape)
library(openxlsx)
library(corrplot)
rm(list = ls())
parse_args <- function(args) {
  res <- list()
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    res[[kv[[1]]]] <- kv[[2]]
  }
  res
}

as_int <- function(x, default) {
  if (is.null(x) || is.na(suppressWarnings(as.integer(x)))) return(as.integer(default))
  as.integer(x)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}

# input Yeo resolution of 7 or 17.
Yeoresolution <- as_int(args$yeo_res, 17L)
if (Yeoresolution == 7) {
  Yeoresolution.delLM <- 6
} else if (Yeoresolution == 17) {
  Yeoresolution.delLM <- 15
} else {
  stop("Invalid Yeoresolution (use 7 or 17): ", Yeoresolution)
}
elementnum <- Yeoresolution.delLM * (Yeoresolution.delLM + 1) / 2
make_plots <- as_int(args$make_plots, 0L) == 1L

is_gpfs <- dir.exists("/ibmgpfs/cuizaixu_lab")
if (is_gpfs) {
  SC_path_Cui <- "/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/defaultatlas"
  SC_path_SNU <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SNU_data/processed/SC"
  SC_path_CCNP <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/SC"

  demopath_Cui <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath/CuiBP"
  demopath_SNU <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath/SNU"
  demopath_CCNP <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/demopath/CCNP"
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ChineseCohort"

  Volume_path_Cui <- "/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/schaefer400_nodevolume"
  Volume_path_SNU <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SNU_data/processed/schaefer400_nodevolume"
  Volume_path_CCNP <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/schaefer400_nodevolume"
} else {
  stop("This script is configured for the GPFS environment. Provide GPFS paths or run on the cluster.")
}

outputFolder <- if (!is.null(args$output_dir)) args$output_dir else file.path(project_root, "outputs", "intermediate", "1st_dataclean", "chinese_cohort")
FigureFolder <- if (!is.null(args$figure_dir)) args$figure_dir else file.path(project_root, "outputs", "figures", "1st_dataclean", "chinese_cohort")
dir.create(outputFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

# CuiBP
Behavior_Cui <- read.xlsx(file.path(project_root, "demopath", "basic_demo_merge_screen.xlsx")) # 152 subjects with complete dMRI & normal anat
Behavior_Cui$Sex <- as.factor(dplyr::recode(as.character(Behavior_Cui$Sex), "M" = "Male", "F" = "Female", .default = as.character(Behavior_Cui$Sex)))
Behavior_Cui <- dplyr::rename(Behavior_Cui, subID = MRI_ID)
Behavior_Cui$Age <- suppressWarnings(as.numeric(Behavior_Cui$Age))
Behavior_Cui <- Behavior_Cui[!is.na(Behavior_Cui$Age) & Behavior_Cui$Age <= 26, ]
# SNU
Behavior_SNU <- read.csv(paste0(demopath_SNU, "/basic_demo_merge_screen.csv")) # 145 subjects with complete dMRI & normal anat
Behavior_SNU <- Behavior_SNU %>% dplyr::rename(all_of(c(subID = "sub_ID", Sex = "Gender_num")))
Behavior_SNU$Sex <- factor(Behavior_SNU$Sex, levels = c(1, 2), labels = c("Male", "Female"))
Behavior_SNU$Age <- suppressWarnings(as.numeric(Behavior_SNU$Age))
Behavior_SNU <- Behavior_SNU[!is.na(Behavior_SNU$Age) & Behavior_SNU$Age <= 26, ]
# CCNP
Behavior_CCNP <- read.csv(paste0(demopath_CCNP, "/basic_demo_devCCNPPEK.csv")) # 323 subjects with complete dMRI & normal anat
Behavior_CCNP$Sex <- as.factor(Behavior_CCNP$Sex)

#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
## qsiprep output matrix is in Yeo 7 order, so reorder schaefer400 index to Yeo 7 order
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
## Rearrange left and right regions.
schaefer400_index$index_7network_LRmixed <- schaefer400_index$index
schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] <- schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] - 200
orderYeo_7<-order(schaefer400_index$index_7network_LRmixed)

schaefer400_index$index_17network_LRmixed <- schaefer400_index$index_17network
schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] <- schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] - 200
orderYeo_17<-order(schaefer400_index$index_17network_LRmixed)

# filter index of P75th and P25th of CV.
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))
deleteindex25 <- readRDS(paste0(interfileFolder, '/CV25_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))

# delete limbic region
# limbicindex: limbic region index in schaefer400_index
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))

# assign each region to Yeo.resolution*Yeo.resolution network.
schaefer400_index.Yeo7 <- schaefer400_index[order(schaefer400_index$index_7network_LRmixed),]
schaefer400_index.Yeo17 <- schaefer400_index[order(schaefer400_index$index_17network_LRmixed),]
schaefer400_index.Yeo7 <- schaefer400_index.Yeo7 %>% mutate(Yeo.resolutionnode = recode_factor(network_label,
                                                                                               "Vis" = 1,
                                                                                               "SomMot" = 2,
                                                                                               "DorsAttn" = 3,
                                                                                               "SalVentAttn" = 4,
                                                                                               "Cont" = 5,
                                                                                               "Default" = 6,
                                                                                               "Limbic" = 0))

summary(schaefer400_index.Yeo7$Yeo.resolutionnode)
schaefer400_index.Yeo17 <- schaefer400_index.Yeo17 %>% 
  mutate(Yeo.resolutionnode = recode_factor(network_label_17network,
                                            "VisCent" = 3,
                                            "VisPeri" = 1,
                                            "SomMotA" =2,
                                            "SomMotB" = 4,
                                            "DorsAttnA" =5,
                                            "DorsAttnB" = 7,
                                            "SalVentAttnA" =9,
                                            "SalVentAttnB" =13,
                                            "ContA" =11,
                                            "ContB" =15,
                                            "ContC" =6,
                                            "DefaultA" = 12,
                                            "DefaultB" = 14,
                                            "DefaultC" = 8,
                                            "TempPar" = 10))

summary(schaefer400_index.Yeo17$Yeo.resolutionnode)
if (Yeoresolution == 7){
  Yeo.resolutionnode <- schaefer400_index.Yeo7$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else if (Yeoresolution == 17){
  Yeo.resolutionnode <- schaefer400_index.Yeo17$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else{
  print("Invalid Yeoresolution!")
}
# SC 376*376 --> Yeoresolution.delLM*Yeoresolution.delLM
matrixYeo.resolution <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
matrixYeo.resolution[lower.tri(matrixYeo.resolution, diag = T)] <- c(1:elementnum)
matrixYeo.resolution[upper.tri(matrixYeo.resolution)] <- t(matrixYeo.resolution)[upper.tri(matrixYeo.resolution)]
matrix_SCYeo.resolution <- matrix(NA, 376, 376)
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    xindex <- which(Yeo.resolutionnode==x)
    yindex <- which(Yeo.resolutionnode==y)
    matrix_SCYeo.resolution[xindex, yindex] <- matrixYeo.resolution[x,y]
  }
}
# an index telling how 376*376 map to 12*12
Yeo.resolution.index <- matrix_SCYeo.resolution[lower.tri(matrix_SCYeo.resolution, diag = T)]
#################################################

#### import SC data
#### Yeoresolution.delLM regions, (Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum SCs
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

# rbind demographic data (align to reference merge script)
behavior_cols <- c("subID", "Age", "Sex", "Handedness", "ICV", "mean_fd")
Behavior <- Behavior_Cui %>% select(all_of(behavior_cols))
Behavior$study <- "CuiBP"
Behavior2 <- Behavior_SNU %>% select(all_of(behavior_cols))
Behavior2$study <- "SNU"
Behavior3 <- Behavior_CCNP %>%
  select(scanID, ScanAge, Sex, handedness, ICV, mean_fd) %>%
  dplyr::rename(all_of(c(subID = "scanID", Handedness = "handedness", Age = "ScanAge")))
Behavior3$Age <- suppressWarnings(as.numeric(Behavior3$Age))
Behavior3 <- Behavior3[!is.na(Behavior3$Age) & Behavior3$Age > 6 & Behavior3$Age <= 26, ]
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
  
  sc_file <- paste0(SC_path, '/', SCname)
  volumefile <- paste0(Volume_path, '/', subID, '_Volume7.txt')
  # all the T1 parcellation for HCPD succeed.
  if (file.exists(sc_file)){
    nodevolume <- read_table(volumefile, col_names = FALSE)
    SCmat <- readMat(sc_file) 
    # load steamline counts matrix & fiber length matrix
    SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    length_raw <- SCmat$schaefer400.radius2.meanlength.connectivity[schaefer376_delLM, schaefer376_delLM]
    if (Yeoresolution == 7){
      SCmat_raw <- SCmat_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
      length_raw <- length_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
    }else if (Yeoresolution == 17){
      SCmat_raw <- SCmat_raw[orderYeo_17, orderYeo_17]
      length_raw <- length_raw[orderYeo_17, orderYeo_17]
    }else{
      print("Invalid Yeoresolution!")
    }
    
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
      group = Yeo.resolution.index,
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
    mean_length75 <- (result$sum_length75 / result$sum_value75)[1:elementnum]
    mean_length25 <- (result$sum_length25 / result$sum_value25)[1:elementnum]
    sumSC.raw75 <- result$sum_value75[1:elementnum]
    sumSC.raw25 <- result$sum_value25[1:elementnum]
    ## node volume
    nodevolume <- as.numeric(nodevolume$X1)
    nodevolume <- nodevolume[schaefer376_delLM] # delete limbic regions
    if (Yeoresolution == 7){
      nodevolume <- nodevolume[orderYeo_7] # sorted by Yeo index
    }else if (Yeoresolution == 17){
      nodevolume <- nodevolume[orderYeo_17]
    }else{
      print("Invalid Yeoresolution!")
    }
    df2 <- data.frame(
      group = Yeo.resolutionnode,
      value = nodevolume
    )
    result2 <- df2 %>%
      group_by(group) %>%
      summarise(sum_value = sum(value))
    nodevolume_sum <- result2$sum_value[1:Yeoresolution.delLM] # sum of nodes' volume for each node fraction (Yeo.resolution).
    
    ### Yeo.resolution*Yeo.resolution
    volumeSC <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
    for (x in 1:Yeoresolution.delLM){
      for (y in 1:Yeoresolution.delLM){
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

saveRDS(SCdata.sum75.merge, file.path(outputFolder, paste0("SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.merge.rds")))
saveRDS(SCdata.sum25.merge, file.path(outputFolder, paste0("SCdata_Yeo", Yeoresolution, "_CV25_sumSCinvnode.sum.msmtcsd.merge.rds")))
saveRDS(SCdata.sum75_noInvNode.merge, file.path(outputFolder, paste0("SCdata_Yeo", Yeoresolution, "_CV75_sumSC.sum.msmtcsd.merge.rds")))
saveRDS(SCdata.sum25_noInvNode.merge, file.path(outputFolder, paste0("SCdata_Yeo", Yeoresolution, "_CV25_sumSC.sum.msmtcsd.merge.rds")))

if (!make_plots) quit(save = "no", status = 0)

### plot matrix
#Yeo.resolution
SCdata.sum75.merge <- readRDS(file.path(outputFolder, paste0("SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.merge.rds")))
SCdata.sum25.merge <- readRDS(file.path(outputFolder, paste0("SCdata_Yeo", Yeoresolution, "_CV25_sumSCinvnode.sum.msmtcsd.merge.rds")))
SCdata.sum75.merge <- SCdata.sum75.merge %>% drop_na("SC.1")
SCdata.sum25.merge <- SCdata.sum25.merge %>% drop_na("SC.1")

meanSC.75<-colMeans(SCdata.sum75.merge[,2:(elementnum+1)])
meanSC.25<-colMeans(SCdata.sum25.merge[,2:(elementnum+1)])
Matsize<-Yeoresolution.delLM
Matrix.Yeoresolution.delLM <- matrix(NA, nrow=Matsize, ncol =Matsize)
indexup <- upper.tri(Matrix.Yeoresolution.delLM)
indexsave <- !indexup ###keep lower triangle and diagonal
index <- as.numeric(meanSC.75)
Matrix.Yeoresolution.delLM[indexsave] <- index
Matrix.Yeoresolution.delLM[indexup] <- t(Matrix.Yeoresolution.delLM)[indexup]
colnames(Matrix.Yeoresolution.delLM) <-seq(1, Matsize)
rownames(Matrix.Yeoresolution.delLM) <-seq(1, Matsize)

dir.create(file.path(FigureFolder, "SCmatrix"), showWarnings = FALSE, recursive = TRUE)
tiff( 
  filename = file.path(FigureFolder, "SCmatrix", paste0("SC_Yeo", Yeoresolution, "_CV75_sumSCinvnode.tiff")),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(Matrix.Yeoresolution.delLM, col=rev(COL2(diverging = "RdBu", n=200)), axes = TRUE)
dev.off()

index <- as.numeric(meanSC.25)
Matrix.Yeoresolution.delLM[indexsave] <- index
Matrix.Yeoresolution.delLM[indexup] <- t(Matrix.Yeoresolution.delLM)[indexup]
colnames(Matrix.Yeoresolution.delLM) <-seq(1, Matsize)
rownames(Matrix.Yeoresolution.delLM) <-seq(1, Matsize)

tiff( 
  filename = file.path(FigureFolder, "SCmatrix", paste0("SC_Yeo", Yeoresolution, "_CV25_sumSCinvnode.tiff")),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(Matrix.Yeoresolution.delLM, col=rev(COL2(diverging = "RdBu", n=200)), axes = TRUE)
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


