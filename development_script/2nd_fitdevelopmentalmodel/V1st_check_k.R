# This script is for check the k value.
rm(list=ls())
library(mgcv)
library(parallel)
interfileFolder_HCP <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
interfileFolder_ABCD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
interfileFolder_ChineseCohort <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ChineseCohort'
functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolderABCD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ABCD_final/SA12/CV75'
FigureFolderHCPD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12/CV75'
FigureFolderChineseCohort<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ChineseCohort_final/SA12/CV75'
resultFolderABCD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ABCD'
resultFolderHCPD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HCPD'
resultFolder_ChineseCohort <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ChineseCohort'

source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/gammsmooth.R'))
SCdata.hcp <- readRDS(paste0(interfileFolder_HCP, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
SCdata.ABCD <- readRDS(paste0(interfileFolder_ABCD, '/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds'))
SCdata.ChineseCohort <- readRDS(paste0(interfileFolder_ChineseCohort, '/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds'))

# load models
#gammodel.HCP.k3 <- readRDS(paste0(interfileFolder_HCP, '/gammodel78_sumSCinvnode_over8_CV75.rds'))
#gammodel.ABCD.k3 <- readRDS(paste0(interfileFolder_ABCD, '/gammodel78_sumSCinvnode_over8_siteall_CV75.rds'))

# HCPD k=3~6
####################################
AIC_k3_6 <- data.frame(region=paste0("SC.", 1:78), AIC.k3=rep(NA,78), AIC.k4=rep(NA,78),
                       AIC.k5=rep(NA,78),AIC.k6=rep(NA,78))
SCdata.hcp$sex <- as.factor(SCdata.hcp$sex)
covariates<-"sex+mean_fd"
dataname<-"SCdata.hcp"
smooth_var<-"age"
for (k in c(4,5,6)){
  for (i in 1:78){
    model.k3.tmp <- gammodel.HCP.k3[[i]]
    AIC_k3_6$region[i] <- as.character(model.k3.tmp$terms[[2]])
    region <- AIC_k3_6$region[i]
    model.kN.tmp<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=k, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    df.tmp <- AIC(model.k3.tmp, model.kN.tmp)
    AIC_k3_6$AIC.k3[i] <- df.tmp$AIC[1]
    AIC_k3_6[i, paste0('AIC.k', k)] <- df.tmp$AIC[2]
  }
}
min.index <- apply(AIC_k3_6, 1, function(x) which.min(x))
AIC_k3_6$min.index <- min.index+1
table(AIC_k3_6$min.index) # 52/78 edges have the minimal AIC when k=3.

## bootstrap
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(mgcv)
  library(tidyr)
  
})

SCdata.hcp$sex <- as.factor(SCdata.hcp$sex)
covariates<-"sex+mean_fd"
smooth_var<-"age"
stratify.var <- SCdata.hcp$site

choice.df <- data.frame(boottime=1:1000)
for (i in 1:78){
  region <- paste0("SC.", i, "_h")
  clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
  
  bootchoice<-parLapply(cl, 1:1000, function(x){
    
    # bootstrap
    set.seed( seed=925 + x )
    SPLIT <- split(1:NROW(SCdata.hcp), stratify.var)
    LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X),replace=TRUE)})
    INDEX <- unsplit(LAPPLY, stratify.var)
    
    SCdata.hcp.subset <- SCdata.hcp[INDEX, ]
    
    # fit models
    AIC_df <- data.frame(knots=rep(NA,4), AIC=rep(NA,4))
    for (knots in c(3:6)){
      modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, "TRUE", covariates))
      gam.model <- gam(modelformula, method = "REML", data = SCdata.hcp.subset)
      AIC_df$knots[knots-2] <- knots
      AIC_df$AIC[knots-2] <- AIC(gam.model)
    }
    
    bestK <- AIC_df$knots[which.min(AIC_df$AIC)]
    
    return(bestK)
    
  })
  
  bootchoice.df <- unlist(bootchoice)
  choice.df[[region]] <- bootchoice.df
}

choice.df$MostFrequent <- apply(choice.df, 1, function(row) {
  freq <- table(row)
  as.numeric(names(freq)[which.max(freq)])
})

table(choice.df$MostFrequent)
# 3   4   5   6 
# 777  72  75  76

choice.df$modelnum_k3 <- rowSums(choice.df[,2:79] == 3)
write.csv(choice.df, paste0(resultFolderHCPD, "/select_k_gam.csv"), row.names = F)
## Best k
ggplot(data = choice.df, aes(MostFrequent, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", position = position_dodge(width = 1), breaks=c(2.5,3.5,4.5,5.5, 6.5), center=0) +
  scale_x_continuous(breaks=c(3,4,5,6), limits = c(2.5, 6.5))+
  labs(x = expression("The optimal "*italic("k")*" value"), y = "Frequency", title = "HCPD") +
  labs(y= "Frequency")+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderHCPD, '/GAM_k_choose.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderHCPD, '/GAM_k_choose.svg', sep = ''), width = 12, height = 10, units = "cm")

## Number of models select k=3
ggplot(data = choice.df, aes(modelnum_k3, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", position = position_dodge(width = 1)) +
  #geom_vline(xintercept = median(choice.df$modelnum_k3), color="red", linetype="dashed", linewidth=1)+
  labs(x = expression("Number of models choosing "*italic("k")*"=3"), y = "Frequency", title = "HCP-D") +
  labs(y= "Frequency")+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderHCPD, '/GAM_k_choose_modelnum.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderHCPD, '/GAM_k_choose_modelnum.svg', sep = ''), width = 12, height = 10, units = "cm")
######################################################


# ABCD k=3~6
####################################################
AIC_k3_6.abcd <- data.frame(region=paste0("SC.", 1:78), AIC.k3=rep(NA,78), AIC.k4=rep(NA,78),
                       AIC.k5=rep(NA,78),AIC.k6=rep(NA,78))
SCdata.ABCD$sex <- as.factor(SCdata.ABCD$sex)
covariates<-"sex+mean_fd"
dataname<-"SCdata.ABCD"
smooth_var<-"age"
for (k in c(4,5,6)){
  for (i in 1:78){
    model.k3.tmp <- gammodel.ABCD.k3[[i]]
    AIC_k3_6.abcd$region[i] <- as.character(model.k3.tmp$gam$terms[[2]])
    region <- AIC_k3_6.abcd$region[i]
    model.kN.tmp<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=k, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    df.tmp <- AIC(model.k3.tmp$mer, model.kN.tmp$mer)
    AIC_k3_6.abcd$AIC.k3[i] <- df.tmp$AIC[1]
    AIC_k3_6.abcd[i, paste0('AIC.k', k)] <- df.tmp$AIC[2]
  }
}

min.index <- apply(AIC_k3_6.abcd, 1, function(x) which.min(x))
AIC_k3_6.abcd$min.index <- min.index+1
table(AIC_k3_6.abcd$min.index) # 71/78 edges have the minimal AIC when k=3.

## bootstrap
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(mgcv)
  library(tidyr)
  library(gamm4)
})

SCdata.ABCD$sex <- as.factor(SCdata.ABCD$sex)
covariates<-"sex+mean_fd"
smooth_var<-"age"
stratify.var <- SCdata.ABCD$siteID

choice.df <- data.frame(boottime=1:1000)
for (i in 1:78){
  region <- paste0("SC.", i, "_h")
  clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
  
  bootchoice<-parLapply(cl, 1:1000, function(x){
    
    # bootstrap
    set.seed( seed=925 + x )
    SPLIT <- split(1:NROW(SCdata.ABCD), stratify.var)
    LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X),replace=TRUE)})
    INDEX <- unsplit(LAPPLY, stratify.var)
    
    SCdata.ABCD.subset <- SCdata.ABCD[INDEX, ]
    
    # fit models
    AIC_df <- data.frame(knots=rep(NA,4), AIC=rep(NA,4))
    for (knots in c(3:6)){
      modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, "TRUE", covariates))
      gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = SCdata.ABCD.subset)
      
      AIC_df$knots[knots-2] <- knots
      AIC_df$AIC[knots-2] <- AIC(gamm.model$mer)
    }
    
    bestK <- AIC_df$knots[which.min(AIC_df$AIC)]
    
    return(bestK)
    
  })
  
  bootchoice.df <- unlist(bootchoice)
  choice.df[[region]] <- bootchoice.df
}

choice.df$MostFrequent <- apply(choice.df, 1, function(row) {
  freq <- table(row)
  as.numeric(names(freq)[which.max(freq)])
})

table(choice.df$MostFrequent)
# 3 
# 1000

choice.df$modelnum_k3 <- rowSums(choice.df[,2:79] == 3)
write.csv(choice.df, paste0(resultFolderABCD, "/select_k_gam.csv"), row.names = F)
## Best k
ggplot(data = choice.df, aes(MostFrequent, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", position = position_dodge(width = 1), breaks=c(2.5,3.5,4.5,5.5, 6.5), center=0) +
  scale_x_continuous(breaks=c(3,4,5,6), limits = c(2.5, 6.5))+
  labs(x = expression("The optimal "*italic("k")*" value"), y = "Frequency", title = "ABCD") +
  labs(y= "Frequency")+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderABCD, '/GAM_k_choose.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderABCD, '/GAM_k_choose.svg', sep = ''), width = 12, height = 10, units = "cm")

## Number of models select k=3
ggplot(data = choice.df, aes(modelnum_k3, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", position = position_dodge(width = 1)) +
  #geom_vline(xintercept = median(choice.df$modelnum_k3), color="red", linetype="dashed", linewidth=1)+
  labs(x = expression("Number of models choosing "*italic("k")*"=3"), y = "Frequency", title = "ABCD") +
  labs(y= "Frequency")+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderABCD, '/GAM_k_choose_modelnum.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderABCD, '/GAM_k_choose_modelnum.svg', sep = ''), width = 12, height = 10, units = "cm")
################################################

# ChineseCohort k=3~6
#########################################
AIC_k3_6 <- data.frame(region=paste0("SC.", 1:78), AIC.k3=rep(NA,78), AIC.k4=rep(NA,78),
                       AIC.k5=rep(NA,78),AIC.k6=rep(NA,78))
SCdata.ChineseCohort$Sex <- as.factor(SCdata.ChineseCohort$Sex)
covariates<-"Sex+mean_fd"
dataname<-"SCdata.ChineseCohort"
smooth_var<-"Age"

## bootstrap
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(mgcv)
  library(tidyr)
  
})

SCdata.ChineseCohort$Sex <- as.factor(SCdata.ChineseCohort$Sex)
covariates<-"Sex+mean_fd"
smooth_var<-"Age"
stratify.var <- SCdata.ChineseCohort$study

choice.df <- data.frame(boottime=1:1000)
for (i in 1:78){
  region <- paste0("SC.", i, "_h")
  clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
  
  bootchoice<-parLapply(cl, 1:1000, function(x){
    
    # bootstrap
    set.seed( seed=925 + x )
    SPLIT <- split(1:NROW(SCdata.ChineseCohort), stratify.var)
    LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X),replace=TRUE)})
    INDEX <- unsplit(LAPPLY, stratify.var)
    
    SCdata.ChineseCohort.subset <- SCdata.ChineseCohort[INDEX, ]
    
    # fit models
    AIC_df <- data.frame(knots=rep(NA,4), AIC=rep(NA,4))
    for (knots in c(3:6)){
      modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, "TRUE", covariates))
      gam.model <- gam(modelformula, method = "REML", data = SCdata.ChineseCohort.subset)
      AIC_df$knots[knots-2] <- knots
      AIC_df$AIC[knots-2] <- AIC(gam.model)
    }
    
    bestK <- AIC_df$knots[which.min(AIC_df$AIC)]
    
    return(bestK)
    
  })
  
  bootchoice.df <- unlist(bootchoice)
  choice.df[[region]] <- bootchoice.df
}

choice.df$MostFrequent <- apply(choice.df, 1, function(row) {
  freq <- table(row)
  as.numeric(names(freq)[which.max(freq)])
})

table(choice.df$MostFrequent)
# 3   4   5   6 
# 814 129  27  30 

choice.df$modelnum_k3 <- rowSums(choice.df[,2:79] == 3)
write.csv(choice.df, paste0(resultFolder_ChineseCohort, "/select_k_gam.csv"), row.names = F)
## Best k
ggplot(data = choice.df, aes(MostFrequent, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", position = position_dodge(width = 1), breaks=c(2.5,3.5,4.5,5.5, 6.5), center=0) +
  scale_x_continuous(breaks=c(3,4,5,6), limits = c(2.5, 6.5))+
  labs(x = expression("The optimal "*italic("k")*" value"), y = "Frequency", title = "Chinese Cohort") +
  labs(y= "Frequency")+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderChineseCohort, '/GAM_k_choose.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderChineseCohort, '/GAM_k_choose.svg', sep = ''), width = 12, height = 10, units = "cm")

## Number of models select k=3
ggplot(data = choice.df, aes(modelnum_k3, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", position = position_dodge(width = 1)) +
  #geom_vline(xintercept = median(choice.df$modelnum_k3), color="red", linetype="dashed", linewidth=1)+
  labs(x = expression("Number of models choosing "*italic("k")*"=3"), y = "Frequency", title = "Chinese Cohort") +
  labs(y= "Frequency")+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderChineseCohort, '/GAM_k_choose_modelnum.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderChineseCohort, '/GAM_k_choose_modelnum.svg', sep = ''), width = 12, height = 10, units = "cm")
###############################################






