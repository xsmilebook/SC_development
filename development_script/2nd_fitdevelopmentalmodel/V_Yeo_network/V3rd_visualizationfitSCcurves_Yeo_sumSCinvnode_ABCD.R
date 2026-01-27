## This script is to generate fitted values from gam models.
## The predicted SC strength will be generated as the age was sampled from 8.9 to 13.8, 
## with a total of 1000 data points, covariables will be set as median or mode.
## The data will be used to draw developmental trajectories.
library(R.matlab)
library(mgcv)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(gratia)
library(RColorBrewer)
library(visreg)
rm(list = ls())
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ABCD_final/Yeo', Yeoresolution)

#### load data
CVthr = 75
gamresultsum.Yeo.delLM<-readRDS(paste0(interfileFolder, '/gamresults_Yeo',elementnum,'_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel_Yeo',elementnum,'_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution,'_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df_Yeo',elementnum,'_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/plotdata_derivatives.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
source(paste0(functionFolder, '/gammsmooth.R'))
#### generate function data
#### Yeoresolution.delLM index & SC rank
MatrixYeoresolution.delLM<-matrix(NA, nrow=Yeoresolution.delLM, ncol=Yeoresolution.delLM)
indexupYeoresolution.delLM <- upper.tri(MatrixYeoresolution.delLM)
indexsaveYeoresolution.delLM <- !indexupYeoresolution.delLM
MatrixYeoresolution.delLM.index<-MatrixYeoresolution.delLM
#(Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum
MatrixYeoresolution.delLM.index[indexsaveYeoresolution.delLM]<-c(1:elementnum)
#SC rank
MatrixYeoresolution.delLM.SCrank<-MatrixYeoresolution.delLM
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    MatrixYeoresolution.delLM.SCrank[x,y]<-(x+y)^2+(x-y)^2
  }
}
MatrixYeoresolution.delLM.SCrank[indexupYeoresolution.delLM]<-NA
MatrixYeoresolution.delLM.SCrank[indexsaveYeoresolution.delLM]<-rank(MatrixYeoresolution.delLM.SCrank[indexsaveYeoresolution.delLM], ties.method = "average")
gamresultsum.Yeo.delLM$SCrank <- MatrixYeoresolution.delLM.SCrank[indexsaveYeoresolution.delLM]
gamresultsum.Yeo.delLM$partialRsq <- as.numeric(gamresultsum.Yeo.delLM$partialRsq)
# plot data
plotdatasum <- list()
for (x in 1:elementnum){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, smooth_var="age")
  plotdata$SC_label <- names(plotdata)[14]
  plotdata[,14] <- NULL
  plotdata$SCrank <- MatrixYeoresolution.delLM.SCrank[indexsaveYeoresolution.delLM][x]
  plotdata$PartialRsq <- gamresultsum.Yeo.delLM$partialRsq[x]
  plotdata$meanderv2 <- gamresultsum.Yeo.delLM$meanderv2[x]
  
  plotdatasum[[x]] <- plotdata
}
saveRDS(plotdatasum, paste0(interfileFolder, "/plotdatasum_scale_TRUE_Yeo", Yeoresolution,".rds"))
plotdatasum.df <- do.call(rbind, lapply(plotdatasum, as.data.frame))
summary(plotdatasum.df$age)
boxplot(plotdatasum.df$PartialRsq)

## plot
maxth = max(abs(gamresultsum.Yeo.delLM$partialRsq))
ggplot()+
  geom_line(data=plotdatasum.df, aes(x=age, y=fit.ratio, group=SC_label, color=PartialRsq), size=0.8, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu",limit=c(-maxth, maxth),na.value="#053061", direction = -1)+
  labs(x="Age (years)", y="SC strength (ratio)")+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"),aspect.ratio = 0.85,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=20, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/Yeo', Yeoresolution.delLM,'_sumSCinvnode_fit/devcurve_SCrank_fit.ratio.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/Yeo', Yeoresolution.delLM,'_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.svg'),  dpi=600, width=14, height =13, units = "cm")

