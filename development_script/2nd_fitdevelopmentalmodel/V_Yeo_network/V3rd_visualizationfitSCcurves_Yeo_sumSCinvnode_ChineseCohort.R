## Validation: parcellate 376 regions to ds.resolution fractions.
## This script is to generate fitted values from gam models.
## The predicted SC strength will be generated as the age was sampled from 6 to 22, 
## with a total of 1000 data points, covariables will be set as mean or mode.
## The data will be used to draw developmental trajectories.
library(R.matlab)
library(mgcv)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(openxlsx)
library(gratia)
library(RColorBrewer)
library(paletteer)
rm(list = ls())
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

# in PC
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ChineseCohort'
functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ChineseCohort'
FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ChineseCohort_final/Yeo', Yeoresolution)

#### load data
CVthr = 75
gamresultsum.Yeo.delLM<-readRDS(paste0(interfileFolder, '/gamresults_Yeo',elementnum,'_sumSCinvnode_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel_Yeo',elementnum,'_sumSCinvnode_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution,'_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatAge.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df_Yeo',elementnum,'_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/plotdata_derivatives.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
#### generate function data
Agevector=SCdata$Age
# set cluster
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
clusterEvalQ(cl, {
  library(R.matlab)
  library(mgcv)
  library(psych)
  library(tidyverse)
  library(parallel)
  library(scales)
  library(openxlsx)
  library(gratia)
  library(RColorBrewer)
  library(paletteer)
  source(paste0(functionFolder, '/plotdata_generate.R'))
  source(paste0(functionFolder, '/plotdata_derivatives.R'))
  source(paste0(functionFolder, '/colorbarvalue.R'))
})

plotdatasum<-parLapply(cl, 1:elementnum, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, dataname=NA, "Age", Agevector)
  return(plotdata)
})
saveRDS(plotdatasum, paste0(interfileFolder, "/plotdatasum_scale_TRUE_Yeo", Yeoresolution,".rds"))
plotdatasum.df<-as.data.frame(matrix(NA, nrow = 1, ncol=17))
names(plotdatasum.df)<-c(names(plotdatasum[[2]])[1:13],"SC_label",  "SCrank", 
                         "PartialRsq", "meanderiv2")
#### SAds.resolution index & SC rank
Matrixds.resolution<-matrix(NA, nrow=Yeoresolution.delLM, ncol=Yeoresolution.delLM)
indexupds.resolution <- upper.tri(Matrixds.resolution)
indexsaveds.resolution <- !indexupds.resolution
Matrixds.resolution.index<-Matrixds.resolution
Matrixds.resolution.index[indexsaveds.resolution]<-c(1:elementnum)
#SC rank
Matrixds.resolution.SCrank<-Matrixds.resolution
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    Matrixds.resolution.SCrank[x,y]<-(x+y)^2+(x-y)^2
  }
}
Matrixds.resolution.SCrank[indexupds.resolution]<-NA
Matrixds.resolution.SCrank[indexsaveds.resolution]<-rank(Matrixds.resolution.SCrank[indexsaveds.resolution], ties.method = "average")
# rbind plotdata
for (i in 1:elementnum){
  tmp<-plotdatasum[[i]][,-14]
  tmp$SC_label<-names(plotdatasum[[i]])[14]
  tmp$SCrank<-Matrixds.resolution.SCrank[indexsaveds.resolution][i]
  tmp$PartialRsq<-gamresultsum.Yeo.delLM$partialRsq[i]
  tmp$meanderiv2<-gamresultsum.Yeo.delLM$meanderv2[i]
  plotdatasum.df<-rbind(plotdatasum.df, tmp)
}
plotdatasum.df<-plotdatasum.df[-1, ]
unique(plotdatasum.df$SC_label)
unique(plotdatasum.df$SCrank)
unique(plotdatasum.df$PartialRsq)

## plot
# colored based on partial R square
lmthr <- max(abs(gamresultsum.Yeo.delLM$partialRsq))
ggplot()+
  geom_line(data=plotdatasum.df, aes(x=Age, y=fit.ratio, group=SC_label, color=PartialRsq), size=1.5, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(-lmthr,lmthr))+
  labs(x="Age (years)", y="SC strength (ratio)")+
  #scale_color_manual(values = rev(brewer.pal(10, "RdBu")))+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=22, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/Yeo', Yeoresolution, '_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/Yeo', Yeoresolution, '_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.svg'), dpi=600, width=15, height =15, units = "cm")
# colored based on mean 2nd derivative
SC_label_derv2_order <- gamresultsum.Yeo.delLM$parcel[order(gamresultsum.Yeo.delLM$meanderv2)]
plotdatasum.df$SC_label2 <- factor(plotdatasum.df$SC_label, levels=SC_label_derv2_order)
lmthr <- max(abs(gamresultsum.Yeo.delLM$meanderv2))
ggplot()+
  geom_line(data=plotdatasum.df, aes(x=Age, y=fit.Z, group=SC_label2, color=meanderiv2), size=0.8, alpha=0.8)+
  #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.meanderiv2, oob = squish) +
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, limit=c(-lmthr,lmthr))+
  labs(x="Age (years)", y="SC strength (z-score)")+
  #scale_color_manual(values = rev(brewer.pal(10, "RdBu")))+
  scale_y_continuous(breaks = c(-1.5, 0.0, 1.5))+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.9,
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/Yeo', Yeoresolution, '_sumSCinvnode_fit/devcurve_SCrank_fit.Z.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/Yeo', Yeoresolution, '_sumSCinvnode_fit/devcurve_SCrank_fit.Z.svg'), dpi=600, width=16, height =16, units = "cm")

