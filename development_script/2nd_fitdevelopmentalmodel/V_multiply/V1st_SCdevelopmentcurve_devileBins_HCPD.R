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
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2

resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA', ds.resolution)

#### load data
CVthr = 75
gamresultsum.SAorder.delLM<-readRDS(paste0(interfileFolder, '/gamresults',elementnum,'_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel',elementnum,'_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution,'_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df',elementnum,'_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/plotdata_derivatives.R'))
#### generate function data
plotdatasum<-mclapply(1:elementnum, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, "age")
  return(plotdata)
}, mc.cores = 2)
plotdatasum.df<-as.data.frame(matrix(NA, nrow = 1, ncol=17))
names(plotdatasum.df)<-c(names(plotdatasum[[2]])[1:13],"SC_label",  "SCrank", 
                         "PartialRsq", "meanderiv2")
#### SAds.resolution index & SC rank
Matrixds.resolution<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
indexupds.resolution <- upper.tri(Matrixds.resolution)
indexsaveds.resolution <- !indexupds.resolution
Matrixds.resolution.index<-Matrixds.resolution
Matrixds.resolution.index[indexsaveds.resolution]<-c(1:elementnum)
#SC rank
Matrixds.resolution.SCrank<-Matrixds.resolution
for (x in 1:ds.resolution){
  for (y in 1:ds.resolution){
    Matrixds.resolution.SCrank[x,y]<-x*y
  }
}
Matrixds.resolution.SCrank[indexupds.resolution]<-NA
Matrixds.resolution.SCrank[indexsaveds.resolution]<-rank(Matrixds.resolution.SCrank[indexsaveds.resolution], ties.method = "average")
# rbind plotdata
for (i in 1:elementnum){
  tmp<-plotdatasum[[i]][,-14]
  
  tmp$SC_label<-names(plotdatasum[[i]])[14]
  tmp$SCrank<-Matrixds.resolution.SCrank[indexsaveds.resolution][i]
  tmp$PartialRsq<-gamresultsum.SAorder.delLM$partialRsq[i]
  tmp$meanderiv2<-gamresultsum.SAorder.delLM$meanderv2[i]
  plotdatasum.df<-rbind(plotdatasum.df, tmp)
}
plotdatasum.df<-plotdatasum.df[-1, ]
unique(plotdatasum.df$SC_label)
unique(plotdatasum.df$SCrank)
unique(plotdatasum.df$PartialRsq)

## Average fitted values for 10 deciles of connectional axis
SAds.resolution_10 <- data.frame(SCrank=Matrixds.resolution.SCrank[indexsaveds.resolution])
SAds.resolution_10 <- SAds.resolution_10 %>%
  mutate(decile = ntile(SCrank, 10))
#table(SAds.resolution_10$decile)
SAds.resolution_10$SC_label <- paste0("SC.", c(1:78), "_h")
write.csv(SAds.resolution_10, paste0(interfileFolder, '/SAds.resolution_10.csv'), row.names = F)
plotdatasum.df.label <- merge(plotdatasum.df, SAds.resolution_10, by="SC_label", all.x=T)
#names(plotdatasum.df.label)
plotdatasum.df.decile <- plotdatasum.df.label %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(fit), SCranktype_order=mean(decile))

#names(plotdatasum.df.decile)
plotdatasum.df.decile <- plotdatasum.df.decile %>%
  group_by(decile) %>%
  mutate(fit.Z = scale(fit.avg), fit.ratio=fit.avg/fit.avg[1])

ggplot(data=plotdatasum.df.decile, aes(x=age, y=fit.Z, group=decile, color=decile))+
  geom_line(size=1.5, alpha=0.8)+
  #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.meanderiv2, oob = squish) +
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="Age (years)", y="SC strength (z-score)")+
  #scale_color_manual(values = rev(brewer.pal(6, "RdBu")))+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),legend.position = "none",
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5))
ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SAds.resolution_decile_sumSCinvnode_fit/devcurve_SCrank_multi_i*j_fit.Z_SCtype10.tiff'), dpi=600, width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SAds.resolution_decile_sumSCinvnode_fit/devcurve_SCrank_multi_i*j_fit.Z_SCtype10.svg'), dpi=600, width=15, height =13, units = "cm")





