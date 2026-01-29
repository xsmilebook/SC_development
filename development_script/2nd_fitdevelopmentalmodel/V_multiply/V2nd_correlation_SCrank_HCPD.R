library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
rm(list = ls())
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2

demopath<-'D:/xuxiaoyu/DMRI_network_development/SC_development/demopath'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA', ds.resolution)

#### load data
CVthr = 75
gamresult<-readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
#### source function
source(paste0(functionFolder, '/SCrankcorr_mult.R'))
#### convert critical ages of insignificantly developmental edges to NA
#### convert critical ages equal to age boundaries to NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method="fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)
gamresult$increase.onset[gamresult$sig==FALSE]<-NA ; gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2,2)==8.08] <- NA
gamresult$increase.offset[gamresult$sig==FALSE]<-NA ; gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2,2)==21.92] <- NA
gamresult$peak.change[gamresult$sig==FALSE]<-NA
gamresult$peak.increase.change[gamresult$sig==FALSE]<-NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")

#### 1. compute correlations to SC rank
computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change",
                     "meanderv2")
for (x in 1:5){
  computevar <- computevar.list[x]
  ds.resolution<-ds.resolution
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution,  dsdata=FALSE)
  if (x==1){
    SCrank_correlation <- correlation.df
  }else{
    SCrank_correlation <- rbind(SCrank_correlation, correlation.df)
  }
}
SCrank_correlation
# ds.resolution         Interest.var   r.spearman   p.spearman
# 1            12           partialRsq  0.006880284 9.523271e-01
# 2            12      increase.onset2 -0.008550088 9.530078e-01
# 3            12     increase.offset2  0.471121920 3.249029e-03
# 4            12 peak.increase.change  0.614498642 1.527675e-08
# 5            12            meanderv2  0.551927815 1.624245e-07

### 2. scatter plots
############################################
# meanderv2
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2

ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=SCrank), size=5)+
  geom_smooth(aes(x=SCrank, y=meanderv2),linewidth=2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.003, 0, 0.003), labels = c(-3, 0, 3))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 0.9,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '_SAmult.tiff'), width=13, height =12, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '_SAmult.svg'), dpi=600, width=20, height =14, units = "cm")







