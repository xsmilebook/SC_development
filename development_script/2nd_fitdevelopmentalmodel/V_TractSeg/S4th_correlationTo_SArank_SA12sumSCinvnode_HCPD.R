#### This script is to conduct correlation analysis 
#### between gam statistical indexes to connectional axis rank.
#### And draw scatter plots & matrix graphs.
#### Figure 2. (a, c); Figure 3. (b, c); Supplementary Figure 2 (a,b); 
#### Supplementary Figure 4; Supplementary Figure 6 (a,b,d,e,g,h)
#### Spearman correlations were conducted.
#### The functions of spin tests came from Váša, F. et al. (2018, https://github.com/frantisekvasa/rotate_parcellation).
#### Network spin tests refer to Hansen et al., 2022, Nature Communications (https://github.com/netneurolab/hansen_crossdisorder_vulnerability/blob/main/code/03_disorder_similarity.py)
#### line 185 ~ 201
library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
rm(list = ls())
demopath<-'D:/xuxiaoyu/DMRI_network_development/SC_development/demopath'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'

#### load data
CVthr = 75
gamresult<-readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE_TractSeg.rds'))
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method="fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.TractSeg.merge.rds'))
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
source('D:/xuxiaoyu/DMRI_network_development/Normative_model/functions/plotmatrix.R')
#### description
meanSC <- colMeans(SCdata[,which(str_detect(names(SCdata), "SC."))])
corr.test(meanSC, gamresult$partialRsq) # r=0.55, p=0
#### convert critical ages of insignificantly developmental edges to NA
#### convert critical ages equal to age boundaries to NA
gamresult$increase.onset[gamresult$sig==FALSE]<-NA ; gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2,2)==8.08] <- NA
gamresult$increase.offset[gamresult$sig==FALSE]<-NA ; gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2,2)==21.92] <- NA
gamresult$peak.change[gamresult$sig==FALSE]<-NA
gamresult$peak.increase.change[gamresult$sig==FALSE]<-NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
summary(gamresult)


#### 1. compute correlations to SC rank
computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change",
                     "meanderv2")
for (x in 1:5){
  computevar <- computevar.list[x]
  ds.resolution<-12
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=FALSE)
  if (x==1){
    SCrank_correlation <- correlation.df
  }else{
    SCrank_correlation <- rbind(SCrank_correlation, correlation.df)
  }
}
SCrank_correlation
# CV 75th
#   ds.resolution         Interest.var r.spearman   p.spearman
# 1            12           partialRsq -0.2619885 2.049864e-02
# 2            12      increase.onset2  0.5045372 3.798861e-03
# 3            12     increase.offset2  0.4544212 3.659540e-03
# 4            12 peak.increase.change  0.8644552 1.501852e-16
# 5            12            meanderv2  0.7914159 6.495174e-18



### 2. scatter plots
############################################
ds.resolution <- 12

# meanderv2
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)


mytheme <- theme(axis.text=element_text(size=23, color="black"), 
                 axis.title =element_text(size=23),aspect.ratio = 0.9,
                 axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
                 plot.title = element_text(size=20, hjust = 0.5, vjust=2),
                 plot.background=element_rect(fill="transparent"),
                 panel.background=element_rect(fill="transparent"),
                 legend.position = "none")

ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=SCrank), size=5)+
  geom_smooth(aes(x=SCrank, y=meanderv2),linewidth=2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.004, -0.002, 0, 0.002), labels = c(-4, -2, 0, 2))+
  theme_classic()+mytheme
  
ggsave(paste0(FigureFolder,'/CV', CVthr, '/TractSeg/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=13, height =12, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/TractSeg/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=width, height =height, units = "cm")


### 3. matrix graphs for resolution of 12
############################################
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2", "meanderv2_control_length")
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))
SCrank_correlation.df<-lapply(1:5, function(x){
  computevar <- computevar.list[x]
  ds.resolution<-12
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution,dsdata=TRUE)
  return(correlation.df)
})
colorbar.prob <- c(0.5, 0.4, 0.6, 0.5, 0.5)


for (i in 1:5){
  computevar <- computevar.list[i]
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[[i]][,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, 12)
  rownames(Matrix.tmp) <-seq(1, 12)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, 12)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  colorbarvalues.tmp <- colorbarvalues(SCrank_correlation.df[[i]][,2], colorbar.prob[i])
  if (computevar=="partialRsq"){
    lmthr <- max(abs(gamresult$partialRsq))
    Matrix.tmp.sig <- matrix(NA, nrow = 12, ncol=12)
    Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresult$pfdr<0.05)
    Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
    colnames(Matrix.tmp.sig) <-seq(1, 12)
    rownames(Matrix.tmp.sig) <-seq(1, 12)
    matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
    matrixtmp.df.sig$nodeid <- seq(1, 12)
    matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
    matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
    matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
    matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
    matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value==0),]
    
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
      geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
      ggtitle(label = computevar)+labs(x=NULL, y=NULL)+
      scale_y_continuous(breaks=NULL, labels = NULL)+
      scale_x_continuous(breaks=NULL, labels = NULL)+
      theme(axis.line = element_blank(), 
            #axis.ticks=element_line(linewidth = 0),
            axis.text.x=element_text(size=12, angle=45, hjust=1), 
            axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
            axis.title =element_text(size=18),
            plot.title = element_text(size=18, hjust = 0.5),
            legend.title=element_text(size=18),
            legend.text=element_text(size=18), 
            panel.background=element_rect(fill=NA),
            panel.grid.major=element_line(linewidth = 0), 
            panel.grid.minor=element_line(linewidth = 1))
    
  }else{
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      #paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.tmp,oob = squish) +
      #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1,values=colorbarvalues.tmp, oob = squish) +
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "grey")+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
      ggtitle(label = computevar)+labs(x=NULL, y=NULL)+
      scale_y_continuous(breaks=NULL, labels = NULL)+
      scale_x_continuous(breaks=NULL, labels = NULL)+
      theme(axis.line = element_blank(), 
            #axis.ticks=element_line(linewidth = 0),
            axis.text.x=element_text(size=12, angle=45, hjust=1), 
            axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
            axis.title =element_text(size=18),
            plot.title = element_text(size=18, hjust = 0.5),
            legend.title=element_text(size=18),
            legend.text=element_text(size=18), 
            panel.background=element_rect(fill=NA),
            panel.grid.major=element_line(linewidth = 0), 
            panel.grid.minor=element_line(linewidth = 1))
  }
  
  Fig
  filename<-paste0(FigureFolder,"/CV", CVthr,  "/TractSeg/Matrix12_sumSCinvnode_gamstats_Age8_22/", computevar, "_12net_delLM_CV75.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}



