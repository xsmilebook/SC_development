## Validation: large-scale matrix of Yeo 17 or 7
#### This script is to conduct correlation analysis 
#### between gam statistical indexes to S-A connectional axis rank.
#### And draw scatter plots & matrix graphs.
library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)

rm(list = ls())
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

demopath<-'D:/xuxiaoyu/DMRI_network_development/SC_development/demopath'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ABCD_final/Yeo', Yeoresolution)

#### load data
CVthr = 75
gamresult_noscale<-readRDS(paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'.rds'))
gamresult<-readRDS(paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gamresult$pfdr <- p.adjust(gamresult$bootstrap_pvalue, method="fdr")
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### 1. calculate correlation to SC rank
boxplot(gamresult$partialRsq)
gamresult <- within(gamresult, 
                    {partialRsq2 <- partialRsq
                    partialRsq2[which(partialRsq2>mean(partialRsq)+3*sd(partialRsq) | partialRsq2<mean(partialRsq)-3*sd(partialRsq))] <- NA
                    })
summary(gamresult$partialRsq2)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.0210804 -0.0002204  0.0006102  0.0020659  0.0043751  0.0263147          3
computevar <- "partialRsq2"
ds.resolution<-Yeoresolution.delLM
SCrank_correlation <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=FALSE)
computevar <- "meanderv2"
SCrank_correlation <- rbind(SCrank_correlation, SCrankcorr(gamresult, computevar, ds.resolution, dsdata=FALSE))

# matrix  of Yeo7
# ds.resolution Interest.var r.spearman   p.spearman
# 1             6  partialRsq2 -0.6935065 4.898974e-04
# 2             6    meanderv2  0.7701299 4.437647e-05

# matrix  of Yeo17
# ds.resolution Interest.var r.spearman   p.spearman
# 1            15  partialRsq2 -0.4235960 1.948253e-06
# 2            15    meanderv2  0.5270858 6.227975e-10

### 2. scatter plots
############################################
ds.resolution <- ds.resolution

## partial Rsq
computevar <- "partialRsq2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution,dsdata=TRUE)
summary(correlation.df$partialRsq2)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq2
maxth <- max(abs(correlation.df$partialRsq2))
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq2, color=partialRsq2), size=5)+
  geom_smooth(aes(x=SCrank, y=partialRsq2), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, limit=c(-maxth, maxth))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120))+
  labs(x="S-A connectional axis rank", y=expression("Age effect (partial "*italic("R")^"2"*")"))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23, color="black"),aspect.ratio = 1,
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=13, height =12, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=20, height =14, units = "cm")

## mean 2nd derivatives
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
summary(gamresult$meanderv2)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2
maxth <- max(abs(correlation.df$meanderv2))

ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=meanderv2), size=5)+
  geom_smooth(aes(x=SCrank, y=meanderv2), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, limit=c(-maxth, maxth))+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.015, -0.010, -0.005,0, 0.005, 0.010), labels = c(-15,-10, -5, 0, 5, 10))+
  #scale_y_continuous(breaks = c(-0.003,0, 0.003, 0.006), labels = c(-3, 0, 3, 6))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 0.9,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=12, height =11.5, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=17.5, height =15, units = "cm")


### 3. matrix graphs for resolution at ds.resolution
############################################
Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)

linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))

computevarlist <- c("partialRsq2", "meanderv2")
maxthr <- c(max(abs(gamresult$partialRsq2)), max(abs(gamresult$meanderv2)))
n=0
for (computevar in computevarlist){
  SCrank_correlation.df<-SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
  
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, ds.resolution)
  rownames(Matrix.tmp) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  n=n+1
  
  Matrix.tmp.sig <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
  Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresult$pfdr<0.05)
  Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
  colnames(Matrix.tmp.sig) <-seq(1, ds.resolution)
  rownames(Matrix.tmp.sig) <-seq(1, ds.resolution)
  matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
  matrixtmp.df.sig$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
  matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
  matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
  matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
  matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value==0),]
  if (computevar=="partialRsq2"){
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu", limit=c(-maxthr[1], maxthr[1]), na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu", limit=c(-maxthr[1], maxthr[1]), na.value = "#053061")+
      geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=10)+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
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
      scale_fill_distiller(type="seq", palette = "RdBu", limit=c(-maxthr[2], maxthr[2]), na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu", limit=c(-maxthr[2], maxthr[2]), na.value = "#053061")+
      #geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
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
  filename<-paste0(FigureFolder,'/CV', CVthr, "/Matrix_Yeo", Yeoresolution, "_sumSCinvnode_gamstats_Age8_22/", computevar, "_delLM_CV", CVthr,"_siteall.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}









