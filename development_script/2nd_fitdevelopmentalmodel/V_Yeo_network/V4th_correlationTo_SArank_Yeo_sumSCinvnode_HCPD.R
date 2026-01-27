## Validation: large-scale matrix of Yeo network
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
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/Yeo', Yeoresolution)

#### load data
CVthr = 75
gamresult<-readRDS(paste0(interfileFolder, '/gamresults_Yeo', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))

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

gamresult <- gamresult %>% mutate(
  meanderv2_2 <- meanderv2
)
gamresult$meanderv2_2[(gamresult$meanderv2<(mean(gamresult$meanderv2)-3*sd(gamresult$meanderv2)))| 
                        (gamresult$meanderv2>(mean(gamresult$meanderv2)+3*sd(gamresult$meanderv2)))] <- NA

#### 1. compute correlations to SC rank
computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change",
                     "meanderv2_2")
for (x in 1:5){
  computevar <- computevar.list[x]
  ds.resolution<-Yeoresolution.delLM
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution,  dsdata=FALSE)
  if (x==1){
    SCrank_correlation <- correlation.df
  }else{
    SCrank_correlation <- rbind(SCrank_correlation, correlation.df)
  }
}
SCrank_correlation
# Yeo 7
# ds.resolution         Interest.var r.spearman   p.spearman
# 1             6           partialRsq  0.1116883 6.298112e-01
# 2             6      increase.onset2  0.4230769 1.497485e-01
# 3             6     increase.offset2  0.3214286 4.820720e-01
# 4             6 peak.increase.change  0.8483073 3.286117e-05
# 5             6            meanderv2  0.8194805 5.497437e-06

# Yeo 17
# ds.resolution         Interest.var r.spearman   p.spearman
# 1            15           partialRsq  0.1804330 4.860129e-02
# 2            15      increase.onset2  0.1955456 1.343092e-01
# 3            15     increase.offset2  0.2922160 6.729510e-02
# 4            15 peak.increase.change  0.5804978 4.786538e-08
# 5            15          meanderv2_2  0.4261272 1.352728e-06


### 2. scatter plots
############################################
ds.resolution <- Yeoresolution.delLM

## partial Rsq
computevar <- "partialRsq"

correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
summary(gamresult$partialRsq)
lmthr <- max(abs(gamresult$partialRsq))
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq, color=partialRsq, size=0.5))+
  geom_smooth(aes(x=SCrank, y=partialRsq), method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(-lmthr,lmthr))+  labs(x="S-A connectional axis rank", y="Partial R2")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder, '/CV', CVthr,'/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV', CVthr,'/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=17, height =14, units = "cm")

# meanderv2
computevar <- "meanderv2_2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2_2
lmthr <- max(abs(gamresult$meanderv2_2), na.rm=T)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2_2, color=SCrank), size=5.5)+
  geom_smooth(aes(x=SCrank, y=meanderv2_2),linewidth=2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120))+
  scale_y_continuous(breaks = c(-0.002, 0, 0.002), labels = c(-2, 0, 2))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 0.9,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=13, height =12, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=17.5, height =15, units = "cm")


### 3. matrix graphs for resolution of ds.resolution
############################################
Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2_2")
linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))

SCrank_correlation.df <- list()
for (x in 1:5){
  computevar <- computevar.list[x]
  ds.resolution<-ds.resolution
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
  SCrank_correlation.df[[x]] <- correlation.df
}

for (i in 1:5){
  computevar <- computevar.list[i]
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[[i]][,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, ds.resolution)
  rownames(Matrix.tmp) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  if (computevar=="partialRsq"){
    lmthr <- max(abs(gamresult$partialRsq))
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
    
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
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
    lmthr <- max(abs(gamresult[[computevar.list[i]]]), na.rm=T)
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "#053061")+
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
  filename<-paste0(FigureFolder,'/CV', CVthr,  "/Matrix", ds.resolution, "_sumSCinvnode_gamstats_Age8_22/", computevar, "_", ds.resolution, "net_delLM_CV75.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}



