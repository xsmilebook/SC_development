#### This script is to conduct correlation analysis in ABCD data
#### between gam statistical indexes to connectional axis rank.
#### And draw scatter plots & matrix graphs.
#### Figure 5. c
#### Spearman correlations were conducted.
#### The functions of spin tests came from Váša, F. et al. (2018, https://github.com/frantisekvasa/rotate_parcellation).

library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
rm(list = ls())
demopath<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/demopath'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA12'

#### load data
CVthr = 75
gamresult<-readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gamresult$pfdr <- p.adjust(gamresult$bootstrap_pvalue, method="fdr")
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
mean_length <- readRDS(paste0(interfileFolder, '/mean_length_SA12.rds'))
meanSC <- colMeans(SCdata[,str_detect(names(SCdata), "SC.")], na.rm = T)
if (CVthr==75){
  mean_length.df <- mean_length[[1]]
  meanlength <- colMeans(mean_length.df[,1:78], na.rm=T)
}else{mean_length.df <- mean_length[[2]]; meanlength <- colMeans(mean_length.df[,1:78], na.rm=T)}
summary(meanlength)
perm.id.full<-readRDS("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA12_sphericalrotations_N10000.rds")

#### description
corr.test(meanSC, gamresult$partialRsq) # r=0.16, p=0.17
boxplot(gamresult$partialRsq)
summary(gamresult)
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### 1. calculate correlation to SC rank
gamresult$partialRsq2 <- gamresult$partialRsq
boxplot(gamresult$partialRsq2)
gamresult$partialRsq2[which.min(gamresult$partialRsq2)] <- NA # remove outlier
computevar <- "partialRsq2"
ds.resolution<-12
SCrank_correlation <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE)
computevar <- "meanderv2"
SCrank_correlation <- rbind(SCrank_correlation, SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE))
# site all
# CV75
#   ds.resolution  Interest.var r.spearman  p.spin
# 1            12    partialRsq2 -0.6825812 0.014275
# 2            12    meanderv2  0.8093480 0.001700
# CV25
#   ds.resolution  Interest.var r.spearman  p.spin
# 1            12    partialRsq2 -0.6629027 0.017925
# 2            12    meanderv2  0.7735849 0.001925

# 1.1 control for fiber length
gamresult$meanlength <- meanlength
gamresult$partialRsq_control_length[which(!is.na(gamresult$partialRsq2))] <- residuals(lm(partialRsq2~meanlength, data=gamresult))
corr.test(gamresult$partialRsq_control_length, gamresult$meanlength, method = "pearson") # r=0
SCrankcorr(gamresult, "partialRsq_control_length", 12, perm.id.full, dsdata=FALSE)
# CV75
#    ds.resolution matsize              Interest.var r.spearman  p.spin
# 1            12     376 partialRsq_control_length -0.6277786 0.012175
# CV25
# 1            12 partialRsq_control_length  -0.601396 0.0178

gamresult$meanderv2_control_length <- residuals(lm(meanderv2~meanlength, data=gamresult))
corr.test(gamresult$meanderv2_control_length, gamresult$meanlength, method = "pearson") # r=0
SCrankcorr(gamresult, "meanderv2_control_length", 12, perm.id.full, dsdata=FALSE)
# CV75
#      ds.resolution             Interest.var r.spearman  p.spin
# 1            12 meanderv2_control_length  0.8376119 0.00105
# CV25
# 1            12 meanderv2_control_length  0.7804391  5e-04

### 2. scatter plots
############################################
ds.resolution <- 12

## partial Rsq
computevar <- "partialRsq2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(gamresult$partialRsq2)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq2
colorbarvalues.Rsq <- colorbarvalues(correlation.df$partialRsq2, 0.4)
if (CVthr == 75){
  mytheme <- theme(axis.text=element_text(size=24.3, color="black"), 
                   axis.title =element_text(size=24.3),aspect.ratio = 0.8,
                   axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
                   plot.title = element_text(size=20, hjust = 0.5, vjust=2),
                   plot.background=element_rect(fill="transparent"),
                   panel.background=element_rect(fill="transparent"),
                   legend.position = "none")
  mywidth = 17.5; myheight = 15
}else{
  mytheme <- theme(axis.text=element_text(size=24, color="black"), 
                   axis.title =element_text(size=24),aspect.ratio = 1.05,axis.line = element_line(linewidth = 0.6),
                   axis.ticks = element_line(linewidth = 0.6),
                   plot.title = element_text(size=15, hjust = 0.5, vjust=0),
                   plot.subtitle = element_text(size=21, hjust = 0.9, vjust=-6),
                   plot.background=element_rect(fill="transparent"),
                   panel.background=element_rect(fill="transparent"),
                   legend.position = "none")
  mywidth = 17; myheight = 14
}


ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq2, color=partialRsq2), size=5)+
  geom_smooth(aes(x=SCrank, y=partialRsq2), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.Rsq)+
  labs(x="S-A connectional axis rank", y=expression("Age effect (partial "*italic("R")^"2"*")"))+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+mytheme
  
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=mywidth, height =myheight, units = "cm")

## mean 2nd derivatives
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(gamresult$meanderv2)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2
colorbarvalues.meanderv2 <- colorbarvalues(correlation.df$meanderv2, 0.46)

if (CVthr==25){
  mytheme <- theme(axis.text=element_text(size=23.2, color="black"), 
                   axis.title =element_text(size=23.2),aspect.ratio = 1,axis.line = element_line(linewidth = 0.6),
                   axis.ticks = element_line(linewidth = 0.6),
                   plot.title = element_text(size=15, hjust = 0.5, vjust=0),
                   plot.subtitle = element_text(size=15, hjust = 0.1, vjust=-6),
                   plot.background=element_rect(fill="transparent"),
                   panel.background=element_rect(fill="transparent"),
                   legend.position = "none")
  mywidth = myheight = 13
}else{
  mytheme <- theme(axis.text=element_text(size=26.4, color="black"), 
                   axis.title =element_text(size=26.4),aspect.ratio = 0.74,
                   axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
                   plot.title = element_text(size=20, hjust = 0.5, vjust=2),
                   plot.background=element_rect(fill="transparent"),
                   panel.background=element_rect(fill="transparent"),
                   legend.position = "none")
  mywidth = 17.5; myheight = 15
}

ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=meanderv2), size=5)+
  geom_smooth(aes(x=SCrank, y=meanderv2), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.meanderv2)+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.005,0,0.005,0.010), labels=c(-5,0,5,10))+
  theme_classic()+mytheme
  
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=mywidth, height =myheight, units = "cm")

### 3. matrix graphs for resolution at 12
############################################
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)

linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))

computevarlist <- c("partialRsq2", "meanderv2", "partialRsq_control_length", "meanderv2_control_length")
colorprob <- c(0.4, 0.46, 0.4, 0.53)
n=0
for (computevar in computevarlist){
  SCrank_correlation.df<-SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
  
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, 12)
  rownames(Matrix.tmp) <-seq(1, 12)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, 12)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  n=n+1
  colorbarvalues.tmp <- colorbarvalues(matrixtmp.df.melt$value, colorprob[n])
  
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
  if (computevar=="partialRsq2"){
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
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
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      #geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
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
  filename<-paste0(FigureFolder,'/CV', CVthr, "/Matrix12_sumSCinvnode_gamstats_Age8_22/", computevar, "_delLM_CV", CVthr,"_siteall.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}


## Plot for the condition regressed out mean fiber length
###########################
## partial Rsq
computevar <- "partialRsq_control_length"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(gamresult$partialRsq_control_length)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq_control_length
colorbarvalues.Rsq <- colorbarvalues(correlation.df$partialRsq_control_length, 0.4)
SCrank_correlation <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE)
rho <- round(SCrank_correlation$r.spearman[1], 2)
pspin <- round(SCrank_correlation$p.spin[1], 3)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq_control_length, color=partialRsq_control_length, size=0.5))+
  geom_smooth(aes(x=SCrank, y=partialRsq_control_length), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.Rsq)+
  #scale_colour_gradientn(colours = RdBu2,values=colorbarvalues.Rsq, space="Lab")+
  labs(x="S-A connectional axis rank", y=expression("Age effect (partial "*italic("R")^"2"*")"),
       subtitle=paste0("rho=", rho, "\n Pspin=", pspin))+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22),aspect.ratio = 0.9,axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size=15, hjust = 0.5, vjust=0),
        plot.subtitle = element_text(size=15, hjust = 0.9, vjust=-6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=17, height =14, units = "cm")

## meanderv2
computevar <- "meanderv2_control_length"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(gamresult$partialRsq_control_length)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2_control_length
colorbarvalues.meanderv2 <- colorbarvalues(correlation.df$meanderv2_control_length, 0.53)
SCrank_correlation <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE)
rho <- round(SCrank_correlation$r.spearman[1], 2)
pspin <- round(SCrank_correlation$p.spin[1], 3)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2_control_length, color=meanderv2_control_length, size=0.5))+
  geom_smooth(aes(x=SCrank, y=meanderv2_control_length), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.Rsq)+
  #scale_colour_gradientn(colours = RdBu2,values=colorbarvalues.Rsq, space="Lab")+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22),aspect.ratio = 0.9,axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size=15, hjust = 0.5, vjust=0),
        plot.subtitle = element_text(size=15, hjust = 0.9, vjust=-6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=17, height =14, units = "cm")



