## Validation for different matrix resolution.
## This script is to generate the alignment between developmental effect and S-A connectional axis.
library(tidyverse)
library(R.matlab)
library(psych)
library(gratia)
library(mgcv)
library(parallel)
library(ggplot2)
library(RColorBrewer)
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
# set path
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/results_ABCD'
  
}else{
  resultFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ABCD_final/Yeo', Yeoresolution)
}

#### load data
CVthr = 75
derivative.posterior.df <- readRDS(paste0(resultFolder, '/derivative.posterior.df.Yeo', Yeoresolution, '_CV', CVthr,'.rds'))
derivative.df <- readRDS(paste0(resultFolder, '/derivative.df_Yeo', elementnum, '_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, "/SCrankcorr.R"))
source(paste0(functionFolder, '/colorbarvalue.R'))
#### calculate S-A connectional axis in Yeoresolution.delLM*Yeoresolution.delLM
########################
MatrixYeoresolution.delLM<-matrix(NA, nrow=Yeoresolution.delLM, ncol=Yeoresolution.delLM)
indexupYeoresolution.delLM <- upper.tri(MatrixYeoresolution.delLM)
indexsaveYeoresolution.delLM <- !indexupYeoresolution.delLM
MatrixYeoresolution.delLM.index<-MatrixYeoresolution.delLM
#(Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum
MatrixYeoresolution.delLM.index[indexsaveYeoresolution.delLM]<-c(1:elementnum)
MatrixYeoresolution.delLM.SCrank<-MatrixYeoresolution.delLM
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    MatrixYeoresolution.delLM.SCrank[x,y]<-x^2+y^2
  }
}
MatrixYeoresolution.delLM.SCrank[indexupYeoresolution.delLM]<-NA
SCrankYeoresolution.delLM<-rank(MatrixYeoresolution.delLM.SCrank[indexsaveYeoresolution.delLM], ties.method = "average")

#### calculate correlation between posterior derivatives and S-A connectional axis
deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
# drawtime is the time of draws
compute.SC.corr <- function(drawtime){
  deriv.Yeoresolution.delLM.drawtmp <- data.frame(age=rep(NA, elementnum*1000), deri.pos=rep(NA, elementnum*1000),
                                              SClabel=rep(NA, elementnum*1000))
  for (i in 1:elementnum){
    df.tmp <- derivative.posterior.df[[i]]
    df.tmp <- df.tmp[df.tmp$draw==paste0("draw", drawtime),]
    lwth <- (i-1)*1000 +1
    upth <- i*1000
    deriv.Yeoresolution.delLM.drawtmp$age[lwth:upth]<-df.tmp$age
    deriv.Yeoresolution.delLM.drawtmp$deri.pos[lwth:upth]<-df.tmp$posterior.derivative
    deriv.Yeoresolution.delLM.drawtmp$SClabel[lwth:upth]<-paste0("SC.", i)
  }
  agerange <- deriv.Yeoresolution.delLM.drawtmp$age[1:1000]
  corr.df <- data.frame(corr.pos.tmp=rep(NA,1000))
  # estimate rho at 1,000 age points
  for (j in 1:1000){
    deri.pos.tmp <- deriv.Yeoresolution.delLM.drawtmp$deri.pos[deriv.Yeoresolution.delLM.drawtmp$age==agerange[j]]
    corr.pos.tmp <- corr.test(deri.pos.tmp, SCrankYeoresolution.delLM, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}

# compute correlation coefficients between S-A connectional axis and 1,000 posterior derivatives at 1,000 age points.
if (str_detect(wdpath, "cuizaixu_lab")){
  deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
    corr.df.tmp <- compute.SC.corr(x)
    return(corr.df.tmp)
  }, mc.cores = 40)
  
  deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
  deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
  write.csv(deri.SCrank.posterior.corr, paste0(resultFolder, '/deri.SCrank_Yeo', Yeoresolution, '_CV', CVthr,'.posterior.diw.corr.csv'), row.names = F)
}

###### extract age of maximal / zero S-A connectional axis correlation: posterior median value + 95% CI
deri.SCrank.posterior.corr <- read.csv(paste0(resultFolder, '/deri.SCrank_Yeo', Yeoresolution, '_CV', CVthr,'.posterior.diw.corr.csv'))
agerange <- unique(derivative.posterior.df[[1]]$age)
#### median corr and 95% CI
posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))

##### plot alignment with S-A connectional axis correlation
#############################################
df.poscorr.diw <- data.frame(age=agerange, median=posterior.corr.diw.median, up.95CI=posterior.corr.diw.CI$X97.5.,
                             lw.95CI=posterior.corr.diw.CI$X2.5.)
loess.median <- loess(median~age, data=df.poscorr.diw, span=0.2)
loess.lw <- loess(lw.95CI~age, data=df.poscorr.diw, span=0.2)
loess.up <- loess(up.95CI~age, data=df.poscorr.diw, span=0.2)
df.poscorr.diw$median.loess <- loess.median$fitted
df.poscorr.diw$lw.95CI.loess <- loess.lw$fitted
df.poscorr.diw$up.95CI.loess <- loess.up$fitted

ggplot(data=df.poscorr.diw)+
  geom_ribbon(aes(x=age, ymin=lw.95CI.loess, ymax=up.95CI.loess), alpha=0.3)+
  geom_line(aes(x=age, y=median.loess), size=1)+
  #scale_y_continuous(breaks=c(-0.5, -0.3, 0.0))+
  theme_classic()+
  labs(x="Age (years)", y="Alignment with S-A connectional axis (rho)")+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.line = element_line(linewidth = 0.55),
        axis.ticks = element_line(linewidth = 0.55),
        aspect.ratio =1.1,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),plot.margin = margin(t=10, r=0, b=0, l=0),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/Yeo', Yeoresolution, '_posDeriv_divweight_corr.tiff'), width=14, height=14, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Yeo', Yeoresolution, '_posDeriv_divweight_corr.svg'), width=13.5, height=15, units="cm")


### correlation plot at different age
#############################################
agerange <- unique(derivative.df$age)
df.poscorr.diw$median[df.poscorr.diw$age==min(agerange)]
# age = 8.916667
min(agerange)
df.age8 <-as.data.frame(matrix(NA, nrow=1, ncol=9))
names(df.age8)<-names(derivative.df)
for (i in 1:elementnum){
  df.tmp <- derivative.df[derivative.df$label_ID==paste0("SC.", i, "_h"),]
  df.tmp <- df.tmp[df.tmp$age==min(agerange), ]
  df.age8 <- rbind(df.age8, df.tmp)
}
df.age8 <- df.age8[-1,]

SCrankcorr(df.age8, "derivative", Yeoresolution.delLM) 
# ds.resolution Interest.var r.spearman   p.spearman
# 1         15   derivative -0.5106294 2.548003e-09

# age=13.7
ntmp<-which.min(abs(agerange-max(agerange))) 
agerange[ntmp]
df.age13 <-as.data.frame(matrix(NA, nrow=1, ncol=9))
names(df.age13)<-names(derivative.df)
for (i in 1:elementnum){
  df.tmp <- derivative.df[derivative.df$label_ID==paste0("SC.", i, "_h"),]
  df.tmp <- df.tmp[df.tmp$age==agerange[ntmp], ]
  df.age13 <- rbind(df.age13, df.tmp)
}
df.age13 <- df.age13[-1,]
SCrankcorr(df.age13, "derivative", Yeoresolution.delLM) 
# ds.resolution Interest.var r.spearman p.spearman
# 1       15   derivative -0.1351129  0.1411961

## display 2 lines in one plot
df.agemerge <- rbind(df.age8, df.age13)
df.agemerge$age <- as.factor(df.agemerge$age)
ggplot(data=df.agemerge)+
  geom_smooth(aes(x=rep(SCrankYeoresolution.delLM, 2), y=derivative, group=age), linewidth=1.2,color="black", method ="lm", 
              se=T)+
  #scale_x_continuous(breaks=c(0,40,80))+
  labs(x="S-A connectional axis rank", y="SC change rate")+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22),aspect.ratio = 0.9,
        plot.title = element_text(size=15, hjust = 0.1, vjust=-5),
        legend.position = "none", plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/deri.diw_corr_SCrank_Yeo', Yeoresolution, 'ageAll.tiff'), width = 20, height = 14, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/deri.diw_corr_SCrank_Yeo', Yeoresolution, 'ageAll.svg'), dpi=600, width = 14, height = 12, units="cm")
#############################################






