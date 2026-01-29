## Validation for different matrix resolution in the HCP-D dataset.
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
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/results_HCPD'
  
}else{
  resultFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HCPD'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
  functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/Yeo', Yeoresolution)
}

#### load data
CVthr = 75
derivative.posterior.df <- readRDS(paste0(resultFolder, '/derivative.posterior.df.Yeo', Yeoresolution, '_CV', CVthr,'.rds'))
derivative.df <- readRDS(paste0(resultFolder, '/derivative.df_Yeo', elementnum, '_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))
source(paste0(functionFolder, "/SCrankcorr.R"))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### calculate S-A connectional axis in Yeoresolution.delLM*Yeoresolution.delLM
########################
Matrixds.resolution<-matrix(NA, nrow=Yeoresolution.delLM, ncol=Yeoresolution.delLM)
indexupds.resolution <- upper.tri(Matrixds.resolution)
indexsaveds.resolution <- !indexupds.resolution
Matrixds.resolution.index<-Matrixds.resolution
#(Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum
Matrixds.resolution.index[indexsaveds.resolution]<-c(1:elementnum)
Matrixds.resolution.SCrank<-Matrixds.resolution
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    Matrixds.resolution.SCrank[x,y]<-x^2+y^2
  }
}
Matrixds.resolution.SCrank[indexupds.resolution]<-NA
SCrankds.resolution<-rank(Matrixds.resolution.SCrank[indexsaveds.resolution], ties.method = "average")

#### calculate correlation between posterior derivatives and S-A connectional axis
deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
# drawtime is the time of draws
compute.SC.corr <- function(drawtime){
  deriv.Yeo.resolution.drawtmp <- data.frame(age=rep(NA, elementnum*1000), deri.pos=rep(NA, elementnum*1000),
                                              SClabel=rep(NA, elementnum*1000))
  for (i in 1:elementnum){
    df.tmp <- derivative.posterior.df[[i]]
    df.tmp <- df.tmp[df.tmp$draw==paste0("draw", drawtime),]
    lwth <- (i-1)*1000 +1
    upth <- i*1000
    deriv.Yeo.resolution.drawtmp$age[lwth:upth]<-df.tmp$age
    deriv.Yeo.resolution.drawtmp$deri.pos[lwth:upth]<-df.tmp$posterior.derivative
    deriv.Yeo.resolution.drawtmp$SClabel[lwth:upth]<-paste0("SC.", i)
  }
  agerange <- deriv.Yeo.resolution.drawtmp$age[1:1000]
  corr.df <- data.frame(corr.pos.tmp=rep(NA,1000))
  for (j in 1:1000){
    deri.pos.tmp <- deriv.Yeo.resolution.drawtmp$deri.pos[deriv.Yeo.resolution.drawtmp$age==agerange[j]]
    corr.pos.tmp <- corr.test(deri.pos.tmp, SCrankds.resolution, method = "spearman")$r
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

age.0.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.min(abs(round(deri.SCrank.posterior.corr[x,], 4)-0)))])
age.0.corr.diw <- as.numeric(unlist(age.0.corr.diw))
age.0.corr.diw.median <- median(age.0.corr.diw) #median age 
age.0.corr.diw.CI <- quantile(age.0.corr.diw, probs = c(0.025, 0.975)) #compute the credible interval based on all draws
age.0.corr.diw.lower <- age.0.corr.diw.CI[[1]]
age.0.corr.diw.upper <- age.0.corr.diw.CI[[2]]
age.0.corr.diw.median # 16.07316 (Yeo7) 14.63 (Yeo17)
age.0.corr.diw.CI # 15.27002 17.25017 (Yeo7) 14.06462 15.18728 (Yeo17)
#### median corr and 95% CI
posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))

##### plot alignment with S-A connectional axis correlation
#############################################
df.poscorr.diw <- data.frame(age=agerange, median=posterior.corr.diw.median, up.95CI=posterior.corr.diw.CI$X97.5.,
                             lw.95CI=posterior.corr.diw.CI$X2.5.)
df.poscorr.diw$zero.corr.CI <- (df.poscorr.diw$age > age.0.corr.diw.lower & df.poscorr.diw$age < age.0.corr.diw.upper)
df.poscorr.diw$zero.corr.window <-df.poscorr.diw$age * df.poscorr.diw$zero.corr.CI
df.poscorr.diw$zero.corr.window[df.poscorr.diw$zero.corr.window==0] <-NA
loess.median <- loess(median~age, data=df.poscorr.diw, span=0.2)
loess.lw <- loess(lw.95CI~age, data=df.poscorr.diw, span=0.2)
loess.up <- loess(up.95CI~age, data=df.poscorr.diw, span=0.2)
df.poscorr.diw$median.loess <- loess.median$fitted
df.poscorr.diw$lw.95CI.loess <- loess.lw$fitted
df.poscorr.diw$up.95CI.loess <- loess.up$fitted

ggplot(data=df.poscorr.diw)+
  geom_ribbon(aes(x=age, ymin=lw.95CI.loess, ymax=up.95CI.loess), alpha=0.3)+
  geom_line(aes(x=age, y=median.loess), linewidth=1.5)+
  #geom_ribbon(aes(x=max.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#595959")+
  geom_ribbon(aes(x=zero.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#F8B01B", alpha=1)+
  geom_ribbon(aes(x=zero.corr.window, ymin=median.loess-0.01, ymax=median.loess+0.01), fill="#B2182B", alpha=1)+
  #geom_vline(aes(xintercept=age.0.corr.diw.median), color="black",linetype="dashed")+
  theme_classic()+
  labs(x="Age (years)", y="rho")+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 1,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")

ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/Yeo', Yeoresolution, '_posDeriv_divweight_corr.tiff'), width=13.5, height=13.5, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Yeo', Yeoresolution, '_posDeriv_divweight_corr.svg'), dpi=600, width=17.5, height =15, units="cm")

### plot age distribution with 0 corr
ggplot() +
  geom_histogram(aes(age.0.corr.diw, y = ..count..),binwidth = 0.5, linewidth=1,
                 color = "black", fill = "white") +
  geom_vline(aes(xintercept = age.0.corr.diw.median), colour = "red", linetype="solid", linewidth=1)+
  labs(x = NULL, y = NULL) +theme_classic()+
  scale_y_discrete(breaks=NULL)+scale_x_continuous(breaks=NULL)+
  theme(axis.line = element_blank(),
        aspect.ratio = 0.5,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.title = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 20))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Agedistribution_0corr.tiff'), width=6, height = 6, units="cm")
ggsave(paste0(FigureFolder, '/CV', CVthr, '/Alignment_development/Agedistribution_0corr.svg'), width=6, height =6, units="cm")
