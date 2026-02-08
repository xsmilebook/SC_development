library(tidyverse)
library(R.matlab)
library(psych)
library(gratia)
library(mgcv)
library(parallel)
library(reshape)

rm(list = ls())
CVthr = 75
elementnum = 78
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/results_HCPD'
  FigureFolder<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Figure_HCPD_final/SA12'
}else{
  resultFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HCPD'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
  functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'
}

# load data
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"))
nrow(SCdata.sum.merge)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
SCdata.diw <- readRDS(paste0(interfileFolder, "/SCdata.diw_SA12CV75.rds"))
SCdata.diw$sex <- as.factor(SCdata.diw$sex)
gamresultsum.df <- readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds'))


# source function
source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/gamderivatives.R'))
source(paste0(functionFolder, '/SCrankcorr.R'))

SCrank.df <- SCrankcorr(gamresultsum.df, "meanderv2", 12, T)
SCrank12<-SCrank.df$SCrank

# define function
compute.SC.corr <- function(drawtime){
  deriv.SA12.drawtmp <- data.frame(age=rep(NA, 78*1000), deri.pos=rep(NA, 78*1000),
                                   SClabel=rep(NA, 78*1000))
  for (i in 1:78){
    df.tmp <- derivative.posterior.df[[i]]
    df.tmp <- df.tmp[df.tmp$draw==paste0("draw", drawtime),]
    lwth <- (i-1)*1000 +1
    upth <- i*1000
    deriv.SA12.drawtmp$age[lwth:upth]<-df.tmp$age
    deriv.SA12.drawtmp$deri.pos[lwth:upth]<-df.tmp$posterior.derivative
    deriv.SA12.drawtmp$SClabel[lwth:upth]<-paste0("SC.", i)
  }
  agerange <- deriv.SA12.drawtmp$age[1:1000]
  corr.df <- data.frame(corr.pos.tmp=rep(NA,1000))
  for (j in 1:1000){
    deri.pos.tmp <- deriv.SA12.drawtmp$deri.pos[deriv.SA12.drawtmp$age==agerange[j]]
    corr.pos.tmp <- corr.test(deri.pos.tmp, SCrank12, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}


covariates<-"sex+mean_fd"
dataname<-"SCdata.diw"
smooth_var<-"age"

AcrossZeroAge <- data.frame(kvalue = c(3:6), age.0.corr.diw.median=rep(NA, 4), 
                            age.0.corr.diw.lower=rep(NA, 4), age.0.corr.diw.upper=rep(NA, 4))
for (k in c(3,4,5,6)){
  # (1) fit models
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-grep("SC.", names(SCdata.diw), value=T)[x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=k, set_fx=F, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 50)
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_over8_CV75_k', k, '_scale_TRUE.rds'))
  
  gammodelsum <- resultsum
  
  # (2) Calculate derivatives
  derivative.posterior.df<-mclapply(1:elementnum, function(x){
    SClabel.tmp <- grep("SC.", names(SCdata.diw), value=T)[x]
    modobj<-gammodelsum[[x]]
    draws<-1000
    increments<-1000
    derivdata<- gam.derivatives(modobj, "age",smoothvector=NA, draws, increments, return_posterior_derivatives = TRUE)
    derivdata$SCrank<-SCrank12[x]
    meanSC <- mean(modobj$model[,SClabel.tmp],na.rm=T)
    derivdata$meanSC<-meanSC
    maxSC <- max(modobj$model[,SClabel.tmp],na.rm=T)
    derivdata$maxSC<-maxSC
    return(derivdata)
  })
  
  saveRDS(derivative.posterior.df, paste0(resultFolder, '/derivative.posterior.df.SA12_CV75_k', k, '.rds'))
  
  # (3) Estimate age window
  deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
  rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
  colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
  
  deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
    corr.df.tmp <- compute.SC.corr(x)
    return(corr.df.tmp)
  }, mc.cores = 40)
  deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
  deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
  write.csv(deri.SCrank.posterior.corr, paste0(resultFolder, '/deri.SCrank12_CV75.posterior.diw.corr_k', k,'.csv'), row.names = F)
  
  # (4) Compute across zero age
  agerange <- unique(derivative.posterior.df[[1]]$age)
  
  # 1. age window when alignment is equal to zero
  age.0.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.min(abs(round(deri.SCrank.posterior.corr[x,], 4)-0)))])
  age.0.corr.diw <- as.numeric(unlist(age.0.corr.diw))
  age.0.corr.diw.median <- median(age.0.corr.diw) #median age #bayes
  age.0.corr.diw.CI <- quantile(age.0.corr.diw, probs = c(0.025, 0.975))
  age.0.corr.diw.lower <- age.0.corr.diw.CI[[1]]
  age.0.corr.diw.upper <- age.0.corr.diw.CI[[2]]
  AcrossZeroAge$age.0.corr.diw.median[k-2] <- age.0.corr.diw.median
  AcrossZeroAge$age.0.corr.diw.lower[k-2] <- age.0.corr.diw.lower
  AcrossZeroAge$age.0.corr.diw.upper[k-2] <- age.0.corr.diw.upper
  
  # diw
  posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
  posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
  # diw 95%CI
  posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
  posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))
  
  ##### plot alignment with S-A connectional axis correlation (Figure 3 (B))
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
  
  saveRDS(df.poscorr.diw, paste0(resultFolder, "/Alignment_plotdf_k", k,".rds"))
}

print(AcrossZeroAge)
#    kvalue age.0.corr.diw.median age.0.corr.diw.lower age.0.corr.diw.upper
# 1      3              15.68544             15.36695             16.03161
# 2      4              14.41149             13.96839             15.02077
# 3      5              14.43919             13.85761             14.99308
# 4      6              14.28687             13.82991             14.89615


