library(R.matlab)
library(mgcv)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(openxlsx)
library(gratia)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(patchwork)
rm(list = ls())

resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'

#### load data
CVthr = 75
gamresultsum.SAorder.delLM<-readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
#### generate fitted values for developmental trajectories

#     Cont Default DorsAttn Limbic SalVentAttn SomMot Vis
# 1     0       0        0      0           0     16  16
# 2     0       0        0      0           1     17  14
# 3     1       0        2      0           0     15  14
# 4     2       0       11      0           2     10   7
# 5     3       1        7      1           4      9   6
# 6     1       6        9      1           4      7   3
# 7     6       8        6      0           8      2   1
# 8     4       7        7      0          12      1   0
# 9     5      13        4      0           9      0   0
# 10   11      14        0      1           5      0   0
# 11   14      16        0      0           1      0   0
# 12    5      26        0      0           0      0   0

# set cluster
num_cores <- detectCores() - 8
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
clusterEvalQ(cl, {
  library(tidyverse)
  library(R.matlab)
  library(psych)
  library(gratia)
  library(mgcv)
  library(parallel)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape)
  source(paste0(functionFolder, '/plotdata_generate.R'))
  source(paste0(functionFolder, '/colorbarvalue.R'))
})

plotdatasum<-parLapply(cl, 1:78, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, dataname = NA, "age")
  return(plotdata)
})
plotdatasum.df<-as.data.frame(matrix(NA, nrow = 1, ncol=17))
names(plotdatasum.df)<-c(names(plotdatasum[[2]])[1:13],"SC_label",  "SCrank", 
                         "PartialRsq", "meanderiv2")
#### SA12 index & SC rank
Matrix12<-matrix(NA, nrow=12, ncol=12)
indexup12 <- upper.tri(Matrix12)
indexsave12 <- !indexup12
Matrix12.index<-Matrix12
#13*12/2=78
Matrix12.index[indexsave12]<-c(1:78)
#S-A connectional axis rank
Matrix12.SCrank<-Matrix12
for (x in 1:12){
  for (y in 1:12){
    Matrix12.SCrank[x,y]<-x^2+y^2
  }
}
Matrix12.SCrank[indexup12]<-NA
Matrix12.SCrank[indexsave12]<-rank(Matrix12.SCrank[indexsave12], ties.method = "average")
# rbind plotdata
for (i in 1:78){
  tmp<-plotdatasum[[i]][,-14]
  tmp$SC_label<-names(plotdatasum[[i]])[14]
  tmp$SCrank<-Matrix12.SCrank[indexsave12][i]
  tmp$PartialRsq<-gamresultsum.SAorder.delLM$partialRsq[i]
  tmp$meanderiv2<-gamresultsum.SAorder.delLM$meanderv2[i]
  plotdatasum.df<-rbind(plotdatasum.df, tmp)
}
plotdatasum.df<-plotdatasum.df[-1, ]

# Define SC labels
Matrix_SClabel <- Matrix12
Matrix_SClabel[1:4, 1:4] <- "S-S"; Matrix_SClabel[5:12, 1:4] <- "S-A"
Matrix_SClabel[5:12, 5:12] <- "A-A"

# Matrix_SClabel[1:4, 1:4] <- "S-S"; Matrix_SClabel[5:8, 1:4] <- "S-M"; Matrix_SClabel[9:12, 1:4] <- "S-A"
# Matrix_SClabel[5:8, 5:8] <- "M-M"; Matrix_SClabel[9:12, 5:8] <- "M-A"; Matrix_SClabel[9:12, 9:12] <- "A-A"
gamresultsum.SAorder.delLM$SCtype <- Matrix_SClabel[indexsave12]

# Developmental trajectories of each SC type.
plotdatasum.df <- plotdatasum.df %>% left_join(select(gamresultsum.SAorder.delLM, c("parcel", "SCtype")),join_by(SC_label==parcel))
lwthr = min(plotdatasum.df$fit.ratio)
upthr = max(plotdatasum.df$fit.ratio)
SCtype.vec <- c("S-S", "S-A", "A-A")

for (i in c(1:3)){
   plotdf.tmp <- plotdatasum.df[which(plotdatasum.df$SCtype==SCtype.vec[i]), ]
   # if (mean(plotdf.tmp$SCrank) <= 26){
   #   palettename <- "Blues"; direction=-1
   # }else if (mean(plotdf.tmp$SCrank) > 26 & mean(plotdf.tmp$SCrank) <= 52){
   #   palettename <- "RdBu"; direction=-1
   # }else{
   #   palettename <- "Reds"; direction=1
   # }
   
   ggplot()+
     geom_line(data=plotdf.tmp, aes(x=age, y=fit.ratio, group=SC_label, color=SCrank), size=0.8, alpha=0.8)+
     scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(1,78))+
     labs(x="Age (years)", y="SC strength (ratio)", color="rank")+
     theme_classic()+
     theme(axis.text=element_text(size=18, color="black"), 
           axis.title =element_text(size=18, color="black"),aspect.ratio = 0.9,
           plot.background=element_rect(fill="transparent"),
           panel.background=element_rect(fill="transparent"),
           plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
   
  ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/SCtype/SCdev_ratio_', SCtype.vec[i],'_3Cat.tiff'), width = 10, height=9, units="cm")
  
}

# correlation within each type
gamresultsum.SAorder.delLM$SCrank <- Matrix12.SCrank[indexsave12]
gamresultsum.SAorder.delLM <- gamresultsum.SAorder.delLM %>% mutate(
  SCtype2 = case_when(str_detect(SCtype, "A")~ "A",
                      .default = "S-S")
  
)

corr.result <- data.frame(SCtype = c("S-S", "A"), corr_2nd_deriv=rep(NA, 2), corr_2nd_deriv_P=rep(NA, 2))
SCtype.vec <- c("S-S", "A")
lmthr <- max(abs(gamresultsum.SAorder.delLM$meanderv2))
for (i in c(1:2)){
  gamresultsum.tmp <- gamresultsum.SAorder.delLM[which(gamresultsum.SAorder.delLM$SCtype2==SCtype.vec[i]), ]
  corr.result$corr_2nd_deriv[i] <- corr.test(gamresultsum.tmp$meanderv2, gamresultsum.tmp$SCrank, method = "spearman")$r
  corr.result$corr_2nd_deriv_P[i] <- corr.test(gamresultsum.tmp$meanderv2, gamresultsum.tmp$SCrank, method = "spearman")$p
  
  ggplot(data=gamresultsum.tmp)+
    geom_point(aes(x=SCrank, y=meanderv2, color=meanderv2), size=5)+
    geom_smooth(aes(x=SCrank, y=meanderv2),linewidth=2, method ="lm", color="black")+
    scale_color_distiller(type="seq", palette = "RdBu", limit=c(-lmthr, lmthr), direction = -1)+
    labs(x="S-A connectional axis rank", y="Second derivative", color="rank")+
    scale_y_continuous(breaks = c(-0.004, -0.002, 0, 0.002, 0.004), labels=c(-4, -2, 0, 2, 4))+
    #scale_x_continuous(breaks = seq(from=1, to=78, by=3))+
    theme_classic()+
    theme(axis.text=element_text(size=18, color="black"), 
          axis.title =element_text(size=18, color="black"),aspect.ratio = 0.9,
          axis.ticks.y = element_blank(),
          plot.background=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
  
  ggsave(paste0(FigureFolder, '/CV',CVthr, '/correlation_sumSCinvnode_SCrank/SCtype/mean_2ndderiv_SCrank_', SCtype.vec[i],'_3Cat.tiff'), width = 11, height=9, units="cm")
}
# SCtype corr_2nd_deriv corr_2nd_deriv_P
# 1    S-S     0.80606061     0.0048620611
# 2    S-M     0.60780002     0.0125027953
# 3    S-A     0.44411765     0.0848227442
# 4    M-M     0.89090909     0.0005421442
# 5    M-A     0.63823529     0.0078011386
# 6    A-A     0.05454545     0.8810361812

#   SCtype corr_2nd_deriv corr_2nd_deriv_P
# 1    S-S      0.8060606     4.862061e-03
# 2    S-A      0.8626157     2.189237e-10
# 3    A-A      0.6419562     2.441883e-05

# SCtype corr_2nd_deriv corr_2nd_deriv_P
# 1    S-S      0.8060606     4.862061e-03
# 2      A      0.6951003     4.814100e-11


## ANOVA and post-hoc test across the 3 categories
gamresultsum.SAorder.delLM$SCrank <- Matrix12.SCrank[lower.tri(Matrix12.SCrank, diag = T)]
anova_results <- aov(meanderv2 ~ SCtype2, data = gamresultsum.SAorder.delLM)
# F = 28.53, P=6.1e-10 Df = 2, 75
# F = 
# S-S VS A

t.test(meanderv2~SCtype2, data = gamresultsum.SAorder.delLM)
# t=6.9482, df=11.26, p=2.152e-05

pairwise_result <- pairwise.t.test(gamresultsum.SAorder.delLM$meanderv2, gamresultsum.SAorder.delLM$SCtype, p.adjust.method = "fdr")
print(pairwise_result)

#     S-S     S-A 
# S-A 5.2e-09 -   
# A-A 5.1e-10 0.37

anova_results <- aov(SCrank ~ SCtype, data = gamresultsum.SAorder.delLM)
# F = 28.53, P=6.1e-10 Df = 2, 75
pairwise_result <- pairwise.t.test(gamresultsum.SAorder.delLM$SCrank, gamresultsum.SAorder.delLM$SCtype, p.adjust.method = "fdr")
print(pairwise_result)

gamresultsum.SAorder.delLM$SCtype2 <- factor(gamresultsum.SAorder.delLM$SCtype2, levels=c("S-S", "A"))
ggplot(data=gamresultsum.SAorder.delLM)+
  geom_boxplot(aes(x=SCtype2, y=meanderv2), fill="transparent", color="black")+
  labs(x=NULL, y="Second derivative")+
  theme_classic()+
  theme(axis.text=element_text(size=18, color="black"), 
        axis.title =element_text(size=18, color="black"),aspect.ratio = 0.9,
        axis.ticks.y = element_blank(),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/correlation_sumSCinvnode_SCrank/SCtype/mean_2ndderiv_2Cat.tiff'), width = 11, height=9, units="cm")


## first derivative at age of 8 and 22
derivative$SCrank <- as.numeric(derivative$SCrank)
derivative_8 <- derivative %>% filter(round(age,4)==8.0833) %>% left_join(select(gamresultsum.SAorder.delLM, c("parcel", "SCtype")), join_by(label_ID==parcel))
derivative_21 <- derivative %>% filter(round(age,4)==21.9167) %>% left_join(select(gamresultsum.SAorder.delLM, c("parcel", "SCtype")), join_by(label_ID==parcel))

corr.result.8 <- data.frame(SCtype = SCtype.vec, corr_1st_deriv8=rep(NA, 6), corr_1st_deriv8_P=rep(NA, 6))
corr.result.21 <- data.frame(SCtype = SCtype.vec, corr_1st_deriv21=rep(NA, 6), corr_1st_deriv21_P=rep(NA, 6))
lmthr8 <- max(abs(derivative_8$derivative))
lmthr21 <- max(abs(derivative_21$derivative))

for (i in c(1:6)){
  derivative_8.tmp <- derivative_8[which(derivative_8$SCtype==SCtype.vec[i]), ]
  derivative_21.tmp <- derivative_21[which(derivative_21$SCtype==SCtype.vec[i]), ]
  
  corr.result.8$corr_1st_deriv8[i] <- corr.test(derivative_8.tmp$derivative, derivative_8.tmp$SCrank, method = "spearman")$r
  corr.result.8$corr_1st_deriv8_P[i] <- corr.test(derivative_8.tmp$derivative, derivative_8.tmp$SCrank, method = "spearman")$p
  
  corr.result.21$corr_1st_deriv21[i] <- corr.test(derivative_21.tmp$derivative, derivative_21.tmp$SCrank, method = "spearman")$r
  corr.result.21$corr_1st_deriv21_P[i] <- corr.test(derivative_21.tmp$derivative, derivative_21.tmp$SCrank, method = "spearman")$p
  
  # 8 years
  ggplot(data=derivative_8.tmp)+
    geom_point(aes(x=SCrank, y=derivative, color=derivative), size=5)+
    geom_smooth(aes(x=SCrank, y=derivative),linewidth=2, method ="lm", color="black")+
    scale_color_distiller(type="seq", palette = "RdBu", limit=c(-lmthr8, lmthr8), direction = -1)+
    labs(x="S-A connectional axis rank", y="SC change rate (8 years)", color="rank")+
    scale_y_continuous(breaks = waiver(), labels=NULL)+
    theme_classic()+
    theme(axis.text=element_text(size=18, color="black"), 
          axis.title =element_text(size=18, color="black"),aspect.ratio = 0.9,
          axis.ticks.y = element_blank(),
          plot.background=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
  
  ggsave(paste0(FigureFolder, '/CV',CVthr, '/correlation_sumSCinvnode_SCrank/SCtype/mean_1stderiv8year_SCrank_', SCtype.vec[i],'.tiff'), width = 10, height=9, units="cm")
  
  # 21 years
  ggplot(data=derivative_21.tmp)+
    geom_point(aes(x=SCrank, y=derivative, color=derivative), size=5)+
    geom_smooth(aes(x=SCrank, y=derivative),linewidth=2, method ="lm", color="black")+
    scale_color_distiller(type="seq", palette = "RdBu", limit=c(-lmthr21, lmthr21), direction = -1)+
    labs(x="S-A connectional axis rank", y="SC change rate (22 years)", color="rank")+
    scale_y_continuous(breaks = waiver(), labels=NULL)+
    theme_classic()+
    theme(axis.text=element_text(size=18, color="black"), 
          axis.title =element_text(size=18, color="black"),aspect.ratio = 0.9,
          axis.ticks.y = element_blank(),
          plot.background=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
  
  ggsave(paste0(FigureFolder, '/CV',CVthr, '/correlation_sumSCinvnode_SCrank/SCtype/mean_1stderiv21year_SCrank_', SCtype.vec[i],'.tiff'), width = 10, height=9, units="cm")
  
}

corr.result.8
# SCtype corr_1st_deriv8 corr_1st_deriv8_P
# 1    S-S     -0.84242424       0.002220031
# 2    S-M     -0.47647059       0.062059174
# 3    S-A     -0.01764706       0.948280150
# 4    M-M     -0.75757576       0.011143447
# 5    M-A     -0.24705882       0.356275042
# 6    A-A      0.13939394       0.700931885

corr.result.21
# SCtype corr_1st_deriv21 corr_1st_deriv21_P
# 1    S-S       0.83030303        0.002940227
# 2    S-M       0.30588235        0.249254299
# 3    S-A       0.70000000        0.002535095
# 4    M-M       0.36969697        0.293050075
# 5    M-A       0.65000000        0.006415707
# 6    A-A       0.05454545        0.881036181







