# Fig. S5
library(tidyverse)
library(parallel)
library(scales)
library(openxlsx)
library(gratia)
library(RColorBrewer)
library(paletteer)
rm(list = ls())

resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'

#### load data
CVthr = 75
derivative <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### Connectional axis rank
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
Matrix12.SCrank[indexsave12]<-rank(Matrix12.SCrank[indexsave12], ties.method = "first")
SCrank <- data.frame(label_ID=paste0("SC.",1:78,"_h"), SCrank=Matrix12.SCrank[indexsave12])
derivative <- derivative %>% left_join(SCrank, by="label_ID")
derivative$se_CI <- (derivative$upper-derivative$lower) / (2*1.96)
derivative$p.value <- 2 * (1 - pnorm(abs(derivative$derivative / derivative$se_CI)))
summary(derivative$p.value)
derivative$significance_pvalue <- (derivative$p.value < 0.05)
derivative <- derivative %>% group_by(age) %>%
  mutate(fdr_pvalue = p.adjust(p.value, method="fdr")) %>% ungroup()
summary(derivative$fdr_pvalue)
derivative$significance_pvalue_fdr <- (derivative$fdr_pvalue < 0.05)
derivative$significant.derivative[which(derivative$significant.derivative==0)] <- NA
derivative$significant.derivative_fdr <- derivative$derivative
derivative$significant.derivative_fdr[which(derivative$significance_pvalue_fdr==F)] <- NA
#### plot significant derivatives
derivative <- derivative[order(derivative$SCrank),]
derivative$SCrank <- factor(derivative$SCrank, levels=as.character(1:78))
derivative$SCrank_label <- paste0("label", derivative$SCrank) 
derivative$SCrank_label <- factor(derivative$SCrank_label, levels=paste0("label", 78:1))
saveRDS(derivative, paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
ggplot(data=derivative)+
  geom_bar(aes(x=age, y=c(rep(1,78000)), fill = significant.derivative_fdr, group=SCrank_label,
               color=significant.derivative_fdr),stat = "identity", position = "stack")+
  scale_fill_distiller(type="seq", palette = "Reds", direction = 1,na.value ="white") +
  scale_color_distiller(type="seq", palette = "Reds", direction = 1,na.value ="white") +
  #scale_y_continuous(breaks = NULL)+
  #geom_hline(yintercept = 78)+
  ylab("S-A connectional axis rank")+xlab("Age (years)")+ggtitle(label="HCP-D")+
  labs(color="SC change rate", fill="SC change rate")+
  #scale_x_continuous(breaks = NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=20, color='black'),
        axis.title = element_text(size = 20),aspect.ratio = 1,
        plot.title=element_text(size=20, color='black', hjust=0.5),
        legend.text = element_text(size=20, color='black'), 
        legend.title = element_text(size=20, color='black'))
ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SA12_sumSCinvnode_fit/significant_derivative.tiff'), width=20, height =14, units = "cm")





