# generate age distribution plots and demographic tables.
library(ggplot2)
library(RColorBrewer)
library('reshape2')
library(gdata)
library(tableone)
display.brewer.all()

rm(list = ls())
interfileFolderABCD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
interfileFolderHCPD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
interfileFolderChineseCohort <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ChineseCohort'
interfileFolderHBN <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HBN'

FigureFolderABCD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ABCD_final'
FigureFolderHCPD <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final'
FigureFolderChineseCohort <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ChineseCohort_final'
FigureFolderHBN <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HBN_final'

dataHCPD <- readRDS(paste0(interfileFolderHCPD, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
dataABCD <- readRDS(paste0(interfileFolderABCD, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
dataABCD$age <- dataABCD$age / 12
dataChineseCohort <- readRDS(paste0(interfileFolderChineseCohort, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatAge.rds"))
dataHBN <- readRDS(paste0(interfileFolderHBN, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
Behavior_HBN <- read.csv("D:/xuxiaoyu/open_dataset_information/HBN/demography/basic_demo_merge_screen_MRIQC.csv")

# HCPD
ggplot(data = dataHCPD, aes(age, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", boundary = 5, position = position_dodge(width = 0.8)) +
  labs(x = "Age (years)", y = "Frequency", title = paste0("HCP-D, N=", nrow(dataHCPD))) +
  scale_x_continuous(limits = c(8, 23), breaks = c(8,11,14,17,20,23)) +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1.2,
    plot.title = element_text(color = "black", size = 18, hjust = 0.5),
    axis.title = element_text(color = "black", size = 18),axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(color = "black", size = 18),
    legend.position = "none")

ggsave(paste(FigureFolderHCPD, '/Fig1_Age_distribution_count_HCPD.tiff', sep = ''), width = 11, height = 13, units = "cm")
ggsave(paste(FigureFolderHCPD, '/Fig1_Age_distribution_count_HCPD.svg', sep = ''), width = 12, height = 12, units = "cm")

dataHCPD <- dataHCPD %>% mutate(age_group = case_when(
  floor(age) < ceiling(age) ~ paste0(floor(age), "_", ceiling(age)),
  floor(age) == ceiling(age) ~ paste0(floor(age), "_", (ceiling(age)+1))
))
dataHCPD$age_group <- as.factor(dataHCPD$age_group)
countstats_HCPD <- as.data.frame(table(dataHCPD$age_group))
names(countstats_HCPD) <- c("age_group", "n")
write.csv(countstats_HCPD, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_HCPD.csv", row.names=F)

countstats_HCPD_F <- countstats_HCPD %>% filter(sex == "F")
countstats_HCPD_M <- countstats_HCPD %>% filter(sex == "M")
countstats_HCPD_F$sex <- NULL
countstats_HCPD_M$sex <- NULL
write.csv(countstats_HCPD_F, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_HCPD_F.csv", row.names=F)
write.csv(countstats_HCPD_M, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_HCPD_M.csv", row.names=F)

# ABCD
fillcolor = c("#83B7D7", "#B4D3E7")
ggplot(data = dataABCD, aes(age, y = ..count.., fill = eventname)) +
  geom_histogram(binwidth = 0.5, color = "black", boundary = 5, position = "stack",linewidth=0.5) +
  labs(x = "Age (years)", y = NULL, title = "ABCD, N=7,104") +
  #scale_x_continuous(limits = c(9, 16), breaks = c(9,10,11,12,13,14,15,16)) +
  #scale_y_continuous(limits = c(0, 50), breaks = c(0,10, 20,30,40, 50)) +
  scale_fill_manual(values = fillcolor)+
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1.2,
        plot.title = element_text(color = "black", size = 18, hjust = 0.5),
        axis.title = element_text(color = "black", size = 18),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 18),
        legend.position = "none")

ggsave(paste(FigureFolderABCD, '/Fig1_Age_distribution_count_ABCD.tiff', sep = ''),  width = 11, height = 13, units = "cm")
ggsave(paste(FigureFolderABCD, '/Fig1_Age_distribution_count_ABCD.svg', sep = ''), width = 12, height = 10, units = "cm")

length(which(table(dataABCD$subID)==2)) # 2570 participants have 2 visits

dataABCD <- dataABCD %>% mutate(age_group = 
                                  cut(age, breaks = seq(floor(min(age)), ceiling(max(age)), by = 0.5), right = FALSE)
)

countstats_ABCD <- as.data.frame(table(dataABCD$age_group, dataABCD$eventname))
names(countstats_ABCD) <- c("age_group", "eventname", "n")
countstats_ABCD <- countstats_ABCD %>% pivot_wider(names_from = eventname,
                                                   values_from = n,
                                                   values_fill = list(n = 0))

write.csv(countstats_ABCD, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_ABCD.csv", row.names=F)

countstats_ABCD_F <- countstats_ABCD %>% filter(sex == 2)
countstats_ABCD_M <- countstats_ABCD %>% filter(sex == 1)
countstats_ABCD_F$sex <- NULL
countstats_ABCD_M$sex <- NULL
names(countstats_ABCD_F)[2:3] <- c("2-Year Follow-up", "Baseline")
names(countstats_ABCD_M)[2:3] <- c("2-Year Follow-up", "Baseline")
write.csv(countstats_ABCD_F, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_ABCD_F.csv", row.names=F)
write.csv(countstats_ABCD_M, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_ABCD_M.csv", row.names=F)

# Chinese Cohort
dataChineseCohort$study[dataChineseCohort$study=="CCNP"] <- paste0(dataChineseCohort$study[dataChineseCohort$study=="CCNP"], "_", dataChineseCohort$ses[dataChineseCohort$study=="CCNP"])

fillcolor <- c("#F7FCF5", "#E5F5E0","#C7E9C0", "#B4D3E7", "#83B7D7")
ggplot(data = dataChineseCohort, aes(Age, y = ..count.., fill = study)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 5, position = "stack",linewidth=0.5) +
  scale_fill_manual(values = fillcolor)+
  labs(x = "Age (years)", y = "Frequency", title = paste0("Chinese Cohort, N=", nrow(dataChineseCohort))) +
  #scale_x_continuous(limits = c(8, 23), breaks = c(8,11,14,17,20,23)) +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  labs(y=NULL)+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1.2,
        plot.title = element_text(color = "black", size = 18, hjust = 0.5),
        axis.title = element_text(color = "black", size = 18),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 18),
        legend.position = "none")

ggsave(paste(FigureFolderChineseCohort, '/Fig1_Age_distribution_count_ChineseCohort.tiff', sep = ''), width = 11, height = 13, units = "cm")
ggsave(paste(FigureFolderChineseCohort, '/Fig1_Age_distribution_count_ChineseCohort.svg', sep = ''), width = 12, height = 10, units = "cm")

dataChineseCohort <- dataChineseCohort %>% mutate(age_group = case_when(
  floor(Age) < ceiling(Age) ~ paste0(floor(Age), "_", ceiling(Age)),
  floor(Age) == ceiling(Age) ~ paste0(floor(Age), "_", (ceiling(Age)+1))
))

countstats_ChineseCohort <- as.data.frame(table(dataChineseCohort$age_group, dataChineseCohort$study))
names(countstats_ChineseCohort) <- c("age_group", "study", "n")
countstats_ChineseCohort <- countstats_ChineseCohort %>% pivot_wider(names_from = study,
                                                   values_from = n,
                                                   values_fill = list(n = 0))

write.csv(countstats_ChineseCohort, "D:/xuxiaoyu/DMRI_network_development/SC_development/website/demographydata/countstats_ChineseCohort.csv", row.names=F)



## Description for demographic information
dataHCPD <- within(dataHCPD, {
  nih_fluidcogcomp_unadjusted[nih_fluidcogcomp_unadjusted<mean(nih_fluidcogcomp_unadjusted, na.rm=T)-3*sd(nih_fluidcogcomp_unadjusted, na.rm=T) | nih_fluidcogcomp_unadjusted>mean(nih_fluidcogcomp_unadjusted, na.rm=T)+3*sd(nih_fluidcogcomp_unadjusted, na.rm=T)] <- NA
})
summary(dataHCPD$nih_fluidcogcomp_unadjusted)
HCPD_vars <- c("age", "sex", "handnessfactor", "race_ethnicity", "mean_fd", "nih_fluidcogcomp_unadjusted", "site", "ICV", "income.adj")
HCPD_tableone <- CreateTableOne(HCPD_vars, data = dataHCPD, factorVars=c("sex", "handnessfactor", "race_ethnicity", "site"))
HCPD_tableone <- as.data.frame(print(HCPD_tableone))
write.csv(HCPD_tableone, paste0(FigureFolderHCPD, '/demographic_info.csv'))

ABCD_vars <- c("age", "sex","eventname", "handness", "race_ethnicity", "mean_fd", "nihtbx_fluidcomp_uncorrected","GENERAL", "siteID", "smri_vol_scs_intracranialv", "income.adj")
ABCD_tableone <- CreateTableOne(ABCD_vars,strata="eventname", data = dataABCD, factorVars=c("sex", "handness", "race_ethnicity", "siteID"), test=F)
ABCD_tableone <- as.data.frame(print(ABCD_tableone))
write.csv(ABCD_tableone, paste0(FigureFolderABCD, '/demographic_info.csv'))

ChineseCohort_vars <- c("Age", "Sex", "Handedness", "mean_fd", "ICV", "CBCLtotalproblem", "FluidComp")
#dataChineseCohort$study <- str_split_i(dataChineseCohort$study, "_", 1)
ChineseCohort_tableone <- CreateTableOne(ChineseCohort_vars,strata="study", data = dataChineseCohort, factorVars=c("Sex", "Handedness"), test=F)
ChineseCohort_tableone <- as.data.frame(print(ChineseCohort_tableone))
write.csv(ChineseCohort_tableone, paste0(FigureFolderChineseCohort, '/demographic_info.csv'))

