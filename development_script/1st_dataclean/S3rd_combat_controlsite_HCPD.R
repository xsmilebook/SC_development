library(R.matlab)
library(psych)
library(mgcv)
library(tidyverse)
rm(list = ls())

functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results'
demopath<-'D:/xuxiaoyu/open_dataset_information/HCP/HCPD_behavior'
source(paste0(functionFolder, '/combat.R'))
CVthr = 75
resolutionds <- 12
edgenum <- (resolutionds+1)*resolutionds / 2
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', resolutionds,'_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
Behavior <- read.csv(paste0(demopath, '/HCPD_demo_behav.csv'))
Behavior$subID <- paste0("sub-", Behavior$subID)
####combat for SC strength###########################################################
#################################################################################
SC_vars <- grep("SC.", names(SCdata), value=T)
model_terms <- c(SC_vars,"subID", "age", "site","sex", "mean_fd")
comtable <- SCdata %>% select(model_terms) %>%
  drop_na()
sitetab <- table(comtable$site)
removesite <- names(sitetab)[which(sitetab < 100)]
comtable <- comtable[which(! comtable$site %in% removesite),]

batch <- as.character(comtable$site)
harmonized_data_age <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data_age) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data_age$subID=comtable$subID
for (i in 1:edgenum){
  ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
  smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+mean_fd"
  region <- SC_vars[i]
  # age model
  modelformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula1, method="REML", data = comtable)
  mod<- model.matrix(gam.model)
  # combat. adjust site effect
  combatdata <- combat(ctab,batch,mod=mod, eb=F,verbose=TRUE)
  harmonized_data_age[,i]<-t(combatdata$dat.combat)
}
dataTable<- merge(harmonized_data_age, Behavior, by="subID")
describe(comtable$SC.1)
describe(dataTable$SC.1_h)
corr.test(comtable$SC.78, dataTable$SC.78_h)
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds,"_CV", CVthr,"_sumSCinvnode.sum.msmtcsd.combatage.rds"))

# cognition model
model_terms <- c(SC_vars,"subID", "age", "site","sex", "mean_fd", "nih_fluidcogcomp_unadjusted")
comtable <- SCdata %>% select(model_terms) %>%
  drop_na() %>% filter(nih_fluidcogcomp_unadjusted < mean(nih_fluidcogcomp_unadjusted)+3*sd(nih_fluidcogcomp_unadjusted) & nih_fluidcogcomp_unadjusted > mean(nih_fluidcogcomp_unadjusted)-3*sd(nih_fluidcogcomp_unadjusted))
sitetab <- table(comtable$site)

batch <- as.character(comtable$site)
harmonized_data_cognition <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data_cognition) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data_cognition$subID=comtable$subID
for (i in 1:edgenum){
  ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
  smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+mean_fd"
  region <- SC_vars[i]; Cogvar="nih_fluidcogcomp_unadjusted"
  # age model
  modelformula2 <- as.formula(sprintf("%s ~ %s + s(%s, k = %s, fx = %s) + %s",region, Cogvar, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula2, method="REML", data = comtable)
  mod<- model.matrix(gam.model)
  # combat. adjust site effect
  combatdata <- combat(ctab,batch,mod=mod, eb=F,verbose=TRUE)
  harmonized_data_cognition[,i]<-t(combatdata$dat.combat)
}
dataTable<- merge(harmonized_data_cognition, Behavior, by="subID")
describe(comtable$SC.1)
describe(dataTable$SC.1_h)
corr.test(comtable$SC.1, dataTable$SC.1_h)
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds,"_CV", CVthr,"_sumSCinvnode.sum.msmtcsd.combatNTBfluid_cog2SC.rds"))

