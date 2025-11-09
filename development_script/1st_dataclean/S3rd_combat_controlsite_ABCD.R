library(R.matlab);
library(psych)
library(mgcv);
library(tidyverse)
library(lme4)
library(gamm4)
rm(list = ls())

functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results'
demopath<-"D:/xuxiaoyu/open_dataset_information/ABCD/info"

source(paste0(functionFolder, '/combat.R'))
resolutionds <- 12
edgenum <- (resolutionds+1)*resolutionds / 2
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', resolutionds, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
Behavior <- read.csv("D:/xuxiaoyu/open_dataset_information/ABCD/info/DemodfScreenFinal.csv")
SCdata$ACSSCORE <- NULL
SC_vars <- grep("SC.", names(SCdata), value=T)

####combat for SC strength###########################################################
#################################################################################
model_terms <- c(SC_vars,"subID", "age", "siteID","sex", "mean_fd", 
                 "scanID")
comtable <- SCdata %>% select(model_terms) %>%
  drop_na() 
sitetab <- table(comtable$siteID) #13 sites

batch <- as.character(comtable$siteID)
harmonized_data_age <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data_age) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data_age$scanID=comtable$scanID
for (i in 1:edgenum){
  ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
  smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+mean_fd"
  region <- SC_vars[i]
  # age model
  modelformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gamm.model <- gamm4(modelformula1, random=~(1|subID), REML=TRUE, data = comtable)
  mod<- model.matrix(gamm.model$gam)
  # combat. adjust site effect
  combatdata <- combat(ctab,batch,mod=mod, eb=F,verbose=TRUE)
  harmonized_data_age[,i]<-t(combatdata$dat.combat)
}
dataTable<- merge(harmonized_data_age, Behavior, by="scanID")
describe(comtable$SC.1)
describe(dataTable$SC.1_h)
corr.test(comtable$SC.1, dataTable$SC.1_h)
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds, "_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))

# pFactor-GENERAL model
model_terms <- c(SC_vars,"subID", "age", "siteID", "sex", "mean_fd", 
                 "scanID", "GENERAL")
comtable <- SCdata %>% select(model_terms) %>%
  drop_na() 
batch <- as.character(comtable$siteID)
harmonized_data_age <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data_age) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data_age$scanID=comtable$scanID
for (i in 1:edgenum){
  ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
  smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+mean_fd"
  pFactor_var <- "GENERAL"
  region <- SC_vars[i]
  # age model
  modelformula1 <- as.formula(sprintf("%s ~ %s+s(%s, k = %s, fx = %s) + %s", region,pFactor_var, smooth_var, knots, set_fx, covariates))
  gamm.model <- gamm4(modelformula1, random=~(1|subID), REML=TRUE, data = comtable)
  mod<- model.matrix(gamm.model$gam)
  # combat. adjust site effect
  combatdata <- combat(ctab,batch,mod=mod, eb=F,verbose=TRUE)
  harmonized_data_age[,i]<-t(combatdata$dat.combat)
}
dataTable<- merge(harmonized_data_age, Behavior, by="scanID")
describe(comtable$SC.50)
describe(dataTable$SC.50_h)
corr.test(comtable$SC.50, dataTable$SC.50_h)
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds, "_CV75_sumSCinvnode.sum.msmtcsd.combatPFactorGeneral.rds"))


# cognition model
model_terms <- c(SC_vars,"subID", "age", "siteID", "sex", "mean_fd", 
                 "scanID", "nihtbx_fluidcomp_uncorrected", "eventname")
comtable <- SCdata %>% select(model_terms) %>%
  drop_na() %>% 
  filter(eventname=="baseline_year_1_arm_1")
sitetab <- table(comtable$siteID)

batch <- as.character(comtable$siteID)
harmonized_data_cognition <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data_cognition) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data_cognition$scanID=comtable$scanID
for (i in 1:edgenum){
  ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
  smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+mean_fd"
  region <- SC_vars[i]; Cogvar="nihtbx_fluidcomp_uncorrected"
  # age model
  modelformula2 <- as.formula(sprintf("%s ~ %s + s(%s, k = %s, fx = %s) + %s",region, Cogvar, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula2, method="REML", data = comtable)
  mod<- model.matrix(gam.model)
  # combat. adjust site effect
  combatdata <- combat(ctab,batch,mod=mod, eb=F,verbose=TRUE)
  harmonized_data_cognition[,i]<-t(combatdata$dat.combat)
}
dataTable<- merge(harmonized_data_cognition, Behavior, by="scanID")
describe(comtable$SC.5)
describe(dataTable$SC.5_h)
corr.test(comtable$SC.50, dataTable$SC.50_h)
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds, "_CV75_sumSCinvnode.sum.msmtcsd.combatNTBfluid_cog2SC.rds"))

