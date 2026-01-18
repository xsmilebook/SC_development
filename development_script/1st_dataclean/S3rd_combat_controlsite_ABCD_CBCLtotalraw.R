library(R.matlab)
library(psych)
library(mgcv)
library(tidyverse)
library(lme4)
library(gamm4)
rm(list = ls())

wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")) {
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD"
  demopath <- file.path(wdpath, "demopath")
  functionFolder <- file.path(wdpath, "gamfunction")
} else {
  interfileFolder <- "D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD"
  demopath <- "D:/xuxiaoyu/open_dataset_information/ABCD/info"
  functionFolder <- "D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction"
}

source(paste0(functionFolder, "/combat.R"))

resolutionds <- 12
edgenum <- (resolutionds + 1) * resolutionds / 2
SCdata <- readRDS(paste0(interfileFolder, "/SCdata_SA", resolutionds, "_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"))
Behavior <- read.csv(file.path(demopath, "DemodfScreenFinal.csv"))

SCdata$ACSSCORE <- NULL
SC_vars <- grep("SC\\.", names(SCdata), value = TRUE)

# CBCL total raw model (sensitivity analysis for p-factor)
model_terms <- c(SC_vars, "subID", "age", "siteID", "sex", "mean_fd",
                 "scanID", "cbcl_scr_syn_totprob_r")
comtable <- SCdata %>% select(model_terms) %>%
  drop_na()
sitetab <- table(comtable$siteID)

batch <- as.character(comtable$siteID)
harmonized_data_cbcl <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data_cbcl) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data_cbcl$scanID <- comtable$scanID

for (i in 1:edgenum) {
  ctab <- t(data.matrix(comtable %>% select(SC_vars[i])))
  smooth_var <- "age"
  knots <- 3
  set_fx <- TRUE
  covariates <- "sex+mean_fd"
  cbcl_var <- "cbcl_scr_syn_totprob_r"
  region <- SC_vars[i]
  modelformula <- as.formula(sprintf("%s ~ %s + s(%s, k = %s, fx = %s) + %s",
                                     region, cbcl_var, smooth_var, knots, set_fx, covariates))
  gamm.model <- gamm4(modelformula, random = ~(1 | subID), REML = TRUE, data = comtable)
  mod <- model.matrix(gamm.model$gam)
  combatdata <- combat(ctab, batch, mod = mod, eb = FALSE, verbose = TRUE)
  harmonized_data_cbcl[, i] <- t(combatdata$dat.combat)
}

dataTable <- merge(harmonized_data_cbcl, Behavior, by = "scanID")
describe(comtable$SC.50)
describe(dataTable$SC.50_h)
corr.test(comtable$SC.50, dataTable$SC.50_h)
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds,
                          "_CV75_sumSCinvnode.sum.msmtcsd.combatCBCLtotalraw.rds"))
