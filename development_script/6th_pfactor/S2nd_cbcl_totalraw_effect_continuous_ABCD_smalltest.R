library(mgcv)
library(parallel)
library(tidyverse)
library(reshape)
library(RColorBrewer)
rm(list = ls())

CVthr <- 75
set.seed(42)

wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")) {
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD"
} else {
  interfileFolder <- file.path(wdpath, "interdataFolder_ABCD")
}
functionFolder <- file.path(wdpath, "gamfunction")
outputFolder <- file.path(wdpath, "outputs", "results", "cbcl_totprob")
combatFolder <- file.path(wdpath, "outputs", "results", "combat_cbcl")

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

source(paste0(functionFolder, "/gamminteraction.R"))
source(paste0(functionFolder, "/SCrankcorr.R"))

input_rds <- file.path(combatFolder, "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatCBCLtotalraw.rds")
SCdata <- readRDS(input_rds)
if (is.data.frame(SCdata$age) || !is.numeric(SCdata$age)) {
  demodf <- read.csv(file.path(wdpath, "demopath", "DemodfScreenFinal.csv"))
  SCdata$age <- demodf$age[match(SCdata$scanID, demodf$scanID)]
}
SCdata[, c("sex", "handness", "race_ethnicity")] <- lapply(SCdata[, c("sex", "handness", "race_ethnicity")], as.factor)
SCdata$age <- SCdata$age / 12
SCdata$cbcl_scr_syn_totprob_r <- as.numeric(SCdata$cbcl_scr_syn_totprob_r)
SCdata$totalstrength <- rowMeans(SCdata[, str_detect(names(SCdata), "SC\\.") & str_detect(names(SCdata), "_h")])

meandistance <- read.csv(paste0(interfileFolder, "/average_EuclideanDistance_12.csv"))
meandistance <- meandistance$Edistance

SCdata <- SCdata[sample(seq_len(nrow(SCdata)), min(300, nrow(SCdata))), , drop = FALSE]
sc_edges <- grep("SC\\.", names(SCdata), value = TRUE)
sc_edges <- sc_edges[seq_len(min(3, length(sc_edges)))]

## 1. CBCL total raw association
smooth_var <- "age"
int_var.predict.percentile <- 0.1
covariates <- "sex+mean_fd"
knots <- 3
set_fx <- TRUE
increments <- 200
stats_only <- FALSE
int_var <- "cbcl_scr_syn_totprob_r"

resultsum <- mclapply(seq_along(sc_edges), function(x) {
  region <- sc_edges[x]
  gamresult <- gamm.smooth.predict.covariateinteraction(region, "SCdata", smooth_var, int_var,
                                                        int_var.predict.percentile, covariates,
                                                        knots, set_fx, increments, stats_only)
  gamresult <- as.data.frame(gamresult)
  return(gamresult)
}, mc.cores = 2)

gamresult.tmp <- do.call(rbind, resultsum)
gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

saveRDS(gamresult.tmp, file.path(outputFolder, paste0("gamresult_Int_age_cbcl_totprob_raw_CV", CVthr, "_smalltest.rds")))

if (length(sc_edges) < 78) {
  quit(save = "no")
}

## 2. Correlation to S-A connectional axis
if (nrow(gamresult.tmp) >= 78) {
  SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = FALSE)
  SCrank.df.cbcl <- SCrankcorr(gamresult.tmp, "T.disease", 12, dsdata = FALSE)
  SCrank.df <- rbind(SCrank.df.age, SCrank.df.cbcl)
  SCrank.df$int_var <- int_var
}

if (nrow(gamresult.tmp) >= 78) {
  SCrank.tmp <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = TRUE)
  SCrank <- SCrank.tmp$SCrank
}

if (nrow(gamresult.tmp) >= 78) {
  gamresult.tmp$meandistance <- meandistance[seq_len(nrow(gamresult.tmp))]
  SCrankresult.whole.controldistance <- SCrankcorr(gamresult.tmp, "T.disease", 12, dsdata = FALSE)
}

if (nrow(gamresult.tmp) >= 78) {
  saveRDS(SCrankresult.whole.controldistance,
          file.path(outputFolder, paste0("SCrankcorr_cbcl_totprob_raw_CV", CVthr, "_smalltest.rds")))
}
