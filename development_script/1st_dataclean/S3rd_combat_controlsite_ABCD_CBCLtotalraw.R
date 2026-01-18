library(R.matlab)
library(psych)
library(mgcv)
library(lme4)
library(gamm4)
rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
test_n <- if (length(args) >= 1) as.integer(args[[1]]) else 0
test_edges <- if (length(args) >= 2) as.integer(args[[2]]) else 0

wdpath <- getwd()
if (grepl("cuizaixu_lab", wdpath)) {
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
Behavior_cbcl <- Behavior[, c("scanID", "cbcl_scr_syn_totprob_r")]
SCdata <- merge(SCdata, Behavior_cbcl, by = "scanID", all.x = TRUE)
SC_vars <- grep("SC\\.", names(SCdata), value = TRUE)
edgenum_run <- if (test_edges > 0) min(edgenum, test_edges) else edgenum
SC_vars_run <- SC_vars[seq_len(edgenum_run)]

# CBCL total raw model (sensitivity analysis for p-factor)
model_terms <- c(SC_vars_run, "subID", "age", "siteID", "sex", "mean_fd",
                 "scanID", "cbcl_scr_syn_totprob_r")
comtable <- SCdata[, model_terms, drop = FALSE]
comtable <- comtable[complete.cases(comtable), , drop = FALSE]
if (test_n > 0) {
  comtable <- do.call(rbind, lapply(split(comtable, comtable$siteID), function(df) {
    utils::head(df, test_n)
  }))
}
sitetab <- table(comtable$siteID)

batch <- as.character(comtable$siteID)
harmonized_data_cbcl <- data.frame(matrix(NA, nrow(comtable), edgenum_run))
names(harmonized_data_cbcl) <- paste0("SC.", c(1:edgenum_run), "_h")
harmonized_data_cbcl$scanID <- comtable$scanID

for (i in 1:edgenum_run) {
  ctab <- t(data.matrix(comtable[, SC_vars_run[i], drop = FALSE]))
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
qc_var <- if ("SC.50" %in% names(comtable)) "SC.50" else SC_vars_run[1]
qc_var_h <- paste0(qc_var, "_h")
describe(comtable[[qc_var]])
if (qc_var_h %in% names(dataTable)) {
  describe(dataTable[[qc_var_h]])
  corr.test(comtable[[qc_var]], dataTable[[qc_var_h]])
}
suffix <- if (test_n > 0 || test_edges > 0) {
  paste0(".test_n", test_n, "_edges", edgenum_run)
} else {
  ""
}
saveRDS(dataTable, paste0(interfileFolder, "/SCdata_SA", resolutionds,
                          "_CV75_sumSCinvnode.sum.msmtcsd.combatCBCLtotalraw", suffix, ".rds"))
