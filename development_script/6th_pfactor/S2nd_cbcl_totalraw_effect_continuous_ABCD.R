Sys.setenv(R_LIBS_USER = "", R_LIBS = "")
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_prefix)) {
  .libPaths(file.path(conda_prefix, "lib", "R", "library"))
}

library(mgcv)
library(parallel)
rm(list = ls())
CVthr <- 75

args <- commandArgs(trailingOnly = TRUE)
test_n <- if (length(args) >= 1) as.integer(args[[1]]) else 0
test_edges <- if (length(args) >= 2) as.integer(args[[2]]) else 0

wdpath <- getwd()
if (grepl("cuizaixu_lab", wdpath)) {
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
SCdata[, c("sex", "handness", "race_ethnicity")] <- lapply(SCdata[, c("sex", "handness", "race_ethnicity")], as.factor)
SCdata$age <- SCdata$age / 12
SCdata$cbcl_scr_syn_totprob_r <- as.numeric(SCdata$cbcl_scr_syn_totprob_r)
SCdata$totalstrength <- rowMeans(SCdata[, grepl("^SC\\.", names(SCdata)) & grepl("_h$", names(SCdata))])

if (test_n > 0) {
  SCdata <- do.call(rbind, lapply(split(SCdata, SCdata$siteID), function(df) {
    utils::head(df, test_n)
  }))
}

meandistance <- read.csv(paste0(interfileFolder, "/average_EuclideanDistance_12.csv"))
meandistance <- meandistance$Edistance

## 1. CBCL total raw association
dataname <- "SCdata"
smooth_var <- "age"
int_var.predict.percentile <- 0.1
covariates <- "sex+mean_fd"
knots <- 3
set_fx <- TRUE
increments <- 1000
stats_only <- if (test_n > 0 || test_edges > 0) FALSE else TRUE
int_var <- "cbcl_scr_syn_totprob_r"

edge_count <- if (test_edges > 0) min(78, test_edges) else 78
mc_cores <- if (test_n > 0 || test_edges > 0) 1 else 50
resultsum <- mclapply(1:edge_count, function(x) {
  region <- grep("^SC\\.", names(SCdata), value = TRUE)[x]
  gamresult <- gamm.smooth.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                        int_var.predict.percentile, covariates,
                                                        knots, set_fx, increments, stats_only)
  gamresult <- as.data.frame(gamresult)
  return(gamresult)
}, mc.cores = mc_cores)

gamresult.tmp <- do.call(rbind, resultsum)
gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

print(paste0(sum(gamresult.tmp$bootstrap_pvalue.fdr < 0.05), " edges have significant age by ", int_var, " effect."))
print(paste0(sum(gamresult.tmp$bootstrap.P.disease.fdr < 0.05), " edges have significant ", int_var, " effect."))

saveRDS(gamresult.tmp, file.path(outputFolder, paste0("gamresult_Int_age_cbcl_totprob_raw_CV", CVthr, ".rds")))

if (test_n > 0 || test_edges > 0) {
  message("Test mode: skip S-A correlation and distance-controlled correlation.")
  quit(save = "no", status = 0)
}

## 2. Correlation to S-A connectional axis
gamresult.tmp <- readRDS(file.path(outputFolder, paste0("gamresult_Int_age_cbcl_totprob_raw_CV", CVthr, ".rds")))
gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = FALSE)
SCrank.df.cbcl <- SCrankcorr(gamresult.tmp, "T.disease", 12, dsdata = FALSE)
SCrank.df <- rbind(SCrank.df.age, SCrank.df.cbcl)
SCrank.df$int_var <- int_var
SCrank.df

SCrank.tmp <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = TRUE)
SCrank <- SCrank.tmp$SCrank

print("Next, correlation between CBCL associations and connectional axis is tested while controlling for Euclidean distance.")
gamresult.tmp$meandistance <- meandistance
print(stats::cor.test(gamresult.tmp$T.disease, gamresult.tmp$meandistance, method = "pearson"))
gamresult.tmp$T.disease_control_distance[which(!is.na(gamresult.tmp$T.disease))] <- residuals(lm(T.disease ~ meandistance, data = gamresult.tmp))
print(stats::cor.test(gamresult.tmp$T.disease_control_distance, gamresult.tmp$meandistance, method = "pearson"))
SCrankresult.whole.controldistance <- SCrankcorr(gamresult.tmp, "T.disease_control_distance", 12, dsdata = FALSE)
print(paste("Correlation coefficient between CBCL associations regressing out fiber distance and connectional axis is",
            round(SCrankresult.whole.controldistance$r.spearman, 2), "with a P value of",
            round(SCrankresult.whole.controldistance$p.spearman, 3)))
print(SCrankresult.whole.controldistance)

saveRDS(SCrankresult.whole.controldistance,
        file.path(outputFolder, paste0("SCrankcorr_cbcl_totprob_raw_CV", CVthr, ".rds")))
