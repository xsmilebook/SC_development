#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(tidyverse)
})

rm(list = ls())

CVthr <- 75
int_var <- "GENERAL"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "pfactor")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
}

euclid_csv <- Sys.getenv(
  "ABCD_EUCLID_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "average_EuclideanDistance_12.csv")
)
if (!file.exists(euclid_csv)) {
  stop("Missing ABCD_EUCLID_CSV: ", euclid_csv)
}
meandistance <- read.csv(euclid_csv)$Edistance

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "gamminteraction.R"))

SCdata <- readRDS(input_rds)
SCdata$age <- as.numeric(SCdata$age) / 12

needed_base <- c("subID", "scanID", "siteID", "age", "sex", "mean_fd", int_var)
missing_base <- setdiff(needed_base, names(SCdata))
if (length(missing_base) > 0) {
  stop("Missing required columns: ", paste(missing_base, collapse = ", "))
}

SCdata[, c("sex")] <- lapply(SCdata[, c("sex"), drop = FALSE], as.factor)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
sc_cols <- sc_cols[str_detect(sc_cols, "_h$")]
if (length(sc_cols) < 78) {
  stop("Expected >=78 ComBat SC.*_h columns, got: ", length(sc_cols))
}

dataname <- "SCdata"
smooth_var <- "age"
covariates <- "sex+mean_fd"
knots <- 3
set_fx <- TRUE
increments <- 1000
int_var_predict_percentile <- 0.1
stats_only <- TRUE

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "50"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 50
num_cores <- min(num_cores, 50L)

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("gamresult_Int_age_pFactor_", int_var, "_CV", CVthr, ".rds"))

if (force || !file.exists(out_rds)) {
  message("[INFO] Fitting p-factor GAMM interaction models (n_edges=78, mc.cores=", num_cores, ")")
  resultsum <- parallel::mclapply(seq_len(78), function(x) {
    region <- sc_cols[[x]]
    gamresult <- gamm.smooth.predict.covariateinteraction(
      region,
      dataname,
      smooth_var,
      int_var,
      int_var_predict_percentile,
      covariates,
      knots,
      set_fx,
      increments,
      stats_only = stats_only
    )
    as.data.frame(gamresult)
  }, mc.cores = num_cores)

  gamresult.tmp <- do.call(rbind, resultsum)
  gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
  gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
  gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

  message(sum(gamresult.tmp$bootstrap_pvalue.fdr < 0.05), " edges have significant age by ", int_var, " effect.")
  message(sum(gamresult.tmp$bootstrap.P.disease.fdr < 0.05), " edges have significant ", int_var, " effect.")

  saveRDS(gamresult.tmp, out_rds)
} else {
  gamresult.tmp <- readRDS(out_rds)
}

gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

message("[INFO] Correlation to connectional axis")
SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = FALSE)
SCrank.df.general <- SCrankcorr(gamresult.tmp, "T.disease", 12, dsdata = FALSE)
SCrank.df <- rbind(SCrank.df.age, SCrank.df.general)
SCrank.df$int_var <- int_var

message("[INFO] Control Euclidean distance for T.disease")
gamresult.tmp$meandistance <- meandistance
gamresult.tmp$T.disease_control_distance[which(!is.na(gamresult.tmp$T.disease))] <-
  residuals(lm(T.disease ~ meandistance, data = gamresult.tmp))
SCrank.df.general.controldistance <- SCrankcorr(gamresult.tmp, "T.disease_control_distance", 12, dsdata = FALSE)

saveRDS(
  list(
    SCrank.df = SCrank.df,
    SCrank.df.general.controldistance = SCrank.df.general.controldistance
  ),
  file.path(resultFolder, paste0("SCrankcorr_summary_Int_age_pFactor_", int_var, "_CV", CVthr, ".rds"))
)

