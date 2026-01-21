#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
output_rds <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_comp_agecorrected_baseline.rds"
test_n <- if (length(args) >= 3) as.integer(args[[3]]) else 0

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

source(file.path("combat_gam", "longitudinal", "neuroCombat.R"))
source(file.path("combat_gam", "longitudinal", "nonlinearlongcombat.R"))

set.seed(42)
n_cores_env <- Sys.getenv("SLURM_CPUS_PER_TASK")
n_cores <- if (nzchar(n_cores_env)) as.integer(n_cores_env) else 1
if (is.na(n_cores) || n_cores < 1) {
  n_cores <- 1
}

scdata <- readRDS(input_rds)
sc_cols <- grep("^SC\\.", names(scdata), value = TRUE)
if (length(sc_cols) == 0) {
  stop("No SC.* columns found in input.")
}

needed <- c(sc_cols, "subID", "scanID", "siteID", "age", "sex", "mean_fd", "eventname", "nihtbx_totalcomp_agecorrected")
missing <- setdiff(needed, names(scdata))
if (length(missing) > 0) {
  stop(paste("Missing columns:", paste(missing, collapse = ", ")))
}

df <- scdata[, needed]
df <- df %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  drop_na()

if (test_n > 0) {
  df <- df %>%
    group_by(siteID) %>%
    filter(n() >= 2) %>%
    slice_head(n = min(2, n())) %>%
    ungroup()
}

normalize_numeric <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.numeric(factor(x)))
  }
  as.numeric(x)
}

features <- as.matrix(df[, sc_cols])
ranef <- df$subID
X_nlin <- as.matrix(df$age)
colnames(X_nlin) <- "age"
X_lin <- df[, c("sex", "mean_fd", "nihtbx_totalcomp_agecorrected"), drop = FALSE]
X_lin <- as.data.frame(lapply(X_lin, normalize_numeric))
colnames(X_lin) <- c("sex", "mean_fd", "nihtbx_totalcomp_agecorrected")
batch <- df$siteID

use_random <- length(unique(ranef)) < nrow(df)
result <- nlongcombat(
  features,
  knot = 3,
  ranef = ranef,
  X_nlin = X_nlin,
  X_lin = X_lin,
  batch = batch,
  use_random = use_random,
  n_cores = n_cores
)

harmonized <- as.data.frame(result$harmonized)
colnames(harmonized) <- paste0(sc_cols, "_h")

out <- data.frame(scanID = df$scanID)
out$subID <- df$subID
out$siteID <- df$siteID
out$age <- df$age
out$sex <- df$sex
out$mean_fd <- df$mean_fd
out$nihtbx_totalcomp_agecorrected <- df$nihtbx_totalcomp_agecorrected
out <- cbind(out, harmonized)

out_dir <- dirname(output_rds)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
saveRDS(out, output_rds)

