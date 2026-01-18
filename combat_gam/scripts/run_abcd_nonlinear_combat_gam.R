#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_abcd_nonlinear_combat_gam.R <input_rds> <output_dir> [test_n]")
}

input_rds <- args[[1]]
output_dir <- args[[2]]
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

base_cols <- c(
  "subID", "scanID", "siteID", "age", "sex", "mean_fd",
  "nihtbx_fluidcomp_uncorrected", "GENERAL", "eventname"
)

keep_cols <- unique(c(sc_cols, base_cols))
missing <- setdiff(keep_cols, names(scdata))
if (length(missing) > 0) {
  stop(paste("Missing columns:", paste(missing, collapse = ", ")))
}

scdata <- scdata[, keep_cols]

normalize_numeric <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.numeric(factor(x)))
  }
  as.numeric(x)
}

run_variant <- function(df, variant, extra_covars, baseline_only = FALSE) {
  if (baseline_only) {
    df <- df %>% filter(eventname == "baseline_year_1_arm_1")
  }
  needed <- c(sc_cols, "subID", "scanID", "siteID", "age", "sex", "mean_fd", extra_covars)
  df <- df[, needed]
  df <- df %>% tidyr::drop_na()
  if (test_n > 0) {
    df <- df %>%
      group_by(siteID) %>%
      filter(n() >= 2) %>%
      slice_head(n = 2) %>%
      ungroup()
  }

  features <- as.matrix(df[, sc_cols])
  ranef <- df$subID
  X_nlin <- as.matrix(df$age)
  colnames(X_nlin) <- "age"
  X_lin <- df[, c("sex", "mean_fd", extra_covars), drop = FALSE]
  X_lin <- as.data.frame(lapply(X_lin, normalize_numeric))
  colnames(X_lin) <- c("sex", "mean_fd", extra_covars)
  batch <- df$siteID

  use_random <- length(unique(ranef)) < nrow(df)
  result <- nlongcombat(features, knot = 3, ranef = ranef, X_nlin = X_nlin, X_lin = X_lin, batch = batch, use_random = use_random, n_cores = n_cores)
  harmonized <- as.data.frame(result$harmonized)
  colnames(harmonized) <- paste0(sc_cols, "_h")

  out <- data.frame(scanID = df$scanID)
  out$subID <- df$subID
  out$siteID <- df$siteID
  out$age <- df$age
  out$sex <- df$sex
  out$mean_fd <- df$mean_fd
  for (col in extra_covars) {
    out[[col]] <- df[[col]]
  }
  out <- cbind(out, harmonized)

  out_path <- file.path(output_dir, paste0("SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_", variant, ".rds"))
  saveRDS(out, out_path)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

run_variant(scdata, "age_sex_meanfd", character(0), baseline_only = FALSE)
run_variant(scdata, "cognition", "nihtbx_fluidcomp_uncorrected", baseline_only = TRUE)
run_variant(scdata, "pfactor", "GENERAL", baseline_only = FALSE)
