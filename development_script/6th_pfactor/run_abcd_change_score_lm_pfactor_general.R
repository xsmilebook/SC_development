#!/usr/bin/env Rscript

rm(list = ls())

# Prefer the active conda environment libraries, and avoid accidental user-library pollution on clusters.
Sys.unsetenv(c("R_LIBS_USER", "R_LIBS"))
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (nzchar(conda_prefix)) {
  .libPaths(c(file.path(conda_prefix, "lib", "R", "library"), .libPaths()))
}
.libPaths(.libPaths()[!grepl("/GPFS/.*/R/packages", .libPaths())])

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
})

CVthr <- 75
int_var <- "GENERAL"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "change_score_lm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
}

plotdatasum_rds <- Sys.getenv(
  "ABCD_PLOTDATASUM_RDS",
  unset = file.path(
    project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
    "abcd", "combat_gam", "CV75", "plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds"
  )
)
if (!file.exists(plotdatasum_rds)) stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)

source(file.path(functionFolder, "lmfunction.R"))
source(file.path(functionFolder, "SCrankcorr.R"))

scanid_to_eventname <- function(scanID) {
  sess <- sub("^.*_ses-", "", as.character(scanID))
  sess <- gsub("([a-z])([A-Z])", "\\1_\\2", sess)
  sess <- gsub("([A-Za-z])([0-9])", "\\1_\\2", sess)
  sess <- gsub("([0-9])([A-Za-z])", "\\1_\\2", sess)
  tolower(sess)
}

SCdata <- readRDS(input_rds)
if (!("eventname" %in% names(SCdata)) && ("scanID" %in% names(SCdata))) {
  SCdata$eventname <- scanid_to_eventname(SCdata$scanID)
}
if (!("eventname" %in% names(SCdata))) {
  stop("Missing eventname (required to construct baseline age): input has no eventname/scanID")
}

needed <- c("subID", "age", "sex", "mean_fd", int_var)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$sex <- as.factor(SCdata$sex)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
sc_cols <- sc_cols[seq_len(78)]

plotdata <- readRDS(plotdatasum_rds)
if (!all(c("SC_label", "fit") %in% names(plotdata))) stop("plotdata missing SC_label/fit: ", plotdatasum_rds)
plot_fit <- plotdata$fit
names(plot_fit) <- as.character(plotdata$SC_label)
missing_fit <- setdiff(sc_cols, names(plot_fit))
if (length(missing_fit) > 0) stop("plotdata missing fits for edges: ", paste(head(missing_fit, 10), collapse = ", "))

# Normalize SC strength to ratio (divide by initial fit)
for (edge in sc_cols) {
  f0 <- as.numeric(plot_fit[[edge]])
  if (is.na(f0) || !is.finite(f0) || f0 == 0) stop("Invalid plotdata fit for edge: ", edge)
  SCdata[[edge]] <- as.numeric(SCdata[[edge]]) / f0
}

dataname <- "SCdata"
covariates <- "age_baseline + sex + mean_fd_mean"

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "40"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 40
num_cores <- min(num_cores, 40L)

n_edges <- as.integer(Sys.getenv("N_EDGES", unset = "78"))
if (is.na(n_edges) || n_edges < 1) n_edges <- 78
n_edges <- min(n_edges, 78L)

out_rds <- file.path(resultFolder, paste0("lmresult_change_score_pfactor_", int_var, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)

message("[INFO] Fitting change-score LM (n_edges=", n_edges, ", mc.cores=", num_cores, ")")
res_list <- parallel::mclapply(seq_len(n_edges), function(i) {
  region <- sc_cols[[i]]
  tryCatch(
    {
      lmresult <- lm.fit.change(
        region = region,
        dataname = dataname,
        int_var = int_var,
        covariates = covariates,
        stats_only = TRUE
      )
      list(ok = TRUE, row = lmresult, region = region, err = NA_character_)
    },
    error = function(e) list(ok = FALSE, row = NULL, region = region, err = conditionMessage(e))
  )
}, mc.cores = num_cores)

ok_mask <- vapply(res_list, function(z) isTRUE(z$ok), logical(1))
if (!all(ok_mask)) {
  first_bad <- which(!ok_mask)[[1]]
  stop("Edge failed: ", res_list[[first_bad]]$region, "\n", res_list[[first_bad]]$err)
}

lmresult <- do.call(rbind, lapply(res_list, `[[`, "row"))
lmresult$p_int_fdr <- p.adjust(lmresult$p_int, method = "fdr")
saveRDS(lmresult, out_rds)
write.csv(lmresult, out_csv, row.names = FALSE)
message("[INFO] Saved: ", out_rds)

message("[INFO] Correlation to connectional axis (beta_int)")
SCrank.df <- SCrankcorr(lmresult, "beta_int", 12, dsdata = FALSE)
saveRDS(SCrank.df, file.path(resultFolder, paste0("SCrankcorr_change_score_pfactor_", int_var, "_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr r=", round(SCrank.df$r.spearman, 3), " p=", signif(SCrank.df$p.spearman, 3))

message("[INFO] Done.")

