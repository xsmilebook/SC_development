## Chinese Cohort (ComBat-GAM) | Fit GAM models per edge (raw + scaled)
##
## Purpose:
## - Fit edge-wise GAM to estimate SC trajectories.
## - Generate plotdata (raw) for trajectory visualizations.
## - Scale each edge by its fitted value at the youngest age (ratio-to-baseline),
##   then refit models on the scaled data ("_scale_TRUE").
##
## Input (default, project-relative):
## - outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds
##
## Outputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/gammodel78_sumSCinvnode_over8_CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/plotdatasum.df_SA12_sumSCinvnode_CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/SCdata.diw_SA12CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/gammodel78_sumSCinvnode_over8_CV75_scale_TRUE.rds
##
## Usage:
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ChineseCohort_combatgam.R
##   Rscript --vanilla .../S1st_fitgammodels_SA_ds_sumSCinvnode_ChineseCohort_combatgam.R --cvthr=75 --n_edges=10 --force=1

rm(list = ls())

library(mgcv)
library(parallel)
library(tidyverse)

parse_args <- function(args) {
  res <- list()
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    res[[kv[[1]]]] <- kv[[2]]
  }
  res
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}
force <- as.integer(if (!is.null(args$force)) args$force else 0L) == 1L

ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else elementnum)
n_edges <- min(elementnum, max(1L, n_edges))

dataset <- "chinese"

functionFolder <- file.path(project_root, "gamfunction")
interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  dataset, "combat_gam", paste0("CV", CVthr)
)
dir.create(interfileFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "chinese",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

message("[INFO] project_root=", project_root)
message("[INFO] dataset=", dataset)
message("[INFO] CVthr=", CVthr)
message("[INFO] input_rds=", input_rds)
message("[INFO] interfileFolder=", interfileFolder)
message("[INFO] n_edges=", n_edges, " n_cores=", n_cores)

SCdata.sum.merge <- readRDS(input_rds)
if (!is.data.frame(SCdata.sum.merge)) stop("Expected a data.frame in input_rds: ", input_rds)

# Normalize covariate column names (Chinese output preserves original names from ComBat-GAM):
# Age/Sex -> age/sex for downstream functions.
if (!"age" %in% names(SCdata.sum.merge) && "Age" %in% names(SCdata.sum.merge)) {
  SCdata.sum.merge$age <- SCdata.sum.merge$Age
}
if (!"sex" %in% names(SCdata.sum.merge) && "Sex" %in% names(SCdata.sum.merge)) {
  SCdata.sum.merge$sex <- SCdata.sum.merge$Sex
}
if (!"mean_fd" %in% names(SCdata.sum.merge) && "meanFD" %in% names(SCdata.sum.merge)) {
  SCdata.sum.merge$mean_fd <- SCdata.sum.merge$meanFD
}

required_covars <- c("age", "sex", "mean_fd")
missing_required <- setdiff(required_covars, names(SCdata.sum.merge))
if (length(missing_required) > 0) stop("Missing required columns in Chinese SCdata: ", paste(missing_required, collapse = ", "))

SCdata.sum.merge$age <- as.numeric(SCdata.sum.merge$age)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
SCdata.sum.merge$mean_fd <- as.numeric(SCdata.sum.merge$mean_fd)

sc_cols_all <- grep("^SC\\.", names(SCdata.sum.merge), value = TRUE)
if (length(sc_cols_all) == 0) stop("No SC.* columns found in input_rds: ", input_rds)
sc_cols <- sc_cols_all
sc_cols_h <- sc_cols_all[grepl("_h$", sc_cols_all)]
if (length(sc_cols_h) > 0) sc_cols <- sc_cols_h
if (length(sc_cols) != elementnum) {
  stop("Expected ", elementnum, " SC edge columns, got ", length(sc_cols), ". Example cols: ", paste(head(sc_cols, 3), collapse = ", "))
}

keep_cols <- c(required_covars, sc_cols)
SCdata.sum.merge <- SCdata.sum.merge[, keep_cols, drop = FALSE]
SCdata.sum.merge <- SCdata.sum.merge[stats::complete.cases(SCdata.sum.merge), , drop = FALSE]
if (nrow(SCdata.sum.merge) < 50) stop("Too few complete rows after filtering: ", nrow(SCdata.sum.merge))

source(file.path(functionFolder, "gamsmooth.R"))
source(file.path(functionFolder, "plotdata_generate.R"))

out_gamresults_raw <- file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds"))
out_gammodel_raw <- file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds"))
out_plotdatasum_df <- file.path(interfileFolder, paste0("plotdatasum.df_SA", ds.resolution, "_sumSCinvnode_CV", CVthr, ".rds"))
out_scdata_diw <- file.path(interfileFolder, paste0("SCdata.diw_SA", ds.resolution, "CV", CVthr, ".rds"))
out_gammodel_scaled <- file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
out_gamresults_scaled <- file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))

should_run <- function(path) force || !file.exists(path)

safe_fit_row <- function(SClabel, dataname, smooth_var, covariates) {
  tryCatch(
    {
      gamresult <- gam.fit.smooth(
        SClabel, dataname, smooth_var, covariates,
        knots = 3, set_fx = TRUE, stats_only = FALSE, mod_only = FALSE
      )
      as.data.frame(gamresult)
    },
    error = function(e) {
      message("[WARN] gam.fit.smooth failed for ", SClabel, ": ", conditionMessage(e))
      NULL
    }
  )
}

covariates <- "sex+mean_fd"
dataname <- "SCdata.sum.merge"
smooth_var <- "age"

## Raw GAM results/models
if (should_run(out_gamresults_raw) || should_run(out_gammodel_raw)) {
  result_rows <- mclapply(sc_cols[seq_len(n_edges)], function(SClabel) {
    safe_fit_row(SClabel, dataname, smooth_var, covariates)
  }, mc.cores = n_cores)

  model_objs <- mclapply(sc_cols[seq_len(n_edges)], function(SClabel) {
    tryCatch(
      gam.fit.smooth(
        SClabel, dataname, smooth_var, covariates,
        knots = 3, set_fx = TRUE, stats_only = TRUE, mod_only = TRUE
      ),
      error = function(e) {
        message("[WARN] gam model fit failed for ", SClabel, ": ", conditionMessage(e))
        NULL
      }
    )
  }, mc.cores = n_cores)

  ok_idx <- which(!vapply(result_rows, is.null, logical(1)) & !vapply(model_objs, is.null, logical(1)))
  sc_cols_ok <- sc_cols[ok_idx]
  if (length(ok_idx) == 0) stop("All GAM fits failed; see warnings above.")
  if (length(ok_idx) < n_edges) {
    failed <- setdiff(sc_cols[seq_len(n_edges)], sc_cols_ok)
    writeLines(failed, file.path(interfileFolder, paste0("failed_edges_sumSCinvnode_CV", CVthr, ".txt")))
    saveRDS(ok_idx, file.path(interfileFolder, paste0("ok_edge_index_sumSCinvnode_CV", CVthr, ".rds")))
    message("[INFO] Successful edges: ", length(ok_idx), " / ", n_edges)
    stop("GAM fits failed for ", length(failed), " edges; see failed_edges_sumSCinvnode_CV", CVthr, ".txt")
  }

  gamresultsum.df <- dplyr::bind_rows(result_rows[ok_idx])
  num_cols <- setdiff(names(gamresultsum.df), "parcel")
  gamresultsum.df[num_cols] <- lapply(gamresultsum.df[num_cols], as.numeric)
  gamresultsum.df$pfdr <- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
  gamresultsum.df$sig <- (gamresultsum.df$pfdr < 0.05)
  gammodelsum <- model_objs[ok_idx]

  saveRDS(gamresultsum.df, out_gamresults_raw)
  saveRDS(gammodelsum, out_gammodel_raw)
} else {
  gamresultsum.df <- readRDS(out_gamresults_raw)
  gammodelsum <- readRDS(out_gammodel_raw)
  sc_cols_ok <- as.character(gamresultsum.df$parcel)
}

## plot data raw (fitted values)
if (should_run(out_plotdatasum_df)) {
  plotdatasum <- mclapply(seq_len(length(gammodelsum)), function(x) {
    modobj <- gammodelsum[[x]]
    plotdata <- plotdata_generate(modobj, "age")
    plotdata$SC_label <- sc_cols_ok[[x]]
    plotdata
  }, mc.cores = n_cores)
  plotdatasum.df <- dplyr::bind_rows(plotdatasum)
  saveRDS(plotdatasum.df, out_plotdatasum_df)
} else {
  plotdatasum.df <- readRDS(out_plotdatasum_df)
}

## Scale SC strength by the fitted value at the youngest age point.
if (should_run(out_scdata_diw)) {
  SCdata.diw <- SCdata.sum.merge
  for (SClabel in sc_cols_ok) {
    plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label == SClabel, ]
    SCdata.diw[, SClabel] <- SCdata.sum.merge[, SClabel] / plotdata.tmp$fit[[1]]
  }
  saveRDS(SCdata.diw, out_scdata_diw)
} else {
  SCdata.diw <- readRDS(out_scdata_diw)
}

## scaled gam models + results
covariates <- "sex+mean_fd"
dataname <- "SCdata.diw"
smooth_var <- "age"
if (should_run(out_gammodel_scaled) || should_run(out_gamresults_scaled)) {
  scaled_models <- mclapply(sc_cols_ok, function(SClabel) {
    tryCatch(
      gam.fit.smooth(
        SClabel, dataname, smooth_var, covariates,
        knots = 3, set_fx = TRUE, stats_only = TRUE, mod_only = TRUE
      ),
      error = function(e) {
        message("[WARN] scaled gam model fit failed for ", SClabel, ": ", conditionMessage(e))
        NULL
      }
    )
  }, mc.cores = n_cores)
  scaled_ok <- which(!vapply(scaled_models, is.null, logical(1)))
  if (length(scaled_ok) == 0) stop("All scaled GAM fits failed; see warnings above.")
  sc_cols_scaled_ok <- sc_cols_ok[scaled_ok]
  scaled_models <- scaled_models[scaled_ok]

  scaled_rows <- mclapply(sc_cols_scaled_ok, function(SClabel) {
    safe_fit_row(SClabel, dataname, smooth_var, covariates)
  }, mc.cores = n_cores)

  scaled_ok2 <- which(!vapply(scaled_rows, is.null, logical(1)))
  if (length(scaled_ok2) == 0) stop("All scaled GAM result extraction failed; see warnings above.")
  if (length(scaled_ok2) < length(sc_cols_scaled_ok)) {
    failed_scaled <- sc_cols_scaled_ok[setdiff(seq_along(sc_cols_scaled_ok), scaled_ok2)]
    writeLines(failed_scaled, file.path(interfileFolder, paste0("failed_edges_sumSCinvnode_CV", CVthr, "_scale_TRUE.txt")))
    message("[INFO] Successful scaled edges: ", length(scaled_ok2), " / ", length(sc_cols_scaled_ok))
  }

  sc_cols_scaled_ok <- sc_cols_scaled_ok[scaled_ok2]
  scaled_models <- scaled_models[scaled_ok2]
  scaled_rows <- scaled_rows[scaled_ok2]

  gamresultsum.df <- dplyr::bind_rows(scaled_rows)
  num_cols <- setdiff(names(gamresultsum.df), "parcel")
  gamresultsum.df[num_cols] <- lapply(gamresultsum.df[num_cols], as.numeric)

  saveRDS(gamresultsum.df, out_gamresults_scaled)
  saveRDS(scaled_models, out_gammodel_scaled)
}

