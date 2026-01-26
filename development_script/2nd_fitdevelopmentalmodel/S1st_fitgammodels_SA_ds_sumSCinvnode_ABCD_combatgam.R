## ABCD (ComBat-GAM) | Fit GAMM models per edge (raw + scaled)
##
## Purpose:
## - Fit edge-wise GAMM to estimate SC trajectories.
## - Generate plotdata (raw) for trajectory visualizations.
## - Scale each edge by its fitted value at the youngest age (ratio-to-baseline),
##   then refit models on the scaled data ("_scale_TRUE").
##
## Inputs:
## - outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds
##
## Outputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gamresults78_sumSCinvnode_over8_siteall_CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gammodel78_sumSCinvnode_over8_siteall_CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/SCdata.diw_SA12CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gammodel78_sumSCinvnode_over8_CV75_scale_TRUE.rds
##
## Usage:
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R
##   Rscript --vanilla .../S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R --cvthr=75 --n_edges=10 --force=1

rm(list = ls())

library(parallel)
library(mgcv)
library(gamm4)
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
CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2
n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else elementnum)
n_edges <- min(elementnum, max(1L, n_edges))

dataset <- "abcd"

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  dataset, "combat_gam", paste0("CV", CVthr)
)
dir.create(interfileFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "abcd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)

out_gamresults_raw <- file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_siteall_CV", CVthr, ".rds"))
out_gammodel_raw <- file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_siteall_CV", CVthr, ".rds"))
out_plotdatasum_df <- file.path(interfileFolder, paste0("plotdatasum.df_SA", ds.resolution, "_sumSCinvnode_siteall_CV", CVthr, ".rds"))
out_scdata_diw <- file.path(interfileFolder, paste0("SCdata.diw_SA", ds.resolution, "CV", CVthr, ".rds"))
out_gamresults_scaled <- file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
out_gammodel_scaled <- file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
out_baseline_csv <- file.path(interfileFolder, paste0("baseline_fits_sumSCinvnode_CV", CVthr, ".csv"))

should_run <- function(path) force || !file.exists(path)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

message("[INFO] project_root=", project_root)
message("[INFO] dataset=", dataset)
message("[INFO] CVthr=", CVthr)
message("[INFO] input_rds=", input_rds)
message("[INFO] interfileFolder=", interfileFolder)
message("[INFO] n_edges=", n_edges, " n_cores=", n_cores)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "gammsmooth.R"))
source(file.path(functionFolder, "plotdata_generate.R"))

SCdata <- readRDS(input_rds)
if (!is.data.frame(SCdata)) stop("Expected a data.frame in input_rds: ", input_rds)

required <- c("age", "sex", "mean_fd")
missing_required <- setdiff(required, names(SCdata))
if (length(missing_required) > 0) stop("Missing required columns in ABCD SCdata: ", paste(missing_required, collapse = ", "))

if (!"subID" %in% names(SCdata)) {
  candidates <- c("src_subject_id", "subjectkey", "subject_id", "participant_id", "id")
  found <- candidates[candidates %in% names(SCdata)][1]
  if (is.na(found) || is.null(found)) {
    stop("Missing subject id column. Provide `subID` or one of: ", paste(candidates, collapse = ", "))
  }
  SCdata$subID <- as.character(SCdata[[found]])
  message("[INFO] subID mapped from column: ", found)
}
SCdata$subID <- as.factor(SCdata$subID)

age_max <- suppressWarnings(max(as.numeric(SCdata$age), na.rm = TRUE))
if (is.finite(age_max) && age_max > 30) {
  SCdata$age <- as.numeric(SCdata$age) / 12
  message("[INFO] Detected age likely in months (max=", sprintf("%.2f", age_max), "); converted to years by /12.")
} else {
  SCdata$age <- as.numeric(SCdata$age)
}

SCdata$sex <- as.factor(SCdata$sex)
SCdata$mean_fd <- as.numeric(SCdata$mean_fd)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
sc_cols <- sc_cols[str_detect(sc_cols, "_h$")]
if (length(sc_cols) != elementnum) {
  stop("Expected ", elementnum, " SC edge columns (SC.*_h), got ", length(sc_cols))
}

keep_cols <- c("age", "sex", "mean_fd", "subID", sc_cols)
SCdata <- SCdata[, keep_cols, drop = FALSE]
SCdata <- SCdata[stats::complete.cases(SCdata), , drop = FALSE]
if (nrow(SCdata) < 50) stop("Too few complete rows after filtering: ", nrow(SCdata))

covariates <- "sex+mean_fd"
dataname <- "SCdata"
smooth_var <- "age"

safe_fit_row <- function(SClabel) {
  tryCatch(
    {
      gammresult <- gamm.fit.smooth(
        SClabel, dataname, smooth_var, covariates,
        knots = 3, set_fx = TRUE, stats_only = FALSE, mod_only = FALSE
      )
      as.data.frame(gammresult)
    },
    error = function(e) {
      message("[WARN] gamm.fit.smooth failed for ", SClabel, ": ", conditionMessage(e))
      NULL
    }
  )
}

## Raw GAMM results/models
if (should_run(out_gamresults_raw) || should_run(out_gammodel_raw)) {
  result_rows <- mclapply(sc_cols[seq_len(n_edges)], safe_fit_row, mc.cores = n_cores)
  model_objs <- mclapply(sc_cols[seq_len(n_edges)], function(SClabel) {
    tryCatch(
      gamm.fit.smooth(
        SClabel, dataname, smooth_var, covariates,
        knots = 3, set_fx = TRUE, stats_only = TRUE, mod_only = TRUE
      ),
      error = function(e) {
        message("[WARN] gamm model fit failed for ", SClabel, ": ", conditionMessage(e))
        NULL
      }
    )
  }, mc.cores = n_cores)

  ok_idx <- which(!vapply(result_rows, is.null, logical(1)) & !vapply(model_objs, is.null, logical(1)))
  sc_cols_ok <- sc_cols[ok_idx]
  if (length(ok_idx) == 0) stop("All GAMM fits failed; see warnings above.")
  if (length(ok_idx) < n_edges) {
    failed <- setdiff(sc_cols[seq_len(n_edges)], sc_cols_ok)
    writeLines(failed, file.path(interfileFolder, paste0("failed_edges_sumSCinvnode_CV", CVthr, ".txt")))
    saveRDS(ok_idx, file.path(interfileFolder, paste0("ok_edge_index_sumSCinvnode_CV", CVthr, ".rds")))
    message("[INFO] Successful edges: ", length(ok_idx), " / ", n_edges)
  }

  gamresultsum.df <- dplyr::bind_rows(result_rows[ok_idx])
  num_cols <- setdiff(names(gamresultsum.df), "parcel")
  gamresultsum.df[num_cols] <- lapply(gamresultsum.df[num_cols], as.numeric)
  gamresultsum.df$pfdr <- p.adjust(gamresultsum.df$bootstrap_pvalue, method = "fdr")
  gamresultsum.df$sig <- (gamresultsum.df$pfdr < 0.05)
  gammodelsum <- model_objs[ok_idx]

  saveRDS(gamresultsum.df, out_gamresults_raw)
  saveRDS(gammodelsum, out_gammodel_raw)
} else {
  gamresultsum.df <- readRDS(out_gamresults_raw)
  gammodelsum <- readRDS(out_gammodel_raw)
  sc_cols_ok <- as.character(gamresultsum.df$parcel)
}

## Plot data (raw fitted values)
if (should_run(out_plotdatasum_df)) {
  plotdatasum <- mclapply(seq_len(length(gammodelsum)), function(i) {
    modobj <- gammodelsum[[i]]
    plotdata <- plotdata_generate(modobj, "age")
    plotdata$SC_label <- sc_cols_ok[[i]]
    plotdata
  }, mc.cores = n_cores)
  plotdatasum.df <- dplyr::bind_rows(plotdatasum)
  saveRDS(plotdatasum.df, out_plotdatasum_df)
} else {
  plotdatasum.df <- readRDS(out_plotdatasum_df)
}

baseline_df <- data.frame(parcel = sc_cols_ok, baseline_fit = NA_real_)
if (nrow(plotdatasum.df) > 0) {
  baseline_df$baseline_fit <- vapply(sc_cols_ok, function(SClabel) {
    tmp <- plotdatasum.df[plotdatasum.df$SC_label == SClabel, , drop = FALSE]
    if (nrow(tmp) == 0) return(NA_real_)
    as.numeric(tmp$fit[[1]])
  }, numeric(1))
}
write.csv(baseline_df, out_baseline_csv, row.names = FALSE)

## Scale by baseline fit (ratio-to-baseline)
if (should_run(out_scdata_diw)) {
  SCdata.diw <- SCdata
  for (SClabel in sc_cols_ok) {
    plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label == SClabel, ]
    SCdata.diw[, SClabel] <- SCdata[, SClabel] / plotdata.tmp$fit[[1]]
  }
  saveRDS(SCdata.diw, out_scdata_diw)
} else {
  SCdata.diw <- readRDS(out_scdata_diw)
}

## Scaled GAMM models + results
dataname <- "SCdata.diw"
if (should_run(out_gammodel_scaled) || should_run(out_gamresults_scaled)) {
  scaled_models <- mclapply(sc_cols_ok, function(SClabel) {
    tryCatch(
      gamm.fit.smooth(
        SClabel, dataname, smooth_var, covariates,
        knots = 3, set_fx = TRUE, stats_only = TRUE, mod_only = TRUE
      ),
      error = function(e) {
        message("[WARN] scaled GAMM model fit failed for ", SClabel, ": ", conditionMessage(e))
        NULL
      }
    )
  }, mc.cores = n_cores)

  scaled_ok <- which(!vapply(scaled_models, is.null, logical(1)))
  if (length(scaled_ok) == 0) stop("All scaled GAMM fits failed; see warnings above.")
  sc_cols_scaled_ok <- sc_cols_ok[scaled_ok]
  scaled_models <- scaled_models[scaled_ok]

  scaled_rows <- mclapply(sc_cols_scaled_ok, function(SClabel) {
    safe_fit_row(SClabel)
  }, mc.cores = n_cores)

  scaled_ok2 <- which(!vapply(scaled_rows, is.null, logical(1)))
  if (length(scaled_ok2) == 0) stop("All scaled GAMM result extraction failed; see warnings above.")
  if (length(scaled_ok2) < length(sc_cols_scaled_ok)) {
    failed_scaled <- sc_cols_scaled_ok[setdiff(seq_along(sc_cols_scaled_ok), scaled_ok2)]
    writeLines(failed_scaled, file.path(interfileFolder, paste0("failed_edges_sumSCinvnode_CV", CVthr, "_scale_TRUE.txt")))
    message("[INFO] Successful scaled edges: ", length(scaled_ok2), " / ", length(sc_cols_scaled_ok))
  }

  sc_cols_scaled_ok <- sc_cols_scaled_ok[scaled_ok2]
  scaled_models <- scaled_models[scaled_ok2]
  scaled_rows <- scaled_rows[scaled_ok2]

  gamresultsum_scaled <- dplyr::bind_rows(scaled_rows)
  num_cols <- setdiff(names(gamresultsum_scaled), "parcel")
  gamresultsum_scaled[num_cols] <- lapply(gamresultsum_scaled[num_cols], as.numeric)
  saveRDS(gamresultsum_scaled, out_gamresults_scaled)
  saveRDS(scaled_models, out_gammodel_scaled)
}
