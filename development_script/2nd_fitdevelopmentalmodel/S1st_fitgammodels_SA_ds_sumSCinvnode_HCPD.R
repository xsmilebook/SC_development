## This script is to fit gam models for each edge (HCP-D).
## Statistical indexes and gam model files will be generated.
##
## Changes vs historical version:
## - All inputs/outputs are under the current SCDevelopment project root.
## - Default input uses the ComBat-GAM HCP-D output:
##   outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds
##
## Usage (examples):
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R
##   Rscript --vanilla .../S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R --n_edges=10
##   Rscript --vanilla .../S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R --input_rds=/path/to/*.rds --cvthr=75

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

# set resolution
ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2

# config
CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else elementnum)
n_edges <- min(elementnum, max(1L, n_edges))

functionFolder <- file.path(project_root, "gamfunction")
interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
dir.create(interfileFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "hcpd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

# load data
SCdata.sum.merge <- readRDS(input_rds)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
sc_cols <- grep("^SC\\.", names(SCdata.sum.merge), value = TRUE)
if (length(sc_cols) != elementnum) {
  stop("Expected ", elementnum, " SC edge columns (SC.*), got ", length(sc_cols))
}

# source function
source(file.path(functionFolder, "gamsmooth.R"))
source(file.path(functionFolder, "plotdata_generate.R"))

## calculate gam results
covariates <- "sex+mean_fd"
dataname <- "SCdata.sum.merge"
smooth_var <- "age"

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

result_rows <- mclapply(sc_cols[seq_len(n_edges)], function(SClabel) {
  safe_fit_row(SClabel, dataname, smooth_var, covariates)
}, mc.cores = n_cores)

## Models are fitted separately (historical behavior) but we must keep alignment:
## only keep edges that succeeded in BOTH the stats pass and the model pass.
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
}

gamresultsum.df <- dplyr::bind_rows(result_rows[ok_idx])
num_cols <- setdiff(names(gamresultsum.df), "parcel")
gamresultsum.df[num_cols] <- lapply(gamresultsum.df[num_cols], as.numeric)
gamresultsum.df$pfdr <- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
gamresultsum.df$sig <- (gamresultsum.df$pfdr < 0.05)
saveRDS(gamresultsum.df, file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds")))

## calculate gam models (only successful edges)
gammodelsum <- model_objs[ok_idx]
saveRDS(gammodelsum, file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds")))

## plot data raw (fitted values)
plotdatasum <- mclapply(seq_len(length(gammodelsum)), function(x) {
  modobj <- gammodelsum[[x]]
  plotdata <- plotdata_generate(modobj, "age")
  plotdata$SC_label <- sc_cols_ok[[x]]
  plotdata
}, mc.cores = n_cores)
plotdatasum.df <- dplyr::bind_rows(plotdatasum)
saveRDS(plotdatasum.df, file.path(interfileFolder, paste0("plotdatasum.df_SA", ds.resolution, "_sumSCinvnode_CV", CVthr, ".rds")))

# Scale SC strength by the fitted value at the youngest age point (used for derivatives and visualization).
SCdata.diw <- SCdata.sum.merge
for (SClabel in sc_cols_ok) {
  plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label == SClabel, ]
  SCdata.diw[, SClabel] <- SCdata.sum.merge[, SClabel] / plotdata.tmp$fit[[1]]
}
saveRDS(SCdata.diw, file.path(interfileFolder, paste0("SCdata.diw_SA", ds.resolution, "CV", CVthr, ".rds")))

## scaled gam models
covariates <- "sex+mean_fd"
dataname <- "SCdata.diw"
smooth_var <- "age"
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

## scaled gam results
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
saveRDS(gamresultsum.df, file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
saveRDS(scaled_models, file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
