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

resultsum <- mclapply(seq_len(n_edges), function(x) {
  SClabel <- sc_cols[[x]]
  region <- SClabel
  gamresult <- gam.fit.smooth(
    region, dataname, smooth_var, covariates,
    knots = 3, set_fx = TRUE, stats_only = FALSE, mod_only = FALSE
  )
  as.data.frame(gamresult)
}, mc.cores = n_cores)

gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df[, c(2:18)] <- lapply(gamresultsum.df[, c(2:18)], as.numeric)
gamresultsum.df$pfdr <- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
gamresultsum.df$sig <- (gamresultsum.df$pfdr < 0.05)
saveRDS(gamresultsum.df, file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds")))

## calculate gam models
resultsum <- mclapply(seq_len(n_edges), function(x) {
  SClabel <- sc_cols[[x]]
  region <- SClabel
  gam.fit.smooth(
    region, dataname, smooth_var, covariates,
    knots = 3, set_fx = TRUE, stats_only = TRUE, mod_only = TRUE
  )
}, mc.cores = n_cores)
saveRDS(resultsum, file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds")))

## plot data raw (fitted values)
gammodelsum <- resultsum
plotdatasum <- mclapply(seq_len(n_edges), function(x) {
  modobj <- gammodelsum[[x]]
  plotdata <- plotdata_generate(modobj, "age")
  plotdata$SC_label <- sc_cols[[x]]
  plotdata
}, mc.cores = n_cores)
plotdatasum.df <- do.call(rbind, lapply(plotdatasum, function(x) data.frame(x)))
saveRDS(plotdatasum.df, file.path(interfileFolder, paste0("plotdatasum.df_SA", ds.resolution, "_sumSCinvnode_CV", CVthr, ".rds")))

# Scale SC strength by the fitted value at the youngest age point (used for derivatives and visualization).
SCdata.diw <- SCdata.sum.merge
for (i in seq_len(n_edges)) {
  SClabel <- sc_cols[[i]]
  plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label == SClabel, ]
  SCdata.diw[, SClabel] <- SCdata.sum.merge[, SClabel] / plotdata.tmp$fit[[1]]
}
saveRDS(SCdata.diw, file.path(interfileFolder, paste0("SCdata.diw_SA", ds.resolution, "CV", CVthr, ".rds")))

## scaled gam models
covariates <- "sex+mean_fd"
dataname <- "SCdata.diw"
smooth_var <- "age"
resultsum <- mclapply(seq_len(n_edges), function(x) {
  SClabel <- sc_cols[[x]]
  region <- SClabel
  gam.fit.smooth(
    region, dataname, smooth_var, covariates,
    knots = 3, set_fx = TRUE, stats_only = TRUE, mod_only = TRUE
  )
}, mc.cores = n_cores)
saveRDS(resultsum, file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))

## scaled gam results
resultsum <- mclapply(seq_len(n_edges), function(x) {
  SClabel <- sc_cols[[x]]
  region <- SClabel
  gamresult <- gam.fit.smooth(
    region, dataname, smooth_var, covariates,
    knots = 3, set_fx = TRUE, stats_only = FALSE, mod_only = FALSE
  )
  as.data.frame(gamresult)
}, mc.cores = n_cores)

gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df[, c(2:18)] <- lapply(gamresultsum.df[, c(2:18)], as.numeric)
saveRDS(gamresultsum.df, file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))

