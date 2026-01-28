## This script is to generate derivative and posterior derivative values for HCP-D GAM models.
## Age will be sampled from 8 to 22, with a total of 1000 data points.
## Covariables will be set as median or mode.
##
## Inputs are expected from:
##   outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/
## and outputs are written to:
##   outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/
##
## Usage (examples):
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_HCPD.R
##   Rscript --vanilla .../S2nd_calculatederivative_HCPD.R --skip_posterior=1 --n_edges=10

library(tidyverse)
library(R.matlab)
library(psych)
library(gratia)
library(mgcv)
library(parallel)

rm(list = ls())

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

# set resolution (default SA12; override for SA7/SA17 via --ds_res=7/17)
ds.resolution <- as.integer(if (!is.null(args$ds_res)) args$ds_res else 12L)
elementnum <- ds.resolution * (ds.resolution + 1) / 2

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else elementnum)
n_edges <- min(elementnum, max(1L, n_edges))
skip_posterior <- as.integer(if (!is.null(args$skip_posterior)) args$skip_posterior else 0L) == 1L

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

#### load data (scaled models)
gamresultsum <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
gammodelsum <- readRDS(file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
if (nrow(gamresultsum) == 0 || length(gammodelsum) == 0) {
  stop("Empty inputs: gamresults or gammodels. Check upstream outputs under: ", interfileFolder)
}

SCrank <- matrix(NA, ds.resolution, ds.resolution)
for (x in 1:ds.resolution) {
  for (y in 1:ds.resolution) {
    SCrank[x, y] <- x^2 + y^2
  }
}
SCrank <- SCrank[lower.tri(SCrank, diag = TRUE)]
SCrank.df <- data.frame(parcel = paste0("SC.", 1:elementnum, "_h"), SCrank = SCrank)

#### source function
source(file.path(functionFolder, "gamderivatives.R"))

## derivative
out_derivative <- file.path(resultFolder, paste0("derivative.df", elementnum, "_CV", CVthr, ".rds"))
out_derivative_posterior <- file.path(resultFolder, paste0("derivative.posterior.df.SA", ds.resolution, "_CV", CVthr, ".rds"))

if (!force && file.exists(out_derivative)) {
  message("[INFO] derivative.df exists, skipping: ", out_derivative)
} else {
derivative.sum <- mclapply(seq_len(min(n_edges, nrow(gamresultsum))), function(x) {
  SClabel.tmp <- gamresultsum$parcel[[x]]
  modobj <- gammodelsum[[x]]
  draws <- 1
  increments <- 1000
  derivdata <- gam.derivatives(modobj, "age", draws, increments, return_posterior_derivatives = FALSE)
  derivdata$label_ID <- gamresultsum$parcel[[x]]
  meanSC <- mean(modobj$model[, SClabel.tmp], na.rm = TRUE)
  derivdata$meanSC <- meanSC
  derivdata
}, mc.cores = n_cores)

derivative.df <- do.call(rbind, lapply(derivative.sum, function(x) data.frame(x)))
saveRDS(derivative.df, out_derivative)
}

## posterior derivative (optional; heavy)
if (!skip_posterior) {
  if (!force && file.exists(out_derivative_posterior)) {
    message("[INFO] posterior derivative exists, skipping: ", out_derivative_posterior)
  } else {
  derivative.posterior.sum <- mclapply(seq_len(min(n_edges, nrow(gamresultsum))), function(x) {
    modobj <- gammodelsum[[x]]
    draws <- 1000
    increments <- 1000
    SClabel <- gamresultsum$parcel[[x]]
    derivdata <- gam.derivatives(modobj, "age", draws, increments, return_posterior_derivatives = TRUE)
    derivdata$SCrank <- SCrank.df$SCrank[SCrank.df$parcel == SClabel]
    meanSC <- mean(modobj$model[, SClabel], na.rm = TRUE)
    derivdata$meanSC <- meanSC
    maxSC <- max(modobj$model[, SClabel], na.rm = TRUE)
    derivdata$maxSC <- maxSC
    derivdata
  }, mc.cores = n_cores)
  saveRDS(derivative.posterior.sum, out_derivative_posterior)
  }
}
