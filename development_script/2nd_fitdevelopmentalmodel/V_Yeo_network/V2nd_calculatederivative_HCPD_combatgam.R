## HCP-D (ComBat-GAM) | Yeo7/Yeo17 | Derivatives from scaled GAM models
##
## Runnable version adapted from:
##   V2nd_calculatederivative_HCPD.R
##
## Inputs:
## - outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/yeo/Yeo{7,17}/combat_gam/CV75/gammodel_Yeo*_scale_TRUE.rds
##
## Outputs:
## - outputs/results/2nd_fitdevelopmentalmodel/hcpd/yeo/Yeo{7,17}/combat_gam/CV75/*
##
## Usage:
##   Rscript --vanilla .../V2nd_calculatederivative_HCPD_combatgam.R --yeo=7
##   Rscript --vanilla .../V2nd_calculatederivative_HCPD_combatgam.R --yeo=17 --skip_posterior=1

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

infer_resolution_from_edges <- function(n_edges) {
  n <- (sqrt(8 * n_edges + 1) - 1) / 2
  n_int <- as.integer(round(n))
  if (n_int < 1 || n_int * (n_int + 1) / 2 != n_edges) {
    stop("Cannot infer resolution from n_edges=", n_edges, " (expected triangular number).")
  }
  n_int
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}
force <- as.integer(if (!is.null(args$force)) args$force else 0L) == 1L

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
yeo <- as.integer(if (!is.null(args$yeo)) args$yeo else 17L)
if (!yeo %in% c(7L, 17L)) stop("Unsupported yeo resolution: ", yeo, " (expected 7 or 17)")

skip_posterior <- as.integer(if (!is.null(args$skip_posterior)) args$skip_posterior else 0L) == 1L

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr)
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "gamderivatives.R"))

gamresult_files <- list.files(interfileFolder, pattern = paste0("^gamresults_Yeo\\d+_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE\\.rds$"), full.names = TRUE)
if (length(gamresult_files) != 1) stop("Expected exactly one scaled gamresults file in ", interfileFolder)
gamresultsum <- readRDS(gamresult_files[[1]])
if (!is.data.frame(gamresultsum) || nrow(gamresultsum) == 0) {
  stop("Empty gamresults inputs under: ", interfileFolder)
}

elementnum <- nrow(gamresultsum)
Yeoresolution.delLM <- infer_resolution_from_edges(elementnum)

gammodel_files <- list.files(interfileFolder, pattern = paste0("^gammodel_Yeo", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE\\.rds$"), full.names = TRUE)
if (length(gammodel_files) != 1) stop("Expected exactly one scaled gammodel file in ", interfileFolder)
gammodelsum <- readRDS(gammodel_files[[1]])

n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else elementnum)
n_edges <- min(elementnum, max(1L, n_edges))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

message("[INFO] yeo=", yeo, " Yeoresolution.delLM=", Yeoresolution.delLM, " elementnum=", elementnum)
message("[INFO] interfileFolder=", interfileFolder)
message("[INFO] resultFolder=", resultFolder)
message("[INFO] n_edges=", n_edges, " n_cores=", n_cores, " skip_posterior=", skip_posterior)

SCrank <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
for (x in 1:Yeoresolution.delLM) {
  for (y in 1:Yeoresolution.delLM) {
    SCrank[x, y] <- x^2 + y^2
  }
}
SCrank <- SCrank[lower.tri(SCrank, diag = TRUE)]
SCrank.df <- data.frame(parcel = paste0("SC.", 1:elementnum, "_h"), SCrank = SCrank)

out_derivative <- file.path(resultFolder, paste0("derivative.df_Yeo", elementnum, "_CV", CVthr, ".rds"))
out_posterior <- file.path(resultFolder, paste0("derivative.posterior.df.Yeo", yeo, "_CV", CVthr, ".rds"))

if (!force && file.exists(out_derivative)) {
  message("[INFO] derivative exists, skipping: ", out_derivative)
} else {
  derivative.sum <- mclapply(seq_len(min(n_edges, elementnum)), function(i) {
    SClabel.tmp <- gamresultsum$parcel[[i]]
    modobj <- gammodelsum[[i]]
    derivdata <- gam.derivatives(modobj, "age", draws = 1, increments = 1000, return_posterior_derivatives = FALSE)
    derivdata$label_ID <- gamresultsum$parcel[[i]]
    meanSC <- mean(modobj$model[, SClabel.tmp], na.rm = TRUE)
    derivdata$meanSC <- meanSC
    derivdata
  }, mc.cores = n_cores)

  derivative.df <- do.call(rbind, lapply(derivative.sum, function(x) data.frame(x)))
  saveRDS(derivative.df, out_derivative)
}

if (!skip_posterior) {
  if (!force && file.exists(out_posterior)) {
    message("[INFO] posterior derivative exists, skipping: ", out_posterior)
  } else {
    derivative.posterior.sum <- mclapply(seq_len(min(n_edges, elementnum)), function(i) {
      modobj <- gammodelsum[[i]]
      SClabel <- gamresultsum$parcel[[i]]
      derivdata <- gam.derivatives(modobj, "age", draws = 1000, increments = 1000, return_posterior_derivatives = TRUE)
      derivdata$SCrank <- SCrank.df$SCrank[SCrank.df$parcel == SClabel]
      meanSC <- mean(modobj$model[, SClabel], na.rm = TRUE)
      derivdata$meanSC <- meanSC
      maxSC <- max(modobj$model[, SClabel], na.rm = TRUE)
      derivdata$maxSC <- maxSC
      derivdata
    }, mc.cores = n_cores)
    saveRDS(derivative.posterior.sum, out_posterior)
  }
}
