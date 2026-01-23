## ABCD (ComBat-GAM) | Compute derivatives + posterior derivatives for scaled GAMM models
##
## Inputs are expected from:
##   outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/
## and outputs are written to:
##   outputs/results/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/
##
## Usage:
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_ABCD_combatgam.R
##   Rscript --vanilla .../S2nd_calculatederivative_ABCD_combatgam.R --skip_posterior=1 --n_edges=10 --force=1

library(parallel)
library(tidyverse)

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

dataset <- "abcd"
CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2
n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else elementnum)
n_edges <- min(elementnum, max(1L, n_edges))
skip_posterior <- as.integer(if (!is.null(args$skip_posterior)) args$skip_posterior else 0L) == 1L

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  dataset, "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  dataset, "combat_gam", paste0("CV", CVthr)
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "gamderivatives.R"))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

models_rds <- file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
results_rds <- file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
if (!file.exists(models_rds)) stop("Missing scaled models: ", models_rds)
if (!file.exists(results_rds)) stop("Missing scaled gamresults: ", results_rds)

gammodelsum <- readRDS(models_rds)
gamresultsum <- readRDS(results_rds)
if (!is.list(gammodelsum) || length(gammodelsum) == 0) stop("Empty gammodel list: ", models_rds)
if (!is.data.frame(gamresultsum) || !"parcel" %in% names(gamresultsum)) stop("Invalid gamresults: ", results_rds)

parcels <- as.character(gamresultsum$parcel)
if (length(parcels) != length(gammodelsum)) {
  message("[WARN] gamresults parcels and model list length differ; using intersection by names when available.")
  if (!is.null(names(gammodelsum)) && all(nzchar(names(gammodelsum)))) {
    parcels <- names(gammodelsum)
  } else {
    parcels <- parcels[seq_len(min(length(parcels), length(gammodelsum)))]
  }
}

# S-A axis value mapping for SC.1_h ... SC.78_h (lower triangle incl diag)
SCrank_mat <- matrix(NA_real_, nrow = ds.resolution, ncol = ds.resolution)
for (x in 1:ds.resolution) for (y in 1:ds.resolution) SCrank_mat[x, y] <- x^2 + y^2
SCrank_vec <- SCrank_mat[lower.tri(SCrank_mat, diag = TRUE)]
parcel_all <- paste0("SC.", seq_len(elementnum), "_h")
SCrank_map <- setNames(SCrank_vec, parcel_all)

out_derivative <- file.path(resultFolder, paste0("derivative.df", elementnum, "_CV", CVthr, ".rds"))
out_posterior <- file.path(resultFolder, paste0("derivative.posterior.df.SA", ds.resolution, "_CV", CVthr, ".rds"))

if (!force && file.exists(out_derivative)) {
  message("[INFO] derivative.df exists, skipping: ", out_derivative)
} else {
  derivative.sum <- mclapply(seq_len(min(n_edges, length(gammodelsum))), function(i) {
    modobj <- gammodelsum[[i]]
    SClabel <- parcels[[i]]
    derivdata <- gam.derivatives(modobj, "age", draws = 1, increments = 1000, return_posterior_derivatives = FALSE)
    derivdata$label_ID <- SClabel
    SClabel2 <- gsub("`", "", SClabel)
    meanSC <- tryCatch(mean(modobj$gam$model[, SClabel2], na.rm = TRUE), error = function(e) NA_real_)
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
    posterior.sum <- mclapply(seq_len(min(n_edges, length(gammodelsum))), function(i) {
      modobj <- gammodelsum[[i]]
      SClabel <- parcels[[i]]
      derivdata <- gam.derivatives(modobj, "age", draws = 1000, increments = 1000, return_posterior_derivatives = TRUE)
      SClabel2 <- gsub("`", "", SClabel)
      derivdata$SCrank <- unname(SCrank_map[[SClabel2]])
      meanSC <- tryCatch(mean(modobj$gam$model[, SClabel2], na.rm = TRUE), error = function(e) NA_real_)
      maxSC <- tryCatch(max(modobj$gam$model[, SClabel2], na.rm = TRUE), error = function(e) NA_real_)
      derivdata$meanSC <- meanSC
      derivdata$maxSC <- maxSC
      derivdata
    }, mc.cores = n_cores)
    saveRDS(posterior.sum, out_posterior)
  }
}
