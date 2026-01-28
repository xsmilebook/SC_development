## HCP-D (ComBat-GAM SA12) | Sensitivity analysis: add SES / ICV covariate
##
## Generate derivative and posterior derivative values for the scaled GAM models.
##
## Inputs:
##   outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/covariates/<cov_tag>/combat_gam/CV75/
## Outputs:
##   outputs/results/2nd_fitdevelopmentalmodel/hcpd/covariates/<cov_tag>/combat_gam/CV75/

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

backfill_covariate_from_demo <- function(scdata, covariate, demo_csv) {
  if (!"subID" %in% names(scdata)) stop("SCdata missing subID")
  if (!file.exists(demo_csv)) stop("demo_csv not found: ", demo_csv)
  demo <- read.csv(demo_csv)
  if (!"subID" %in% names(demo)) stop("demo_csv missing subID: ", demo_csv)
  if (!covariate %in% names(demo)) stop("demo_csv missing covariate '", covariate, "': ", demo_csv)

  demo2 <- demo %>%
    dplyr::select(subID, !!covariate) %>%
    dplyr::mutate(subID = as.character(subID)) %>%
    dplyr::group_by(subID) %>%
    dplyr::summarise(!!covariate := median(.data[[covariate]], na.rm = TRUE), .groups = "drop")

  sc2 <- scdata %>%
    dplyr::mutate(subID = as.character(subID)) %>%
    dplyr::left_join(demo2, by = "subID", suffix = c("", ".demo"))

  if (!covariate %in% names(sc2)) return(sc2)
  demo_col <- paste0(covariate, ".demo")
  if (demo_col %in% names(sc2)) {
    sc2[[covariate]] <- dplyr::coalesce(sc2[[covariate]], sc2[[demo_col]])
    sc2[[demo_col]] <- NULL
  }
  sc2
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
skip_posterior <- as.integer(if (!is.null(args$skip_posterior)) args$skip_posterior else 0L) == 1L

cov_tag <- if (!is.null(args$cov_tag)) args$cov_tag else "SES"
add_covariate <- if (!is.null(args$add_covariate)) {
  args$add_covariate
} else {
  if (cov_tag == "ICV") "ICV" else "income.adj"
}
demo_csv <- if (!is.null(args$demo_csv)) args$demo_csv else file.path(project_root, "demopath", "HCPD_demo_behav.csv")
input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "hcpd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "covariates", cov_tag, "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "covariates", cov_tag, "combat_gam", paste0("CV", CVthr)
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

## If the original input data is missing the requested covariate column, backfill
## from demopath/HCPD_demo_behav.csv and save a project-local copy for auditing.
if (file.exists(input_rds)) {
  sc_tmp <- readRDS(input_rds)
  if (!add_covariate %in% names(sc_tmp) || any(is.na(sc_tmp[[add_covariate]]))) {
    n_missing_before <- if (add_covariate %in% names(sc_tmp)) sum(is.na(sc_tmp[[add_covariate]])) else nrow(sc_tmp)
    message("[INFO] Backfilling covariate '", add_covariate, "' from: ", demo_csv, " (missing before=", n_missing_before, ")")
    sc_tmp <- backfill_covariate_from_demo(sc_tmp, add_covariate, demo_csv)
    n_missing_after <- if (add_covariate %in% names(sc_tmp)) sum(is.na(sc_tmp[[add_covariate]])) else nrow(sc_tmp)
    message("[INFO] Covariate '", add_covariate, "' missing after backfill=", n_missing_after)
    out_sc_backfilled <- file.path(interfileFolder, paste0("SCdata_input_backfilled_", add_covariate, ".rds"))
    saveRDS(sc_tmp, out_sc_backfilled)
    message("[INFO] Saved backfilled input RDS (for reference): ", out_sc_backfilled)
  }
}

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

source(file.path(functionFolder, "gamderivatives.R"))

out_derivative <- file.path(resultFolder, paste0("derivative.df", elementnum, "_CV", CVthr, ".rds"))
out_derivative_posterior <- file.path(resultFolder, paste0("derivative.posterior.df.SA", ds.resolution, "_CV", CVthr, ".rds"))

if (!force && file.exists(out_derivative)) {
  message("[INFO] derivative.df exists; skipping: ", out_derivative)
} else {
  derivative.sum <- mclapply(seq_len(min(n_edges, nrow(gamresultsum))), function(x) {
    SClabel.tmp <- gamresultsum$parcel[[x]]
    modobj <- gammodelsum[[x]]
    draws <- 1
    increments <- 1000
    derivdata <- gam.derivatives(modobj, "age", draws, increments, return_posterior_derivatives = FALSE)
    derivdata$label_ID <- SClabel.tmp
    derivdata$meanSC <- mean(modobj$model[, SClabel.tmp], na.rm = TRUE)
    derivdata
  }, mc.cores = n_cores)

  derivative.df <- do.call(rbind, lapply(derivative.sum, function(x) data.frame(x)))
  saveRDS(derivative.df, out_derivative)
}

if (!skip_posterior) {
  if (!force && file.exists(out_derivative_posterior)) {
    message("[INFO] posterior derivative exists; skipping: ", out_derivative_posterior)
  } else {
    derivative.posterior.sum <- mclapply(seq_len(min(n_edges, nrow(gamresultsum))), function(x) {
      modobj <- gammodelsum[[x]]
      draws <- 1000
      increments <- 1000
      SClabel <- gamresultsum$parcel[[x]]
      derivdata <- gam.derivatives(modobj, "age", draws, increments, return_posterior_derivatives = TRUE)
      derivdata$SCrank <- SCrank.df$SCrank[SCrank.df$parcel == SClabel]
      derivdata$meanSC <- mean(modobj$model[, SClabel], na.rm = TRUE)
      derivdata$maxSC <- max(modobj$model[, SClabel], na.rm = TRUE)
      derivdata
    }, mc.cores = n_cores)
    saveRDS(derivative.posterior.sum, out_derivative_posterior)
  }
}
