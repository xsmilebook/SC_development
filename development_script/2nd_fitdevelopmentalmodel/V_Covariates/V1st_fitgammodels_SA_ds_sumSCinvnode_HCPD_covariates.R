## HCP-D (ComBat-GAM SA12) | Sensitivity analysis: add SES / ICV covariate
##
## Fits GAM models for each edge with covariates:
##   sex + mean_fd + <add_covariate>
##
## Default input:
##   outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds
## If add_covariate is missing from the input RDS, it will be backfilled from:
##   demopath/HCPD_demo_behav.csv (by subID)
##
## Outputs are written under:
##   outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/covariates/<cov_tag>/combat_gam/CV75/
##
## Usage:
##   Rscript --vanilla .../V1st_fitgammodels_SA_ds_sumSCinvnode_HCPD_covariates.R --add_covariate=income.adj --cov_tag=SES
##   Rscript --vanilla .../V1st_fitgammodels_SA_ds_sumSCinvnode_HCPD_covariates.R --add_covariate=ICV --cov_tag=ICV

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

normalize_subid <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^sub[-_]", "", x)
  x
}

backfill_covariate_from_demo <- function(scdata, covariate, demo_csv) {
  if (!"subID" %in% names(scdata)) stop("SCdata missing subID")
  if (!file.exists(demo_csv)) stop("demo_csv not found: ", demo_csv)

  demo <- read.csv(demo_csv)
  if (!"subID" %in% names(demo)) stop("demo_csv missing subID: ", demo_csv)
  if (!covariate %in% names(demo)) stop("demo_csv missing covariate '", covariate, "': ", demo_csv)

  sc2 <- scdata %>%
    dplyr::mutate(subID_key = normalize_subid(subID))
  demo2 <- demo %>%
    dplyr::mutate(subID_key = normalize_subid(subID)) %>%
    dplyr::select(subID_key, !!covariate) %>%
    dplyr::group_by(subID_key) %>%
    dplyr::summarise(!!covariate := median(.data[[covariate]], na.rm = TRUE), .groups = "drop")

  sc2 <- sc2 %>%
    dplyr::left_join(demo2, by = "subID_key", suffix = c("", ".demo"))

  demo_col <- paste0(covariate, ".demo")
  if (demo_col %in% names(sc2)) {
    sc2[[covariate]] <- dplyr::coalesce(sc2[[covariate]], sc2[[demo_col]])
    sc2[[demo_col]] <- NULL
  }
  sc2$subID_key <- NULL
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

add_covariate <- if (!is.null(args$add_covariate)) args$add_covariate else "income.adj"
cov_tag <- if (!is.null(args$cov_tag)) args$cov_tag else if (add_covariate == "ICV") "ICV" else "SES"

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "hcpd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)

demo_csv <- if (!is.null(args$demo_csv)) args$demo_csv else file.path(project_root, "demopath", "HCPD_demo_behav.csv")

functionFolder <- file.path(project_root, "gamfunction")
interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "covariates", cov_tag, "combat_gam", paste0("CV", CVthr)
)
dir.create(interfileFolder, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

source(file.path(functionFolder, "gamsmooth.R"))
source(file.path(functionFolder, "plotdata_generate.R"))

SCdata.sum.merge <- readRDS(input_rds)
if (!"subID" %in% names(SCdata.sum.merge)) stop("input_rds missing column subID: ", input_rds)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)

sc_cols <- grep("^SC\\.", names(SCdata.sum.merge), value = TRUE)
if (length(sc_cols) != elementnum) {
  stop("Expected ", elementnum, " SC edge columns (SC.*), got ", length(sc_cols))
}

if (!add_covariate %in% names(SCdata.sum.merge) || any(is.na(SCdata.sum.merge[[add_covariate]]))) {
  n_missing_before <- if (add_covariate %in% names(SCdata.sum.merge)) sum(is.na(SCdata.sum.merge[[add_covariate]])) else nrow(SCdata.sum.merge)
  message("[INFO] Backfilling covariate '", add_covariate, "' from: ", demo_csv, " (missing before=", n_missing_before, ")")
  SCdata.sum.merge <- backfill_covariate_from_demo(SCdata.sum.merge, add_covariate, demo_csv)
  if (!add_covariate %in% names(SCdata.sum.merge)) {
    stop("Failed to provide add_covariate '", add_covariate, "' in analysis data after backfill.")
  }
  n_missing_after <- sum(is.na(SCdata.sum.merge[[add_covariate]]))
  message("[INFO] Covariate '", add_covariate, "' missing after backfill=", n_missing_after)
}

required_cov <- c("age", "sex", "mean_fd", add_covariate)
missing_required <- setdiff(required_cov, names(SCdata.sum.merge))
if (length(missing_required) > 0) stop("Missing required columns: ", paste(missing_required, collapse = ", "))

complete_mask <- complete.cases(SCdata.sum.merge[, required_cov, drop = FALSE])
message("[INFO] Complete cases for covariates (age, sex, mean_fd, ", add_covariate, "): ", sum(complete_mask), " / ", nrow(SCdata.sum.merge))
if (sum(complete_mask) < 30) {
  bad_ids <- head(unique(SCdata.sum.merge$subID[!complete_mask]), 10)
  stop("Too few complete cases (", sum(complete_mask), ") after covariate backfill; check subID alignment and missingness. Example subIDs with missing covariates: ",
       paste(bad_ids, collapse = ", "))
}
SCdata.sum.merge <- SCdata.sum.merge[complete_mask, , drop = FALSE]

covariates <- paste0("sex+mean_fd+", add_covariate)
dataname <- "SCdata.sum.merge"
smooth_var <- "age"

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

message("[INFO] cov_tag=", cov_tag, " add_covariate=", add_covariate, " covariates=", covariates)

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

## scaled gam models/results (for derivatives + visualization)
dataname <- "SCdata.diw"
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

  scaled_stats <- mclapply(sc_cols_scaled_ok, function(SClabel) {
    safe_fit_row(SClabel, dataname, smooth_var, covariates)
  }, mc.cores = n_cores)
  stats_ok <- which(!vapply(scaled_stats, is.null, logical(1)))
  if (length(stats_ok) == 0) stop("All scaled GAM stats fits failed; see warnings above.")
  sc_cols_final <- sc_cols_scaled_ok[stats_ok]
  scaled_models <- scaled_models[scaled_ok][stats_ok]

  gamresultsum_scaled.df <- dplyr::bind_rows(scaled_stats[stats_ok])
  num_cols <- setdiff(names(gamresultsum_scaled.df), "parcel")
  gamresultsum_scaled.df[num_cols] <- lapply(gamresultsum_scaled.df[num_cols], as.numeric)
  gamresultsum_scaled.df$pfdr <- p.adjust(gamresultsum_scaled.df$anova.smooth.pvalue, method = "fdr")
  gamresultsum_scaled.df$sig <- (gamresultsum_scaled.df$pfdr < 0.05)

  saveRDS(scaled_models, out_gammodel_scaled)
  saveRDS(gamresultsum_scaled.df, out_gamresults_scaled)
  message("[INFO] Scaled models: ", length(scaled_models), " edges; significant edges (pfdr<0.05): ", sum(gamresultsum_scaled.df$sig, na.rm = TRUE))
} else {
  message("[INFO] Scaled outputs exist; skipping scaled fit under: ", interfileFolder)
}
