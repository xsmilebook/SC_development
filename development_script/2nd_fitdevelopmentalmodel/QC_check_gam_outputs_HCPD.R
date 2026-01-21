## QC script for HCP-D developmental GAM outputs (ComBat-GAM input).
## Goal: validate that gamresults/gammodels are consistent and that mgcv::predict.gam
## works for each edge (including se.fit path when available).
##
## Usage:
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/QC_check_gam_outputs_HCPD.R --project_root=/path/to/SCDevelopment --cvthr=75
##   Rscript --vanilla .../QC_check_gam_outputs_HCPD.R --cvthr=75 --force=1

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

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
force <- as.integer(if (!is.null(args$force)) args$force else 0L) == 1L

ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2
expected_parcels <- paste0("SC.", seq_len(elementnum), "_h")

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr), "qc"
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  gamresults_raw = file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds")),
  gammodel_raw = file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, ".rds")),
  gamresults_scaled = file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")),
  gammodel_scaled = file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")),
  plotdatasum_scaled = file.path(interfileFolder, "plotdatasum_scale_TRUE_SA12.rds")
)

for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) stop("Missing required file: ", nm, " => ", paths[[nm]])
}

qc_out_summary <- file.path(resultFolder, paste0("qc_summary_CV", CVthr, ".csv"))
qc_out_edges <- file.path(resultFolder, paste0("qc_edges_CV", CVthr, ".csv"))
qc_out_rds <- file.path(resultFolder, paste0("qc_details_CV", CVthr, ".rds"))

if (!force && file.exists(qc_out_summary) && file.exists(qc_out_edges) && file.exists(qc_out_rds)) {
  message("[INFO] QC outputs exist, skipping: ", resultFolder, " (set --force=1 to re-run)")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(mgcv))

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "plotdata_generate.R"))

read_obj <- function(path) {
  tryCatch(readRDS(path), error = function(e) structure(list(message = conditionMessage(e)), class = "qc_read_error"))
}

gamresults_raw <- read_obj(paths$gamresults_raw)
gammodel_raw <- read_obj(paths$gammodel_raw)
gamresults_scaled <- read_obj(paths$gamresults_scaled)
gammodel_scaled <- read_obj(paths$gammodel_scaled)
plotdatasum_scaled <- read_obj(paths$plotdatasum_scaled)

is_df <- function(x) inherits(x, "data.frame")
is_list <- function(x) is.list(x) && !inherits(x, "qc_read_error")

check_gamresults <- function(df) {
  out <- list()
  out$is_df <- is_df(df)
  if (!out$is_df) return(out)
  out$nrow <- nrow(df)
  out$has_parcel <- "parcel" %in% names(df)
  out$n_unique_parcel <- if (out$has_parcel) length(unique(df$parcel)) else NA_integer_
  out$n_missing_parcel <- if (out$has_parcel) sum(is.na(df$parcel) | df$parcel == "") else NA_integer_
  out$n_not_expected <- if (out$has_parcel) sum(!df$parcel %in% expected_parcels) else NA_integer_
  out$dupe_parcels <- if (out$has_parcel) any(duplicated(df$parcel)) else NA
  out
}

check_gammodels <- function(models) {
  out <- list()
  out$is_list <- is_list(models)
  if (!out$is_list) return(out)
  out$length <- length(models)
  out$n_null <- sum(vapply(models, is.null, logical(1)))
  out$n_gam <- sum(vapply(models, function(m) {
    if (is.null(m)) return(FALSE)
    inherits(m, "gam") || (!is.null(m$gam) && inherits(m$gam, "gam"))
  }, logical(1)))
  out
}

res_raw <- check_gamresults(gamresults_raw)
res_scaled <- check_gamresults(gamresults_scaled)
mod_raw <- check_gammodels(gammodel_raw)
mod_scaled <- check_gammodels(gammodel_scaled)

plot_scaled_type <- if (inherits(plotdatasum_scaled, "try-error")) "try-error" else if (is_list(plotdatasum_scaled)) "list" else class(plotdatasum_scaled)[1]
plot_scaled_len <- if (is_list(plotdatasum_scaled)) length(plotdatasum_scaled) else NA_integer_

edge_qc <- function(gamresults, gammodels, label) {
  if (!is_df(gamresults) || !is_list(gammodels)) {
    return(data.frame())
  }
  n <- min(nrow(gamresults), length(gammodels))
  if (n == 0) return(data.frame())

  parcels <- as.character(gamresults$parcel[seq_len(n)])
  details <- lapply(seq_len(n), function(i) {
    m <- gammodels[[i]]
    parcel <- parcels[[i]]
    cls <- if (is.null(m)) NA_character_ else paste(class(m), collapse = ",")
    gm <- if (inherits(m, "gam")) m else if (!is.null(m$gam) && inherits(m$gam, "gam")) m$gam else NULL
    ok_gam <- !is.null(gm)
    formula_chr <- if (ok_gam) paste0(deparse(formula(gm)), collapse = "") else NA_character_
    ok_terms <- if (!ok_gam) FALSE else {
      all(c("age", "sex", "mean_fd") %in% names(gm$model))
    }
    pred_ok <- FALSE
    se_ok <- FALSE
    pred_msg <- NA_character_
    if (ok_gam) {
      pd <- tryCatch(plotdata_generate(gm, "age"), error = function(e) structure(list(msg = conditionMessage(e)), class = "qc_pred_error"))
      if (!inherits(pd, "qc_pred_error")) {
        pred_ok <- is_df(pd) && all(c("fit", "se.fit") %in% names(pd))
        se_ok <- pred_ok && any(!is.na(pd$se.fit))
      } else {
        pred_msg <- pd$msg
      }
    } else {
      pred_msg <- "not a gam object"
    }
    data.frame(
      stage = label,
      idx = i,
      parcel = parcel,
      model_class = cls,
      is_gam = ok_gam,
      has_terms_age_sex_meanfd = ok_terms,
      formula = formula_chr,
      predict_ok = pred_ok,
      se_fit_available = se_ok,
      predict_error = pred_msg,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, details)
}

edges_raw <- edge_qc(gamresults_raw, gammodel_raw, "raw")
edges_scaled <- edge_qc(gamresults_scaled, gammodel_scaled, "scaled")
edges_all <- rbind(edges_raw, edges_scaled)

summary_df <- data.frame(
  CVthr = CVthr,
  elementnum = elementnum,
  gamresults_raw_n = if (!is.null(res_raw$nrow)) res_raw$nrow else NA_integer_,
  gamresults_raw_dupe_parcel = if (!is.null(res_raw$dupe_parcels)) res_raw$dupe_parcels else NA,
  gammodels_raw_len = if (!is.null(mod_raw$length)) mod_raw$length else NA_integer_,
  gammodels_raw_null = if (!is.null(mod_raw$n_null)) mod_raw$n_null else NA_integer_,
  gamresults_scaled_n = if (!is.null(res_scaled$nrow)) res_scaled$nrow else NA_integer_,
  gamresults_scaled_dupe_parcel = if (!is.null(res_scaled$dupe_parcels)) res_scaled$dupe_parcels else NA,
  gammodels_scaled_len = if (!is.null(mod_scaled$length)) mod_scaled$length else NA_integer_,
  gammodels_scaled_null = if (!is.null(mod_scaled$n_null)) mod_scaled$n_null else NA_integer_,
  plotdatasum_scaled_type = plot_scaled_type,
  plotdatasum_scaled_len = plot_scaled_len,
  edges_predict_fail = sum(!edges_all$predict_ok, na.rm = TRUE),
  edges_se_fit_missing = sum(edges_all$predict_ok & !edges_all$se_fit_available, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(summary_df, qc_out_summary, row.names = FALSE)
write.csv(edges_all, qc_out_edges, row.names = FALSE)
saveRDS(
  list(
    paths = paths,
    gamresults_raw_check = res_raw,
    gamresults_scaled_check = res_scaled,
    gammodel_raw_check = mod_raw,
    gammodel_scaled_check = mod_scaled,
    edges = edges_all
  ),
  qc_out_rds
)

message("[QC] Wrote: ", qc_out_summary)
message("[QC] Wrote: ", qc_out_edges)
message("[QC] Wrote: ", qc_out_rds)
