## ABCD (ComBat-GAM) | Fit GAMM models per edge (scaled)
##
## Purpose (for downstream derivative + alignment analyses):
## - Fit edge-wise GAMM (random intercept for subject) to estimate SC trajectories.
## - Scale each edge by its fitted value at the youngest age (ratio-to-baseline),
##   and refit models on the scaled data ("_scale_TRUE").
##
## Inputs:
## - outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds
##
## Outputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/SCdata.diw_SA12CV75.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gammodel78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/baseline_fits_sumSCinvnode_CV75.csv
##
## Usage:
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R
##   Rscript --vanilla .../S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R --cvthr=75 --n_edges=10 --force=1

rm(list = ls())

library(parallel)
library(mgcv)
library(gamm4)

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

out_scdata_diw <- file.path(interfileFolder, paste0("SCdata.diw_SA", ds.resolution, "CV", CVthr, ".rds"))
out_models_scaled <- file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
out_results_scaled <- file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds"))
out_baseline_csv <- file.path(interfileFolder, paste0("baseline_fits_sumSCinvnode_CV", CVthr, ".csv"))

should_run <- function(path) force || !file.exists(path)
if (!should_run(out_scdata_diw) && !should_run(out_models_scaled) && !should_run(out_results_scaled) && !should_run(out_baseline_csv)) {
  message("[INFO] Outputs exist; skipping (set --force=1 to re-run): ", interfileFolder)
  quit(save = "no", status = 0)
}

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

message("[INFO] project_root=", project_root)
message("[INFO] dataset=", dataset)
message("[INFO] CVthr=", CVthr)
message("[INFO] input_rds=", input_rds)
message("[INFO] interfileFolder=", interfileFolder)
message("[INFO] n_edges=", n_edges, " n_cores=", n_cores)

SCdata <- readRDS(input_rds)
if (!is.data.frame(SCdata)) stop("Expected a data.frame in input_rds: ", input_rds)

required <- c("age", "sex", "mean_fd")
missing_required <- setdiff(required, names(SCdata))
if (length(missing_required) > 0) stop("Missing required columns in ABCD SCdata: ", paste(missing_required, collapse = ", "))

# Map subject id column to `subID` for gamm4 random intercept
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

# Normalize age units to years (ABCD is often stored in months)
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
if (length(sc_cols) != elementnum) {
  stop("Expected ", elementnum, " SC edge columns (SC.*), got ", length(sc_cols))
}

# Drop rows with missing covariates or any SC edge missing (to keep model fitting stable)
keep_cols <- c("age", "sex", "mean_fd", "subID", sc_cols)
SCdata <- SCdata[, keep_cols, drop = FALSE]
SCdata <- SCdata[stats::complete.cases(SCdata), , drop = FALSE]
if (nrow(SCdata) < 50) stop("Too few complete rows after filtering: ", nrow(SCdata))

fit_one_edge <- function(SClabel, dat) {
  tryCatch(
    {
      f <- as.formula(sprintf("`%s` ~ s(age, k = 3, fx = TRUE) + sex + mean_fd", SClabel))
      gamm4::gamm4(f, random = ~(1 | subID), REML = TRUE, data = dat)
    },
    error = function(e) {
      message("[WARN] gamm4 failed for ", SClabel, ": ", conditionMessage(e))
      NULL
    }
  )
}

baseline_fit_one <- function(modobj) {
  if (is.null(modobj) || is.null(modobj$gam) || !inherits(modobj$gam, "gam")) return(NA_real_)
  gam <- modobj$gam
  df <- gam$model
  smooth_var <- "age"
  theseVars <- attr(gam$terms, "term.labels")
  varClasses <- attr(gam$terms, "dataClasses")
  pred <- data.frame(init = 0)
  for (v in theseVars) {
    thisClass <- varClasses[[v]]
    if (v == smooth_var) {
      pred[[v]] <- min(as.numeric(df[[v]]), na.rm = TRUE)
    } else {
      if (thisClass %in% c("factor", "ordered") || is.factor(df[[v]]) || is.ordered(df[[v]])) {
        tab <- table(df[[v]])
        lvl <- names(tab)[which.max(tab)][1]
        pred[[v]] <- factor(lvl, levels = levels(df[[v]]), ordered = is.ordered(df[[v]]))
      } else {
        pred[[v]] <- median(as.numeric(df[[v]]), na.rm = TRUE)
      }
    }
  }
  as.numeric(mgcv::predict.gam(gam, newdata = pred, se.fit = FALSE))[1]
}

message("[INFO] Fitting unscaled GAMM models (for baseline scaling)...")
raw_models <- mclapply(sc_cols[seq_len(n_edges)], function(SClabel) fit_one_edge(SClabel, SCdata), mc.cores = n_cores)
ok_raw <- which(!vapply(raw_models, is.null, logical(1)))
if (length(ok_raw) == 0) stop("All unscaled GAMM fits failed; see warnings above.")
sc_cols_ok <- sc_cols[ok_raw]
raw_models <- raw_models[ok_raw]
message("[INFO] Unscaled successful edges: ", length(sc_cols_ok), " / ", n_edges)

baseline_fits <- vapply(raw_models, baseline_fit_one, numeric(1))
baseline_df <- data.frame(parcel = sc_cols_ok, baseline_fit = baseline_fits)
write.csv(baseline_df, out_baseline_csv, row.names = FALSE)

bad_baseline <- which(!is.finite(baseline_fits) | baseline_fits <= 0)
if (length(bad_baseline) > 0) {
  message("[WARN] Non-positive/invalid baseline fits for ", length(bad_baseline), " edges; these edges will be dropped in scaled pipeline.")
  keep2 <- setdiff(seq_along(sc_cols_ok), bad_baseline)
  sc_cols_ok <- sc_cols_ok[keep2]
  baseline_fits <- baseline_fits[keep2]
}
if (length(sc_cols_ok) == 0) stop("No edges left after removing invalid baseline fits.")

message("[INFO] Creating scaled SCdata (ratio-to-baseline) ...")
SCdata_diw <- SCdata[, c("age", "sex", "mean_fd", "subID", sc_cols_ok), drop = FALSE]
for (i in seq_along(sc_cols_ok)) {
  SCdata_diw[[sc_cols_ok[[i]]]] <- SCdata_diw[[sc_cols_ok[[i]]]] / baseline_fits[[i]]
}
saveRDS(SCdata_diw, out_scdata_diw)

message("[INFO] Fitting scaled GAMM models ...")
scaled_models <- mclapply(sc_cols_ok, function(SClabel) fit_one_edge(SClabel, SCdata_diw), mc.cores = n_cores)
ok_scaled <- which(!vapply(scaled_models, is.null, logical(1)))
if (length(ok_scaled) == 0) stop("All scaled GAMM fits failed; see warnings above.")

scaled_models <- scaled_models[ok_scaled]
sc_cols_scaled_ok <- sc_cols_ok[ok_scaled]
names(scaled_models) <- sc_cols_scaled_ok

# Minimal results object for downstream steps (derivative scripts need `$parcel`).
gamresults_scaled <- data.frame(parcel = sc_cols_scaled_ok, stringsAsFactors = FALSE)

saveRDS(scaled_models, out_models_scaled)
saveRDS(gamresults_scaled, out_results_scaled)
message("[INFO] Saved scaled models: ", out_models_scaled)

