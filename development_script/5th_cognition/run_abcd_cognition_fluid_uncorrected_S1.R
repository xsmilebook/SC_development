#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(psych)
  library(reshape)
  library(RColorBrewer)
  library(tidyverse)
})

rm(list = ls())

CVthr <- 75
Cogvar <- "nihtbx_fluidcomp_uncorrected"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "cognition")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "cognition")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (cognition variant)")
}

euclid_csv <- Sys.getenv(
  "ABCD_EUCLID_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "average_EuclideanDistance_12.csv")
)
if (!file.exists(euclid_csv)) {
  stop("Missing ABCD_EUCLID_CSV: ", euclid_csv)
}

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "gamcog.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

SCdata <- readRDS(input_rds)
meandistance <- read.csv(euclid_csv)$Edistance
SCdata$age <- as.numeric(SCdata$age) / 12

if (!Cogvar %in% names(SCdata)) {
  stop("Missing cognition variable in input: ", Cogvar)
}

nonna_index <- which(!is.na(SCdata[, Cogvar]))
SCdata.cog <- SCdata[nonna_index, , drop = FALSE]
if ("eventname" %in% names(SCdata.cog)) {
  SCdata.cog <- SCdata.cog[SCdata.cog$eventname == "baseline_year_1_arm_1", , drop = FALSE]
}

cogagemodel <- gam(stats::as.formula(paste0(Cogvar, "~ s(age,k=3, fx=TRUE)+sex+mean_fd")), data = SCdata.cog)
t <- summary(cogagemodel)
message("age, sex, mean_fd can explain ", round(t$r.sq, 3), " variance of cognition.")

dataname <- "SCdata.cog"
smooth_var <- "age"
covariates <- "sex+mean_fd"
knots <- 3
corrmethod <- "pearson"

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "60"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 60
num_cores <- min(num_cores, 60L)

run_mclapply_with_fallback <- function(x, fun, n_workers) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1) n_workers <- 1L
  tries <- c(n_workers, 60L, 50L, 40L, 30L, 20L, 16L, 8L, 4L, 2L, 1L)
  tries <- unique(tries[tries <= n_workers])
  for (k in tries) {
    message("[INFO] Trying mclapply with mc.cores=", k)
    out <- try(parallel::mclapply(x, fun, mc.cores = k), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    msg <- as.character(out)
    message("[WARN] mclapply failed: ", msg)
    message("[WARN] Retrying with fewer workers.")
  }
  stop("mclapply failed for all worker settings; likely process/memory limits.")
}

make_error_row <- function(parcel, cognition_var, err) {
  data.frame(
    parcel = parcel,
    cognition_var = cognition_var,
    gam.smooth.t = NA_real_,
    gam.smooth.pvalue = NA_real_,
    anova.cov.pvalue = NA_real_,
    anova.cov.int.pvalue = NA_real_,
    partialRsq.int = NA_real_,
    partialRsq = NA_real_,
    correstimate = NA_real_,
    corrp = NA_real_,
    error = as.character(err),
    stringsAsFactors = FALSE
  )
}

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("SC_Cog_results_", Cogvar, "_CV", CVthr, "_cognition.rds"))

if (force || !file.exists(out_rds)) {
  sc_labels <- grep("SC\\.", names(SCdata), value = TRUE)
  if (length(sc_labels) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_labels))

  resultsum <- run_mclapply_with_fallback(seq_len(78), function(i) {
    SClabel <- sc_labels[[i]]
    tryCatch(
      {
        gamresult <- gam.fit.cognition(
          SClabel,
          dataname,
          Cogvar,
          smooth_var,
          covariates,
          knots,
          corrmethod,
          set_fx = TRUE,
          stats_only = TRUE
        )
        list(ok = TRUE, row = as.data.frame(gamresult), err = NA_character_)
      },
      error = function(e) {
        list(ok = FALSE, row = NULL, err = conditionMessage(e))
      }
    )
  }, num_cores)

  ok_mask <- vapply(resultsum, function(z) isTRUE(z$ok), logical(1))
  if (!any(ok_mask)) {
    errs <- unique(vapply(resultsum, function(z) as.character(z$err), character(1)))
    stop("All edges failed inside mclapply. First errors:\n", paste(head(errs, 10), collapse = "\n"))
  }

  rows <- vector("list", length = 78)
  for (i in seq_len(78)) {
    if (isTRUE(resultsum[[i]]$ok)) {
      rows[[i]] <- resultsum[[i]]$row
    } else {
      rows[[i]] <- make_error_row(sc_labels[[i]], Cogvar, resultsum[[i]]$err)
    }
  }
  SC_Cog_results.df <- do.call(rbind, rows)

  if ("error" %in% names(SC_Cog_results.df)) {
    failed_n <- sum(!is.na(SC_Cog_results.df$error) & nzchar(SC_Cog_results.df$error))
    message("[INFO] Edge-level failures: ", failed_n, " / 78")
  }

  SC_Cog_results.df[, c(3:10)] <- lapply(SC_Cog_results.df[, c(3:10)], as.numeric)

  if (!"corrp" %in% names(SC_Cog_results.df)) stop("Missing column: corrp")
  if (!"anova.cov.pvalue" %in% names(SC_Cog_results.df)) stop("Missing column: anova.cov.pvalue")
  if (!"gam.smooth.pvalue" %in% names(SC_Cog_results.df)) stop("Missing column: gam.smooth.pvalue")

  SC_Cog_results.df$corr.p.fdr <- p.adjust(as.numeric(SC_Cog_results.df$corrp), method = "fdr")
  SC_Cog_results.df$anova.cov.p.fdr <- p.adjust(as.numeric(SC_Cog_results.df$anova.cov.pvalue), method = "fdr")
  SC_Cog_results.df$gam.smooth.p.fdr <- p.adjust(as.numeric(SC_Cog_results.df$gam.smooth.pvalue), method = "fdr")

  saveRDS(SC_Cog_results.df, out_rds)
} else {
  SC_Cog_results.df <- readRDS(out_rds)
}

message(sum(SC_Cog_results.df$anova.cov.p.fdr < 0.05), " edges have significant associations with ", Cogvar, ".")

SC_Cog_results.df.whole <- SC_Cog_results.df
SCrankresult.whole <- SCrankcorr(SC_Cog_results.df.whole, "gam.smooth.t", 12)
message(
  "Correlation coefficient between cognitive associations and connectional axis is ",
  round(SCrankresult.whole$r.spearman, 2), " with P=", round(SCrankresult.whole$p.spearman, 3)
)

message("Control Euclidean distance.")
SC_Cog_results.df.whole$meandistance <- meandistance
SC_Cog_results.df.whole$gam.smooth.t_control_distance[which(!is.na(SC_Cog_results.df.whole$gam.smooth.t))] <-
  residuals(lm(gam.smooth.t ~ meandistance, data = SC_Cog_results.df.whole))
SCrankresult.whole.controllength <- SCrankcorr(SC_Cog_results.df.whole, "gam.smooth.t_control_distance", 12, dsdata = FALSE)

saveRDS(
  list(
    SCrankresult.whole = SCrankresult.whole,
    SCrankresult.whole.controllength = SCrankresult.whole.controllength
  ),
  file.path(resultFolder, paste0("SCrankcorr_summary_", Cogvar, "_CV", CVthr, "_cognition.rds"))
)

