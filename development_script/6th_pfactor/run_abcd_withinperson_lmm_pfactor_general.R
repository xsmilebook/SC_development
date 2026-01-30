#!/usr/bin/env Rscript

rm(list = ls())

# Prefer the active conda environment libraries, and avoid accidental user-library pollution on clusters.
Sys.unsetenv(c("R_LIBS_USER", "R_LIBS"))
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (nzchar(conda_prefix)) {
  .libPaths(c(file.path(conda_prefix, "lib", "R", "library"), .libPaths()))
}
.libPaths(.libPaths()[!grepl("/GPFS/.*/R/packages", .libPaths())])

suppressPackageStartupMessages({
  library(lme4)
  library(pbkrtest)
  library(parallel)
  library(dplyr)
  library(ggplot2)
})

CVthr <- 75
int_var <- "GENERAL"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "withinperson_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "6th_pfactor", "abcd", "withinperson_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
}

source(file.path(functionFolder, "lmminteraction.R"))

SCdata <- readRDS(input_rds)
needed <- c("subID", "age", "sex", "mean_fd", int_var)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$sex <- as.factor(SCdata$sex)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))

n_edges <- as.integer(Sys.getenv("N_EDGES", unset = "78"))
if (is.na(n_edges) || n_edges < 1) n_edges <- 78
n_edges <- min(n_edges, 78L)

pb_method <- Sys.getenv("PB_METHOD", unset = "KR")
pb_nsim <- as.integer(Sys.getenv("PB_NSIM", unset = "1000"))
if (is.na(pb_nsim) || pb_nsim < 1) pb_nsim <- 1000

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "40"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 40
num_cores <- min(num_cores, 40L)

out_rds <- file.path(resultFolder, paste0("lmmresult_time_by_pfactor_", int_var, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)

message("[INFO] Fitting LMM interaction models (n_edges=", n_edges, ", mc.cores=", num_cores, ", PB_METHOD=", pb_method, ", PB_NSIM=", pb_nsim, ")")
res_list <- parallel::mclapply(seq_len(n_edges), function(i) {
  region <- sc_cols[[i]]
  lmm.time.predict.covariateinteraction(
    region = region,
    dataname = "SCdata",
    age_var = "age",
    int_var = int_var,
    covariates = "sex+mean_fd",
    int_var_mode = "as_is",
    pb_method = pb_method,
    pb_nsim = pb_nsim,
    stats_only = TRUE
  )
}, mc.cores = num_cores)

if (any(vapply(res_list, inherits, logical(1), what = "try-error"))) {
  stop("At least one edge failed (mclapply returned try-error). Set N_EDGES small and rerun to locate the failing edge.")
}

lmmresult <- do.call(rbind, res_list)
lmmresult$p_time_int_fdr <- p.adjust(lmmresult$p_time_int, method = "fdr")
saveRDS(lmmresult, out_rds)
write.csv(lmmresult, out_csv, row.names = FALSE)
message("[INFO] Saved: ", out_rds)

SCdata$totalstrength <- rowMeans(SCdata[, sc_cols[seq_len(78)], drop = FALSE], na.rm = TRUE)
out_total <- lmm.time.predict.covariateinteraction(
  region = "totalstrength",
  dataname = "SCdata",
  age_var = "age",
  int_var = int_var,
  covariates = "sex+mean_fd",
  int_var_mode = "as_is",
  pb_method = pb_method,
  pb_nsim = pb_nsim,
  stats_only = FALSE
)
pred <- out_total$predicted
pred$lower <- pred$.fitted - 1.96 * pred$.se
pred$upper <- pred$.fitted + 1.96 * pred$.se

Fig <- ggplot(pred, aes(x = time, y = .fitted, group = int_level, linetype = int_level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = int_level), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9, color = "black") +
  geom_point(size = 2, color = "black") +
  theme_classic() +
  scale_linetype_manual(values = c(low10 = "dashed", high90 = "solid")) +
  scale_fill_manual(values = c(low10 = "grey60", high90 = "grey20")) +
  labs(x = "Years from baseline", y = "Predicted totalstrength (fixed effects only)", linetype = NULL, fill = NULL)

ggsave(file.path(FigureFolder, paste0("pred_totalstrength_time_by_", int_var, "_low10_high90.pdf")), Fig, width = 12, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("pred_totalstrength_time_by_", int_var, "_low10_high90.tiff")), Fig, width = 12, height = 10, units = "cm", bg = "transparent", dpi = 600)

message("[INFO] Done.")
