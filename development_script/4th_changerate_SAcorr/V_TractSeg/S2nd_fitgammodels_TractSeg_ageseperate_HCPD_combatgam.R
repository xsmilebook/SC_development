## HCP-D (ComBat-GAM) | TractSeg (major-bundle) | Age-separated GAM analysis at flip age
##
## Purpose:
## - Split scaled SC edges into younger/older groups at the flip age (from S1).
## - Fit edge-wise GAM models in each subgroup and compute partial-Rsq vs S-A rank
##   correlations; generate scatter + matrix figures.
## - Print key numeric results to stdout (captured by SLURM logs).
##
## Inputs (project-relative defaults):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/tractseg/combat_gam/CV75/SCdata.diw_SA12CV75_TractSeg.rds
## - outputs/results/4th_changerate_SAcorr/hcpd/tractseg/combat_gam/CV75/alignment_summary_flip_age.csv
##
## Outputs:
## - outputs/intermediate/4th_changerate_SAcorr/hcpd/tractseg/combat_gam/CV75/gamresults*_TractSeg_*.rds
## - outputs/figures/4th_changerate_SAcorr/hcpd/tractseg/combat_gam/CV75/correlation_sumSCinvnode_SCrank_{younger,older}/*.tiff + *.pdf
##
## Usage:
##   Rscript --vanilla development_script/4th_changerate_SAcorr/V_TractSeg/S2nd_fitgammodels_TractSeg_ageseperate_HCPD_combatgam.R
##   Rscript --vanilla .../S2nd_fitgammodels_TractSeg_ageseperate_HCPD_combatgam.R --n_edges=10 --force=1

rm(list = ls())

library(mgcv)
library(parallel)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

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

as_int <- function(x, default) {
  if (is.null(x) || is.na(suppressWarnings(as.integer(x)))) return(as.integer(default))
  as.integer(x)
}

as_num <- function(x, default) {
  if (is.null(x) || is.na(suppressWarnings(as.numeric(x)))) return(as.numeric(default))
  as.numeric(x)
}

melt_df <- function(x, id.vars) {
  if (requireNamespace("reshape2", quietly = TRUE)) {
    return(reshape2::melt(x, id.vars = id.vars))
  }
  if (requireNamespace("reshape", quietly = TRUE)) {
    return(reshape::melt(x, id.vars = id.vars))
  }
  stop("Missing package for melt(): install reshape2 or reshape in the runtime environment.")
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

dataset <- "hcpd"
CVthr <- as_num(args$cvthr, 75)
force <- as_int(args$force, 0L) == 1L

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

inter2_dir <- file.path(project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel", dataset, "tractseg", "combat_gam", paste0("CV", CVthr))
scdata_diw_rds <- if (!is.null(args$scdata_diw_rds)) {
  args$scdata_diw_rds
} else {
  file.path(inter2_dir, paste0("SCdata.diw_SA12CV", CVthr, "_TractSeg.rds"))
}
if (!file.exists(scdata_diw_rds)) stop("Missing scdata_diw_rds: ", scdata_diw_rds)

result_dir <- file.path(project_root, "outputs", "results", "4th_changerate_SAcorr", dataset, "tractseg", "combat_gam", paste0("CV", CVthr))
flip_csv <- if (!is.null(args$flip_age_csv)) {
  args$flip_age_csv
} else {
  file.path(result_dir, "alignment_summary_flip_age.csv")
}

sepage <- if (!is.null(args$sepage)) {
  as_num(args$sepage, NA_real_)
} else if (file.exists(flip_csv)) {
  flip_df <- read.csv(flip_csv)
  if (!"flip_age_median" %in% names(flip_df)) stop("flip_age_csv missing flip_age_median: ", flip_csv)
  as.numeric(flip_df$flip_age_median[[1]])
} else {
  NA_real_
}
if (is.na(sepage)) stop("Missing sepage: pass --sepage=... or ensure flip CSV exists: ", flip_csv)

inter4_dir <- file.path(project_root, "outputs", "intermediate", "4th_changerate_SAcorr", dataset, "tractseg", "combat_gam", paste0("CV", CVthr))
figure_root <- file.path(project_root, "outputs", "figures", "4th_changerate_SAcorr", dataset, "tractseg", "combat_gam", paste0("CV", CVthr))
dir.create(inter4_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_root, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "gamsmooth.R"))
source(file.path(functionFolder, "SCrankcorr.R"))

SCdata.diw <- readRDS(scdata_diw_rds)
if (!is.data.frame(SCdata.diw)) stop("Expected data.frame in scdata_diw_rds: ", scdata_diw_rds)

required <- c("age", "sex", "mean_fd")
missing_required <- setdiff(required, names(SCdata.diw))
if (length(missing_required) > 0) stop("Missing required columns in SCdata.diw: ", paste(missing_required, collapse = ", "))

SCdata.diw$age <- as.numeric(SCdata.diw$age)
SCdata.diw$sex <- as.factor(SCdata.diw$sex)
SCdata.diw$mean_fd <- as.numeric(SCdata.diw$mean_fd)

sc_cols_all <- grep("^SC\\.", names(SCdata.diw), value = TRUE)
if (length(sc_cols_all) == 0) stop("No SC.* columns found in scdata_diw_rds: ", scdata_diw_rds)
sc_cols <- sc_cols_all
sc_cols_h <- sc_cols_all[grepl("_h$", sc_cols_all)]
if (length(sc_cols_h) > 0) sc_cols <- sc_cols_h

elementnum <- length(sc_cols)
matsize <- infer_resolution_from_edges(elementnum)

n_edges <- as_int(args$n_edges, elementnum)
n_edges <- min(elementnum, max(1L, n_edges))

message("[INFO] dataset=", dataset, " tractseg matsize=", matsize, " elementnum=", elementnum, " n_edges=", n_edges)
message("[INFO] scdata_diw_rds=", scdata_diw_rds)
message("[INFO] sepage=", sepage, " n_cores=", n_cores)

SCdata.younger <- SCdata.diw[SCdata.diw$age <= sepage, , drop = FALSE]
SCdata.older <- SCdata.diw[SCdata.diw$age > sepage, , drop = FALSE]
if (nrow(SCdata.younger) < 50 || nrow(SCdata.older) < 50) {
  message("[WARN] Small subgroup sizes: younger=", nrow(SCdata.younger), " older=", nrow(SCdata.older))
}

covariates <- "sex+mean_fd"
smooth_var <- "age"

should_run <- function(path) force || !file.exists(path)

fit_group <- function(dataname, data_df, set_fx, out_results_rds, out_models_rds) {
  if (should_run(out_results_rds)) {
    res_list <- mclapply(seq_len(n_edges), function(i) {
      SClabel <- sc_cols[[i]]
      tryCatch(
        as.data.frame(gam.fit.smooth(SClabel, dataname, smooth_var, covariates, knots = 3, set_fx = set_fx, stats_only = FALSE, mod_only = FALSE)),
        error = function(e) {
          message("[WARN] gam.fit.smooth failed for ", dataname, " ", SClabel, ": ", conditionMessage(e))
          NULL
        }
      )
    }, mc.cores = n_cores)
    res_df <- dplyr::bind_rows(res_list)
    if (nrow(res_df) == 0) stop("Empty GAM results for group: ", dataname)
    num_cols <- setdiff(names(res_df), "parcel")
    res_df[num_cols] <- lapply(res_df[num_cols], as.numeric)
    saveRDS(res_df, out_results_rds)
  }
  res_df <- readRDS(out_results_rds)

  if (should_run(out_models_rds)) {
    mod_list <- mclapply(seq_len(n_edges), function(i) {
      SClabel <- sc_cols[[i]]
      tryCatch(
        gam.fit.smooth(SClabel, dataname, smooth_var, covariates, knots = 3, set_fx = set_fx, stats_only = TRUE, mod_only = TRUE),
        error = function(e) {
          message("[WARN] gam model failed for ", dataname, " ", SClabel, ": ", conditionMessage(e))
          NULL
        }
      )
    }, mc.cores = n_cores)
    saveRDS(mod_list, out_models_rds)
  }

  res_df
}

out_young_results <- file.path(inter4_dir, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_younger_TractSeg.rds"))
out_old_results <- file.path(inter4_dir, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_older_TractSeg.rds"))
out_young_models <- file.path(inter4_dir, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_younger_TractSeg.rds"))
out_old_models <- file.path(inter4_dir, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_older_TractSeg.rds"))

gamresult_y <- fit_group("SCdata.younger", SCdata.younger, set_fx = TRUE, out_young_results, out_young_models)
gamresult_o <- fit_group("SCdata.older", SCdata.older, set_fx = FALSE, out_old_results, out_old_models)

gamresult_y$partialRsq1 <- as.numeric(gamresult_y$partialRsq)
mu <- mean(gamresult_y$partialRsq1, na.rm = TRUE)
sdv <- sd(gamresult_y$partialRsq1, na.rm = TRUE)
gamresult_y$partialRsq1[gamresult_y$partialRsq1 > mu + 3 * sdv | gamresult_y$partialRsq1 < mu - 3 * sdv] <- NA_real_

computevar_y <- "partialRsq1"
ct_y <- SCrankcorr(gamresult_y, computevar_y, matsize, dsdata = FALSE)
message(sprintf("[RESULT] dataset=%s tractseg younger %s: r=%.5f p=%.3g (n=%d)", dataset, computevar_y, ct_y$r.spearman[[1]], ct_y$p.spearman[[1]], nrow(gamresult_y)))

FigCorrYoung <- file.path(figure_root, "correlation_sumSCinvnode_SCrank_younger")
dir.create(FigCorrYoung, showWarnings = FALSE, recursive = TRUE)

correlation.df_y <- SCrankcorr(gamresult_y, computevar_y, matsize, dsdata = TRUE)
maxthr_y <- max(abs(correlation.df_y[[computevar_y]]), na.rm = TRUE)
p_sc_y <- ggplot(correlation.df_y) +
  geom_point(aes(x = SCrank, y = .data[[computevar_y]], color = .data[[computevar_y]]), size = 3) +
  geom_smooth(aes(x = SCrank, y = .data[[computevar_y]]), method = "lm", color = "black") +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-maxthr_y, maxthr_y)) +
  labs(x = "S-A connectional axis rank", y = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15),
    aspect.ratio = 0.8,
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )
ggsave(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_n", matsize, ".tiff")), p_sc_y, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_n", matsize, ".pdf")), p_sc_y, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")

make_matrix_plot <- function(values, matsize, maxthr) {
  mtrixplot <- matrix(NA_real_, matsize, matsize)
  mtrixplot[lower.tri(mtrixplot, diag = TRUE)] <- values[seq_len(min(length(values), matsize * (matsize + 1) / 2))]
  Matrix.tmp <- mtrixplot
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  mat_df <- as.data.frame(Matrix.tmp)
  mat_df$nodeid <- seq_len(matsize)
  mat_m <- melt_df(mat_df, id.vars = c("nodeid"))
  mat_m$variable <- as.numeric(mat_m$variable)
  mat_m$nodeid <- 0 - mat_m$nodeid
  mat_m$value <- as.numeric(mat_m$value)
  linerange_frame <- data.frame(
    x = c(0.5, matsize + 0.5), ymin = rep(-matsize - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -matsize - 0.5), xmin = rep(0.5, times = 2), xmax = rep(matsize + 0.5, times = 2)
  )
  ggplot(mat_m) +
    geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
    scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-maxthr, maxthr)) +
    scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-maxthr, maxthr)) +
    geom_segment(aes(x = 0.5, y = -0.5, xend = matsize + 0.5, yend = -0.5 - matsize), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
    labs(x = "", y = "") +
    scale_y_continuous(breaks = NULL, labels = NULL) +
    scale_x_continuous(breaks = NULL, labels = NULL) +
    theme(
      axis.line = element_line(linewidth = 0),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid.major = element_line(linewidth = 0),
      panel.grid.minor = element_line(linewidth = 1)
    )
}

maxthr_mat_y <- max(abs(gamresult_y$partialRsq1), na.rm = TRUE)
p_mat_y <- make_matrix_plot(gamresult_y$partialRsq1, matsize, maxthr_mat_y)
ggsave(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_matrix_n", matsize, ".tiff")), p_mat_y, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_matrix_n", matsize, ".pdf")), p_mat_y, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")

computevar_o <- "partialRsq"
ct_o <- SCrankcorr(gamresult_o, computevar_o, matsize, dsdata = FALSE)
message(sprintf("[RESULT] dataset=%s tractseg older %s: r=%.5f p=%.3g (n=%d)", dataset, computevar_o, ct_o$r.spearman[[1]], ct_o$p.spearman[[1]], nrow(gamresult_o)))

FigCorrOld <- file.path(figure_root, "correlation_sumSCinvnode_SCrank_older")
dir.create(FigCorrOld, showWarnings = FALSE, recursive = TRUE)

correlation.df_o <- SCrankcorr(gamresult_o, computevar_o, matsize, dsdata = TRUE)
maxthr_o <- max(abs(correlation.df_o[[computevar_o]]), na.rm = TRUE)
p_sc_o <- ggplot(correlation.df_o) +
  geom_point(aes(x = SCrank, y = .data[[computevar_o]], color = .data[[computevar_o]]), size = 3) +
  geom_smooth(aes(x = SCrank, y = .data[[computevar_o]]), method = "lm", color = "black") +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-maxthr_o, maxthr_o)) +
  labs(x = "S-A connectional axis rank", y = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15),
    aspect.ratio = 0.8,
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )
ggsave(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_n", matsize, ".tiff")), p_sc_o, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_n", matsize, ".pdf")), p_sc_o, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")

maxthr_mat_o <- max(abs(as.numeric(gamresult_o$partialRsq)), na.rm = TRUE)
p_mat_o <- make_matrix_plot(as.numeric(gamresult_o$partialRsq), matsize, maxthr_mat_o)
ggsave(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_matrix_n", matsize, ".tiff")), p_mat_o, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_matrix_n", matsize, ".pdf")), p_mat_o, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")

