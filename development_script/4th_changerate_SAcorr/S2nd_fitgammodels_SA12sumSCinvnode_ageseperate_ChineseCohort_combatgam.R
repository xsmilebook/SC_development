## Rscript version of:
##   development_script/4th_changerate_SAcorr/S2nd_fitgammodels_SA12sumSCinvnode_ageseperate_HCPD.Rmd
##
## Purpose:
## - Split Chinese Cohort (scaled) SC edges into younger/older groups at the flip age.
## - Fit edge-wise GAM models in each subgroup and generate partial-Rsq vs S-A rank
##   scatter + 12x12 matrix figures (matching the original Rmd outputs).
## - Print numeric correlation results to stdout for SLURM logs.
##
## Inputs (project-relative defaults):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/SCdata.diw_SA12CV75.rds
## - outputs/results/4th_changerate_SAcorr/hcpd/combat_gam/CV75/alignment_summary_flip_age.csv
##
## Outputs:
## - outputs/intermediate/4th_changerate_SAcorr/chinese/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_{young,old}.rds
## - outputs/intermediate/4th_changerate_SAcorr/chinese/combat_gam/CV75/gammodel78_sumSCinvnode_over8_CV75_{younger,older}.rds
## - outputs/figures/4th_changerate_SAcorr/chinese/combat_gam/CV75/{correlation_sumSCinvnode_SCrank_younger,...}/*.tiff + *.pdf
##
## Usage:
##   Rscript --vanilla development_script/4th_changerate_SAcorr/S2nd_fitgammodels_SA12sumSCinvnode_ageseperate_ChineseCohort_combatgam.R
##   Rscript --vanilla .../S2nd_fitgammodels...R --sepage=15.5 --n_edges=10 --force=1

rm(list = ls())

library(mgcv)
library(parallel)
library(ggplot2)
library(RColorBrewer)

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

save_svg_if_available <- function(filename, plot, width_cm, height_cm, dpi = 600) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("[INFO] svglite not available; skip svg: ", filename)
    return(invisible(FALSE))
  }
  ggsave(filename, plot, dpi = dpi, width = width_cm, height = height_cm, units = "cm", bg = "transparent")
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}

dataset <- "chinese"
CVthr <- as_num(args$cvthr, 75)
ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2
force <- as_int(args$force, 0L) == 1L
save_svg <- as_int(args$save_svg, 0L) == 1L

n_edges <- as_int(args$n_edges, elementnum)
n_edges <- min(elementnum, max(1L, n_edges))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

inter2_dir <- file.path(project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel", dataset, "combat_gam", paste0("CV", CVthr))
scdata_diw_rds <- if (!is.null(args$scdata_diw_rds)) {
  args$scdata_diw_rds
} else {
  file.path(inter2_dir, paste0("SCdata.diw_SA", ds.resolution, "CV", CVthr, ".rds"))
}
if (!file.exists(scdata_diw_rds)) stop("Missing scdata_diw_rds: ", scdata_diw_rds)

flip_csv <- if (!is.null(args$flip_age_csv)) {
  args$flip_age_csv
} else {
  file.path(project_root, "outputs", "results", "4th_changerate_SAcorr", dataset, "combat_gam", paste0("CV", CVthr), "alignment_summary_flip_age.csv")
}

sepage <- if (!is.null(args$sepage)) {
  as_num(args$sepage, NA_real_)
} else if (file.exists(flip_csv)) {
  flip_df <- read.csv(flip_csv)
  if (!"flip_age_median" %in% names(flip_df)) stop("flip_age_csv missing flip_age_median: ", flip_csv)
  as_num(flip_df$flip_age_median[[1]], NA_real_)
} else {
  NA_real_
}
if (!is.finite(sepage)) stop("Missing/invalid sepage. Provide --sepage=... or a valid --flip_age_csv=... (default: ", flip_csv, ")")

inter_out <- file.path(project_root, "outputs", "intermediate", "4th_changerate_SAcorr", dataset, "combat_gam", paste0("CV", CVthr))
dir.create(inter_out, showWarnings = FALSE, recursive = TRUE)

FigureRoot <- file.path(project_root, "outputs", "figures", "4th_changerate_SAcorr", dataset, "combat_gam", paste0("CV", CVthr))
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)

FigCorrYoung <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank_younger")
FigCorrOld <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank_older")
FigMatYoung <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats_younger")
FigMatOld <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats_older")
dir.create(FigCorrYoung, showWarnings = FALSE, recursive = TRUE)
dir.create(FigCorrOld, showWarnings = FALSE, recursive = TRUE)
dir.create(FigMatYoung, showWarnings = FALSE, recursive = TRUE)
dir.create(FigMatOld, showWarnings = FALSE, recursive = TRUE)

out_young_results <- file.path(inter_out, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_young.rds"))
out_old_results <- file.path(inter_out, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_old.rds"))
out_young_models <- file.path(inter_out, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_younger.rds"))
out_old_models <- file.path(inter_out, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_older.rds"))

message("[INFO] dataset=", dataset, " CVthr=", CVthr, " sepage=", sprintf("%.5f", sepage))
message("[INFO] scdata_diw_rds=", scdata_diw_rds)
message("[INFO] flip_age_csv=", flip_csv, " (exists=", file.exists(flip_csv), ")")
message("[INFO] n_edges=", n_edges, " n_cores=", n_cores, " force=", force)
message("[INFO] save_svg=", save_svg)

SCdata.sum.merge <- readRDS(scdata_diw_rds)
if (!is.data.frame(SCdata.sum.merge)) stop("SCdata.diw is not a data.frame: ", scdata_diw_rds)

required_cols <- c("age", "sex", "mean_fd")
missing_cols <- setdiff(required_cols, names(SCdata.sum.merge))
if (length(missing_cols) > 0) stop("SCdata.diw missing required columns: ", paste(missing_cols, collapse = ", "))

sc_cols <- grep("^SC\\.", names(SCdata.sum.merge), value = TRUE)
if (length(sc_cols) != elementnum) stop("Expected ", elementnum, " SC edge columns; got ", length(sc_cols))

SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
SCdata.sum.merge$mean_fd <- as.numeric(SCdata.sum.merge$mean_fd)
SCdata.sum.merge$age <- as.numeric(SCdata.sum.merge$age)

SCdata.sum.merge_younger <- SCdata.sum.merge[SCdata.sum.merge$age < sepage, , drop = FALSE]
SCdata.sum.merge_older <- SCdata.sum.merge[SCdata.sum.merge$age >= sepage, , drop = FALSE]
message("[INFO] n_younger=", nrow(SCdata.sum.merge_younger), " n_older=", nrow(SCdata.sum.merge_older))

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "gamsmooth.R"))
source(file.path(functionFolder, "SCrankcorr.R"))

should_run <- function(path) force || !file.exists(path)

covariates <- "sex+mean_fd"
smooth_var <- "age"

fit_group <- function(dataname, data_obj, set_fx, out_results_rds, out_models_rds) {
  assign(dataname, data_obj, envir = .GlobalEnv)

  if (should_run(out_results_rds)) {
    res_list <- mclapply(seq_len(min(n_edges, length(sc_cols))), function(i) {
      SClabel <- sc_cols[[i]]
      x <- tryCatch(
        {
          as.data.frame(gam.fit.smooth(SClabel, dataname, smooth_var, covariates, knots = 3, set_fx = set_fx, stats_only = FALSE, mod_only = FALSE))
        },
        error = function(e) {
          message("[WARN] gam.fit.smooth failed for ", dataname, " ", SClabel, ": ", conditionMessage(e))
          NULL
        }
      )
      x
    }, mc.cores = n_cores)

    res_df <- dplyr::bind_rows(res_list)
    if (nrow(res_df) == 0) stop("Empty gam results for group: ", dataname)
    num_cols <- setdiff(names(res_df), "parcel")
    res_df[num_cols] <- lapply(res_df[num_cols], as.numeric)
    saveRDS(res_df, out_results_rds)
  }

  res_df <- readRDS(out_results_rds)

  if (should_run(out_models_rds)) {
    mod_list <- mclapply(seq_len(min(n_edges, length(sc_cols))), function(i) {
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

  list(results = res_df)
}

young <- fit_group("SCdata.sum.merge_younger", SCdata.sum.merge_younger, set_fx = TRUE, out_young_results, out_young_models)
old <- fit_group("SCdata.sum.merge_older", SCdata.sum.merge_older, set_fx = FALSE, out_old_results, out_old_models)

# Younger: outlier-filtered partialRsq1
gamresultsum.df_younger <- young$results
gamresultsum.df_younger$partialRsq1 <- gamresultsum.df_younger$partialRsq
mu <- mean(gamresultsum.df_younger$partialRsq1, na.rm = TRUE)
sdv <- sd(gamresultsum.df_younger$partialRsq1, na.rm = TRUE)
gamresultsum.df_younger$partialRsq1[gamresultsum.df_younger$partialRsq1 > mu + 3 * sdv | gamresultsum.df_younger$partialRsq1 < mu - 3 * sdv] <- NA

computevar_y <- "partialRsq1"
ct_y <- SCrankcorr(gamresultsum.df_younger, computevar_y, ds.resolution, dsdata = FALSE)
message(sprintf("[RESULT] younger %s: r=%.5f p=%.3g (n=%d)", computevar_y, ct_y$r.spearman[[1]], ct_y$p.spearman[[1]], nrow(gamresultsum.df_younger)))

correlation.df <- SCrankcorr(gamresultsum.df_younger, computevar_y, ds.resolution, dsdata = TRUE)
maxthr <- max(abs(gamresultsum.df_younger$partialRsq1), na.rm = TRUE)

p_sc_y <- ggplot(correlation.df) +
  geom_point(aes(x = SCrank, y = partialRsq1, color = partialRsq1), size = 3) +
  geom_smooth(aes(x = SCrank, y = partialRsq1), method = "lm", color = "black") +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-maxthr, maxthr)) +
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
ggsave(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_n", ds.resolution, ".tiff")), p_sc_y, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_n", ds.resolution, ".pdf")), p_sc_y, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
if (save_svg) save_svg_if_available(file.path(FigCorrYoung, paste0(computevar_y, "_SCrankcorr_n", ds.resolution, ".svg")), p_sc_y, width_cm = 10, height_cm = 8)

# Younger matrix (partialRsq1 placed into lower triangle, symmetric)
mtrixplot <- matrix(NA_real_, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = TRUE)] <- gamresultsum.df_younger$partialRsq1[seq_len(min(length(gamresultsum.df_younger$partialRsq1), elementnum))]
Matrix.tmp <- mtrixplot
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
mat_df <- as.data.frame(Matrix.tmp)
mat_df$nodeid <- seq_len(ds.resolution)
mat_m <- melt_df(mat_df, id.vars = c("nodeid"))
mat_m$variable <- as.numeric(mat_m$variable)
mat_m$nodeid <- 0 - mat_m$nodeid
mat_m$value <- as.numeric(mat_m$value)
linerange_frame <- data.frame(
  x = c(0.5, 12 + 0.5), ymin = rep(-12 - 0.5, times = 2), ymax = rep(-0.5, times = 2),
  y = c(-0.5, -12 - 0.5), xmin = rep(0.5, times = 2), xmax = rep(12 + 0.5, times = 2)
)
RdBucol <- rev(brewer.pal(11, "RdBu"))
p_mat_y <- ggplot(mat_m) +
  geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
  scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-maxthr, maxthr), na.value = RdBucol[11]) +
  scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-maxthr, maxthr), na.value = RdBucol[11]) +
  geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
  geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
  geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
  ggtitle(label = "partialRsq") +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, angle = 315, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid.major = element_line(linewidth = 0),
    panel.grid.minor = element_line(linewidth = 1)
  )
ggsave(file.path(FigMatYoung, paste0("partialRsq_12net_delLM_CV", CVthr, ".tiff")), p_mat_y, dpi = 600, height = 18, width = 20, units = "cm", bg = "transparent")
ggsave(file.path(FigMatYoung, paste0("partialRsq_12net_delLM_CV", CVthr, ".pdf")), p_mat_y, dpi = 600, height = 18, width = 20, units = "cm", bg = "transparent")

# Older: partialRsq
gamresultsum.df_older <- old$results
computevar_o <- "partialRsq"
ct_o <- SCrankcorr(gamresultsum.df_older, computevar_o, ds.resolution, dsdata = FALSE)
message(sprintf("[RESULT] older %s: r=%.5f p=%.3g (n=%d)", computevar_o, ct_o$r.spearman[[1]], ct_o$p.spearman[[1]], nrow(gamresultsum.df_older)))

correlation.df2 <- SCrankcorr(gamresultsum.df_older, computevar_o, ds.resolution, dsdata = TRUE)
maxthr2 <- max(abs(correlation.df2$partialRsq), na.rm = TRUE)
p_sc_o <- ggplot(correlation.df2) +
  geom_point(aes(x = SCrank, y = partialRsq, color = partialRsq), size = 3) +
  geom_smooth(aes(x = SCrank, y = partialRsq), method = "lm", color = "black") +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-maxthr2, maxthr2)) +
  labs(x = "S-A connectional axis rank", y = NULL) +
  scale_y_continuous(breaks = c(0.00, 0.02, 0.04)) +
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
ggsave(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_n", ds.resolution, ".tiff")), p_sc_o, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_n", ds.resolution, ".pdf")), p_sc_o, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
if (save_svg) save_svg_if_available(file.path(FigCorrOld, paste0(computevar_o, "_SCrankcorr_n", ds.resolution, ".svg")), p_sc_o, width_cm = 10, height_cm = 8)

Matrix.tmp2 <- matrix(NA_real_, nrow = 12, ncol = 12)
Matrix.tmp2[lower.tri(Matrix.tmp2, diag = TRUE)] <- correlation.df2[[2]]
Matrix.tmp2[upper.tri(Matrix.tmp2)] <- t(Matrix.tmp2)[upper.tri(Matrix.tmp2)]
mat_df2 <- as.data.frame(Matrix.tmp2)
mat_df2$nodeid <- seq_len(12)
mat_m2 <- melt_df(mat_df2, id.vars = c("nodeid"))
mat_m2$variable <- as.numeric(mat_m2$variable)
mat_m2$nodeid <- 0 - mat_m2$nodeid
mat_m2$value <- as.numeric(mat_m2$value)
p_mat_o <- ggplot(mat_m2) +
  geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
  scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-maxthr2, maxthr2), na.value = "grey") +
  scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-maxthr2, maxthr2), na.value = "grey") +
  geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
  geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
  geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
  ggtitle(label = "partialRsq") +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, angle = 315, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid.major = element_line(linewidth = 0),
    panel.grid.minor = element_line(linewidth = 1)
  )
ggsave(file.path(FigMatOld, paste0("partialRsq_12net_delLM_CV", CVthr, ".tiff")), p_mat_o, dpi = 600, height = 18, width = 20, units = "cm", bg = "transparent")
ggsave(file.path(FigMatOld, paste0("partialRsq_12net_delLM_CV", CVthr, ".pdf")), p_mat_o, dpi = 600, height = 18, width = 20, units = "cm", bg = "transparent")
