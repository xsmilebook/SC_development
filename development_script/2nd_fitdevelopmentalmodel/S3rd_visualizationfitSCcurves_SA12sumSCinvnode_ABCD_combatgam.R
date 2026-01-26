## Rscript version of:
##   S3rd_visualizationfitSCcurves_SA12sumSCinvnode_ABCD.Rmd
##
## This script generates fitted values from scaled GAMM models and produces
## developmental trajectory plots for ABCD (SA12, sumSCinvnode), using outputs
## from the ComBat-GAM pipeline.
##
## Default inputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gammodel78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/results/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/derivative.df78_CV75.rds
##
## Outputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/plotdatasum_scale_TRUE_SA12.rds
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/SA12_10.csv
## - outputs/figures/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/** (tiff+pdf)

rm(list = ls())

library(parallel)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(patchwork)

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

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "abcd", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "abcd", "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "abcd", "combat_gam")

dir.create(interfileFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "plotdata_generate.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

## load data
gamresultsum.SAorder.delLM <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
gammodelsum <- readRDS(file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
derivative <- readRDS(file.path(resultFolder, paste0("derivative.df", elementnum, "_CV", CVthr, ".rds")))

ensure_sig_derivative_fdr <- function(derivative_df) {
  if (!is.data.frame(derivative_df) || nrow(derivative_df) == 0) return(derivative_df)
  if ("significant.derivative_fdr" %in% names(derivative_df) &&
    length(derivative_df$significant.derivative_fdr) == nrow(derivative_df)) {
    return(derivative_df)
  }

  if ("derivative" %in% names(derivative_df) && "significance_pvalue_fdr" %in% names(derivative_df)) {
    sig_val <- as.numeric(derivative_df$derivative)
    sig_mask <- as.logical(derivative_df$significance_pvalue_fdr)
    sig_val[which(is.na(sig_mask) | !sig_mask)] <- NA_real_
    derivative_df$significant.derivative_fdr <- sig_val
    return(derivative_df)
  }

  if ("derivative" %in% names(derivative_df) && "significant" %in% names(derivative_df)) {
    sig_val <- as.numeric(derivative_df$derivative)
    sig_mask <- as.logical(derivative_df$significant)
    sig_val[which(is.na(sig_mask) | !sig_mask)] <- NA_real_
    derivative_df$significant.derivative_fdr <- sig_val
    return(derivative_df)
  }

  if ("significant.derivative" %in% names(derivative_df)) {
    sig_val <- as.numeric(derivative_df$significant.derivative)
    derivative_df$significant.derivative_fdr <- sig_val
    return(derivative_df)
  }

  warning("Missing significance columns in derivative.df; setting significant.derivative_fdr to NA.")
  derivative_df$significant.derivative_fdr <- rep(NA_real_, nrow(derivative_df))
  derivative_df
}

derivative <- ensure_sig_derivative_fdr(derivative)

FigureFolder <- file.path(FigureRoot, paste0("CV", CVthr))
FigureFolder_SCfit <- file.path(FigureFolder, "SA12_sumSCinvnode_fit")
FigureFolder_SCdecile <- file.path(FigureFolder, "SA12_decile_sumSCinvnode_fit")
dir.create(FigureFolder_SCfit, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder_SCdecile, showWarnings = FALSE, recursive = TRUE)

out_plotdatasum_list <- file.path(interfileFolder, "plotdatasum_scale_TRUE_SA12.rds")
out_sa12_csv <- file.path(interfileFolder, "SA12_10.csv")
out_fig1 <- file.path(FigureFolder_SCfit, "devcurve_Rsq_fit.ratio.tiff")
out_fig2 <- file.path(FigureFolder_SCfit, "devcurve_meanderv2_fit.Z.tiff")

if (!force && file.exists(out_plotdatasum_list) && file.exists(out_fig1) && file.exists(out_fig2)) {
  message("[INFO] S3 outputs exist; skipping S3. Set --force=1 to re-run.")
  quit(save = "no", status = 0)
}

## generate fitted values for developmental trajectories
if (nrow(gamresultsum.SAorder.delLM) == 0 || length(gammodelsum) == 0) {
  stop("Empty inputs: gamresults or gammodels. Check upstream S1/S2 outputs under: ", interfileFolder)
}
if (length(gammodelsum) != nrow(gamresultsum.SAorder.delLM)) {
  message(
    "[WARN] gammodels length (", length(gammodelsum),
    ") != gamresults rows (", nrow(gamresultsum.SAorder.delLM),
    "). Using min() and skipping missing models."
  )
}
n_edges <- min(length(gammodelsum), nrow(gamresultsum.SAorder.delLM))

plot_one <- function(idx) {
  modobj <- gammodelsum[[idx]]
  if (is.null(modobj)) return(NULL)
  tryCatch(
    plotdata_generate(modobj, "age"),
    error = function(e) {
      message("[WARN] plotdata_generate failed for idx=", idx, ": ", conditionMessage(e))
      NULL
    }
  )
}

plotdatasum <- mclapply(seq_len(n_edges), plot_one, mc.cores = n_cores)
saveRDS(plotdatasum, out_plotdatasum_list)

## SA12 index & SC rank
Matrix12 <- matrix(NA, nrow = 12, ncol = 12)
indexup12 <- upper.tri(Matrix12)
indexsave12 <- !indexup12
Matrix12.SCrank <- Matrix12
for (x in 1:12) {
  for (y in 1:12) {
    Matrix12.SCrank[x, y] <- x^2 + y^2
  }
}
Matrix12.SCrank[indexup12] <- NA
Matrix12.SCrank[indexsave12] <- rank(Matrix12.SCrank[indexsave12], ties.method = "average")

parcel_all <- paste0("SC.", seq_len(elementnum), "_h")
SCrank_map <- setNames(Matrix12.SCrank[indexsave12], parcel_all)

ok_plot_idx <- which(!vapply(plotdatasum, is.null, logical(1)))
if (length(ok_plot_idx) == 0) stop("All plotdata_generate() calls failed; see warnings above.")

plotdatasum.df <- dplyr::bind_rows(lapply(ok_plot_idx, function(i) {
  tmp <- as.data.frame(plotdatasum[[i]])
  sc_label <- gamresultsum.SAorder.delLM$parcel[[i]]
  if (!is.null(sc_label) && sc_label %in% names(tmp)) {
    tmp <- tmp[, setdiff(names(tmp), sc_label), drop = FALSE]
  }
  tmp$SC_label <- sc_label
  tmp$SCrank <- unname(SCrank_map[[sc_label]])
  tmp$PartialRsq <- gamresultsum.SAorder.delLM$partialRsq[[i]]
  tmp$meanderv2 <- gamresultsum.SAorder.delLM$meanderv2[[i]]
  tmp
}))

## Save SA12 decile labels
SA12_10 <- data.frame(SCrank = Matrix12.SCrank[indexsave12]) %>%
  mutate(decile = ntile(SCrank, 10))
SA12_10$SC_label <- parcel_all
write.csv(SA12_10, out_sa12_csv, row.names = FALSE)

## Plots: 78 developmental trajectories (fit.ratio colored by partial R^2)
lmthr <- max(abs(gamresultsum.SAorder.delLM$partialRsq))
p1 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.ratio, group = SC_label, color = PartialRsq), linewidth = 0.8, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-lmthr, lmthr), guide = "none") +
  labs(x = "Age (years)", y = "SC strength (ratio)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20.5, color = "black"),
    axis.title = element_text(size = 20.5, color = "black"),
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "none"
  )
ggsave(file.path(FigureFolder_SCfit, "devcurve_Rsq_fit.ratio.tiff"), p1, width = 20, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder_SCfit, "devcurve_Rsq_fit.ratio.pdf"), p1, dpi = 600, width = 15, height = 15, units = "cm", bg = "transparent")

## Plots: 78 developmental trajectories (fit.Z colored by mean 2nd derivative)
colorbarvalues.meanderv2 <- colorbarvalues(
  plotdatasum.df$meanderv2,
  abs(min(plotdatasum.df$meanderv2)) / (max(plotdatasum.df$meanderv2) - min(plotdatasum.df$meanderv2))
)
SC_label_derv2_order <- gamresultsum.SAorder.delLM$parcel[order(gamresultsum.SAorder.delLM$meanderv2)]
plotdatasum.df$SC_label2 <- factor(plotdatasum.df$SC_label, levels = SC_label_derv2_order)
p2 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.Z, group = SC_label2, color = meanderv2), linewidth = 0.8, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", values = colorbarvalues.meanderv2, direction = -1, guide = "none") +
  labs(x = "Age (years)", y = "SC strength (z-score)") +
  scale_y_continuous(breaks = c(-1.5, 0.0, 1.5)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20.5, color = "black"),
    axis.title = element_text(size = 20.5, color = "black"),
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "none"
  )
ggsave(file.path(FigureFolder_SCfit, "devcurve_meanderv2_fit.Z.tiff"), p2, width = 20, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder_SCfit, "devcurve_meanderv2_fit.Z.pdf"), p2, dpi = 600, width = 15, height = 15, units = "cm", bg = "transparent")

## Decile trajectories
plotdatasum.df.label <- merge(plotdatasum.df, SA12_10, by = "SC_label", all.x = TRUE)
plotdatasum.df.decile <- plotdatasum.df.label %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(fit), SCranktype_order = mean(decile), .groups = "drop")
plotdatasum.df.decile <- plotdatasum.df.decile %>%
  group_by(decile) %>%
  mutate(fit.Z = scale(fit.avg))

p3 <- ggplot(data = plotdatasum.df.decile, aes(x = age, y = fit.Z, group = decile, color = decile)) +
  geom_line(linewidth = 1.3, alpha = 1) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1) +
  labs(x = "Age (years)", y = "SC strength (z-score)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22.5, color = "black"),
    axis.title = element_text(size = 22.5, color = "black"),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )

ggsave(file.path(FigureFolder_SCdecile, "devcurve_SCrank_fit.Z_SCtype10.tiff"), p3, dpi = 600, width = 20, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder_SCdecile, "devcurve_SCrank_fit.Z_SCtype10.pdf"), p3, dpi = 600, width = 15, height = 13, units = "cm", bg = "transparent")
