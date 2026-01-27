## HCP-D (ComBat-GAM) | TractSeg (major-bundle) | Visualization of fitted SC curves
##
## Adapted to be runnable (project-root relative) and to match the original
## TractSeg figure outputs (two trajectory plots; tiff+svg).

rm(list = ls())

library(parallel)
library(tidyverse)
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
  "hcpd", "tractseg", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "tractseg", "combat_gam", paste0("CV", CVthr)
)

FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "hcpd", "tractseg", "combat_gam")
FigureFolder <- file.path(FigureRoot, paste0("CV", CVthr), "TractSeg", "SA12_sumSCinvnode_fit")
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "plotdata_generate.R"))

gamresult <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE_TractSeg.rds")))
gammodelsum <- readRDS(file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE_TractSeg.rds")))
if (nrow(gamresult) == 0 || length(gammodelsum) == 0) {
  stop("Empty inputs: gamresults or gammodels. Check upstream outputs under: ", interfileFolder)
}

out_fig1_tiff <- file.path(FigureFolder, "devcurve_Rsq_fit.ratio.tiff")
out_fig1_svg <- file.path(FigureFolder, "devcurve_Rsq_fit.ratio.svg")
out_fig2_tiff <- file.path(FigureFolder, "devcurve_meanderv2_fit.Z.tiff")
out_fig2_svg <- file.path(FigureFolder, "devcurve_meanderv2_fit.Z.svg")

if (!force && file.exists(out_fig1_tiff) && file.exists(out_fig1_svg) && file.exists(out_fig2_tiff) && file.exists(out_fig2_svg)) {
  message("[INFO] TractSeg S3 outputs exist; skipping. Set --force=1 to re-run.")
  quit(save = "no", status = 0)
}

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

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

plotdatasum <- mclapply(seq_len(min(length(gammodelsum), nrow(gamresult))), plot_one, mc.cores = n_cores)
ok_plot_idx <- which(!vapply(plotdatasum, is.null, logical(1)))
if (length(ok_plot_idx) == 0) stop("All plotdata_generate() calls failed; see warnings above.")

plotdatasum.df <- dplyr::bind_rows(lapply(ok_plot_idx, function(i) {
  tmp <- as.data.frame(plotdatasum[[i]])
  sc_label <- gamresult$parcel[[i]]
  if (!is.null(sc_label) && sc_label %in% names(tmp)) {
    tmp <- tmp[, setdiff(names(tmp), sc_label), drop = FALSE]
  }
  tmp$SC_label <- sc_label
  tmp$PartialRsq <- as.numeric(gamresult$partialRsq[[i]])
  tmp$meanderiv2 <- as.numeric(gamresult$meanderv2[[i]])
  tmp
}))

## Plot 1: ratio trajectories colored by partial R^2
lmthr <- max(abs(gamresult$partialRsq), na.rm = TRUE)
p1 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.ratio, group = SC_label, color = PartialRsq), linewidth = 0.8, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-lmthr, lmthr), guide = "none") +
  labs(x = "Age (years)", y = "SC strength (ratio)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    aspect.ratio = 0.9,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "none"
  )
ggsave(out_fig1_tiff, p1, width = 20, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig1_svg, p1, dpi = 600, width = 15, height = 15, units = "cm", bg = "transparent")

## Plot 2: z-scored trajectories colored by mean 2nd derivative
SC_label_derv2_order <- gamresult$parcel[order(gamresult$meanderv2)]
plotdatasum.df$SC_label2 <- factor(plotdatasum.df$SC_label, levels = SC_label_derv2_order)
lmthr2 <- max(abs(gamresult$meanderv2), na.rm = TRUE)
p2 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.Z, group = SC_label2, color = meanderiv2), linewidth = 0.8, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-lmthr2, lmthr2), guide = "none") +
  labs(x = "Age (years)", y = "SC strength (z-score)") +
  scale_y_continuous(breaks = c(-1.5, 0.0, 1.5)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    aspect.ratio = 0.9,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "none"
  )
ggsave(out_fig2_tiff, p2, width = 20, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig2_svg, p2, dpi = 600, width = 15, height = 15, units = "cm", bg = "transparent")

