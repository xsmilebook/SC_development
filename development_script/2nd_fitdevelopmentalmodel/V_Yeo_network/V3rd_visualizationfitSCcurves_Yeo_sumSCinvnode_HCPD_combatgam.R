## HCP-D (ComBat-GAM) | Yeo7/Yeo17 | Visualization of fitted SC curves
##
## Runnable version adapted from:
##   V3rd_visualizationfitSCcurves_Yeo_sumSCinvnode_HCPD.R
##
## Outputs (match original figure naming):
## - devcurve_Rsq_fit.ratio.{tiff,svg}
## - devcurve_SCrank_fit.Z.{tiff,svg}

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
force <- as.integer(if (!is.null(args$force)) args$force else 0L) == 1L

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
yeo <- as.integer(if (!is.null(args$yeo)) args$yeo else 17L)
if (!yeo %in% c(7L, 17L)) stop("Unsupported yeo resolution: ", yeo, " (expected 7 or 17)")

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam")
FigureFolder <- file.path(FigureRoot, paste0("CV", CVthr), paste0("Yeo", yeo, "_sumSCinvnode_fit"))
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "plotdata_generate.R"))

gamresult_files <- list.files(interfileFolder, pattern = paste0("^gamresults_Yeo\\d+_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE\\.rds$"), full.names = TRUE)
if (length(gamresult_files) != 1) stop("Expected one scaled gamresults file in: ", interfileFolder)
gamresultsum <- readRDS(gamresult_files[[1]])
elementnum <- nrow(gamresultsum)
Yeoresolution.delLM <- infer_resolution_from_edges(elementnum)

gammodel_files <- list.files(interfileFolder, pattern = paste0("^gammodel_Yeo", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE\\.rds$"), full.names = TRUE)
if (length(gammodel_files) != 1) stop("Expected one scaled gammodel file in: ", interfileFolder)
gammodelsum <- readRDS(gammodel_files[[1]])

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

out_fig1_tiff <- file.path(FigureFolder, "devcurve_Rsq_fit.ratio.tiff")
out_fig2_tiff <- file.path(FigureFolder, "devcurve_SCrank_fit.Z.tiff")
out_fig1_svg <- file.path(FigureFolder, "devcurve_Rsq_fit.ratio.svg")
out_fig2_svg <- file.path(FigureFolder, "devcurve_SCrank_fit.Z.svg")
out_plotdatasum <- file.path(interfileFolder, paste0("plotdatasum_scale_TRUE_Yeo", yeo, ".rds"))

if (!force && file.exists(out_fig1_tiff) && file.exists(out_fig2_tiff) && file.exists(out_fig1_svg) && file.exists(out_fig2_svg)) {
  message("[INFO] Yeo S3 outputs exist; skipping. Set --force=1 to re-run.")
  quit(save = "no", status = 0)
}

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

plotdatasum <- mclapply(seq_len(min(length(gammodelsum), elementnum)), plot_one, mc.cores = n_cores)
saveRDS(plotdatasum, out_plotdatasum)

Matrixds <- matrix(NA, nrow = Yeoresolution.delLM, ncol = Yeoresolution.delLM)
indexup <- upper.tri(Matrixds)
indexsave <- !indexup
Matrixds.SCrank <- Matrixds
for (x in 1:Yeoresolution.delLM) {
  for (y in 1:Yeoresolution.delLM) {
    Matrixds.SCrank[x, y] <- (x + y)^2 + (x - y)^2
  }
}
Matrixds.SCrank[indexup] <- NA
Matrixds.SCrank[indexsave] <- rank(Matrixds.SCrank[indexsave], ties.method = "average")
parcel_all <- paste0("SC.", seq_len(elementnum), "_h")
SCrank_map <- setNames(Matrixds.SCrank[indexsave], parcel_all)

ok_plot_idx <- which(!vapply(plotdatasum, is.null, logical(1)))
if (length(ok_plot_idx) == 0) stop("All plotdata_generate() calls failed; see warnings above.")

plotdatasum.df <- dplyr::bind_rows(lapply(ok_plot_idx, function(i) {
  tmp <- as.data.frame(plotdatasum[[i]])
  sc_label <- gamresultsum$parcel[[i]]
  if (!is.null(sc_label) && sc_label %in% names(tmp)) {
    tmp <- tmp[, setdiff(names(tmp), sc_label), drop = FALSE]
  }
  tmp$SC_label <- sc_label
  tmp$SCrank <- unname(SCrank_map[[sc_label]])
  tmp$PartialRsq <- as.numeric(gamresultsum$partialRsq[[i]])
  tmp$meanderiv2 <- as.numeric(gamresultsum$meanderv2[[i]])
  tmp
}))

## Plot 1: ratio trajectories colored by partial R^2
lmthr <- max(abs(gamresultsum$partialRsq), na.rm = TRUE)
p1 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.ratio, group = SC_label, color = PartialRsq), linewidth = 1.5, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-lmthr, lmthr), guide = "none") +
  labs(x = "Age (years)", y = "SC strength (ratio)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 22, hjust = 0.5),
    legend.position = "none"
  )
ggsave(out_fig1_tiff, p1, width = 15, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig1_svg, p1, dpi = 600, width = 15, height = 12, units = "cm", bg = "transparent")

## Plot 2: z-scored trajectories colored by S-A rank (SCrank)
SC_label_derv2_order <- gamresultsum$parcel[order(gamresultsum$meanderv2)]
plotdatasum.df$SC_label2 <- factor(plotdatasum.df$SC_label, levels = SC_label_derv2_order)
p2 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.Z, group = SC_label2, color = SCrank), linewidth = 0.8, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, guide = "none") +
  labs(x = "Age (years)", y = "SC strength (z-score)") +
  scale_y_continuous(breaks = c(-1.5, 0.0, 1.5)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 23, color = "black"),
    axis.title = element_text(size = 23),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )
ggsave(out_fig2_tiff, p2, width = 15, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig2_svg, p2, dpi = 600, width = 16, height = 16, units = "cm", bg = "transparent")

