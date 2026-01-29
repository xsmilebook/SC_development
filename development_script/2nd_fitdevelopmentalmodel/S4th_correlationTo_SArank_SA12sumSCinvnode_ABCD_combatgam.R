## Rscript version of:
##   S4th_correlationTo_SArank_SA12sumSCinvnode_ABCD.Rmd
##
## This script performs Spearman correlations between GAMM-derived developmental
## statistics and the S-A connectional axis rank, and generates scatter plots.
##
## Default inputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds
## - wd/interdataFolder_ABCD/average_EuclideanDistance_12.csv
##
## Outputs:
## - outputs/results/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/SCrank_correlation_summary.csv
## - outputs/figures/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/correlation_sumSCinvnode_SCrank/*.tiff + *.pdf
## - (optional) outputs/figures/.../Matrix12_sumSCinvnode_gamstats_Age8_22/*.tiff

rm(list = ls())

library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
library(ggplot2)

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
make_matrix_graphs <- as.integer(if (!is.null(args$make_matrix_graphs)) args$make_matrix_graphs else 0L) == 1L

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "abcd", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "abcd", "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "abcd", "combat_gam", paste0("CV", CVthr))

dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

gamresult <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
if (!is.data.frame(gamresult) || !"parcel" %in% names(gamresult)) {
  stop("Invalid gamresults input: ", interfileFolder)
}
if (!"bootstrap_pvalue" %in% names(gamresult)) {
  stop("Missing bootstrap_pvalue in gamresults; check S1 outputs.")
}
gamresult$pfdr <- p.adjust(gamresult$bootstrap_pvalue, method = "fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "abcd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)
SCdata <- readRDS(input_rds)

euclid_csv <- if (!is.null(args$euclid_csv)) {
  args$euclid_csv
} else {
  file.path(project_root, "wd", "interdataFolder_ABCD", "average_EuclideanDistance_12.csv")
}
if (!file.exists(euclid_csv)) stop("Missing euclid_csv: ", euclid_csv)
EucDistance <- read.csv(euclid_csv)

message(sum(gamresult$sig), " edges have significant developmental effects.")

out_summary <- file.path(resultFolder, "SCrank_correlation_summary.csv")
if (!force && file.exists(out_summary)) {
  message("[INFO] S4 summary exists, skipping: ", out_summary, " (set --force=1 to re-run)")
  quit(save = "no", status = 0)
}

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (length(sc_cols) == 0) stop("No SC.* columns found in input_rds: ", input_rds)
meanSC <- colMeans(SCdata[, sc_cols, drop = FALSE], na.rm = TRUE)

parcel_all <- paste0("SC.", seq_len(elementnum), "_h")
meanSC_map <- setNames(meanSC, parcel_all)

if (!"Edistance" %in% names(EucDistance)) stop("Euclidean distance CSV missing column 'Edistance': ", euclid_csv)
if (nrow(EucDistance) == elementnum) {
  Edist_map <- setNames(EucDistance$Edistance, parcel_all)
} else if ("parcel" %in% names(EucDistance)) {
  Edist_map <- setNames(EucDistance$Edistance, as.character(EucDistance$parcel))
} else {
  stop("Euclidean distance rows do not match elementnum and no 'parcel' column to align: ", euclid_csv)
}

meanSC_aligned <- meanSC_map[gamresult$parcel]
corr.test(meanSC_aligned, gamresult$partialRsq)

gamresult <- within(gamresult, {
  partialRsq2 <- partialRsq
  partialRsq2[partialRsq2 > mean(partialRsq2, na.rm = TRUE) + 3 * sd(partialRsq2, na.rm = TRUE) |
    partialRsq2 < mean(partialRsq2, na.rm = TRUE) - 3 * sd(partialRsq2, na.rm = TRUE)] <- NA
})

## compute correlations to SC rank
computevar.list <- c("partialRsq2", "meanderv2")
SCrank_correlation <- do.call(
  rbind,
  lapply(computevar.list, function(computevar) SCrankcorr(gamresult, computevar, ds.resolution, dsdata = FALSE))
)

## Validation: Control for Euclidean distance
gamresult$EucDistance <- unname(Edist_map[gamresult$parcel])
fit_pr <- lm(partialRsq2 ~ EucDistance, data = gamresult, na.action = na.exclude)
fit_md <- lm(meanderv2 ~ EucDistance, data = gamresult, na.action = na.exclude)
gamresult$partialRsq_control_distance <- as.numeric(residuals(fit_pr))
gamresult$meanderv2_control_distance <- as.numeric(residuals(fit_md))
SCrank_correlation <- rbind(
  SCrank_correlation,
  SCrankcorr(gamresult, "partialRsq_control_distance", ds.resolution, dsdata = FALSE),
  SCrankcorr(gamresult, "meanderv2_control_distance", ds.resolution, dsdata = FALSE)
)
write.csv(SCrank_correlation, out_summary, row.names = FALSE)

## scatter plots
FigCorrFolder <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank")
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

get_scatter_style <- function(computevar, cvthr) {
  if (computevar %in% c("partialRsq2", "partialRsq_control_distance")) {
    if (cvthr == 75) {
      return(list(
        theme = theme(
          axis.text = element_text(size = 24.3, color = "black"),
          axis.title = element_text(size = 24.3),
          aspect.ratio = 0.8,
          axis.line = element_line(linewidth = 0.6),
          axis.ticks = element_line(linewidth = 0.6),
          plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "none"
        ),
        svg_width = 17.5,
        svg_height = 15
      ))
    }
    return(list(
      theme = theme(
        axis.text = element_text(size = 24, color = "black"),
        axis.title = element_text(size = 24),
        aspect.ratio = 1.05,
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size = 15, hjust = 0.5, vjust = 0),
        plot.subtitle = element_text(size = 21, hjust = 0.9, vjust = -6),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
      ),
      svg_width = 17,
      svg_height = 14
    ))
  }

  if (cvthr == 25) {
    return(list(
      theme = theme(
        axis.text = element_text(size = 23.2, color = "black"),
        axis.title = element_text(size = 23.2),
        aspect.ratio = 1,
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size = 15, hjust = 0.5, vjust = 0),
        plot.subtitle = element_text(size = 15, hjust = 0.1, vjust = -6),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
      ),
      svg_width = 13,
      svg_height = 13
    ))
  }

  list(
    theme = theme(
      axis.text = element_text(size = 26.4, color = "black"),
      axis.title = element_text(size = 26.4),
      aspect.ratio = 0.74,
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),
      plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none"
    ),
    svg_width = 17.5,
    svg_height = 15
  )
}

get_colorbar_prob <- function(computevar) {
  if (computevar == "partialRsq2") return(0.4)
  if (computevar == "meanderv2") return(0.46)
  if (computevar == "partialRsq_control_distance") return(0.4)
  if (computevar == "meanderv2_control_distance") return(0.53)
  0.4
}

plot_one_scatter <- function(computevar, ylab) {
  df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  names(df) <- c("SCrank", computevar)
  style <- get_scatter_style(computevar, CVthr)
  colorbar_prob <- get_colorbar_prob(computevar)
  colorbar_vals <- colorbarvalues(df[[computevar]], colorbar_prob)

  p <- ggplot(df) +
    geom_point(aes(x = SCrank, y = .data[[computevar]], color = .data[[computevar]]), size = 5) +
    geom_smooth(aes(x = SCrank, y = .data[[computevar]]), method = "lm", color = "black", linewidth = 1.2) +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, values = colorbar_vals) +
    labs(x = "S-A connectional axis rank", y = ylab) +
    theme_classic() + style$theme

  if (computevar == "meanderv2") {
    p <- p + scale_y_continuous(breaks = c(-0.005, 0, 0.005, 0.010), labels = c(-5, 0, 5, 10))
  }

  attr(p, "svg_width") <- style$svg_width
  attr(p, "svg_height") <- style$svg_height
  p
}

scatter_targets <- list(
  partialRsq2 = expression("Age effect (partial "*italic("R")^"2"*")"),
  meanderv2 = "Second derivative",
  partialRsq_control_distance = expression("Age effect (partial "*italic("R")^"2"*")"),
  meanderv2_control_distance = "Second derivative"
)

for (nm in names(scatter_targets)) {
  if (!nm %in% names(gamresult)) next
  p <- plot_one_scatter(nm, scatter_targets[[nm]])
  pdf_path <- file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".pdf"))
  ggsave(pdf_path, p, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")

  if (.Platform$OS.type == "windows") {
    svg_width <- attr(p, "svg_width")
    svg_height <- attr(p, "svg_height")
    svg_path <- file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".svg"))
    ggsave(svg_path, p, dpi = 600, width = svg_width, height = svg_height, units = "cm", bg = "transparent")
  }
}

## matrix graphs (optional; heavy)
if (make_matrix_graphs) {
  Matrix.tmp <- matrix(NA, nrow = 12, ncol = 12)
  computevar.list2 <- c("partialRsq2", "meanderv2", "partialRsq_control_distance", "meanderv2_control_distance")

  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5), ymin = rep(-12 - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5), xmin = rep(0.5, times = 2), xmax = rep(12 + 0.5, times = 2)
  )

  SCrank_correlation.df <- mclapply(seq_along(computevar.list2), function(i) {
    computevar <- computevar.list2[[i]]
    SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  }, mc.cores = min(6L, n_cores))

  colorbar.prob <- c(
    0.5,
    abs(min(gamresult$meanderv2, na.rm = TRUE)) / (max(gamresult$meanderv2, na.rm = TRUE) - min(gamresult$meanderv2, na.rm = TRUE)),
    0.5,
    abs(min(gamresult$meanderv2_control_distance, na.rm = TRUE)) /
      (max(gamresult$meanderv2_control_distance, na.rm = TRUE) - min(gamresult$meanderv2_control_distance, na.rm = TRUE))
  )

  FigMatrixFolder <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats_Age8_22")
  dir.create(FigMatrixFolder, showWarnings = FALSE, recursive = TRUE)

  for (i in seq_along(computevar.list2)) {
    computevar <- computevar.list2[[i]]
    Matrix.tmp[lower.tri(Matrix.tmp, diag = TRUE)] <- SCrank_correlation.df[[i]][, 2]
    Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
    colnames(Matrix.tmp) <- seq(1, 12)
    rownames(Matrix.tmp) <- seq(1, 12)

    corrplot(Matrix.tmp, method = "color", type = "lower", tl.col = "black", tl.cex = 0.8)
    ggsave(file.path(FigMatrixFolder, paste0("Matrix12_sumSCinvnode_", computevar, ".tiff")), width = 12, height = 10, units = "cm", dpi = 600)
  }
}
