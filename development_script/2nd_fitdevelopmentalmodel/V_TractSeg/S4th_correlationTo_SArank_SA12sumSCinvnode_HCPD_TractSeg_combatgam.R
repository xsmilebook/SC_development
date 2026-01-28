## HCP-D (ComBat-GAM) | TractSeg (major-bundle) | Correlation to S-A rank + matrix graphs
##
## Runnable version adapted from:
##   V_TractSeg/S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD.R
##
## Figures (match original naming):
## - TractSeg/correlation_sumSCinvnode_SCrank/meanmeanderv2_SCrankcorr_n12.(tiff,pdf)
## - TractSeg/Matrix12_sumSCinvnode_gamstats_Age8_22/{computevar}_12net_delLM_CV75.tiff

rm(list = ls())

library(tidyverse)
library(parallel)
library(psych)
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

print_spearman_summary <- function(summary_df, prefix = "[INFO]") {
  required <- c("ds.resolution", "Interest.var", "r.spearman", "p.spearman")
  missing <- setdiff(required, names(summary_df))
  if (length(missing) > 0) {
    message(prefix, " SCrank summary missing columns: ", paste(missing, collapse = ", "))
    return(invisible(NULL))
  }

  summary_df <- summary_df %>%
    mutate(
      ds.resolution = as.integer(ds.resolution),
      Interest.var = as.character(Interest.var),
      r.spearman = as.numeric(r.spearman),
      p.spearman = as.numeric(p.spearman)
    )

  ds <- unique(summary_df$ds.resolution)
  ds_label <- if (length(ds) == 1) as.character(ds[[1]]) else paste(ds, collapse = ",")
  message(prefix, " SCrank Spearman summary (ds.resolution=", ds_label, "):")

  for (i in seq_len(nrow(summary_df))) {
    message(prefix, "  ", summary_df$Interest.var[[i]], ": r=", sprintf("%.4f", summary_df$r.spearman[[i]]), ", p=", format(summary_df$p.spearman[[i]], digits = 3, scientific = TRUE))
  }
  invisible(NULL)
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

make_matrix_graphs <- as.integer(if (!is.null(args$make_matrix_graphs)) args$make_matrix_graphs else 1L) == 1L

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "tractseg", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "tractseg", "combat_gam", paste0("CV", CVthr)
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "hcpd", "tractseg", "combat_gam", paste0("CV", CVthr), "TractSeg")
FigCorrFolder <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank")
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

gamresult <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE_TractSeg.rds")))
if (!is.data.frame(gamresult) || !"parcel" %in% names(gamresult)) stop("Invalid gamresults under: ", interfileFolder)

gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)

## boundary handling (match original HCP-D constants)
gamresult$increase.onset[gamresult$sig == FALSE] <- NA
gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2, 2) == 8.08] <- NA
gamresult$increase.offset[gamresult$sig == FALSE] <- NA
gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2, 2) == 21.92] <- NA
gamresult$peak.change[gamresult$sig == FALSE] <- NA
gamresult$peak.increase.change[gamresult$sig == FALSE] <- NA

out_summary <- file.path(resultFolder, "SCrank_correlation_summary.csv")
if (force || !file.exists(out_summary)) {
  computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change", "meanderv2")
  SCrank_correlation <- do.call(
    rbind,
    lapply(computevar.list, function(computevar) SCrankcorr(gamresult, computevar, ds.resolution, dsdata = FALSE))
  )
  write.csv(SCrank_correlation, out_summary, row.names = FALSE)
}
tryCatch({
  SCrank_correlation <- read.csv(out_summary, stringsAsFactors = FALSE)
  print_spearman_summary(SCrank_correlation)
}, error = function(e) {
  message("[WARN] Failed to read/print S4 summary: ", out_summary, " | ", conditionMessage(e))
})

## scatter plot: meanderv2 (match original)
out_meanderv2_tiff <- file.path(FigCorrFolder, paste0("meanmeanderv2_SCrankcorr_n", ds.resolution, ".tiff"))
out_meanderv2_pdf <- file.path(FigCorrFolder, paste0("meanmeanderv2_SCrankcorr_n", ds.resolution, ".pdf"))

if (force || !file.exists(out_meanderv2_tiff) || !file.exists(out_meanderv2_pdf)) {
  correlation.df <- SCrankcorr(gamresult, "meanderv2", ds.resolution, dsdata = TRUE)
  mytheme <- theme(
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
  p <- ggplot(data = correlation.df) +
    geom_point(aes(x = SCrank, y = meanderv2, color = SCrank), size = 5) +
    geom_smooth(aes(x = SCrank, y = meanderv2), linewidth = 2, method = "lm", color = "black") +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, guide = "none") +
    labs(x = "S-A connectional axis rank", y = "Second derivative") +
    scale_y_continuous(breaks = c(-0.004, -0.002, 0, 0.002), labels = c(-4, -2, 0, 2)) +
    theme_classic() + mytheme
  ggsave(out_meanderv2_tiff, p, width = 13, height = 12, units = "cm", bg = "transparent")
  ggsave(out_meanderv2_pdf, p, dpi = 600, width = 13, height = 12, units = "cm", bg = "transparent")
}

## matrix graphs (match original naming)
if (make_matrix_graphs) {
  FigMatrixFolder <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats_Age8_22")
  dir.create(FigMatrixFolder, showWarnings = FALSE, recursive = TRUE)

  computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2")
  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5), ymin = rep(-12 - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5), xmin = rep(0.5, times = 2), xmax = rep(12 + 0.5, times = 2)
  )

  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
  n_cores <- max(1L, n_cores)

  SCrank_correlation.df <- mclapply(seq_along(computevar.list), function(i) {
    computevar <- computevar.list[[i]]
    SCrankcorr(gamresult, computevar, 12, dsdata = TRUE)
  }, mc.cores = min(6L, n_cores))

  colorbar.prob <- c(0.5, 0.4, 0.6, 0.5, 0.5)

  for (i in seq_along(computevar.list)) {
    computevar <- computevar.list[[i]]
    filename <- file.path(FigMatrixFolder, paste0(computevar, "_12net_delLM_CV", CVthr, ".tiff"))
    if (!force && file.exists(filename)) next

    Matrix.tmp <- matrix(NA, nrow = 12, ncol = 12)
    Matrix.tmp[lower.tri(Matrix.tmp, diag = TRUE)] <- SCrank_correlation.df[[i]][, 2]
    Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
    colnames(Matrix.tmp) <- seq(1, 12)
    rownames(Matrix.tmp) <- seq(1, 12)
    matrixtmp.df <- as.data.frame(Matrix.tmp)
    matrixtmp.df$nodeid <- seq(1, 12)
    matrixtmp.df.melt <- melt(matrixtmp.df, id.vars = c("nodeid"))
    matrixtmp.df.melt$variable <- as.numeric(matrixtmp.df.melt$variable)
    matrixtmp.df.melt$nodeid <- 0 - matrixtmp.df.melt$nodeid
    matrixtmp.df.melt$value <- as.numeric(matrixtmp.df.melt$value)

    colorbarvalues.tmp <- colorbarvalues(SCrank_correlation.df[[i]][, 2], colorbar.prob[[i]])

    if (computevar == "partialRsq") {
      lmthr <- max(abs(gamresult$partialRsq), na.rm = TRUE)
      Matrix.tmp.sig <- matrix(NA, nrow = 12, ncol = 12)
      Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = TRUE)] <- (gamresult$pfdr < 0.05)
      Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
      colnames(Matrix.tmp.sig) <- seq(1, 12)
      rownames(Matrix.tmp.sig) <- seq(1, 12)
      matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
      matrixtmp.df.sig$nodeid <- seq(1, 12)
      matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig, id.vars = c("nodeid"))
      matrixtmp.df.sig.melt$variable <- as.numeric(matrixtmp.df.sig.melt$variable)
      matrixtmp.df.sig.melt$nodeid <- 0 - matrixtmp.df.sig.melt$nodeid
      matrixtmp.df.sig.melt$value <- as.numeric(matrixtmp.df.sig.melt$value)
      matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value == 0), ]

      Fig <- ggplot(data = matrixtmp.df.melt) +
        geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
        scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "grey") +
        scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "grey") +
        geom_text(data = matrixtmp.df.sig.melt, aes(x = variable, y = nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size = 6) +
        geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
        geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
        geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
        ggtitle(label = computevar) + labs(x = NULL, y = NULL) +
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
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 1),
          plot.background = element_rect(fill = "transparent", color = NA)
        )
    } else {
      Fig <- ggplot(data = matrixtmp.df.melt) +
        geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
        scale_fill_distiller(type = "seq", palette = "RdBu", values = colorbarvalues.tmp, na.value = "grey") +
        scale_color_distiller(type = "seq", palette = "RdBu", values = colorbarvalues.tmp, na.value = "grey") +
        geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
        geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
        geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
        ggtitle(label = computevar) + labs(x = NULL, y = NULL) +
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
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 1),
          plot.background = element_rect(fill = "transparent", color = NA)
        )
    }

    ggsave(filename, Fig, height = 18, width = 20, units = "cm", dpi = 600, bg = "transparent")
  }
}
