## HCP-D (ComBat-GAM) | Yeo7/Yeo17 | Correlation to S-A rank + matrix graphs
##
## Runnable version adapted from:
##   V4th_correlationTo_SArank_Yeo_sumSCinvnode_HCPD.R
##
## Figures (match original naming):
## - correlation_sumSCinvnode_SCrank/meanpartialRsq_SCrankcorr_n{ds}.(tiff,pdf)
## - correlation_sumSCinvnode_SCrank/meanmeanderv2_2_SCrankcorr_n{ds}.(tiff,pdf)
## - Matrix{ds}_sumSCinvnode_gamstats_Age8_22/{computevar}_{ds}net_delLM_CV75.tiff

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
is_windows <- .Platform$OS.type == "windows"
is_interactive <- interactive()
skip_compute_on_windows <- as.integer(if (!is.null(args$skip_compute_on_windows)) args$skip_compute_on_windows else 1L) == 1L
skip_compute <- is_windows && skip_compute_on_windows

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
yeo <- as.integer(if (!is.null(args$yeo)) args$yeo else 17L)
if (!yeo %in% c(7L, 17L)) stop("Unsupported yeo resolution: ", yeo, " (expected 7 or 17)")

make_matrix_graphs <- as.integer(if (!is.null(args$make_matrix_graphs)) args$make_matrix_graphs else 1L) == 1L
if (is_windows && is_interactive && is.null(args$make_matrix_graphs)) {
  make_matrix_graphs <- FALSE
}

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "hcpd", "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr))
FigCorrFolder <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

gamresult_files <- list.files(interfileFolder, pattern = paste0("^gamresults_Yeo\\d+_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE\\.rds$"), full.names = TRUE)
if (length(gamresult_files) != 1) stop("Expected one scaled gamresults file in: ", interfileFolder)
gamresult <- readRDS(gamresult_files[[1]])
if (!is.data.frame(gamresult) || !"parcel" %in% names(gamresult)) stop("Invalid gamresults: ", gamresult_files[[1]])

elementnum <- nrow(gamresult)
ds.resolution <- infer_resolution_from_edges(elementnum)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)
if (is_windows) {
  n_cores <- 1L
}

gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)

## match original boundary handling (HCP-D age range)
gamresult$increase.onset[gamresult$sig == FALSE] <- NA
gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2, 2) == 8.08] <- NA
gamresult$increase.offset[gamresult$sig == FALSE] <- NA
gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2, 2) == 21.92] <- NA
gamresult$peak.change[gamresult$sig == FALSE] <- NA
gamresult$peak.increase.change[gamresult$sig == FALSE] <- NA

gamresult <- gamresult %>% mutate(meanderv2_2 = meanderv2)
gamresult$meanderv2_2[(gamresult$meanderv2 < (mean(gamresult$meanderv2, na.rm = TRUE) - 3 * sd(gamresult$meanderv2, na.rm = TRUE))) |
  (gamresult$meanderv2 > (mean(gamresult$meanderv2, na.rm = TRUE) + 3 * sd(gamresult$meanderv2, na.rm = TRUE)))] <- NA

out_summary <- file.path(resultFolder, "SCrank_correlation_summary.csv")
if (!force && file.exists(out_summary) && !skip_compute) {
  message("[INFO] S4 summary exists; skipping summary recompute: ", out_summary)
  tryCatch({
    SCrank_correlation <- read.csv(out_summary, stringsAsFactors = FALSE)
    print_spearman_summary(SCrank_correlation)
  }, error = function(e) {
    message("[WARN] Failed to read/print S4 summary: ", out_summary, " | ", conditionMessage(e))
  })
} else if (!skip_compute) {
  computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change", "meanderv2_2")
  SCrank_correlation <- do.call(
    rbind,
    lapply(computevar.list, function(computevar) SCrankcorr(gamresult, computevar, ds.resolution, dsdata = FALSE))
  )
  write.csv(SCrank_correlation, out_summary, row.names = FALSE)
  print_spearman_summary(SCrank_correlation)
} else if (!file.exists(out_summary)) {
  message("[WARN] Skip summary computation on Windows; correlation summary not written: ", out_summary)
}

## scatter plots (match original)
out_partial_tiff <- file.path(FigCorrFolder, paste0("meanpartialRsq_SCrankcorr_n", ds.resolution, ".tiff"))
out_partial_pdf <- file.path(FigCorrFolder, paste0("meanpartialRsq_SCrankcorr_n", ds.resolution, ".pdf"))
out_partial_svg <- file.path(FigCorrFolder, paste0("meanpartialRsq_SCrankcorr_n", ds.resolution, ".svg"))
out_meanderv2_tiff <- file.path(FigCorrFolder, paste0("meanmeanderv2_2_SCrankcorr_n", ds.resolution, ".tiff"))
out_meanderv2_pdf <- file.path(FigCorrFolder, paste0("meanmeanderv2_2_SCrankcorr_n", ds.resolution, ".pdf"))
out_meanderv2_svg <- file.path(FigCorrFolder, paste0("meanmeanderv2_2_SCrankcorr_n", ds.resolution, ".svg"))

require_svg <- is_windows && requireNamespace("svglite", quietly = TRUE)
if (!force &&
  file.exists(out_partial_tiff) && file.exists(out_partial_pdf) &&
  file.exists(out_meanderv2_tiff) && file.exists(out_meanderv2_pdf) &&
  (!require_svg || (file.exists(out_partial_svg) && file.exists(out_meanderv2_svg)))) {
  message("[INFO] Yeo S4 scatter outputs exist; skipping. Set --force=1 to re-run.")
} else {
  mytheme <- theme(
    axis.text = element_text(size = 23, color = "black"),
    axis.title = element_text(size = 23),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )

  x_breaks <- if (ds.resolution == 7) {
    seq(0, 20, by = 5)
  } else if (ds.resolution == 17) {
    seq(0, 120, by = 20)
  } else {
    pretty(c(0, elementnum), n = 6)
  }
  x_limits <- if (ds.resolution == 7) {
    c(0, 20)
  } else if (ds.resolution == 17) {
    c(0, 120)
  } else {
    c(0, elementnum)
  }
  scr_min <- if (ds.resolution == 7) 0 else if (ds.resolution == 17) 0 else 0
  scr_max <- if (ds.resolution == 7) 20 else if (ds.resolution == 17) 120 else elementnum

  ## partial Rsq
  correlation.df <- SCrankcorr(gamresult, "partialRsq", ds.resolution, dsdata = TRUE)
  p1 <- ggplot(data = correlation.df) +
    geom_point(aes(x = SCrank, y = partialRsq, color = SCrank), size = 3.5, alpha = 0.9) +
    geom_smooth(aes(x = SCrank, y = partialRsq), method = "lm", color = "black") +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(scr_min, scr_max), guide = "none") +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    labs(x = "S-A connectional axis rank", y = "Partial R2") +
    theme_classic() + mytheme
  ggsave(out_partial_tiff, p1, width = 17, height = 14, units = "cm", bg = "transparent")
  ggsave(out_partial_pdf, p1, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
  if (require_svg) {
    ggsave(out_partial_svg, p1, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
  } else if (is_windows) {
    message("[WARN] svglite not available; skip svg output on Windows.")
  }

  ## meanderv2_2
  correlation.df <- SCrankcorr(gamresult, "meanderv2_2", ds.resolution, dsdata = TRUE)
  p2 <- ggplot(data = correlation.df) +
    geom_point(aes(x = SCrank, y = meanderv2_2, color = SCrank), size = 5.5, alpha = 0.9) +
    geom_smooth(aes(x = SCrank, y = meanderv2_2), linewidth = 2, method = "lm", color = "black") +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(scr_min, scr_max), guide = "none") +
    labs(x = "S-A connectional axis rank", y = "Second derivative") +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(breaks = c(-0.002, 0, 0.002), labels = c(-2, 0, 2)) +
    theme_classic() + mytheme
  ggsave(out_meanderv2_tiff, p2, width = 13, height = 12, units = "cm", bg = "transparent")
  ggsave(out_meanderv2_pdf, p2, dpi = 600, width = 17.5, height = 15, units = "cm", bg = "transparent")
  if (require_svg) {
    ggsave(out_meanderv2_svg, p2, dpi = 600, width = 17.5, height = 15, units = "cm", bg = "transparent")
  } else if (is_windows) {
    message("[WARN] svglite not available; skip svg output on Windows.")
  }
}

## matrix graphs (match original naming)
if (make_matrix_graphs) {
  FigMatrixFolder <- file.path(FigureRoot, paste0("Matrix", ds.resolution, "_sumSCinvnode_gamstats_Age8_22"))
  dir.create(FigMatrixFolder, showWarnings = FALSE, recursive = TRUE)

  computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2_2")
  linerange_frame <- data.frame(
    x = c(0.5, ds.resolution + 0.5), ymin = rep(-ds.resolution - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -ds.resolution - 0.5), xmin = rep(0.5, times = 2), xmax = rep(ds.resolution + 0.5, times = 2)
  )

  if (is_windows) {
    SCrank_correlation.df <- lapply(seq_along(computevar.list), function(i) {
      computevar <- computevar.list[[i]]
      SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
    })
  } else {
    SCrank_correlation.df <- mclapply(seq_along(computevar.list), function(i) {
      computevar <- computevar.list[[i]]
      SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
    }, mc.cores = min(6L, n_cores))
  }

  for (i in seq_along(computevar.list)) {
    computevar <- computevar.list[[i]]
    filename <- file.path(FigMatrixFolder, paste0(computevar, "_", ds.resolution, "net_delLM_CV", CVthr, ".tiff"))
    if (!force && file.exists(filename)) next

    Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol = ds.resolution)
    Matrix.tmp[lower.tri(Matrix.tmp, diag = TRUE)] <- SCrank_correlation.df[[i]][, 2]
    Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
    colnames(Matrix.tmp) <- seq(1, ds.resolution)
    rownames(Matrix.tmp) <- seq(1, ds.resolution)
    matrixtmp.df <- as.data.frame(Matrix.tmp)
    matrixtmp.df$nodeid <- seq(1, ds.resolution)
    matrixtmp.df.melt <- melt(matrixtmp.df, id.vars = c("nodeid"))
    matrixtmp.df.melt$variable <- as.numeric(matrixtmp.df.melt$variable)
    matrixtmp.df.melt$nodeid <- 0 - matrixtmp.df.melt$nodeid
    matrixtmp.df.melt$value <- as.numeric(matrixtmp.df.melt$value)

    if (computevar == "partialRsq") {
      lmthr <- max(abs(gamresult$partialRsq), na.rm = TRUE)
      Matrix.tmp.sig <- matrix(NA, nrow = ds.resolution, ncol = ds.resolution)
      Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = TRUE)] <- (gamresult$pfdr < 0.05)
      Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
      colnames(Matrix.tmp.sig) <- seq(1, ds.resolution)
      rownames(Matrix.tmp.sig) <- seq(1, ds.resolution)
      matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
      matrixtmp.df.sig$nodeid <- seq(1, ds.resolution)
      matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig, id.vars = c("nodeid"))
      matrixtmp.df.sig.melt$variable <- as.numeric(matrixtmp.df.sig.melt$variable)
      matrixtmp.df.sig.melt$nodeid <- 0 - matrixtmp.df.sig.melt$nodeid
      matrixtmp.df.sig.melt$value <- as.numeric(matrixtmp.df.sig.melt$value)
      matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value == 0), ]

      Fig <- ggplot(data = matrixtmp.df.melt) +
        geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
        scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "grey") +
        scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "grey") +
        geom_text(data = matrixtmp.df.sig.melt, aes(x = variable, y = nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size = 10) +
        geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
        geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
        geom_segment(aes(x = 0.5, y = -0.5, xend = ds.resolution + 0.5, yend = -ds.resolution - 0.5), color = "black", linewidth = 0.5) +
        labs(x = NULL, y = NULL) +
        scale_y_continuous(breaks = NULL, labels = NULL) +
        scale_x_continuous(breaks = NULL, labels = NULL) +
        theme(
          axis.line = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12, angle = 315, hjust = 1, vjust = 1),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 1),
          plot.background = element_rect(fill = "transparent", color = NA)
        )
    } else {
      lmthr <- max(abs(gamresult[[computevar]]), na.rm = TRUE)
      Fig <- ggplot(data = matrixtmp.df.melt) +
        geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
        scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "#053061") +
        scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "#053061") +
        geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
        geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
        geom_segment(aes(x = 0.5, y = -0.5, xend = ds.resolution + 0.5, yend = -ds.resolution - 0.5), color = "black", linewidth = 0.5) +
        labs(x = NULL, y = NULL) +
        scale_y_continuous(breaks = NULL, labels = NULL) +
        scale_x_continuous(breaks = NULL, labels = NULL) +
        theme(
          axis.line = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12, angle = 315, hjust = 1, vjust = 1),
          axis.title = element_text(size = 18),
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
