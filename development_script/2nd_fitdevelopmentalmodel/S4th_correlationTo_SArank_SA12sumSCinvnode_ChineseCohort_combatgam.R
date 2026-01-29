## Rscript version of:
##   S4th_correlationTo_SArank_SA12sumSCinvnode_ChineseCohort.Rmd
##
## This script performs Spearman correlations between GAM-derived developmental
## statistics and the S-A connectional axis rank, and generates scatter plots.
##
## Default inputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds
## - wd/interdataFolder_HCPD/average_EuclideanDistance_12.csv
##
## Outputs:
## - outputs/results/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/SCrank_correlation_summary.csv
## - outputs/figures/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/correlation_sumSCinvnode_SCrank/*.tiff (+ one *.pdf)
## - (optional) outputs/figures/.../Matrix12_sumSCinvnode_gamstats_Age8_22/*.tiff

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
is_windows <- .Platform$OS.type == "windows"
skip_compute_on_windows <- as.integer(if (!is.null(args$skip_compute_on_windows)) args$skip_compute_on_windows else 1L) == 1L
skip_compute <- is_windows && skip_compute_on_windows

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2
make_matrix_graphs <- as.integer(if (!is.null(args$make_matrix_graphs)) args$make_matrix_graphs else 0L) == 1L

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "chinese", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "chinese", "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "chinese", "combat_gam", paste0("CV", CVthr))

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
if (!"anova.smooth.pvalue" %in% names(gamresult)) {
  stop("Missing anova.smooth.pvalue in gamresults; check S1 outputs.")
}
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "chinese",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)
SCdata <- readRDS(input_rds)

euclid_csv <- if (!is.null(args$euclid_csv)) {
  args$euclid_csv
} else {
  file.path(project_root, "wd", "interdataFolder_HCPD", "average_EuclideanDistance_12.csv")
}
if (!file.exists(euclid_csv)) stop("Missing euclid_csv: ", euclid_csv)
EucDistance <- read.csv(euclid_csv)

message(sum(gamresult$sig), " edges have significant developmental effects.")

out_summary <- file.path(resultFolder, "SCrank_correlation_summary.csv")
FigCorrFolder <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank")
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

out_fig_meanderv2 <- file.path(FigCorrFolder, paste0("meanmeanderv2_c_SCrankcorr_n", ds.resolution, ".tiff"))
out_fig_meanderv2_pdf <- file.path(FigCorrFolder, paste0("meanmeanderv2_c_SCrankcorr_n", ds.resolution, ".pdf"))
out_fig_meanderv2_ctrl_tiff <- file.path(FigCorrFolder, paste0("meanmeanderv2_c_control_distance_SCrankcorr_n", ds.resolution, ".tiff"))
out_fig_meanderv2_ctrl_pdf <- file.path(FigCorrFolder, paste0("meanmeanderv2_c_control_distance_SCrankcorr_n", ds.resolution, ".pdf"))

if (!force && !skip_compute &&
  file.exists(out_summary) &&
  file.exists(out_fig_meanderv2) &&
  file.exists(out_fig_meanderv2_pdf) &&
  file.exists(out_fig_meanderv2_ctrl_tiff) &&
  file.exists(out_fig_meanderv2_ctrl_pdf)) {
  message("[INFO] S4 outputs exist; skipping S4. Set --force=1 to re-run.")
  tryCatch({
    SCrank_correlation <- read.csv(out_summary, stringsAsFactors = FALSE)
    print_spearman_summary(SCrank_correlation)
  }, error = function(e) {
    message("[WARN] Failed to read/print S4 summary: ", out_summary, " | ", conditionMessage(e))
  })
  quit(save = "no", status = 0)
}

placeholder_plot <- function(title, subtitle) {
  ggplot() +
    annotate("text", x = 0, y = 0.2, label = title, size = 6) +
    annotate("text", x = 0, y = -0.2, label = subtitle, size = 4.5) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", color = NA))
}

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (length(sc_cols) == 0) stop("No SC.* columns found in input_rds: ", input_rds)
meanSC <- colMeans(SCdata[, sc_cols, drop = FALSE], na.rm = TRUE)
meanSC_map <- setNames(meanSC, sc_cols)

parcel_all <- paste0("SC.", seq_len(elementnum), "_h")

if (!"Edistance" %in% names(EucDistance)) stop("Euclidean distance CSV missing column 'Edistance': ", euclid_csv)
if ("SC_label" %in% names(EucDistance)) {
  Edist_map <- setNames(EucDistance$Edistance, as.character(EucDistance$SC_label))
} else if ("parcel" %in% names(EucDistance)) {
  Edist_map <- setNames(EucDistance$Edistance, as.character(EucDistance$parcel))
} else if (nrow(EucDistance) == elementnum) {
  Edist_map <- setNames(EucDistance$Edistance, parcel_all)
} else {
  stop("Euclidean distance CSV cannot be aligned (need SC_label/parcel or nrow==78): ", euclid_csv)
}

meanSC_aligned <- meanSC_map[gamresult$parcel]
gamresult$EucDistance <- unname(Edist_map[gamresult$parcel])

## convert critical ages of insignificantly developmental edges to NA
## convert critical ages equal to age boundaries to NA (age boundaries depend on each edge/model).
for (nm in c("increase.onset", "increase.offset", "decrease.onset", "decrease.offset", "change.onset", "change.offset", "peak.change", "peak.increase.change")) {
  if (nm %in% names(gamresult)) gamresult[[nm]][gamresult$sig == FALSE] <- NA
}

## Robust meanderv2 outlier handling (match historical Rmd behavior)
if ("meanderv2" %in% names(gamresult)) {
  gamresult <- gamresult %>% mutate(
    meanderv2_c = case_when(
      meanderv2 > mean(meanderv2, na.rm = TRUE) + 3 * sd(meanderv2, na.rm = TRUE) ~ NA,
      meanderv2 < mean(meanderv2, na.rm = TRUE) - 3 * sd(meanderv2, na.rm = TRUE) ~ NA,
      .default = meanderv2
    )
  )
} else {
  gamresult$meanderv2_c <- NA_real_
}

## Validation: Control for Euclidean distance
if ("partialRsq" %in% names(gamresult)) {
  fit_pr <- lm(partialRsq ~ EucDistance, data = gamresult, na.action = na.exclude)
  gamresult$partialRsq_control_distance <- as.numeric(residuals(fit_pr))
}

if ("meanderv2_c" %in% names(gamresult)) {
  fit_md <- lm(meanderv2_c ~ EucDistance, data = gamresult, na.action = na.exclude)
  gamresult$meanderv2_c_control_distance <- as.numeric(residuals(fit_md))
}

if (!skip_compute) {
  ## correlations to SC rank (summary table)
  computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2_c")
  SCrank_correlation <- do.call(
    rbind,
    lapply(computevar.list, function(computevar) SCrankcorr(gamresult, computevar, ds.resolution, dsdata = FALSE))
  )

  if ("partialRsq_control_distance" %in% names(gamresult)) {
    SCrank_correlation <- rbind(
      SCrank_correlation,
      SCrankcorr(gamresult, "partialRsq_control_distance", ds.resolution, dsdata = FALSE)
    )
  }

  if ("meanderv2_c_control_distance" %in% names(gamresult)) {
    SCrank_correlation <- rbind(
      SCrank_correlation,
      SCrankcorr(gamresult, "meanderv2_c_control_distance", ds.resolution, dsdata = FALSE)
    )
  }

  write.csv(SCrank_correlation, out_summary, row.names = FALSE)
  print_spearman_summary(SCrank_correlation)
} else {
  if (!file.exists(out_summary)) {
    message("[WARN] Skip summary computation on Windows; correlation summary not written: ", out_summary)
  } else {
    message("[INFO] Skip summary computation on Windows; using existing summary: ", out_summary)
  }
}

## scatter plots (match historical Chinese Rmd output count/names)
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

## 1) meanderv2_c: tiff + pdf (+ svg on Windows)
if ("meanderv2_c" %in% names(gamresult)) {
  computevar <- "meanderv2_c"
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  if (!any(is.finite(correlation.df$meanderv2_c))) {
    p_meanderv2 <- placeholder_plot("meanderv2_c", "No finite values available.")
  } else {
    p_meanderv2 <- ggplot(data = correlation.df) +
      geom_point(aes(x = SCrank, y = meanderv2_c, color = SCrank), size = 5) +
      geom_smooth(aes(x = SCrank, y = meanderv2_c), linewidth = 2, method = "lm", color = "black") +
      scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, guide = "none") +
      labs(x = "S-A connectional axis rank", y = "Second derivative") +
      scale_y_continuous(breaks = c(-0.002, -0.001, 0, 0.001), labels = c(-2, -1, 0, 1)) +
      theme_classic() + mytheme
  }
  ggsave(out_fig_meanderv2, p_meanderv2, width = 13, height = 12, units = "cm", bg = "transparent")
  ggsave(out_fig_meanderv2_pdf, p_meanderv2, dpi = 600, width = 13, height = 12, units = "cm", bg = "transparent")
  if (is_windows) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      message("[WARN] svglite not available; skip svg output on Windows.")
    } else {
      out_fig_meanderv2_svg <- file.path(FigCorrFolder, paste0("meanmeanderv2_c_SCrankcorr_n", ds.resolution, ".svg"))
      ggsave(out_fig_meanderv2_svg, p_meanderv2, dpi = 600, width = 13, height = 12, units = "cm", bg = "transparent")
    }
  }
}

## 2) meanderv2_c_control_distance: tiff + pdf (+ svg on Windows)
if ("meanderv2_c_control_distance" %in% names(gamresult)) {
  computevar <- "meanderv2_c_control_distance"
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  if (!any(is.finite(correlation.df$meanderv2_c_control_distance))) {
    p_meanderv2_ctrl <- placeholder_plot("meanderv2_c_control_distance", "No finite values available.")
  } else {
    p_meanderv2_ctrl <- ggplot(data = correlation.df) +
      geom_point(aes(x = SCrank, y = meanderv2_c_control_distance, color = SCrank), size = 5.5) +
      geom_smooth(aes(x = SCrank, y = meanderv2_c_control_distance), linewidth = 1.2, method = "lm", color = "black") +
      scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, guide = "none") +
      labs(x = "S-A connectional axis rank", y = "Second derivative") +
      scale_y_continuous(breaks = c(-0.002, -0.001, 0, 0.001), labels = c(-2, -1, 0, 1)) +
      theme_classic() + mytheme
  }
  ggsave(out_fig_meanderv2_ctrl_tiff, p_meanderv2_ctrl, width = 16, height = 14, units = "cm", bg = "transparent")
  ggsave(out_fig_meanderv2_ctrl_pdf, p_meanderv2_ctrl, dpi = 600, width = 16, height = 14, units = "cm", bg = "transparent")
  if (is_windows) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      message("[WARN] svglite not available; skip svg output on Windows.")
    } else {
      out_fig_meanderv2_ctrl_svg <- file.path(FigCorrFolder, paste0("meanmeanderv2_c_control_distance_SCrankcorr_n", ds.resolution, ".svg"))
      ggsave(out_fig_meanderv2_ctrl_svg, p_meanderv2_ctrl, dpi = 600, width = 16, height = 14, units = "cm", bg = "transparent")
    }
  }
}

## matrix graphs (optional)
if (make_matrix_graphs) {
  Matrix.tmp <- matrix(NA, nrow = 12, ncol = 12)
  computevar.list2 <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2_c")

  FigMatrixFolder <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats_Age8_22")
  dir.create(FigMatrixFolder, showWarnings = FALSE, recursive = TRUE)

  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5), ymin = rep(-12 - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5), xmin = rep(0.5, times = 2), xmax = rep(12 + 0.5, times = 2)
  )

  SCrank_correlation.df <- mclapply(seq_along(computevar.list2), function(i) {
    computevar <- computevar.list2[[i]]
    SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  }, mc.cores = min(6L, n_cores))

  colorbar.prob <- c(
    0.5, 0.4, 0.6, 0.5,
    abs(min(gamresult$meanderv2, na.rm = TRUE)) / (max(gamresult$meanderv2, na.rm = TRUE) - min(gamresult$meanderv2, na.rm = TRUE))
  )

  for (i in seq_along(computevar.list2)) {
    computevar <- computevar.list2[[i]]
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

    values_vec <- SCrank_correlation.df[[i]][, 2]
    if (!any(is.finite(values_vec))) {
      Fig <- placeholder_plot(computevar, "No finite values available.")
      filename <- file.path(FigMatrixFolder, paste0(computevar, "_12net_delLM_CV", CVthr, ".tiff"))
      ggsave(filename, Fig, height = 18, width = 20, units = "cm", dpi = 600, bg = "transparent")
      next
    }

    colorbarvalues.tmp <- colorbarvalues(values_vec, colorbar.prob[[i]])
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
        scale_fill_distiller(type = "seq", palette = "RdBu", limit = c(-lmthr, lmthr), na.value = "grey") +
        scale_color_distiller(type = "seq", palette = "RdBu", limit = c(-lmthr, lmthr), na.value = "grey") +
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
    } else if (computevar == "meanderv2_c") {
      lmthr <- max(abs(gamresult$meanderv2_c), na.rm = TRUE)
      Fig <- ggplot(data = matrixtmp.df.melt) +
        geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
        scale_fill_distiller(type = "seq", palette = "RdBu", limit = c(-lmthr, lmthr), na.value = "#053061") +
        scale_color_distiller(type = "seq", palette = "RdBu", limit = c(-lmthr, lmthr), na.value = "#053061") +
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

    filename <- file.path(FigMatrixFolder, paste0(computevar, "_12net_delLM_CV", CVthr, ".tiff"))
    ggsave(filename, Fig, height = 18, width = 20, units = "cm", dpi = 600, bg = "transparent")
  }
}
