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
## - outputs/figures/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/correlation_sumSCinvnode_SCrank/*.tiff + *.pdf
## - (optional) outputs/figures/.../Matrix12_sumSCinvnode_gamstats_*.tiff

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
if (!force && file.exists(out_summary)) {
  message("[INFO] S4 summary exists, skipping: ", out_summary, " (set --force=1 to re-run)")
  quit(save = "no", status = 0)
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

## correlations to SC rank (summary table)
computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2_c")
SCrank_correlation <- do.call(
  rbind,
  lapply(computevar.list, function(computevar) SCrankcorr(gamresult, computevar, ds.resolution, dsdata = FALSE))
)

## Validation: Control for Euclidean distance
if ("partialRsq" %in% names(gamresult)) {
  fit_pr <- lm(partialRsq ~ EucDistance, data = gamresult, na.action = na.exclude)
  gamresult$partialRsq_control_distance <- as.numeric(residuals(fit_pr))
  SCrank_correlation <- rbind(
    SCrank_correlation,
    SCrankcorr(gamresult, "partialRsq_control_distance", ds.resolution, dsdata = FALSE)
  )
}

if ("meanderv2_c" %in% names(gamresult)) {
  fit_md <- lm(meanderv2_c ~ EucDistance, data = gamresult, na.action = na.exclude)
  gamresult$meanderv2_c_control_distance <- as.numeric(residuals(fit_md))
  SCrank_correlation <- rbind(
    SCrank_correlation,
    SCrankcorr(gamresult, "meanderv2_c_control_distance", ds.resolution, dsdata = FALSE)
  )
}

write.csv(SCrank_correlation, out_summary, row.names = FALSE)

## scatter plots
FigCorrFolder <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank")
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

plot_one_scatter <- function(computevar, ylab) {
  df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  names(df) <- c("SCrank", computevar)
  ct <- suppressWarnings(corr.test(df$SCrank, df[[computevar]], method = "spearman"))
  extract_corr <- function(ct_obj) {
    r_obj <- ct_obj$r
    p_obj <- ct_obj$p
    if (is.matrix(r_obj) || is.data.frame(r_obj)) {
      r_val <- if (ncol(r_obj) >= 2) r_obj[1, 2] else r_obj[1, 1]
      p_val <- if (is.matrix(p_obj) || is.data.frame(p_obj)) {
        if (ncol(p_obj) >= 2) p_obj[1, 2] else p_obj[1, 1]
      } else {
        as.numeric(p_obj)[1]
      }
      return(list(r = as.numeric(r_val), p = as.numeric(p_val)))
    }
    list(r = as.numeric(r_obj)[1], p = as.numeric(p_obj)[1])
  }
  rp <- extract_corr(ct)
  r <- rp$r
  p <- rp$p
  ggplot(df, aes(x = SCrank, y = .data[[computevar]])) +
    geom_point(alpha = 0.9, size = 3.2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.0, color = "black") +
    labs(
      x = "S-A rank",
      y = ylab,
      title = sprintf("%s (Spearman r=%.3f, p=%.3g)", computevar, r, p)
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
}

scatter_targets <- list(
  partialRsq = "Age effect (partial R^2)",
  meanderv2_c = "Mean 2nd derivative",
  partialRsq_control_distance = "Partial R^2 (residualized by Euclidean distance)",
  meanderv2_c_control_distance = "Mean 2nd derivative (residualized by Euclidean distance)"
)

for (nm in names(scatter_targets)) {
  if (!nm %in% names(gamresult)) next
  p <- plot_one_scatter(nm, scatter_targets[[nm]])
  ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".tiff")), p, dpi = 600, width = 20, height = 12, units = "cm", bg = "transparent")
  ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".pdf")), p, dpi = 600, width = 16, height = 14, units = "cm", bg = "transparent")
}

## matrix graphs (optional)
if (make_matrix_graphs) {
  Matrix.tmp <- matrix(NA, nrow = 12, ncol = 12)
  computevar.list2 <- c("partialRsq", "meanderv2_c", "partialRsq_control_distance", "meanderv2_c_control_distance")

  FigMatrixFolder <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats")
  dir.create(FigMatrixFolder, showWarnings = FALSE, recursive = TRUE)

  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5), ymin = rep(-12 - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5), xmin = rep(0.5, times = 2), xmax = rep(12 + 0.5, times = 2)
  )

  SCrank_correlation.df <- mclapply(seq_along(computevar.list2), function(i) {
    computevar <- computevar.list2[[i]]
    SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
  }, mc.cores = min(6L, n_cores))

  for (i in seq_along(computevar.list2)) {
    computevar <- computevar.list2[[i]]
    Matrix.tmp[,] <- NA
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

    lmthr <- max(abs(matrixtmp.df.melt$value), na.rm = TRUE)
    p <- ggplot(data = matrixtmp.df.melt) +
      geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
      scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "grey") +
      scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-lmthr, lmthr), na.value = "grey") +
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

    ggsave(file.path(FigMatrixFolder, paste0("Matrix12_sumSCinvnode_", computevar, ".tiff")), p, width = 12, height = 10, units = "cm", dpi = 600, bg = "transparent")
    ggsave(file.path(FigMatrixFolder, paste0("Matrix12_sumSCinvnode_", computevar, ".pdf")), p, width = 12, height = 10, units = "cm", dpi = 600, bg = "transparent")
  }
}
