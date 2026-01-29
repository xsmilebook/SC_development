## Rscript version of:
##   S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD.Rmd
##
## This script performs Spearman correlations between GAM-derived developmental
## statistics and the S-A connectional axis rank, and generates scatter plots.
##
## Default inputs (project-relative):
## - outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds
## - outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds
## - wd/interdataFolder_HCPD/average_EuclideanDistance_12.csv
##
## Outputs:
## - outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/SCrank_correlation_summary.csv
## - outputs/figures/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/correlation_sumSCinvnode_SCrank/*.tiff + *.pdf
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
ds.resolution <- as.integer(if (!is.null(args$ds_res)) args$ds_res else 12L)
elementnum <- ds.resolution * (ds.resolution + 1) / 2
make_matrix_graphs <- as.integer(if (!is.null(args$make_matrix_graphs)) args$make_matrix_graphs else 0L) == 1L
if (is_windows && is_interactive && is.null(args$make_matrix_graphs)) {
  make_matrix_graphs <- FALSE
}
sa_axis_mode <- tolower(if (!is.null(args$sa_axis_mode)) args$sa_axis_mode else "addsquare")
out_tag <- if (!is.null(args$out_tag)) args$out_tag else ""
tag_suffix <- if (nzchar(out_tag)) paste0("_", out_tag) else ""

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
resultFolder <- file.path(
  project_root, "outputs", "results", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "hcpd", "combat_gam", paste0("CV", CVthr))

dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "colorbarvalue.R"))

infer_resolution_from_edges <- function(n_edges) {
  n <- (sqrt(8 * n_edges + 1) - 1) / 2
  n_int <- as.integer(round(n))
  if (n_int < 1 || n_int * (n_int + 1) / 2 != n_edges) {
    stop("Cannot infer resolution from n_edges=", n_edges, " (expected triangular number).")
  }
  n_int
}

build_sa_rank <- function(ds, mode = c("addsquare", "multiply")) {
  mode <- match.arg(mode)
  mat <- matrix(NA_real_, nrow = ds, ncol = ds)
  for (x in 1:ds) {
    for (y in 1:ds) {
      mat[x, y] <- if (mode == "multiply") x * y else x^2 + y^2
    }
  }
  mat <- mat[lower.tri(mat, diag = TRUE)]
  rank(mat, ties.method = "average")
}

sa_axis_mode <- if (sa_axis_mode %in% c("addsquare", "x2+y2", "sum_sq")) "addsquare" else if (sa_axis_mode %in% c("multiply", "x*y", "mult")) "multiply" else sa_axis_mode
if (!sa_axis_mode %in% c("addsquare", "multiply")) stop("Invalid --sa_axis_mode: ", sa_axis_mode, " (use addsquare|multiply)")

SCrankcorr_custom <- function(gamresult, computevar, sa_rank_map, ds.resolution, dsdata = FALSE) {
  if (!"parcel" %in% names(gamresult)) stop("gamresult missing 'parcel' column.")
  if (!computevar %in% names(gamresult)) stop("gamresult missing computevar: ", computevar)
  scr <- unname(sa_rank_map[as.character(gamresult$parcel)])
  df <- data.frame(SCrank = scr, computevar = as.numeric(gamresult[[computevar]]))

  ct <- suppressWarnings(corr.test(df, method = "spearman"))
  correstimate <- ct$r[2, 1]
  p.spearman <- ct$p[2, 1]
  SCrankdf <- data.frame(
    ds.resolution = ds.resolution,
    sa_axis_mode = sa_axis_mode,
    Interest.var = computevar,
    r.spearman = correstimate,
    p.spearman = p.spearman
  )
  if (dsdata) {
    names(df) <- c("SCrank", computevar)
    return(df)
  }
  SCrankdf
}

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)
if (is_windows) {
  n_cores <- 1L
}

## load data
gamresult <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "hcpd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)
SCdata <- readRDS(input_rds)

euclid_csv <- if (!is.null(args$euclid_csv)) {
  args$euclid_csv
} else {
  file.path(project_root, "wd", "interdataFolder_HCPD", paste0("average_EuclideanDistance_", ds.resolution, ".csv"))
}
do_euclid <- is.character(euclid_csv) && nzchar(euclid_csv) && file.exists(euclid_csv)
EuricDistance <- if (do_euclid) read.csv(euclid_csv) else NULL

message(sum(gamresult$sig), " edges have significant developmental effects.")

out_summary <- file.path(resultFolder, paste0("SCrank_correlation_summary", tag_suffix, ".csv"))
should_recompute_summary <- !skip_compute && (force || !file.exists(out_summary))
sum_df_existing <- NULL
if (!should_recompute_summary) {
  message("[INFO] S4 summary exists; will reuse unless --force=1: ", out_summary)
  sum_df_existing <- tryCatch(read.csv(out_summary), error = function(e) NULL)
  if (is.data.frame(sum_df_existing) && all(c("Interest.var", "r.spearman", "p.spearman") %in% names(sum_df_existing))) {
    message(sprintf("[INFO] SCrank Spearman summary (ds.resolution=%d sa_axis_mode=%s tag=%s)", ds.resolution, sa_axis_mode, ifelse(nzchar(out_tag), out_tag, "<none>")))
    for (i in seq_len(nrow(sum_df_existing))) {
      message(sprintf("[INFO]  %s: r=%.5f, p=%.3g", as.character(sum_df_existing$Interest.var[[i]]), as.numeric(sum_df_existing$r.spearman[[i]]), as.numeric(sum_df_existing$p.spearman[[i]])))
    }
  }
}

## description
sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (length(sc_cols) == 0) stop("No SC.* columns found in input_rds: ", input_rds)
meanSC <- colMeans(SCdata[, sc_cols, drop = FALSE])

parcel_all <- paste0("SC.", seq_len(elementnum), "_h")
meanSC_map <- setNames(meanSC, sc_cols)

Edist_map <- NULL
if (do_euclid) {
  if (!"Edistance" %in% names(EuricDistance)) stop("Euclidean distance CSV missing column 'Edistance': ", euclid_csv)
  if ("SC_label" %in% names(EuricDistance)) {
    Edist_map <- setNames(as.numeric(EuricDistance$Edistance), as.character(EuricDistance$SC_label))
  } else if ("parcel" %in% names(EuricDistance)) {
    Edist_map <- setNames(as.numeric(EuricDistance$Edistance), as.character(EuricDistance$parcel))
  } else if (nrow(EuricDistance) == elementnum) {
    Edist_map <- setNames(as.numeric(EuricDistance$Edistance), parcel_all)
  } else {
    stop("Euclidean distance CSV cannot be aligned (need SC_label/parcel or rows==elementnum): ", euclid_csv)
  }

  # Sanity check: require full coverage for the current resolution; otherwise skip to avoid misleading results.
  ed_mapped <- unname(Edist_map[as.character(gamresult$parcel)])
  n_mapped <- sum(!is.na(ed_mapped))
  if (n_mapped != elementnum) {
    message(sprintf(
      "[WARN] Euclidean distance map coverage %d/%d (ds.resolution=%d); skip control-distance correlations. (euclid_csv=%s)",
      n_mapped, elementnum, ds.resolution, euclid_csv
    ))
    do_euclid <- FALSE
    EuricDistance <- NULL
    Edist_map <- NULL
  }
}

meanSC_aligned <- meanSC_map[gamresult$parcel]
if (!skip_compute) {
  if (sum(!is.na(meanSC_aligned)) == elementnum) {
    corr.test(meanSC_aligned, gamresult$partialRsq)
    if (do_euclid) {
      corr.test(meanSC_aligned, Edist_map[gamresult$parcel])
    }
  } else {
    message(sprintf(
      "[WARN] meanSC columns (%d) do not cover ds.resolution=%d (edges=%d); skip meanSC correlation checks.",
      length(sc_cols), ds.resolution, elementnum
    ))
  }
}

## convert critical ages of insignificantly developmental edges to NA
if (do_euclid) {
  gamresult$EuricDistance <- unname(Edist_map[gamresult$parcel])
} else {
  gamresult$EuricDistance <- NA_real_
}
gamresult$increase.onset[gamresult$sig == FALSE] <- NA
gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2, 2) == 8.08] <- NA
gamresult$increase.offset[gamresult$sig == FALSE] <- NA
gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2, 2) == 21.92] <- NA
gamresult$peak.change[gamresult$sig == FALSE] <- NA
gamresult$peak.increase.change[gamresult$sig == FALSE] <- NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")

## compute correlations to SC rank
computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change", "meanderv2")
sa_rank_all <- build_sa_rank(ds.resolution, mode = sa_axis_mode)
sa_rank_map <- setNames(sa_rank_all, parcel_all)
SCrank_correlation <- NULL
if (should_recompute_summary) {
  SCrank_correlation <- do.call(
    rbind,
    lapply(computevar.list, function(computevar) SCrankcorr_custom(gamresult, computevar, sa_rank_map, ds.resolution, dsdata = FALSE))
  )
}

## Validation: Control for Euclidean distance
if (do_euclid) {
  gamresult$meanderv2_control_distance <- residuals(lm(meanderv2 ~ EuricDistance, data = gamresult, na.action = na.exclude))
  gamresult$partialRsq_control_distance <- residuals(lm(partialRsq ~ EuricDistance, data = gamresult, na.action = na.exclude))
  if (should_recompute_summary) {
    SCrank_correlation <- rbind(
      SCrank_correlation,
      SCrankcorr_custom(gamresult, "meanderv2_control_distance", sa_rank_map, ds.resolution, dsdata = FALSE),
      SCrankcorr_custom(gamresult, "partialRsq_control_distance", sa_rank_map, ds.resolution, dsdata = FALSE)
    )
  }
} else {
  message("[INFO] euclid_csv not provided or missing; skip control-distance correlations.")
}
if (should_recompute_summary) {
  write.csv(SCrank_correlation, out_summary, row.names = FALSE)
} else if (skip_compute && !file.exists(out_summary)) {
  message("[WARN] Skip summary computation on Windows; correlation summary not written: ", out_summary)
}

## scatter plots
FigCorrFolder <- file.path(FigureRoot, paste0("correlation_sumSCinvnode_SCrank", tag_suffix))
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

plots_complete <- function() {
  scatter_targets <- c(
    "partialRsq",
    "meanderv2",
    "meanderv2_control_distance",
    "partialRsq_control_distance",
    "increase.onset2",
    "increase.offset2",
    "peak.increase.change"
  )
  present <- scatter_targets[scatter_targets %in% names(gamresult)]
  if (length(present) == 0) return(FALSE)
  required <- c()
  for (nm in present) {
    required <- c(
      required,
      file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".tiff")),
      file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".pdf"))
    )
    if (is_windows && requireNamespace("svglite", quietly = TRUE)) {
      required <- c(
        required,
        file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".svg"))
      )
    }
  }
  all(file.exists(required))
}

require_svg <- is_windows && requireNamespace("svglite", quietly = TRUE)
if (!force && !should_recompute_summary && plots_complete()) {
  message("[INFO] S4 plots already exist; skipping plotting (set --force=1 to re-draw).")
  if (!is_interactive) quit(save = "no", status = 0)
}

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

plot_one_scatter <- function(computevar, ylab) {
  df <- SCrankcorr_custom(gamresult, computevar, sa_rank_map, ds.resolution, dsdata = TRUE)
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

  point_size <- if (computevar == "partialRsq") 3.5 else 5.5
  x_breaks <- if (ds.resolution == 7) {
    seq(0, 20, by = 5)
  } else if (ds.resolution == 17) {
    seq(0, 120, by = 20)
  } else {
    seq(0, 80, by = 20)
  }
  x_limits <- if (ds.resolution == 7) {
    c(0, 20)
  } else if (ds.resolution == 17) {
    c(0, 120)
  } else {
    c(0, 80)
  }

  scr_min <- if (ds.resolution == 7) 0 else if (ds.resolution == 17) 0 else min(df$SCrank, na.rm = TRUE)
  scr_max <- if (ds.resolution == 7) 20 else if (ds.resolution == 17) 120 else max(df$SCrank, na.rm = TRUE)

  p <- ggplot(df) +
    geom_point(aes(x = SCrank, y = .data[[computevar]], color = SCrank), size = point_size, alpha = 0.9) +
    geom_smooth(aes(x = SCrank, y = .data[[computevar]]), linewidth = 2, method = "lm", color = "black") +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(scr_min, scr_max), guide = "none") +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    labs(
      x = "S-A connectional axis rank",
      y = ylab
    ) +
    theme_classic() + mytheme

  if (computevar == "meanderv2" || computevar == "meanderv2_control_distance") {
    p <- p + scale_y_continuous(breaks = c(-0.003, 0, 0.003), labels = c(-3, 0, 3))
  }

  p
}

scatter_targets <- list(
  partialRsq = "Partial R^2 (signed)",
  meanderv2 = "Second derivative",
  meanderv2_control_distance = "Second derivative",
  partialRsq_control_distance = "Partial R^2 (residualized by Euclidean distance)",
  increase.onset2 = "Increase onset age",
  increase.offset2 = "Increase offset age",
  peak.increase.change = "Peak increasing change age"
)

for (nm in names(scatter_targets)) {
  if (!nm %in% names(gamresult)) next
  p <- plot_one_scatter(nm, scatter_targets[[nm]])
  if (nm == "partialRsq") {
    ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".tiff")), p, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
    ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".pdf")), p, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
    if (require_svg) {
      ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".svg")), p, dpi = 600, width = 17, height = 14, units = "cm", bg = "transparent")
    } else if (is_windows) {
      message("[WARN] svglite not available; skip svg output on Windows.")
    }
  } else {
    ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".tiff")), p, dpi = 600, width = 13, height = 12, units = "cm", bg = "transparent")
    ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".pdf")), p, dpi = 600, width = 17.5, height = 15, units = "cm", bg = "transparent")
    if (require_svg) {
      ggsave(file.path(FigCorrFolder, paste0("mean", nm, "_SCrankcorr_n", ds.resolution, ".svg")), p, dpi = 600, width = 17.5, height = 15, units = "cm", bg = "transparent")
    } else if (is_windows) {
      message("[WARN] svglite not available; skip svg output on Windows.")
    }
  }
}

## matrix graphs for resolution of 12 (optional; heavy)
if (make_matrix_graphs) {
  Matrix.tmp <- matrix(NA, nrow = 12, ncol = 12)
  computevar.list2 <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change", "meanderv2", "meanderv2_control_distance")

  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5), ymin = rep(-12 - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5), xmin = rep(0.5, times = 2), xmax = rep(12 + 0.5, times = 2)
  )

  if (is_windows) {
    SCrank_correlation.df <- lapply(seq_along(computevar.list2), function(i) {
      computevar <- computevar.list2[[i]]
      SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
    })
  } else {
    SCrank_correlation.df <- mclapply(seq_along(computevar.list2), function(i) {
      computevar <- computevar.list2[[i]]
      SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)
    }, mc.cores = min(6L, n_cores))
  }

  colorbar.prob <- c(
    0.5, 0.4, 0.6, 0.5,
    abs(min(gamresult$meanderv2, na.rm = TRUE)) / (max(gamresult$meanderv2, na.rm = TRUE) - min(gamresult$meanderv2, na.rm = TRUE)),
    0.5
  )

  FigMatrixFolder <- file.path(FigureRoot, "Matrix12_sumSCinvnode_gamstats_Age8_22")
  dir.create(FigMatrixFolder, showWarnings = FALSE, recursive = TRUE)

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
      Fig <- ggplot(data = matrixtmp.df.melt) +
        geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
        scale_fill_distiller(type = "seq", palette = "RdBu", values = colorbarvalues.tmp, na.value = "grey") +
        scale_color_distiller(type = "seq", palette = "RdBu", values = colorbarvalues.tmp, na.value = "grey") +
        geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
        geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
        geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
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

    ggsave(file.path(FigMatrixFolder, paste0(computevar, "_12net_delLM_CV", CVthr, ".tiff")), Fig, height = 18, width = 20, units = "cm", bg = "transparent")
  }
}
