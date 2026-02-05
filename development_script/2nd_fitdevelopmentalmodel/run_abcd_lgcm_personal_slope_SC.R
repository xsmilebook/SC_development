#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
})

rm(list = ls())

CVthr <- 75

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

resultFolder <- file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", "abcd", "lgcm_personal_slope")
FigureFolder <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "abcd", "lgcm_personal_slope")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

out_rds <- file.path(resultFolder, paste0("lgcm_personal_slope_SC_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)
force <- identical(Sys.getenv("FORCE", unset = "0"), "1")

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (age/sex/mean_fd variant; longitudinal)")
}

plotdatasum_rds <- Sys.getenv(
  "ABCD_PLOTDATASUM_RDS",
  unset = file.path(
    project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
    "abcd", "combat_gam", "CV75", "plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds"
  )
)
if (!file.exists(plotdatasum_rds)) stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)
plotdata <- readRDS(plotdatasum_rds)

sa12_csv <- Sys.getenv(
  "ABCD_SA12_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "SA12_10.csv")
)
if (!file.exists(sa12_csv)) stop("Missing ABCD_SA12_CSV: ", sa12_csv)
SA12_10 <- read.csv(sa12_csv, stringsAsFactors = FALSE)

source(file.path(project_root, "gamfunction", "SCrankcorr.R"))

scanid_to_eventname <- function(scanID) {
  sess <- sub("^.*_ses-", "", as.character(scanID))
  sess <- gsub("([a-z])([A-Z])", "\\1_\\2", sess)
  sess <- gsub("([A-Za-z])([0-9])", "\\1_\\2", sess)
  sess <- gsub("([0-9])([A-Za-z])", "\\1_\\2", sess)
  tolower(sess)
}

vec_to_mat <- function(vec, ds = 12) {
  mat <- matrix(NA, ds, ds)
  idx <- which(lower.tri(mat, diag = TRUE))
  if (length(vec) > length(idx)) {
    stop("vec length exceeds lower-triangle size: ", length(vec), " > ", length(idx))
  }
  mat[idx[seq_along(vec)]] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  dimnames(mat) <- list(as.character(seq_len(ds)), as.character(seq_len(ds)))
  mat
}

plot_matrix <- function(mat, title, out_base) {
  df_melt <- as.data.frame(as.table(mat))
  names(df_melt) <- c("nodeid", "variable", "value")
  node_raw <- df_melt$nodeid
  var_raw <- df_melt$variable
  df_melt$nodeid <- suppressWarnings(as.numeric(as.character(node_raw)))
  df_melt$variable <- suppressWarnings(as.numeric(as.character(var_raw)))
  if (all(is.na(df_melt$nodeid))) df_melt$nodeid <- as.integer(node_raw)
  if (all(is.na(df_melt$variable))) df_melt$variable <- as.integer(var_raw)
  df_melt$nodeid <- -df_melt$nodeid
  df_melt$value <- as.numeric(df_melt$value)

  limthr <- if (all(is.na(df_melt$value))) 1 else max(abs(df_melt$value), na.rm = TRUE)
  if (!is.finite(limthr) || limthr == 0) limthr <- 1

  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5),
    ymin = rep(-12 - 0.5, times = 2),
    ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5),
    xmin = rep(0.5, times = 2),
    xmax = rep(12 + 0.5, times = 2)
  )

  p <- ggplot(data = df_melt) +
    geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
    scale_fill_distiller(type = "seq", palette = "RdBu", direction = -1, na.value = "grey", limits = c(-limthr, limthr)) +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, na.value = "grey", limits = c(-limthr, limthr)) +
    geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
    annotate("segment", x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5, color = "black", linewidth = 0.5) +
    ggtitle(label = title) +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(breaks = NULL, labels = NULL) +
    scale_x_continuous(breaks = NULL, labels = NULL) +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, angle = 315, hjust = 1, vjust = 1),
      axis.title = element_text(size = 18),
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0),
      panel.grid.minor = element_line(linewidth = 1)
    )

  ggsave(paste0(out_base, ".tiff"), p, height = 18, width = 20, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), p, height = 18, width = 20, units = "cm", bg = "transparent")
}

if (!force && file.exists(out_rds)) {
  message("[INFO] Found existing results, loading (set FORCE=1 to recompute): ", out_rds)
  res_df <- readRDS(out_rds)
} else {
  SCdata <- readRDS(input_rds)
  if (!("eventname" %in% names(SCdata)) && ("scanID" %in% names(SCdata))) {
    SCdata$eventname <- scanid_to_eventname(SCdata$scanID)
  }
  if (!("eventname" %in% names(SCdata))) {
    stop("Missing eventname (required to define baseline timepoint): input has no eventname/scanID")
  }

  SCdata$age <- as.numeric(SCdata$age)
  message(
    "[INFO] SCdata age range (years): ",
    round(min(SCdata$age, na.rm = TRUE), 3), "–", round(max(SCdata$age, na.rm = TRUE), 3)
  )

  needed <- c("subID", "age", "sex", "mean_fd")
  missing <- setdiff(needed, names(SCdata))
  if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
  SCdata$sex <- as.factor(SCdata$sex)

  sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
  if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
  if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
  sc_cols <- sc_cols[seq_len(78)]

  plot_fit <- plotdata$fit
  names(plot_fit) <- as.character(plotdata$SC_label)
  missing_fit <- setdiff(sc_cols, names(plot_fit))
  if (length(missing_fit) > 0) stop("plotdata missing fits for edges: ", paste(head(missing_fit, 10), collapse = ", "))

  # Normalize SC strength to ratio (divide by initial fit), consistent with LGCM cognition/pfactor scripts.
  for (edge in sc_cols) {
    f0 <- as.numeric(plot_fit[[edge]])
    if (is.na(f0) || !is.finite(f0) || f0 == 0) stop("Invalid plotdata fit for edge: ", edge)
    SCdata[[edge]] <- as.numeric(SCdata[[edge]]) / f0
  }

  # Build per-subject t0/t1 rows (min vs max age) and compute slope_per_year per edge.
  idx_by_sub <- split(seq_len(nrow(SCdata)), as.character(SCdata$subID))
  rows <- vector("list", length(idx_by_sub))
  for (k in seq_along(idx_by_sub)) {
    ii <- idx_by_sub[[k]]
    tsub <- SCdata$age[ii]
    if (length(unique(round(tsub, 6))) < 2) next
    i0 <- ii[which.min(tsub)]
    i1 <- ii[which.max(tsub)]
    dt <- SCdata$age[i1] - SCdata$age[i0]
    if (!is.finite(dt) || dt <= 0) next

    out <- data.frame(
      subID = SCdata$subID[i0],
      age_t0 = SCdata$age[i0],
      delta_age = dt,
      stringsAsFactors = FALSE
    )
    for (edge in sc_cols) {
      out[[paste0(edge, "_t0")]] <- SCdata[[edge]][i0]
      out[[paste0(edge, "_t1")]] <- SCdata[[edge]][i1]
    }
    rows[[k]] <- out
  }
  dat_sub <- dplyr::bind_rows(rows)
  if (nrow(dat_sub) < 10) stop("Too few subjects after two-timepoint selection: ", nrow(dat_sub))

  num_cores <- as.integer(Sys.getenv("SLOPE_CORES", unset = "16"))
  if (is.na(num_cores) || num_cores < 1) num_cores <- 16
  num_cores <- min(num_cores, parallel::detectCores())

  fit_edge <- function(i, data_sub, edges) {
    edge <- edges[[i]]
    t0_col <- paste0(edge, "_t0")
    t1_col <- paste0(edge, "_t1")
    slope <- (data_sub[[t1_col]] - data_sub[[t0_col]]) / data_sub$delta_age
    slope <- slope[is.finite(slope)]
    data.frame(
      edge = edge,
      n_sub = length(slope),
      personal_slope_mean = if (length(slope) > 0) mean(slope) else NA_real_,
      personal_slope_sd = if (length(slope) > 1) stats::sd(slope) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  message("[INFO] Computing SC personal slope per edge (LGCM-style slope_per_year; n_edges=78)")
  if (.Platform$OS.type == "windows") {
    message("[INFO] Windows parallel: ", num_cores, " workers")
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("dat_sub", "sc_cols", "fit_edge"),
      envir = environment()
    )
    res_list <- parallel::parLapply(cl, seq_along(sc_cols), fit_edge, data_sub = dat_sub, edges = sc_cols)
  } else {
    res_list <- lapply(seq_along(sc_cols), fit_edge, data_sub = dat_sub, edges = sc_cols)
  }

  res_df <- dplyr::bind_rows(res_list)
  saveRDS(res_df, out_rds)
  write.csv(res_df, out_csv, row.names = FALSE)
}

if (!all(c("edge", "personal_slope_mean") %in% names(res_df))) {
  stop("Invalid res_df loaded from: ", out_rds, "\nMissing required columns.")
}

message("[INFO] Matrix + S-A axis correlation (personal_slope_mean)")
slope_mat <- vec_to_mat(res_df$personal_slope_mean, ds = 12)
plot_matrix(
  slope_mat,
  "SC personal slope (mean slope_per_year)",
  file.path(FigureFolder, paste0("matrix_personal_slope_mean_SC_CV", CVthr))
)

SCrank.df <- SCrankcorr(res_df, "personal_slope_mean", 12, dsdata = FALSE)
saveRDS(SCrank.df, file.path(resultFolder, paste0("SCrankcorr_personal_slope_mean_SC_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr (personal slope mean) r=", round(SCrank.df$r.spearman, 3), " p=", signif(SCrank.df$p.spearman, 3))

SCrank.data <- SCrankcorr(res_df, "personal_slope_mean", 12, dsdata = TRUE)
limthr <- max(abs(SCrank.data$personal_slope_mean), na.rm = TRUE)
if (!is.finite(limthr) || limthr == 0) limthr <- 1

scatterFig <- ggplot(data = SCrank.data) +
  geom_point(aes(x = SCrank, y = personal_slope_mean, color = personal_slope_mean), size = 5) +
  geom_smooth(aes(x = SCrank, y = personal_slope_mean), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr, limthr)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 23, color = "black"),
    axis.title = element_text(size = 23),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "S-A axis rank", y = "Mean personal slope (per year)")

ggsave(
  file.path(FigureFolder, paste0("scatter_personal_slope_mean_vs_SCrank_SC_CV", CVthr, ".tiff")),
  scatterFig, height = 13, width = 13, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("scatter_personal_slope_mean_vs_SCrank_SC_CV", CVthr, ".pdf")),
  scatterFig, height = 13, width = 13, units = "cm", bg = "transparent"
)

