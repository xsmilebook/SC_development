#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(tidyverse)
  library(reshape)
  library(RColorBrewer)
})

rm(list = ls())

CVthr <- 75
int_var <- "GENERAL"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "pfactor")
FigureFolder <- file.path(project_root, "outputs", "figures", "6th_pfactor", "abcd", "pfactor")
intermediateFolder <- file.path(project_root, "outputs", "intermediate", "6th_pfactor", "abcd", "pfactor")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(intermediateFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
}

euclid_csv <- Sys.getenv(
  "ABCD_EUCLID_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "average_EuclideanDistance_12.csv")
)
if (!file.exists(euclid_csv)) {
  stop("Missing ABCD_EUCLID_CSV: ", euclid_csv)
}
meandistance <- read.csv(euclid_csv)$Edistance

sa12_csv <- Sys.getenv(
  "ABCD_SA12_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "SA12_10.csv")
)
if (!file.exists(sa12_csv)) {
  stop("Missing ABCD_SA12_CSV: ", sa12_csv)
}
SA12_10 <- read.csv(sa12_csv, stringsAsFactors = FALSE)

plotdatasum_rds <- Sys.getenv(
  "ABCD_PLOTDATASUM_RDS",
  unset = "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds"
)
if (!file.exists(plotdatasum_rds)) {
  stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)
}
plotdata <- readRDS(plotdatasum_rds)

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "gamminteraction.R"))

SCdata <- readRDS(input_rds)
age_to_years <- function(age_raw) {
  age_num <- as.numeric(age_raw)
  mx <- suppressWarnings(max(age_num, na.rm = TRUE))
  if (is.finite(mx) && mx > 24) {
    age_num / 12
  } else {
    age_num
  }
}
SCdata$age <- age_to_years(SCdata$age)

needed_base <- c("subID", "scanID", "siteID", "age", "sex", "mean_fd", int_var)
missing_base <- setdiff(needed_base, names(SCdata))
if (length(missing_base) > 0) {
  stop("Missing required columns: ", paste(missing_base, collapse = ", "))
}

SCdata[, c("sex")] <- lapply(SCdata[, c("sex"), drop = FALSE], as.factor)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
sc_cols <- sc_cols[str_detect(sc_cols, "_h$")]
if (length(sc_cols) < 78) {
  stop("Expected >=78 ComBat SC.*_h columns, got: ", length(sc_cols))
}

SCdata$totalstrength <- rowMeans(SCdata[, sc_cols, drop = FALSE])

dataname <- "SCdata"
smooth_var <- "age"
covariates <- "sex+mean_fd"
knots <- 3
set_fx <- TRUE
increments <- 1000
int_var_predict_percentile <- 0.1
stats_only <- TRUE

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "50"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 50
num_cores <- min(num_cores, 50L)

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("gamresult_Int_age_pFactor_", int_var, "_CV", CVthr, ".rds"))

if (force || !file.exists(out_rds)) {
  message("[INFO] Fitting p-factor GAMM interaction models (n_edges=78, mc.cores=", num_cores, ")")
  resultsum <- parallel::mclapply(seq_len(78), function(x) {
    region <- sc_cols[[x]]
    gamresult <- gamm.smooth.predict.covariateinteraction(
      region,
      dataname,
      smooth_var,
      int_var,
      int_var_predict_percentile,
      covariates,
      knots,
      set_fx,
      increments,
      stats_only = stats_only
    )
    as.data.frame(gamresult)
  }, mc.cores = num_cores)

  gamresult.tmp <- do.call(rbind, resultsum)
  gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
  gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
  gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

  message(sum(gamresult.tmp$bootstrap_pvalue.fdr < 0.05), " edges have significant age by ", int_var, " effect.")
  message(sum(gamresult.tmp$bootstrap.P.disease.fdr < 0.05), " edges have significant ", int_var, " effect.")

  saveRDS(gamresult.tmp, out_rds)
} else {
  gamresult.tmp <- readRDS(out_rds)
}

gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

message("[INFO] Correlation to connectional axis")
SCrank.tmp <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = TRUE)
SCrank <- SCrank.tmp$SCrank
gamresult.tmp$SCrank <- SCrank
SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = FALSE)
SCrank.df.general <- SCrankcorr(gamresult.tmp, "T.disease", 12, dsdata = FALSE)
SCrank.df <- rbind(SCrank.df.age, SCrank.df.general)
SCrank.df$int_var <- int_var

message("[INFO] Control Euclidean distance for T.disease")
gamresult.tmp$meandistance <- meandistance
gamresult.tmp$T.disease_control_distance[which(!is.na(gamresult.tmp$T.disease))] <-
  residuals(lm(T.disease ~ meandistance, data = gamresult.tmp))
SCrank.df.general.controldistance <- SCrankcorr(gamresult.tmp, "T.disease_control_distance", 12, dsdata = FALSE)

saveRDS(
  list(
    SCrank.df = SCrank.df,
    SCrank.df.general.controldistance = SCrank.df.general.controldistance
  ),
  file.path(resultFolder, paste0("SCrankcorr_summary_Int_age_pFactor_", int_var, "_CV", CVthr, ".rds"))
)

save_svg_or_pdf <- function(filename_svg, plot_obj, width, height, units = "cm") {
  ok <- TRUE
  tryCatch(
    {
      ggsave(filename_svg, plot_obj, dpi = 600, width = width, height = height, units = units, bg = "transparent")
    },
    error = function(e) {
      ok <<- FALSE
      message("[WARN] svg save failed: ", filename_svg, " (", conditionMessage(e), "); saving pdf instead")
    }
  )
  if (!ok) {
    filename_pdf <- sub("\\.svg$", ".pdf", filename_svg)
    ggsave(filename_pdf, plot_obj, dpi = 600, width = width, height = height, units = units, bg = "transparent")
  }
}

## ---- Figures: scatter plots + matrix graphs (match the original Rmd) ----
dir.create(file.path(FigureFolder, "Disease", "pFactor"), showWarnings = FALSE, recursive = TRUE)
for (Interest.var in c("IntpartialRsq", "T.disease", "T.disease_control_distance")) {
  tmpvar <- gamresult.tmp[, Interest.var]
  limthr <- max(abs(gamresult.tmp[, Interest.var]), na.rm = TRUE)
  if (str_detect(Interest.var, "T\\.disease") && int_var == "GENERAL") {
    ytitle <- expression(italic("p") * "-factor association (" * italic("T") * " value)")
  } else {
    ytitle <- Interest.var
  }

  if (str_detect(Interest.var, "_control_distance")) {
    mytheme <- theme(
      axis.text = element_text(size = 23, color = "black"),
      axis.title = element_text(size = 23),
      aspect.ratio = 0.97,
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),
      plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none"
    )
    mywidth <- 13
    myheight <- 13
  } else {
    mytheme <- theme(
      axis.text = element_text(size = 23.4, color = "black"),
      axis.title = element_text(size = 23.4),
      aspect.ratio = 1,
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),
      plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none"
    )
    mywidth <- 15
    myheight <- 15
  }

  scatterFig <- ggplot(data = gamresult.tmp) +
    geom_point(aes(x = SCrank, y = tmpvar, color = tmpvar), size = 5) +
    geom_smooth(aes(x = SCrank, y = tmpvar), linewidth = 1.4, method = "lm", color = "black") +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr, limthr)) +
    labs(x = "S-A connectional axis rank", y = ytitle) +
    theme_classic() +
    mytheme

  out_base <- file.path(FigureFolder, "Disease", "pFactor", paste0(Interest.var, "_", int_var, "_SCrankcorr"))
  ggsave(paste0(out_base, ".tiff"), scatterFig, width = mywidth, height = myheight, units = "cm", bg = "transparent")
  save_svg_or_pdf(paste0(out_base, ".svg"), scatterFig, width = mywidth, height = myheight, units = "cm")

  Matrix.tmp.T <- matrix(NA, 12, 12)
  Matrix.tmp.T[lower.tri(Matrix.tmp.T, diag = TRUE)] <- tmpvar
  Matrix.tmp.T[upper.tri(Matrix.tmp.T)] <- t(Matrix.tmp.T)[upper.tri(Matrix.tmp.T)]
  colnames(Matrix.tmp.T) <- seq(1, 12)
  rownames(Matrix.tmp.T) <- seq(1, 12)
  matrixtmp.df <- as.data.frame(Matrix.tmp.T)
  matrixtmp.df$nodeid <- seq(1, 12)
  matrixtmp.df.melt <- melt(matrixtmp.df, id.vars = c("nodeid"))
  matrixtmp.df.melt$variable <- as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid <- 0 - matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value <- as.numeric(matrixtmp.df.melt$value)

  Matrix.tmp.sig <- matrix(NA, 12, 12)
  if (Interest.var == "IntpartialRsq") {
    Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = TRUE)] <- (gamresult.tmp$bootstrap_pvalue.fdr < 0.05)
  } else {
    Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = TRUE)] <- (gamresult.tmp$bootstrap.P.disease.fdr < 0.05)
  }
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

  titlematrix <- paste0(int_var, "_", Interest.var)
  linerange_frame <- data.frame(
    x = c(0.5, 12 + 0.5),
    ymin = rep(-12 - 0.5, times = 2),
    ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12 - 0.5),
    xmin = rep(0.5, times = 2),
    xmax = rep(12 + 0.5, times = 2)
  )

  MatFig <- ggplot(data = matrixtmp.df.melt) +
    geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
    scale_fill_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    scale_color_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    geom_text(data = matrixtmp.df.sig.melt, aes(x = variable, y = nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size = 8) +
    geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
    ggtitle(label = titlematrix) +
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

  filename <- file.path(FigureFolder, "Disease", "pFactor", paste0(Interest.var, "_", int_var, "_Matrix12.tiff"))
  ggsave(filename, MatFig, height = 18, width = 20, units = "cm", bg = "transparent")
}

## ---- Figures: interaction plot (development curves by p-factor level) ----
message("[INFO] Building interaction-curve figures (high vs low p-factor)")
dir.create(file.path(FigureFolder, "Disease", "pFactor", "Interaction"), showWarnings = FALSE, recursive = TRUE)

# Scale SC strength by their initial strength in the age span (as in the original Rmd).
SCdata.diw <- SCdata
for (region in sc_cols[seq_len(78)]) {
  plotdata.tmp <- plotdata[plotdata$SC_label == region, , drop = FALSE]
  if (!("fit" %in% names(plotdata.tmp)) || nrow(plotdata.tmp) < 1 || is.na(plotdata.tmp$fit[[1]]) || plotdata.tmp$fit[[1]] == 0) {
    message("[WARN] Missing/invalid plotdata fit for ", region, "; skip scaling for this edge")
    next
  }
  SCdata.diw[[region]] <- SCdata[[region]] / plotdata.tmp$fit[[1]]
}
SCdata.diw[, sc_cols] <- lapply(SCdata.diw[, sc_cols, drop = FALSE], as.numeric)

dataname <- "SCdata.diw"
stats_only <- FALSE
cache_rds <- file.path(intermediateFolder, paste0("plotdata_high90_low10_pFactor_", int_var, "_develop_CV", CVthr, ".rds"))
if (force || !file.exists(cache_rds)) {
  resultsum <- vector("list", length = 78)
  for (x in 1:78) {
    region <- sc_cols[[x]]

    int_var_predict_percentile <- 0.1
    result.all <- gamm.smooth.predict.covariateinteraction(
      region, dataname, smooth_var, int_var, int_var_predict_percentile,
      covariates, knots, set_fx, increments, stats_only = stats_only
    )
    plotdata.low <- result.all[[2]]
    plotdata.low$pFactor <- int_var
    plotdata.low$label <- "low"

    int_var_predict_percentile <- 0.9
    result.all <- gamm.smooth.predict.covariateinteraction(
      region, dataname, smooth_var, int_var, int_var_predict_percentile,
      covariates, knots, set_fx, increments, stats_only = stats_only
    )
    plotdata.high <- result.all[[2]]
    plotdata.high$pFactor <- int_var
    plotdata.high$label <- "high"

    plotdata.edge <- rbind(plotdata.low, plotdata.high)
    plotdata.edge$SC_label <- region
    resultsum[[x]] <- plotdata.edge
  }
  saveRDS(resultsum, cache_rds)
}

resultsum <- readRDS(cache_rds)
plotdata.all <- do.call(rbind, resultsum)
plotdata.all <- merge(plotdata.all, SA12_10, by = "SC_label")

plotdf.decile.low <- plotdata.all %>%
  filter(label == "low") %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(.fitted), decile = mean(decile), .groups = "drop")
plotdf.decile.low$label <- "low"

plotdf.decile.high <- plotdata.all %>%
  filter(label == "high") %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(.fitted), decile = mean(decile), .groups = "drop")
plotdf.decile.high$label <- "high"

plotdf.decile <- rbind(plotdf.decile.low, plotdf.decile.high)

colorid <- rev(brewer.pal(10, "RdBu"))
for (i in 1:10) {
  plotdf.tmp <- plotdf.decile[plotdf.decile$decile == i, , drop = FALSE]
  colorindex <- colorid[i]
  age_label_mult <- if (max(plotdf.tmp$age, na.rm = TRUE) <= 2) 10 else 1
  if (i == 1 || i == 6) {
    mytheme <- theme(
      axis.text = element_text(size = 21, color = "black"),
      axis.title = element_text(size = 21),
      aspect.ratio = 1,
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.border = element_rect(fill = NA, color = "transparent"),
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "none"
    )
  } else {
    mytheme <- theme(
      axis.text.x = element_text(size = 21, color = "black"),
      axis.text.y = element_text(size = 21, color = "transparent"),
      axis.title.x = element_text(size = 21),
      axis.title.y = element_text(size = 21, colour = "transparent"),
      aspect.ratio = 1,
      axis.line.x = element_line(linewidth = 0.5),
      axis.line.y = element_line(linewidth = 0.5, colour = "transparent"),
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.ticks.y = element_line(linewidth = 0.5, colour = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_line(linewidth = 0.5, colour = "transparent"),
      panel.border = element_rect(fill = NA, color = "transparent"),
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "none"
    )
  }

  Fig <- ggplot(data = plotdf.tmp) +
    geom_line(aes(x = age, y = fit.avg, group = label, linetype = label), linewidth = 1.2, color = colorindex) +
    scale_x_continuous(labels = function(x) x * age_label_mult) +
    scale_y_continuous(breaks = c(0.9, 1.0), limits = c(0.90, 1.1)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = NULL, y = "SC strength (ratio)") +
    mytheme

  out_base <- file.path(FigureFolder, "Disease", "pFactor", "Interaction", paste0("developmentcurve_decile", i))
  ggsave(paste0(out_base, ".tiff"), Fig, width = 10, height = 10, units = "cm", bg = "transparent")
  save_svg_or_pdf(paste0(out_base, ".svg"), Fig, width = 10, height = 10, units = "cm")
}
