library(mgcv)
library(parallel)
library(ggplot2)
library(reshape)
library(RColorBrewer)
rm(list = ls())
CVthr <- 75

wdpath <- getwd()
if (grepl("cuizaixu_lab", wdpath, fixed = TRUE)) {
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD"
} else {
  interfileFolder <- file.path(wdpath, "wd", "interdataFolder_ABCD")
}
functionFolder <- file.path(wdpath, "gamfunction")
outputFolder <- file.path(wdpath, "outputs", "results", "cbcl_totprob")
combatFolder <- file.path(wdpath, "outputs", "results", "combat_cbcl")
figureFolder <- file.path(wdpath, "outputs", "figures", "cbcl_totprob")
figureFolderAssociation <- file.path(figureFolder, "Association")
figureFolderInteraction <- file.path(figureFolder, "Interaction")

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}
dir.create(figureFolderAssociation, showWarnings = FALSE, recursive = TRUE)
dir.create(figureFolderInteraction, showWarnings = FALSE, recursive = TRUE)

source(paste0(functionFolder, "/gamminteraction.R"))
source(paste0(functionFolder, "/SCrankcorr.R"))

input_rds <- file.path(
  wdpath, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cbcl.rds"
)
SCdata <- readRDS(input_rds)

need_demo <- (is.data.frame(SCdata$age) || !is.numeric(SCdata$age)) ||
  !all(c("handness", "race_ethnicity") %in% names(SCdata))
if (need_demo) {
  demodf <- read.csv(file.path(wdpath, "demopath", "DemodfScreenFinal.csv"))
  key <- match(SCdata$scanID, demodf$scanID)
  if (is.data.frame(SCdata$age) || !is.numeric(SCdata$age)) {
    SCdata$age <- demodf$age[key]
  }
  if (!"handness" %in% names(SCdata)) {
    SCdata$handness <- demodf$handness[key]
  }
  if (!"race_ethnicity" %in% names(SCdata)) {
    SCdata$race_ethnicity <- demodf$race_ethnicity[key]
  }
}
SCdata[, c("sex", "handness", "race_ethnicity")] <- lapply(SCdata[, c("sex", "handness", "race_ethnicity")], as.factor)
SCdata$age <- SCdata$age / 12
SCdata$cbcl_scr_syn_totprob_r <- as.numeric(SCdata$cbcl_scr_syn_totprob_r)
sc_cols <- grepl("SC\\.", names(SCdata)) & grepl("_h", names(SCdata))
SCdata$totalstrength <- rowMeans(SCdata[, sc_cols, drop = FALSE])

meandistance <- read.csv(paste0(interfileFolder, "/average_EuclideanDistance_12.csv"))
meandistance <- meandistance$Edistance

# Detect Windows system
is_windows <- .Platform$OS.type == "windows"

if (is_windows) {
  message("Windows system detected - skipping parallel processing and using existing files")
  cores <- 1
} else {
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (is.na(cores) || cores < 1) {
    cores <- 1
  }
  cores <- min(cores, 50)
}

## 1. CBCL total raw association
dataname <- "SCdata"
smooth_var <- "age"
int_var.predict.percentile <- 0.1
covariates <- "sex+mean_fd"
knots <- 3
set_fx <- TRUE
increments <- 1000
stats_only <- TRUE
int_var <- "cbcl_scr_syn_totprob_r"

# Check if result file exists (for Windows)
gamresult_file <- file.path(outputFolder, paste0("gamresult_Int_age_cbcl_totprob_raw_CV", CVthr, ".rds"))

if (is_windows && file.exists(gamresult_file)) {
  message("Using existing CBCL association results file on Windows")
  gamresult.tmp <- readRDS(gamresult_file)
} else {
  resultsum <- mclapply(1:78, function(x) {
    region <- grep("SC\\.", names(SCdata), value = TRUE)[x]
    gamresult <- gamm.smooth.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                          int_var.predict.percentile, covariates,
                                                          knots, set_fx, increments, stats_only)
    gamresult <- as.data.frame(gamresult)
    return(gamresult)
  }, mc.cores = cores)
  
  gamresult.tmp <- do.call(rbind, resultsum)
  gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
  gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
  gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")
  
  print(paste0(sum(gamresult.tmp$bootstrap_pvalue.fdr < 0.05), " edges have significant age by ", int_var, " effect."))
  print(paste0(sum(gamresult.tmp$bootstrap.P.disease.fdr < 0.05), " edges have significant ", int_var, " effect."))
  
  saveRDS(gamresult.tmp, gamresult_file)
}

## 2. Correlation to S-A connectional axis
# Check if gamresult file exists, if not try to read from the variable created above
if (exists("gamresult.tmp")) {
  message("Using gamresult.tmp from previous step")
} else {
  gamresult_file <- file.path(outputFolder, paste0("gamresult_Int_age_cbcl_totprob_raw_CV", CVthr, ".rds"))
  if (file.exists(gamresult_file)) {
    gamresult.tmp <- readRDS(gamresult_file)
    message("Reading existing gamresult file")
  } else {
    stop("No gamresult data available. Please run the analysis on a Linux system first.")
  }
}
gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")

SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = FALSE)
SCrank.df.cbcl <- SCrankcorr(gamresult.tmp, "T.disease", 12, dsdata = FALSE)
SCrank.df <- rbind(SCrank.df.age, SCrank.df.cbcl)
SCrank.df$int_var <- int_var

#  print correlation
cat("\n=== correlation results ===\n")
print(SCrank.df)
cat("\n")


SCrank.tmp <- SCrankcorr(gamresult.tmp, "IntpartialRsq", 12, dsdata = TRUE)
gamresult.tmp$SCrank <- SCrank.tmp$SCrank

print("Next, correlation between CBCL associations and connectional axis is tested while controlling for Euclidean distance.")
gamresult.tmp$meandistance <- meandistance
corr.test(gamresult.tmp$T.disease, gamresult.tmp$meandistance, method = "pearson")
gamresult.tmp$T.disease_control_distance[which(!is.na(gamresult.tmp$T.disease))] <- residuals(lm(T.disease ~ meandistance, data = gamresult.tmp))
corr.test(gamresult.tmp$T.disease_control_distance, gamresult.tmp$meandistance, method = "pearson")
SCrankresult.whole.controldistance <- SCrankcorr(gamresult.tmp, "T.disease_control_distance", 12, dsdata = FALSE)
print(paste("Correlation coefficient between CBCL associations regressing out fiber distance and connectional axis is",
            round(SCrankresult.whole.controldistance$r.spearman, 2), "with a P value of",
            round(SCrankresult.whole.controldistance$p.spearman, 3)))
print(SCrankresult.whole.controldistance)

SCrank.df.controldistance <- SCrankresult.whole.controldistance
SCrank.df.controldistance$Interest.var <- "T.disease_control_distance"
SCrank.df.controldistance$int_var <- "cbcl_scr_syn_totprob_r"
SCrank.df.all <- rbind(SCrank.df, SCrank.df.controldistance)

print(SCrank.df.all)
cat("\n")

saveRDS(SCrankresult.whole.controldistance,
        file.path(outputFolder, paste0("SCrankcorr_cbcl_totprob_raw_CV", CVthr, ".rds")))

## 3. Plots: t-value matrix and S-A rank correlation
plot_vars <- c("IntpartialRsq", "T.disease", "T.disease_control_distance")
for (Interest.var in plot_vars) {
  message(paste("Processing variable:", Interest.var))
  tmpvar <- gamresult.tmp[, Interest.var]
  message(paste("Variable range:", min(tmpvar, na.rm = TRUE), "to", max(tmpvar, na.rm = TRUE)))
  limthr <- max(abs(tmpvar), na.rm = TRUE)
  ytitle <- Interest.var
  if (grepl("T.disease", Interest.var, fixed = TRUE)) {
    ytitle <- expression("CBCL score association (" * italic("T") * " value)")
  }

  if (grepl("_control_distance", Interest.var, fixed = TRUE)) {
    mytheme <- theme(axis.text = element_text(size = 23, color = "black"),
                     axis.title = element_text(size = 23),
                     aspect.ratio = 0.97,
                     axis.line = element_line(linewidth = 0.6),
                     axis.ticks = element_line(linewidth = 0.6),
                     plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
                     plot.background = element_rect(fill = "transparent"),
                     panel.background = element_rect(fill = "transparent"),
                     legend.position = "none")
    mywidth <- 15
    myheight <- 15
  } else {
    mytheme <- theme(axis.text = element_text(size = 23.4, color = "black"),
                     axis.title = element_text(size = 23.4),
                     aspect.ratio = 1,
                     axis.line = element_line(linewidth = 0.6),
                     axis.ticks = element_line(linewidth = 0.6),
                     plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
                     plot.background = element_rect(fill = "transparent"),
                     panel.background = element_rect(fill = "transparent"),
                     legend.position = "none")
    mywidth <- 18
    myheight <- 18
  }

  scatterFig <- ggplot(data = gamresult.tmp) +
    geom_point(aes(x = SCrank, y = tmpvar, color = tmpvar), size = 5) +
    geom_smooth(aes(x = SCrank, y = tmpvar), linewidth = 1.4, method = "lm", color = "black") +
    scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr, limthr)) +
    labs(x = "S-A connectional axis rank", y = ytitle) +
    theme_classic() +
    mytheme +
    theme(axis.title.y = element_text(margin = margin(r = 10)))
  ggsave(file.path(figureFolderAssociation, paste0(Interest.var, "_", int_var, "_SCrankcorr.tiff")),
         scatterFig, width = mywidth, height = myheight, units = "cm")
  scatter_filename <- file.path(figureFolderAssociation, paste0(Interest.var, "_", int_var, "_SCrankcorr.svg"))
  message(paste("Saving scatter plot:", scatter_filename))
  ggsave(scatter_filename, scatterFig, device = grDevices::svg, width = mywidth, height = myheight, units = "cm")

  ## Matrix plot generation for each variable
  Matrix.tmp.T <- matrix(NA, 12, 12)
  Matrix.tmp.T[lower.tri(Matrix.tmp.T, diag = TRUE)] <- tmpvar
  Matrix.tmp.T[upper.tri(Matrix.tmp.T)] <- t(Matrix.tmp.T)[upper.tri(Matrix.tmp.T)]
  colnames(Matrix.tmp.T) <- seq_len(12)
  rownames(Matrix.tmp.T) <- seq_len(12)
  matrixtmp.df <- as.data.frame(Matrix.tmp.T)
  matrixtmp.df$nodeid <- seq_len(12)
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
  matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
  matrixtmp.df.sig$nodeid <- seq_len(12)
  matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig, id.vars = c("nodeid"))
  matrixtmp.df.sig.melt$variable <- as.numeric(matrixtmp.df.sig.melt$variable)
  matrixtmp.df.sig.melt$nodeid <- 0 - matrixtmp.df.sig.melt$nodeid
  matrixtmp.df.sig.melt$value <- as.logical(matrixtmp.df.sig.melt$value)
  matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[which(matrixtmp.df.sig.melt$value), ]

  linerange_frame <- data.frame(
    x = c(0.5, 12.5),
    ymin = rep(-12.5, times = 2),
    ymax = rep(-0.5, times = 2),
    y = c(-0.5, -12.5),
    xmin = rep(0.5, times = 2),
    xmax = rep(12.5, times = 2)
  )
  MatFig <- ggplot(data = matrixtmp.df.melt) +
    geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
    scale_fill_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    scale_color_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    geom_text(data = matrixtmp.df.sig.melt, aes(x = variable, y = nodeid, label = "*"),
              vjust = 0.7, hjust = 0.5, size = 8) +
    geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0.5, y = -0.5, xend = 12.5, yend = -12.5), color = "black", linewidth = 0.5) +
    ggtitle(paste0(int_var, "_", Interest.var)) +
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
  matrix_filename <- file.path(figureFolderAssociation, paste0(Interest.var, "_", int_var, "_Matrix12.tiff"))
  message(paste("Saving matrix plot:", matrix_filename))
  ggsave(matrix_filename, MatFig, height = 18, width = 20, units = "cm")
  
} # End of for loop for plot_vars

## 4. Developmental trajectories: high (90th) vs low (10th) CBCL across five S-A deciles
sa_path <- file.path(interfileFolder, "SA12_10.csv")
plotdata_path <- file.path(interfileFolder, paste0("plotdatasum.df_SA12_sumSCinvnode_siteall_CV", CVthr, ".rds"))
has_scaling <- file.exists(plotdata_path)
SCdata.diw <- SCdata
if (has_scaling) {
  plotdata.ref <- readRDS(plotdata_path)
  for (x in 1:78) {
    region <- grep("SC\\.", names(SCdata), value = TRUE)[x]
    plotdata.tmp <- plotdata.ref[plotdata.ref$SC_label == region, ]
    if (nrow(plotdata.tmp) > 0 && is.numeric(plotdata.tmp$fit[1]) && !is.na(plotdata.tmp$fit[1]) && plotdata.tmp$fit[1] != 0) {
      SCdata.diw[, region] <- SCdata[, region] / plotdata.tmp$fit[1]
    }
  }
}
SCdata.diw[, grep("SC\\.", names(SCdata.diw), value = TRUE)] <- lapply(SCdata.diw[, grep("SC\\.", names(SCdata.diw), value = TRUE)], as.numeric)
assign("SCdata.diw", SCdata.diw, envir = .GlobalEnv)

trajectory_cache <- file.path(outputFolder, paste0("plotdata_high90_low10_cbcl_totprob_develop_CV", CVthr, ".rds"))

if (is_windows && file.exists(trajectory_cache)) {
  message("Using existing developmental trajectory results file on Windows")
} else if (!file.exists(trajectory_cache)) {
  stats_only <- FALSE
  increments <- 200
  dataname <- "SCdata.diw"
  resultsum <- mclapply(1:78, function(x) {
    region <- grep("SC\\.", names(SCdata.diw), value = TRUE)[x]
    int_var.predict.percentile <- 0.1
    result.all <- gamm.smooth.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                           int_var.predict.percentile, covariates,
                                                           knots, set_fx, increments, stats_only)
    plotdata.low <- result.all[[2]]
    plotdata.low$label <- "low"
    plotdata.low$SC_label <- region

    int_var.predict.percentile <- 0.9
    result.all <- gamm.smooth.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                           int_var.predict.percentile, covariates,
                                                           knots, set_fx, increments, stats_only)
    plotdata.high <- result.all[[2]]
    plotdata.high$label <- "high"
    plotdata.high$SC_label <- region

    rbind(plotdata.low, plotdata.high)
  }, mc.cores = min(cores, 20))
  saveRDS(resultsum, trajectory_cache)
}

SA12_10 <- read.csv(sa_path)

# Check if trajectory cache exists
if (!file.exists(trajectory_cache)) {
  stop("Trajectory data file not found. Please run the analysis on a Linux system first to generate the trajectory data.")
}

plotdata <- do.call(rbind, readRDS(trajectory_cache))
plotdata <- merge(plotdata, SA12_10[, c("SC_label", "decile")], by = "SC_label")

fit_col <- if (".fitted" %in% names(plotdata)) ".fitted" else "fitted.centered"
plotdata$fit_value <- plotdata[, fit_col]
plotdf.decile <- aggregate(fit_value ~ decile + age + label, data = plotdata, FUN = mean)

colorid <- rev(brewer.pal(10, "RdBu"))
for (i in 1:10) {
  plotdf.tmp <- plotdf.decile[plotdf.decile$decile == i, ]
  colorindex <- colorid[i]

  if (i == 1 || i == 6) {
    mytheme <- theme(axis.text = element_text(size = 21, color = "black"),
                     axis.title = element_text(size = 21),
                     aspect.ratio = 1,
                     axis.line = element_line(linewidth = 0.5),
                     axis.ticks = element_line(linewidth = 0.5),
                     plot.background = element_rect(fill = "transparent"),
                     panel.border = element_rect(fill = NA, color = "transparent"),
                     panel.background = element_rect(fill = "transparent", color = "transparent"),
                     legend.position = "none")
  } else {
    mytheme <- theme(axis.text.x = element_text(size = 21, color = "black"),
                     axis.text.y = element_text(size = 21, color = "transparent"),
                     axis.title.x = element_text(size = 21),
                     axis.title.y = element_text(size = 21, colour = "transparent"),
                     aspect.ratio = 1,
                     axis.line.x = element_line(linewidth = 0.5),
                     axis.line.y = element_line(linewidth = 0.5, colour = "transparent"),
                     axis.ticks.x = element_line(linewidth = 0.5),
                     axis.ticks.y = element_line(linewidth = 0.5, colour = "transparent"),
                     plot.background = element_rect(fill = "transparent"),
                     panel.grid = element_line(linewidth = 0.5, colour = "transparent"),
                     panel.border = element_rect(fill = NA, color = "transparent"),
                     panel.background = element_rect(fill = "transparent", color = "transparent"),
                     legend.position = "none")
  }

  Fig <- ggplot(data = plotdf.tmp) +
    geom_line(aes(x = age, y = fit_value, group = label, linetype = label), linewidth = 1.2, color = colorindex) +
    scale_y_continuous(breaks = c(0.9, 1.0), limits = c(0.90, 1.1)) +
    scale_linetype_manual(values = c(high = "dashed", low = "solid")) +
    labs(x = NULL, y = "SC strength (ratio)") +
    mytheme

  ggsave(file.path(figureFolderInteraction, paste0("developmentcurve_decile", i, ".tiff")),
         Fig, width = 10, height = 10, units = "cm")
  ggsave(file.path(figureFolderInteraction, paste0("developmentcurve_decile", i, ".svg")),
         Fig, device = grDevices::svg, width = 10, height = 10, units = "cm")
}
