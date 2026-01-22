#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(psych)
  library(reshape)
  library(RColorBrewer)
  library(tidyverse)
})

rm(list = ls())

CVthr <- 75
Cogvar <- "nihtbx_fluidcomp_uncorrected"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "cognition")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "cognition")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (cognition variant)")
}

euclid_csv <- Sys.getenv(
  "ABCD_EUCLID_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "average_EuclideanDistance_12.csv")
)
if (!file.exists(euclid_csv)) {
  stop("Missing ABCD_EUCLID_CSV: ", euclid_csv)
}

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "gamcog.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

# `gamcog.R` loads `ecostats` -> `mvabund`, which provides a `theme()` function
# that can mask `ggplot2::theme()` and break plotting. Pin ggplot2 functions.
theme <- ggplot2::theme
theme_classic <- ggplot2::theme_classic

SCdata <- readRDS(input_rds)
meandistance <- read.csv(euclid_csv)$Edistance
SCdata$age <- as.numeric(SCdata$age) / 12

if (!Cogvar %in% names(SCdata)) {
  stop("Missing cognition variable in input: ", Cogvar)
}

nonna_index <- which(!is.na(SCdata[, Cogvar]))
SCdata.cog <- SCdata[nonna_index, , drop = FALSE]
if ("eventname" %in% names(SCdata.cog)) {
  SCdata.cog <- SCdata.cog[SCdata.cog$eventname == "baseline_year_1_arm_1", , drop = FALSE]
}

cogagemodel <- gam(stats::as.formula(paste0(Cogvar, "~ s(age,k=3, fx=TRUE)+sex+mean_fd")), data = SCdata.cog)
t <- summary(cogagemodel)
message("age, sex, mean_fd can explain ", round(t$r.sq, 3), " variance of cognition.")

dataname <- "SCdata.cog"
smooth_var <- "age"
covariates <- "sex+mean_fd"
knots <- 3
corrmethod <- "pearson"

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "60"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 60
num_cores <- min(num_cores, 60L)

run_mclapply_with_fallback <- function(x, fun, n_workers) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1) n_workers <- 1L
  tries <- c(n_workers, 60L, 50L, 40L, 30L, 20L, 16L, 8L, 4L, 2L, 1L)
  tries <- unique(tries[tries <= n_workers])
  for (k in tries) {
    message("[INFO] Trying mclapply with mc.cores=", k)
    out <- try(parallel::mclapply(x, fun, mc.cores = k), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    msg <- as.character(out)
    message("[WARN] mclapply failed: ", msg)
    message("[WARN] Retrying with fewer workers.")
  }
  stop("mclapply failed for all worker settings; likely process/memory limits.")
}

make_error_row <- function(parcel, cognition_var, err) {
  data.frame(
    parcel = parcel,
    cognition_var = cognition_var,
    gam.smooth.t = NA_real_,
    gam.smooth.pvalue = NA_real_,
    anova.cov.pvalue = NA_real_,
    anova.cov.int.pvalue = NA_real_,
    partialRsq.int = NA_real_,
    partialRsq = NA_real_,
    correstimate = NA_real_,
    corrp = NA_real_,
    error = as.character(err),
    stringsAsFactors = FALSE
  )
}

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("SC_Cog_results_", Cogvar, "_CV", CVthr, "_cognition.rds"))

if (force || !file.exists(out_rds)) {
  sc_labels <- grep("SC\\.", names(SCdata), value = TRUE)
  if (length(sc_labels) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_labels))

  resultsum <- run_mclapply_with_fallback(seq_len(78), function(i) {
    SClabel <- sc_labels[[i]]
    tryCatch(
      {
        gamresult <- gam.fit.cognition(
          SClabel,
          dataname,
          Cogvar,
          smooth_var,
          covariates,
          knots,
          corrmethod,
          set_fx = TRUE,
          stats_only = TRUE
        )
        list(ok = TRUE, row = as.data.frame(gamresult), err = NA_character_)
      },
      error = function(e) {
        list(ok = FALSE, row = NULL, err = conditionMessage(e))
      }
    )
  }, num_cores)

  ok_mask <- vapply(resultsum, function(z) isTRUE(z$ok), logical(1))
  if (!any(ok_mask)) {
    errs <- unique(vapply(resultsum, function(z) as.character(z$err), character(1)))
    stop("All edges failed inside mclapply. First errors:\n", paste(head(errs, 10), collapse = "\n"))
  }

  rows <- vector("list", length = 78)
  for (i in seq_len(78)) {
    if (isTRUE(resultsum[[i]]$ok)) {
      rows[[i]] <- resultsum[[i]]$row
    } else {
      rows[[i]] <- make_error_row(sc_labels[[i]], Cogvar, resultsum[[i]]$err)
    }
  }
  SC_Cog_results.df <- do.call(rbind, rows)

  if ("error" %in% names(SC_Cog_results.df)) {
    failed_n <- sum(!is.na(SC_Cog_results.df$error) & nzchar(SC_Cog_results.df$error))
    message("[INFO] Edge-level failures: ", failed_n, " / 78")
  }

  SC_Cog_results.df[, c(3:10)] <- lapply(SC_Cog_results.df[, c(3:10)], as.numeric)

  if (!"corrp" %in% names(SC_Cog_results.df)) stop("Missing column: corrp")
  if (!"anova.cov.pvalue" %in% names(SC_Cog_results.df)) stop("Missing column: anova.cov.pvalue")
  if (!"gam.smooth.pvalue" %in% names(SC_Cog_results.df)) stop("Missing column: gam.smooth.pvalue")

  SC_Cog_results.df$corr.p.fdr <- p.adjust(as.numeric(SC_Cog_results.df$corrp), method = "fdr")
  SC_Cog_results.df$anova.cov.p.fdr <- p.adjust(as.numeric(SC_Cog_results.df$anova.cov.pvalue), method = "fdr")
  SC_Cog_results.df$gam.smooth.p.fdr <- p.adjust(as.numeric(SC_Cog_results.df$gam.smooth.pvalue), method = "fdr")

  saveRDS(SC_Cog_results.df, out_rds)
} else {
  SC_Cog_results.df <- readRDS(out_rds)
}

message(sum(SC_Cog_results.df$anova.cov.p.fdr < 0.05), " edges have significant associations with ", Cogvar, ".")

SC_Cog_results.df.whole <- SC_Cog_results.df
SCrankresult.whole <- SCrankcorr(SC_Cog_results.df.whole, "gam.smooth.t", 12)
message(
  "Correlation coefficient between cognitive associations and connectional axis is ",
  round(SCrankresult.whole$r.spearman, 2), " with P=", round(SCrankresult.whole$p.spearman, 3)
)

message("Control Euclidean distance.")
SC_Cog_results.df.whole$meandistance <- meandistance
SC_Cog_results.df.whole$gam.smooth.t_control_distance[which(!is.na(SC_Cog_results.df.whole$gam.smooth.t))] <-
  residuals(lm(gam.smooth.t ~ meandistance, data = SC_Cog_results.df.whole))
SCrankresult.whole.controllength <- SCrankcorr(SC_Cog_results.df.whole, "gam.smooth.t_control_distance", 12, dsdata = FALSE)

saveRDS(
  list(
    SCrankresult.whole = SCrankresult.whole,
    SCrankresult.whole.controllength = SCrankresult.whole.controllength
  ),
  file.path(resultFolder, paste0("SCrankcorr_summary_", Cogvar, "_CV", CVthr, "_cognition.rds"))
)

fig_dir <- file.path(FigureFolder, Cogvar)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

correlation.df <- SCrankcorr(SC_Cog_results.df, "gam.smooth.t", 12, dsdata = TRUE)
correlation.df$sig <- (SC_Cog_results.df$anova.cov.p.fdr < 0.05)
Matrix.tmp <- matrix(NA, nrow = 12, ncol = 12)
linerange_frame <- data.frame(
  x = c(0.5, 12 + 0.5),
  ymin = rep(-12 - 0.5, times = 2),
  ymax = rep(-0.5, times = 2),
  y = c(-0.5, -12 - 0.5),
  xmin = rep(0.5, times = 2),
  xmax = rep(12 + 0.5, times = 2)
)
SC_Cog_results.tmp <- SC_Cog_results.df.whole
SC_Cog_results.tmp$SCrank <- correlation.df$SCrank

lwth <- min(SC_Cog_results.tmp$gam.smooth.t, na.rm = TRUE)
mytheme <- theme(
  axis.text = element_text(size = 23.2, color = "black"),
  axis.title = element_text(size = 23.2),
  aspect.ratio = 0.88,
  axis.line = element_line(linewidth = 0.6),
  axis.ticks = element_line(linewidth = 0.6),
  plot.title = element_text(size = 15, hjust = 0.5, vjust = 0),
  plot.subtitle = element_text(size = 15, hjust = 0.9, vjust = -6),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  plot.margin = margin(t = 10, r = 5, b = 5, l = 5, unit = "pt"),
  legend.position = "none"
)
width <- height <- 16.5

p_scatter <- ggplot(data = SC_Cog_results.tmp) +
  geom_point(aes(x = SCrank, y = gam.smooth.t, color = gam.smooth.t), size = 5.5) +
  geom_smooth(aes(x = SCrank, y = gam.smooth.t), method = "lm", color = "black", linewidth = 1.4) +
  scale_colour_distiller(type = "seq", palette = "RdBu", limits = c(lwth, -lwth), direction = -1) +
  labs(x = "S-A connectional axis rank", y = expression("Cognitive association (" * italic("T") * " value)")) +
  theme_classic() + mytheme

ggsave(file.path(fig_dir, "CorrTvalue_SCrankcorr_n12_siteall.tiff"), p_scatter, width = 17, height = 14, units = "cm", bg = "transparent")
ggsave(file.path(fig_dir, "CorrTvalue_SCrankcorr_n12_siteall.pdf"), p_scatter, dpi = 600, width = width, height = height, units = "cm", bg = "transparent")

Matrix.tmp.T <- Matrix.tmp
Matrix.tmp.T[lower.tri(Matrix.tmp.T, diag = TRUE)] <- SC_Cog_results.tmp$gam.smooth.t
Matrix.tmp.T[upper.tri(Matrix.tmp.T)] <- t(Matrix.tmp.T)[upper.tri(Matrix.tmp.T)]
colnames(Matrix.tmp.T) <- seq(1, 12)
rownames(Matrix.tmp.T) <- seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp.T)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df, id.vars = c("nodeid"))
matrixtmp.df.melt$variable <- as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid <- 0 - matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value <- as.numeric(matrixtmp.df.melt$value)

Matrix.tmp.sig <- Matrix.tmp
Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = TRUE)] <- correlation.df$sig
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
titlematrix <- Cogvar

p_matrix <- ggplot(data = matrixtmp.df.melt) +
  geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
  scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(lwth, -lwth), na.value = "grey") +
  scale_color_distiller(type = "seq", palette = "RdBu", limits = c(lwth, -lwth), na.value = "grey") +
  geom_text(data = matrixtmp.df.sig.melt, aes(x = variable, y = nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size = 6) +
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

ggsave(file.path(fig_dir, "CorrTvalue_Matrix_n12_siteall.tiff"), p_matrix, height = 15, width = 16, units = "cm", bg = "transparent")
ggsave(file.path(fig_dir, "CorrTvalue_Matrix_n12_siteall.pdf"), p_matrix, dpi = 600, height = 15, width = 16, units = "cm", bg = "transparent")

correlation2.df <- SCrankcorr(SC_Cog_results.df.whole, "gam.smooth.t_control_distance", 12, dsdata = TRUE)
lwth2 <- abs(min(SC_Cog_results.df.whole$gam.smooth.t_control_distance, na.rm = TRUE))
p_scatter2 <- ggplot(data = correlation2.df) +
  geom_point(aes(x = SCrank, y = gam.smooth.t_control_distance, color = gam.smooth.t_control_distance), size = 5) +
  geom_smooth(aes(x = SCrank, y = gam.smooth.t_control_distance), method = "lm", color = "black", linewidth = 1.4) +
  scale_colour_distiller(type = "seq", palette = "RdBu", limits = c(-lwth2, lwth2), direction = -1) +
  labs(x = "S-A connectional axis rank", y = expression("Cognitive association (" * italic("T") * " value)")) +
  theme_classic() +
  theme(
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

ggsave(file.path(fig_dir, "CorrTvalue_SCrankcorr_n12_siteall_control_distance.tiff"), p_scatter2, width = 13, height = 13, units = "cm", bg = "transparent")
ggsave(file.path(fig_dir, "CorrTvalue_SCrankcorr_n12_siteall_control_distance.pdf"), p_scatter2, dpi = 600, width = 15, height = 13.5, units = "cm", bg = "transparent")
