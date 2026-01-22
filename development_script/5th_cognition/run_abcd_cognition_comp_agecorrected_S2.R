#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(psych)
  library(reshape)
  library(RColorBrewer)
  library(tidyverse)
  library(corrplot)
  library(ggcorrplot)
  library(geomtextpath)
})

rm(list = ls())

CVthr <- 75
Cogvar <- "nihtbx_totalcomp_agecorrected"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "comp_agecorrected")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "comp_agecorrected")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "gamcog.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_comp_agecorrected_baseline.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam_comp_agecorrected_baseline.sbatch")
}

SCdata <- readRDS(input_rds)
SCdata$age <- as.numeric(SCdata$age) / 12

in_rds <- file.path(resultFolder, paste0("SC_Cog_results_", Cogvar, "_CV", CVthr, "_comp_agecorrected.rds"))
if (!file.exists(in_rds)) {
  stop("Missing S1 output: ", in_rds, "\nRun first: Rscript development_script/5th_cognition/run_abcd_cognition_comp_agecorrected_S1.R")
}
SC_Cog_results.df <- readRDS(in_rds)

SC_Cog_results.df$SClabelorder <- c(1:78)
SCrank_df <- SCrankcorr(SC_Cog_results.df, "gam.cog.t", 12, dsdata = TRUE)
SC_Cog_results.df$SCrank <- SCrank_df$SCrank

## 1) "component" plot (for this analysis, cognition measure is the age-corrected total score)
comp_plot <- data.frame(category = "NIH Toolbox total cognition (age-corrected)", value = 1)
p_comp <- ggplot(comp_plot, aes(x = category, y = value)) +
  geom_col(fill = brewer.pal(3, "Set2")[[1]]) +
  coord_flip() +
  labs(x = NULL, y = NULL, title = "Cognition measure used in this analysis") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
ggsave(file.path(FigureFolder, "NTB_totalcognition_agecorrected_component.tiff"), p_comp, width = 16, height = 6, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, "NTB_totalcognition_agecorrected_component.pdf"), p_comp, dpi = 600, width = 16, height = 6, units = "cm", bg = "transparent")

## 2) Scatter plots (3 example connections at ~10/50/90% SCrank)
nonna_index <- which(!is.na(SCdata[, Cogvar]))
SCdata.cog <- SCdata[nonna_index, , drop = FALSE]
if ("eventname" %in% names(SCdata.cog)) {
  SCdata.cog <- SCdata.cog[SCdata.cog$eventname == "baseline_year_1_arm_1", , drop = FALSE]
}

SCdata.cog[, str_detect(names(SCdata.cog), "SC.")] <- lapply(SCdata.cog[, str_detect(names(SCdata.cog), "SC.")], as.numeric)
dataname <- "SCdata.cog"
smooth_var <- "age"
covariates <- "sex+mean_fd"
knots <- 3
corrmethod <- "pearson"

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "72"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 72

make_cluster_fallback <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 1) n <- 1L
  stepdowns <- c(60L, 50L, 40L, 30L, 20L)
  stepdowns <- stepdowns[stepdowns < n]
  tries <- unique(pmax(1L, c(n, stepdowns, 16L, 8L, 4L, 2L, 1L)))
  for (k in tries) {
    cl <- try(parallel::makeCluster(k), silent = TRUE)
    if (!inherits(cl, "try-error")) {
      message("[INFO] Using PSOCK workers: ", k)
      return(cl)
    }
    message("[WARN] makeCluster(", k, ") failed; retrying with fewer workers.")
  }
  stop("Failed to create any PSOCK cluster (process limit or memory pressure).")
}

pick_nearest <- function(target) {
  idx <- which.min(abs(SC_Cog_results.df$SCrank - target))
  idx[[1]]
}
targets <- as.numeric(stats::quantile(SC_Cog_results.df$SCrank, probs = c(0.1, 0.5, 0.9)))
edge_idx <- unique(vapply(targets, pick_nearest, integer(1)))
edge_idx <- edge_idx[seq_len(min(length(edge_idx), 3))]

cl <- make_cluster_fallback(num_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(mgcv)
    library(parallel)
    library(psych)
    library(gratia)
    library(tidyverse)
    library(reshape)
  })
  source(file.path(functionFolder, "gaminteraction.R"))
  source(file.path(functionFolder, "gamcog.R"))
}))

plotdata_list <- parLapply(cl, edge_idx, function(x) {
  SClabel <- grep("SC.", names(SCdata.cog), value = TRUE)[x]
  gamresult.df <- gam.fit.cognition(
    SClabel,
    dataname,
    Cogvar,
    smooth_var,
    covariates,
    knots,
    corrmethod,
    set_fx = TRUE,
    stats_only = FALSE
  )
  out <- as.data.frame(gamresult.df[[2]])
  out$SC_label <- SClabel
  out
})

sc_fig_dir <- file.path(FigureFolder, Cogvar)
dir.create(sc_fig_dir, showWarnings = FALSE, recursive = TRUE)

for (k in seq_along(edge_idx)) {
  N <- edge_idx[[k]]
  plotdata_N <- plotdata_list[[k]]
  SCrank <- SC_Cog_results.df$SCrank[N]
  message("Scatterplot for connection ", N, " with SCrank ", SCrank)

  p_sc <- ggplot(data = plotdata_N) +
    geom_point(aes(x = SCres, y = cogres), color = "grey", size = 0.8) +
    geom_smooth(aes(x = SCres, y = cogres), linewidth = 1.4, method = "lm", color = "black") +
    labs(x = "SC strength (residual)", y = "Cognition (residual)") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 22, color = "black"),
      axis.title = element_text(size = 22),
      aspect.ratio = 0.9,
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none"
    )

  ggsave(file.path(sc_fig_dir, paste0("SC", N, "_scatterplot.tiff")), p_sc, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
  ggsave(file.path(sc_fig_dir, paste0("SC", N, "_scatterplot.pdf")), p_sc, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
}

## 3) Distribution of significant correlations
SC_Cog_results.df.sig <- SC_Cog_results.df %>% filter(anova.cov.p.fdr < 0.05)
p_hist <- ggplot(data = SC_Cog_results.df.sig, aes(correstimate, y = ..count..)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "gray", position = position_dodge()) +
  labs(x = expression("Correlation (" * italic("r") * ")"), y = "Frequency", title = "Cognitive association") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 0),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )
ggsave(file.path(sc_fig_dir, "SigCorrestimateDistribution.tiff"), p_hist, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
ggsave(file.path(sc_fig_dir, "SigCorrestimateDistribution.pdf"), p_hist, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
