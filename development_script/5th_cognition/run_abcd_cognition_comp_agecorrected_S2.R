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
Cogvar <- "nihtbx_fluidcomp_agecorrected"
variant_tag <- Sys.getenv("COG_ASSOC_TAG", unset = "")
variant_suffix <- if (nzchar(variant_tag)) paste0("_", variant_tag) else ""

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

in_rds <- file.path(
  resultFolder,
  paste0("SC_Cog_results_", Cogvar, "_CV", CVthr, "_comp_agecorrected", variant_suffix, ".rds")
)
if (!file.exists(in_rds)) {
  stop("Missing S1 output: ", in_rds, "\nRun first: Rscript development_script/5th_cognition/run_abcd_cognition_comp_agecorrected_S1.R")
}
SC_Cog_results.df <- readRDS(in_rds)

SC_Cog_results.df$SClabelorder <- c(1:78)
SCrank_df <- SCrankcorr(SC_Cog_results.df, "gam.smooth.t", 12, dsdata = TRUE)
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
ggsave(
  file.path(FigureFolder, paste0("NTB_totalcognition_agecorrected_component", variant_suffix, ".tiff")),
  p_comp,
  width = 16,
  height = 6,
  units = "cm",
  bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("NTB_totalcognition_agecorrected_component", variant_suffix, ".pdf")),
  p_comp,
  dpi = 600,
  width = 16,
  height = 6,
  units = "cm",
  bg = "transparent"
)

## 2) Scatter plots (3 example connections at ~10/50/90% SCrank)
nonna_index <- which(!is.na(SCdata[, Cogvar]))
SCdata.cog <- SCdata[nonna_index, , drop = FALSE]
if ("eventname" %in% names(SCdata.cog)) {
  SCdata.cog <- SCdata.cog[SCdata.cog$eventname == "baseline_year_1_arm_1", , drop = FALSE]
}

SCdata.cog[, str_detect(names(SCdata.cog), "SC.")] <- lapply(SCdata.cog[, str_detect(names(SCdata.cog), "SC.")], as.numeric)
dataname <- "SCdata.cog"
smooth_var <- ""
covariates <- "mean_fd"
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

pick_nearest <- function(target) {
  idx <- which.min(abs(SC_Cog_results.df$SCrank - target))
  idx[[1]]
}
targets <- as.numeric(stats::quantile(SC_Cog_results.df$SCrank, probs = c(0.1, 0.5, 0.9)))
edge_idx <- unique(vapply(targets, pick_nearest, integer(1)))
edge_idx <- edge_idx[seq_len(min(length(edge_idx), 3))]

sc_labels <- grep("SC\\.", names(SCdata.cog), value = TRUE)
get_plotdata_for_edge <- function(edge_index) {
  SClabel <- sc_labels[[edge_index]]
  tryCatch(
    {
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
      list(ok = TRUE, idx = edge_index, label = SClabel, df = out, err = NA_character_)
    },
    error = function(e) {
      list(ok = FALSE, idx = edge_index, label = SClabel, df = NULL, err = conditionMessage(e))
    }
  )
}

pick_working_edge <- function(target_rank, max_try = 20) {
  ord <- order(abs(SC_Cog_results.df$SCrank - target_rank))
  ord <- ord[seq_len(min(length(ord), max_try))]
  for (idx in ord) {
    res <- get_plotdata_for_edge(idx)
    if (isTRUE(res$ok)) return(res)
    message("[WARN] Plotdata failed for edge ", idx, " (", res$label, "): ", res$err)
  }
  list(ok = FALSE, idx = NA_integer_, label = NA_character_, df = NULL, err = "no working edge found")
}

sc_fig_dir <- file.path(FigureFolder, Cogvar)
dir.create(sc_fig_dir, showWarnings = FALSE, recursive = TRUE)

for (k in seq_along(targets)) {
  target_rank <- targets[[k]]
  res <- pick_working_edge(target_rank, max_try = 30)
  if (!isTRUE(res$ok)) {
    message("[WARN] No working edge found for quantile ", k, "; saving placeholder plot")
    out_base <- file.path(sc_fig_dir, paste0("SC_missing_q", k, "_scatterplot", variant_suffix))
    p_empty <- ggplot() +
      geom_text(aes(x = 0, y = 0), label = paste0("No valid edge for q", k, " (", variant_tag, ")")) +
      theme_void()
    ggsave(paste0(out_base, ".tiff"), p_empty, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
    ggsave(paste0(out_base, ".pdf"), p_empty, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
    next
  }

  plotdata_N <- res$df
  N <- res$idx
  SCrank <- SC_Cog_results.df$SCrank[N]
  message("Scatterplot for connection ", N, " with SCrank ", SCrank, " (target=", target_rank, ")")
  out_base <- file.path(sc_fig_dir, paste0("SC", N, "_scatterplot", variant_suffix))

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

  ggsave(paste0(out_base, ".tiff"), p_sc, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), p_sc, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
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
ggsave(
  file.path(sc_fig_dir, paste0("SigCorrestimateDistribution", variant_suffix, ".tiff")),
  p_hist,
  width = 13.5,
  height = 13.5,
  units = "cm",
  bg = "transparent"
)
ggsave(
  file.path(sc_fig_dir, paste0("SigCorrestimateDistribution", variant_suffix, ".pdf")),
  p_hist,
  dpi = 600,
  width = 13.5,
  height = 13.5,
  units = "cm",
  bg = "transparent"
)
