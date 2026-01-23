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
Cogvar <- "nihtbx_fluidcomp_agecorrected"
Cogvar_base <- "nihtbx_fluidcomp_agecorrected_base"
mode <- Sys.getenv("COG_ASSOC_MODE", unset = "original")
variant_tag <- Sys.getenv("COG_ASSOC_TAG", unset = "")
variant_suffix <- if (nzchar(variant_tag)) paste0("_", variant_tag) else if (mode == "meanfd_only") "_meanfd_only" else ""

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "comp_agecorrected")
intermediateFolder <- file.path(project_root, "outputs", "intermediate", "5th_cognition", "abcd", "comp_agecorrected")
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(intermediateFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch")
}

sa12_csv <- Sys.getenv(
  "ABCD_SA12_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "SA12_10.csv")
)
if (!file.exists(sa12_csv)) stop("Missing ABCD_SA12_CSV: ", sa12_csv)
SA12_10 <- read.csv(sa12_csv, stringsAsFactors = FALSE)

plotdatasum_rds <- Sys.getenv(
  "ABCD_PLOTDATASUM_RDS",
  unset = "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds"
)
if (!file.exists(plotdatasum_rds)) stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)
plotdata <- readRDS(plotdatasum_rds)

source(file.path(functionFolder, "gamminteraction.R"))

SCdata <- readRDS(input_rds)
SCdata$age <- as.numeric(SCdata$age) / 12

if (!all(c("subID", "eventname", "age", "mean_fd") %in% names(SCdata))) {
  stop("Missing required columns in SCdata: subID/age/mean_fd")
}
if (!("sex" %in% names(SCdata))) stop("Missing required column: sex")
SCdata$sex <- as.factor(SCdata$sex)

if (!Cogvar %in% names(SCdata)) {
  message(
    "[WARN] Missing phenotype column in longitudinal SCdata: ", Cogvar, "\n",
    "[WARN] SCdata input: ", input_rds, "\n",
    "[WARN] Will backfill baseline cognition from demopath/DemodfScreenFinal.csv and join by subID."
  )
  demopath_csv <- file.path(project_root, "demopath", "DemodfScreenFinal.csv")
  if (!file.exists(demopath_csv)) {
    stop("Missing demopath/DemodfScreenFinal.csv (git-ignored, required to backfill baseline cognition): ", demopath_csv)
  }
  Demodf <- read.csv(demopath_csv, stringsAsFactors = FALSE)
  missing_demo <- setdiff(c("subID", "eventname", Cogvar), names(Demodf))
  if (length(missing_demo) > 0) {
    stop("Missing required columns in demopath/DemodfScreenFinal.csv: ", paste(missing_demo, collapse = ", "))
  }
  Cogdf <- Demodf %>%
    select(subID, eventname, all_of(Cogvar)) %>%
    drop_na() %>%
    filter(str_detect(eventname, "base")) %>%
    select(subID, !!Cogvar_base := all_of(Cogvar)) %>%
    distinct()
  SCdata <- SCdata %>% left_join(Cogdf, by = "subID")
  if (!Cogvar_base %in% names(SCdata)) stop("Baseline cognition join failed, missing: ", Cogvar_base)
}

if (Cogvar %in% names(SCdata) && !Cogvar_base %in% names(SCdata)) {
  Cogdf <- SCdata %>%
    select(subID, eventname, all_of(Cogvar)) %>%
    drop_na() %>%
    filter(str_detect(eventname, "base")) %>%
    select(subID, !!Cogvar_base := all_of(Cogvar)) %>%
    distinct()
  SCdata <- SCdata %>% left_join(Cogdf, by = "subID")
  if (!Cogvar_base %in% names(SCdata)) stop("Baseline cognition join failed, missing: ", Cogvar_base)
}

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))

# Scale SC strength by initial fit (ratio) for visualization.
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
smooth_var <- "age"
int_var <- Cogvar_base
covariates <- "sex+mean_fd"
knots <- 3
set_fx <- TRUE
increments <- 1000
stats_only <- FALSE

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "40"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 40
num_cores <- min(num_cores, 40L)

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
trajectory_cache <- file.path(
  intermediateFolder,
  paste0("plotdata_high90_low10_", Cogvar_base, "_develop_CV", CVthr, variant_suffix, ".rds")
)

if (force || !file.exists(trajectory_cache)) {
  message("[INFO] Generating cognition-by-age interaction plot data (high90 vs low10), mc.cores=", num_cores, " suffix=", variant_suffix)
  resultsum <- parallel::mclapply(seq_len(78), function(i) {
    region <- sc_cols[[i]]

    tryCatch(
      {
        int_var_predict_percentile <- 0.1
        out_low <- gamm.smooth.predict.covariateinteraction(
          region, dataname, smooth_var, int_var, int_var_predict_percentile,
          covariates, knots, set_fx, increments, stats_only = stats_only
        )[[2]]
        out_low$SC_label <- region
        out_low$cognitionlevel <- "low"

        int_var_predict_percentile <- 0.9
        out_high <- gamm.smooth.predict.covariateinteraction(
          region, dataname, smooth_var, int_var, int_var_predict_percentile,
          covariates, knots, set_fx, increments, stats_only = stats_only
        )[[2]]
        out_high$SC_label <- region
        out_high$cognitionlevel <- "high"

        list(ok = TRUE, df = rbind(out_low, out_high), err = NA_character_, region = region)
      },
      error = function(e) list(ok = FALSE, df = NULL, err = conditionMessage(e), region = region)
    )
  }, mc.cores = num_cores)

  ok_mask <- vapply(resultsum, function(z) isTRUE(z$ok), logical(1))
  if (!any(ok_mask)) {
    errs <- unique(vapply(resultsum, function(z) paste0(z$region, ": ", z$err), character(1)))
    stop("All edges failed in interaction fitting. First errors:\n", paste(head(errs, 10), collapse = "\n"))
  }

  saveRDS(resultsum, trajectory_cache)
}

resultsum <- readRDS(trajectory_cache)
plotdf <- do.call(rbind, lapply(resultsum, `[[`, "df"))
plotdf <- merge(plotdf, SA12_10, by.x = "SC_label", by.y = "SC_label")

plotdf.decile.low <- plotdf %>%
  filter(cognitionlevel == "low") %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(.fitted), decile = mean(decile), .groups = "drop")
plotdf.decile.low$cognitionlevel <- "low"

plotdf.decile.high <- plotdf %>%
  filter(cognitionlevel == "high") %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(.fitted), decile = mean(decile), .groups = "drop")
plotdf.decile.high$cognitionlevel <- "high"

plotdf.decile <- rbind(plotdf.decile.low, plotdf.decile.high)

interaction_dir <- file.path(FigureFolder, Cogvar, "Interaction")
dir.create(interaction_dir, showWarnings = FALSE, recursive = TRUE)

colorid <- rev(brewer.pal(10, "RdBu"))
for (i in 1:10) {
  plotdf.tmp <- plotdf.decile[plotdf.decile$decile == i, , drop = FALSE]
  colorindex <- colorid[i]
  age_label_mult <- if (max(plotdf.tmp$age, na.rm = TRUE) <= 2) 10 else 1

  if (i == 1) {
    mytheme <- theme(
      axis.text = element_text(size = 21, color = "black"),
      axis.title = element_text(size = 21),
      aspect.ratio = 1.2,
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
    geom_line(aes(x = age, y = fit.avg, group = cognitionlevel, linetype = cognitionlevel), linewidth = 1.2, color = colorindex) +
    scale_x_continuous(labels = function(x) x * age_label_mult) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_y_continuous(breaks = c(0.9, 1.0, 1.1), limits = c(0.89, 1.1)) +
    labs(x = NULL, y = "SC strength (ratio)") +
    mytheme

  out_base <- file.path(interaction_dir, paste0("developmentcurve_decile", i, variant_suffix))
  ggsave(paste0(out_base, ".tiff"), Fig, width = 10, height = 10, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), Fig, dpi = 600, width = 10, height = 10, units = "cm", bg = "transparent")
}
