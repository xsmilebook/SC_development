#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(parallel)
})

rm(list = ls())

CVthr <- 75
Pvar <- "GENERAL"
Pvar_base <- "GENERAL_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "age_wp_bp_lmm_tvalue")
FigureFolder <- file.path(project_root, "outputs", "figures", "6th_pfactor", "abcd", "age_wp_bp_lmm_tvalue")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
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

scanid_to_eventname <- function(scanID) {
  sess <- sub("^.*_ses-", "", as.character(scanID))
  sess <- gsub("([a-z])([A-Z])", "\\1_\\2", sess)
  sess <- gsub("([A-Za-z])([0-9])", "\\1_\\2", sess)
  sess <- gsub("([0-9])([A-Za-z])", "\\1_\\2", sess)
  tolower(sess)
}

SCdata <- readRDS(input_rds)
if (!("eventname" %in% names(SCdata)) && ("scanID" %in% names(SCdata))) {
  SCdata$eventname <- scanid_to_eventname(SCdata$scanID)
}
if (!("eventname" %in% names(SCdata))) {
  stop("Missing eventname (required to construct baseline pfactor): input has no eventname/scanID")
}

needed <- c("subID", "age", "sex", "mean_fd", Pvar)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)

Pvardf <- SCdata %>%
  select(subID, eventname, all_of(Pvar)) %>%
  filter(!is.na(.data[[Pvar]])) %>%
  filter(grepl("base", eventname, ignore.case = TRUE)) %>%
  select(subID, !!Pvar_base := all_of(Pvar)) %>%
  distinct()
if (nrow(Pvardf) < 1) stop("No baseline p-factor rows found in: ", input_rds)

SCdata <- SCdata %>% left_join(Pvardf, by = "subID")
if (!Pvar_base %in% names(SCdata)) stop("Baseline p-factor join failed, missing: ", Pvar_base)
SCdata <- SCdata[!is.na(SCdata[[Pvar_base]]), , drop = FALSE]

sub_n <- table(SCdata$subID)
keep_sub <- names(sub_n[sub_n >= 2])
SCdata <- SCdata[SCdata$subID %in% keep_sub, , drop = FALSE]

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
sc_cols <- sc_cols[seq_len(78)]

plot_fit <- plotdata$fit
names(plot_fit) <- as.character(plotdata$SC_label)
missing_fit <- setdiff(sc_cols, names(plot_fit))
if (length(missing_fit) > 0) stop("plotdata missing fits for edges: ", paste(head(missing_fit, 10), collapse = ", "))

for (edge in sc_cols) {
  f0 <- as.numeric(plot_fit[[edge]])
  if (is.na(f0) || !is.finite(f0) || f0 == 0) stop("Invalid plotdata fit for edge: ", edge)
  SCdata[[edge]] <- as.numeric(SCdata[[edge]]) / f0
}

SCdata$age_bp <- ave(SCdata$age, SCdata$subID, FUN = mean)
SCdata$age_wp <- SCdata$age - SCdata$age_bp

base_by_sub <- tapply(SCdata[[Pvar_base]], SCdata$subID, function(x) x[which.max(!is.na(x))][1])
base_by_sub <- base_by_sub[!is.na(base_by_sub)]
q10 <- as.numeric(quantile(base_by_sub, 0.1, na.rm = TRUE))
q90 <- as.numeric(quantile(base_by_sub, 0.9, na.rm = TRUE))

age_min <- min(SCdata$age, na.rm = TRUE)
age_max <- max(SCdata$age, na.rm = TRUE)
age_seq <- seq(age_min, age_max, length.out = 100)

age_bp_mean <- mean(SCdata$age_bp, na.rm = TRUE)
mean_fd_mean <- mean(SCdata$mean_fd, na.rm = TRUE)
sex_tab <- table(SCdata$sex)
sex_ref <- names(sex_tab)[which.max(sex_tab)]
subid_ref <- as.character(SCdata$subID[1])

num_cores <- as.integer(Sys.getenv("LMM_CORES", unset = "16"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 16
num_cores <- min(num_cores, parallel::detectCores())

fit_edge <- function(i, data_all, edges, cov_name, q10, q90, age_seq, age_bp_mean, sex_ref, mean_fd_mean, subid_ref) {
  edge <- edges[[i]]
  df <- data_all[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", cov_name, edge)]
  names(df)[ncol(df)] <- "y"
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    return(list(
      row = data.frame(edge = edge, n_sub = nrow(df), stringsAsFactors = FALSE),
      pred_low = data.frame(),
      pred_high = data.frame()
    ))
  }

  df$cov <- df[[cov_name]]
  form <- stats::as.formula("y ~ age_wp * cov + age_bp + sex + mean_fd + (1 + age_wp || subID)")
  fit <- tryCatch(lme4::lmer(form, data = df, REML = FALSE), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(
      row = data.frame(edge = edge, n_sub = nrow(df), stringsAsFactors = FALSE),
      pred_low = data.frame(),
      pred_high = data.frame()
    ))
  }

  new_base <- data.frame(
    subID = subid_ref,
    age_wp = age_seq - age_bp_mean,
    age_bp = age_bp_mean,
    sex = sex_ref,
    mean_fd = mean_fd_mean,
    cov = q10,
    age = age_seq,
    stringsAsFactors = FALSE
  )
  if (is.factor(df$sex)) new_base$sex <- factor(sex_ref, levels = levels(df$sex))
  if (is.factor(df$subID)) new_base$subID <- factor(subid_ref, levels = levels(df$subID))

  pred_low_val <- tryCatch(
    stats::predict(fit, newdata = new_base, re.form = NA, allow.new.levels = TRUE),
    error = function(e) rep(NA_real_, length(age_seq))
  )
  pred_low <- data.frame(
    age = age_seq,
    .fitted = as.numeric(pred_low_val),
    label = "low",
    SC_label = edge,
    stringsAsFactors = FALSE
  )

  new_base$cov <- q90
  pred_high_val <- tryCatch(
    stats::predict(fit, newdata = new_base, re.form = NA, allow.new.levels = TRUE),
    error = function(e) rep(NA_real_, length(age_seq))
  )
  pred_high <- data.frame(
    age = age_seq,
    .fitted = as.numeric(pred_high_val),
    label = "high",
    SC_label = edge,
    stringsAsFactors = FALSE
  )

  list(
    row = data.frame(edge = edge, n_sub = nrow(df), stringsAsFactors = FALSE),
    pred_low = pred_low,
    pred_high = pred_high
  )
}

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_pfactor_tvalue_", Pvar_base, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)
pred_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_pfactor_pred_", Pvar_base, "_CV", CVthr, ".rds"))

if (!force && file.exists(pred_rds)) {
  message("[INFO] Found existing prediction results, loading (set FORCE=1 to recompute)")
  pred_list <- readRDS(pred_rds)
  res_df <- if (file.exists(out_rds)) readRDS(out_rds) else data.frame()
} else {
  message("[INFO] Fitting LMM (pfactor) with age_wp + age_bp and interaction")
  if (.Platform$OS.type == "windows") {
    message("[INFO] Windows parallel: ", num_cores, " workers")
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("SCdata", "sc_cols", "fit_edge", "q10", "q90", "age_seq", "age_bp_mean", "sex_ref", "mean_fd_mean", "subid_ref"),
      envir = environment()
    )
    res_list <- parallel::parLapply(
      cl, seq_along(sc_cols), fit_edge,
      data_all = SCdata,
      edges = sc_cols,
      cov_name = Pvar_base,
      q10 = q10,
      q90 = q90,
      age_seq = age_seq,
      age_bp_mean = age_bp_mean,
      sex_ref = sex_ref,
      mean_fd_mean = mean_fd_mean,
      subid_ref = subid_ref
    )
  } else {
    res_list <- lapply(
      seq_along(sc_cols), fit_edge,
      data_all = SCdata,
      edges = sc_cols,
      cov_name = Pvar_base,
      q10 = q10,
      q90 = q90,
      age_seq = age_seq,
      age_bp_mean = age_bp_mean,
      sex_ref = sex_ref,
      mean_fd_mean = mean_fd_mean,
      subid_ref = subid_ref
    )
  }

  res_rows <- lapply(res_list, `[[`, "row")
  pred_low_list <- lapply(res_list, `[[`, "pred_low")
  pred_high_list <- lapply(res_list, `[[`, "pred_high")

  res_df <- dplyr::bind_rows(res_rows)
  pred_list <- list(
    low = dplyr::bind_rows(pred_low_list),
    high = dplyr::bind_rows(pred_high_list)
  )

  saveRDS(res_df, out_rds)
  write.csv(res_df, out_csv, row.names = FALSE)
  saveRDS(pred_list, pred_rds)
}

plotdata.all <- dplyr::bind_rows(pred_list$low, pred_list$high)
plotdata.all <- merge(plotdata.all, SA12_10, by = "SC_label")

plotdf.decile.low <- plotdata.all %>%
  filter(label == "low") %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(.fitted, na.rm = TRUE), decile = mean(decile), .groups = "drop")
plotdf.decile.low$label <- "low"

plotdf.decile.high <- plotdata.all %>%
  filter(label == "high") %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(.fitted, na.rm = TRUE), decile = mean(decile), .groups = "drop")
plotdf.decile.high$label <- "high"

plotdf.decile <- rbind(plotdf.decile.low, plotdf.decile.high)
plotdf.decile$label <- factor(plotdf.decile$label, levels = c("low", "high"))

colorid <- rev(brewer.pal(10, "RdBu"))

interactionFolder <- file.path(FigureFolder, "Interaction")
dir.create(interactionFolder, showWarnings = FALSE, recursive = TRUE)

plot_list <- vector("list", length = 10)
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
    scale_y_continuous(breaks = c(0.9, 1.0, 1.1), limits = c(0.90, 1.1)) +
    scale_linetype_manual(values = c(low = "dashed", high = "solid")) +
    labs(x = NULL, y = "SC strength (ratio)") +
    mytheme

  out_base <- file.path(interactionFolder, paste0("developmentcurve_decile", i))
  ggsave(paste0(out_base, ".tiff"), Fig, width = 10, height = 10, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), Fig, width = 10, height = 10, units = "cm", bg = "transparent")
  plot_list[[i]] <- Fig
}

combined_fig <- patchwork::wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_base <- file.path(interactionFolder, "developmentcurve_decile_all_2x5")
ggsave(paste0(combined_base, ".tiff"), combined_fig, width = 50, height = 20, units = "cm", bg = "transparent")
ggsave(paste0(combined_base, ".pdf"), combined_fig, width = 50, height = 20, units = "cm", bg = "transparent")

message("[INFO] Done.")
