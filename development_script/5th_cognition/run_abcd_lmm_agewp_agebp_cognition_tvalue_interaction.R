#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(parallel)
})

rm(list = ls())

CVthr <- 75
Cogvar <- "nihtbx_fluidcomp_uncorrected"
Cogvar_base <- "nihtbx_fluidcomp_uncorrected_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "age_wp_bp_lmm_tvalue")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "age_wp_bp_lmm_tvalue")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (age/sex/mean_fd variant; longitudinal)")
}

cog_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
)
if (!file.exists(cog_rds)) {
  stop("Missing cognition baseline input: ", cog_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (cognition variant)")
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

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "pb_lmm_anova.R"))

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
  stop("Missing eventname (required to construct baseline cognition): input has no eventname/scanID")
}

needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)

SCcog <- readRDS(cog_rds)
if (!("eventname" %in% names(SCcog)) && ("scanID" %in% names(SCcog))) {
  SCcog$eventname <- scanid_to_eventname(SCcog$scanID)
}
missing_cog <- setdiff(c("subID", "eventname", Cogvar), names(SCcog))
if (length(missing_cog) > 0) stop("Missing required columns in cognition SCdata: ", paste(missing_cog, collapse = ", "))

Cogdf <- SCcog %>%
  select(subID, eventname, all_of(Cogvar)) %>%
  filter(!is.na(.data[[Cogvar]])) %>%
  filter(grepl("base", eventname, ignore.case = TRUE)) %>%
  select(subID, !!Cogvar_base := all_of(Cogvar)) %>%
  distinct()
if (nrow(Cogdf) < 1) stop("No baseline cognition rows found in: ", cog_rds)

SCdata <- SCdata %>% left_join(Cogdf, by = "subID")
if (!Cogvar_base %in% names(SCdata)) stop("Baseline cognition join failed, missing: ", Cogvar_base)
SCdata <- SCdata[!is.na(SCdata[[Cogvar_base]]), , drop = FALSE]

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

base_by_sub <- tapply(SCdata[[Cogvar_base]], SCdata$subID, function(x) x[which.max(!is.na(x))][1])
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

  if (all(is.na(df_melt$value))) {
    limthr <- 1
  } else {
    limthr <- max(abs(df_melt$value), na.rm = TRUE)
  }
  if (!is.finite(limthr) || limthr == 0) {
    message("[WARN] Matrix values are all NA/0 for: ", title, "; set limthr=1 for plotting")
    limthr <- 1
  }

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
    scale_fill_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    scale_color_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
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

plot_matrix_sig <- function(mat, sig_mat, title, out_base) {
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

  sig_df <- as.data.frame(as.table(sig_mat))
  names(sig_df) <- c("nodeid", "variable", "sig")
  sig_df$nodeid <- suppressWarnings(as.numeric(as.character(sig_df$nodeid)))
  sig_df$variable <- suppressWarnings(as.numeric(as.character(sig_df$variable)))
  if (all(is.na(sig_df$nodeid))) sig_df$nodeid <- as.integer(sig_df$nodeid)
  if (all(is.na(sig_df$variable))) sig_df$variable <- as.integer(sig_df$variable)
  sig_df$nodeid <- -sig_df$nodeid
  sig_df <- sig_df[!is.na(sig_df$sig) & sig_df$sig, , drop = FALSE]

  if (all(is.na(df_melt$value))) {
    limthr <- 1
  } else {
    limthr <- max(abs(df_melt$value), na.rm = TRUE)
  }
  if (!is.finite(limthr) || limthr == 0) {
    message("[WARN] Matrix values are all NA/0 for: ", title, "; set limthr=1 for plotting")
    limthr <- 1
  }

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
    scale_fill_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    scale_color_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(-limthr, limthr)) +
    geom_text(data = sig_df, aes(x = variable, y = nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size = 8) +
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
num_cores <- as.integer(Sys.getenv("LMM_CORES", unset = "16"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 16
num_cores <- min(num_cores, parallel::detectCores())

pb_nsim <- as.integer(Sys.getenv("PB_NSIM", unset = "1000"))
if (is.na(pb_nsim) || pb_nsim < 1) pb_nsim <- 1000
pb_seed <- as.integer(Sys.getenv("PB_SEED", unset = "925"))
if (is.na(pb_seed) || pb_seed < 1) pb_seed <- 925

fit_edge <- function(i, data_all, edges, cov_name, q10, q90, age_seq, age_bp_mean, sex_ref, mean_fd_mean, subid_ref, pb_nsim, pb_seed) {
  edge <- edges[[i]]
  df <- data_all[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", cov_name, edge)]
  names(df)[ncol(df)] <- "y"
  names(df)[ncol(df) - 1] <- "cog_base"
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    return(list(
      row = data.frame(edge = edge, n_sub = nrow(df), t_cog = NA_real_,
                       bootstrap.P.cognition = NA_real_,
                       bootstrap_pvalue = NA_real_),
      pred_low = NULL,
      pred_high = NULL
    ))
  }

  df$subID <- as.factor(df$subID)
  df$sex <- as.factor(df$sex)

  t_cog <- NA_real_
  p_boot_main <- NA_real_
  p_boot_int <- NA_real_
  main_fit <- tryCatch(
    lme4::lmer(y ~ age_wp + age_bp + sex + mean_fd + (1 + age_wp || subID) + cog_base, data = df, REML = FALSE),
    error = function(e) NULL
  )
  null_main <- tryCatch(
    lme4::lmer(y ~ age_wp + age_bp + sex + mean_fd + (1 + age_wp || subID), data = df, REML = FALSE),
    error = function(e) NULL
  )
  if (!is.null(main_fit)) {
    sm <- summary(main_fit)
    if ("cog_base" %in% rownames(sm$coefficients)) {
      t_cog <- sm$coefficients["cog_base", "t value"]
    }
    if (!is.null(null_main)) {
      p_boot_main <- pb_lmm_anova(main_fit, null_main, nsim = pb_nsim, seed = pb_seed + i)
    }
  }

  int_fit <- tryCatch(
    lme4::lmer(y ~ age_wp * cog_base + age_bp + sex + mean_fd + (1 + age_wp || subID), data = df, REML = FALSE),
    error = function(e) NULL
  )
  if (!is.null(int_fit) && !is.null(main_fit)) {
    p_boot_int <- pb_lmm_anova(int_fit, main_fit, nsim = pb_nsim, seed = pb_seed + 1000L + i)
  }

  pred_low <- NULL
  pred_high <- NULL
  if (!is.null(int_fit)) {
    sub_levels <- levels(df$subID)
    sex_levels <- levels(df$sex)
    sub_ref <- if (length(sub_levels) > 0) sub_levels[[1]] else subid_ref
    sex_ref_use <- if (sex_ref %in% sex_levels) sex_ref else sex_levels[[1]]

    newdata_low <- data.frame(
      subID = factor(sub_ref, levels = sub_levels),
      age_bp = age_bp_mean,
      age_wp = age_seq - age_bp_mean,
      sex = factor(sex_ref_use, levels = sex_levels),
      mean_fd = mean_fd_mean,
      cog_base = q10
    )
    newdata_high <- newdata_low
    newdata_high$cog_base <- q90

    pred_low_vals <- tryCatch(predict(int_fit, newdata = newdata_low, re.form = NA), error = function(e) rep(NA_real_, length(age_seq)))
    pred_high_vals <- tryCatch(predict(int_fit, newdata = newdata_high, re.form = NA), error = function(e) rep(NA_real_, length(age_seq)))

    pred_low <- data.frame(age = age_seq, .fitted = as.numeric(pred_low_vals), label = "low", SC_label = edge)
    pred_high <- data.frame(age = age_seq, .fitted = as.numeric(pred_high_vals), label = "high", SC_label = edge)
  }

  list(
    row = data.frame(edge = edge, n_sub = nrow(df),
                     t_cog = as.numeric(t_cog),
                     bootstrap.P.cognition = as.numeric(p_boot_main),
                     bootstrap_pvalue = as.numeric(p_boot_int)),
    pred_low = pred_low,
    pred_high = pred_high
  )
}

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_cognition_tvalue_", Cogvar_base, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)
pred_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_cognition_pred_", Cogvar_base, "_CV", CVthr, ".rds"))

if (!force && file.exists(out_rds) && file.exists(pred_rds)) {
  message("[INFO] Found existing results, loading (set FORCE=1 to recompute)")
  res_df <- readRDS(out_rds)
  pred_list <- readRDS(pred_rds)
  need_cols <- c("t_cog", "bootstrap.P.cognition", "bootstrap_pvalue")
  missing_cols <- setdiff(need_cols, names(res_df))
  if (length(missing_cols) > 0) {
    stop(
      "Existing results are missing new columns: ", paste(missing_cols, collapse = ", "),
      "\nSet FORCE=1 to recompute: ", out_rds
    )
  }
} else {
  message("[INFO] Fitting LMM (cognition) with age_wp + age_bp and interaction")
  if (.Platform$OS.type == "windows") {
    message("[INFO] Windows parallel: ", num_cores, " workers")
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("SCdata", "sc_cols", "fit_edge", "q10", "q90", "age_seq", "age_bp_mean", "sex_ref", "mean_fd_mean", "subid_ref", "pb_lmm_anova"),
      envir = environment()
    )
    res_list <- parallel::parLapply(
      cl, seq_along(sc_cols), fit_edge,
      data_all = SCdata,
      edges = sc_cols,
      cov_name = Cogvar_base,
      q10 = q10,
      q90 = q90,
      age_seq = age_seq,
      age_bp_mean = age_bp_mean,
      sex_ref = sex_ref,
      mean_fd_mean = mean_fd_mean,
      subid_ref = subid_ref,
      pb_nsim = pb_nsim,
      pb_seed = pb_seed
    )
  } else {
    res_list <- lapply(
      seq_along(sc_cols), fit_edge,
      data_all = SCdata,
      edges = sc_cols,
      cov_name = Cogvar_base,
      q10 = q10,
      q90 = q90,
      age_seq = age_seq,
      age_bp_mean = age_bp_mean,
      sex_ref = sex_ref,
      mean_fd_mean = mean_fd_mean,
      subid_ref = subid_ref,
      pb_nsim = pb_nsim,
      pb_seed = pb_seed
    )
  }

  res_rows <- lapply(res_list, `[[`, "row")
  pred_low_list <- lapply(res_list, `[[`, "pred_low")
  pred_high_list <- lapply(res_list, `[[`, "pred_high")

  res_df <- dplyr::bind_rows(res_rows)
  res_df$bootstrap_pvalue.fdr <- p.adjust(res_df$bootstrap_pvalue, method = "fdr")
  res_df$bootstrap.P.cognition.fdr <- p.adjust(res_df$bootstrap.P.cognition, method = "fdr")
  pred_list <- list(
    low = dplyr::bind_rows(pred_low_list),
    high = dplyr::bind_rows(pred_high_list)
  )

  saveRDS(res_df, out_rds)
  write.csv(res_df, out_csv, row.names = FALSE)
  saveRDS(pred_list, pred_rds)
}

if (!("bootstrap_pvalue.fdr" %in% names(res_df))) {
  res_df$bootstrap_pvalue.fdr <- p.adjust(res_df$bootstrap_pvalue, method = "fdr")
}
if (!("bootstrap.P.cognition.fdr" %in% names(res_df))) {
  res_df$bootstrap.P.cognition.fdr <- p.adjust(res_df$bootstrap.P.cognition, method = "fdr")
}

message("[INFO] cognition t value matrix + S-A axis correlation")
mat_cog_t <- vec_to_mat(res_df$t_cog, ds = 12)
sig_cog <- vec_to_mat(res_df$bootstrap.P.cognition.fdr < 0.05, ds = 12)
plot_matrix_sig(
  mat_cog_t,
  sig_cog,
  "Cognition t value (age_wp/bp LMM)",
  file.path(FigureFolder, paste0("matrix_cognition_tvalue_", Cogvar_base, "_CV", CVthr))
)

SCrank.df.cog <- SCrankcorr(res_df, "t_cog", 12, dsdata = FALSE)
saveRDS(SCrank.df.cog, file.path(resultFolder, paste0("SCrankcorr_cognition_tvalue_", Cogvar_base, "_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr (cognition t value) r=", round(SCrank.df.cog$r.spearman, 3), " p=", signif(SCrank.df.cog$p.spearman, 3))

SCrank.data.cog <- SCrankcorr(res_df, "t_cog", 12, dsdata = TRUE)
limthr <- max(abs(SCrank.data.cog$t_cog), na.rm = TRUE)
if (!is.finite(limthr) || limthr == 0) limthr <- 1
scatterFig <- ggplot(SCrank.data.cog) +
  geom_point(aes(x = SCrank, y = t_cog, color = t_cog), size = 5) +
  geom_smooth(aes(x = SCrank, y = t_cog), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr, limthr)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 23, color = "black"),
    axis.title = element_text(size = 23),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  ) +
  labs(x = "S-A connectional axis rank", y = "Cognition t value")

ggsave(file.path(FigureFolder, paste0("scatter_cognition_tvalue_vs_SCrank_", Cogvar_base, "_CV", CVthr, ".tiff")),
       scatterFig, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_cognition_tvalue_vs_SCrank_", Cogvar_base, "_CV", CVthr, ".pdf")),
       scatterFig, width = 15, height = 15, units = "cm", bg = "transparent")

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

colorid <- rev(brewer.pal(10, "RdBu"))

interactionFolder <- file.path(FigureFolder, "Interaction")
dir.create(interactionFolder, showWarnings = FALSE, recursive = TRUE)

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
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = NULL, y = "SC strength (ratio)") +
    mytheme

  out_base <- file.path(interactionFolder, paste0("developmentcurve_decile", i))
  ggsave(paste0(out_base, ".tiff"), Fig, width = 10, height = 10, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), Fig, width = 10, height = 10, units = "cm", bg = "transparent")
}

message("[INFO] Done.")
