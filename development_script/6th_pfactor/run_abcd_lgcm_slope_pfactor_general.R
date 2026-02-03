#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpattern)
})

rm(list = ls())

CVthr <- 75
Pvar <- "GENERAL"
Pvar_base <- "GENERAL_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "lgcm_slope")
FigureFolder <- file.path(project_root, "outputs", "figures", "6th_pfactor", "abcd", "lgcm_slope")
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

source(file.path(project_root, "gamfunction", "SCrankcorr.R"))

scanid_to_eventname <- function(scanID) {
  sess <- sub("^.*_ses-", "", as.character(scanID))
  sess <- gsub("([a-z])([A-Z])", "\\1_\\2", sess)
  sess <- gsub("([A-Za-z])([0-9])", "\\1_\\2", sess)
  sess <- gsub("([0-9])([A-Za-z])", "\\1_\\2", sess)
  tolower(sess)
}

age_to_years <- function(age_raw) {
  age_num <- as.numeric(age_raw)
  mx <- suppressWarnings(max(age_num, na.rm = TRUE))
  if (is.finite(mx) && mx > 24) {
    return(age_num / 12)
  }
  if (is.finite(mx) && mx > 0 && mx <= 2) {
    mx12 <- mx * 12
    mn <- suppressWarnings(min(age_num, na.rm = TRUE))
    mn12 <- mn * 12
    if (is.finite(mx12) && mx12 >= 6 && mx12 <= 30 && is.finite(mn12) && mn12 >= 4) {
      return(age_num * 12)
    }
  }
  age_num
}

SCdata <- readRDS(input_rds)
if (!("eventname" %in% names(SCdata)) && ("scanID" %in% names(SCdata))) {
  SCdata$eventname <- scanid_to_eventname(SCdata$scanID)
}
if (!("eventname" %in% names(SCdata))) {
  stop("Missing eventname (required to define baseline timepoint): input has no eventname/scanID")
}

SCdata$age <- age_to_years(SCdata$age)
message(
  "[INFO] SCdata age range (years): ",
  round(min(SCdata$age, na.rm = TRUE), 3), "–", round(max(SCdata$age, na.rm = TRUE), 3)
)

needed <- c("subID", "age", "sex", "mean_fd", Pvar)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
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

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
sc_cols <- sc_cols[seq_len(78)]

plot_fit <- plotdata$fit
names(plot_fit) <- as.character(plotdata$SC_label)
missing_fit <- setdiff(sc_cols, names(plot_fit))
if (length(missing_fit) > 0) stop("plotdata missing fits for edges: ", paste(head(missing_fit, 10), collapse = ", "))

# Normalize SC strength to ratio (divide by initial fit)
for (edge in sc_cols) {
  f0 <- as.numeric(plot_fit[[edge]])
  if (is.na(f0) || !is.finite(f0) || f0 == 0) stop("Invalid plotdata fit for edge: ", edge)
  SCdata[[edge]] <- as.numeric(SCdata[[edge]]) / f0
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

  limthr <- max(abs(df_melt$value), na.rm = TRUE)
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
    geom_segment(aes(x = 0.5, y = -0.5, xend = 12 + 0.5, yend = -12 - 0.5), color = "black", linewidth = 0.5) +
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

# Build per-subject t0/t1 rows
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
    sex = SCdata$sex[i0],
    mean_fd_t0 = SCdata$mean_fd[i0],
    mean_fd_t1 = SCdata$mean_fd[i1],
    pfactor_base = SCdata[[Pvar_base]][i0],
    delta_age = dt,
    stringsAsFactors = FALSE
  )
  for (edge in sc_cols) {
    v0 <- SCdata[[edge]][i0]
    v1 <- SCdata[[edge]][i1]
    out[[paste0(edge, "_t0")]] <- v0
    out[[paste0(edge, "_t1")]] <- v1
  }
  rows[[k]] <- out
}

dat_sub <- dplyr::bind_rows(rows)
if (nrow(dat_sub) < 10) stop("Too few subjects after two-timepoint selection: ", nrow(dat_sub))
dat_sub$sex <- as.factor(dat_sub$sex)

message("[INFO] Fitting per-edge LGCM-style slope LM (n_edges=78)")
res_rows <- vector("list", length(sc_cols))
pred_low <- numeric(length(sc_cols))
pred_high <- numeric(length(sc_cols))
for (i in seq_along(sc_cols)) {
  edge <- sc_cols[[i]]
  t0_col <- paste0(edge, "_t0")
  t1_col <- paste0(edge, "_t1")
  slope <- (dat_sub[[t1_col]] - dat_sub[[t0_col]]) / dat_sub$delta_age

  df <- data.frame(
    slope_per_year = slope,
    age_t0 = dat_sub$age_t0,
    SC_t0 = dat_sub[[t0_col]],
    pfactor_base = dat_sub$pfactor_base,
    sex = dat_sub$sex,
    mean_fd_t0 = dat_sub$mean_fd_t0,
    mean_fd_t1 = dat_sub$mean_fd_t1
  )
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    res_rows[[i]] <- data.frame(
      edge = edge,
      n_sub = nrow(df),
      beta_pfactor = NA_real_,
      t_pfactor = NA_real_,
      p_pfactor = NA_real_,
      stringsAsFactors = FALSE
    )
    pred_low[i] <- NA_real_
    pred_high[i] <- NA_real_
    next
  }

  lm_slope <- lm(slope_per_year ~ age_t0 + SC_t0 + pfactor_base + sex + mean_fd_t0 + mean_fd_t1, data = df)
  sm <- summary(lm_slope)
  beta_pf <- sm$coefficients["pfactor_base", "Estimate"]
  t_pf <- sm$coefficients["pfactor_base", "t value"]
  p_pf <- sm$coefficients["pfactor_base", "Pr(>|t|)"]

  res_rows[[i]] <- data.frame(
    edge = edge,
    n_sub = nrow(df),
    beta_pfactor = as.numeric(beta_pf),
    t_pfactor = as.numeric(t_pf),
    p_pfactor = as.numeric(p_pf),
    stringsAsFactors = FALSE
  )

  q10 <- quantile(df$pfactor_base, 0.1, na.rm = TRUE)
  q90 <- quantile(df$pfactor_base, 0.9, na.rm = TRUE)
  pred <- predict(lm_slope, newdata = df)
  pred_low[i] <- mean(pred[df$pfactor_base <= q10], na.rm = TRUE)
  pred_high[i] <- mean(pred[df$pfactor_base >= q90], na.rm = TRUE)
}

res_df <- do.call(rbind, res_rows)
res_df$p_pfactor_fdr <- p.adjust(res_df$p_pfactor, method = "fdr")
saveRDS(res_df, file.path(resultFolder, paste0("lgcm_slope_results_pfactor_", Pvar, "_CV", CVthr, ".rds")))
write.csv(res_df, file.path(resultFolder, paste0("lgcm_slope_results_pfactor_", Pvar, "_CV", CVthr, ".csv")), row.names = FALSE)

message("[INFO] Effect matrix + S-A axis correlation (beta_pfactor)")
beta_mat <- vec_to_mat(res_df$beta_pfactor, ds = 12)
sig_mat <- vec_to_mat(res_df$p_pfactor_fdr < 0.05, ds = 12)
plot_matrix_sig(
  beta_mat,
  sig_mat,
  "LGCM slope: pfactor effect (beta)",
  file.path(FigureFolder, paste0("matrix_lgcm_slope_beta_pfactor_", Pvar, "_CV", CVthr))
)

SCrank.df.beta <- SCrankcorr(res_df, "beta_pfactor", 12, dsdata = FALSE)
saveRDS(SCrank.df.beta, file.path(resultFolder, paste0("SCrankcorr_lgcm_slope_pfactor_", Pvar, "_CV", CVthr, "_beta.rds")))
message("[INFO] SCrankcorr (beta) r=", round(SCrank.df.beta$r.spearman, 3), " p=", signif(SCrank.df.beta$p.spearman, 3))

SCrank.data.beta <- SCrankcorr(res_df, "beta_pfactor", 12, dsdata = TRUE)
limthr.beta <- max(abs(SCrank.data.beta$beta_pfactor), na.rm = TRUE)
if (!is.finite(limthr.beta) || limthr.beta == 0) limthr.beta <- 1
scatterFig.beta <- ggplot(data = SCrank.data.beta) +
  geom_point(aes(x = SCrank, y = beta_pfactor, color = beta_pfactor), size = 5) +
  geom_smooth(aes(x = SCrank, y = beta_pfactor), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr.beta, limthr.beta)) +
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
  labs(x = "S-A connectional axis rank", y = "slope beta (pfactor)")

ggsave(
  file.path(FigureFolder, paste0("scatter_beta_vs_SCrank_lgcm_slope_pfactor_", Pvar, "_CV", CVthr, ".pdf")),
  scatterFig.beta, width = 13, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("scatter_beta_vs_SCrank_lgcm_slope_pfactor_", Pvar, "_CV", CVthr, ".tiff")),
  scatterFig.beta, width = 13, height = 12, units = "cm", bg = "transparent", dpi = 600
)

message("[INFO] Effect matrix + S-A axis correlation (t_pfactor)")
t_mat <- vec_to_mat(res_df$t_pfactor, ds = 12)
sig_mat_t <- vec_to_mat(res_df$p_pfactor_fdr < 0.05, ds = 12)
plot_matrix_sig(
  t_mat,
  sig_mat_t,
  "LGCM slope: pfactor effect (t value)",
  file.path(FigureFolder, paste0("matrix_lgcm_slope_tvalue_pfactor_", Pvar, "_CV", CVthr))
)

SCrank.df.t <- SCrankcorr(res_df, "t_pfactor", 12, dsdata = FALSE)
saveRDS(SCrank.df.t, file.path(resultFolder, paste0("SCrankcorr_lgcm_slope_pfactor_", Pvar, "_CV", CVthr, "_tvalue.rds")))
message("[INFO] SCrankcorr (t) r=", round(SCrank.df.t$r.spearman, 3), " p=", signif(SCrank.df.t$p.spearman, 3))

SCrank.data.t <- SCrankcorr(res_df, "t_pfactor", 12, dsdata = TRUE)
limthr.t <- max(abs(SCrank.data.t$t_pfactor), na.rm = TRUE)
if (!is.finite(limthr.t) || limthr.t == 0) limthr.t <- 1
scatterFig.t <- ggplot(data = SCrank.data.t) +
  geom_point(aes(x = SCrank, y = t_pfactor, color = t_pfactor), size = 5) +
  geom_smooth(aes(x = SCrank, y = t_pfactor), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr.t, limthr.t)) +
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
  labs(x = "S-A connectional axis rank", y = "slope t value (pfactor)")

ggsave(
  file.path(FigureFolder, paste0("scatter_tvalue_vs_SCrank_lgcm_slope_pfactor_", Pvar, "_CV", CVthr, ".pdf")),
  scatterFig.t, width = 13, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("scatter_tvalue_vs_SCrank_lgcm_slope_pfactor_", Pvar, "_CV", CVthr, ".tiff")),
  scatterFig.t, width = 13, height = 12, units = "cm", bg = "transparent", dpi = 600
)

# Decile summary for low/high predicted slopes


decile_df <- data.frame(
  SC_label = sc_cols,
  pred_low = pred_low,
  pred_high = pred_high,
  stringsAsFactors = FALSE
)
decile_df <- merge(decile_df, SA12_10, by.x = "SC_label", by.y = "SC_label", all.x = TRUE)
decile_summary <- decile_df %>%
  group_by(decile) %>%
  summarise(
    low_mean = mean(pred_low, na.rm = TRUE),
    high_mean = mean(pred_high, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  decile_summary,
  file.path(resultFolder, paste0("lgcm_slope_decile_low_high_pfactor_", Pvar, "_CV", CVthr, ".csv")),
  row.names = FALSE
)

plotdf <- data.frame(
  decile = decile_summary$decile,
  high = decile_summary$high_mean,
  low = decile_summary$low_mean
)
plotdf_long <- reshape(
  plotdf,
  varying = c("high", "low"),
  v.names = "mean",
  timevar = "group",
  times = c("high", "low"),
  direction = "long"
)
plotdf_long$decile <- as.integer(plotdf_long$decile)
plotdf_long$group <- factor(plotdf_long$group, levels = c("high", "low"))
colorid <- rev(brewer.pal(10, "RdBu"))
names(colorid) <- as.character(1:10)

bar_fig <- ggplot(plotdf_long, aes(x = factor(decile), y = mean, fill = factor(decile), group = group)) +
  geom_col_pattern(
    aes(pattern = group, alpha = group, linetype = group),
    color = "black",
    pattern_color = "black",
    pattern_fill = "transparent",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.6,
    width = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = colorid, guide = "none") +
  scale_pattern_manual(values = c(high = "none", low = "stripe"), guide = "none") +
  scale_alpha_manual(values = c(high = 1, low = 0.35), name = "Group", labels = c(high = "High level", low = "Low level")) +
  scale_linetype_manual(values = c(high = "solid", low = "dashed"), name = "Group", labels = c(high = "High level", low = "Low level")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  guides(
    alpha = guide_legend(override.aes = list(fill = "grey70", color = "black", pattern = c("none", "stripe"))),
    linetype = guide_legend(override.aes = list(fill = "grey70", color = "black"))
  ) +
  labs(x = "S-A decile", y = "Predicted slope per year")

ggsave(
  file.path(FigureFolder, paste0("bar_lgcm_slope_deciles_low_high_pfactor_", Pvar, "_CV", CVthr, ".tiff")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("bar_lgcm_slope_deciles_low_high_pfactor_", Pvar, "_CV", CVthr, ".pdf")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)

message("[INFO] Done.")
