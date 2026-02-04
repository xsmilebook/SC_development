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

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", "abcd", "age_wp_bp_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "abcd", "age_wp_bp_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

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

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "pb_lmm_anova.R"))

SCdata <- readRDS(input_rds)
needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)

message(
  "[INFO] SCdata age range (years): ",
  round(min(SCdata$age, na.rm = TRUE), 3), "–", round(max(SCdata$age, na.rm = TRUE), 3)
)

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

pb_nsim <- as.integer(Sys.getenv("PB_NSIM", unset = "1000"))
if (is.na(pb_nsim) || pb_nsim < 1) pb_nsim <- 1000
pb_seed <- as.integer(Sys.getenv("PB_SEED", unset = "925"))
if (is.na(pb_seed) || pb_seed < 1) pb_seed <- 925
num_cores <- as.integer(Sys.getenv("LMM_CORES", unset = "16"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 16
num_cores <- min(num_cores, parallel::detectCores())

fit_edge <- function(i, data_all, edges, pb_nsim, pb_seed) {
  edge <- edges[[i]]
  df <- data_all[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", edge)]
  names(df)[ncol(df)] <- "y"
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    return(data.frame(edge = edge, n_sub = nrow(df), beta_wp = NA_real_, beta_bp = NA_real_,
                      t_wp = NA_real_, t_bp = NA_real_, partial_r2_bp = NA_real_, p_bp = NA_real_))
  }

  full <- lme4::lmer(y ~ age_wp + age_bp + sex + mean_fd + (1 + age_wp || subID), data = df, REML = FALSE)
  null <- lme4::lmer(y ~ age_wp + sex + mean_fd + (1 + age_wp || subID), data = df, REML = FALSE)
  sm <- summary(full)
  beta_wp <- sm$coefficients["age_wp", "Estimate"]
  beta_bp <- sm$coefficients["age_bp", "Estimate"]
  t_wp <- sm$coefficients["age_wp", "t value"]
  t_bp <- sm$coefficients["age_bp", "t value"]

  sse_full <- sum(residuals(full)^2)
  sse_null <- sum(residuals(null)^2)
  partial_r2 <- if (is.finite(sse_null) && sse_null > 0) (sse_null - sse_full) / sse_null else NA_real_
  p_bp <- pb_lmm_anova(full, null, nsim = pb_nsim, seed = pb_seed + i)

  data.frame(
    edge = edge,
    n_sub = nrow(df),
    beta_wp = as.numeric(beta_wp),
    beta_bp = as.numeric(beta_bp),
    t_wp = as.numeric(t_wp),
    t_bp = as.numeric(t_bp),
    partial_r2_bp = as.numeric(partial_r2),
    p_bp = as.numeric(p_bp),
    stringsAsFactors = FALSE
  )
}

message("[INFO] Fitting LMM (SC) with age_wp + age_bp")
if (.Platform$OS.type == "windows") {
  message("[INFO] Windows parallel: ", num_cores, " workers")
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(
    cl,
    varlist = c("SCdata", "sc_cols", "pb_nsim", "pb_seed", "fit_edge", "pb_lmm_anova"),
    envir = environment()
  )
  res_list <- parallel::parLapply(cl, seq_along(sc_cols), fit_edge, data_all = SCdata, edges = sc_cols, pb_nsim = pb_nsim, pb_seed = pb_seed)
} else {
  res_list <- lapply(seq_along(sc_cols), fit_edge, data_all = SCdata, edges = sc_cols, pb_nsim = pb_nsim, pb_seed = pb_seed)
}

res_df <- do.call(rbind, res_list)
res_df$p_bp_fdr <- p.adjust(res_df$p_bp, method = "fdr")
out_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_SC_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)
saveRDS(res_df, out_rds)
write.csv(res_df, out_csv, row.names = FALSE)

message("[INFO] age_bp partial R2 matrix + S-A axis correlation")
mat_bp <- vec_to_mat(res_df$partial_r2_bp, ds = 12)
sig_bp <- vec_to_mat(res_df$p_bp_fdr < 0.05, ds = 12)
plot_matrix_sig(
  mat_bp,
  sig_bp,
  "SC age_bp partial R2",
  file.path(FigureFolder, paste0("matrix_age_bp_partialR2_SC_CV", CVthr))
)

SCrank.df.bp <- SCrankcorr(res_df, "partial_r2_bp", 12, dsdata = FALSE)
saveRDS(SCrank.df.bp, file.path(resultFolder, paste0("SCrankcorr_age_bp_partialR2_SC_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr (age_bp partial R2) r=", round(SCrank.df.bp$r.spearman, 3), " p=", signif(SCrank.df.bp$p.spearman, 3))

SCrank.data.bp <- SCrankcorr(res_df, "partial_r2_bp", 12, dsdata = TRUE)
limthr <- max(abs(SCrank.data.bp$partial_r2_bp), na.rm = TRUE)
if (!is.finite(limthr) || limthr == 0) limthr <- 1
scatterFig <- ggplot(SCrank.data.bp) +
  geom_point(aes(x = SCrank, y = partial_r2_bp, color = partial_r2_bp), size = 5) +
  geom_smooth(aes(x = SCrank, y = partial_r2_bp), method = "lm", color = "black", linewidth = 1.4) +
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
  labs(x = "S-A connectional axis rank", y = "age_bp partial R2")
ggsave(file.path(FigureFolder, paste0("scatter_age_bp_partialR2_vs_SCrank_SC_CV", CVthr, ".tiff")),
       scatterFig, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_age_bp_partialR2_vs_SCrank_SC_CV", CVthr, ".pdf")),
       scatterFig, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] Done.")
