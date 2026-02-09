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
resultFolder <- file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", "abcd", "age_wp_bp_lmm_baselineage")
FigureFolder <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "abcd", "age_wp_bp_lmm_baselineage")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_tag <- "_baselineage_2tp"
out_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_SC_CV", CVthr, out_tag, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)

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

source(file.path(functionFolder, "SCrankcorr.R"))

scanid_to_eventname <- function(scanID) {
  sess <- sub("^.*_ses-", "", as.character(scanID))
  sess <- gsub("([a-z])([A-Z])", "\\1_\\2", sess)
  sess <- gsub("([A-Za-z])([0-9])", "\\1_\\2", sess)
  sess <- gsub("([0-9])([A-Za-z])", "\\1_\\2", sess)
  tolower(sess)
}

compute_baseline_age <- function(df, sub_col = "subID", age_col = "age", event_col = "eventname") {
  subs <- unique(as.character(df[[sub_col]]))
  out <- rep(NA_real_, length(subs))
  for (i in seq_along(subs)) {
    sid <- subs[[i]]
    idx <- which(as.character(df[[sub_col]]) == sid & !is.na(df[[age_col]]))
    if (length(idx) < 1) next
    age_i <- as.numeric(df[[age_col]][idx])
    event_i <- as.character(df[[event_col]][idx])
    idx_base <- which(grepl("base", event_i, ignore.case = TRUE))
    if (length(idx_base) > 0) {
      out[[i]] <- age_i[idx_base][which.min(age_i[idx_base])]
    } else {
      out[[i]] <- age_i[which.min(age_i)]
    }
  }
  data.frame(subID = subs, age_bp = out, stringsAsFactors = FALSE)
}

SCdata <- readRDS(input_rds)
if (!("eventname" %in% names(SCdata)) && ("scanID" %in% names(SCdata))) {
  SCdata$eventname <- scanid_to_eventname(SCdata$scanID)
}
if (!("eventname" %in% names(SCdata))) {
  stop("Missing eventname (required to construct baseline age): input has no eventname/scanID")
}
needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)

message(
  "[INFO] SCdata age range (years): ",
  round(min(SCdata$age, na.rm = TRUE), 3), "-", round(max(SCdata$age, na.rm = TRUE), 3)
)

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

base_age_df <- compute_baseline_age(SCdata)
if (any(is.na(base_age_df$age_bp))) {
  stop("Baseline age extraction failed for subjects: ", paste(base_age_df$subID[is.na(base_age_df$age_bp)], collapse = ", "))
}
SCdata <- SCdata %>% left_join(base_age_df, by = "subID")
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

plot_matrix <- function(mat, title, out_base, sig_mat = NULL) {
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

  sig_df <- data.frame()
  if (!is.null(sig_mat)) {
    sig_df <- as.data.frame(as.table(sig_mat))
    names(sig_df) <- c("nodeid", "variable", "sig")
    node_sig_raw <- sig_df$nodeid
    var_sig_raw <- sig_df$variable
    sig_df$nodeid <- suppressWarnings(as.numeric(as.character(node_sig_raw)))
    sig_df$variable <- suppressWarnings(as.numeric(as.character(var_sig_raw)))
    if (all(is.na(sig_df$nodeid))) sig_df$nodeid <- as.integer(node_sig_raw)
    if (all(is.na(sig_df$variable))) sig_df$variable <- as.integer(var_sig_raw)
    sig_df$nodeid <- -sig_df$nodeid
    sig_df <- sig_df[!is.na(sig_df$sig) & sig_df$sig, , drop = FALSE]
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

  if (nrow(sig_df) > 0) {
    p <- p + geom_text(data = sig_df, aes(x = variable, y = nodeid, label = "*"), vjust = 0.65, hjust = 0.5, size = 8)
  }

  ggsave(paste0(out_base, ".tiff"), p, height = 18, width = 20, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), p, height = 18, width = 20, units = "cm", bg = "transparent")
}

save_colorbar <- function(limthr, out_base) {
  if (!is.finite(limthr) || limthr == 0) limthr <- 1
  cb <- data.frame(
    x = seq(-limthr, limthr, length.out = 600),
    y = 1,
    z = seq(-limthr, limthr, length.out = 600)
  )

  p <- ggplot(cb, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "#F7F7F7",
      high = "#D7191C",
      midpoint = 0,
      limits = c(-limthr, limthr),
      guide = "none"
    ) +
    coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(0.1,0.1, 0.1, 0.1, "mm")
    )

  ggsave(paste0(out_base, ".tiff"), p, width = 12, height = 1.5, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".pdf"), p, width = 12, height = 1.5, units = "cm", bg = "transparent")
}

num_cores <- as.integer(Sys.getenv("LMM_CORES", unset = "16"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 16
num_cores <- min(num_cores, parallel::detectCores())

fit_edge <- function(i, data_all, edges) {
  edge <- edges[[i]]
  df <- data_all[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", edge)]
  names(df)[ncol(df)] <- "y"
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    return(data.frame(edge = edge, n_sub = nrow(df), beta_wp = NA_real_, beta_bp = NA_real_,
                      t_wp = NA_real_, t_bp = NA_real_, p_wp = NA_real_, p_bp = NA_real_))
  }

  full <- tryCatch(
    lme4::lmer(y ~ age_wp + age_bp + sex + mean_fd + (1 | subID), data = df, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(full)) {
    return(data.frame(edge = edge, n_sub = nrow(df), beta_wp = NA_real_, beta_bp = NA_real_,
                      t_wp = NA_real_, t_bp = NA_real_, p_wp = NA_real_, p_bp = NA_real_))
  }
  sm <- summary(full)
  beta_wp <- sm$coefficients["age_wp", "Estimate"]
  beta_bp <- sm$coefficients["age_bp", "Estimate"]
  t_wp <- sm$coefficients["age_wp", "t value"]
  t_bp <- sm$coefficients["age_bp", "t value"]
  red_wp <- tryCatch(
    lme4::lmer(y ~ age_bp + sex + mean_fd + (1 | subID), data = df, REML = FALSE),
    error = function(e) NULL
  )
  red_bp <- tryCatch(
    lme4::lmer(y ~ age_wp + sex + mean_fd + (1 | subID), data = df, REML = FALSE),
    error = function(e) NULL
  )
  get_lrt_p <- function(red, full_model) {
    if (is.null(red) || is.null(full_model)) return(NA_real_)
    atb <- tryCatch(stats::anova(red, full_model), error = function(e) NULL)
    if (is.null(atb) || nrow(atb) < 2 || !("Pr(>Chisq)" %in% names(atb))) return(NA_real_)
    as.numeric(atb$`Pr(>Chisq)`[2])
  }
  p_wp <- get_lrt_p(red_wp, full)
  p_bp <- get_lrt_p(red_bp, full)

  data.frame(
    edge = edge,
    n_sub = nrow(df),
    beta_wp = as.numeric(beta_wp),
    beta_bp = as.numeric(beta_bp),
    t_wp = as.numeric(t_wp),
    t_bp = as.numeric(t_bp),
    p_wp = as.numeric(p_wp),
    p_bp = as.numeric(p_bp),
    stringsAsFactors = FALSE
  )
}

need_refit <- force
if (!need_refit && file.exists(out_rds)) {
  message("[INFO] Found existing results, loading (set FORCE=1 to recompute): ", out_rds)
  res_df <- readRDS(out_rds)
  need_cols <- c("t_wp", "t_bp", "beta_wp", "beta_bp", "p_wp", "p_bp")
  missing_cols <- setdiff(need_cols, names(res_df))
  if (length(missing_cols) > 0) {
    message("[INFO] Existing results are missing new columns: ", paste(missing_cols, collapse = ", "), "; recomputing")
    need_refit <- TRUE
  }
}

if (need_refit || !file.exists(out_rds)) {
  message("[INFO] Fitting LMM (SC baseline-age decomposition) with age_wp + age_bp")
  if (.Platform$OS.type == "windows") {
    message("[INFO] Windows parallel: ", num_cores, " workers")
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("SCdata", "sc_cols", "fit_edge"),
      envir = environment()
    )
    res_list <- parallel::parLapply(cl, seq_along(sc_cols), fit_edge, data_all = SCdata, edges = sc_cols)
  } else {
    res_list <- lapply(seq_along(sc_cols), fit_edge, data_all = SCdata, edges = sc_cols)
  }

  res_df <- do.call(rbind, res_list)
}

res_df$p_wp_fdr <- stats::p.adjust(res_df$p_wp, method = "fdr")
res_df$p_bp_fdr <- stats::p.adjust(res_df$p_bp, method = "fdr")
saveRDS(res_df, out_rds)
write.csv(res_df, out_csv, row.names = FALSE)

message(
  "[INFO] FDR significant edges: age_wp=", sum(!is.na(res_df$p_wp_fdr) & res_df$p_wp_fdr < 0.05),
  ", age_bp=", sum(!is.na(res_df$p_bp_fdr) & res_df$p_bp_fdr < 0.05)
)

sig_wp_mat <- vec_to_mat(res_df$p_wp_fdr < 0.05, ds = 12)
sig_bp_mat <- vec_to_mat(res_df$p_bp_fdr < 0.05, ds = 12)

message("[INFO] age_wp t value matrix + S-A axis correlation")
mat_wp_t <- vec_to_mat(res_df$t_wp, ds = 12)
plot_matrix(
  mat_wp_t,
  "SC age_wp t value",
  file.path(FigureFolder, paste0("matrix_age_wp_tvalue_SC_CV", CVthr, out_tag)),
  sig_mat = sig_wp_mat
)
limthr_wp_t_mat <- 19.85
if (!is.finite(limthr_wp_t_mat) || limthr_wp_t_mat == 0) limthr_wp_t_mat <- 1
save_colorbar(
  limthr_wp_t_mat,
  file.path(FigureFolder, paste0("matrix_age_wp_tvalue_SC_CV", CVthr, out_tag, "_colorbar"))
)

SCrank.df.wp_t <- SCrankcorr(res_df, "t_wp", 12, dsdata = FALSE)
saveRDS(SCrank.df.wp_t, file.path(resultFolder, paste0("SCrankcorr_age_wp_tvalue_SC_CV", CVthr, out_tag, ".rds")))
message("[INFO] SCrankcorr (age_wp t value) r=", round(SCrank.df.wp_t$r.spearman, 3), " p=", signif(SCrank.df.wp_t$p.spearman, 3))

SCrank.data.wp_t <- SCrankcorr(res_df, "t_wp", 12, dsdata = TRUE)
limthr_wp_t <- max(abs(SCrank.data.wp_t$t_wp), na.rm = TRUE)
if (!is.finite(limthr_wp_t) || limthr_wp_t == 0) limthr_wp_t <- 1
scatterFig.wp_t <- ggplot(SCrank.data.wp_t) +
  geom_point(aes(x = SCrank, y = t_wp, color = t_wp), size = 5) +
  geom_smooth(aes(x = SCrank, y = t_wp), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr_wp_t, limthr_wp_t)) +
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
  labs(x = "S-A connectional axis rank", y = "Within-person age effect (T value)")

ggsave(file.path(FigureFolder, paste0("scatter_age_wp_tvalue_vs_SCrank_SC_CV", CVthr, out_tag, ".tiff")),
       scatterFig.wp_t, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_age_wp_tvalue_vs_SCrank_SC_CV", CVthr, out_tag, ".pdf")),
       scatterFig.wp_t, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] age_bp t value matrix + S-A axis correlation")
mat_bp_t <- vec_to_mat(res_df$t_bp, ds = 12)
plot_matrix(
  mat_bp_t,
  "SC age_bp t value",
  file.path(FigureFolder, paste0("matrix_age_bp_tvalue_SC_CV", CVthr, out_tag)),
  sig_mat = sig_bp_mat
)

SCrank.df.bp_t <- SCrankcorr(res_df, "t_bp", 12, dsdata = FALSE)
saveRDS(SCrank.df.bp_t, file.path(resultFolder, paste0("SCrankcorr_age_bp_tvalue_SC_CV", CVthr, out_tag, ".rds")))
message("[INFO] SCrankcorr (age_bp t value) r=", round(SCrank.df.bp_t$r.spearman, 3), " p=", signif(SCrank.df.bp_t$p.spearman, 3))

SCrank.data.bp_t <- SCrankcorr(res_df, "t_bp", 12, dsdata = TRUE)
limthr_bp_t <- max(abs(SCrank.data.bp_t$t_bp), na.rm = TRUE)
if (!is.finite(limthr_bp_t) || limthr_bp_t == 0) limthr_bp_t <- 1
scatterFig.bp_t <- ggplot(SCrank.data.bp_t) +
  geom_point(aes(x = SCrank, y = t_bp, color = t_bp), size = 5) +
  geom_smooth(aes(x = SCrank, y = t_bp), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr_bp_t, limthr_bp_t)) +
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
  labs(x = "S-A connectional axis rank", y = "age_bp t value")

ggsave(file.path(FigureFolder, paste0("scatter_age_bp_tvalue_vs_SCrank_SC_CV", CVthr, out_tag, ".tiff")),
       scatterFig.bp_t, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_age_bp_tvalue_vs_SCrank_SC_CV", CVthr, out_tag, ".pdf")),
       scatterFig.bp_t, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] age_wp fixed effect (beta_wp) matrix + S-A axis correlation")
mat_wp <- vec_to_mat(res_df$beta_wp, ds = 12)
plot_matrix(
  mat_wp,
  "SC age_wp fixed effect (beta)",
  file.path(FigureFolder, paste0("matrix_age_wp_beta_SC_CV", CVthr, out_tag)),
  sig_mat = sig_wp_mat
)

SCrank.df.wp <- SCrankcorr(res_df, "beta_wp", 12, dsdata = FALSE)
saveRDS(SCrank.df.wp, file.path(resultFolder, paste0("SCrankcorr_age_wp_beta_SC_CV", CVthr, out_tag, ".rds")))
message("[INFO] SCrankcorr (age_wp beta) r=", round(SCrank.df.wp$r.spearman, 3), " p=", signif(SCrank.df.wp$p.spearman, 3))

SCrank.data.wp <- SCrankcorr(res_df, "beta_wp", 12, dsdata = TRUE)
limthr_wp <- max(abs(SCrank.data.wp$beta_wp), na.rm = TRUE)
if (!is.finite(limthr_wp) || limthr_wp == 0) limthr_wp <- 1
scatterFig.wp <- ggplot(SCrank.data.wp) +
  geom_point(aes(x = SCrank, y = beta_wp, color = beta_wp), size = 5) +
  geom_smooth(aes(x = SCrank, y = beta_wp), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr_wp, limthr_wp)) +
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
  labs(x = "S-A connectional axis rank", y = "age_wp fixed effect (beta)")

ggsave(file.path(FigureFolder, paste0("scatter_age_wp_beta_vs_SCrank_SC_CV", CVthr, out_tag, ".tiff")),
       scatterFig.wp, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_age_wp_beta_vs_SCrank_SC_CV", CVthr, out_tag, ".pdf")),
       scatterFig.wp, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] age_bp fixed effect (beta_bp) matrix + S-A axis correlation")
mat_bp_beta <- vec_to_mat(res_df$beta_bp, ds = 12)
plot_matrix(
  mat_bp_beta,
  "SC age_bp fixed effect (beta)",
  file.path(FigureFolder, paste0("matrix_age_bp_beta_SC_CV", CVthr, out_tag)),
  sig_mat = sig_bp_mat
)

SCrank.df.bp_beta <- SCrankcorr(res_df, "beta_bp", 12, dsdata = FALSE)
saveRDS(SCrank.df.bp_beta, file.path(resultFolder, paste0("SCrankcorr_age_bp_beta_SC_CV", CVthr, out_tag, ".rds")))
message("[INFO] SCrankcorr (age_bp beta) r=", round(SCrank.df.bp_beta$r.spearman, 3), " p=", signif(SCrank.df.bp_beta$p.spearman, 3))

SCrank.data.bp_beta <- SCrankcorr(res_df, "beta_bp", 12, dsdata = TRUE)
limthr_bp_beta <- max(abs(SCrank.data.bp_beta$beta_bp), na.rm = TRUE)
if (!is.finite(limthr_bp_beta) || limthr_bp_beta == 0) limthr_bp_beta <- 1
scatterFig.bp_beta <- ggplot(SCrank.data.bp_beta) +
  geom_point(aes(x = SCrank, y = beta_bp, color = beta_bp), size = 5) +
  geom_smooth(aes(x = SCrank, y = beta_bp), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr_bp_beta, limthr_bp_beta)) +
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
  labs(x = "S-A connectional axis rank", y = "age_bp fixed effect (beta)")

ggsave(file.path(FigureFolder, paste0("scatter_age_bp_beta_vs_SCrank_SC_CV", CVthr, out_tag, ".tiff")),
       scatterFig.bp_beta, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_age_bp_beta_vs_SCrank_SC_CV", CVthr, out_tag, ".pdf")),
       scatterFig.bp_beta, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] Done.")
