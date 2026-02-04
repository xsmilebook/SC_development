#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpattern)
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
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "age_wp_bp_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "age_wp_bp_lmm")
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

fit_edge <- function(i, data_all, edges, base_by_sub, q10, q90, pb_nsim, pb_seed) {
  edge <- edges[[i]]
  df <- data_all[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", edge)]
  names(df)[ncol(df)] <- "y"
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    return(list(
      row = data.frame(edge = edge, n_sub = nrow(df), beta_wp = NA_real_, beta_bp = NA_real_,
                       t_wp = NA_real_, t_bp = NA_real_, partial_r2_bp = NA_real_, p_bp = NA_real_,
                       personal_t_low_high = NA_real_, personal_p_low_high = NA_real_,
                       personal_cor = NA_real_, personal_cor_p = NA_real_,
                       personal_low_mean = NA_real_, personal_high_mean = NA_real_),
      edge = edge
    ))
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

  re <- tryCatch(ranef(full)[["subID"]], error = function(e) NULL)
  if (is.null(re) || !"age_wp" %in% names(re)) {
    return(list(
      row = data.frame(edge = edge, n_sub = nrow(df), beta_wp = as.numeric(beta_wp), beta_bp = as.numeric(beta_bp),
                       t_wp = as.numeric(t_wp), t_bp = as.numeric(t_bp), partial_r2_bp = as.numeric(partial_r2),
                       p_bp = as.numeric(p_bp), personal_t_low_high = NA_real_, personal_p_low_high = NA_real_,
                       personal_cor = NA_real_, personal_cor_p = NA_real_,
                       personal_low_mean = NA_real_, personal_high_mean = NA_real_),
      edge = edge
    ))
  }

  personal <- as.numeric(fixef(full)["age_wp"]) + as.numeric(re[,"age_wp"])
  slope_df <- data.frame(subID = rownames(re), personal = personal, stringsAsFactors = FALSE)
  base_df <- data.frame(subID = names(base_by_sub), base = as.numeric(base_by_sub), stringsAsFactors = FALSE)
  slope_df <- merge(slope_df, base_df, by = "subID")

  low_idx <- slope_df$base <= q10
  high_idx <- slope_df$base >= q90
  low_mean <- mean(slope_df$personal[low_idx], na.rm = TRUE)
  high_mean <- mean(slope_df$personal[high_idx], na.rm = TRUE)

  t_out <- tryCatch(stats::t.test(slope_df$personal[low_idx], slope_df$personal[high_idx]), error = function(e) NULL)
  t_val <- if (is.null(t_out) || !is.finite(t_out$statistic)) 0 else as.numeric(t_out$statistic)
  p_val <- if (is.null(t_out) || !is.finite(t_out$p.value)) 1 else as.numeric(t_out$p.value)

  cor_out <- tryCatch(stats::cor.test(slope_df$personal, slope_df$base), error = function(e) NULL)
  cor_val <- if (is.null(cor_out) || !is.finite(cor_out$estimate)) NA_real_ else as.numeric(cor_out$estimate)
  cor_p <- if (is.null(cor_out) || !is.finite(cor_out$p.value)) NA_real_ else as.numeric(cor_out$p.value)

  list(
    row = data.frame(
      edge = edge,
      n_sub = nrow(df),
      beta_wp = as.numeric(beta_wp),
      beta_bp = as.numeric(beta_bp),
      t_wp = as.numeric(t_wp),
      t_bp = as.numeric(t_bp),
      partial_r2_bp = as.numeric(partial_r2),
      p_bp = as.numeric(p_bp),
      personal_t_low_high = as.numeric(t_val),
      personal_p_low_high = as.numeric(p_val),
      personal_cor = as.numeric(cor_val),
      personal_cor_p = as.numeric(cor_p),
      personal_low_mean = as.numeric(low_mean),
      personal_high_mean = as.numeric(high_mean),
      stringsAsFactors = FALSE
    ),
    edge = edge
  )
}

message("[INFO] Fitting LMM (cognition) with age_wp + age_bp")
if (.Platform$OS.type == "windows") {
  message("[INFO] Windows parallel: ", num_cores, " workers")
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(
    cl,
    varlist = c("SCdata", "sc_cols", "base_by_sub", "q10", "q90", "pb_nsim", "pb_seed", "fit_edge", "pb_lmm_anova"),
    envir = environment()
  )
  res_list <- parallel::parLapply(
    cl, seq_along(sc_cols), fit_edge,
    data_all = SCdata, edges = sc_cols, base_by_sub = base_by_sub, q10 = q10, q90 = q90, pb_nsim = pb_nsim, pb_seed = pb_seed
  )
} else {
  res_list <- lapply(
    seq_along(sc_cols), fit_edge,
    data_all = SCdata, edges = sc_cols, base_by_sub = base_by_sub, q10 = q10, q90 = q90, pb_nsim = pb_nsim, pb_seed = pb_seed
  )
}

res_df <- do.call(rbind, lapply(res_list, `[[`, "row"))
res_df$p_bp_fdr <- p.adjust(res_df$p_bp, method = "fdr")
res_df$personal_p_low_high_fdr <- p.adjust(res_df$personal_p_low_high, method = "fdr")
res_df$personal_cor_p_fdr <- p.adjust(res_df$personal_cor_p, method = "fdr")
res_df$partial_r2_bp_3sd <- res_df$partial_r2_bp
pr_mean <- mean(res_df$partial_r2_bp, na.rm = TRUE)
pr_sd <- stats::sd(res_df$partial_r2_bp, na.rm = TRUE)
if (is.finite(pr_mean) && is.finite(pr_sd) && pr_sd > 0) {
  out_idx <- res_df$partial_r2_bp > pr_mean + 3 * pr_sd | res_df$partial_r2_bp < pr_mean - 3 * pr_sd
  res_df$partial_r2_bp_3sd[out_idx] <- NA_real_
}

out_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_cognition_", Cogvar_base, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)
saveRDS(res_df, out_rds)
write.csv(res_df, out_csv, row.names = FALSE)

message("[INFO] age_bp partial R2 matrix + S-A axis correlation")
mat_bp <- vec_to_mat(res_df$partial_r2_bp_3sd, ds = 12)
sig_bp <- vec_to_mat(res_df$p_bp_fdr < 0.05, ds = 12)
plot_matrix_sig(
  mat_bp,
  sig_bp,
  "Cognition age_bp partial R2",
  file.path(FigureFolder, paste0("matrix_age_bp_partialR2_", Cogvar_base, "_CV", CVthr))
)

SCrank.df.bp <- SCrankcorr(res_df, "partial_r2_bp_3sd", 12, dsdata = FALSE)
saveRDS(SCrank.df.bp, file.path(resultFolder, paste0("SCrankcorr_age_bp_partialR2_", Cogvar_base, "_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr (age_bp partial R2) r=", round(SCrank.df.bp$r.spearman, 3), " p=", signif(SCrank.df.bp$p.spearman, 3))

SCrank.data.bp <- SCrankcorr(res_df, "partial_r2_bp_3sd", 12, dsdata = TRUE)
limthr <- max(abs(SCrank.data.bp$partial_r2_bp_3sd), na.rm = TRUE)
if (!is.finite(limthr) || limthr == 0) limthr <- 1
scatterFig.bp <- ggplot(SCrank.data.bp) +
  geom_point(aes(x = SCrank, y = partial_r2_bp_3sd, color = partial_r2_bp_3sd), size = 5) +
  geom_smooth(aes(x = SCrank, y = partial_r2_bp_3sd), method = "lm", color = "black", linewidth = 1.4) +
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
ggsave(file.path(FigureFolder, paste0("scatter_age_bp_partialR2_vs_SCrank_", Cogvar_base, "_CV", CVthr, ".tiff")),
       scatterFig.bp, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_age_bp_partialR2_vs_SCrank_", Cogvar_base, "_CV", CVthr, ".pdf")),
       scatterFig.bp, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] age_wp personal slope t-value matrix (low10 vs high10)")
mat_t <- vec_to_mat(res_df$personal_t_low_high, ds = 12)
sig_t <- vec_to_mat(res_df$personal_p_low_high_fdr < 0.05, ds = 12)
plot_matrix_sig(
  mat_t,
  sig_t,
  "Personal slope t-value (low10 vs high10)",
  file.path(FigureFolder, paste0("matrix_personal_agewp_tvalue_low10_high10_", Cogvar_base, "_CV", CVthr))
)

message("[INFO] age_wp personal slope vs cognition correlation matrix + S-A axis")
mat_cor <- vec_to_mat(res_df$personal_cor, ds = 12)
sig_cor <- vec_to_mat(res_df$personal_cor_p_fdr < 0.05, ds = 12)
plot_matrix_sig(
  mat_cor,
  sig_cor,
  "Personal slope vs cognition (r)",
  file.path(FigureFolder, paste0("matrix_personal_agewp_corr_", Cogvar_base, "_CV", CVthr))
)

SCrank.df.cor <- SCrankcorr(res_df, "personal_cor", 12, dsdata = FALSE)
saveRDS(SCrank.df.cor, file.path(resultFolder, paste0("SCrankcorr_personal_agewp_corr_", Cogvar_base, "_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr (personal slope r) r=", round(SCrank.df.cor$r.spearman, 3), " p=", signif(SCrank.df.cor$p.spearman, 3))

SCrank.data.cor <- SCrankcorr(res_df, "personal_cor", 12, dsdata = TRUE)
limthr.cor <- max(abs(SCrank.data.cor$personal_cor), na.rm = TRUE)
if (!is.finite(limthr.cor) || limthr.cor == 0) limthr.cor <- 1
scatterFig.cor <- ggplot(SCrank.data.cor) +
  geom_point(aes(x = SCrank, y = personal_cor, color = personal_cor), size = 5) +
  geom_smooth(aes(x = SCrank, y = personal_cor), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr.cor, limthr.cor)) +
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
  labs(x = "S-A connectional axis rank", y = "Personal slope vs cognition (r)")
ggsave(file.path(FigureFolder, paste0("scatter_personal_agewp_corr_vs_SCrank_", Cogvar_base, "_CV", CVthr, ".tiff")),
       scatterFig.cor, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_personal_agewp_corr_vs_SCrank_", Cogvar_base, "_CV", CVthr, ".pdf")),
       scatterFig.cor, width = 15, height = 15, units = "cm", bg = "transparent")

message("[INFO] Decile bar plot (low/high personal slope)")
decile_df <- data.frame(
  SC_label = sc_cols,
  low_mean = res_df$personal_low_mean,
  high_mean = res_df$personal_high_mean,
  stringsAsFactors = FALSE
)
decile_df <- merge(decile_df, SA12_10, by.x = "SC_label", by.y = "SC_label", all.x = TRUE)
decile_summary <- decile_df %>%
  group_by(decile) %>%
  summarise(
    low_mean = mean(low_mean, na.rm = TRUE),
    high_mean = mean(high_mean, na.rm = TRUE),
    .groups = "drop"
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
  labs(x = "S-A decile", y = "Mean personal slope (age_wp)")

ggsave(
  file.path(FigureFolder, paste0("bar_personal_agewp_deciles_low_high_", Cogvar_base, "_CV", CVthr, ".tiff")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("bar_personal_agewp_deciles_low_high_", Cogvar_base, "_CV", CVthr, ".pdf")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)

message("[INFO] Done.")
