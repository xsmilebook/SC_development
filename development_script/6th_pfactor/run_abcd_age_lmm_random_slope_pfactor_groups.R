#!/usr/bin/env Rscript

rm(list = ls())

# Prefer the active conda environment libraries, and avoid accidental user-library pollution on clusters.
Sys.unsetenv(c("R_LIBS_USER", "R_LIBS"))
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (nzchar(conda_prefix)) {
  .libPaths(c(file.path(conda_prefix, "lib", "R", "library"), .libPaths()))
}
.libPaths(.libPaths()[!grepl("/GPFS/.*/R/packages", .libPaths())])

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
})

CVthr <- 75
Pvar <- "GENERAL"
Pvar_base <- "GENERAL_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "age_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "6th_pfactor", "abcd", "age_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
}

source(file.path(functionFolder, "SCrankcorr.R"))

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
  if (is.finite(mx) && mx > 24) return(age_num / 12)
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
  stop("Missing eventname (required for baseline grouping): input has no eventname/scanID")
}

needed <- c("subID", "age", "sex", "mean_fd", Pvar)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$sex <- as.factor(SCdata$sex)

# Baseline p-factor values for grouping
Pfactor_df <- SCdata %>%
  select(subID, eventname, all_of(Pvar)) %>%
  filter(!is.na(.data[[Pvar]])) %>%
  filter(grepl("base", eventname, ignore.case = TRUE)) %>%
  select(subID, !!Pvar_base := all_of(Pvar)) %>%
  distinct()
if (nrow(Pfactor_df) < 1) stop("No baseline p-factor rows found in: ", input_rds)
SCdata <- SCdata %>% left_join(Pfactor_df, by = "subID")
if (!Pvar_base %in% names(SCdata)) stop("Baseline p-factor join failed, missing: ", Pvar_base)

# Prepare age (years) and baseline age
SCdata$age <- age_to_years(SCdata$age)
base_age <- tapply(SCdata$age, SCdata$subID, function(x) min(x, na.rm = TRUE))
SCdata$age_baseline <- base_age[match(as.character(SCdata$subID), names(base_age))]

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
n_edges <- as.integer(Sys.getenv("N_EDGES", unset = "78"))
if (is.na(n_edges) || n_edges < 1) n_edges <- 78
n_edges <- min(n_edges, 78L)
sc_cols <- sc_cols[seq_len(n_edges)]

# Subject groups based on baseline p-factor (lowest/highest 10%)
sub_pf <- SCdata %>%
  select(subID, all_of(Pvar_base)) %>%
  distinct() %>%
  filter(!is.na(.data[[Pvar_base]]))
q10 <- quantile(sub_pf[[Pvar_base]], 0.1, na.rm = TRUE)
q90 <- quantile(sub_pf[[Pvar_base]], 0.9, na.rm = TRUE)
sub_low <- sub_pf$subID[sub_pf[[Pvar_base]] <= q10]
sub_high <- sub_pf$subID[sub_pf[[Pvar_base]] >= q90]

groups <- list(
  all = unique(SCdata$subID),
  low10 = sub_low,
  high10 = sub_high
)

fit_edge <- function(df, edge_col) {
  # Drop rows with missing covariates or response
  df <- df[!is.na(df$age) & !is.na(df$sex) & !is.na(df$mean_fd) & !is.na(df[[edge_col]]), , drop = FALSE]

  # Remove outliers (3 SD) in edge
  y <- df[[edge_col]]
  out_y <- which(y < mean(y, na.rm = TRUE) - 3 * sd(y, na.rm = TRUE) |
                   y > mean(y, na.rm = TRUE) + 3 * sd(y, na.rm = TRUE))
  if (length(out_y) > 0) df <- df[-out_y, , drop = FALSE]

  # Keep subjects with >=2 timepoints
  ok <- tapply(df$age, df$subID, function(x) length(unique(round(x, 6))) >= 2)
  df <- df[df$subID %in% names(ok)[ok], , drop = FALSE]
  if (nrow(df) < 10) return(list(ok = FALSE, beta = NA, t = NA, rand_mean = NA))

  df$sex <- as.factor(df$sex)
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    check.nobs.vs.nRE = "ignore",
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev = "ignore"
  )
  fml <- as.formula(sprintf("%s ~ age + sex + mean_fd + (1 + age | subID)", edge_col))
  mod <- suppressWarnings(lmer(fml, data = df, REML = FALSE, control = ctrl))
  sm <- summary(mod)
  beta_age <- sm$coefficients["age", "Estimate"]
  t_age <- sm$coefficients["age", "t value"]
  re <- ranef(mod)$subID
  rand_age_mean <- mean(re[,"age"], na.rm = TRUE)
  list(ok = TRUE, beta = beta_age, t = t_age, rand_mean = rand_age_mean)
}

vec_to_mat <- function(vec, ds = 12) {
  mat <- matrix(NA, ds, ds)
  mat[lower.tri(mat, diag = TRUE)] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}

plot_matrix <- function(mat, title, out_base) {
  df_melt <- as.data.frame(as.table(mat))
  names(df_melt) <- c("nodeid", "variable", "value")
  df_melt$nodeid <- as.numeric(as.character(df_melt$nodeid))
  df_melt$variable <- as.numeric(as.character(df_melt$variable))
  df_melt$nodeid <- -df_melt$nodeid
  df_melt$value <- as.numeric(df_melt$value)

  limthr <- max(abs(df_melt$value), na.rm = TRUE)
  if (!is.finite(limthr) || limthr == 0) {
    message("[WARN] Matrix values are all NA/zero for: ", title, "; set limthr=1 for plotting")
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

run_group <- function(group_name, sub_ids) {
  df <- SCdata[SCdata$subID %in% sub_ids, , drop = FALSE]
  results <- lapply(sc_cols, function(edge) {
    out <- fit_edge(df, edge)
    data.frame(edge = edge, ok = out$ok, beta_age = out$beta, t_age = out$t, rand_age_mean = out$rand_mean)
  })
  do.call(rbind, results)
}

message("[INFO] Fitting LMM per edge for all/low10/high10 groups")
res_all <- run_group("all", groups$all)
res_low <- run_group("low10", groups$low10)
res_high <- run_group("high10", groups$high10)

saveRDS(list(all = res_all, low10 = res_low, high10 = res_high),
        file.path(resultFolder, paste0("age_lmm_random_slope_results_", Pvar_base, "_CV", CVthr, ".rds")))
write.csv(res_all, file.path(resultFolder, paste0("age_lmm_random_slope_results_all_", Pvar_base, "_CV", CVthr, ".csv")), row.names = FALSE)

# Correlation with S-A axis (full sample)
message("[INFO] Correlation to connectional axis (fixed age beta)")
SCrank.fixed <- SCrankcorr(res_all, "beta_age", 12, dsdata = FALSE)
saveRDS(SCrank.fixed, file.path(resultFolder, paste0("SCrankcorr_age_fixed_", Pvar_base, "_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr fixed r=", round(SCrank.fixed$r.spearman, 3), " p=", signif(SCrank.fixed$p.spearman, 3))

message("[INFO] Correlation to connectional axis (mean random slope)")
SCrank.rand <- SCrankcorr(res_all, "rand_age_mean", 12, dsdata = FALSE)
saveRDS(SCrank.rand, file.path(resultFolder, paste0("SCrankcorr_age_random_", Pvar_base, "_CV", CVthr, ".rds")))
message("[INFO] SCrankcorr random r=", round(SCrank.rand$r.spearman, 3), " p=", signif(SCrank.rand$p.spearman, 3))

# Scatter plots (full sample)
SCrank.fixed.df <- SCrankcorr(res_all, "beta_age", 12, dsdata = TRUE)
limthr <- max(abs(SCrank.fixed.df$beta_age), na.rm = TRUE)
p_fixed <- ggplot(SCrank.fixed.df) +
  geom_point(aes(x = SCrank, y = beta_age, color = beta_age), size = 5) +
  geom_smooth(aes(x = SCrank, y = beta_age), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr, limthr)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 23.4, color = "black"),
    axis.title = element_text(size = 23.4),
    aspect.ratio = 1,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  ) +
  labs(x = "S-A connectional axis rank", y = "Fixed age effect (beta)")

ggsave(file.path(FigureFolder, paste0("scatter_fixed_age_vs_SCrank_", Pvar_base, "_CV", CVthr, ".tiff")), p_fixed, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_fixed_age_vs_SCrank_", Pvar_base, "_CV", CVthr, ".pdf")), p_fixed, width = 15, height = 15, units = "cm", bg = "transparent")

SCrank.rand.df <- SCrankcorr(res_all, "rand_age_mean", 12, dsdata = TRUE)
limthr2 <- max(abs(SCrank.rand.df$rand_age_mean), na.rm = TRUE)
p_rand <- ggplot(SCrank.rand.df) +
  geom_point(aes(x = SCrank, y = rand_age_mean, color = rand_age_mean), size = 5) +
  geom_smooth(aes(x = SCrank, y = rand_age_mean), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr2, limthr2)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 23.4, color = "black"),
    axis.title = element_text(size = 23.4),
    aspect.ratio = 1,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  ) +
  labs(x = "S-A connectional axis rank", y = "Mean random age effect")

ggsave(file.path(FigureFolder, paste0("scatter_random_age_vs_SCrank_", Pvar_base, "_CV", CVthr, ".tiff")), p_rand, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_random_age_vs_SCrank_", Pvar_base, "_CV", CVthr, ".pdf")), p_rand, width = 15, height = 15, units = "cm", bg = "transparent")

# Matrices (full sample, low/high groups)
mat_fixed_all <- vec_to_mat(res_all$beta_age)
mat_rand_all <- vec_to_mat(res_all$rand_age_mean)
mat_fixed_low <- vec_to_mat(res_low$beta_age)
mat_rand_low <- vec_to_mat(res_low$rand_age_mean)
mat_fixed_high <- vec_to_mat(res_high$beta_age)
mat_rand_high <- vec_to_mat(res_high$rand_age_mean)

saveRDS(
  list(
    fixed_all = mat_fixed_all,
    random_all = mat_rand_all,
    fixed_low10 = mat_fixed_low,
    random_low10 = mat_rand_low,
    fixed_high10 = mat_fixed_high,
    random_high10 = mat_rand_high
  ),
  file.path(resultFolder, paste0("age_lmm_matrices_", Pvar_base, "_CV", CVthr, ".rds"))
)

plot_matrix(mat_fixed_all, "Fixed age effect (all)", file.path(FigureFolder, paste0("matrix_fixed_age_all_", Pvar_base, "_CV", CVthr)))
plot_matrix(mat_rand_all, "Random age effect (all)", file.path(FigureFolder, paste0("matrix_random_age_all_", Pvar_base, "_CV", CVthr)))
plot_matrix(mat_fixed_low, "Fixed age effect (low10)", file.path(FigureFolder, paste0("matrix_fixed_age_low10_", Pvar_base, "_CV", CVthr)))
plot_matrix(mat_rand_low, "Random age effect (low10)", file.path(FigureFolder, paste0("matrix_random_age_low10_", Pvar_base, "_CV", CVthr)))
plot_matrix(mat_fixed_high, "Fixed age effect (high10)", file.path(FigureFolder, paste0("matrix_fixed_age_high10_", Pvar_base, "_CV", CVthr)))
plot_matrix(mat_rand_high, "Random age effect (high10)", file.path(FigureFolder, paste0("matrix_random_age_high10_", Pvar_base, "_CV", CVthr)))

message("[INFO] Done.")
