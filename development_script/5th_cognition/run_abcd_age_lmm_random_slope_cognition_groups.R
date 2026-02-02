#!/usr/bin/env Rscript

rm(list = ls())

# Prefer the active conda environment libraries, and avoid accidental user-library pollution on clusters.
Sys.unsetenv(c("R_LIBS_USER", "R_LIBS"))
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (nzchar(conda_prefix)) {
  .libPaths(c(file.path(conda_prefix, "lib", "R", "library"), .libPaths()))
}
.libPaths(.libPaths()[!grepl("/GPFS/.*/R/packages", .libPaths())])

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
if ("--help" %in% args) {
  cat("Usage: Rscript development_script/5th_cognition/run_abcd_age_lmm_random_slope_cognition_groups.R [options]\n",
      "Options:\n",
      "  --project-root PATH     Project root (default: auto-detect via script location)\n",
      "  --input-rds PATH        Longitudinal SCdata RDS (default: outputs/results/combat_gam/abcd/*age_sex_meanfd.rds)\n",
      "  --cog-rds PATH          Baseline cognition RDS (default: outputs/results/combat_gam/abcd/*combatgam_cognition.rds)\n",
      "  --result-dir PATH       Output results dir (default: outputs/results/5th_cognition/abcd/age_lmm)\n",
      "  --figure-dir PATH       Output figure dir (default: outputs/figures/5th_cognition/abcd/age_lmm)\n",
      "  --cvthr INT             CV threshold (default: 75)\n",
      "  --cogvar NAME           Cognition variable (default: nihtbx_fluidcomp_uncorrected)\n",
      "  --cogvar-base NAME      Baseline cognition variable name (default: nihtbx_fluidcomp_uncorrected_base)\n",
      "  --n-edges INT           Number of edges for quick tests (default: env N_EDGES or 78)\n",
      "\n",
      "Environment overrides (if args not provided): PROJECT_ROOT, INPUT_RDS, COG_RDS, RESULT_DIR, FIGURE_DIR, CVTHR, COGVAR, COGVAR_BASE, N_EDGES\n",
      sep = "")
  quit(save = "no", status = 0)
}

get_arg <- function(args, key, default = NULL) {
  hit <- which(args == key)
  if (length(hit) == 1 && length(args) >= hit + 1) return(args[hit + 1])
  default
}

is_abs_path <- function(path) {
  grepl("^[A-Za-z]:[\\\\/]", path) || grepl("^/", path)
}

resolve_path <- function(path, base_dir) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (is_abs_path(path)) return(normalizePath(path, mustWork = FALSE))
  normalizePath(file.path(base_dir, path), mustWork = FALSE)
}

get_script_path <- function() {
  arg_file <- sub("^--file=", "", full_args[grepl("^--file=", full_args)][1])
  if (length(arg_file) == 0 || !nzchar(arg_file)) return(NULL)
  normalizePath(arg_file, mustWork = FALSE)
}

find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, mustWork = FALSE)
  for (i in seq_len(8)) {
    if (file.exists(file.path(cur, "ARCHITECTURE.md"))) return(cur)
    parent <- dirname(cur)
    if (parent == cur) break
    cur <- parent
  }
  NULL
}

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
})

CVthr <- as.integer(get_arg(args, "--cvthr", Sys.getenv("CVTHR", unset = "75")))
if (is.na(CVthr)) CVthr <- 75
Cogvar <- get_arg(args, "--cogvar", Sys.getenv("COGVAR", unset = "nihtbx_fluidcomp_uncorrected"))
Cogvar_base <- get_arg(args, "--cogvar-base", Sys.getenv("COGVAR_BASE", unset = "nihtbx_fluidcomp_uncorrected_base"))

project_root_arg <- get_arg(args, "--project-root", Sys.getenv("PROJECT_ROOT", unset = ""))
script_path <- get_script_path()
script_dir <- if (!is.null(script_path)) dirname(script_path) else getwd()
project_root <- if (nzchar(project_root_arg)) {
  resolve_path(project_root_arg, getwd())
} else {
  find_project_root(script_dir)
}
if (is.null(project_root) || !file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Cannot find project root. Provide --project-root or PROJECT_ROOT env, or run under SCDevelopment.")
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- resolve_path(get_arg(args, "--result-dir", Sys.getenv("RESULT_DIR", unset = "")),
                             project_root)
FigureFolder <- resolve_path(get_arg(args, "--figure-dir", Sys.getenv("FIGURE_DIR", unset = "")),
                             project_root)
if (is.null(resultFolder)) {
  resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "age_lmm")
}
if (is.null(FigureFolder)) {
  FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "age_lmm")
}
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- resolve_path(get_arg(args, "--input-rds", Sys.getenv("INPUT_RDS", unset = "")),
                          project_root)
if (is.null(input_rds)) {
  input_rds <- file.path(
    project_root, "outputs", "results", "combat_gam", "abcd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
  )
}
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds,
       "\nProvide --input-rds or INPUT_RDS env, or generate via combat_gam/sbatch/abcd_combat_gam.sbatch (age/sex/mean_fd variant; longitudinal).")
}

cog_rds <- resolve_path(get_arg(args, "--cog-rds", Sys.getenv("COG_RDS", unset = "")),
                        project_root)
if (is.null(cog_rds)) {
  cog_rds <- file.path(
    project_root, "outputs", "results", "combat_gam", "abcd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
  )
}
if (!file.exists(cog_rds)) {
  stop("Missing cognition baseline input: ", cog_rds,
       "\nProvide --cog-rds or COG_RDS env, or generate via combat_gam/sbatch/abcd_combat_gam.sbatch (cognition variant).")
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
  stop("Missing eventname (required for baseline age): input has no eventname/scanID")
}

SCcog <- readRDS(cog_rds)
if (!("eventname" %in% names(SCcog)) && ("scanID" %in% names(SCcog))) {
  SCcog$eventname <- scanid_to_eventname(SCcog$scanID)
}

needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$sex <- as.factor(SCdata$sex)

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

# Prepare age (years) and baseline age
SCdata$age <- age_to_years(SCdata$age)
base_age <- tapply(SCdata$age, SCdata$subID, function(x) min(x, na.rm = TRUE))
SCdata$age_baseline <- base_age[match(as.character(SCdata$subID), names(base_age))]

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
n_edges <- as.integer(get_arg(args, "--n-edges", Sys.getenv("N_EDGES", unset = "78")))
if (is.na(n_edges) || n_edges < 1) n_edges <- 78
n_edges <- min(n_edges, 78L)
sc_cols <- sc_cols[seq_len(n_edges)]
out_suffix <- if (n_edges < 78) paste0("_N", n_edges) else ""
if (n_edges < 78) {
  message("[WARN] N_EDGES=", n_edges, " < 78: outputs are partial; matrices will be mostly NA (grey) except computed edges.")
}

# Subject groups based on baseline cognition (lowest/highest 10%)
sub_cog <- SCdata %>%
  select(subID, all_of(Cogvar_base)) %>%
  distinct() %>%
  filter(!is.na(.data[[Cogvar_base]]))
q10 <- quantile(sub_cog[[Cogvar_base]], 0.1, na.rm = TRUE)
q90 <- quantile(sub_cog[[Cogvar_base]], 0.9, na.rm = TRUE)
sub_low <- sub_cog$subID[sub_cog[[Cogvar_base]] <= q10]
sub_high <- sub_cog$subID[sub_cog[[Cogvar_base]] >= q90]

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
  # Random slopes are centered at 0; use mean absolute slope to summarize magnitude.
  rand_age_mean_abs <- mean(abs(re[, "age"]), na.rm = TRUE)
  list(ok = TRUE, beta = beta_age, t = t_age, rand_mean = rand_age_mean_abs)
}

vec_to_mat <- function(vec, ds = 12) {
  mat <- matrix(NA, ds, ds)
  idx <- which(lower.tri(mat, diag = TRUE))
  if (length(vec) > length(idx)) {
    stop("vec length exceeds lower-triangle size: ", length(vec), " > ", length(idx))
  }
  # Fill without recycling (supports N_EDGES < full matrix for quick tests).
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

  lwth <- min(df_melt$value, na.rm = TRUE)
  if (!is.finite(lwth)) {
    message("[WARN] Matrix values are all NA for: ", title, "; set lwth=-1 for plotting")
    lwth <- -1
  }
  if (lwth > 0) {
    # Ensure symmetric limits around 0 for positive-only matrices.
    lwth <- -lwth
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
    scale_fill_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(lwth, -lwth)) +
    scale_color_distiller(type = "seq", palette = "RdBu", na.value = "grey", limits = c(lwth, -lwth)) +
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
        file.path(resultFolder, paste0("age_lmm_random_slope_results_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
write.csv(res_all, file.path(resultFolder, paste0("age_lmm_random_slope_results_all_", Cogvar_base, "_CV", CVthr, out_suffix, ".csv")), row.names = FALSE)

# Correlation with S-A axis (full sample)
message("[INFO] Correlation to connectional axis (fixed age beta)")
SCrank.fixed <- SCrankcorr(res_all, "beta_age", 12, dsdata = FALSE)
saveRDS(SCrank.fixed, file.path(resultFolder, paste0("SCrankcorr_age_fixed_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
message("[INFO] SCrankcorr fixed r=", round(SCrank.fixed$r.spearman, 3), " p=", signif(SCrank.fixed$p.spearman, 3))

message("[INFO] Correlation to connectional axis (mean random slope)")
SCrank.rand <- SCrankcorr(res_all, "rand_age_mean", 12, dsdata = FALSE)
saveRDS(SCrank.rand, file.path(resultFolder, paste0("SCrankcorr_age_random_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
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
ggsave(file.path(FigureFolder, paste0("scatter_fixed_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")), p_fixed, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_fixed_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")), p_fixed, width = 15, height = 15, units = "cm", bg = "transparent")

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
ggsave(file.path(FigureFolder, paste0("scatter_random_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")), p_rand, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_random_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")), p_rand, width = 15, height = 15, units = "cm", bg = "transparent")

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
  file.path(resultFolder, paste0("age_lmm_matrices_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds"))
)

plot_matrix(mat_fixed_all, "Fixed age effect (all)", file.path(FigureFolder, paste0("matrix_fixed_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_rand_all, "Random age effect (all)", file.path(FigureFolder, paste0("matrix_random_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_rand_low, "Random age effect (low10)", file.path(FigureFolder, paste0("matrix_random_age_low10_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_rand_high, "Random age effect (high10)", file.path(FigureFolder, paste0("matrix_random_age_high10_", Cogvar_base, "_CV", CVthr, out_suffix)))

message("[INFO] Done.")
