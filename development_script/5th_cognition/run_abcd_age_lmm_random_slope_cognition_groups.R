#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
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
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "age_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "age_lmm")
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
if (!file.exists(plotdatasum_rds)) {
  stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)
}
plotdata <- readRDS(plotdatasum_rds)

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "lmm_age_random_slope.R"))

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
if ("eventname" %in% names(SCdata)) {
  message("[INFO] SCdata eventname table:\n", paste(capture.output(print(table(SCdata$eventname))), collapse = "\n"))
}

SCdata$age <- age_to_years(SCdata$age)
message(
  "[INFO] SCdata age range (years): ",
  round(min(SCdata$age, na.rm = TRUE), 3), "–", round(max(SCdata$age, na.rm = TRUE), 3)
)

needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
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

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))
n_edges <- as.integer(Sys.getenv("N_EDGES", unset = "78"))
if (is.na(n_edges) || n_edges < 1) n_edges <- 78
n_edges <- min(n_edges, 78L)
sc_cols <- sc_cols[seq_len(n_edges)]
out_suffix <- if (n_edges < 78) paste0("_N", n_edges) else ""
if (n_edges < 78) {
  message("[WARN] N_EDGES=", n_edges, " < 78: outputs are partial; matrices will be mostly NA (grey) except computed edges.")
}

# Scale SC strength by initial fit (ratio) before LMM.
SCdata.diw <- SCdata
for (region in sc_cols[seq_len(n_edges)]) {
  plotdata.tmp <- plotdata[plotdata$SC_label == region, , drop = FALSE]
  if (!("fit" %in% names(plotdata.tmp)) || nrow(plotdata.tmp) < 1 || is.na(plotdata.tmp$fit[[1]]) || plotdata.tmp$fit[[1]] == 0) {
    message("[WARN] Missing/invalid plotdata fit for ", region, "; skip scaling for this edge")
    next
  }
  SCdata.diw[[region]] <- SCdata[[region]] / plotdata.tmp$fit[[1]]
}
SCdata.diw[, sc_cols] <- lapply(SCdata.diw[, sc_cols, drop = FALSE], as.numeric)
SCdata <- SCdata.diw

run_all <- function() {
  dataname <- "SCdata_all"
  assign(dataname, SCdata, envir = .GlobalEnv)
  results <- lapply(sc_cols, function(edge) {
    lmm.age.random.slope(edge, dataname, return_model = TRUE)
  })
  stats <- dplyr::bind_rows(lapply(results, `[[`, "stats"))
  models <- lapply(results, `[[`, "model")
  names(models) <- sc_cols
  list(stats = stats, models = models)
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
  if (lwth > 0) lwth <- -lwth

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

message("[INFO] Fitting LMM per edge (random slope: age || subID)")
res_all_out <- run_all()
res_all <- res_all_out$stats
model_list <- res_all_out$models

saveRDS(res_all,
        file.path(resultFolder, paste0("age_lmm_random_slope_results_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
saveRDS(model_list,
        file.path(resultFolder, paste0("age_lmm_random_slope_models_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
write.csv(
  res_all,
  file.path(resultFolder, paste0("age_lmm_random_slope_results_all_", Cogvar_base, "_CV", CVthr, out_suffix, ".csv")),
  row.names = FALSE
)

message("[INFO] Correlation to connectional axis (fixed age beta)")
SCrank.fixed <- SCrankcorr(res_all, "beta_age", 12, dsdata = FALSE)
saveRDS(SCrank.fixed, file.path(resultFolder, paste0("SCrankcorr_age_fixed_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
message("[INFO] SCrankcorr fixed r=", round(SCrank.fixed$r.spearman, 3), " p=", signif(SCrank.fixed$p.spearman, 3))

message("[INFO] Correlation to connectional axis (mean random slope)")
SCrank.rand <- SCrankcorr(res_all, "rand_age_mean", 12, dsdata = FALSE)
saveRDS(SCrank.rand, file.path(resultFolder, paste0("SCrankcorr_age_random_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
message("[INFO] SCrankcorr random r=", round(SCrank.rand$r.spearman, 3), " p=", signif(SCrank.rand$p.spearman, 3))

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
ggsave(file.path(FigureFolder, paste0("scatter_fixed_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")),
       p_fixed, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_fixed_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")),
       p_fixed, width = 15, height = 15, units = "cm", bg = "transparent")

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
ggsave(file.path(FigureFolder, paste0("scatter_random_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")),
       p_rand, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_random_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")),
       p_rand, width = 15, height = 15, units = "cm", bg = "transparent")

mat_fixed_all <- vec_to_mat(res_all$beta_age)
mat_rand_all <- vec_to_mat(res_all$rand_age_mean)
saveRDS(
  list(
    fixed_all = mat_fixed_all,
    random_all = mat_rand_all
  ),
  file.path(resultFolder, paste0("age_lmm_matrices_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds"))
)

plot_matrix(mat_fixed_all, "Fixed age effect (all)", file.path(FigureFolder, paste0("matrix_fixed_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_rand_all, "Random age effect (all)", file.path(FigureFolder, paste0("matrix_random_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
message("[INFO] Done.")
