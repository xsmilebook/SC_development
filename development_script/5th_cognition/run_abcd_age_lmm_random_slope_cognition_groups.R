#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpattern)
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

sa12_csv <- Sys.getenv(
  "ABCD_SA12_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "SA12_10.csv")
)
if (!file.exists(sa12_csv)) {
  stop("Missing ABCD_SA12_CSV: ", sa12_csv)
}
SA12_10 <- read.csv(sa12_csv, stringsAsFactors = FALSE)

source(file.path(functionFolder, "SCrankcorr.R"))
source(file.path(functionFolder, "lmm_age_random_slope.R"))

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
if ("eventname" %in% names(SCdata)) {
  message("[INFO] SCdata eventname table:\n", paste(capture.output(print(table(SCdata$eventname))), collapse = "\n"))
}

SCdata$age <- as.numeric(SCdata$age)
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

# Keep subjects with >=2 unique timepoints before fitting.
sub_time_ok <- tapply(SCdata$age, SCdata$subID, function(x) length(unique(round(x, 6))) >= 2)
SCdata <- SCdata[SCdata$subID %in% names(sub_time_ok)[sub_time_ok], , drop = FALSE]
if (nrow(SCdata) < 10) stop("Too few rows after >=2 timepoint filter: ", nrow(SCdata))

run_all <- function() {
  results <- lapply(sc_cols, function(edge) {
    df_edge <- SCdata
    y <- df_edge[[edge]]
    y_mean <- mean(y, na.rm = TRUE)
    y_sd <- stats::sd(y, na.rm = TRUE)
    if (is.finite(y_mean) && is.finite(y_sd) && y_sd > 0) {
      out_sub <- unique(df_edge$subID[which(y < y_mean - 2.69 * y_sd | y > y_mean + 2.69 * y_sd)])
      if (length(out_sub) > 0) df_edge <- df_edge[!df_edge$subID %in% out_sub, , drop = FALSE]
    }
    dataname <- paste0("SCdata_all_", edge)
    assign(dataname, df_edge, envir = .GlobalEnv)
    lmm.age.random.slope(edge, dataname, return_model = TRUE, return_slopes = TRUE)
  })
  stats <- dplyr::bind_rows(lapply(results, `[[`, "stats"))
  models <- lapply(results, `[[`, "model")
  slopes <- lapply(results, `[[`, "rand_slopes")
  names(models) <- sc_cols
  list(stats = stats, models = models, slopes = slopes)
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

message("[INFO] Fitting LMM per edge (random slope: age || subID)")
res_all_out <- run_all()
res_all <- res_all_out$stats
model_list <- res_all_out$models
slopes_list <- res_all_out$slopes

rand_slope_var <- vapply(seq_along(sc_cols), function(i) {
  m <- model_list[[i]]
  if (is.null(m) || !inherits(m, "lmerMod")) return(NA_real_)
  vc_df <- tryCatch(as.data.frame(lme4::VarCorr(m)), error = function(e) NULL)
  if (!is.data.frame(vc_df) || nrow(vc_df) < 1) return(NA_real_)
  idx <- which(vc_df$grp != "Residual" & vc_df$var1 == "age" & (is.na(vc_df$var2) | vc_df$var2 == ""))
  if (length(idx) < 1) return(NA_real_)
  as.numeric(vc_df$vcov[idx[[1]]])
}, numeric(1))
rand_slope_all_zero <- vapply(slopes_list, function(df) {
  if (!is.data.frame(df) || nrow(df) == 0 || !"random_slope" %in% names(df)) return(NA)
  rs <- as.numeric(df$random_slope)
  isTRUE(all(is.finite(rs) & abs(rs) < 1e-12))
}, logical(1))
res_all$rand_slope_var <- rand_slope_var
res_all$rand_slope_all_zero <- rand_slope_all_zero
zero_edges <- sc_cols[!is.na(rand_slope_all_zero) & rand_slope_all_zero]
message("[INFO] ranef_i(age)==0 edges: ", length(zero_edges), "/", length(sc_cols))
if (length(zero_edges) > 0) {
  writeLines(
    zero_edges,
    con = file.path(resultFolder, paste0("edges_rand_age_all_zero_", Cogvar_base, "_CV", CVthr, out_suffix, ".txt"))
  )
}

personal_slope <- vapply(slopes_list, function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(NA_real_)
  mean(df$fixed_slope + df$random_slope, na.rm = TRUE)
}, numeric(1))
res_all$personal_slope <- personal_slope

# low10/high10 groups for personal slope comparison
sub_cog <- SCdata %>%
  select(subID, all_of(Cogvar_base)) %>%
  distinct() %>%
  filter(!is.na(.data[[Cogvar_base]]))
q10 <- quantile(sub_cog[[Cogvar_base]], 0.1, na.rm = TRUE)
q90 <- quantile(sub_cog[[Cogvar_base]], 0.9, na.rm = TRUE)
sub_low <- sub_cog$subID[sub_cog[[Cogvar_base]] <= q10]
sub_high <- sub_cog$subID[sub_cog[[Cogvar_base]] >= q90]

personal_low_mean <- numeric(length(sc_cols))
personal_high_mean <- numeric(length(sc_cols))
personal_t <- numeric(length(sc_cols))
personal_p <- numeric(length(sc_cols))
for (i in seq_along(sc_cols)) {
  df <- slopes_list[[i]]
  if (!is.data.frame(df) || nrow(df) == 0) {
    personal_low_mean[i] <- NA_real_
    personal_high_mean[i] <- NA_real_
    personal_t[i] <- NA_real_
    personal_p[i] <- NA_real_
    next
  }
  df$personal <- df$fixed_slope + df$random_slope
  low_vals <- df$personal[df$subID %in% sub_low]
  high_vals <- df$personal[df$subID %in% sub_high]
  personal_low_mean[i] <- if (length(low_vals) > 0) mean(low_vals, na.rm = TRUE) else NA_real_
  personal_high_mean[i] <- if (length(high_vals) > 0) mean(high_vals, na.rm = TRUE) else NA_real_
  if (length(low_vals) > 1 && length(high_vals) > 1) {
    tt <- tryCatch(t.test(low_vals, high_vals), error = function(e) NULL)
    if (!is.null(tt)) {
      personal_t[i] <- unname(tt$statistic)
      personal_p[i] <- tt$p.value
    } else {
      personal_t[i] <- 0
      personal_p[i] <- 1
    }
  } else {
    personal_t[i] <- 0
    personal_p[i] <- 1
  }
}
res_all$personal_low10_mean <- personal_low_mean
res_all$personal_high10_mean <- personal_high_mean
res_all$personal_t_low_high <- personal_t
res_all$personal_p_low_high <- personal_p
res_all$personal_p_low_high_fdr <- p.adjust(res_all$personal_p_low_high, method = "fdr")

decile_df <- data.frame(
  SC_label = sc_cols,
  personal_low = res_all$personal_low10_mean,
  personal_high = res_all$personal_high10_mean,
  stringsAsFactors = FALSE
)
decile_df <- merge(decile_df, SA12_10, by.x = "SC_label", by.y = "SC_label", all.x = TRUE)
if (any(is.na(decile_df$decile))) {
  message("[WARN] Missing decile mapping for ", sum(is.na(decile_df$decile)), " edges in SA12_10.")
}
decile_summary <- decile_df %>%
  group_by(decile) %>%
  summarise(
    low_mean = mean(personal_low, na.rm = TRUE),
    high_mean = mean(personal_high, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  decile_summary,
  file.path(resultFolder, paste0("personal_slope_decile_low_high_", Cogvar_base, "_CV", CVthr, out_suffix, ".csv")),
  row.names = FALSE
)

saveRDS(res_all,
        file.path(resultFolder, paste0("age_lmm_random_slope_results_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
saveRDS(model_list,
        file.path(resultFolder, paste0("age_lmm_random_slope_models_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
saveRDS(slopes_list,
        file.path(resultFolder, paste0("age_lmm_random_slope_personal_slopes_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
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

message("[INFO] Correlation to connectional axis (mean personal slope)")
SCrank.personal <- SCrankcorr(res_all, "personal_slope", 12, dsdata = FALSE)
saveRDS(SCrank.personal, file.path(resultFolder, paste0("SCrankcorr_age_personal_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds")))
message("[INFO] SCrankcorr personal r=", round(SCrank.personal$r.spearman, 3), " p=", signif(SCrank.personal$p.spearman, 3))

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

SCrank.personal.df <- SCrankcorr(res_all, "personal_slope", 12, dsdata = TRUE)
limthr3 <- max(abs(SCrank.personal.df$personal_slope), na.rm = TRUE)
p_personal <- ggplot(SCrank.personal.df) +
  geom_point(aes(x = SCrank, y = personal_slope, color = personal_slope), size = 5) +
  geom_smooth(aes(x = SCrank, y = personal_slope), method = "lm", color = "black", linewidth = 1.4) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr3, limthr3)) +
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
  labs(x = "S-A connectional axis rank", y = "Mean personal age effect")
ggsave(file.path(FigureFolder, paste0("scatter_personal_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")),
       p_personal, width = 15, height = 15, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_personal_age_vs_SCrank_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")),
       p_personal, width = 15, height = 15, units = "cm", bg = "transparent")

mat_fixed_all <- vec_to_mat(res_all$beta_age)
mat_rand_all <- vec_to_mat(res_all$rand_age_mean)
mat_personal_all <- vec_to_mat(res_all$personal_slope)
mat_personal_low <- vec_to_mat(res_all$personal_low10_mean)
mat_personal_high <- vec_to_mat(res_all$personal_high10_mean)
mat_t_low_high <- vec_to_mat(res_all$personal_t_low_high)
sig_low_high <- vec_to_mat(res_all$personal_p_low_high_fdr < 0.05)
saveRDS(
  list(
    fixed_all = mat_fixed_all,
    random_all = mat_rand_all,
    personal_all = mat_personal_all,
    personal_low10 = mat_personal_low,
    personal_high10 = mat_personal_high,
    personal_t_low_high = mat_t_low_high,
    personal_sig_low_high = sig_low_high
  ),
  file.path(resultFolder, paste0("age_lmm_matrices_", Cogvar_base, "_CV", CVthr, out_suffix, ".rds"))
)

plot_matrix(mat_fixed_all, "Fixed age effect (all)", file.path(FigureFolder, paste0("matrix_fixed_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_rand_all, "Random age effect (all)", file.path(FigureFolder, paste0("matrix_random_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_personal_all, "Personal age effect (all)", file.path(FigureFolder, paste0("matrix_personal_age_all_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_personal_low, "Personal age effect (low10)", file.path(FigureFolder, paste0("matrix_personal_age_low10_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix(mat_personal_high, "Personal age effect (high10)", file.path(FigureFolder, paste0("matrix_personal_age_high10_", Cogvar_base, "_CV", CVthr, out_suffix)))
plot_matrix_sig(mat_t_low_high, sig_low_high, "Personal slope t-value (low10 vs high10)", file.path(FigureFolder, paste0("matrix_personal_age_tvalue_low10_high10_", Cogvar_base, "_CV", CVthr, out_suffix)))

# Decile bar plot (low vs high personal slope)
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

dodge_width <- 0.8
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
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  guides(
    alpha = guide_legend(override.aes = list(fill = "grey70", color = "black", pattern = c("none", "stripe"))),
    linetype = guide_legend(override.aes = list(fill = "grey70", color = "black"))
  ) +
  labs(x = "S-A decile", y = "Mean personal age effect")

ggsave(
  file.path(FigureFolder, paste0("bar_personal_age_deciles_low_high_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("bar_personal_age_deciles_low_high_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)

# Decile bar plot for random slope (low vs high)
decile_df_rand <- data.frame(
  SC_label = sc_cols,
  random_low = vapply(slopes_list, function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(NA_real_)
    mean(df$random_slope[df$subID %in% sub_low], na.rm = TRUE)
  }, numeric(1)),
  random_high = vapply(slopes_list, function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(NA_real_)
    mean(df$random_slope[df$subID %in% sub_high], na.rm = TRUE)
  }, numeric(1)),
  stringsAsFactors = FALSE
)
decile_df_rand <- decile_df_rand %>%
  filter(!(is.finite(random_low) & is.finite(random_high) & random_low == 0 & random_high == 0))
decile_df_rand <- merge(decile_df_rand, SA12_10, by.x = "SC_label", by.y = "SC_label", all.x = TRUE)
decile_summary_rand <- decile_df_rand %>%
  group_by(decile) %>%
  summarise(
    low_mean = mean(random_low, na.rm = TRUE),
    high_mean = mean(random_high, na.rm = TRUE),
    .groups = "drop"
  )
plotdf_r <- data.frame(
  decile = decile_summary_rand$decile,
  high = decile_summary_rand$high_mean,
  low = decile_summary_rand$low_mean
)
plotdf_r_long <- reshape(
  plotdf_r,
  varying = c("high", "low"),
  v.names = "mean",
  timevar = "group",
  times = c("high", "low"),
  direction = "long"
)
plotdf_r_long$decile <- as.integer(plotdf_r_long$decile)
plotdf_r_long$group <- factor(plotdf_r_long$group, levels = c("high", "low"))

bar_fig_rand <- ggplot(plotdf_r_long, aes(x = factor(decile), y = mean, fill = factor(decile), group = group)) +
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
  labs(x = "S-A decile", y = "Mean random age effect")

ggsave(
  file.path(FigureFolder, paste0("bar_random_age_deciles_low_high_", Cogvar_base, "_CV", CVthr, out_suffix, ".tiff")),
  bar_fig_rand, width = 18, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("bar_random_age_deciles_low_high_", Cogvar_base, "_CV", CVthr, out_suffix, ".pdf")),
  bar_fig_rand, width = 18, height = 12, units = "cm", bg = "transparent"
)
message("[INFO] Done.")
