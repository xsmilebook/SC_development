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
  library(pbkrtest)
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

CVthr <- 75
int_var <- "GENERAL"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "6th_pfactor", "abcd", "withinperson_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "6th_pfactor", "abcd", "withinperson_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)")
}

source(file.path(functionFolder, "lmminteraction.R"))
source(file.path(functionFolder, "SCrankcorr.R"))

SCdata <- readRDS(input_rds)
needed <- c("subID", "age", "sex", "mean_fd", int_var)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$sex <- as.factor(SCdata$sex)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (length(sc_cols) < 78) stop("Expected >=78 SC.* columns, got: ", length(sc_cols))

n_edges <- as.integer(Sys.getenv("N_EDGES", unset = "78"))
if (is.na(n_edges) || n_edges < 1) n_edges <- 78
n_edges <- min(n_edges, 78L)

pb_method <- Sys.getenv("PB_METHOD", unset = "KR")
pb_nsim <- as.integer(Sys.getenv("PB_NSIM", unset = "1000"))
if (is.na(pb_nsim) || pb_nsim < 1) pb_nsim <- 1000

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "40"))
if (is.na(num_cores) || num_cores < 1) num_cores <- 40
num_cores <- min(num_cores, 40L)

out_rds <- file.path(resultFolder, paste0("lmmresult_time_by_pfactor_", int_var, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)

message("[INFO] Fitting LMM interaction models (n_edges=", n_edges, ", mc.cores=", num_cores, ", PB_METHOD=", pb_method, ", PB_NSIM=", pb_nsim, ")")
res_list <- parallel::mclapply(seq_len(n_edges), function(i) {
  region <- sc_cols[[i]]
  lmm.time.predict.covariateinteraction(
    region = region,
    dataname = "SCdata",
    age_var = "age",
    int_var = int_var,
    covariates = "sex+mean_fd",
    int_var_mode = "as_is",
    pb_method = pb_method,
    pb_nsim = pb_nsim,
    stats_only = TRUE
  )
}, mc.cores = num_cores)

if (any(vapply(res_list, inherits, logical(1), what = "try-error"))) {
  stop("At least one edge failed (mclapply returned try-error). Set N_EDGES small and rerun to locate the failing edge.")
}

lmmresult <- do.call(rbind, res_list)
lmmresult$p_time_int_fdr <- p.adjust(lmmresult$p_time_int, method = "fdr")
saveRDS(lmmresult, out_rds)
write.csv(lmmresult, out_csv, row.names = FALSE)
message("[INFO] Saved: ", out_rds)

message("[INFO] Correlation to connectional axis (t_time_int)")
SCrank.df.t <- SCrankcorr(lmmresult, "t_time_int", 12, dsdata = FALSE)
saveRDS(SCrank.df.t, file.path(resultFolder, paste0("SCrankcorr_lmm_time_by_pfactor_", int_var, "_CV", CVthr, "_tvalue.rds")))
message("[INFO] SCrankcorr (t) r=", round(SCrank.df.t$r.spearman, 3), " p=", signif(SCrank.df.t$p.spearman, 3))

message("[INFO] Scatter plot: t_time_int vs S-A rank")
SCrank.data.t <- SCrankcorr(lmmresult, "t_time_int", 12, dsdata = TRUE)
limthr.t <- max(abs(SCrank.data.t$t_time_int), na.rm = TRUE)
scatterFig.t <- ggplot(data = SCrank.data.t) +
  geom_point(aes(x = SCrank, y = t_time_int, color = t_time_int), size = 3) +
  geom_smooth(aes(x = SCrank, y = t_time_int), method = "lm", color = "black", linewidth = 0.9) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-limthr.t, limthr.t)) +
  theme_classic() +
  labs(x = "S-A connectional axis rank", y = "LMM t value")
ggsave(file.path(FigureFolder, paste0("scatter_tvalue_vs_SCrank_lmm_pfactor_", int_var, "_CV", CVthr, ".pdf")), scatterFig.t, width = 12, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("scatter_tvalue_vs_SCrank_lmm_pfactor_", int_var, "_CV", CVthr, ".tiff")), scatterFig.t, width = 12, height = 10, units = "cm", bg = "transparent", dpi = 600)

SCdata$totalstrength <- rowMeans(SCdata[, sc_cols[seq_len(78)], drop = FALSE], na.rm = TRUE)
out_total <- lmm.time.predict.covariateinteraction(
  region = "totalstrength",
  dataname = "SCdata",
  age_var = "age",
  int_var = int_var,
  covariates = "sex+mean_fd",
  int_var_mode = "as_is",
  pb_method = pb_method,
  pb_nsim = pb_nsim,
  stats_only = FALSE
)
pred <- out_total$predicted
pred$lower <- pred$.fitted - 1.96 * pred$.se
pred$upper <- pred$.fitted + 1.96 * pred$.se

Fig <- ggplot(pred, aes(x = time, y = .fitted, group = int_level, linetype = int_level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = int_level), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9, color = "black") +
  geom_point(size = 2, color = "black") +
  theme_classic() +
  scale_linetype_manual(values = c(low10 = "dashed", high90 = "solid")) +
  scale_fill_manual(values = c(low10 = "grey60", high90 = "grey20")) +
  labs(x = "Years from baseline", y = "Predicted totalstrength (fixed effects only)", linetype = NULL, fill = NULL)

ggsave(file.path(FigureFolder, paste0("pred_totalstrength_time_by_", int_var, "_low10_high90.pdf")), Fig, width = 12, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("pred_totalstrength_time_by_", int_var, "_low10_high90.tiff")), Fig, width = 12, height = 10, units = "cm", bg = "transparent", dpi = 600)

message("[INFO] Decile-wise within-person SC change vs pfactor (S-A axis deciles)")
sa12_csv <- Sys.getenv(
  "ABCD_SA12_CSV",
  unset = file.path(project_root, "wd", "interdataFolder_ABCD", "SA12_10.csv")
)
if (!file.exists(sa12_csv)) stop("Missing ABCD_SA12_CSV: ", sa12_csv)
SA12_10 <- read.csv(sa12_csv, stringsAsFactors = FALSE)
needed_sa <- c("SC_label", "decile")
missing_sa <- setdiff(needed_sa, names(SA12_10))
if (length(missing_sa) > 0) stop("Missing columns in SA12_10: ", paste(missing_sa, collapse = ", "))

edge_map <- SA12_10 %>%
  select(SC_label, decile) %>%
  filter(SC_label %in% sc_cols) %>%
  distinct()
if (nrow(edge_map) < 78) {
  missing_edges <- setdiff(sc_cols, edge_map$SC_label)
  stop("SA12_10 missing edges: ", paste(head(missing_edges, 10), collapse = ", "))
}
edges_by_decile <- split(edge_map$SC_label, edge_map$decile)
deciles <- sort(unique(as.integer(names(edges_by_decile))))
if (length(deciles) != 10) stop("Expected 10 deciles, got: ", paste(deciles, collapse = ", "))

# Normalize SC strength to ratio (divide by initial fit) to match previous figures.
plotdatasum_rds <- Sys.getenv(
  "ABCD_PLOTDATASUM_RDS",
  unset = file.path(
    project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
    "abcd", "combat_gam", "CV75", "plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds"
  )
)
if (!file.exists(plotdatasum_rds)) stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)
plotdata <- readRDS(plotdatasum_rds)
if (!all(c("SC_label", "fit") %in% names(plotdata))) stop("plotdata missing SC_label/fit: ", plotdatasum_rds)
plot_fit <- plotdata$fit
names(plot_fit) <- as.character(plotdata$SC_label)
missing_fit <- setdiff(sc_cols, names(plot_fit))
if (length(missing_fit) > 0) stop("plotdata missing fits for edges: ", paste(head(missing_fit, 10), collapse = ", "))

scanid_to_eventname <- function(scanID) {
  sess <- sub("^.*_ses-", "", as.character(scanID))
  sess <- gsub("([a-z])([A-Z])", "\\1_\\2", sess)
  sess <- gsub("([A-Za-z])([0-9])", "\\1_\\2", sess)
  sess <- gsub("([0-9])([A-Za-z])", "\\1_\\2", sess)
  tolower(sess)
}
if (!("eventname" %in% names(SCdata)) && ("scanID" %in% names(SCdata))) {
  SCdata$eventname <- scanid_to_eventname(SCdata$scanID)
}

SCdata_dec <- SCdata
for (edge in sc_cols[seq_len(78)]) {
  f0 <- as.numeric(plot_fit[[edge]])
  if (is.na(f0) || !is.finite(f0) || f0 == 0) stop("Invalid plotdata fit for edge: ", edge)
  SCdata_dec[[edge]] <- as.numeric(SCdata_dec[[edge]]) / f0
}

base_age <- get_baseline_age(SCdata_dec, subid_var = "subID", age_var = "age", event_var = "eventname")
age_years <- age_to_years(SCdata_dec$age)
SCdata_dec$time <- age_years - base_age[match(as.character(SCdata_dec$subID), names(base_age))]

for (d in deciles) {
  d_edges <- edges_by_decile[[as.character(d)]]
  SCdata_dec[[paste0("SC_decile", d)]] <- rowMeans(SCdata_dec[, d_edges, drop = FALSE], na.rm = TRUE)
}

decile_cols <- paste0("SC_decile", deciles)
needed_cols <- c("subID", "time", "sex", "mean_fd", int_var, decile_cols)
SCdata_dec <- SCdata_dec[, needed_cols, drop = FALSE]
SCdata_dec <- SCdata_dec[is.finite(SCdata_dec$time) & !is.na(SCdata_dec[[int_var]]), , drop = FALSE]

idx_by_sub <- split(seq_len(nrow(SCdata_dec)), as.character(SCdata_dec$subID))
rows <- vector("list", length(idx_by_sub))
for (k in seq_along(idx_by_sub)) {
  ii <- idx_by_sub[[k]]
  tsub <- SCdata_dec$time[ii]
  if (length(unique(round(tsub, 6))) < 2) next
  i0 <- ii[which.min(tsub)]
  i1 <- ii[which.max(tsub)]
  dt <- SCdata_dec$time[i1] - SCdata_dec$time[i0]
  if (!is.finite(dt) || dt <= 0) next

  out <- data.frame(
    subID = SCdata_dec$subID[i0],
    sex = SCdata_dec$sex[i0],
    mean_fd_mean = mean(SCdata_dec$mean_fd[ii], na.rm = TRUE),
    pfactor_base = SCdata_dec[[int_var]][i0],
    delta_t = dt,
    stringsAsFactors = FALSE
  )
  for (d in deciles) {
    v0 <- SCdata_dec[[paste0("SC_decile", d)]][i0]
    v1 <- SCdata_dec[[paste0("SC_decile", d)]][i1]
    out[[paste0("delta_decile", d, "_per_year")]] <- (v1 - v0) / dt
  }
  rows[[k]] <- out
}
delta_dec <- dplyr::bind_rows(rows)
delta_dec$sex <- as.factor(delta_dec$sex)

delta_long <- delta_dec %>%
  pivot_longer(
    cols = starts_with("delta_decile") & ends_with("_per_year"),
    names_to = "decile",
    values_to = "delta_per_year"
  )
delta_long$decile <- as.integer(gsub("^delta_decile([0-9]+)_per_year$", "\\1", delta_long$decile))

delta_long$res_delta <- NA_real_
for (d in deciles) {
  idx <- which(delta_long$decile == d)
  dd <- delta_long[idx, , drop = FALSE]
  delta_long$res_delta[idx] <- residuals(lm(delta_per_year ~ sex + mean_fd_mean, data = dd))
}

rp_rows <- vector("list", length(deciles))
for (j in seq_along(deciles)) {
  d <- deciles[[j]]
  dd <- delta_long[delta_long$decile == d, , drop = FALSE]
  dd <- dd[is.finite(dd$pfactor_base) & is.finite(dd$res_delta), , drop = FALSE]
  ct <- suppressWarnings(stats::cor.test(dd$pfactor_base, dd$res_delta))
  message(
    "[INFO] Decile ", d,
    " r=", signif(unname(ct$estimate), 4),
    " p=", signif(ct$p.value, 4),
    " n=", nrow(dd)
  )
  rp_rows[[j]] <- data.frame(
    decile = d,
    r = unname(ct$estimate),
    p = ct$p.value,
    n = nrow(dd)
  )
}
rp_df <- do.call(rbind, rp_rows)
write.csv(rp_df, file.path(resultFolder, paste0("delta_SC_deciles_vs_pfactor_", int_var, "_residualized_rp.csv")), row.names = FALSE)

Fig_dec <- ggplot(delta_long, aes(x = pfactor_base, y = res_delta)) +
  geom_point(size = 0.9, alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, color = "black") +
  facet_wrap(~decile, ncol = 5) +
  theme_classic() +
  labs(
    x = "p-factor",
    y = "Within-person Δ SC ratio/year residual"
  )

out_base <- file.path(FigureFolder, paste0("delta_SC_deciles_vs_pfactor_", int_var, "_residualized"))
ggsave(paste0(out_base, ".pdf"), Fig_dec, width = 22, height = 12, units = "cm", bg = "transparent")
ggsave(paste0(out_base, ".tiff"), Fig_dec, width = 22, height = 12, units = "cm", bg = "transparent", dpi = 600)

message("[INFO] Done.")
