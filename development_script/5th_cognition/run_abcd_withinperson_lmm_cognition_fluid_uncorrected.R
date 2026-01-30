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
Cogvar <- "nihtbx_fluidcomp_uncorrected"
Cogvar_base <- "nihtbx_fluidcomp_uncorrected_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "withinperson_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "withinperson_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (age/sex/mean_fd variant; longitudinal)")
}

source(file.path(functionFolder, "lmminteraction.R"))

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
SCdata$sex <- as.factor(SCdata$sex)

# Baseline cognition comes from the cognition-specific ComBat output (baseline-only),
# then merged into the longitudinal SCdata by subID to keep "baseline fill" logic consistent.
cog_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
)
if (!file.exists(cog_rds)) {
  stop("Missing cognition baseline input: ", cog_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (cognition variant)")
}
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
if (all(is.na(SCdata[[Cogvar_base]]))) stop("Baseline cognition is all NA after join: ", Cogvar_base)

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

out_rds <- file.path(resultFolder, paste0("lmmresult_time_by_", Cogvar_base, "_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)

message("[INFO] Fitting LMM interaction models (n_edges=", n_edges, ", mc.cores=", num_cores, ", PB_METHOD=", pb_method, ", PB_NSIM=", pb_nsim, ")")
res_list <- parallel::mclapply(seq_len(n_edges), function(i) {
  region <- sc_cols[[i]]
  tryCatch(
    {
      out <- lmm.time.predict.covariateinteraction(
        region = region,
        dataname = "SCdata",
        age_var = "age",
        int_var = Cogvar_base,
        covariates = "sex+mean_fd",
        int_var_mode = "as_is",
        pb_method = pb_method,
        pb_nsim = pb_nsim,
        stats_only = TRUE
      )
      list(ok = TRUE, row = out, region = region, err = NA_character_)
    },
    error = function(e) list(ok = FALSE, row = NULL, region = region, err = conditionMessage(e))
  )
}, mc.cores = num_cores)

ok_mask <- vapply(res_list, function(z) isTRUE(z$ok), logical(1))
if (!all(ok_mask)) {
  first_bad <- which(!ok_mask)[[1]]
  stop("Edge failed: ", res_list[[first_bad]]$region, "\n", res_list[[first_bad]]$err)
}

lmmresult <- do.call(rbind, lapply(res_list, `[[`, "row"))
lmmresult$p_time_int_fdr <- p.adjust(lmmresult$p_time_int, method = "fdr")
saveRDS(lmmresult, out_rds)
write.csv(lmmresult, out_csv, row.names = FALSE)
message("[INFO] Saved: ", out_rds)

message("[INFO] totalstrength within-person change vs cognition (residualized for sex + mean_fd)")
SCdata$totalstrength <- rowMeans(SCdata[, sc_cols[seq_len(78)], drop = FALSE], na.rm = TRUE)

SCdata_time <- SCdata
base_age <- get_baseline_age(SCdata_time, subid_var = "subID", age_var = "age", event_var = "eventname")
age_years <- age_to_years(SCdata_time$age)
SCdata_time$time <- age_years - base_age[match(as.character(SCdata_time$subID), names(base_age))]

delta_df <- SCdata_time %>%
  select(subID, time, totalstrength, sex, mean_fd, all_of(Cogvar_base)) %>%
  filter(!is.na(.data[[Cogvar_base]]) & !is.na(totalstrength) & !is.na(time)) %>%
  group_by(subID) %>%
  summarise(
    time0 = min(time, na.rm = TRUE),
    time1 = max(time, na.rm = TRUE),
    ts0 = totalstrength[which.min(time)],
    ts1 = totalstrength[which.max(time)],
    mean_fd_mean = mean(mean_fd, na.rm = TRUE),
    sex = sex[which.max(!is.na(sex))][1],
    cognition_base = .data[[Cogvar_base]][which.max(!is.na(.data[[Cogvar_base]]))][1],
    .groups = "drop"
  ) %>%
  filter(is.finite(time1 - time0) & (time1 - time0) > 0)

delta_df$delta_totalstrength_per_year <- (delta_df$ts1 - delta_df$ts0) / (delta_df$time1 - delta_df$time0)
delta_df$sex <- as.factor(delta_df$sex)

lm_delta <- lm(delta_totalstrength_per_year ~ cognition_base + sex + mean_fd_mean, data = delta_df)
lm_delta_sum <- summary(lm_delta)
message(
  "[INFO] lm(delta) cognition beta=",
  signif(lm_delta_sum$coefficients["cognition_base", "Estimate"], 4),
  ", p=",
  signif(lm_delta_sum$coefficients["cognition_base", "Pr(>|t|)"], 4)
)

res_delta <- residuals(lm(delta_totalstrength_per_year ~ sex + mean_fd_mean, data = delta_df))
plot_df <- data.frame(
  delta_res = res_delta,
  cognition_base = delta_df$cognition_base
)

Fig <- ggplot(plot_df, aes(x = cognition_base, y = delta_res)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
  theme_classic() +
  labs(
    x = "Baseline cognition",
    y = "Within-person Δ totalstrength/year residual"
  )

ggsave(file.path(FigureFolder, paste0("delta_totalstrength_vs_", Cogvar_base, "_residualized.pdf")), Fig, width = 12, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("delta_totalstrength_vs_", Cogvar_base, "_residualized.tiff")), Fig, width = 12, height = 10, units = "cm", bg = "transparent", dpi = 600)

message("[INFO] Decile-wise within-person SC change vs cognition (S-A axis deciles)")
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
if (length(deciles) != 10) {
  stop("Expected 10 deciles, got: ", paste(deciles, collapse = ", "))
}

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

SCdata_dec <- SCdata_time
for (edge in sc_cols[seq_len(78)]) {
  f0 <- as.numeric(plot_fit[[edge]])
  if (is.na(f0) || !is.finite(f0) || f0 == 0) stop("Invalid plotdata fit for edge: ", edge)
  SCdata_dec[[edge]] <- as.numeric(SCdata_dec[[edge]]) / f0
}

for (d in deciles) {
  d_edges <- edges_by_decile[[as.character(d)]]
  SCdata_dec[[paste0("SC_decile", d)]] <- rowMeans(SCdata_dec[, d_edges, drop = FALSE], na.rm = TRUE)
}

decile_cols <- paste0("SC_decile", deciles)
needed_cols <- c("subID", "time", "sex", "mean_fd", Cogvar_base, decile_cols)
SCdata_dec <- SCdata_dec[, needed_cols, drop = FALSE]
SCdata_dec <- SCdata_dec[is.finite(SCdata_dec$time) & !is.na(SCdata_dec[[Cogvar_base]]), , drop = FALSE]

idx_by_sub <- split(seq_len(nrow(SCdata_dec)), as.character(SCdata_dec$subID))
rows <- vector("list", length(idx_by_sub))
names(rows) <- names(idx_by_sub)
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
    cognition_base = SCdata_dec[[Cogvar_base]][i0],
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

for (d in deciles) {
  dd <- delta_long[delta_long$decile == d, , drop = FALSE]
  dd <- dd[is.finite(dd$cognition_base) & is.finite(dd$res_delta), , drop = FALSE]
  ct <- suppressWarnings(stats::cor.test(dd$cognition_base, dd$res_delta))
  message(
    "[INFO] Decile ", d,
    " r=", signif(unname(ct$estimate), 4),
    " p=", signif(ct$p.value, 4),
    " n=", nrow(dd)
  )

  Fig_d <- ggplot(dd, aes(x = cognition_base, y = res_delta)) +
    geom_point(size = 1.0, alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
    theme_classic() +
    labs(
      x = "Baseline cognition",
      y = "Within-person Δ SC ratio/year residual"
    )

  out_base <- file.path(FigureFolder, sprintf("delta_SC_decile%02d_vs_%s_residualized", d, Cogvar_base))
  ggsave(paste0(out_base, ".pdf"), Fig_d, width = 12, height = 10, units = "cm", bg = "transparent")
  ggsave(paste0(out_base, ".tiff"), Fig_d, width = 12, height = 10, units = "cm", bg = "transparent", dpi = 600)
}

message("[INFO] Done.")
