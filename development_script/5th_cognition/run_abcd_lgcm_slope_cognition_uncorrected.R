#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
})

rm(list = ls())

CVthr <- 75
Cogvar <- "nihtbx_fluidcomp_uncorrected"
Cogvar_base <- "nihtbx_fluidcomp_uncorrected_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "lgcm_slope")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "lgcm_slope")
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

needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
if (!("sex" %in% names(SCdata))) stop("Missing required column: sex")
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
    cog_base = SCdata[[Cogvar_base]][i0],
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
    cog_base = dat_sub$cog_base,
    sex = dat_sub$sex,
    mean_fd_t0 = dat_sub$mean_fd_t0,
    mean_fd_t1 = dat_sub$mean_fd_t1
  )
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10) {
    res_rows[[i]] <- data.frame(
      edge = edge,
      n_sub = nrow(df),
      beta_cog = NA_real_,
      t_cog = NA_real_,
      p_cog = NA_real_,
      stringsAsFactors = FALSE
    )
    pred_low[i] <- NA_real_
    pred_high[i] <- NA_real_
    next
  }

  lm_slope <- lm(slope_per_year ~ age_t0 + SC_t0 + cog_base + sex + mean_fd_t0 + mean_fd_t1, data = df)
  sm <- summary(lm_slope)
  beta_cog <- sm$coefficients["cog_base", "Estimate"]
  t_cog <- sm$coefficients["cog_base", "t value"]
  p_cog <- sm$coefficients["cog_base", "Pr(>|t|)"]

  res_rows[[i]] <- data.frame(
    edge = edge,
    n_sub = nrow(df),
    beta_cog = as.numeric(beta_cog),
    t_cog = as.numeric(t_cog),
    p_cog = as.numeric(p_cog),
    stringsAsFactors = FALSE
  )

  q10 <- quantile(df$cog_base, 0.1, na.rm = TRUE)
  q90 <- quantile(df$cog_base, 0.9, na.rm = TRUE)
  pred <- predict(lm_slope, newdata = df)
  pred_low[i] <- mean(pred[df$cog_base <= q10], na.rm = TRUE)
  pred_high[i] <- mean(pred[df$cog_base >= q90], na.rm = TRUE)
}

res_df <- do.call(rbind, res_rows)
res_df$p_cog_fdr <- p.adjust(res_df$p_cog, method = "fdr")
saveRDS(res_df, file.path(resultFolder, paste0("lgcm_slope_results_", Cogvar_base, "_CV", CVthr, ".rds")))
write.csv(res_df, file.path(resultFolder, paste0("lgcm_slope_results_", Cogvar_base, "_CV", CVthr, ".csv")), row.names = FALSE)

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
  file.path(resultFolder, paste0("lgcm_slope_decile_low_high_", Cogvar_base, "_CV", CVthr, ".csv")),
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
  geom_col(aes(alpha = group, linetype = group), color = "black", width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colorid) +
  scale_alpha_manual(values = c(high = 1, low = 0.35)) +
  scale_linetype_manual(values = c(high = "solid", low = "dashed")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "S-A decile", y = "Predicted slope per year")

ggsave(
  file.path(FigureFolder, paste0("bar_lgcm_slope_deciles_low_high_", Cogvar_base, "_CV", CVthr, ".tiff")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)
ggsave(
  file.path(FigureFolder, paste0("bar_lgcm_slope_deciles_low_high_", Cogvar_base, "_CV", CVthr, ".pdf")),
  bar_fig, width = 18, height = 12, units = "cm", bg = "transparent"
)

message("[INFO] Done.")
