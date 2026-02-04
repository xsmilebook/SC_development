#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(lcmm)
  library(dplyr)
})

rm(list = ls())

CVthr <- 75
edge_name <- Sys.getenv("EDGE_NAME", unset = "SC.1_h")

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "age_lmm", "tmp_lcmm_lmm_compare")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)

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
if (!file.exists(plotdatasum_rds)) {
  stop("Missing ABCD_PLOTDATASUM_RDS: ", plotdatasum_rds)
}
plotdata <- readRDS(plotdatasum_rds)

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
SCdata$age <- as.numeric(SCdata$age)
message(
  "[INFO] SCdata age range (years): ",
  round(min(SCdata$age, na.rm = TRUE), 3), "–", round(max(SCdata$age, na.rm = TRUE), 3)
)

needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$sex <- as.factor(SCdata$sex)

sc_cols <- grep("^SC\\.", names(SCdata), value = TRUE)
if (any(grepl("_h$", sc_cols))) sc_cols <- sc_cols[grepl("_h$", sc_cols)]
if (!edge_name %in% sc_cols) stop("EDGE_NAME not found in SCdata: ", edge_name)

# Scale SC strength by initial fit (ratio), same as S3.
plotdata.tmp <- plotdata[plotdata$SC_label == edge_name, , drop = FALSE]
if (!("fit" %in% names(plotdata.tmp)) || nrow(plotdata.tmp) < 1 || is.na(plotdata.tmp$fit[[1]]) || plotdata.tmp$fit[[1]] == 0) {
  stop("Missing/invalid plotdata fit for edge: ", edge_name)
}
SCdata[[edge_name]] <- SCdata[[edge_name]] / plotdata.tmp$fit[[1]]
SCdata[[edge_name]] <- as.numeric(SCdata[[edge_name]])

# Keep subjects with >=2 unique timepoints.
sub_time_ok <- tapply(SCdata$age, SCdata$subID, function(x) length(unique(round(x, 6))) >= 2)
SCdata <- SCdata[SCdata$subID %in% names(sub_time_ok)[sub_time_ok], , drop = FALSE]
if (nrow(SCdata) < 10) stop("Too few rows after >=2 timepoint filter: ", nrow(SCdata))

# lcmm ID encoding (as specified).
SCdata$ID <- factor(SCdata$subID, levels = unique(SCdata$subID), labels = seq_len(length(unique(SCdata$subID))))
SCdata$ID <- as.numeric(SCdata$ID)

message("[INFO] Fitting lcmm: ", edge_name)
lcmm_fit <- lcmm::lcmm(
  as.formula(sprintf("%s ~ age + sex + mean_fd", edge_name)),
  random = ~ age,
  subject = "ID",
  data = SCdata,
  ng = 1
)

message("[INFO] Fitting lmm: ", edge_name)
lmm_fit <- lme4::lmer(
  as.formula(sprintf("%s ~ age + sex + mean_fd + (1 + age || subID)", edge_name)),
  data = SCdata,
  REML = TRUE
)

lcmm_re_age <- lcmm_fit[["predRE"]][["age"]]
n_sub <- length(unique(SCdata$subID))
lmm_re_age_ranef <- as.numeric(lme4::ranef(lmm_fit)[["subID"]][["age"]])

# Align by subject ID (lcmm uses numeric ID).
id_map <- data.frame(
  subID = unique(SCdata$subID),
  ID = as.numeric(factor(unique(SCdata$subID), levels = unique(SCdata$subID))),
  stringsAsFactors = FALSE
)

id_order <- sort(unique(SCdata$ID))
if (length(id_order) != length(lcmm_re_age)) {
  warning(
    "Length mismatch between lcmm predRE$age (", length(lcmm_re_age),
    ") and unique(SCdata$ID) (", length(id_order), "). Will align by min length."
  )
  min_len <- min(length(id_order), length(lcmm_re_age))
  id_order <- id_order[seq_len(min_len)]
  lcmm_re_age <- lcmm_re_age[seq_len(min_len)]
}

lcmm_df <- data.frame(ID = id_order, lcmm_re_age = as.numeric(lcmm_re_age), stringsAsFactors = FALSE)
lcmm_df <- merge(lcmm_df, id_map, by = "ID", all.x = TRUE)

lmm_re_df <- data.frame(
  subID = rownames(lme4::ranef(lmm_fit)[["subID"]]),
  lmm_re_age_ranef = lmm_re_age_ranef,
  stringsAsFactors = FALSE
)

cmp_df <- merge(lcmm_df, lmm_re_df, by = "subID", all = FALSE)

# Fixed effects
lmm_beta_age <- as.numeric(lme4::fixef(lmm_fit)["age"])
lcmm_beta_age <- NA_real_
if (!is.null(names(lcmm_fit$best)) && "age" %in% names(lcmm_fit$best)) {
  lcmm_beta_age <- as.numeric(lcmm_fit$best[["age"]])
}

# Personal slopes = fixed + random
cmp_df$lmm_slope_personal <- lmm_beta_age + cmp_df$lmm_re_age_ranef
cmp_df$lcmm_slope_personal <- lcmm_beta_age + cmp_df$lcmm_re_age

corr_ranef <- suppressWarnings(cor(cmp_df$lcmm_re_age, cmp_df$lmm_re_age_ranef, use = "pairwise.complete.obs"))
corr_personal <- suppressWarnings(cor(cmp_df$lcmm_slope_personal, cmp_df$lmm_slope_personal, use = "pairwise.complete.obs"))
max_diff_ranef <- max(abs(cmp_df$lcmm_re_age - cmp_df$lmm_re_age_ranef), na.rm = TRUE)
max_diff_personal <- max(abs(cmp_df$lcmm_slope_personal - cmp_df$lmm_slope_personal), na.rm = TRUE)

# ---- lcmm scale transform (latent -> observed) ----
lin2 <- as.numeric(lcmm_fit$best["Linear 2 (std err)"])  # ≈ residual SD on observed scale

# fixed effect on observed scale
lcmm_beta_age_obs <- lcmm_beta_age * lin2

# random slope (predRE) on observed scale
cmp_df$lcmm_re_age_obs <- cmp_df$lcmm_re_age * lin2

# personal slope on observed scale
cmp_df$lcmm_slope_personal_obs <- (lcmm_beta_age + cmp_df$lcmm_re_age) * lin2

# lmer personal slope already on observed scale
cmp_df$lmm_slope_personal_obs <- lmm_beta_age + cmp_df$lmm_re_age_ranef

# now compare on same scale
corr_personal_obs <- suppressWarnings(cor(
  cmp_df$lcmm_slope_personal_obs, cmp_df$lmm_slope_personal_obs, use="pairwise.complete.obs"
))
max_diff_personal_obs <- max(abs(
  cmp_df$lcmm_slope_personal_obs - cmp_df$lmm_slope_personal_obs
), na.rm=TRUE)


summary_df <- data.frame(
  edge = edge_name,
  n_sub = n_sub,
  n_sub_overlap = nrow(cmp_df),
  len_lcmm = length(lcmm_re_age),
  beta_age_lmm = lmm_beta_age,
  beta_age_lcmm = lcmm_beta_age,
  corr_lcmm_vs_lmm_ranef = corr_ranef,
  corr_personal_slope = corr_personal,
  max_abs_diff_lcmm_vs_lmm_ranef = max_diff_ranef,
  max_abs_diff_personal_slope = max_diff_personal,
  stringsAsFactors = FALSE
)

out_base <- file.path(resultFolder, paste0("lcmm_lmm_compare_", edge_name, "_CV", CVthr))
write.csv(summary_df, paste0(out_base, "_summary.csv"), row.names = FALSE)
saveRDS(
  list(
    summary = summary_df,
    lcmm_fit = lcmm_fit,
    lmm_fit = lmm_fit,
    lcmm_re_age = lcmm_re_age,
    lmm_re_age_ranef = lmm_re_age_ranef,
    beta_age_lmm = lmm_beta_age,
    beta_age_lcmm = lcmm_beta_age,
    compare = cmp_df
  ),
  paste0(out_base, "_detail.rds")
)

message("[INFO] Done. Summary saved: ", paste0(out_base, "_summary.csv"))
