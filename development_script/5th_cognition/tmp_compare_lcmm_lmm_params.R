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
SCdata$age <- age_to_years(SCdata$age)
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
  REML = FALSE
)

lcmm_re_age <- lcmm_fit[["predRE"]][["age"]]
lmm_u <- lmm_fit@u
n_sub <- length(unique(SCdata$subID))
if (length(lmm_u) != n_sub * 2) {
  warning("Unexpected lmm @u length: ", length(lmm_u), " (n_sub=", n_sub, ")")
}

# For (1 + age || subID), u is concatenated (intercept, age) by subject.
lmm_re_age_u <- tail(lmm_u, n_sub)
lmm_re_age_ranef <- as.numeric(lme4::ranef(lmm_fit)[["subID"]][["age"]])

min_len <- min(length(lcmm_re_age), length(lmm_re_age_u), length(lmm_re_age_ranef))
cmp_df <- data.frame(
  lcmm_re_age = lcmm_re_age[seq_len(min_len)],
  lmm_re_age_u = lmm_re_age_u[seq_len(min_len)],
  lmm_re_age_ranef = lmm_re_age_ranef[seq_len(min_len)],
  stringsAsFactors = FALSE
)

corr_u <- suppressWarnings(cor(cmp_df$lcmm_re_age, cmp_df$lmm_re_age_u, use = "pairwise.complete.obs"))
corr_ranef <- suppressWarnings(cor(cmp_df$lcmm_re_age, cmp_df$lmm_re_age_ranef, use = "pairwise.complete.obs"))
max_diff_u <- max(abs(cmp_df$lcmm_re_age - cmp_df$lmm_re_age_u), na.rm = TRUE)
max_diff_ranef <- max(abs(cmp_df$lcmm_re_age - cmp_df$lmm_re_age_ranef), na.rm = TRUE)

summary_df <- data.frame(
  edge = edge_name,
  n_sub = n_sub,
  len_lcmm = length(lcmm_re_age),
  len_lmm_u = length(lmm_u),
  corr_lcmm_vs_lmm_u = corr_u,
  corr_lcmm_vs_lmm_ranef = corr_ranef,
  max_abs_diff_lcmm_vs_lmm_u = max_diff_u,
  max_abs_diff_lcmm_vs_lmm_ranef = max_diff_ranef,
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
    lmm_re_age_u = lmm_re_age_u,
    lmm_re_age_ranef = lmm_re_age_ranef,
    compare = cmp_df
  ),
  paste0(out_base, "_detail.rds")
)

message("[INFO] Done. Summary saved: ", paste0(out_base, "_summary.csv"))
