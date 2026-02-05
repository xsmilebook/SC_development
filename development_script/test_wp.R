suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(parallel)
})

rm(list = ls())

CVthr <- 75

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

functionFolder <- file.path(project_root, "gamfunction")
resultFolder <- file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", "abcd", "age_wp_bp_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "abcd", "age_wp_bp_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(resultFolder, paste0("lmm_agewp_bp_SC_CV", CVthr, ".rds"))
out_csv <- sub("\\.rds$", ".csv", out_rds)

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

SCdata <- readRDS(input_rds)
needed <- c("subID", "age", "sex", "mean_fd")
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)

message(
  "[INFO] SCdata age range (years): ",
  round(min(SCdata$age, na.rm = TRUE), 3), "–", round(max(SCdata$age, na.rm = TRUE), 3)
)

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

edge <- "SC.2_h"
df <- SCdata[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", edge)]
lmm_model <- lme4::lmer(SC.2_h ~ age_wp + age_bp + sex + mean_fd + (1 + age_wp || subID),
                        data = df, REML = FALSE)


