#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(pbkrtest)
  library(parallel)
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
resultFolder <- file.path(project_root, "outputs", "results", "5th_cognition", "abcd", "withinperson_lmm")
FigureFolder <- file.path(project_root, "outputs", "figures", "5th_cognition", "abcd", "withinperson_lmm")
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
)
if (!file.exists(input_rds)) {
  stop("Missing input_rds: ", input_rds, "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (cognition variant)")
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

if (!Cogvar %in% names(SCdata)) {
  stop("Missing cognition variable in input: ", Cogvar, "\nInput: ", input_rds)
}
Cogdf <- SCdata %>%
  select(subID, eventname, all_of(Cogvar)) %>%
  filter(!is.na(.data[[Cogvar]])) %>%
  filter(grepl("base", eventname, ignore.case = TRUE)) %>%
  select(subID, !!Cogvar_base := all_of(Cogvar)) %>%
  distinct()
SCdata <- SCdata %>% left_join(Cogdf, by = "subID")
if (!Cogvar_base %in% names(SCdata)) stop("Baseline cognition join failed, missing: ", Cogvar_base)

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
res_cog <- residuals(lm(cognition_base ~ sex + mean_fd_mean, data = delta_df))
plot_df <- data.frame(res_delta = res_delta, res_cog = res_cog)

Fig <- ggplot(plot_df, aes(x = res_cog, y = res_delta)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "black") +
  theme_classic() +
  labs(
    x = "Baseline cognition (residualized: sex + mean_fd)",
    y = "Within-person Δ totalstrength / year (residualized: sex + mean_fd)"
  )

ggsave(file.path(FigureFolder, paste0("delta_totalstrength_vs_", Cogvar_base, "_residualized.pdf")), Fig, width = 12, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(FigureFolder, paste0("delta_totalstrength_vs_", Cogvar_base, "_residualized.tiff")), Fig, width = 12, height = 10, units = "cm", bg = "transparent", dpi = 600)

message("[INFO] Done.")
