## Check k values for GAM smooth using ComBat-GAM data.
##
## Usage:
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/V_check_k/V1st_check_k.R \
##     --project_root=/path/to/SCDevelopment \
##     --dataset=chinese \
##     --k_values=3,4,5,6
##
## Outputs (project-relative defaults):
## - outputs/results/2nd_fitdevelopmentalmodel/<dataset>/check_k/aic_by_edge.csv
## - outputs/results/2nd_fitdevelopmentalmodel/<dataset>/check_k/aic_summary.csv
## - outputs/figures/2nd_fitdevelopmentalmodel/<dataset>/check_k/aic_compare_by_k.{tiff,pdf}

rm(list = ls())

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

parse_args <- function(args) {
  res <- list()
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    res[[kv[[1]]]] <- kv[[2]]
  }
  res
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}

dataset <- if (!is.null(args$dataset)) tolower(args$dataset) else "chinese"
if (!dataset %in% c("hcpd", "abcd", "chinese")) {
  stop("Unsupported dataset: ", dataset, " (expected hcpd|abcd|chinese)")
}

k_values <- if (!is.null(args$k_values)) {
  as.integer(strsplit(args$k_values, ",", fixed = TRUE)[[1]])
} else {
  c(3L, 4L, 5L, 6L)
}
if (any(is.na(k_values)) || any(k_values < 3)) {
  stop("Invalid k_values: ", args$k_values)
}

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else if (dataset == "hcpd") {
  file.path(project_root, "outputs", "results", "combat_gam", "hcpd",
            "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds")
} else if (dataset == "abcd") {
  file.path(project_root, "outputs", "results", "combat_gam", "abcd",
            "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds")
} else {
  file.path(project_root, "outputs", "results", "combat_gam", "chinese",
            "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds")
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)

output_dir <- if (!is.null(args$output_dir)) {
  args$output_dir
} else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "check_k")
}
figure_dir <- if (!is.null(args$figure_dir)) {
  args$figure_dir
} else {
  file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", dataset, "check_k")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

SCdata <- readRDS(input_rds)
if (!is.data.frame(SCdata)) stop("Expected a data.frame in input_rds: ", input_rds)

# Normalize covariate column names
if (!"age" %in% names(SCdata) && "Age" %in% names(SCdata)) SCdata$age <- SCdata$Age
if (!"sex" %in% names(SCdata) && "Sex" %in% names(SCdata)) SCdata$sex <- SCdata$Sex
if (!"mean_fd" %in% names(SCdata) && "meanFD" %in% names(SCdata)) SCdata$mean_fd <- SCdata$meanFD

required <- c("age", "sex", "mean_fd")
missing_required <- setdiff(required, names(SCdata))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)
SCdata$mean_fd <- as.numeric(SCdata$mean_fd)

sc_cols_all <- grep("^SC\\.", names(SCdata), value = TRUE)
if (length(sc_cols_all) == 0) stop("No SC.* columns found in input_rds: ", input_rds)
sc_cols <- sc_cols_all
sc_cols_h <- sc_cols_all[grepl("_h$", sc_cols_all)]
if (length(sc_cols_h) > 0) sc_cols <- sc_cols_h

n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else length(sc_cols))
n_edges <- min(length(sc_cols), max(1L, n_edges))
sc_cols <- sc_cols[seq_len(n_edges)]

SCdata <- SCdata[, c(required, sc_cols), drop = FALSE]
SCdata <- SCdata[stats::complete.cases(SCdata), , drop = FALSE]
if (nrow(SCdata) < 50) stop("Too few complete rows after filtering: ", nrow(SCdata))

message("[INFO] dataset=", dataset)
message("[INFO] input_rds=", input_rds)
message("[INFO] k_values=", paste(k_values, collapse = ","))
message("[INFO] n_edges=", n_edges, " n_cores=", n_cores)
message("[INFO] output_dir=", output_dir)
message("[INFO] figure_dir=", figure_dir)

fit_aic_one_edge <- function(edge) {
  out <- data.frame(region = edge, k = k_values, AIC = NA_real_)
  for (i in seq_along(k_values)) {
    k <- k_values[[i]]
    modelformula <- as.formula(sprintf("%s ~ s(age, k = %d, fx = TRUE) + sex + mean_fd", edge, k))
    aic_val <- tryCatch({
      model <- mgcv::gam(modelformula, method = "REML", data = SCdata)
      stats::AIC(model)
    }, error = function(e) NA_real_)
    out$AIC[[i]] <- aic_val
  }
  out
}

aic_list <- mclapply(sc_cols, fit_aic_one_edge, mc.cores = n_cores)
aic_df <- dplyr::bind_rows(aic_list)

if (!any(is.finite(aic_df$AIC))) {
  stop("All AIC values are NA; check model inputs.")
}

# Delta AIC per edge (relative to best k)
aic_df <- aic_df %>%
  group_by(region) %>%
  mutate(AIC_min = min(AIC, na.rm = TRUE), delta_AIC = AIC - AIC_min) %>%
  ungroup()

best_k_df <- aic_df %>%
  group_by(region) %>%
  filter(AIC == min(AIC, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  select(region, best_k = k, best_AIC = AIC)

write.csv(aic_df, file.path(output_dir, "aic_by_edge.csv"), row.names = FALSE)
write.csv(best_k_df, file.path(output_dir, "aic_summary.csv"), row.names = FALSE)

# AIC comparison figure
plot_df <- aic_df %>% filter(is.finite(delta_AIC))
plot_df$k <- factor(plot_df$k, levels = sort(unique(plot_df$k)))

p <- ggplot(plot_df, aes(x = k, y = delta_AIC)) +
  geom_boxplot(outlier.shape = NA, fill = "#D9E6F2") +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.2, color = "#2C3E50") +
  labs(x = "k", y = "Delta AIC (vs. best)", title = paste0(toupper(dataset), " ComBat-GAM")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )

ggsave(file.path(figure_dir, "aic_compare_by_k.tiff"), p, dpi = 600, width = 14, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "aic_compare_by_k.pdf"), p, dpi = 600, width = 14, height = 10, units = "cm", bg = "transparent")
