## HCP-D (ComBat-GAM SA12) | 4th: Age-resolved S-A alignment with added covariates
##
## This script computes the age-resolved alignment between posterior derivatives
## (SC change rates) and S-A connectional axis rank (Spearman rho), across draws.
##
## Two recommended variants:
## - SES: cov_tag=SES (derivatives from 2nd devmodel covariate pipeline with income.adj)
## - ICV: cov_tag=ICV (derivatives from 2nd devmodel covariate pipeline with ICV)
##
## Outputs:
## - outputs/results/4th_changerate_SAcorr/hcpd/covariates/<cov_tag>/combat_gam/CV75/*
## - outputs/figures/4th_changerate_SAcorr/hcpd/covariates/<cov_tag>/combat_gam/CV75/Alignment_development/*.tiff + *.pdf
##
## Notes:
## - Figures are tiff+pdf only (no svg).
## - Numeric results are printed as [RESULT] lines for SLURM logs.

rm(list = ls())

library(parallel)
library(ggplot2)

melt_df <- function(x, id.vars) {
  if (requireNamespace("reshape2", quietly = TRUE)) {
    return(reshape2::melt(x, id.vars = id.vars))
  }
  if (requireNamespace("reshape", quietly = TRUE)) {
    return(reshape::melt(x, id.vars = id.vars))
  }
  stop("Missing package for melt(): install reshape2 or reshape in the runtime environment.")
}

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

as_int <- function(x, default) {
  if (is.null(x) || is.na(suppressWarnings(as.integer(x)))) return(as.integer(default))
  as.integer(x)
}

as_num <- function(x, default) {
  if (is.null(x) || is.na(suppressWarnings(as.numeric(x)))) return(as.numeric(default))
  as.numeric(x)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}

dataset <- "hcpd"
CVthr <- as_num(args$cvthr, 75)
ds.resolution <- 12
cov_tag <- if (!is.null(args$cov_tag)) args$cov_tag else "SES"

result_dir <- file.path(project_root, "outputs", "results", "4th_changerate_SAcorr", dataset, "covariates", cov_tag, "combat_gam", paste0("CV", CVthr))
figure_dir <- file.path(project_root, "outputs", "figures", "4th_changerate_SAcorr", dataset, "covariates", cov_tag, "combat_gam", paste0("CV", CVthr), "Alignment_development")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

in_derivative_posterior <- if (!is.null(args$derivative_posterior_rds)) {
  args$derivative_posterior_rds
} else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "covariates", cov_tag, "combat_gam", paste0("CV", CVthr),
            paste0("derivative.posterior.df.SA", ds.resolution, "_CV", CVthr, ".rds"))
}
in_derivative <- if (!is.null(args$derivative_rds)) {
  args$derivative_rds
} else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "covariates", cov_tag, "combat_gam", paste0("CV", CVthr),
            paste0("derivative.df", ds.resolution * (ds.resolution + 1) / 2, "_CV", CVthr, ".rds"))
}

if (!file.exists(in_derivative_posterior)) stop("Missing derivative_posterior_rds: ", in_derivative_posterior)
if (!file.exists(in_derivative)) stop("Missing derivative_rds: ", in_derivative)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

n_draws_req <- as_int(args$n_draws, 1000L)
force <- as_int(args$force, 0L) == 1L

out_corr_gz <- file.path(result_dir, "alignment_posterior_corr_draw_age.csv.gz")
out_curve <- file.path(result_dir, "alignment_summary_age_curve.csv")
out_flip <- file.path(result_dir, "alignment_summary_flip_age.csv")

if (!force && file.exists(out_corr_gz) && file.exists(out_curve) && file.exists(out_flip)) {
  message("[INFO] Outputs exist; skipping (set --force=1 to re-run): ", result_dir)
  quit(save = "no", status = 0)
}

message("[INFO] dataset=", dataset)
message("[INFO] cov_tag=", cov_tag)
message("[INFO] CVthr=", CVthr)
message("[INFO] in_derivative_posterior=", in_derivative_posterior)
message("[INFO] in_derivative=", in_derivative)
message("[INFO] n_cores=", n_cores)
message("[INFO] n_draws_req=", n_draws_req)

posterior_list <- readRDS(in_derivative_posterior)
if (!is.list(posterior_list) || length(posterior_list) == 0) stop("Invalid posterior derivative list: ", in_derivative_posterior)

edge_sarank_value <- vapply(posterior_list, function(df) {
  if (!is.data.frame(df) || !"SCrank" %in% names(df)) return(NA_real_)
  unique(as.numeric(df$SCrank))[1]
}, numeric(1))
if (anyNA(edge_sarank_value)) {
  stop("Missing SCrank in posterior derivative list; check 2nd devmodel derivative generation outputs: ", in_derivative_posterior)
}
sa_rank <- edge_sarank_value

draw_labels_all <- unique(posterior_list[[1]]$draw)
draw_labels_all <- as.character(draw_labels_all)
draw_labels_all <- draw_labels_all[order(draw_labels_all)]
if (length(draw_labels_all) == 0) stop("No draw labels in posterior derivatives: ", in_derivative_posterior)

if (n_draws_req > length(draw_labels_all)) {
  message("[WARN] Requested n_draws=", n_draws_req, " but only ", length(draw_labels_all), " available; using all.")
  n_draws_req <- length(draw_labels_all)
}
draw_keep <- draw_labels_all[seq_len(n_draws_req)]

age_values <- unique(as.numeric(posterior_list[[1]]$age))
age_values <- age_values[order(age_values)]
if (length(age_values) == 0) stop("No age values in posterior derivatives: ", in_derivative_posterior)

compute_corr_one_draw <- function(draw_label) {
  mat <- matrix(NA_real_, nrow = length(posterior_list), ncol = length(age_values))
  for (i in seq_along(posterior_list)) {
    df <- posterior_list[[i]]
    df <- df[df$draw == draw_label, , drop = FALSE]
    if (nrow(df) == 0) next
    idx <- match(age_values, as.numeric(df$age))
    mat[i, ] <- as.numeric(df$posterior.derivative)[idx]
  }
  rho <- rep(NA_real_, ncol(mat))
  for (j in seq_len(ncol(mat))) {
    v <- mat[, j]
    if (all(is.na(v))) next
    rho[[j]] <- suppressWarnings(cor(v, sa_rank, method = "spearman", use = "pairwise.complete.obs"))
  }
  rho
}

corr_rows <- mclapply(draw_keep, compute_corr_one_draw, mc.cores = n_cores)
corr_mat <- do.call(rbind, corr_rows)
colnames(corr_mat) <- paste0("age_", format(age_values, trim = TRUE, scientific = FALSE))
rownames(corr_mat) <- draw_keep

corr_mat_round <- round(corr_mat, 4)

con <- gzfile(out_corr_gz, open = "wt")
on.exit(close(con), add = TRUE)
write.csv(cbind(draw = rownames(corr_mat), as.data.frame(corr_mat)), con, row.names = FALSE)

median_rho <- apply(corr_mat_round, 2, median, na.rm = TRUE)
ci_low <- apply(corr_mat_round, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
ci_high <- apply(corr_mat_round, 2, function(x) as.numeric(quantile(x, probs = 0.975, na.rm = TRUE)))

flip_age_each_draw <- vapply(seq_len(nrow(corr_mat)), function(i) {
  row <- corr_mat_round[i, ]
  idx <- which.min(abs(row - 0))
  age_values[[idx]]
}, numeric(1))
flip_median <- median(flip_age_each_draw, na.rm = TRUE)
flip_ci <- as.numeric(quantile(flip_age_each_draw, probs = c(0.025, 0.975), na.rm = TRUE))

curve_df <- data.frame(
  age = age_values,
  rho_median = median_rho,
  rho_ci_low = ci_low,
  rho_ci_high = ci_high
)
write.csv(curve_df, out_curve, row.names = FALSE)

flip_df <- data.frame(
  cov_tag = cov_tag,
  n_edges = length(posterior_list),
  n_draws = nrow(corr_mat),
  n_age_points = length(age_values),
  flip_age_median = flip_median,
  flip_age_ci_low = flip_ci[[1]],
  flip_age_ci_high = flip_ci[[2]]
)
write.csv(flip_df, out_flip, row.names = FALSE)

min_idx <- which.min(median_rho)
max_idx <- which.max(median_rho)
message(sprintf("[RESULT] cov_tag=%s n_edges=%d n_draws=%d n_age_points=%d", cov_tag, length(posterior_list), nrow(corr_mat), length(age_values)))
message(sprintf("[RESULT] flip_age_median=%.5f CI95=[%.5f, %.5f]", flip_median, flip_ci[[1]], flip_ci[[2]]))
message(sprintf("[RESULT] rho_median_min=%.5f at age=%.5f", median_rho[[min_idx]], age_values[[min_idx]]))
message(sprintf("[RESULT] rho_median_max=%.5f at age=%.5f", median_rho[[max_idx]], age_values[[max_idx]]))
message(sprintf("[RESULT] rho_median_first=%.5f at age=%.5f", median_rho[[1]], age_values[[1]]))
message(sprintf("[RESULT] rho_median_last=%.5f at age=%.5f", median_rho[[length(median_rho)]], age_values[[length(age_values)]]))

curve_df$in_flip_window <- (curve_df$age > flip_ci[[1]] & curve_df$age < flip_ci[[2]])
curve_df$flip_window_x <- ifelse(curve_df$in_flip_window, curve_df$age, NA_real_)

loess_m <- loess(rho_median ~ age, data = curve_df, span = 0.2)
loess_l <- loess(rho_ci_low ~ age, data = curve_df, span = 0.2)
loess_u <- loess(rho_ci_high ~ age, data = curve_df, span = 0.2)
curve_df$rho_median_loess <- loess_m$fitted
curve_df$rho_ci_low_loess <- loess_l$fitted
curve_df$rho_ci_high_loess <- loess_u$fitted

p_align <- ggplot(curve_df) +
  geom_ribbon(aes(x = age, ymin = rho_ci_low_loess, ymax = rho_ci_high_loess), alpha = 0.25) +
  geom_line(aes(x = age, y = rho_median_loess), linewidth = 1.2) +
  geom_ribbon(aes(x = flip_window_x, ymin = rho_ci_low_loess, ymax = rho_ci_high_loess), fill = "#F8B01B", alpha = 0.9) +
  theme_classic() +
  labs(x = "Age (years)", y = "Alignment with S-A axis (Spearman rho)",
       title = sprintf("HCP-D (ComBat-GAM) %s | flip age median=%.2f", cov_tag, flip_median)) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 11, hjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_alignment.tiff"), p_align, dpi = 600, width = 18, height = 12, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_alignment.pdf"), p_align, dpi = 600, width = 18, height = 12, units = "cm", bg = "transparent")

p_hist <- ggplot(data.frame(flip_age = flip_age_each_draw), aes(x = flip_age)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "white", linewidth = 0.6) +
  geom_vline(xintercept = flip_median, color = "red", linewidth = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Count", title = "Flip-age distribution across draws") +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 11, hjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
ggsave(file.path(figure_dir, "Flip_age_histogram.tiff"), p_hist, dpi = 600, width = 14, height = 10, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "Flip_age_histogram.pdf"), p_hist, dpi = 600, width = 14, height = 10, units = "cm", bg = "transparent")

