## Rscript version of:
##   development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_HCPD.Rmd
##
## Purpose:
## - Compute age-resolved alignment between SC change rates (posterior derivatives)
##   and the S-A connectional axis rank (Spearman rho), across posterior draws.
## - Generate figures and print numeric results to stdout (captured by SLURM logs).
##
## Default inputs (project-relative):
## - outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/derivative.posterior.df.SA12_CV75.rds
## - outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV75/derivative.df78_CV75.rds
##
## Outputs:
## - outputs/results/4th_changerate_SAcorr/hcpd/combat_gam/CV75/alignment_posterior_corr_draw_age.csv.gz
## - outputs/results/4th_changerate_SAcorr/hcpd/combat_gam/CV75/alignment_summary_age_curve.csv
## - outputs/results/4th_changerate_SAcorr/hcpd/combat_gam/CV75/alignment_summary_flip_age.csv
## - outputs/figures/4th_changerate_SAcorr/hcpd/combat_gam/CV75/Alignment_development/*.tiff + *.pdf
##
## Usage:
##   Rscript --vanilla development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_HCPD.R
##   Rscript --vanilla .../S1st_SAcorr_alongAge_HCPD.R --cvthr=75 --n_draws=100 --force=1

rm(list = ls())

library(parallel)
library(ggplot2)

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

result_dir <- file.path(project_root, "outputs", "results", "4th_changerate_SAcorr", dataset, "combat_gam", paste0("CV", CVthr))
figure_dir <- file.path(project_root, "outputs", "figures", "4th_changerate_SAcorr", dataset, "combat_gam", paste0("CV", CVthr), "Alignment_development")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

in_derivative_posterior <- if (!is.null(args$derivative_posterior_rds)) {
  args$derivative_posterior_rds
} else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "combat_gam", paste0("CV", CVthr),
            paste0("derivative.posterior.df.SA", ds.resolution, "_CV", CVthr, ".rds"))
}
in_derivative <- if (!is.null(args$derivative_rds)) {
  args$derivative_rds
} else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "combat_gam", paste0("CV", CVthr),
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
message("[INFO] CVthr=", CVthr)
message("[INFO] in_derivative_posterior=", in_derivative_posterior)
message("[INFO] in_derivative=", in_derivative)
message("[INFO] n_cores=", n_cores)
message("[INFO] n_draws_req=", n_draws_req)

posterior_list <- readRDS(in_derivative_posterior)
if (!is.list(posterior_list) || length(posterior_list) == 0) stop("Invalid posterior derivative list: ", in_derivative_posterior)

edge_sarank_value <- vapply(posterior_list, function(df) {
  if (!is.data.frame(df)) return(NA_real_)
  if (!"SCrank" %in% names(df)) return(NA_real_)
  unique(as.numeric(df$SCrank))[1]
}, numeric(1))
if (any(is.na(edge_sarank_value))) {
  stop("Some posterior derivative entries missing numeric 'SCrank'. Upstream derivative script should add it.")
}
sa_rank <- rank(edge_sarank_value, ties.method = "average")

age_values <- sort(unique(as.numeric(posterior_list[[1]]$age)))
if (length(age_values) < 10) stop("Too few unique age points in posterior derivatives: ", length(age_values))

all_draws <- sort(unique(as.character(posterior_list[[1]]$draw)))
if (length(all_draws) == 0) stop("No 'draw' labels found in posterior derivatives.")

draw_keep <- all_draws[seq_len(min(length(all_draws), n_draws_req))]
message("[INFO] age_points=", length(age_values), " draw_count=", length(draw_keep), " edges=", length(posterior_list))

compute_corr_one_draw <- function(draw_label) {
  # Build edge x age matrix for this draw
  mat <- matrix(NA_real_, nrow = length(posterior_list), ncol = length(age_values))
  for (i in seq_along(posterior_list)) {
    df <- posterior_list[[i]]
    df <- df[df$draw == draw_label, , drop = FALSE]
    if (nrow(df) == 0) next
    # Align by age
    idx <- match(age_values, as.numeric(df$age))
    mat[i, ] <- as.numeric(df$posterior.derivative)[idx]
  }
  # Spearman rho at each age across edges
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

# Match the historical Rmd behavior: round correlations to 4 decimals before
# computing medians/CI and flip-age (rho≈0) window.
corr_mat_round <- round(corr_mat, 4)

# Save full draw x age correlations (compressed CSV)
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
  n_edges = length(posterior_list),
  n_draws = nrow(corr_mat),
  n_age_points = length(age_values),
  flip_age_median = flip_median,
  flip_age_ci_low = flip_ci[[1]],
  flip_age_ci_high = flip_ci[[2]]
)
write.csv(flip_df, out_flip, row.names = FALSE)

# Print numeric results to stdout for SLURM logs
min_idx <- which.min(median_rho)
max_idx <- which.max(median_rho)
message(sprintf("[RESULT] n_edges=%d n_draws=%d n_age_points=%d", length(posterior_list), nrow(corr_mat), length(age_values)))
message(sprintf("[RESULT] flip_age_median=%.5f CI95=[%.5f, %.5f]", flip_median, flip_ci[[1]], flip_ci[[2]]))
message(sprintf("[RESULT] rho_median_min=%.5f at age=%.5f", median_rho[[min_idx]], age_values[[min_idx]]))
message(sprintf("[RESULT] rho_median_max=%.5f at age=%.5f", median_rho[[max_idx]], age_values[[max_idx]]))
message(sprintf("[RESULT] rho_median_first=%.5f at age=%.5f", median_rho[[1]], age_values[[1]]))
message(sprintf("[RESULT] rho_median_last=%.5f at age=%.5f", median_rho[[length(median_rho)]], age_values[[length(age_values)]]))

# Figures: alignment curve + flip-age histogram
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
       title = sprintf("HCP-D (ComBat-GAM) CV%.0f | flip age median=%.2f", CVthr, flip_median)) +
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
  geom_vline(xintercept = flip_median, color = "red", linewidth = 0.9) +
  theme_classic() +
  labs(x = "Flip age (years)", y = "Count", title = "Flip age distribution across posterior draws") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 11, hjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_flip_age_hist.tiff"), p_hist, dpi = 600, width = 18, height = 12, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_flip_age_hist.pdf"), p_hist, dpi = 600, width = 18, height = 12, units = "cm", bg = "transparent")
