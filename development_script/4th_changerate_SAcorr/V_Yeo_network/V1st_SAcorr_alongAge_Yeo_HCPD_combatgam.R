## HCP-D (ComBat-GAM) | Yeo7/Yeo17 | Age-resolved alignment with S-A axis
##
## Purpose:
## - Compute age-resolved alignment between SC change rates (posterior derivatives)
##   and the S-A connectional axis rank (Spearman rho) across posterior draws.
## - Generate figures and print key numeric results to stdout (captured by SLURM logs).
##
## Inputs (project-relative defaults):
## - outputs/results/2nd_fitdevelopmentalmodel/hcpd/yeo/Yeo{7,17}/combat_gam/CV75/derivative.posterior.df.Yeo{7,17}_CV75.rds
##
## Outputs:
## - outputs/results/4th_changerate_SAcorr/hcpd/yeo/Yeo{7,17}/combat_gam/CV75/alignment_*.csv*
## - outputs/figures/4th_changerate_SAcorr/hcpd/yeo/Yeo{7,17}/combat_gam/CV75/Alignment_development/*.tiff + *.pdf
##
## Usage:
##   Rscript --vanilla development_script/4th_changerate_SAcorr/V_Yeo_network/V1st_SAcorr_alongAge_Yeo_HCPD_combatgam.R --yeo=7
##   Rscript --vanilla .../V1st_SAcorr_alongAge_Yeo_HCPD_combatgam.R --yeo=17 --n_draws=200 --force=1

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

infer_resolution_from_edges <- function(n_edges) {
  n <- (sqrt(8 * n_edges + 1) - 1) / 2
  n_int <- as.integer(round(n))
  if (n_int < 1 || n_int * (n_int + 1) / 2 != n_edges) {
    stop("Cannot infer resolution from n_edges=", n_edges, " (expected triangular number).")
  }
  n_int
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}

dataset <- "hcpd"
CVthr <- as_num(args$cvthr, 75)
yeo <- as_int(args$yeo, 17L)
if (!yeo %in% c(7L, 17L)) stop("Unsupported yeo resolution: ", yeo, " (expected 7 or 17)")

n_draws_req <- as_int(args$n_draws, 1000L)
force <- as_int(args$force, 0L) == 1L

result_dir <- file.path(project_root, "outputs", "results", "4th_changerate_SAcorr", dataset, "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr))
figure_dir <- file.path(project_root, "outputs", "figures", "4th_changerate_SAcorr", dataset, "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr), "Alignment_development")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

in_derivative_posterior <- if (!is.null(args$derivative_posterior_rds)) {
  args$derivative_posterior_rds
} else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "yeo", paste0("Yeo", yeo), "combat_gam", paste0("CV", CVthr),
            paste0("derivative.posterior.df.Yeo", yeo, "_CV", CVthr, ".rds"))
}
if (!file.exists(in_derivative_posterior)) stop("Missing derivative_posterior_rds: ", in_derivative_posterior)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

out_corr_gz <- file.path(result_dir, "alignment_posterior_corr_draw_age.csv.gz")
out_curve <- file.path(result_dir, "alignment_summary_age_curve.csv")
out_flip <- file.path(result_dir, "alignment_summary_flip_age.csv")

if (!force && file.exists(out_corr_gz) && file.exists(out_curve) && file.exists(out_flip)) {
  message("[INFO] Outputs exist; skipping. Set --force=1 to re-run.")
  quit(save = "no", status = 0)
}

posterior_list <- readRDS(in_derivative_posterior)
if (!is.list(posterior_list) || length(posterior_list) == 0) stop("Expected a non-empty list in: ", in_derivative_posterior)

n_edges <- length(posterior_list)
matsize <- infer_resolution_from_edges(n_edges)

age_values <- sort(unique(as.numeric(posterior_list[[1]]$age)))
if (length(age_values) < 10) stop("Too few unique age points in posterior derivatives: ", length(age_values))

all_draws <- sort(unique(as.character(posterior_list[[1]]$draw)))
if (length(all_draws) == 0) stop("No 'draw' labels found in posterior derivatives.")
draw_keep <- all_draws[seq_len(min(length(all_draws), n_draws_req))]

# S-A axis rank for the inferred matrix size
edge_sa_value <- matrix(NA_real_, matsize, matsize)
for (x in 1:matsize) {
  for (y in 1:matsize) {
    edge_sa_value[x, y] <- x^2 + y^2
  }
}
edge_sa_value <- edge_sa_value[lower.tri(edge_sa_value, diag = TRUE)]
sa_rank <- rank(edge_sa_value, ties.method = "average")

message("[INFO] dataset=", dataset, " yeo=", yeo, " matsize=", matsize, " edges=", n_edges)
message("[INFO] age_points=", length(age_values), " draw_count=", length(draw_keep), " n_cores=", n_cores)

compute_corr_one_draw <- function(draw_label) {
  mat <- matrix(NA_real_, nrow = n_edges, ncol = length(age_values))
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

# Match historical behavior: round to 4 decimals before summary/flip-age.
corr_mat_round <- round(corr_mat, 4)

con <- gzfile(out_corr_gz, open = "wt")
on.exit(close(con), add = TRUE)
write.csv(cbind(draw = rownames(corr_mat), as.data.frame(corr_mat)), con, row.names = FALSE)

median_rho <- apply(corr_mat_round, 2, median, na.rm = TRUE)
ci_low <- apply(corr_mat_round, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
ci_high <- apply(corr_mat_round, 2, function(x) as.numeric(quantile(x, probs = 0.975, na.rm = TRUE)))

flip_age_each_draw <- vapply(seq_len(nrow(corr_mat_round)), function(i) {
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
  n_edges = n_edges,
  n_draws = nrow(corr_mat),
  n_age_points = length(age_values),
  flip_age_median = flip_median,
  flip_age_ci_low = flip_ci[[1]],
  flip_age_ci_high = flip_ci[[2]]
)
write.csv(flip_df, out_flip, row.names = FALSE)

min_idx <- which.min(median_rho)
max_idx <- which.max(median_rho)
message(sprintf("[RESULT] dataset=%s yeo=%d n_edges=%d n_draws=%d n_age_points=%d", dataset, yeo, n_edges, nrow(corr_mat), length(age_values)))
message(sprintf("[RESULT] flip_age_median=%.5f CI95=[%.5f, %.5f]", flip_median, flip_ci[[1]], flip_ci[[2]]))
message(sprintf("[RESULT] rho_median_min=%.5f at age=%.5f", median_rho[[min_idx]], age_values[[min_idx]]))
message(sprintf("[RESULT] rho_median_max=%.5f at age=%.5f", median_rho[[max_idx]], age_values[[max_idx]]))

curve_df$in_flip_window <- (curve_df$age > flip_ci[[1]] & curve_df$age < flip_ci[[2]])
curve_df$flip_window_x <- ifelse(curve_df$in_flip_window, curve_df$age, NA_real_)

loess_m <- loess(rho_median ~ age, data = curve_df, span = 0.2)
loess_l <- loess(rho_ci_low ~ age, data = curve_df, span = 0.2)
loess_u <- loess(rho_ci_high ~ age, data = curve_df, span = 0.2)
curve_df$rho_median_loess <- loess_m$fitted
curve_df$rho_ci_low_loess <- loess_l$fitted
curve_df$rho_ci_high_loess <- loess_u$fitted

p_align <- ggplot(curve_df) +
  geom_ribbon(aes(x = age, ymin = rho_ci_low_loess, ymax = rho_ci_high_loess), alpha = 0.3) +
  geom_line(aes(x = age, y = rho_median_loess), linewidth = 1.5) +
  geom_ribbon(aes(x = flip_window_x, ymin = rho_ci_low_loess, ymax = rho_ci_high_loess), fill = "#F8B01B", alpha = 1) +
  geom_ribbon(aes(x = flip_window_x, ymin = rho_median_loess - 0.01, ymax = rho_median_loess + 0.01), fill = "#B2182B", alpha = 1) +
  theme_classic() +
  labs(x = "Age (years)", y = "rho") +
  theme(
    axis.text = element_text(size = 23, color = "black"),
    axis.title = element_text(size = 23),
    aspect.ratio = 1,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )

ggsave(file.path(figure_dir, paste0("Yeo", yeo, "_posDeriv_divweight_corr.tiff")), p_align, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, paste0("Yeo", yeo, "_posDeriv_divweight_corr.pdf")), p_align, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")

p_hist <- ggplot() +
  geom_histogram(aes(flip_age_each_draw, y = ..count..), binwidth = 0.5, linewidth = 1, color = "black", fill = "white") +
  geom_vline(aes(xintercept = flip_median), colour = "red", linetype = "solid", linewidth = 1) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  theme(
    axis.line = element_blank(),
    aspect.ratio = 0.5,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 20)
  )

ggsave(file.path(figure_dir, paste0("Agedistribution_0corr_Yeo", yeo, ".tiff")), p_hist, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, paste0("Agedistribution_0corr_Yeo", yeo, ".pdf")), p_hist, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")

