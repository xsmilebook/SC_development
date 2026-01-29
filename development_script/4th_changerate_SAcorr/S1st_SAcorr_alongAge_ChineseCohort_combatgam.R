## Rscript version of:
##   development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_ChineseCohort.Rmd
##
## Purpose:
## - Compute age-resolved alignment between SC change rates (posterior derivatives)
##   and the S-A connectional axis rank (Spearman rho), across posterior draws.
## - Generate figures and print numeric results to stdout (captured by SLURM logs).
##
## Default inputs (project-relative):
## - outputs/results/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/derivative.posterior.df.SA12_CV75.rds
## - outputs/results/2nd_fitdevelopmentalmodel/chinese/combat_gam/CV75/derivative.df78_CV75.rds
##
## Outputs:
## - outputs/results/4th_changerate_SAcorr/chinese/combat_gam/CV75/alignment_posterior_corr_draw_age.csv.gz
## - outputs/results/4th_changerate_SAcorr/chinese/combat_gam/CV75/alignment_summary_age_curve.csv
## - outputs/results/4th_changerate_SAcorr/chinese/combat_gam/CV75/alignment_summary_flip_age.csv
## - outputs/figures/4th_changerate_SAcorr/chinese/combat_gam/CV75/Alignment_development/*.tiff + *.pdf
##
## Usage:
##   Rscript --vanilla development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_ChineseCohort_combatgam.R
##   Rscript --vanilla .../S1st_SAcorr_alongAge_ChineseCohort_combatgam.R --cvthr=75 --n_draws=100 --force=1

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

is_windows <- .Platform$OS.type == "windows"
plot_only <- is_windows

dataset <- "chinese"
CVthr <- as_num(args$cvthr, 75)
ds.resolution <- 12
save_svg <- is_windows

in_euclid <- if (!is.null(args$euclid_csv)) {
  args$euclid_csv
} else {
  # Prefer the generic SA12 distance table if present.
  # (We keep distances under existing wd/* folders; Chinese does not have its own copy.)
  candidates <- c(
    file.path(project_root, "wd", "interdataFolder_ABCD", "average_EuclideanDistance_12.csv"),
    file.path(project_root, "wd", "interdataFolder_HCPD", "average_EuclideanDistance_12.csv")
  )
  candidates[file.exists(candidates)][[1]]
}

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

out_corr_gz_cd <- file.path(result_dir, "alignment_posterior_corr_draw_age_control_distance.csv.gz")
out_curve_cd <- file.path(result_dir, "alignment_summary_age_curve_control_distance.csv")
out_flip_cd <- file.path(result_dir, "alignment_summary_flip_age_control_distance.csv")

if (plot_only) {
  message("[INFO] Windows detected: plot-only mode (use existing results).")
  if (!file.exists(out_corr_gz) || !file.exists(out_curve) || !file.exists(out_flip) ||
      !file.exists(out_corr_gz_cd) || !file.exists(out_curve_cd) || !file.exists(out_flip_cd)) {
    stop("Plot-only mode requires existing outputs: ", result_dir)
  }
} else if (!force &&
           file.exists(out_corr_gz) && file.exists(out_curve) && file.exists(out_flip) &&
           file.exists(out_corr_gz_cd) && file.exists(out_curve_cd) && file.exists(out_flip_cd)) {
  message("[INFO] Outputs exist; skipping (set --force=1 to re-run): ", result_dir)
  quit(save = "no", status = 0)
}

message("[INFO] dataset=", dataset)
message("[INFO] CVthr=", CVthr)
message("[INFO] in_derivative_posterior=", in_derivative_posterior)
message("[INFO] in_derivative=", in_derivative)
message("[INFO] in_euclid=", in_euclid)
message("[INFO] n_cores=", n_cores)
message("[INFO] n_draws_req=", n_draws_req)
message("[INFO] save_svg=", save_svg)

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

if (!file.exists(in_euclid)) stop("Missing euclidean distance CSV (--euclid_csv): ", in_euclid)
euclid_df <- read.csv(in_euclid, stringsAsFactors = FALSE)
if (!all(c("Edistance", "SC_label") %in% names(euclid_df))) {
  stop("Euclidean distance CSV missing required columns Edistance/SC_label: ", in_euclid)
}
edge_labels <- paste0("SC.", seq_len(length(posterior_list)), "_h")
dist_match <- match(edge_labels, as.character(euclid_df$SC_label))
if (any(is.na(dist_match))) {
  missing_labs <- edge_labels[is.na(dist_match)]
  stop("Euclidean distance CSV missing labels (first few): ",
       paste(head(missing_labs, 10), collapse = ", "),
       " ... file=", in_euclid)
}
dist_vec <- as.numeric(euclid_df$Edistance[dist_match])
if (any(!is.finite(dist_vec))) stop("Non-finite Edistance found after matching labels in: ", in_euclid)

age_values <- sort(unique(as.numeric(posterior_list[[1]]$age)))
if (length(age_values) < 10) stop("Too few unique age points in posterior derivatives: ", length(age_values))

all_draws <- sort(unique(as.character(posterior_list[[1]]$draw)))
if (length(all_draws) == 0) stop("No 'draw' labels found in posterior derivatives.")

draw_keep <- all_draws[seq_len(min(length(all_draws), n_draws_req))]
message("[INFO] age_points=", length(age_values), " draw_count=", length(draw_keep), " edges=", length(posterior_list))

if (plot_only) {
  curve_df <- read.csv(out_curve)
  curve_df_cd <- read.csv(out_curve_cd)
  flip_df <- read.csv(out_flip)
  flip_df_cd <- read.csv(out_flip_cd)
  corr_df <- read.csv(gzfile(out_corr_gz))
  corr_df_cd <- read.csv(gzfile(out_corr_gz_cd))
  age_values <- as.numeric(sub("^age_", "", names(corr_df)[-1]))
  corr_mat <- as.matrix(corr_df[, -1, drop = FALSE])
  corr_mat_cd <- as.matrix(corr_df_cd[, -1, drop = FALSE])
  corr_mat_round <- round(corr_mat, 4)
  corr_mat_round_cd <- round(corr_mat_cd, 4)
  sum_main <- list(
    median_rho = as.numeric(curve_df$rho_median),
    ci_low = as.numeric(curve_df$rho_ci_low),
    ci_high = as.numeric(curve_df$rho_ci_high),
    flip_age_each_draw = vapply(seq_len(nrow(corr_mat_round)), function(i) {
      row <- corr_mat_round[i, ]
      idx <- which.min(abs(row - 0))
      age_values[[idx]]
    }, numeric(1)),
    flip_median = as.numeric(flip_df$flip_age_median[[1]]),
    flip_ci = c(as.numeric(flip_df$flip_age_ci_low[[1]]), as.numeric(flip_df$flip_age_ci_high[[1]]))
  )
  sum_cd <- list(
    median_rho = as.numeric(curve_df_cd$rho_median),
    ci_low = as.numeric(curve_df_cd$rho_ci_low),
    ci_high = as.numeric(curve_df_cd$rho_ci_high),
    flip_age_each_draw = vapply(seq_len(nrow(corr_mat_round_cd)), function(i) {
      row <- corr_mat_round_cd[i, ]
      idx <- which.min(abs(row - 0))
      age_values[[idx]]
    }, numeric(1)),
    flip_median = as.numeric(flip_df_cd$flip_age_median[[1]]),
    flip_ci = c(as.numeric(flip_df_cd$flip_age_ci_low[[1]]), as.numeric(flip_df_cd$flip_age_ci_high[[1]]))
  )
} else {
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
    rho_cd <- rep(NA_real_, ncol(mat))
    for (j in seq_len(ncol(mat))) {
      v <- mat[, j]
      if (all(is.na(v))) next
      rho[[j]] <- suppressWarnings(cor(v, sa_rank, method = "spearman", use = "pairwise.complete.obs"))

      ok <- is.finite(v) & is.finite(dist_vec) & is.finite(sa_rank)
      if (sum(ok) < 3) next
      resid <- suppressWarnings(residuals(stats::lm(v[ok] ~ dist_vec[ok])))
      rho_cd[[j]] <- suppressWarnings(cor(resid, sa_rank[ok], method = "spearman", use = "pairwise.complete.obs"))
    }
    list(rho = rho, rho_control_distance = rho_cd)
  }

  corr_rows <- mclapply(draw_keep, compute_corr_one_draw, mc.cores = n_cores)
  corr_mat <- do.call(rbind, lapply(corr_rows, `[[`, "rho"))
  corr_mat_cd <- do.call(rbind, lapply(corr_rows, `[[`, "rho_control_distance"))
  colnames(corr_mat) <- paste0("age_", format(age_values, trim = TRUE, scientific = FALSE))
  rownames(corr_mat) <- draw_keep
  colnames(corr_mat_cd) <- colnames(corr_mat)
  rownames(corr_mat_cd) <- draw_keep

  # Match the historical Rmd behavior: round correlations to 4 decimals before
  # computing medians/CI and flip-age (rho≈0) window.
  corr_mat_round <- round(corr_mat, 4)
  corr_mat_round_cd <- round(corr_mat_cd, 4)

  write_corr_gz <- function(path, mat) {
    con <- gzfile(path, open = "wt")
    write.csv(cbind(draw = rownames(mat), as.data.frame(mat)), con, row.names = FALSE)
    close(con)
  }

  # Save full draw x age correlations (compressed CSV)
  write_corr_gz(out_corr_gz, corr_mat)
  write_corr_gz(out_corr_gz_cd, corr_mat_cd)

  summarize_corr <- function(mat_round) {
    median_rho <- apply(mat_round, 2, median, na.rm = TRUE)
    ci_low <- apply(mat_round, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
    ci_high <- apply(mat_round, 2, function(x) as.numeric(quantile(x, probs = 0.975, na.rm = TRUE)))

    flip_age_each_draw <- vapply(seq_len(nrow(mat_round)), function(i) {
      row <- mat_round[i, ]
      idx <- which.min(abs(row - 0))
      age_values[[idx]]
    }, numeric(1))
    flip_median <- median(flip_age_each_draw, na.rm = TRUE)
    flip_ci <- as.numeric(quantile(flip_age_each_draw, probs = c(0.025, 0.975), na.rm = TRUE))

    list(
      median_rho = median_rho,
      ci_low = ci_low,
      ci_high = ci_high,
      flip_age_each_draw = flip_age_each_draw,
      flip_median = flip_median,
      flip_ci = flip_ci
    )
  }

  sum_main <- summarize_corr(corr_mat_round)
  sum_cd <- summarize_corr(corr_mat_round_cd)

  write_summary <- function(curve_path, flip_path, sum_obj, n_edges, n_draws, n_age_points) {
    curve_df <- data.frame(
      age = age_values,
      rho_median = sum_obj$median_rho,
      rho_ci_low = sum_obj$ci_low,
      rho_ci_high = sum_obj$ci_high
    )
    write.csv(curve_df, curve_path, row.names = FALSE)

    flip_df <- data.frame(
      n_edges = n_edges,
      n_draws = n_draws,
      n_age_points = n_age_points,
      flip_age_median = sum_obj$flip_median,
      flip_age_ci_low = sum_obj$flip_ci[[1]],
      flip_age_ci_high = sum_obj$flip_ci[[2]]
    )
    write.csv(flip_df, flip_path, row.names = FALSE)
    curve_df
  }

  curve_df <- write_summary(out_curve, out_flip, sum_main, length(posterior_list), nrow(corr_mat), length(age_values))
  curve_df_cd <- write_summary(out_curve_cd, out_flip_cd, sum_cd, length(posterior_list), nrow(corr_mat_cd), length(age_values))
}

# Print numeric results to stdout for SLURM logs
min_idx <- which.min(sum_main$median_rho)
max_idx <- which.max(sum_main$median_rho)
message(sprintf("[RESULT] n_edges=%d n_draws=%d n_age_points=%d", length(posterior_list), nrow(corr_mat), length(age_values)))
message(sprintf("[RESULT] flip_age_median=%.5f CI95=[%.5f, %.5f]", sum_main$flip_median, sum_main$flip_ci[[1]], sum_main$flip_ci[[2]]))
message(sprintf("[RESULT] rho_median_min=%.5f at age=%.5f", sum_main$median_rho[[min_idx]], age_values[[min_idx]]))
message(sprintf("[RESULT] rho_median_max=%.5f at age=%.5f", sum_main$median_rho[[max_idx]], age_values[[max_idx]]))
message(sprintf("[RESULT] rho_median_first=%.5f at age=%.5f", sum_main$median_rho[[1]], age_values[[1]]))
message(sprintf("[RESULT] rho_median_last=%.5f at age=%.5f", sum_main$median_rho[[length(sum_main$median_rho)]], age_values[[length(age_values)]]))

min_idx_cd <- which.min(sum_cd$median_rho)
max_idx_cd <- which.max(sum_cd$median_rho)
message(sprintf("[RESULT] control_distance flip_age_median=%.5f CI95=[%.5f, %.5f]", sum_cd$flip_median, sum_cd$flip_ci[[1]], sum_cd$flip_ci[[2]]))
message(sprintf("[RESULT] control_distance rho_median_min=%.5f at age=%.5f", sum_cd$median_rho[[min_idx_cd]], age_values[[min_idx_cd]]))
message(sprintf("[RESULT] control_distance rho_median_max=%.5f at age=%.5f", sum_cd$median_rho[[max_idx_cd]], age_values[[max_idx_cd]]))

# Figures: alignment curve + flip-age histogram
prep_curve_for_plot <- function(curve_df, flip_ci) {
  curve_df$in_flip_window <- (curve_df$age > flip_ci[[1]] & curve_df$age < flip_ci[[2]])
  curve_df$flip_window_x <- ifelse(curve_df$in_flip_window, curve_df$age, NA_real_)

  loess_m <- loess(rho_median ~ age, data = curve_df, span = 0.2)
  loess_l <- loess(rho_ci_low ~ age, data = curve_df, span = 0.2)
  loess_u <- loess(rho_ci_high ~ age, data = curve_df, span = 0.2)
  curve_df$rho_median_loess <- loess_m$fitted
  curve_df$rho_ci_low_loess <- loess_l$fitted
  curve_df$rho_ci_high_loess <- loess_u$fitted
  curve_df
}

mytheme <- theme(
  axis.text = element_text(size = 27, color = "black"),
  axis.title = element_text(size = 27),
  axis.line = element_line(linewidth = 0.6),
  axis.ticks = element_line(linewidth = 0.6),
  aspect.ratio = 1,
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.title = element_text(size = 15),
  legend.text = element_text(size = 15)
)

make_alignment_plot <- function(curve_df) {
  ggplot(curve_df) +
    geom_ribbon(aes(x = age, ymin = rho_ci_low_loess, ymax = rho_ci_high_loess), alpha = 0.3) +
    geom_line(aes(x = age, y = rho_median_loess), linewidth = 1.5) +
    geom_ribbon(aes(x = flip_window_x, ymin = rho_ci_low_loess, ymax = rho_ci_high_loess), fill = "#F8B01B", alpha = 1) +
    geom_ribbon(aes(x = flip_window_x, ymin = rho_median_loess - 0.01, ymax = rho_median_loess + 0.01), fill = "#B2182B", alpha = 1) +
    theme_classic() +
    labs(x = "Age (years)", y = "rho") +
    mytheme
}

make_flip_hist_plot <- function(flip_age_each_draw, flip_median) {
  ggplot() +
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
}

curve_df_plot <- prep_curve_for_plot(curve_df, sum_main$flip_ci)
p_align <- make_alignment_plot(curve_df_plot)

ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_alignment.tiff"), p_align, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_alignment.pdf"), p_align, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")

ggsave(file.path(figure_dir, "SA12_posDeriv_divweight_corr.tiff"), p_align, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "SA12_posDeriv_divweight_corr.pdf"), p_align, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")

p_hist_rmd <- make_flip_hist_plot(sum_main$flip_age_each_draw, sum_main$flip_median)
ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_flip_age_hist.tiff"), p_hist_rmd, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "SA12_posDeriv_SAaxis_flip_age_hist.pdf"), p_hist_rmd, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "Agedistribution_0corr.tiff"), p_hist_rmd, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "Agedistribution_0corr.pdf"), p_hist_rmd, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")

curve_df_plot_cd <- prep_curve_for_plot(curve_df_cd, sum_cd$flip_ci)
p_align_cd <- make_alignment_plot(curve_df_plot_cd)
ggsave(file.path(figure_dir, "SA12_posDeriv_divweight_corr_control_distance.tiff"), p_align_cd, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "SA12_posDeriv_divweight_corr_control_distance.pdf"), p_align_cd, dpi = 600, width = 13.5, height = 13.5, units = "cm", bg = "transparent")

p_hist_cd <- make_flip_hist_plot(sum_cd$flip_age_each_draw, sum_cd$flip_median)
ggsave(file.path(figure_dir, "Agedistribution_0corr_control_distance.tiff"), p_hist_cd, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")
ggsave(file.path(figure_dir, "Agedistribution_0corr_control_distance.pdf"), p_hist_cd, dpi = 600, width = 6, height = 6, units = "cm", bg = "transparent")

save_svg_if_available <- function(filename, plot, width_cm, height_cm, dpi = 600) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("[INFO] svglite not available; skip svg: ", filename)
    return(invisible(FALSE))
  }
  ggsave(filename, plot, dpi = dpi, width = width_cm, height = height_cm, units = "cm", bg = "transparent")
  invisible(TRUE)
}

if (save_svg) {
  save_svg_if_available(file.path(figure_dir, "SA12_posDeriv_divweight_corr.svg"), p_align, width_cm = 13.5, height_cm = 13.5)
  save_svg_if_available(file.path(figure_dir, "SA12_posDeriv_SAaxis_alignment.svg"), p_align, width_cm = 13.5, height_cm = 13.5)
  save_svg_if_available(file.path(figure_dir, "SA12_posDeriv_SAaxis_flip_age_hist.svg"), p_hist_rmd, width_cm = 6, height_cm = 6)
  save_svg_if_available(file.path(figure_dir, "Agedistribution_0corr.svg"), p_hist_rmd, width_cm = 6, height_cm = 6)
  save_svg_if_available(file.path(figure_dir, "SA12_posDeriv_divweight_corr_control_distance.svg"), p_align_cd, width_cm = 13.5, height_cm = 13.5)
  save_svg_if_available(file.path(figure_dir, "SA12_posDeriv_SAaxis_alignment_control_distance.svg"), p_align_cd, width_cm = 13.5, height_cm = 13.5)
  save_svg_if_available(file.path(figure_dir, "Agedistribution_0corr_control_distance.svg"), p_hist_cd, width_cm = 6, height_cm = 6)
}

# Additional figures from the original Rmd: age-specific derivative scatter/matrix and derivative line plot.
derivative_df <- readRDS(in_derivative)
if (!is.data.frame(derivative_df)) stop("Invalid derivative.df RDS: ", in_derivative)
if (!all(c("age", "label_ID", "derivative") %in% names(derivative_df))) {
  stop("derivative.df missing required columns age/label_ID/derivative: ", in_derivative)
}

age_grid <- sort(unique(as.numeric(derivative_df$age)))
age_min <- min(age_grid)
age_max <- max(age_grid)
age_mid <- age_grid[[which.min(abs(age_grid - sum_main$flip_median))]]
message(sprintf("[INFO] age_min=%.5f age_mid(nearest_flip)=%.5f age_max=%.5f", age_min, age_mid, age_max))

build_age_deriv_df <- function(target_age) {
  parcel_all <- paste0("SC.", seq_len(length(posterior_list)), "_h")
  df_out <- data.frame(label_ID = parcel_all, derivative = NA_real_)
  for (i in seq_along(parcel_all)) {
    lab <- parcel_all[[i]]
    sub <- derivative_df[derivative_df$label_ID == lab, , drop = FALSE]
    if (nrow(sub) == 0) next
    idx <- which.min(abs(as.numeric(sub$age) - target_age))
    df_out$derivative[[i]] <- as.numeric(sub$derivative[[idx]])
  }
  df_out
}

df_age_min <- build_age_deriv_df(age_min)
df_age_mid <- build_age_deriv_df(age_mid)
df_age_max <- build_age_deriv_df(age_max)

spearman_r <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
message(sprintf("[RESULT] age_min rho(derivative,SArank)=%.5f", spearman_r(df_age_min$derivative, sa_rank)))
message(sprintf("[RESULT] age_mid rho(derivative,SArank)=%.5f", spearman_r(df_age_mid$derivative, sa_rank)))
message(sprintf("[RESULT] age_max rho(derivative,SArank)=%.5f", spearman_r(df_age_max$derivative, sa_rank)))

# 3-line scatter plot
df_scatter <- rbind(
  data.frame(age = sprintf("%.2f", age_min), SCrank = sa_rank, derivative = df_age_min$derivative),
  data.frame(age = sprintf("%.2f", age_mid), SCrank = sa_rank, derivative = df_age_mid$derivative),
  data.frame(age = sprintf("%.2f", age_max), SCrank = sa_rank, derivative = df_age_max$derivative)
)
df_scatter$age <- as.factor(df_scatter$age)

p_scatter <- ggplot(df_scatter) +
  geom_smooth(aes(x = SCrank, y = derivative, group = age), color = "black", method = "lm", se = TRUE) +
  labs(x = "S-A connectional axis rank", y = "SC change rate") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15.4, color = "black"),
    axis.title = element_text(size = 15.4, color = "black"),
    aspect.ratio = 0.82,
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
ggsave(file.path(figure_dir, "deri.diw_corr_SCrank12ageAll.tiff"), p_scatter, dpi = 600, width = 14, height = 12, units = "cm", bg = "transparent")
if (save_svg) save_svg_if_available(file.path(figure_dir, "deri.diw_corr_SCrank12ageAll.svg"), p_scatter, width_cm = 12, height_cm = 10)
ggsave(file.path(figure_dir, "deri.diw_corr_SCrank12ageAll.pdf"), p_scatter, dpi = 600, width = 14, height = 12, units = "cm", bg = "transparent")

make_matrix_plot <- function(values_vec, title, out_tiff) {
  matsize <- ds.resolution
  mat <- matrix(NA_real_, nrow = matsize, ncol = matsize)
  idx_lower <- lower.tri(mat, diag = TRUE)
  idx_upper <- upper.tri(mat)
  mat[idx_lower] <- values_vec
  mat[idx_upper] <- t(mat)[idx_upper]
  colnames(mat) <- seq_len(matsize)
  rownames(mat) <- seq_len(matsize)
  mat_df <- as.data.frame(mat)
  mat_df$nodeid <- seq_len(matsize)
  m <- melt_df(mat_df, id.vars = c("nodeid"))
  m$variable <- as.numeric(m$variable)
  m$nodeid <- 0 - m$nodeid

  linerange_frame <- data.frame(
    x = c(0.5, matsize + 0.5), ymin = rep(-matsize - 0.5, times = 2), ymax = rep(-0.5, times = 2),
    y = c(-0.5, -matsize - 0.5), xmin = rep(0.5, times = 2), xmax = rep(matsize + 0.5, times = 2)
  )

  ggplot(m) +
    geom_tile(aes(x = variable, y = nodeid, fill = value, color = value)) +
    scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(-maxderiv, maxderiv)) +
    scale_color_distiller(type = "seq", palette = "RdBu", limits = c(-maxderiv, maxderiv)) +
    geom_segment(aes(x = 0.5, y = -0.5, xend = matsize + 0.5, yend = -0.5 - matsize), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(y = y, xmin = xmin, xmax = xmax), color = "black", linewidth = 0.5) +
    geom_linerange(data = linerange_frame, aes(x = x, ymin = ymin, ymax = ymax), color = "black", linewidth = 0.5) +
    ggtitle(label = title) +
    labs(x = "", y = "") +
    scale_y_continuous(breaks = NULL, labels = NULL) +
    scale_x_continuous(breaks = NULL, labels = NULL) +
    theme(
      axis.line = element_line(linewidth = 0),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid.major = element_line(linewidth = 0),
      panel.grid.minor = element_line(linewidth = 1)
    )
}

maxderiv <- max(abs(c(df_age_min$derivative, df_age_mid$derivative, df_age_max$derivative)), na.rm = TRUE)

p_mat_min <- make_matrix_plot(df_age_min$derivative, sprintf("Age = %.2f", age_min), "Deri_SA12_diw_age8.tiff")
ggsave(file.path(figure_dir, "Deri_SA12_diw_age8.tiff"), p_mat_min, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")
p_mat_mid <- make_matrix_plot(df_age_mid$derivative, sprintf("Age = %.2f", age_mid), "Deri_SA12_diw_age15.tiff")
ggsave(file.path(figure_dir, "Deri_SA12_diw_age15.tiff"), p_mat_mid, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")
p_mat_max <- make_matrix_plot(df_age_max$derivative, sprintf("Age = %.2f", age_max), "Deri_SA12_diw_age21.tiff")
ggsave(file.path(figure_dir, "Deri_SA12_diw_age21.tiff"), p_mat_max, dpi = 600, height = 13, width = 15, units = "cm", bg = "transparent")

# Derivative line plot across age
sc_rank_first <- rank(edge_sarank_value, ties.method = "first")
rank_df <- data.frame(label_ID = paste0("SC.", seq_len(length(posterior_list)), "_h"), SCrank12 = sc_rank_first)
derivative_merge <- merge(derivative_df, rank_df, by = "label_ID", all.x = TRUE)
p_line <- ggplot(derivative_merge) +
  geom_line(aes(x = age, y = derivative, group = label_ID, color = SCrank12), size = 1.4, alpha = 1) +
  scale_color_distiller(type = "seq", palette = "RdBu") +
  scale_fill_distiller(type = "seq", palette = "RdBu") +
  labs(x = "Age (years)", y = "SC change rate", color = "Axis rank") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 21.2, color = "black"),
    axis.title = element_text(size = 21.2, color = "black"),
    aspect.ratio = 0.65,
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )
ggsave(file.path(figure_dir, "derivative_diw_SA12_changerate.tiff"), p_line, dpi = 600, width = 20, height = 16, units = "cm", bg = "transparent")
if (save_svg) save_svg_if_available(file.path(figure_dir, "derivative_diw_SA12_changerate.svg"), p_line, width_cm = 20, height_cm = 16)
ggsave(file.path(figure_dir, "derivative_diw_SA12_changerate.pdf"), p_line, dpi = 600, width = 20, height = 16, units = "cm", bg = "transparent")
