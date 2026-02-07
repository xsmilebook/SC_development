#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
})

rm(list = ls())

CVthr <- 75
Pvar <- "GENERAL"
Pvar_base <- "GENERAL_base"

project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("Please run from SCDevelopment project root (missing ARCHITECTURE.md): ", project_root)
}

resultFolder <- file.path(
  project_root, "outputs", "results", "6th_pfactor", "abcd",
  "age_wp_bp_lmm_baselineage_interaction_decileavg"
)
FigureFolder <- file.path(
  project_root, "outputs", "figures", "6th_pfactor", "abcd",
  "age_wp_bp_lmm_baselineage_interaction_decileavg"
)
dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(
  project_root, "outputs", "results", "combat_gam", "abcd",
  "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
)
if (!file.exists(input_rds)) {
  stop(
    "Missing input_rds: ", input_rds,
    "\nRun first: sbatch combat_gam/sbatch/abcd_combat_gam.sbatch (p-factor variant)"
  )
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

ensure_eventname <- function(df, data_name) {
  if (!("eventname" %in% names(df)) && ("scanID" %in% names(df))) {
    df$eventname <- scanid_to_eventname(df$scanID)
  }
  if (!("eventname" %in% names(df))) {
    stop("Missing eventname in ", data_name, " (required for baseline extraction)")
  }
  df
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

compute_baseline_age <- function(df, sub_col = "subID", age_col = "age", event_col = "eventname") {
  subs <- unique(as.character(df[[sub_col]]))
  out <- rep(NA_real_, length(subs))
  for (i in seq_along(subs)) {
    sid <- subs[[i]]
    idx <- which(as.character(df[[sub_col]]) == sid & !is.na(df[[age_col]]))
    if (length(idx) < 1) next
    age_i <- as.numeric(df[[age_col]][idx])
    event_i <- as.character(df[[event_col]][idx])
    idx_base <- which(grepl("base", event_i, ignore.case = TRUE))
    if (length(idx_base) > 0) {
      out[[i]] <- age_i[idx_base][which.min(age_i[idx_base])]
    } else {
      out[[i]] <- age_i[which.min(age_i)]
    }
  }
  data.frame(subID = subs, age_bp = out, stringsAsFactors = FALSE)
}

extract_baseline_cov <- function(df, value_col, out_col, sub_col = "subID", event_col = "eventname", age_col = "age") {
  if (!value_col %in% names(df)) stop("Missing covariate column: ", value_col)
  idx <- !is.na(df[[value_col]]) & grepl("base", as.character(df[[event_col]]), ignore.case = TRUE)
  base_df <- df[idx, , drop = FALSE]
  if (nrow(base_df) < 1) stop("No baseline rows found for: ", value_col)
  if (age_col %in% names(base_df)) {
    base_df[[age_col]] <- as.numeric(base_df[[age_col]])
    base_df <- base_df[order(as.character(base_df[[sub_col]]), base_df[[age_col]]), , drop = FALSE]
  } else {
    base_df <- base_df[order(as.character(base_df[[sub_col]])), , drop = FALSE]
  }
  base_df <- base_df[!duplicated(as.character(base_df[[sub_col]])), c(sub_col, value_col), drop = FALSE]
  names(base_df) <- c("subID", out_col)
  base_df
}

predict_low_high <- function(
    fit_obj, df_model, age_wp_vec, age_bp_vec, x_vec, x_name,
    sex_ref, mean_fd_ref, subid_ref, q_low, q_high
) {
  if (is.null(fit_obj)) return(data.frame())

  new_template <- data.frame(
    subID = subid_ref,
    age_wp = as.numeric(age_wp_vec[[1]]),
    age_bp = as.numeric(age_bp_vec[[1]]),
    sex = sex_ref,
    mean_fd = mean_fd_ref,
    cov = q_low,
    stringsAsFactors = FALSE
  )
  if (is.factor(df_model$sex)) {
    new_template$sex <- factor(sex_ref, levels = levels(df_model$sex))
  }
  if (is.factor(df_model$subID)) {
    new_template$subID <- factor(subid_ref, levels = levels(df_model$subID))
  }

  n_pred <- length(x_vec)
  new_low <- new_template[rep(1, n_pred), , drop = FALSE]
  new_low$age_wp <- as.numeric(age_wp_vec)
  new_low$age_bp <- as.numeric(age_bp_vec)
  pred_low <- tryCatch(
    stats::predict(fit_obj, newdata = new_low, re.form = NA, allow.new.levels = TRUE),
    error = function(e) rep(NA_real_, n_pred)
  )

  new_high <- new_low
  new_high$cov <- q_high
  pred_high <- tryCatch(
    stats::predict(fit_obj, newdata = new_high, re.form = NA, allow.new.levels = TRUE),
    error = function(e) rep(NA_real_, n_pred)
  )

  out_low <- data.frame(
    label = "low",
    .fitted = as.numeric(pred_low),
    stringsAsFactors = FALSE
  )
  out_high <- data.frame(
    label = "high",
    .fitted = as.numeric(pred_high),
    stringsAsFactors = FALSE
  )
  out_low[[x_name]] <- as.numeric(x_vec)
  out_high[[x_name]] <- as.numeric(x_vec)
  rbind(out_low, out_high)
}

fit_outcome_models <- function(
    data_all, y_col, cov_col, q_low, q_high, age_wp_seq, age_bp_seq,
    age_wp_ref, age_bp_ref, sex_ref, mean_fd_ref, subid_ref
) {
  df <- data_all[, c("subID", "age_wp", "age_bp", "sex", "mean_fd", cov_col, y_col)]
  names(df)[ncol(df)] <- "y"
  df <- df[complete.cases(df), , drop = FALSE]

  row_info <- data.frame(
    outcome = y_col,
    n_obs = nrow(df),
    beta_agewp_cov = NA_real_,
    beta_agebp_cov = NA_real_,
    p_agewp_cov = NA_real_,
    p_agebp_cov = NA_real_,
    stringsAsFactors = FALSE
  )
  if (nrow(df) < 10) {
    return(list(row = row_info, pred_wp = data.frame(), pred_bp = data.frame()))
  }

  df$cov <- df[[cov_col]]
  form_wp <- stats::as.formula("y ~ age_wp * cov + age_bp + sex + mean_fd + (1 | subID)")
  form_bp <- stats::as.formula("y ~ age_wp + age_bp * cov + sex + mean_fd + (1 | subID)")

  fit_wp <- tryCatch(lme4::lmer(form_wp, data = df, REML = FALSE), error = function(e) NULL)
  fit_bp <- tryCatch(lme4::lmer(form_bp, data = df, REML = FALSE), error = function(e) NULL)
  form_red <- stats::as.formula("y ~ age_wp + cov + age_bp + sex + mean_fd + (1 | subID)")
  fit_red <- tryCatch(lme4::lmer(form_red, data = df, REML = FALSE), error = function(e) NULL)

  get_coef <- function(fit_obj, coef_name) {
    if (is.null(fit_obj)) return(NA_real_)
    cf <- tryCatch(lme4::fixef(fit_obj), error = function(e) NULL)
    if (is.null(cf)) return(NA_real_)
    if (!coef_name %in% names(cf)) return(NA_real_)
    as.numeric(cf[[coef_name]])
  }
  row_info$beta_agewp_cov <- get_coef(fit_wp, "age_wp:cov")
  row_info$beta_agebp_cov <- get_coef(fit_bp, "age_bp:cov")
  get_lrt_p <- function(fit_red_obj, fit_full_obj) {
    if (is.null(fit_red_obj) || is.null(fit_full_obj)) return(NA_real_)
    tb <- tryCatch(stats::anova(fit_red_obj, fit_full_obj), error = function(e) NULL)
    if (is.null(tb) || nrow(tb) < 2 || !("Pr(>Chisq)" %in% names(tb))) return(NA_real_)
    as.numeric(tb$`Pr(>Chisq)`[2])
  }
  row_info$p_agewp_cov <- get_lrt_p(fit_red, fit_wp)
  row_info$p_agebp_cov <- get_lrt_p(fit_red, fit_bp)

  pred_wp <- predict_low_high(
    fit_wp, df_model = df,
    age_wp_vec = age_wp_seq,
    age_bp_vec = rep(age_bp_ref, length(age_wp_seq)),
    x_vec = age_wp_seq,
    x_name = "age_wp",
    sex_ref = sex_ref, mean_fd_ref = mean_fd_ref, subid_ref = subid_ref,
    q_low = q_low, q_high = q_high
  )
  pred_bp <- predict_low_high(
    fit_bp, df_model = df,
    age_wp_vec = rep(age_wp_ref, length(age_bp_seq)),
    age_bp_vec = age_bp_seq,
    x_vec = age_bp_seq,
    x_name = "age_bp",
    sex_ref = sex_ref, mean_fd_ref = mean_fd_ref, subid_ref = subid_ref,
    q_low = q_low, q_high = q_high
  )
  if (nrow(pred_wp) > 0) pred_wp$outcome <- y_col
  if (nrow(pred_bp) > 0) pred_bp$outcome <- y_col

  list(row = row_info, pred_wp = pred_wp, pred_bp = pred_bp)
}

plot_decile_curves <- function(
    plot_df, out_dir, out_prefix, line_types, x_col,
    p_col = NULL, p_df = NULL, y_lab = "SC strength (ratio)"
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  plot_df$label <- factor(plot_df$label, levels = c("low", "high"))
  colorid <- rev(RColorBrewer::brewer.pal(10, "RdBu"))
  x_values <- as.numeric(plot_df[[x_col]])

  x_min <- floor(min(x_values, na.rm = TRUE))
  x_max <- ceiling(max(x_values, na.rm = TRUE))
  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    x_breaks <- pretty(x_values, n = 5)
  } else {
    x_breaks <- seq(x_min, x_max, by = 1)
  }

  plots <- vector("list", length = 10)
  for (i in 1:10) {
    tmp <- plot_df[plot_df$decile == i, , drop = FALSE]
    if (nrow(tmp) < 1) {
      tmp <- data.frame(
        decile = i,
        label = factor(c("low", "high"), levels = c("low", "high")),
        fit.avg = NA_real_,
        stringsAsFactors = FALSE
      )
      tmp[[x_col]] <- mean(x_values, na.rm = TRUE)
    }
    color_i <- colorid[[i]]

    if (i == 1 || i == 6) {
      mytheme <- theme(
        axis.text = element_text(size = 21, color = "black"),
        axis.title = element_text(size = 21),
        aspect.ratio = 1,
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(fill = NA, color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none"
      )
    } else {
      mytheme <- theme(
        axis.text.x = element_text(size = 21, color = "black"),
        axis.text.y = element_text(size = 21, color = "transparent"),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21, colour = "transparent"),
        aspect.ratio = 1,
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5, colour = "transparent"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5, colour = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_line(linewidth = 0.5, colour = "transparent"),
        panel.border = element_rect(fill = NA, color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none"
      )
    }

    fig <- ggplot(tmp) +
      geom_line(
        aes(x = .data[[x_col]], y = fit.avg, group = label, linetype = label),
        linewidth = 1.2, color = color_i, na.rm = TRUE
      ) +
      scale_linetype_manual(values = line_types) +
      scale_x_continuous(
        breaks = x_breaks,
        labels = function(x) sprintf("%d", as.integer(round(x)))
      ) +
      scale_y_continuous(
        breaks = c(0.9, 1.0, 1.1),
        limits = c(0.85, 1.15),
        labels = function(y) sprintf("%.1f", y)
      ) +
      labs(x = "Age", y = y_lab) +
      mytheme
    if (!is.null(p_df) && !is.null(p_col) && p_col %in% names(p_df)) {
        p_row <- p_df[p_df$decile == i, , drop = FALSE]
        if (nrow(p_row) > 0) {
          p_val <- as.numeric(p_row[[p_col]][1])
          p_txt <- if (is.finite(p_val)) sprintf("FDR p=%.3g", p_val) else "FDR p=NA"
        x_pos <- suppressWarnings(min(tmp[[x_col]], na.rm = TRUE))
        if (!is.finite(x_pos)) x_pos <- suppressWarnings(mean(x_values, na.rm = TRUE))
        if (!is.finite(x_pos)) x_pos <- 0
        y_pos <- 1.145
        fig <- fig + annotate("text", x = x_pos, y = y_pos, label = p_txt, hjust = 0, vjust = 1, size = 6)
      }
    }

    out_base <- file.path(out_dir, paste0(out_prefix, "_decile", i))
    ggsave(paste0(out_base, ".tiff"), fig, width = 10, height = 10, units = "cm", bg = "transparent")
    ggsave(paste0(out_base, ".pdf"), fig, width = 10, height = 10, units = "cm", bg = "transparent")
    plots[[i]] <- fig
  }

  combined <- patchwork::wrap_plots(plots, nrow = 2, ncol = 5)
  combined_base <- file.path(out_dir, paste0(out_prefix, "_all_2x5"))
  ggsave(paste0(combined_base, ".tiff"), combined, width = 50, height = 20, units = "cm", bg = "transparent")
  ggsave(paste0(combined_base, ".pdf"), combined, width = 50, height = 20, units = "cm", bg = "transparent")
}

SCdata <- readRDS(input_rds)
SCdata <- ensure_eventname(SCdata, "SCdata")
needed <- c("subID", "age", "sex", "mean_fd", Pvar)
missing <- setdiff(needed, names(SCdata))
if (length(missing) > 0) stop("Missing required columns in SCdata: ", paste(missing, collapse = ", "))
SCdata$age <- as.numeric(SCdata$age)
SCdata$sex <- as.factor(SCdata$sex)

Pvardf <- extract_baseline_cov(SCdata, value_col = Pvar, out_col = Pvar_base)
SCdata <- SCdata %>% left_join(Pvardf, by = "subID")
SCdata <- SCdata[!is.na(SCdata[[Pvar_base]]), , drop = FALSE]

sub_n <- table(SCdata$subID)
keep_sub <- names(sub_n[sub_n >= 2])
SCdata <- SCdata[SCdata$subID %in% keep_sub, , drop = FALSE]

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

base_age_df <- compute_baseline_age(SCdata)
if (any(is.na(base_age_df$age_bp))) {
  stop("Baseline age extraction failed for subjects: ", paste(base_age_df$subID[is.na(base_age_df$age_bp)], collapse = ", "))
}
SCdata <- SCdata %>% left_join(base_age_df, by = "subID")
SCdata$age_wp <- SCdata$age - SCdata$age_bp

base_by_sub <- tapply(SCdata[[Pvar_base]], SCdata$subID, function(x) x[which(!is.na(x))[1]])
base_by_sub <- base_by_sub[!is.na(base_by_sub)]
q10 <- as.numeric(stats::quantile(base_by_sub, 0.1, na.rm = TRUE))
q90 <- as.numeric(stats::quantile(base_by_sub, 0.9, na.rm = TRUE))

age_wp_seq <- seq(min(SCdata$age_wp, na.rm = TRUE), max(SCdata$age_wp, na.rm = TRUE), length.out = 100)
age_bp_seq <- seq(min(SCdata$age_bp, na.rm = TRUE), max(SCdata$age_bp, na.rm = TRUE), length.out = 100)
age_wp_ref <- mean(SCdata$age_wp, na.rm = TRUE)
age_bp_ref <- mean(SCdata$age_bp, na.rm = TRUE)
mean_fd_ref <- mean(SCdata$mean_fd, na.rm = TRUE)
sex_tab <- table(SCdata$sex)
sex_ref <- names(sex_tab)[which.max(sex_tab)]
subid_ref <- as.character(SCdata$subID[1])

edge_map <- SA12_10[, c("SC_label", "decile"), drop = FALSE]
edge_map <- edge_map[!duplicated(edge_map$SC_label), , drop = FALSE]
edge_map$decile <- as.integer(edge_map$decile)
edge_map <- edge_map[!is.na(edge_map$decile) & edge_map$SC_label %in% sc_cols, , drop = FALSE]
if (!all(1:10 %in% unique(edge_map$decile))) {
  stop("SA12 decile mapping must contain all deciles 1..10 for current edges")
}

force <- as.integer(Sys.getenv("FORCE", unset = "0")) == 1
out_rds <- file.path(
  resultFolder,
  paste0("lmm_agewp_agebp_baselineage_interaction_decileavg_", Pvar_base, "_CV", CVthr, ".rds")
)

need_refit <- force
if (!need_refit && file.exists(out_rds)) {
  message("[INFO] Found existing results, loading (set FORCE=1 to recompute)")
  result_obj <- readRDS(out_rds)
  has_axis <- TRUE
  if (!("age_wp" %in% names(result_obj$edge_wp_decile))) has_axis <- FALSE
  if (!("age_bp" %in% names(result_obj$edge_bp_decile))) has_axis <- FALSE
  if (!("age_wp" %in% names(result_obj$decile_wp_plot))) has_axis <- FALSE
  if (!("age_bp" %in% names(result_obj$decile_bp_plot))) has_axis <- FALSE
  has_fdr <- all(c("p_agewp_cov_fdr", "p_agebp_cov_fdr") %in% names(result_obj$decile_summary))
  if (!has_axis || !has_fdr) {
    message("[INFO] Existing cache uses old format; recomputing with age_wp/age_bp x-axis and FDR results")
    need_refit <- TRUE
  } else {
    cache_wp_min <- min(result_obj$edge_wp_decile$age_wp, na.rm = TRUE)
    cache_wp_max <- max(result_obj$edge_wp_decile$age_wp, na.rm = TRUE)
    target_wp_min <- min(SCdata$age_wp, na.rm = TRUE)
    target_wp_max <- max(SCdata$age_wp, na.rm = TRUE)
    cache_bp_min <- min(result_obj$edge_bp_decile$age_bp, na.rm = TRUE)
    cache_bp_max <- max(result_obj$edge_bp_decile$age_bp, na.rm = TRUE)
    target_bp_min <- min(SCdata$age_bp, na.rm = TRUE)
    target_bp_max <- max(SCdata$age_bp, na.rm = TRUE)
    if (
      is.finite(cache_wp_min) && is.finite(cache_wp_max) && is.finite(target_wp_min) && is.finite(target_wp_max) &&
      is.finite(cache_bp_min) && is.finite(cache_bp_max) && is.finite(target_bp_min) && is.finite(target_bp_max)
    ) {
      if (
        cache_wp_min > target_wp_min + 1e-6 || cache_wp_max < target_wp_max - 1e-6 ||
        cache_bp_min > target_bp_min + 1e-6 || cache_bp_max < target_bp_max - 1e-6
      ) {
        message("[INFO] Existing cache does not cover full age_wp/age_bp span; recomputing")
        need_refit <- TRUE
      }
    }
  }
}

if (need_refit || !file.exists(out_rds)) {
  message("[INFO] Fitting edge-wise interaction models: age_wp*cov and age_bp*cov")
  edge_fit_list <- lapply(
    sc_cols,
    function(edge) {
      fit_outcome_models(
        data_all = SCdata, y_col = edge, cov_col = Pvar_base,
        q_low = q10, q_high = q90,
        age_wp_seq = age_wp_seq, age_bp_seq = age_bp_seq,
        age_wp_ref = age_wp_ref, age_bp_ref = age_bp_ref,
        sex_ref = sex_ref, mean_fd_ref = mean_fd_ref, subid_ref = subid_ref
      )
    }
  )

  edge_summary <- dplyr::bind_rows(lapply(edge_fit_list, `[[`, "row"))
  edge_pred_wp <- dplyr::bind_rows(lapply(edge_fit_list, `[[`, "pred_wp"))
  edge_pred_bp <- dplyr::bind_rows(lapply(edge_fit_list, `[[`, "pred_bp"))

  edge_wp_dec <- data.frame()
  edge_bp_dec <- data.frame()
  if (nrow(edge_pred_wp) > 0) {
    edge_pred_wp <- edge_pred_wp %>% rename(SC_label = outcome)
    edge_wp_dec <- edge_pred_wp %>%
      inner_join(edge_map, by = "SC_label") %>%
      group_by(decile, age_wp, label) %>%
      summarise(fit.avg = safe_mean(.fitted), .groups = "drop")
  }
  if (nrow(edge_pred_bp) > 0) {
    edge_pred_bp <- edge_pred_bp %>% rename(SC_label = outcome)
    edge_bp_dec <- edge_pred_bp %>%
      inner_join(edge_map, by = "SC_label") %>%
      group_by(decile, age_bp, label) %>%
      summarise(fit.avg = safe_mean(.fitted), .groups = "drop")
  }

  message("[INFO] Building decile-aggregated SC ratio (1..10) and fitting interaction models")
  SCdata_dec <- SCdata
  for (d in 1:10) {
    edges_d <- edge_map$SC_label[edge_map$decile == d]
    if (length(edges_d) < 1) stop("No edges found for decile ", d)
    SCdata_dec[[paste0("SC_decile", d)]] <- rowMeans(SCdata_dec[, edges_d, drop = FALSE], na.rm = TRUE)
  }

  decile_outcomes <- paste0("SC_decile", 1:10)
  decile_fit_list <- lapply(
    decile_outcomes,
    function(dec_col) {
      fit_outcome_models(
        data_all = SCdata_dec, y_col = dec_col, cov_col = Pvar_base,
        q_low = q10, q_high = q90,
        age_wp_seq = age_wp_seq, age_bp_seq = age_bp_seq,
        age_wp_ref = age_wp_ref, age_bp_ref = age_bp_ref,
        sex_ref = sex_ref, mean_fd_ref = mean_fd_ref, subid_ref = subid_ref
      )
    }
  )

  decile_summary <- dplyr::bind_rows(lapply(decile_fit_list, `[[`, "row"))
  decile_pred_wp <- dplyr::bind_rows(lapply(decile_fit_list, `[[`, "pred_wp"))
  decile_pred_bp <- dplyr::bind_rows(lapply(decile_fit_list, `[[`, "pred_bp"))

  if (nrow(decile_summary) > 0) {
    decile_summary$decile <- as.integer(sub("^SC_decile", "", decile_summary$outcome))
    decile_summary$p_agewp_cov_fdr <- stats::p.adjust(decile_summary$p_agewp_cov, method = "fdr")
    decile_summary$p_agebp_cov_fdr <- stats::p.adjust(decile_summary$p_agebp_cov, method = "fdr")
  } else {
    decile_summary$decile <- integer(0)
    decile_summary$p_agewp_cov_fdr <- numeric(0)
    decile_summary$p_agebp_cov_fdr <- numeric(0)
  }

  decile_wp_plot <- data.frame()
  decile_bp_plot <- data.frame()
  if (nrow(decile_pred_wp) > 0) {
    decile_pred_wp$decile <- as.integer(sub("^SC_decile", "", decile_pred_wp$outcome))
    decile_wp_plot <- decile_pred_wp %>%
      group_by(decile, age_wp, label) %>%
      summarise(fit.avg = safe_mean(.fitted), .groups = "drop")
  }
  if (nrow(decile_pred_bp) > 0) {
    decile_pred_bp$decile <- as.integer(sub("^SC_decile", "", decile_pred_bp$outcome))
    decile_bp_plot <- decile_pred_bp %>%
      group_by(decile, age_bp, label) %>%
      summarise(fit.avg = safe_mean(.fitted), .groups = "drop")
  }

  result_obj <- list(
    q10 = q10,
    q90 = q90,
    edge_summary = edge_summary,
    edge_wp_decile = edge_wp_dec,
    edge_bp_decile = edge_bp_dec,
    decile_summary = decile_summary,
    decile_wp_plot = decile_wp_plot,
    decile_bp_plot = decile_bp_plot
  )
  saveRDS(result_obj, out_rds)

  write.csv(
    edge_summary,
    file.path(resultFolder, paste0("edge_summary_baselineage_interaction_", Pvar_base, "_CV", CVthr, ".csv")),
    row.names = FALSE
  )
  write.csv(
    decile_summary,
    file.path(resultFolder, paste0("decile_avgSC_summary_baselineage_interaction_", Pvar_base, "_CV", CVthr, ".csv")),
    row.names = FALSE
  )
}

line_types <- c(low = "solid", high = "dashed")

if (nrow(result_obj$edge_wp_decile) > 0) {
  plot_decile_curves(
    result_obj$edge_wp_decile,
    out_dir = file.path(FigureFolder, "edge_first", "age_wp_interaction"),
    out_prefix = "developmentcurve_edgefirst_agewp",
    line_types = line_types,
    x_col = "age_wp"
  )
}
if (nrow(result_obj$edge_bp_decile) > 0) {
  plot_decile_curves(
    result_obj$edge_bp_decile,
    out_dir = file.path(FigureFolder, "edge_first", "age_bp_interaction"),
    out_prefix = "developmentcurve_edgefirst_agebp",
    line_types = line_types,
    x_col = "age_bp"
  )
}
if (nrow(result_obj$decile_wp_plot) > 0) {
  plot_decile_curves(
    result_obj$decile_wp_plot,
    out_dir = file.path(FigureFolder, "decile_avg_sc_first", "age_wp_interaction"),
    out_prefix = "developmentcurve_decileavgSC_agewp",
    line_types = line_types,
    x_col = "age_wp",
    p_col = "p_agewp_cov_fdr",
    p_df = result_obj$decile_summary
  )
}
if (nrow(result_obj$decile_bp_plot) > 0) {
  plot_decile_curves(
    result_obj$decile_bp_plot,
    out_dir = file.path(FigureFolder, "decile_avg_sc_first", "age_bp_interaction"),
    out_prefix = "developmentcurve_decileavgSC_agebp",
    line_types = line_types,
    x_col = "age_bp",
    p_col = "p_agebp_cov_fdr",
    p_df = result_obj$decile_summary
  )
}

message("[INFO] Done.")
