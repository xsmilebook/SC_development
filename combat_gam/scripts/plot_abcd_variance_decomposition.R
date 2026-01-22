#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
combat_base <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
combat_cog <- if (length(args) >= 3) args[[3]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
combat_pfactor <- if (length(args) >= 4) args[[4]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
out_dir <- if (length(args) >= 5) args[[5]] else "outputs/figures/combat_gam"
combat_cbcl <- if (length(args) >= 6) args[[6]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cbcl.rds"
combat_comp_agecorrected <- if (length(args) >= 7) args[[7]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_comp_agecorrected_baseline.rds"

conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_prefix)) {
  conda_r_lib <- file.path(conda_prefix, "lib", "R", "library")
  if (dir.exists(conda_r_lib)) {
    Sys.setenv(R_LIBS_USER = conda_r_lib, R_LIBS = conda_r_lib)
    .libPaths(c(conda_r_lib, .libPaths()))
  }
}

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(gamm4)
})

base_predictors <- c("age", "sex", "mean_fd", "siteID")
cog_predictors <- c("age", "sex", "mean_fd", "cognition", "siteID")
pfactor_predictors <- c("age", "sex", "mean_fd", "pfactor", "siteID")
cbcl_predictors <- c("age", "sex", "mean_fd", "cbcl_totprob", "siteID")
comp_agecorrected_predictors <- c("age", "sex", "mean_fd", "fluidcomp_agecorrected", "siteID")

maybe_backfill_abcd_from_demopath <- function(dat, needed_cols, demo_path = file.path("demopath", "DemodfScreenFinal.csv")) {
  if (all(needed_cols %in% names(dat))) {
    return(dat)
  }
  if (!"scanID" %in% names(dat)) {
    return(dat)
  }
  if (!file.exists(demo_path)) {
    return(dat)
  }
  demo <- read.csv(demo_path, stringsAsFactors = FALSE)
  join_cols <- intersect(c("scanID", needed_cols), names(demo))
  if (!("scanID" %in% join_cols) || length(join_cols) < 2) {
    return(dat)
  }
  demo <- demo[, join_cols, drop = FALSE]
  demo <- demo[!duplicated(demo$scanID), , drop = FALSE]
  dat %>% left_join(demo, by = "scanID")
}

prepare_raw <- function(path,
                        include_cognition = FALSE,
                        include_pfactor = FALSE,
                        include_cbcl = FALSE,
                        include_comp_agecorrected = FALSE,
                        baseline_only = FALSE) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "eventname")
  if ("subID" %in% names(dat)) {
    needed <- c(needed, "subID")
  }
  backfill_cols <- character()
  if (include_cognition) {
    needed <- c(needed, "nihtbx_fluidcomp_uncorrected")
    backfill_cols <- c(backfill_cols, "nihtbx_fluidcomp_uncorrected")
  }
  if (include_pfactor) {
    needed <- c(needed, "GENERAL")
    backfill_cols <- c(backfill_cols, "GENERAL")
  }
  if (include_cbcl) {
    needed <- c(needed, "cbcl_scr_syn_totprob_r")
    backfill_cols <- c(backfill_cols, "cbcl_scr_syn_totprob_r")
  }
  if (include_comp_agecorrected) {
    needed <- c(needed, "nihtbx_fluidcomp_agecorrected")
    backfill_cols <- c(backfill_cols, "nihtbx_fluidcomp_agecorrected")
  }
  dat <- maybe_backfill_abcd_from_demopath(dat, unique(backfill_cols))
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in raw ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  if (baseline_only) {
    dat <- dat %>% filter(eventname == "baseline_year_1_arm_1")
  }
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(siteID),
      sex = as.factor(sex)
    )
  if (include_cognition) {
    dat$cognition <- dat$nihtbx_fluidcomp_uncorrected
  }
  if (include_pfactor) {
    dat$pfactor <- dat$GENERAL
  }
  if (include_cbcl) {
    dat$cbcl_totprob <- dat$cbcl_scr_syn_totprob_r
  }
  if (include_comp_agecorrected) {
    dat$fluidcomp_agecorrected <- dat$nihtbx_fluidcomp_agecorrected
  }
  list(df = dat, sc_cols = sc_cols)
}

prepare_combat <- function(path,
                           include_cognition = FALSE,
                           include_pfactor = FALSE,
                           include_cbcl = FALSE,
                           include_comp_agecorrected = FALSE) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\..*_h$", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd")
  if ("subID" %in% names(dat)) {
    needed <- c(needed, "subID")
  }
  if (include_cognition) {
    needed <- c(needed, "nihtbx_fluidcomp_uncorrected")
  }
  if (include_pfactor) {
    needed <- c(needed, "GENERAL")
  }
  if (include_cbcl) {
    needed <- c(needed, "cbcl_scr_syn_totprob_r")
  }
  if (include_comp_agecorrected) {
    needed <- c(needed, "nihtbx_fluidcomp_agecorrected")
  }
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in ComBat ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(siteID),
      sex = as.factor(sex)
    )
  if (include_cognition) {
    dat$cognition <- dat$nihtbx_fluidcomp_uncorrected
  }
  if (include_pfactor) {
    dat$pfactor <- dat$GENERAL
  }
  if (include_cbcl) {
    dat$cbcl_totprob <- dat$cbcl_scr_syn_totprob_r
  }
  if (include_comp_agecorrected) {
    dat$fluidcomp_agecorrected <- dat$nihtbx_fluidcomp_agecorrected
  }
  list(df = dat, sc_cols = sc_cols)
}

r2_from_gam <- function(fit) {
  r2 <- tryCatch(summary(fit)$r.sq, error = function(e) NA_real_)
  if (is.null(r2) || is.na(r2) || !is.finite(r2)) {
    return(0)
  }
  as.numeric(r2)
}

build_gam_terms <- function(vars) {
  terms <- character()
  if ("age" %in% vars) {
    terms <- c(terms, "s(age, k=3, bs='tp', fx=TRUE)")
  }
  others <- setdiff(vars, "age")
  if (length(others) > 0) {
    terms <- c(terms, others)
  }
  if (length(terms) == 0) {
    "1"
  } else {
    paste(terms, collapse = " + ")
  }
}

can_use_random_intercept <- function(df, re_var) {
  if (!re_var %in% names(df)) {
    return(FALSE)
  }
  v <- df[[re_var]]
  if (all(is.na(v))) {
    return(FALSE)
  }
  v <- as.factor(v)
  n_obs <- nrow(df)
  n_levels <- nlevels(v)
  if (n_levels < 2) {
    return(FALSE)
  }
  if (n_levels >= n_obs) {
    return(FALSE)
  }
  if (anyDuplicated(v) == 0) {
    return(FALSE)
  }
  TRUE
}

fit_r2_gamm4_abcd <- function(df, vars, re_var = "subID") {
  terms <- build_gam_terms(vars)
  form <- as.formula(paste0("y ~ ", terms))
  if (!can_use_random_intercept(df, re_var)) {
    fit <- mgcv::gam(form, data = df, method = "REML")
    return(r2_from_gam(fit))
  }

  df[[re_var]] <- as.factor(df[[re_var]])
  fit <- tryCatch(
    gamm4::gamm4(form, random = stats::as.formula(paste0("~(1|", re_var, ")")), data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit2 <- mgcv::gam(form, data = df, method = "REML")
    return(r2_from_gam(fit2))
  }

  r2_from_gam(fit$gam)
}

compute_sequential_r2 <- function(y, df, ordered_predictors, re_var = "subID") {
  keep_cols <- unique(c(ordered_predictors, intersect(re_var, names(df))))
  data <- df[, keep_cols, drop = FALSE]
  data$y <- y
  data <- data %>% drop_na()
  if (nrow(data) == 0) {
    return(NULL)
  }

  included <- character(0)
  r2_prev <- fit_r2_gamm4_abcd(data, included, re_var = re_var)
  contrib <- setNames(numeric(length(ordered_predictors)), ordered_predictors)

  for (x in ordered_predictors) {
    included <- c(included, x)
    r2_curr <- fit_r2_gamm4_abcd(data, included, re_var = re_var)
    contrib[[x]] <- r2_curr - r2_prev
    r2_prev <- r2_curr
  }

  list(r2_full = r2_prev, contrib = contrib)
}

compute_variance_decomp <- function(df, sc_cols, label, predictors, strip_suffix = FALSE) {
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "72"))
  if (is.na(cores) || cores < 1) {
    cores <- 40
  }
  cores <- min(cores, 40)
  results <- mclapply(sc_cols, function(col) {
    y <- df[[col]]
    res <- tryCatch(
      compute_sequential_r2(y, df, predictors, re_var = "subID"),
      error = function(e) {
        message(sprintf("[plot_abcd_variance_decomposition] skip %s: %s", col, e$message))
        NULL
      }
    )
    if (is.null(res)) {
      return(NULL)
    }
    edge_base <- if (strip_suffix) sub("_h$", "", col) else col
    data.frame(
      edge = col,
      edge_base = edge_base,
      condition = label,
      total_r2 = res$r2_full,
      t(res$contrib),
      check.names = FALSE
    )
  }, mc.cores = cores)
  rows <- Filter(Negate(is.null), results)
  if (length(rows) == 0) {
    empty <- data.frame(
      edge = character(),
      edge_base = character(),
      condition = character(),
      total_r2 = numeric(),
      check.names = FALSE
    )
    for (p in predictors) {
      empty[[p]] <- numeric()
    }
    return(empty)
  }
  bind_rows(rows)
}

palette <- c(
  siteID = "#F8766D",
  cognition = "#B79F00",
  fluidcomp_agecorrected = "#E64B35",
  mean_fd = "#00BA38",
  sex = "#00B0F6",
  age = "#C77CFF",
  cbcl_totprob = "#4DBBD5",
  pfactor = "#6A3D9A"
)

log_r2_breakdown <- function(df, predictors, label) {
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("[R2] %s: no results to summarize", label))
    return(invisible(NULL))
  }
  needed <- c("condition", "total_r2", predictors)
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    message(sprintf("[R2] %s: missing columns: %s", label, paste(missing, collapse = ", ")))
    return(invisible(NULL))
  }

  summary_df <- df %>%
    group_by(condition) %>%
    summarise(
      n_edges = n(),
      mean_total_r2 = mean(total_r2, na.rm = TRUE),
      median_total_r2 = median(total_r2, na.rm = TRUE),
      across(all_of(predictors), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      across(all_of(predictors), ~ median(.x, na.rm = TRUE), .names = "median_{.col}"),
      .groups = "drop"
    )

  message(sprintf("[R2] %s (mgcv summary r.sq; sequential contributions may be <0)", label))
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, , drop = FALSE]
    cond <- as.character(row$condition[[1]])
    message(sprintf(
      "[R2] %s | %s | n_edges=%d | total_r2 mean=%.4f median=%.4f",
      label, cond, as.integer(row$n_edges[[1]]), row$mean_total_r2[[1]], row$median_total_r2[[1]]
    ))
    for (p in predictors) {
      m <- row[[paste0("mean_", p)]][[1]]
      md <- row[[paste0("median_", p)]][[1]]
      message(sprintf("[R2]   %s: mean=%.4f median=%.4f", p, m, md))
    }
  }
  invisible(summary_df)
}

plot_variant <- function(label,
                         raw_path,
                         combat_path,
                         predictors,
                         include_cognition,
                         include_pfactor,
                         include_cbcl,
                         include_comp_agecorrected,
                         baseline_only,
                         out_prefix) {
  raw_data <- prepare_raw(
    raw_path,
    include_cognition = include_cognition,
    include_pfactor = include_pfactor,
    include_cbcl = include_cbcl,
    include_comp_agecorrected = include_comp_agecorrected,
    baseline_only = baseline_only
  )
  combat_data <- prepare_combat(
    combat_path,
    include_cognition = include_cognition,
    include_pfactor = include_pfactor,
    include_cbcl = include_cbcl,
    include_comp_agecorrected = include_comp_agecorrected
  )

  raw_results <- compute_variance_decomp(raw_data$df, raw_data$sc_cols, "Raw", predictors, strip_suffix = FALSE)
  combat_results <- compute_variance_decomp(combat_data$df, combat_data$sc_cols, "ComBat", predictors, strip_suffix = TRUE)

  if (nrow(raw_results) > 0) {
    order_edges <- raw_results %>%
      arrange(desc(total_r2)) %>%
      pull(edge_base)
  } else if (nrow(combat_results) > 0) {
    order_edges <- combat_results %>%
      arrange(desc(total_r2)) %>%
      pull(edge_base)
  } else {
    message(sprintf("[plot_abcd_variance_decomposition] no valid edges for %s; skip plotting", label))
    return(invisible(NULL))
  }

  combined <- bind_rows(raw_results, combat_results)
  combined$edge_base <- factor(combined$edge_base, levels = order_edges)

  log_r2_breakdown(combined, predictors, label)

  plot_data <- combined %>%
    pivot_longer(cols = all_of(predictors), names_to = "predictor", values_to = "r2")

  plot_data$predictor <- factor(plot_data$predictor, levels = predictors)

  p <- ggplot(plot_data, aes(x = edge_base, y = r2, fill = predictor)) +
    geom_col(width = 0.9, color = "black", linewidth = 0.15) +
    facet_grid(condition ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = palette, breaks = levels(plot_data$predictor)) +
    labs(
      x = "SC edges",
      y = "R square",
      fill = "type",
      title = label
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  ggsave(file.path(out_dir, paste0(out_prefix, ".png")), p, width = 14, height = 8, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")), p, width = 14, height = 8)
}

plot_variant(
  label = "ABCD (age + sex + meanFD)",
  raw_path = raw_rds,
  combat_path = combat_base,
  predictors = base_predictors,
  include_cognition = FALSE,
  include_pfactor = FALSE,
  include_cbcl = FALSE,
  include_comp_agecorrected = FALSE,
  baseline_only = FALSE,
  out_prefix = "abcd_variance_decomp_base"
)

plot_variant(
  label = "ABCD (cognition)",
  raw_path = raw_rds,
  combat_path = combat_cog,
  predictors = cog_predictors,
  include_cognition = TRUE,
  include_pfactor = FALSE,
  include_cbcl = FALSE,
  include_comp_agecorrected = FALSE,
  baseline_only = TRUE,
  out_prefix = "abcd_variance_decomp_cognition_fluid_uncorrected"
)

plot_variant(
  label = "ABCD (p-factor)",
  raw_path = raw_rds,
  combat_path = combat_pfactor,
  predictors = pfactor_predictors,
  include_cognition = FALSE,
  include_pfactor = TRUE,
  include_cbcl = FALSE,
  include_comp_agecorrected = FALSE,
  baseline_only = FALSE,
  out_prefix = "abcd_variance_decomp_pfactor"
)

plot_variant(
  label = "ABCD (CBCL total problems)",
  raw_path = raw_rds,
  combat_path = combat_cbcl,
  predictors = cbcl_predictors,
  include_cognition = FALSE,
  include_pfactor = FALSE,
  include_cbcl = TRUE,
  include_comp_agecorrected = FALSE,
  baseline_only = FALSE,
  out_prefix = "abcd_variance_decomp_cbcl_totprob"
)

plot_variant(
  label = "ABCD (fluid cognition age-corrected; baseline-only)",
  raw_path = raw_rds,
  combat_path = combat_comp_agecorrected,
  predictors = comp_agecorrected_predictors,
  include_cognition = FALSE,
  include_pfactor = FALSE,
  include_cbcl = FALSE,
  include_comp_agecorrected = TRUE,
  baseline_only = TRUE,
  out_prefix = "abcd_variance_decomp_cognition_fluidcomp_agecorrected"
)
