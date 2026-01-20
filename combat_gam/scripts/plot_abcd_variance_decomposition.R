#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
combat_base <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
combat_cog <- if (length(args) >= 3) args[[3]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
combat_pfactor <- if (length(args) >= 4) args[[4]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_pfactor.rds"
out_dir <- if (length(args) >= 5) args[[5]] else "outputs/figures/combat_gam"

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

prepare_raw <- function(path, include_cognition = FALSE, include_pfactor = FALSE, baseline_only = FALSE) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "eventname")
  if ("subID" %in% names(dat)) {
    needed <- c(needed, "subID")
  }
  if (include_cognition) {
    needed <- c(needed, "nihtbx_fluidcomp_uncorrected")
  }
  if (include_pfactor) {
    needed <- c(needed, "GENERAL")
  }
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
  list(df = dat, sc_cols = sc_cols)
}

prepare_combat <- function(path, include_cognition = FALSE, include_pfactor = FALSE) {
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
  list(df = dat, sc_cols = sc_cols)
}

calc_r2 <- function(y, fitted) {
  y <- as.numeric(y)
  fitted <- as.numeric(fitted)
  tss <- sum((y - mean(y))^2)
  rss <- sum((y - fitted)^2)
  if (tss == 0) {
    return(0)
  }
  1 - rss / tss
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

fit_r2_gamm4_abcd <- function(df, vars, re_var = "subID") {
  terms <- build_gam_terms(vars)
  if (!re_var %in% names(df)) {
    fit <- mgcv::gam(as.formula(paste0("y ~ ", terms)), data = df, method = "REML")
    return(calc_r2(df$y, fitted(fit)))
  }
  df[[re_var]] <- as.factor(df[[re_var]])
  form <- as.formula(paste0("y ~ ", terms))
  fit <- gamm4::gamm4(form, random = stats::as.formula(paste0("~(1|", re_var, ")")), data = df)
  calc_r2(df$y, fitted(fit$mer))
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
    cores <- 72
  }
  cores <- min(cores, 72)
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
  bind_rows(Filter(Negate(is.null), results))
}

palette <- c(
  siteID = "#F8766D",
  cognition = "#B79F00",
  mean_fd = "#00BA38",
  sex = "#00B0F6",
  age = "#C77CFF",
  pfactor = "#6A3D9A"
)

plot_variant <- function(label, raw_path, combat_path, predictors, include_cognition, include_pfactor, baseline_only, out_prefix) {
  raw_data <- prepare_raw(raw_path, include_cognition = include_cognition, include_pfactor = include_pfactor, baseline_only = baseline_only)
  combat_data <- prepare_combat(combat_path, include_cognition = include_cognition, include_pfactor = include_pfactor)

  raw_results <- compute_variance_decomp(raw_data$df, raw_data$sc_cols, "Raw", predictors, strip_suffix = FALSE)
  combat_results <- compute_variance_decomp(combat_data$df, combat_data$sc_cols, "ComBat", predictors, strip_suffix = TRUE)

  order_edges <- raw_results %>%
    arrange(desc(total_r2)) %>%
    pull(edge_base)

  combined <- bind_rows(raw_results, combat_results)
  combined$edge_base <- factor(combined$edge_base, levels = order_edges)

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
  baseline_only = TRUE,
  out_prefix = "abcd_variance_decomp_cognition"
)

plot_variant(
  label = "ABCD (p-factor)",
  raw_path = raw_rds,
  combat_path = combat_pfactor,
  predictors = pfactor_predictors,
  include_cognition = FALSE,
  include_pfactor = TRUE,
  baseline_only = FALSE,
  out_prefix = "abcd_variance_decomp_pfactor"
)
