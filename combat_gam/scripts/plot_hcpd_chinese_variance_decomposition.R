#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
hcpd_raw <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
hcpd_combat <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
ch_raw <- if (length(args) >= 3) args[[3]] else "/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge_new.rds"
ch_combat <- if (length(args) >= 4) args[[4]] else "outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
out_dir <- if (length(args) >= 5) args[[5]] else "outputs/figures/combat_gam"
hcpd_prefix <- if (length(args) >= 6) args[[6]] else "hcpd_variance_decomp"
ch_prefix <- if (length(args) >= 7) args[[7]] else "chinese_variance_decomp"

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
})

predictors <- c("siteID", "age", "sex", "mean_fd")

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

fit_r2_gam <- function(df, vars) {
  formula <- as.formula(paste0("y ~ ", build_gam_terms(vars)))
  fit <- mgcv::gam(formula, data = df, method = "REML")
  calc_r2(df$y, fitted(fit))
}

compute_delta_r2 <- function(y, df, predictors) {
  data <- df[, predictors, drop = FALSE]
  data$y <- y
  data <- data %>% drop_na()
  if (nrow(data) == 0) {
    delta <- setNames(rep(NA_real_, length(predictors)), predictors)
    return(list(r2_full = NA_real_, delta = delta))
  }

  r2_full <- fit_r2_gam(data, predictors)
  delta <- setNames(numeric(length(predictors)), predictors)
  for (x in predictors) {
    reduced <- setdiff(predictors, x)
    r2_reduced <- fit_r2_gam(data, reduced)
    delta[[x]] <- r2_full - r2_reduced
  }
  list(r2_full = r2_full, delta = delta)
}

compute_variance_decomp <- function(df, sc_cols, label, strip_suffix = FALSE) {
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (is.na(cores) || cores < 1) {
    cores <- 1
  }
  results <- mclapply(sc_cols, function(col) {
    y <- df[[col]]
    res <- compute_delta_r2(y, df, predictors)
    edge_base <- if (strip_suffix) sub("_h$", "", col) else col
    data.frame(
      edge = col,
      edge_base = edge_base,
      condition = label,
      total_r2 = res$r2_full,
      t(res$delta),
      check.names = FALSE
    )
  }, mc.cores = cores)
  bind_rows(results)
}

prepare_hcpd_raw <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "subID", "site", "age", "sex", "mean_fd")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in HCPD raw:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(site),
      sex = as.factor(sex)
    )
  list(df = dat, sc_cols = sc_cols)
}

prepare_hcpd_combat <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\..*_h$", names(dat), value = TRUE)
  needed <- c(sc_cols, "subID", "site", "age", "sex", "mean_fd")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in HCPD ComBat:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(site),
      sex = as.factor(sex)
    )
  list(df = dat, sc_cols = sc_cols)
}

prepare_chinese_raw <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "subID", "study", "Age", "Sex", "mean_fd")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in Chinese raw:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(study),
      age = Age,
      sex = as.factor(Sex)
    )
  list(df = dat, sc_cols = sc_cols)
}

prepare_chinese_combat <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\..*_h$", names(dat), value = TRUE)
  needed <- c(sc_cols, "subID", "study", "Age", "Sex", "mean_fd")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in Chinese ComBat:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(study),
      age = Age,
      sex = as.factor(Sex)
    )
  list(df = dat, sc_cols = sc_cols)
}

palette <- c(
  siteID = "#F8766D",
  mean_fd = "#00BA38",
  sex = "#00B0F6",
  age = "#C77CFF"
)

plot_variant <- function(label, raw_data, combat_data, out_prefix) {
  raw_results <- compute_variance_decomp(raw_data$df, raw_data$sc_cols, "Raw", strip_suffix = FALSE)
  combat_results <- compute_variance_decomp(combat_data$df, combat_data$sc_cols, "ComBat", strip_suffix = TRUE)

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

hcpd_raw_data <- prepare_hcpd_raw(hcpd_raw)
hcpd_combat_data <- prepare_hcpd_combat(hcpd_combat)
plot_variant("HCPD", hcpd_raw_data, hcpd_combat_data, hcpd_prefix)

ch_raw_data <- prepare_chinese_raw(ch_raw)
ch_combat_data <- prepare_chinese_combat(ch_combat)
plot_variant("Chinese Cohort", ch_raw_data, ch_combat_data, ch_prefix)
