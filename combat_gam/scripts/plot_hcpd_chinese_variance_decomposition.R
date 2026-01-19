#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
hcpd_raw <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
hcpd_combat <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
ch_raw <- if (length(args) >= 3) args[[3]] else "/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
ch_combat <- if (length(args) >= 4) args[[4]] else "outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

predictors <- c("siteID", "age", "sex", "mean_fd")

subset_bits <- function(mask, p) {
  bitwAnd(mask, bitwShiftL(1L, 0:(p - 1))) > 0
}

r2_for_subset <- function(y, df, vars) {
  if (length(vars) == 0) {
    fit <- lm(y ~ 1, data = df)
  } else {
    formula <- as.formula(paste("y ~", paste(vars, collapse = " + ")))
    fit <- lm(formula, data = df)
  }
  tss <- sum((y - mean(y))^2)
  rss <- sum(residuals(fit)^2)
  if (tss == 0) {
    return(0)
  }
  1 - rss / tss
}

shapley_r2 <- function(y, df, vars) {
  p <- length(vars)
  masks <- 0:(2^p - 1)
  r2_vals <- numeric(length(masks))
  names(r2_vals) <- masks
  for (i in seq_along(masks)) {
    mask <- masks[[i]]
    included <- vars[subset_bits(mask, p)]
    r2_vals[[i]] <- r2_for_subset(y, df, included)
  }
  total_r2 <- r2_vals[[as.character(2^p - 1)]]
  contrib <- numeric(p)
  for (j in seq_len(p)) {
    weights <- numeric(length(masks))
    for (i in seq_along(masks)) {
      mask <- masks[[i]]
      if (bitwAnd(mask, bitwShiftL(1L, j - 1)) == 0) {
        s <- sum(subset_bits(mask, p))
        weight <- factorial(s) * factorial(p - s - 1) / factorial(p)
        with_mask <- mask + bitwShiftL(1L, j - 1)
        weights[[i]] <- weight * (r2_vals[[as.character(with_mask)]] - r2_vals[[as.character(mask)]])
      }
    }
    contrib[[j]] <- sum(weights)
  }
  names(contrib) <- vars
  list(total_r2 = total_r2, contrib = contrib)
}

compute_variance_decomp <- function(df, sc_cols, label, strip_suffix = FALSE) {
  results <- lapply(sc_cols, function(col) {
    y <- df[[col]]
    data <- df[, predictors]
    data$y <- y
    data <- data %>% drop_na()
    res <- shapley_r2(y, data, predictors)
    edge_base <- if (strip_suffix) sub("_h$", "", col) else col
    data.frame(
      edge = col,
      edge_base = edge_base,
      condition = label,
      total_r2 = res$total_r2,
      t(res$contrib),
      check.names = FALSE
    )
  })
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

build_annotation <- function(plot_data, predictors) {
  summary <- plot_data %>%
    group_by(condition, predictor) %>%
    summarise(mean = mean(r2, na.rm = TRUE), max = max(r2, na.rm = TRUE), .groups = "drop")

  summary$line <- sprintf("%s: mean=%.4f, max=%.4f", summary$predictor, summary$mean, summary$max)
  ann <- summary %>%
    arrange(match(predictor, predictors)) %>%
    group_by(condition) %>%
    summarise(label = paste(line, collapse = "\n"), .groups = "drop")

  ann
}

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

  ann <- build_annotation(plot_data, predictors)
  ann$x <- levels(plot_data$edge_base)[[1]]
  ann$y <- Inf

  p <- ggplot(plot_data, aes(x = edge_base, y = r2, fill = predictor)) +
    geom_col(width = 0.9, color = "black", linewidth = 0.15) +
    facet_grid(condition ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = palette, breaks = levels(plot_data$predictor)) +
    geom_text(
      data = ann,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.1,
      size = 3
    ) +
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
plot_variant("HCPD", hcpd_raw_data, hcpd_combat_data, "hcpd_variance_decomp")

ch_raw_data <- prepare_chinese_raw(ch_raw)
ch_combat_data <- prepare_chinese_combat(ch_combat)
plot_variant("Chinese Cohort", ch_raw_data, ch_combat_data, "chinese_variance_decomp")
