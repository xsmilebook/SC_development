#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
combat_rds <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_cognition.rds"
output_png <- if (length(args) >= 3) args[[3]] else "outputs/figures/combat_gam/abcd_variance_decomp.png"
output_pdf <- if (length(args) >= 4) args[[4]] else "outputs/figures/combat_gam/abcd_variance_decomp.pdf"

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

predictors <- c("siteID", "age", "sex", "mean_fd", "cognition")

prepare_raw <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "nihtbx_fluidcomp_uncorrected", "eventname")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in raw ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat %>% filter(eventname == "baseline_year_1_arm_1")
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(siteID),
      sex = as.factor(sex),
      cognition = nihtbx_fluidcomp_uncorrected
    )
  list(df = dat, sc_cols = sc_cols)
}

prepare_combat <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\..*_h$", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "nihtbx_fluidcomp_uncorrected")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in ComBat ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(siteID),
      sex = as.factor(sex),
      cognition = nihtbx_fluidcomp_uncorrected
    )
  list(df = dat, sc_cols = sc_cols)
}

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

raw_data <- prepare_raw(raw_rds)
combat_data <- prepare_combat(combat_rds)

raw_results <- compute_variance_decomp(raw_data$df, raw_data$sc_cols, "Raw", strip_suffix = FALSE)
combat_results <- compute_variance_decomp(combat_data$df, combat_data$sc_cols, "ComBat", strip_suffix = TRUE)

order_edges <- raw_results %>%
  arrange(desc(total_r2)) %>%
  pull(edge_base)

combined <- bind_rows(raw_results, combat_results)
combined$edge_base <- factor(combined$edge_base, levels = order_edges)

plot_data <- combined %>%
  pivot_longer(cols = all_of(predictors), names_to = "predictor", values_to = "r2")

palette <- c(
  siteID = "#F8766D",
  cognition = "#B79F00",
  mean_fd = "#00BA38",
  sex = "#00B0F6",
  age = "#C77CFF"
)

plot_data$predictor <- factor(plot_data$predictor, levels = c("siteID", "cognition", "mean_fd", "sex", "age"))

p <- ggplot(plot_data, aes(x = edge_base, y = r2, fill = predictor)) +
  geom_col(width = 0.9, color = "black", linewidth = 0.15) +
  facet_grid(condition ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = palette, breaks = levels(plot_data$predictor)) +
  labs(
    x = "SC edges",
    y = "R square",
    fill = "type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    legend.position = "right"
  )

out_dir <- dirname(output_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

ggsave(output_png, p, width = 14, height = 8, dpi = 300)

ggsave(output_pdf, p, width = 14, height = 8)
