#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
combat_cbcl <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_cbcl/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatCBCLtotalraw.rds"
demopath <- if (length(args) >= 3) args[[3]] else "demopath/DemodfScreenFinal.csv"
out_dir <- if (length(args) >= 4) args[[4]] else "outputs/figures/combat_gam"

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(gamm4)
})

predictors <- c("age", "sex", "mean_fd", "cbcl", "siteID")

prepare_raw <- function(path, demo_path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "eventname")
  if ("subID" %in% names(dat)) {
    needed <- c(needed, "subID")
  }
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in raw ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  if (!"cbcl_scr_syn_totprob_r" %in% names(dat)) {
    demo <- read.csv(demo_path)
    cbcl <- demo[, c("scanID", "cbcl_scr_syn_totprob_r")]
    dat <- dat %>% left_join(cbcl, by = "scanID")
  }
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(siteID),
      sex = as.factor(sex),
      cbcl = cbcl_scr_syn_totprob_r
    )
  list(df = dat, sc_cols = sc_cols)
}

prepare_combat <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\..*_h$", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "cbcl_scr_syn_totprob_r")
  if ("subID" %in% names(dat)) {
    needed <- c(needed, "subID")
  }
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in ComBat ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat %>% drop_na()
  dat <- dat %>%
    mutate(
      siteID = as.factor(siteID),
      sex = as.factor(sex),
      cbcl = cbcl_scr_syn_totprob_r
    )
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

fit_r2_gamm_abcd <- function(df, vars, re_var = "subID") {
  terms <- build_gam_terms(vars)
  form <- as.formula(paste0("y ~ ", terms))
  if (!can_use_random_intercept(df, re_var)) {
    fit <- mgcv::gam(form, data = df, method = "REML")
    return(calc_r2(df$y, fitted(fit)))
  }

  df[[re_var]] <- as.factor(df[[re_var]])
  fit <- tryCatch(
    gamm4::gamm4(form, random = stats::as.formula(paste0("~(1|", re_var, ")")), data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit2 <- mgcv::gam(form, data = df, method = "REML")
    return(calc_r2(df$y, fitted(fit2)))
  }

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
  r2_prev <- fit_r2_gamm_abcd(data, included, re_var = re_var)
  contrib <- setNames(numeric(length(ordered_predictors)), ordered_predictors)
  for (x in ordered_predictors) {
    included <- c(included, x)
    r2_curr <- fit_r2_gamm_abcd(data, included, re_var = re_var)
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
        message(sprintf("[S4th_plot_abcd_variance_decomp_cbcl] skip %s: %s", col, e$message))
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
  mean_fd = "#00BA38",
  sex = "#00B0F6",
  age = "#C77CFF",
  cbcl = "#6A3D9A"
)

raw_data <- prepare_raw(raw_rds, demopath)
combat_data <- prepare_combat(combat_cbcl)

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
  message("[S4th_plot_abcd_variance_decomp_cbcl] no valid edges; skip plotting")
  quit(status = 0)
}

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
    title = "ABCD (CBCL total raw)"
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

summary_row <- plot_data %>%
  group_by(condition, predictor) %>%
  summarise(mean_r2 = mean(r2, na.rm = TRUE), max_r2 = max(r2, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = paste0(predictor, ": mean=", sprintf("%.4f", mean_r2), ", max=", sprintf("%.4f", max_r2))) %>%
  select(condition, label)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

label_df <- summary_row %>%
  group_by(condition) %>%
  summarise(label = paste(label, collapse = "\n"), .groups = "drop") %>%
  mutate(
    edge_base = levels(combined$edge_base)[1],
    r2 = max(plot_data$r2, na.rm = TRUE)
  )

ggsave(file.path(out_dir, "abcd_variance_decomp_cbcl_totalraw.png"), p, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "abcd_variance_decomp_cbcl_totalraw.pdf"), p, width = 14, height = 8)

ymax <- max(plot_data$r2, na.rm = TRUE)
p_fixed <- ggplot(plot_data, aes(x = edge_base, y = r2, fill = predictor)) +
  geom_col(width = 0.9, color = "black", linewidth = 0.15) +
  facet_grid(condition ~ ., scales = "fixed", switch = "y") +
  scale_fill_manual(values = palette, breaks = levels(plot_data$predictor)) +
  scale_y_continuous(limits = c(0, ymax)) +
  geom_text(
    data = label_df,
    aes(x = edge_base, y = r2, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.2
  ) +
  labs(
    x = "SC edges",
    y = "R square",
    fill = "type",
    title = "ABCD (CBCL total raw, fixed y)"
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

ggsave(file.path(out_dir, "abcd_variance_decomp_cbcl_totalraw_fixed.png"), p_fixed, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "abcd_variance_decomp_cbcl_totalraw_fixed.pdf"), p_fixed, width = 14, height = 8)
