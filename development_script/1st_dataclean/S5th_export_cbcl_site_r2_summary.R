#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
combat_cbcl <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_cbcl/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatCBCLtotalraw.rds"
demopath <- if (length(args) >= 3) args[[3]] else "demopath/DemodfScreenFinal.csv"
out_dir <- if (length(args) >= 4) args[[4]] else "outputs/results/combat_cbcl"

suppressPackageStartupMessages({
  library(dplyr)
})

prepare_raw <- function(path, demo_path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\.", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "eventname")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in raw ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  demo <- read.csv(demo_path)
  cbcl <- demo[, c("scanID", "cbcl_scr_syn_totprob_r")]
  dat <- merge(dat, cbcl, by = "scanID", all.x = TRUE)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat$siteID <- as.factor(dat$siteID)
  dat$sex <- as.factor(dat$sex)
  dat$cbcl <- dat$cbcl_scr_syn_totprob_r
  list(df = dat, sc_cols = sc_cols)
}

prepare_combat <- function(path) {
  dat <- readRDS(path)
  sc_cols <- grep("^SC\\..*_h$", names(dat), value = TRUE)
  needed <- c(sc_cols, "scanID", "siteID", "age", "sex", "mean_fd", "cbcl_scr_syn_totprob_r")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop(paste("Missing columns in ComBat ABCD:", paste(missing, collapse = ", ")))
  }
  dat <- dat[, needed]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat$siteID <- as.factor(dat$siteID)
  dat$sex <- as.factor(dat$sex)
  dat$cbcl <- dat$cbcl_scr_syn_totprob_r
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
  if (tss == 0) return(0)
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
  contrib
}

compute_site_r2 <- function(df, sc_cols) {
  predictors <- c("siteID", "age", "sex", "mean_fd", "cbcl")
  sapply(sc_cols, function(col) {
    y <- df[[col]]
    d <- df[, predictors]
    d$y <- y
    d <- d[complete.cases(d), , drop = FALSE]
    shapley_r2(d$y, d, predictors)[["siteID"]]
  })
}

raw <- prepare_raw(raw_rds, demopath)
combat <- prepare_combat(combat_cbcl)

raw_site <- compute_site_r2(raw$df, raw$sc_cols)
combat_site <- compute_site_r2(combat$df, combat$sc_cols)

raw_df <- data.frame(edge = raw$sc_cols, site_r2_raw = raw_site, stringsAsFactors = FALSE)
combat_df <- data.frame(edge = sub("_h$", "", combat$sc_cols), site_r2_combat = combat_site, stringsAsFactors = FALSE)

merged <- raw_df %>%
  inner_join(combat_df, by = "edge") %>%
  arrange(desc(site_r2_combat))

summary_row <- data.frame(
  edge = "SUMMARY",
  site_r2_raw = mean(raw_site),
  site_r2_combat = mean(combat_site),
  stringsAsFactors = FALSE
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write.csv(merged, file.path(out_dir, "cbcl_site_r2_raw_vs_combat.csv"), row.names = FALSE)
write.csv(summary_row, file.path(out_dir, "cbcl_site_r2_summary.csv"), row.names = FALSE)
