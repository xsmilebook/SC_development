#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_rds <- if (length(args) >= 1) args[[1]] else "/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
combat_cbcl <- if (length(args) >= 2) args[[2]] else "outputs/results/combat_cbcl/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatCBCLtotalraw.rds"
demopath <- if (length(args) >= 3) args[[3]] else "demopath/DemodfScreenFinal.csv"
out_dir <- if (length(args) >= 4) args[[4]] else "outputs/results/combat_cbcl"

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
})

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
  if ("subID" %in% names(dat)) {
    needed <- c(needed, "subID")
  }
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

fit_r2_gamm_abcd <- function(df, vars, re_var = "subID") {
  terms <- build_gam_terms(vars)
  if (!re_var %in% names(df)) {
    fit <- mgcv::gam(as.formula(paste0("y ~ ", terms)), data = df, method = "REML")
    return(calc_r2(df$y, fitted(fit)))
  }
  df[[re_var]] <- as.factor(df[[re_var]])
  formula <- as.formula(paste0("y ~ ", terms, " + s(", re_var, ", bs='re')"))
  fit <- mgcv::gam(formula, data = df, method = "REML")
  calc_r2(df$y, fitted(fit))
}

compute_delta_r2 <- function(y, df, predictors, re_var = "subID") {
  keep_cols <- unique(c(predictors, intersect(re_var, names(df))))
  data <- df[, keep_cols, drop = FALSE]
  data$y <- y
  data <- data[complete.cases(data), , drop = FALSE]
  if (nrow(data) == 0) {
    delta <- setNames(rep(NA_real_, length(predictors)), predictors)
    return(list(r2_full = NA_real_, delta = delta))
  }
  r2_full <- fit_r2_gamm_abcd(data, predictors, re_var = re_var)
  delta <- setNames(numeric(length(predictors)), predictors)
  for (x in predictors) {
    reduced <- setdiff(predictors, x)
    r2_reduced <- fit_r2_gamm_abcd(data, reduced, re_var = re_var)
    delta[[x]] <- r2_full - r2_reduced
  }
  list(r2_full = r2_full, delta = delta)
}

compute_site_r2 <- function(df, sc_cols) {
  predictors <- c("siteID", "age", "sex", "mean_fd", "cbcl")
  keep <- predictors
  if ("subID" %in% names(df)) {
    keep <- c(keep, "subID")
  }
  sapply(sc_cols, function(col) {
    y <- df[[col]]
    d <- df[, keep]
    d$y <- y
    d <- d[complete.cases(d), , drop = FALSE]
    compute_delta_r2(d$y, d, predictors, re_var = "subID")$delta[["siteID"]]
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
