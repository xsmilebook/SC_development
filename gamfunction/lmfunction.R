library(dplyr)
set.seed(925)

age_to_years <- function(age_raw) {
  age_num <- as.numeric(age_raw)
  mx <- suppressWarnings(max(age_num, na.rm = TRUE))
  if (is.finite(mx) && mx > 24) {
    return(age_num / 12)
  }
  # Heuristic: some pipelines stored age in "years/12".
  if (is.finite(mx) && mx > 0 && mx <= 2) {
    mx12 <- mx * 12
    mn <- suppressWarnings(min(age_num, na.rm = TRUE))
    mn12 <- mn * 12
    if (is.finite(mx12) && mx12 >= 6 && mx12 <= 30 && is.finite(mn12) && mn12 >= 4) {
      return(age_num * 12)
    }
  }
  age_num
}

get_baseline_age <- function(df, subid_var = "subID", age_var = "age", event_var = "eventname") {
  if (!(subid_var %in% names(df))) stop("Missing column: ", subid_var)
  if (!(age_var %in% names(df))) stop("Missing column: ", age_var)

  subid <- df[[subid_var]]
  age_years <- age_to_years(df[[age_var]])

  if (event_var %in% names(df)) {
    event <- as.character(df[[event_var]])
    is_base <- grepl("base", event, ignore.case = TRUE)
    base_age <- tapply(age_years[is_base], subid[is_base], min, na.rm = TRUE)
  } else {
    base_age <- NULL
  }

  min_age <- tapply(age_years, subid, min, na.rm = TRUE)
  if (is.null(base_age)) return(min_age)

  base_age2 <- min_age
  common_ids <- intersect(names(base_age), names(min_age))
  base_age2[common_ids] <- base_age[common_ids]
  base_age2
}

#### FIT CHANGE-SCORE LM ####
## Compute delta_SC/year within subject and fit:
## change score ~ int_var + covariates
## covariates expected: baseline age, sex, mean_fd (mean across timepoints)
lm.fit.change <- function(region,
                          dataname,
                          int_var,
                          covariates,
                          age_var = "age",
                          subid_var = "subID",
                          event_var = "eventname",
                          stats_only = TRUE) {

  df <- get(dataname)
  needed <- unique(c(subid_var, age_var, int_var, region))
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) stop("Missing required columns in data: ", paste(missing, collapse = ", "))

  # compute time (years from baseline)
  base_age <- get_baseline_age(df, subid_var = subid_var, age_var = age_var, event_var = event_var)
  age_years <- age_to_years(df[[age_var]])
  df$time <- age_years - base_age[match(as.character(df[[subid_var]]), names(base_age))]

  # remove int_var outliers (3 SD)
  int_vec <- df[[int_var]]
  outlier_int <- which(
    int_vec < mean(int_vec, na.rm = TRUE) - 3 * stats::sd(int_vec, na.rm = TRUE) |
      int_vec > mean(int_vec, na.rm = TRUE) + 3 * stats::sd(int_vec, na.rm = TRUE)
  )
  if (length(outlier_int) > 0) df[outlier_int, int_var] <- NA

  # remove y outliers (3 SD)
  y_vec <- df[[region]]
  outlier_y <- which(
    y_vec < mean(y_vec, na.rm = TRUE) - 3 * stats::sd(y_vec, na.rm = TRUE) |
      y_vec > mean(y_vec, na.rm = TRUE) + 3 * stats::sd(y_vec, na.rm = TRUE)
  )
  if (length(outlier_y) > 0) df[outlier_y, region] <- NA

  # keep valid rows
  keep <- !is.na(df[[region]]) & !is.na(df[[int_var]]) & is.finite(df$time)
  df <- df[keep, , drop = FALSE]

  # build subject-level change-score data
  idx_by_sub <- split(seq_len(nrow(df)), as.character(df[[subid_var]]))
  rows <- vector("list", length(idx_by_sub))
  for (k in seq_along(idx_by_sub)) {
    ii <- idx_by_sub[[k]]
    tsub <- df$time[ii]
    if (length(unique(round(tsub, 6))) < 2) next
    i0 <- ii[which.min(tsub)]
    i1 <- ii[which.max(tsub)]
    dt <- df$time[i1] - df$time[i0]
    if (!is.finite(dt) || dt <= 0) next

    out <- data.frame(
      subID = df[[subid_var]][i0],
      change_per_year = (df[[region]][i1] - df[[region]][i0]) / dt,
      age_baseline = age_years[i0],
      sex = df$sex[i0],
      mean_fd_mean = mean(df$mean_fd[ii], na.rm = TRUE),
      int_base = df[[int_var]][i0],
      stringsAsFactors = FALSE
    )
    rows[[k]] <- out
  }
  dat_sub <- dplyr::bind_rows(rows)
  if (nrow(dat_sub) < 10) stop("Too few subjects after change-score construction: ", nrow(dat_sub))

  # outliers on change score and int_var (3 SD)
  cs <- dat_sub$change_per_year
  out_cs <- which(cs < mean(cs, na.rm = TRUE) - 3 * stats::sd(cs, na.rm = TRUE) |
                    cs > mean(cs, na.rm = TRUE) + 3 * stats::sd(cs, na.rm = TRUE))
  if (length(out_cs) > 0) dat_sub <- dat_sub[-out_cs, , drop = FALSE]

  iv <- dat_sub$int_base
  out_iv <- which(iv < mean(iv, na.rm = TRUE) - 3 * stats::sd(iv, na.rm = TRUE) |
                    iv > mean(iv, na.rm = TRUE) + 3 * stats::sd(iv, na.rm = TRUE))
  if (length(out_iv) > 0) dat_sub <- dat_sub[-out_iv, , drop = FALSE]

  dat_sub$sex <- as.factor(dat_sub$sex)

  covariates <- if (is.null(covariates)) "" else as.character(covariates)
  covariates <- trimws(covariates)
  if (!nzchar(covariates)) {
    covariates <- "age_baseline + sex + mean_fd_mean"
  }

  modelformula <- as.formula(sprintf("change_per_year ~ int_base + %s", covariates))
  modelformula.null <- as.formula(sprintf("change_per_year ~ %s", covariates))
  model <- stats::lm(modelformula, data = dat_sub)
  model_null <- stats::lm(modelformula.null, data = dat_sub)
  sm <- summary(model)

  beta_int <- sm$coefficients["int_base", "Estimate"]
  t_int <- sm$coefficients["int_base", "t value"]
  p_int <- sm$coefficients["int_base", "Pr(>|t|)"]

  sse_full <- sum(residuals(model)^2)
  sse_null <- sum(residuals(model_null)^2)
  partialRsq <- (sse_null - sse_full) / sse_null

  stats.results <- data.frame(
    parcel = as.character(region),
    int_var = as.character(int_var),
    n_sub = nrow(dat_sub),
    beta_int = as.numeric(beta_int),
    t_int = as.numeric(t_int),
    p_int = as.numeric(p_int),
    partialRsq = as.numeric(partialRsq),
    stringsAsFactors = FALSE
  )

  if (stats_only) return(stats.results)
  list(stats = stats.results, data = dat_sub)
}

