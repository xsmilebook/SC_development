library(lme4)
library(pbkrtest)
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

fill_baseline_value <- function(df, var, subid_var = "subID", event_var = "eventname") {
  if (!(var %in% names(df))) stop("Missing column: ", var)
  if (!(subid_var %in% names(df))) stop("Missing column: ", subid_var)
  if (!(event_var %in% names(df))) stop("Missing column: ", event_var)

  is_base <- grepl("base", as.character(df[[event_var]]), ignore.case = TRUE)
  base_vals <- tapply(df[[var]][is_base], df[[subid_var]][is_base], function(x) x[which.max(!is.na(x))])
  out <- df[[var]]
  m <- match(as.character(df[[subid_var]]), names(base_vals))
  out[!is.na(m)] <- base_vals[m[!is.na(m)]]
  out
}

make_constant_pred_data <- function(modelobj, time_var, int_var, time_values, int_values) {
  mf <- model.frame(modelobj)
  vars <- setdiff(names(mf), as.character(formula(modelobj)[[2]]))

  pred <- data.frame(init = rep(0, length(time_values) * length(int_values)))
  for (v in vars) {
    if (v == time_var) next
    if (v == int_var) next
    x <- mf[[v]]
    if (is.numeric(x)) {
      pred[[v]] <- rep(stats::median(x, na.rm = TRUE), nrow(pred))
    } else if (is.factor(x) || is.ordered(x)) {
      tab <- table(x)
      lev <- names(tab)[which.max(tab)]
      pred[[v]] <- factor(rep(lev, nrow(pred)), levels = levels(x), ordered = is.ordered(x))
    } else {
      pred[[v]] <- rep(x[[1]], nrow(pred))
    }
  }

  grid <- expand.grid(
    time_value = time_values,
    int_value = int_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred[[time_var]] <- grid$time_value
  pred[[int_var]] <- grid$int_value
  pred$init <- NULL
  pred
}

predict_fixed_only <- function(modelobj, newdata) {
  fixed_formula <- lme4::nobars(formula(modelobj))
  mm <- model.matrix(fixed_formula, newdata)
  beta <- lme4::fixef(modelobj)
  vc <- as.matrix(stats::vcov(modelobj))
  fit <- as.numeric(mm %*% beta)
  se <- as.numeric(sqrt(diag(mm %*% vc %*% t(mm))))
  data.frame(.fitted = fit, .se = se)
}

#### FIT LMM (random intercept + random slope) ####
## Model: SC ~ time + time*int_var + covariates + (1 + time | subID)
## time is years-from-baseline (baseline=0), computed from age-baseline_age within subject.
lmm.time.predict.covariateinteraction <- function(
  region,
  dataname,
  age_var,
  int_var,
  covariates,
  int_var_mode = c("as_is", "baseline_fill"),
  subid_var = "subID",
  event_var = "eventname",
  pb_method = c("KR", "PB"),
  pb_nsim = 1000,
  optimizer = "bobyqa",
  stats_only = TRUE
) {
  int_var_mode <- match.arg(int_var_mode)
  pb_method <- match.arg(pb_method)

  lmm.data <- get(dataname)
  needed <- unique(c(subid_var, age_var, int_var, region))
  missing <- setdiff(needed, names(lmm.data))
  if (length(missing) > 0) stop("Missing required columns in data: ", paste(missing, collapse = ", "))

  if (int_var_mode == "baseline_fill") {
    if (!(event_var %in% names(lmm.data))) stop("baseline_fill requires event_var column: ", event_var)
    int_used <- paste0(int_var, "_base")
    lmm.data[[int_used]] <- fill_baseline_value(lmm.data, int_var, subid_var = subid_var, event_var = event_var)
  } else {
    int_used <- int_var
  }

  # Compute time from baseline age
  base_age <- get_baseline_age(lmm.data, subid_var = subid_var, age_var = age_var, event_var = event_var)
  age_years <- age_to_years(lmm.data[[age_var]])
  m <- match(as.character(lmm.data[[subid_var]]), names(base_age))
  time <- age_years - base_age[m]
  lmm.data$time <- time

  # Clean data: int_var and region outliers, NA filtering
  int_vec <- lmm.data[[int_used]]
  outlier_int <- which(
    int_vec < mean(int_vec, na.rm = TRUE) - 3 * stats::sd(int_vec, na.rm = TRUE) |
      int_vec > mean(int_vec, na.rm = TRUE) + 3 * stats::sd(int_vec, na.rm = TRUE)
  )
  if (length(outlier_int) > 0) lmm.data[outlier_int, int_used] <- NA

  keep <- which(
    !is.na(lmm.data[[int_used]]) &
      !is.na(lmm.data[[region]]) &
      is.finite(lmm.data$time)
  )
  lmm.data <- lmm.data[keep, , drop = FALSE]

  tmp <- lmm.data[[region]]
  outlier_y <- which(
    tmp < mean(tmp, na.rm = TRUE) - 3 * stats::sd(tmp, na.rm = TRUE) |
      tmp > mean(tmp, na.rm = TRUE) + 3 * stats::sd(tmp, na.rm = TRUE)
  )
  if (length(outlier_y) > 0) lmm.data <- lmm.data[-outlier_y, , drop = FALSE]

  # Keep subjects with >=2 distinct time points
  subj <- as.character(lmm.data[[subid_var]])
  tvals <- lmm.data$time
  ok_sub <- tapply(tvals, subj, function(x) length(unique(round(x, 6))) >= 2)
  ok_ids <- names(ok_sub)[ok_sub]
  lmm.data <- lmm.data[subj %in% ok_ids, , drop = FALSE]

  sd_int <- stats::sd(lmm.data[[int_used]], na.rm = TRUE)
  if (is.na(sd_int) || !is.finite(sd_int) || sd_int == 0) {
    stop("int_var has zero/undefined variance after cleaning: ", int_used, " (sd=", sd_int, ")")
  }
  uniq_time <- unique(round(lmm.data$time, 6))
  uniq_time <- uniq_time[is.finite(uniq_time)]
  if (length(uniq_time) < 2) {
    stop("time has <2 unique finite values after cleaning.")
  }

  # Ensure factor covariates remain factors if present
  if ("sex" %in% names(lmm.data)) lmm.data$sex <- as.factor(lmm.data$sex)

  covariates <- if (is.null(covariates)) "" else as.character(covariates)
  covariates <- trimws(covariates)

  fixed_terms <- paste0("time * ", int_used)
  if (nzchar(covariates)) fixed_terms <- paste0(fixed_terms, " + ", covariates)
  modelformula <- stats::as.formula(sprintf("%s ~ %s + (1 + time | %s)", region, fixed_terms, subid_var))
  modelformula.null <- stats::as.formula(sprintf("%s ~ time + %s%s%s + (1 + time | %s)",
    region,
    int_used,
    if (nzchar(covariates)) " + " else "",
    covariates,
    subid_var
  ))

  control <- lme4::lmerControl(optimizer = optimizer)
  m_full <- lme4::lmer(modelformula, data = lmm.data, REML = FALSE, control = control)
  m_null <- lme4::lmer(modelformula.null, data = lmm.data, REML = FALSE, control = control)

  int_pvalue <- if (pb_method == "KR") {
    kr_test <- pbkrtest::KRmodcomp(m_full, m_null)$test
    if (is.data.frame(kr_test)) {
      as.numeric(kr_test[["p.value"]][[1]])
    } else {
      as.numeric(kr_test[1, "p.value"])
    }
  } else {
    refdist <- pbkrtest::PBrefdist(m_full, m_null, nsim = pb_nsim)
    pbkrtest::PBmodcomp(m_full, m_null, ref = refdist)$test["PBtest", "p.value"]
  }

  cf <- lme4::fixef(m_full)
  int_name <- paste0("time:", int_used)
  if (!(int_name %in% names(cf))) {
    # In case model matrix orders interaction term differently (e.g., int_used:time)
    alt <- paste0(int_used, ":time")
    if (!(alt %in% names(cf))) stop("Cannot find interaction term in fixef: ", int_name)
    int_name <- alt
  }
  se_full <- sqrt(diag(as.matrix(stats::vcov(m_full))))
  beta_int <- unname(cf[[int_name]])
  se_int <- unname(se_full[[int_name]])
  t_int <- beta_int / se_int

  # Effect size: predicted difference in change between q90 and q10 of int_var at median follow-up time
  q10 <- as.numeric(stats::quantile(lmm.data[[int_used]], 0.1, na.rm = TRUE))
  q90 <- as.numeric(stats::quantile(lmm.data[[int_used]], 0.9, na.rm = TRUE))
  fu_time <- suppressWarnings(stats::median(tapply(lmm.data$time, lmm.data[[subid_var]], max, na.rm = TRUE), na.rm = TRUE))
  if (is.na(fu_time) || !is.finite(fu_time)) stop("Invalid follow-up time summary (fu_time_median=", fu_time, ").")
  delta_change_q90_vs_q10 <- beta_int * fu_time * (q90 - q10)

  stats.results <- data.frame(
    parcel = as.character(region),
    int_var = as.character(int_used),
    n_obs = nrow(lmm.data),
    n_sub = length(unique(lmm.data[[subid_var]])),
    fu_time_median = as.numeric(fu_time),
    beta_time_int = as.numeric(beta_int),
    se_time_int = as.numeric(se_int),
    t_time_int = as.numeric(t_int),
    p_time_int = as.numeric(int_pvalue),
    int_q10 = q10,
    int_q90 = q90,
    delta_change_q90_vs_q10 = as.numeric(delta_change_q90_vs_q10),
    stringsAsFactors = FALSE
  )

  if (stats_only) return(stats.results)

  pred_data <- make_constant_pred_data(
    modelobj = m_full,
    time_var = "time",
    int_var = int_used,
    time_values = c(0, fu_time),
    int_values = c(q10, q90)
  )
  pred_fit <- predict_fixed_only(m_full, pred_data)
  pred <- cbind(pred_data[, c("time", int_used), drop = FALSE], pred_fit)
  names(pred)[2] <- "int_value"
  pred$int_level <- ifelse(abs(pred$int_value - q10) < 1e-8, "low10", "high90")
  pred$parcel <- as.character(region)

  list(stats = stats.results, predicted = pred)
}
