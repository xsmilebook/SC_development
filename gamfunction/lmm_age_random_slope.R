library(lme4)

lmm.age.random.slope <- function(region,
                                 dataname,
                                 age_var = "age",
                                 sex_var = "sex",
                                 fd_var = "mean_fd",
                                 subid_var = "subID",
                                 min_obs = 10,
                                 stats_only = TRUE) {
  df <- get(dataname)
  needed <- unique(c(region, age_var, sex_var, fd_var, subid_var))
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

  keep <- !is.na(df[[region]]) & !is.na(df[[age_var]]) & !is.na(df[[sex_var]]) & !is.na(df[[fd_var]])
  df <- df[keep, , drop = FALSE]

  ok <- tapply(df[[age_var]], df[[subid_var]], function(x) length(unique(round(x, 6))) >= 2)
  df <- df[df[[subid_var]] %in% names(ok)[ok], , drop = FALSE]
  if (nrow(df) < min_obs) {
    return(data.frame(
      parcel = as.character(region),
      ok = FALSE,
      beta_age = NA_real_,
      t_age = NA_real_,
      rand_age_mean = NA_real_,
      n_obs = nrow(df),
      n_sub = length(unique(df[[subid_var]])),
      stringsAsFactors = FALSE
    ))
  }

  df[[sex_var]] <- as.factor(df[[sex_var]])
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    check.nobs.vs.nRE = "ignore",
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev = "ignore"
  )
  fml <- as.formula(
    sprintf("%s ~ %s + %s + %s + (1 + %s || %s)",
            region, age_var, sex_var, fd_var, age_var, subid_var)
  )

  mod <- suppressWarnings(lmer(fml, data = df, REML = FALSE, control = ctrl))
  sm <- summary(mod)
  beta_age <- sm$coefficients[age_var, "Estimate"]
  t_age <- sm$coefficients[age_var, "t value"]
  re <- ranef(mod)[[subid_var]]
  rand_age_mean <- mean(abs(re[, age_var]), na.rm = TRUE)
  out <- data.frame(
    parcel = as.character(region),
    ok = TRUE,
    beta_age = as.numeric(beta_age),
    t_age = as.numeric(t_age),
    rand_age_mean = as.numeric(rand_age_mean),
    n_obs = nrow(df),
    n_sub = length(unique(df[[subid_var]])),
    stringsAsFactors = FALSE
  )

  if (stats_only) return(out)
  list(stats = out, data = df)
}
