library(lme4)

lmm.age.random.slope <- function(region,
                                 dataname,
                                 age_var = "age",
                                 sex_var = "sex",
                                 fd_var = "mean_fd",
                                 subid_var = "subID",
                                 min_obs = 10,
                                 stats_only = TRUE,
                                 return_slopes = FALSE,
                                 return_model = FALSE) {
  df <- get(dataname)
  needed <- unique(c(region, age_var, sex_var, fd_var, subid_var))
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

  keep <- !is.na(df[[region]]) & !is.na(df[[age_var]]) & !is.na(df[[sex_var]]) & !is.na(df[[fd_var]])
  df <- df[keep, , drop = FALSE]

  if (nrow(df) < min_obs) {
    out <- data.frame(
      parcel = as.character(region),
      ok = FALSE,
      beta_age = NA_real_,
      t_age = NA_real_,
      rand_age_mean = NA_real_,
      n_obs = nrow(df),
      n_sub = length(unique(df[[subid_var]])),
      stringsAsFactors = FALSE
    )
    if (return_slopes || return_model || !stats_only) {
      return(list(stats = out, rand_slopes = data.frame(), model = NULL))
    }
    return(out)
  }

  df[[sex_var]] <- as.factor(df[[sex_var]])
  fml <- as.formula(
    sprintf("%s ~ %s + %s + %s + (1 + %s || %s)",
            region, age_var, sex_var, fd_var, age_var, subid_var)
  )

  ctrl <- lmerControl(checkConv = list(grad = lme4::.makeCC(action = "warning", tol = 0.02)))
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

  rand_slopes <- data.frame(
    subID = rownames(re),
    fixed_slope = as.numeric(beta_age),
    random_slope = as.numeric(re[, age_var]),
    stringsAsFactors = FALSE
  )

  if (return_slopes || return_model || !stats_only) {
    return(list(stats = out, rand_slopes = rand_slopes, model = mod))
  }
  out
}
