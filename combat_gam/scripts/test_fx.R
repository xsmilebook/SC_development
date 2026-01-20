#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
})

# -----------------------------
# Config
# -----------------------------
rds_path <- "D:/code/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge_new.rds"
age_col  <- "Age"  # <- 如果你的年龄列名不是 Age，请改这里
sc_cols  <- paste0("SC.", 1:78)

# 数值比较容差（根据数据尺度可调）
tol_fit_rmse <- 1e-10
tol_pred_rmse <- 1e-10

# 输出文件
out_csv <- "D:/code/compare_fxTRUE_vs_penalized_k3_tp_SC1toSC78.csv"

# -----------------------------
# Load data
# -----------------------------
df <- readRDS(rds_path)

if (!age_col %in% names(df)) {
  stop(sprintf("Age column '%s' not found in data. Available columns include: %s",
               age_col, paste(head(names(df), 30), collapse = ", ")))
}

missing_sc <- setdiff(sc_cols, names(df))
if (length(missing_sc) > 0) {
  stop(sprintf("Missing SC columns: %s", paste(missing_sc, collapse = ", ")))
}

# Age grid for prediction comparison
age_vec <- df[[age_col]]
age_vec <- age_vec[is.finite(age_vec)]
if (length(age_vec) < 5) stop("Not enough finite Age values to build prediction grid.")

age_grid <- seq(min(age_vec), max(age_vec), length.out = 200)
newdata_grid <- data.frame(tmpAge = age_grid)
names(newdata_grid)[1] <- age_col  # rename to actual age_col

# -----------------------------
# Helper: fit & compare one column
# -----------------------------
fit_compare_one <- function(y_col, df, age_col, newdata_grid,
                            tol_fit_rmse, tol_pred_rmse) {
  
  # Subset to complete cases for this y and age
  keep <- is.finite(df[[age_col]]) & is.finite(df[[y_col]])
  d <- df[keep, c(age_col, y_col)]
  names(d) <- c("Age__", "Y__")  # standardize names inside function
  
  n <- nrow(d)
  if (n < 10) {
    return(data.frame(
      y_col = y_col, n = n,
      edf_fx = NA_real_, edf_pen = NA_real_,
      sp_pen = NA_real_,
      rss_fx = NA_real_, rss_pen = NA_real_,
      aic_fx = NA_real_, aic_pen = NA_real_,
      fit_rmse = NA_real_, fit_maxabs = NA_real_,
      pred_rmse = NA_real_, pred_maxabs = NA_real_,
      identical_fit = NA, identical_pred = NA,
      note = "Too few valid rows"
    ))
  }
  
  # Build formulas (using standardized names)
  f_fx  <- as.formula("Y__ ~ s(Age__, k=3, bs='tp', fx=TRUE)")
  f_pen <- as.formula("Y__ ~ s(Age__, k=3, bs='tp')")
  
  # Fit models
  m_fx  <- gam(f_fx,  data = d, method = "REML")
  m_pen <- gam(f_pen, data = d, method = "REML")
  
  # Extract edf / sp
  # For a single smooth term, s.table[1, "edf"] is fine
  st_fx  <- summary(m_fx)$s.table
  st_pen <- summary(m_pen)$s.table
  
  edf_fx  <- as.numeric(st_fx[1, "edf"])
  edf_pen <- as.numeric(st_pen[1, "edf"])
  
  # Penalized smooth has estimated smoothing parameter (sp); fx=TRUE typically has NA or fixed
  sp_pen <- tryCatch({
    sp <- m_pen$sp
    if (length(sp) >= 1) as.numeric(sp[1]) else NA_real_
  }, error = function(e) NA_real_)
  
  # Compare fitted values on training data
  fit_fx  <- fitted(m_fx)
  fit_pen <- fitted(m_pen)
  
  fit_diff <- fit_fx - fit_pen
  fit_rmse <- sqrt(mean(fit_diff^2))
  fit_maxabs <- max(abs(fit_diff))
  
  # Compare predictions on Age grid
  nd <- newdata_grid
  names(nd) <- "Age__"  # align to standardized name
  
  pred_fx  <- as.numeric(predict(m_fx,  newdata = nd, type = "response"))
  pred_pen <- as.numeric(predict(m_pen, newdata = nd, type = "response"))
  
  pred_diff <- pred_fx - pred_pen
  pred_rmse <- sqrt(mean(pred_diff^2))
  pred_maxabs <- max(abs(pred_diff))
  
  # Basic fit stats
  rss_fx  <- sum(residuals(m_fx)^2)
  rss_pen <- sum(residuals(m_pen)^2)
  
  aic_fx  <- AIC(m_fx)
  aic_pen <- AIC(m_pen)
  
  identical_fit  <- (fit_rmse <= tol_fit_rmse)
  identical_pred <- (pred_rmse <= tol_pred_rmse)
  
  note <- ""
  # Common situation: with k=3, penalized smooth may shrink to ~linear or near-fixed,
  # producing near-identical fits.
  # We'll just annotate when extremely close.
  if (identical_fit && identical_pred) note <- "Practically identical under tolerance"
  
  data.frame(
    y_col = y_col, n = n,
    edf_fx = edf_fx, edf_pen = edf_pen,
    sp_pen = sp_pen,
    rss_fx = rss_fx, rss_pen = rss_pen,
    aic_fx = aic_fx, aic_pen = aic_pen,
    fit_rmse = fit_rmse, fit_maxabs = fit_maxabs,
    pred_rmse = pred_rmse, pred_maxabs = pred_maxabs,
    identical_fit = identical_fit, identical_pred = identical_pred,
    note = note,
    row.names = NULL
  )
}

# -----------------------------
# Run over SC.1 ~ SC.78
# -----------------------------
results_list <- lapply(sc_cols, fit_compare_one,
                       df = df, age_col = age_col, newdata_grid = newdata_grid,
                       tol_fit_rmse = tol_fit_rmse, tol_pred_rmse = tol_pred_rmse)

res <- do.call(rbind, results_list)

# Add a convenience flag: "any difference"
res$any_difference <- !(res$identical_fit & res$identical_pred)

# Sort by biggest prediction discrepancy
res <- res[order(res$pred_rmse, decreasing = TRUE), ]

# Save
write.csv(res, out_csv, row.names = FALSE)

# Print a quick summary
cat("Done.\n")
cat("Output:", out_csv, "\n\n")

cat("How many columns show any difference (under tolerance)?\n")
cat(sum(res$any_difference, na.rm = TRUE), "out of", nrow(res), "\n\n")

cat("Top 10 by pred_rmse:\n")
print(head(res[, c("y_col","n","edf_fx","edf_pen","sp_pen","pred_rmse","pred_maxabs","fit_rmse","fit_maxabs","aic_fx","aic_pen","any_difference")], 10))
