library(pbkrtest)

pb_lmm_anova <- function(full_model, null_model, nsim = 1000, seed = 925) {
  if (!inherits(full_model, "lmerMod") || !inherits(null_model, "lmerMod")) {
    stop("pb_lmm_anova expects lmerMod objects (full_model, null_model).")
  }
  if (is.na(nsim) || nsim < 1) nsim <- 1000
  if (is.na(seed) || seed < 1) seed <- 925

  set.seed(seed)
  pb <- tryCatch(
    pbkrtest::PBmodcomp(full_model, null_model, nsim = nsim),
    error = function(e) NULL
  )
  if (is.null(pb)) return(NA_real_)

  pval <- tryCatch(
    as.numeric(pb$test["PBtest", "p.value"]),
    error = function(e) NA_real_
  )
  if (!is.finite(pval)) return(NA_real_)
  pval
}
