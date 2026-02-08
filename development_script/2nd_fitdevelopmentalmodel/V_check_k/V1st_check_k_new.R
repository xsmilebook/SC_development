## Check GAM k values using ComBat-GAM data.
## ABCD uses GAMM (gamm4) for AIC comparison.
##
## Usage (examples):
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/V_check_k/V1st_check_k_new.R \
##     --dataset=chinese
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/V_check_k/V1st_check_k_new.R \
##     --dataset=abcd --k_values=3,4,5,6 --n_edges=10
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/V_check_k/V1st_check_k_new.R \
##     --input_rds=/path/to/SCdata.rds --output_dir=/path/to/output --figure_dir=/path/to/figs
##   Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/V_check_k/V1st_check_k_new.R \
##     --bootstrap_n=1000

rm(list = ls())

library(mgcv)
library(parallel)

parse_args <- function(args) {
  res <- list()
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    res[[kv[[1]]]] <- kv[[2]]
  }
  res
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(if (!is.null(args$project_root)) args$project_root else getwd(), mustWork = FALSE)
if (!file.exists(file.path(project_root, "ARCHITECTURE.md"))) {
  stop("project_root does not look like SCDevelopment (missing ARCHITECTURE.md): ", project_root)
}

dataset <- tolower(if (!is.null(args$dataset)) args$dataset else "chinese")

resolve_path <- function(path) {
  if (is.null(path)) return(NULL)
  if (!grepl("^/", path)) {
    path <- file.path(project_root, path)
  }
  normalizePath(path, mustWork = FALSE)
}

input_rds <- if (!is.null(args$input_rds)) {
  resolve_path(args$input_rds)
} else {
  switch(
    dataset,
    hcpd = file.path(project_root, "outputs", "results", "combat_gam", "hcpd",
                     "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"),
    abcd = file.path(project_root, "outputs", "results", "combat_gam", "abcd",
                     "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"),
    chinese = file.path(project_root, "outputs", "results", "combat_gam", "chinese",
                        "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"),
    stop("Unsupported dataset: ", dataset, " (expected hcpd|abcd|chinese)")
  )
}
if (!file.exists(input_rds)) stop("Missing input_rds: ", input_rds)

k_values_raw <- if (!is.null(args$k_values)) args$k_values else "3,4,5,6"
k_values <- as.integer(strsplit(k_values_raw, ",", fixed = TRUE)[[1]])
k_values <- k_values[!is.na(k_values)]
if (length(k_values) == 0) stop("No valid k_values: ", k_values_raw)

bootstrap_n <- as.integer(if (!is.null(args$bootstrap_n)) args$bootstrap_n else 1000)
if (is.na(bootstrap_n) || bootstrap_n < 0) stop("Invalid bootstrap_n: ", args$bootstrap_n)

output_dir <- resolve_path(if (!is.null(args$output_dir)) args$output_dir else {
  file.path(project_root, "outputs", "results", "2nd_fitdevelopmentalmodel", dataset, "check_k")
})
figure_dir <- resolve_path(if (!is.null(args$figure_dir)) args$figure_dir else {
  file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", dataset, "check_k")
})

if (!startsWith(output_dir, project_root)) stop("output_dir must be under project_root: ", output_dir)
if (!startsWith(figure_dir, project_root)) stop("figure_dir must be under project_root: ", figure_dir)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

scdata <- readRDS(input_rds)

age_col <- if ("age" %in% names(scdata)) "age" else if ("Age" %in% names(scdata)) "Age" else NA_character_
sex_col <- if ("sex" %in% names(scdata)) "sex" else if ("Sex" %in% names(scdata)) "Sex" else NA_character_
mean_fd_col <- if ("mean_fd" %in% names(scdata)) "mean_fd" else if ("meanFD" %in% names(scdata)) "meanFD" else NA_character_

if (is.na(age_col)) stop("Missing age/Age column in input_rds.")
if (is.na(sex_col)) stop("Missing sex/Sex column in input_rds.")
if (is.na(mean_fd_col)) stop("Missing mean_fd column in input_rds.")

if (age_col != "age") scdata$age <- scdata[[age_col]]
if (sex_col != "sex") scdata$sex <- scdata[[sex_col]]
if (mean_fd_col != "mean_fd") scdata$mean_fd <- scdata[[mean_fd_col]]
scdata$sex <- as.factor(scdata$sex)

if (dataset == "abcd") {
  if (!requireNamespace("gamm4", quietly = TRUE)) {
    stop("Missing package 'gamm4' required for ABCD GAMM fitting.")
  }
  if (!("subID" %in% names(scdata))) {
    if ("subjectkey" %in% names(scdata)) {
      scdata$subID <- scdata$subjectkey
    } else {
      stop("ABCD GAMM requires subID (or subjectkey) column for random effects.")
    }
  }
}

sc_cols <- grep("^SC\\.", names(scdata), value = TRUE)
if (length(sc_cols) == 0) stop("No SC.* columns found in input_rds.")

n_edges <- as.integer(if (!is.null(args$n_edges)) args$n_edges else length(sc_cols))
if (is.na(n_edges) || n_edges < 1) n_edges <- length(sc_cols)
n_edges <- min(length(sc_cols), n_edges)

message("[INFO] dataset=", dataset)
message("[INFO] input_rds=", input_rds)
message("[INFO] k_values=", paste(k_values, collapse = ","))
message("[INFO] n_edges=", n_edges)
message("[INFO] bootstrap_n=", bootstrap_n)
message("[INFO] output_dir=", output_dir)
message("[INFO] figure_dir=", figure_dir)
message("[INFO] n_cores=", n_cores)

compute_aic <- function(data, sc_label, k) {
  modelformula <- as.formula(sprintf("%s ~ s(age, k = %s, fx = TRUE) + sex + mean_fd", sc_label, k))
  if (dataset == "abcd") {
    mod <- tryCatch(
      gamm4::gamm4(modelformula, random = ~(1|subID), REML = TRUE, data = data),
      error = function(e) NULL
    )
    if (is.null(mod)) return(NA_real_)
    return(AIC(mod$mer))
  }
  mod <- tryCatch(
    gam(modelformula, method = "REML", data = data),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NA_real_)
  AIC(mod)
}

fit_edge <- function(sc_label) {
  aic_vals <- sapply(k_values, function(k) compute_aic(scdata, sc_label, k))
  names(aic_vals) <- paste0("AIC_k", k_values)
  data.frame(region = sc_label, t(aic_vals), check.names = FALSE, stringsAsFactors = FALSE)
}

result_rows <- mclapply(sc_cols[seq_len(n_edges)], fit_edge, mc.cores = n_cores)
aic_by_edge <- do.call(rbind, result_rows)

aic_cols <- paste0("AIC_k", k_values)
best_k <- apply(aic_by_edge[, aic_cols, drop = FALSE], 1, function(x) {
  if (all(is.na(x))) return(NA_integer_)
  k_values[which.min(x)]
})
aic_by_edge$best_k <- best_k

best_table <- as.data.frame(table(best_k, useNA = "ifany"), stringsAsFactors = FALSE)
colnames(best_table) <- c("k", "n_edges")
best_table$k <- as.character(best_table$k)
valid_total <- sum(best_table$n_edges[best_table$k != "NA"])
best_table$proportion <- if (valid_total > 0) {
  best_table$n_edges / valid_total
} else {
  NA_real_
}

mean_aic <- sapply(aic_cols, function(col) mean(aic_by_edge[[col]], na.rm = TRUE))
mean_aic_df <- data.frame(k = k_values, mean_aic = as.numeric(mean_aic))

aic_summary <- merge(best_table, mean_aic_df, by.x = "k", by.y = "k", all = TRUE)

write.csv(aic_by_edge, file.path(output_dir, "aic_by_edge.csv"), row.names = FALSE)
write.csv(aic_summary, file.path(output_dir, "aic_summary.csv"), row.names = FALSE)

if (all(is.na(unlist(aic_by_edge[, aic_cols, drop = FALSE])))) {
  warning("All AIC values are NA; skipping AIC figure output.")
} else {
  aic_long <- data.frame(
    k = rep(k_values, each = nrow(aic_by_edge)),
    AIC = as.vector(t(as.matrix(aic_by_edge[, aic_cols, drop = FALSE])))
  )

  pdf(file.path(figure_dir, "aic_compare_by_k.pdf"), width = 6, height = 4)
  boxplot(AIC ~ k, data = aic_long, col = "#B4D3E7",
          xlab = "k", ylab = "AIC",
          main = paste0("AIC compare (", toupper(dataset), ")"))
  dev.off()

  tiff(file.path(figure_dir, "aic_compare_by_k.tiff"), width = 1800, height = 1200, res = 300)
  boxplot(AIC ~ k, data = aic_long, col = "#B4D3E7",
          xlab = "k", ylab = "AIC",
          main = paste0("AIC compare (", toupper(dataset), ")"))
  dev.off()
}

if (bootstrap_n > 0) {
  stratify_candidates <- c("siteID", "site", "study", "site_id", "siteid")
  stratify_var <- stratify_candidates[stratify_candidates %in% names(scdata)]
  stratify_var <- if (length(stratify_var) > 0) stratify_var[1] else NA_character_

  message("[INFO] bootstrap stratify_var=", ifelse(is.na(stratify_var), "<none>", stratify_var))

  sample_indices <- function(n) {
    if (is.na(stratify_var)) {
      sample.int(n, n, replace = TRUE)
    } else {
      split_idx <- split(seq_len(n), scdata[[stratify_var]])
      sampled <- lapply(split_idx, function(x) sample(x, length(x), replace = TRUE))
      unlist(sampled, use.names = FALSE)
    }
  }

  bootstrap_best <- mclapply(seq_len(bootstrap_n), function(i) {
    set.seed(925 + i)
    idx <- sample_indices(nrow(scdata))
    data_sub <- scdata[idx, , drop = FALSE]

    edge_best <- rep(NA_integer_, n_edges)
    for (j in seq_len(n_edges)) {
      sc_label <- sc_cols[j]
      aic_vals <- sapply(k_values, function(k) compute_aic(data_sub, sc_label, k))
      if (!all(is.na(aic_vals))) {
        edge_best[j] <- k_values[which.min(aic_vals)]
      }
    }

    if (all(is.na(edge_best))) return(NA_integer_)
    freq <- table(edge_best)
    as.integer(names(freq)[which.max(freq)])
  }, mc.cores = n_cores)

  bootstrap_best <- as.integer(unlist(bootstrap_best))
  counts <- setNames(rep(0L, length(k_values)), k_values)
  freq_table <- table(bootstrap_best, useNA = "no")
  counts[names(freq_table)] <- as.integer(freq_table)
  total_valid <- sum(counts)
  freq <- if (total_valid > 0) counts / total_valid else rep(NA_real_, length(counts))

  bootstrap_df <- data.frame(k = as.integer(names(counts)), n_boot = counts, frequency = freq)
  write.csv(bootstrap_df, file.path(output_dir, "bootstrap_k_counts.csv"), row.names = FALSE)

  if (total_valid > 0) {
    pdf(file.path(figure_dir, "bootstrap_best_k_hist.pdf"), width = 6, height = 4)
    bar_mid <- barplot(freq, names.arg = names(counts), col = "#B4D3E7",
                       xlab = "k", ylab = "Frequency",
                       main = paste0("Bootstrap best k (", toupper(dataset), ")"))
    text(bar_mid, freq, labels = counts, pos = 3, cex = 0.8)
    dev.off()

    tiff(file.path(figure_dir, "bootstrap_best_k_hist.tiff"), width = 1800, height = 1200, res = 300)
    bar_mid <- barplot(freq, names.arg = names(counts), col = "#B4D3E7",
                       xlab = "k", ylab = "Frequency",
                       main = paste0("Bootstrap best k (", toupper(dataset), ")"))
    text(bar_mid, freq, labels = counts, pos = 3, cex = 0.8)
    dev.off()
  } else {
    warning("No valid bootstrap results; skipping bootstrap figure output.")
  }
}
