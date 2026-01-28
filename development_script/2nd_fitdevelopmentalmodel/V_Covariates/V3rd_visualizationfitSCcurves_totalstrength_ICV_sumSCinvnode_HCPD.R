## HCP-D (ComBat-GAM SA12) | Sensitivity analysis: add SES / ICV covariate
##
## Generate developmental trajectory plots from scaled GAM models.
## Outputs are aligned in style to HCP-D Yeo7 figures (sizes/colors/themes).
##
## Inputs:
##   outputs/intermediate/2nd_fitdevelopmentalmodel/hcpd/covariates/<cov_tag>/combat_gam/CV75/
## Outputs:
##   outputs/figures/2nd_fitdevelopmentalmodel/hcpd/covariates/<cov_tag>/combat_gam/CV75/
##
## NOTE: Despite the historical filename, this script supports both SES and ICV via --cov_tag.

rm(list = ls())

library(parallel)
library(tidyverse)
library(ggplot2)

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
force <- as.integer(if (!is.null(args$force)) args$force else 0L) == 1L

CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
cov_tag <- if (!is.null(args$cov_tag)) args$cov_tag else "SES"

ds.resolution <- 12
elementnum <- ds.resolution * (ds.resolution + 1) / 2

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "covariates", cov_tag, "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(
  project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel",
  "hcpd", "covariates", cov_tag, "combat_gam", paste0("CV", CVthr)
)

FigureFolder_SCfit <- file.path(FigureRoot, "SA12_sumSCinvnode_fit")
FigureFolder_SCdecile <- file.path(FigureRoot, "SA12_decile_sumSCinvnode_fit")
dir.create(FigureFolder_SCfit, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder_SCdecile, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "plotdata_generate.R"))
source(file.path(functionFolder, "colorbarvalue.R"))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores) || n_cores < 1) n_cores <- parallel::detectCores()
n_cores <- max(1L, n_cores)

gamresultsum <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
gammodelsum <- readRDS(file.path(interfileFolder, paste0("gammodel", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
if (nrow(gamresultsum) == 0 || length(gammodelsum) == 0) stop("Empty inputs under: ", interfileFolder)

out_plotdatasum <- file.path(interfileFolder, "plotdatasum_scale_TRUE_SA12.rds")
out_sa12_csv <- file.path(interfileFolder, "SA12_10.csv")

out_fig1_tiff <- file.path(FigureFolder_SCfit, "devcurve_Rsq_fit.ratio.tiff")
out_fig1_pdf <- file.path(FigureFolder_SCfit, "devcurve_Rsq_fit.ratio.pdf")
out_fig2_tiff <- file.path(FigureFolder_SCfit, "devcurve_meanderv2_fit.Z.tiff")
out_fig2_pdf <- file.path(FigureFolder_SCfit, "devcurve_meanderv2_fit.Z.pdf")
out_fig3_tiff <- file.path(FigureFolder_SCdecile, "devcurve_SCrank_fit.Z_SCtype10.tiff")
out_fig3_pdf <- file.path(FigureFolder_SCdecile, "devcurve_SCrank_fit.Z_SCtype10.pdf")

if (!force && file.exists(out_fig1_tiff) && file.exists(out_fig1_pdf) &&
  file.exists(out_fig2_tiff) && file.exists(out_fig2_pdf) &&
  file.exists(out_fig3_tiff) && file.exists(out_fig3_pdf)) {
  message("[INFO] S3 outputs exist; skipping. Set --force=1 to re-run.")
  quit(save = "no", status = 0)
}

plot_one <- function(idx) {
  modobj <- gammodelsum[[idx]]
  if (is.null(modobj)) return(NULL)
  tryCatch(
    plotdata_generate(modobj, "age"),
    error = function(e) {
      message("[WARN] plotdata_generate failed for idx=", idx, ": ", conditionMessage(e))
      NULL
    }
  )
}

plotdatasum <- mclapply(seq_len(min(length(gammodelsum), nrow(gamresultsum))), plot_one, mc.cores = n_cores)
saveRDS(plotdatasum, out_plotdatasum)

Matrix12 <- matrix(NA, nrow = ds.resolution, ncol = ds.resolution)
indexup12 <- upper.tri(Matrix12)
indexsave12 <- !indexup12
Matrix12.SCrank <- Matrix12
for (x in 1:ds.resolution) {
  for (y in 1:ds.resolution) {
    Matrix12.SCrank[x, y] <- x^2 + y^2
  }
}
Matrix12.SCrank[indexup12] <- NA
Matrix12.SCrank[indexsave12] <- rank(Matrix12.SCrank[indexsave12], ties.method = "average")

parcel_all <- paste0("SC.", seq_len(elementnum), "_h")
SCrank_map <- setNames(Matrix12.SCrank[indexsave12], parcel_all)

ok_plot_idx <- which(!vapply(plotdatasum, is.null, logical(1)))
if (length(ok_plot_idx) == 0) stop("All plotdata_generate() calls failed; see warnings above.")

plotdatasum.df <- dplyr::bind_rows(lapply(ok_plot_idx, function(i) {
  tmp <- as.data.frame(plotdatasum[[i]])
  sc_label <- gamresultsum$parcel[[i]]
  if (!is.null(sc_label) && sc_label %in% names(tmp)) {
    tmp <- tmp[, setdiff(names(tmp), sc_label), drop = FALSE]
  }
  tmp$SC_label <- sc_label
  tmp$SCrank <- unname(SCrank_map[[sc_label]])
  tmp$PartialRsq <- as.numeric(gamresultsum$partialRsq[[i]])
  tmp$meanderv2 <- as.numeric(gamresultsum$meanderv2[[i]])
  tmp
}))

## Save SA12 decile labels
SA12_10 <- data.frame(SCrank = Matrix12.SCrank[indexsave12]) %>%
  mutate(decile = ntile(SCrank, 10))
SA12_10$SC_label <- parcel_all
write.csv(SA12_10, out_sa12_csv, row.names = FALSE)

## Plot 1: ratio trajectories colored by partial R^2 (Yeo7 style)
lmthr <- max(abs(plotdatasum.df$PartialRsq), na.rm = TRUE)
p1 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.ratio, group = SC_label, color = PartialRsq), linewidth = 1.5, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-lmthr, lmthr), guide = "none") +
  labs(x = "Age (years)", y = "SC strength (ratio)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 22, hjust = 0.5),
    legend.position = "none"
  )
ggsave(out_fig1_tiff, p1, width = 15, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig1_pdf, p1, dpi = 600, width = 15, height = 12, units = "cm", bg = "transparent")

## Plot 2: z-scored trajectories colored by mean 2nd derivative (Yeo7-like)
colorbarvalues.meanderv2 <- colorbarvalues(
  plotdatasum.df$meanderv2,
  abs(min(plotdatasum.df$meanderv2, na.rm = TRUE)) /
    (max(plotdatasum.df$meanderv2, na.rm = TRUE) - min(plotdatasum.df$meanderv2, na.rm = TRUE))
)
SC_label_derv2_order <- gamresultsum$parcel[order(gamresultsum$meanderv2)]
plotdatasum.df$SC_label2 <- factor(plotdatasum.df$SC_label, levels = SC_label_derv2_order)
p2 <- ggplot() +
  geom_line(data = plotdatasum.df, aes(x = age, y = fit.Z, group = SC_label2, color = meanderv2), linewidth = 0.8, alpha = 0.8) +
  scale_color_distiller(type = "seq", palette = "RdBu", values = colorbarvalues.meanderv2, direction = -1, guide = "none") +
  labs(x = "Age (years)", y = "SC strength (z-score)") +
  scale_y_continuous(breaks = c(-1.5, 0.0, 1.5)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 23, color = "black"),
    axis.title = element_text(size = 23),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )
ggsave(out_fig2_tiff, p2, width = 15, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig2_pdf, p2, dpi = 600, width = 16, height = 16, units = "cm", bg = "transparent")

## Plot 3: decile-average trajectories (fixed output count; no svg)
plotdatasum.df.label <- merge(plotdatasum.df, SA12_10, by = "SC_label", all.x = TRUE)
plotdatasum.df.decile <- plotdatasum.df.label %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(fit, na.rm = TRUE), SCranktype_order = mean(decile, na.rm = TRUE), .groups = "drop")
plotdatasum.df.decile <- plotdatasum.df.decile %>%
  group_by(decile) %>%
  mutate(fit.Z = scale(fit.avg))

p3 <- ggplot(data = plotdatasum.df.decile, aes(x = age, y = fit.Z, group = decile, color = decile)) +
  geom_line(linewidth = 1.3, alpha = 1) +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, guide = "none") +
  labs(x = "Age (years)", y = "SC strength (z-score)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 22.5, color = "black"),
    axis.title = element_text(size = 22.5, color = "black"),
    aspect.ratio = 0.9,
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none"
  )

ggsave(out_fig3_tiff, p3, dpi = 600, width = 15, height = 14, units = "cm", bg = "transparent")
ggsave(out_fig3_pdf, p3, dpi = 600, width = 16, height = 16, units = "cm", bg = "transparent")

