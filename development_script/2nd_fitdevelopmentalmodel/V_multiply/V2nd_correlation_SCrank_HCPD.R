## HCP-D (ComBat-GAM) | S-A rank multiply | Scatter plot (meanderv2)
## Runnable version aligned to project-relative paths.

rm(list = ls())

library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
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

is_windows <- .Platform$OS.type == "windows"
CVthr <- as.numeric(if (!is.null(args$cvthr)) args$cvthr else 75)
ds.resolution <- as.integer(if (!is.null(args$ds_res)) args$ds_res else 12L)
elementnum <- ds.resolution * (ds.resolution + 1) / 2

interfileFolder <- file.path(
  project_root, "outputs", "intermediate", "2nd_fitdevelopmentalmodel",
  "hcpd", "combat_gam", paste0("CV", CVthr)
)
FigureRoot <- file.path(project_root, "outputs", "figures", "2nd_fitdevelopmentalmodel", "hcpd", "combat_gam", paste0("CV", CVthr))
FigCorrFolder <- file.path(FigureRoot, "correlation_sumSCinvnode_SCrank")
dir.create(FigureRoot, showWarnings = FALSE, recursive = TRUE)
dir.create(FigCorrFolder, showWarnings = FALSE, recursive = TRUE)

functionFolder <- file.path(project_root, "gamfunction")
source(file.path(functionFolder, "SCrankcorr_mult.R"))

input_rds <- if (!is.null(args$input_rds)) {
  args$input_rds
} else {
  file.path(
    project_root, "outputs", "results", "combat_gam", "hcpd",
    "SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
  )
}

gamresult <- readRDS(file.path(interfileFolder, paste0("gamresults", elementnum, "_sumSCinvnode_over8_CV", CVthr, "_scale_TRUE.rds")))
SCdata <- readRDS(input_rds)

## convert critical ages of insignificantly developmental edges to NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)
gamresult$increase.onset[gamresult$sig == FALSE] <- NA
gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2, 2) == 8.08] <- NA
gamresult$increase.offset[gamresult$sig == FALSE] <- NA
gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2, 2) == 21.92] <- NA
gamresult$peak.change[gamresult$sig == FALSE] <- NA
gamresult$peak.increase.change[gamresult$sig == FALSE] <- NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")

## scatter plot (meanderv2)
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata = TRUE)

mytheme <- theme(
  axis.text = element_text(size = 23, color = "black"),
  axis.title = element_text(size = 23),
  aspect.ratio = 0.9,
  axis.line = element_line(linewidth = 0.6),
  axis.ticks = element_line(linewidth = 0.6),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.position = "none"
)

p <- ggplot(data = correlation.df) +
  geom_point(aes(x = SCrank, y = meanderv2, color = SCrank), size = 5) +
  geom_smooth(aes(x = SCrank, y = meanderv2), linewidth = 2, method = "lm", color = "black") +
  scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, guide = "none") +
  labs(x = "S-A connectional axis rank", y = "Second derivative") +
  scale_y_continuous(breaks = c(-0.003, 0, 0.003), labels = c(-3, 0, 3)) +
  theme_classic() + mytheme

ggsave(file.path(FigCorrFolder, paste0("mean", computevar, "_SCrankcorr_n", ds.resolution, "_SAmult.tiff")),
  p, width = 13, height = 12, units = "cm", bg = "transparent"
)
ggsave(file.path(FigCorrFolder, paste0("mean", computevar, "_SCrankcorr_n", ds.resolution, "_SAmult.pdf")),
  p, dpi = 600, width = 13, height = 12, units = "cm", bg = "transparent"
)
if (is_windows) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("[WARN] svglite not available; skip svg output on Windows.")
  } else {
    ggsave(file.path(FigCorrFolder, paste0("mean", computevar, "_SCrankcorr_n", ds.resolution, "_SAmult.svg")),
      p, dpi = 600, width = 20, height = 14, units = "cm", bg = "transparent"
    )
  }
}
