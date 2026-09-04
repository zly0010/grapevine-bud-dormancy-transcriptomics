#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
# Package library paths are supplied by the user environment.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: plot_figure4e_v2.R <results_v2_dir> <figure_output_dir>")
}
results_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty figure output directory: ", output_dir)
}

meta <- read.csv(file.path(results_dir, "stouffer_meta_v2.csv"), check.names = FALSE)
loo <- read.csv(file.path(results_dir, "stouffer_leave_one_dataset_out_v2.csv"), check.names = FALSE)
core_modules <- c("blue", "turquoise")

full <- meta %>%
  filter(module_color %in% core_modules) %>%
  transmute(
    module_color,
    omitted_dataset = "All datasets",
    datasets_retained = datasets_combined,
    combined_z,
    combined_one_sided_p,
    fdr_bh_across_modules,
    analysis = "full"
  )
sens <- loo %>%
  filter(module_color %in% core_modules) %>%
  transmute(
    module_color,
    omitted_dataset = paste("Without", omitted_dataset),
    datasets_retained,
    combined_z,
    combined_one_sided_p,
    fdr_bh_across_modules = NA_real_,
    analysis = "leave_one_dataset_out"
  )
plot_data <- bind_rows(full, sens) %>%
  mutate(
    omitted_dataset = factor(
      omitted_dataset,
      levels = rev(c(
        "All datasets", "Without GSE273240", "Without GSE184114",
        "Without GSE277812", "Without GSE337039"
      ))
    ),
    module_color = factor(module_color, levels = core_modules),
    point_label = sprintf("%.2f", combined_z)
  )

write.csv(plot_data, file.path(output_dir, "Figure4e_v2_source_data.csv"), row.names = FALSE)

COL <- c(blue = "#356D9A", turquoise = "#159D98")
p <- ggplot(plot_data, aes(combined_z, omitted_dataset, colour = module_color)) +
  geom_vline(xintercept = qnorm(0.95), linetype = 2, linewidth = 0.35, colour = "#6B7280") +
  geom_segment(
    aes(x = qnorm(0.95), xend = combined_z, yend = omitted_dataset),
    linewidth = 0.45, alpha = 0.42
  ) +
  geom_point(size = 2.4) +
  geom_text(aes(label = point_label), nudge_x = 0.14, hjust = 0, size = 2.2,
            family = "Arial", show.legend = FALSE) +
  facet_wrap(~module_color, nrow = 1) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(limits = c(1.35, 7.9), breaks = c(1.645, 3, 4, 5, 6, 7),
                     labels = c("1.645", "3", "4", "5", "6", "7"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Leave-one-dataset-out coherence",
    subtitle = "B = 10,000 size-matched random gene sets; weights = sqrt(n)",
    x = "Weighted Stouffer Z",
    y = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    text = element_text(colour = "#202124"),
    axis.title = element_text(size = 7.5),
    axis.text = element_text(size = 7.0, colour = "#202124"),
    axis.line = element_line(linewidth = 0.35, colour = "#202124"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#202124"),
    plot.title = element_text(size = 8.0, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 6.7, colour = "#6B7280", hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7.5),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 6)
  )

stem <- file.path(output_dir, "Figure4e_v2")
width_mm <- 120
height_mm <- 72
ggsave(paste0(stem, ".svg"), p, width = width_mm, height = height_mm,
       units = "mm", device = svglite::svglite, bg = "white")
ggsave(paste0(stem, ".pdf"), p, width = width_mm, height = height_mm,
       units = "mm", device = cairo_pdf, bg = "white")
ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm,
              units = "mm", res = 600, background = "white")
print(p)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm,
               units = "mm", res = 600, compression = "lzw", background = "white")
print(p)
dev.off()

writeLines(
  c(
    "Figure contract:",
    "- Core conclusion: blue and turquoise retain strong combined coherence after any one validation dataset is removed.",
    "- Evidence role: robustness/validation.",
    "- Archetype: quantitative dot plot.",
    "- Backend: R only.",
    "- Statistical threshold: one-sided P = 0.05 (Z = 1.645).",
    "- Source data: Figure4e_v2_source_data.csv.",
    "- Exports: editable SVG/PDF, 600 dpi PNG/TIFF."
  ),
  file.path(output_dir, "Figure4e_v2_QA_NOTES.txt"),
  useBytes = TRUE
)

message("FIGURE4E_V2_OK")
