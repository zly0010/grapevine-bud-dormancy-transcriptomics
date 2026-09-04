#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
# Package library paths are supplied by the user environment.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: plot_figure4_full_v2.R <results_v2_dir> <legacy_figure_source_dir> <output_dir>")
}
results_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
source_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

COL <- list(
  ink = "#202124", muted = "#6B7280", grid = "#D9DDE3", pale = "#F5F6F7",
  blue = "#356D9A", turquoise = "#159D98", brown = "#9A6A3A", grey = "#A7ADB4",
  up = "#C4473A", down = "#356D9A", accent = "#D49A2A"
)
module_cols <- c(blue = COL$blue, turquoise = COL$turquoise, brown = COL$brown, grey = COL$grey)

theme_nature <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL$ink),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = COL$ink),
      axis.line = element_line(linewidth = 0.35, colour = COL$ink),
      axis.ticks = element_line(linewidth = 0.35, colour = COL$ink),
      axis.ticks.length = grid::unit(1.5, "mm"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = COL$muted, hjust = 0),
      plot.margin = margin(5, 6, 5, 6),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.8)
    )
}
theme_heatmap <- function(base_size = 7.5) {
  theme_nature(base_size) + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.background = element_blank())
}
panel_tag <- function(p, tag) {
  p + labs(tag = tag) + theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8),
    plot.tag.position = c(0, 1)
  )
}

coh <- read.csv(file.path(results_dir, "module_randomset_validation_v2.csv"), check.names = FALSE) %>%
  filter(module_color %in% c("blue", "turquoise")) %>%
  mutate(
    dataset = factor(dataset, levels = c("GSE273240", "GSE184114", "GSE277812", "GSE337039")),
    module_color = factor(module_color, levels = c("blue", "turquoise"))
  )
p_coh <- ggplot(coh, aes(dataset, module_color, fill = observed_mean_pairwise_correlation)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("r = %.2f\n%.0f%% genes", observed_mean_pairwise_correlation, 100 * gene_coverage)),
            size = 2.2, lineheight = 0.9, family = "Arial") +
  scale_fill_gradientn(colours = c("#F1F3F4", "#A9C4D4", "#159D98"), limits = c(0, 0.6), oob = scales::squish) +
  labs(title = "Whole-module coherence across datasets", x = NULL, y = NULL, fill = "Mean r") +
  theme_heatmap() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

null <- read.csv(file.path(results_dir, "module_randomset_null_full_v2.csv"), check.names = FALSE) %>%
  filter(module_color %in% c("blue", "turquoise")) %>%
  mutate(dataset = factor(dataset, levels = c("GSE273240", "GSE184114", "GSE277812", "GSE337039")))
p_null <- ggplot(null, aes(dataset, null_mean_pairwise_correlation, fill = module_color)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.4, alpha = 0.5,
               position = position_dodge(width = 0.66)) +
  geom_point(data = coh, aes(x = dataset, y = observed_mean_pairwise_correlation, colour = module_color),
             shape = 21, fill = "white", size = 2.2, stroke = 0.75,
             position = position_dodge(width = 0.66), inherit.aes = FALSE) +
  scale_fill_manual(values = module_cols) + scale_colour_manual(values = module_cols) +
  labs(title = "Observed coherence exceeds 10,000 matched null sets", x = NULL,
       y = "Mean pairwise correlation", fill = NULL, colour = NULL) +
  theme_nature() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                         legend.position = "top", legend.justification = "left")

scores <- read.delim(file.path(source_dir, "figure_4c_projected_module_score_heatmap.tsv"), check.names = FALSE) %>%
  mutate(module_color = factor(module_color, levels = c("turquoise", "blue")),
         dataset_short = factor(dataset_short, levels = c("273240", "184114", "277812", "337039")))
p_scores <- ggplot(scores, aes(sample_index, module_color, fill = projected_module_score)) +
  geom_tile() + facet_grid(. ~ dataset_short, scales = "free_x") +
  scale_fill_gradient2(low = "#356D9A", mid = "#F7F7F4", high = "#C4473A", midpoint = 0,
                       limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(title = "Projected module scores in independent samples", x = "Samples ordered within dataset",
       y = NULL, fill = "Score") +
  theme_heatmap() + theme(axis.text.x = element_blank(), strip.text = element_text(size = 6.6), legend.position = "right")

effects <- read.delim(file.path(source_dir, "figure_4d_strongest_context_effects.tsv"), check.names = FALSE) %>%
  mutate(module_color = factor(module_color, levels = c("blue", "turquoise")))
p_effects <- ggplot(effects, aes(dataset, module_color)) +
  geom_tile(fill = "#F3F4F5", colour = "white", linewidth = 0.5) +
  geom_point(aes(size = partial_eta_squared, fill = -log10(fdr)), shape = 21, colour = "white", stroke = 0.35) +
  geom_text(aes(label = effect_short), nudge_y = -0.27, size = 1.8, family = "Arial", colour = COL$ink) +
  scale_fill_gradient(low = "#E7ECEF", high = "#C4473A",
                      guide = guide_colorbar(barheight = grid::unit(12, "mm"),
                                             barwidth = grid::unit(2.4, "mm"), title.position = "top")) +
  scale_size_continuous(range = c(3, 7), limits = c(0, 1), breaks = c(0.5, 1.0)) +
  labs(title = "Strongest significant context effect", x = NULL, y = NULL,
       size = expression(partial~eta^2), fill = expression(-log[10](FDR))) +
  theme_heatmap() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right",
                          legend.box = "vertical", legend.spacing.y = grid::unit(1.2, "mm"),
                          legend.key.height = grid::unit(3.2, "mm"))

meta <- read.csv(file.path(results_dir, "stouffer_meta_v2.csv"), check.names = FALSE)
loo <- read.csv(file.path(results_dir, "stouffer_leave_one_dataset_out_v2.csv"), check.names = FALSE)
lodo <- bind_rows(
  meta %>% filter(module_color %in% c("blue", "turquoise")) %>%
    transmute(module_color, omitted_dataset = "All datasets", weighted_stouffer_z = combined_z),
  loo %>% filter(module_color %in% c("blue", "turquoise")) %>%
    transmute(module_color, omitted_dataset = paste("Without", omitted_dataset), weighted_stouffer_z = combined_z)
) %>%
  mutate(
    omitted_dataset = factor(omitted_dataset, levels = rev(c(
      "All datasets", "Without GSE273240", "Without GSE184114", "Without GSE277812", "Without GSE337039"
    ))),
    module_color = factor(module_color, levels = c("blue", "turquoise"))
  )
p_lodo <- ggplot(lodo, aes(weighted_stouffer_z, omitted_dataset, colour = module_color)) +
  geom_vline(xintercept = 1.645, linetype = 2, linewidth = 0.4, colour = COL$muted) +
  geom_point(size = 2.5) + facet_wrap(~module_color, nrow = 1) +
  scale_colour_manual(values = module_cols) +
  scale_x_continuous(limits = c(1.4, 7.7), breaks = c(1.645, 3, 4, 5, 6, 7),
                     labels = c("1.645", "3", "4", "5", "6", "7")) +
  labs(title = "Leave-one-dataset-out coherence", subtitle = "B = 10,000; weights = sqrt(n)",
       x = "Weighted Stouffer Z", y = NULL, colour = NULL) +
  theme_nature() + theme(legend.position = "none")

fig4 <- (panel_tag(p_coh, "a") | panel_tag(p_null, "b")) /
  (panel_tag(p_scores, "c") | panel_tag(p_effects, "d")) /
  panel_tag(p_lodo, "e") + plot_layout(heights = c(0.86, 0.82, 0.72))

stem <- file.path(output_dir, "Figure_4_cross_dataset_module_validation_v2")
ggsave(paste0(stem, ".svg"), fig4, width = 180, height = 210, units = "mm", device = svglite::svglite, bg = "white")
ggsave(paste0(stem, ".pdf"), fig4, width = 180, height = 210, units = "mm", device = cairo_pdf, bg = "white")
ragg::agg_png(paste0(stem, ".png"), width = 180, height = 210, units = "mm", res = 300, background = "white")
print(fig4)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = 180, height = 210, units = "mm", res = 600,
               compression = "lzw", background = "white")
print(fig4)
dev.off()

message("FIGURE4_FULL_V2_OK")
