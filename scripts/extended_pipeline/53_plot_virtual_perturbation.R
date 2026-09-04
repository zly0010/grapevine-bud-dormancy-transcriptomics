options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: 53_plot_virtual_perturbation.R <advanced_project> [output_dir]")
}
project <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(project, "results", "26_virtual_perturbation_2026-08-06", "source_data")
if (!dir.exists(analysis_dir)) analysis_dir <- file.path(project, "04_分析结果", "26_virtual_perturbation", "source_data")
out_root <- if (length(args) == 2) args[[2]] else file.path(project, "figure_build")
fig_dir <- file.path(out_root, "figures")
src_dir <- file.path(out_root, "source_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

COL <- list(
  ink = "#202124", muted = "#6B7280", grid = "#D9DDE3", pale = "#F5F6F7",
  blue = "#356D9A", turquoise = "#159D98", accent = "#D49A2A",
  forest = "#356D9A", sparse = "#B05A47", stable = "#2A8C68", unstable = "#A7ADB4"
)
module_cols <- c(blue = COL$blue, turquoise = COL$turquoise)
method_cols <- c(forest = COL$forest, sparse = COL$sparse)

theme_nature <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL$ink),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = COL$ink),
      axis.line = element_line(linewidth = 0.35, colour = COL$ink),
      axis.ticks = element_line(linewidth = 0.35, colour = COL$ink),
      axis.ticks.length = unit(1.5, "mm"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = COL$muted, hjust = 0),
      plot.margin = margin(5, 6, 5, 6),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.8),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(4, "mm")
    )
}

panel_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 8),
      plot.tag.position = c(0, 1)
    )
}

save_figure <- function(plot, stem, width_mm = 180, height_mm = 170) {
  ggsave(file.path(fig_dir, paste0(stem, ".svg")), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, bg = "white")
  ragg::agg_png(file.path(fig_dir, paste0(stem, ".png")), width = width_mm, height = height_mm,
                units = "mm", res = 300, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = width_mm, height = height_mm,
                 units = "mm", res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

read_tsv <- function(name) {
  read.delim(file.path(analysis_dir, name), sep = "\t", check.names = FALSE, quote = "", comment.char = "")
}
write_source <- function(x, name) {
  write.table(x, file.path(src_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

summary_df <- read_tsv("04_perturbation_consensus_summary.tsv")
fold_df <- read_tsv("09_leave_one_variety_perturbation_summary.tsv")
validation <- read_tsv("06_leave_one_variety_prediction.tsv")
fold_all <- read_tsv("08_leave_one_variety_perturbation.tsv")

as_flag <- function(x) tolower(as.character(x)) %in% c("true", "1", "t")
summary_df <- summary_df %>% mutate(
  robust_computational_hit = as_flag(robust_computational_hit),
  seed_direction_consistent = as_flag(seed_direction_consistent),
  model_direction_concordant = as_flag(model_direction_concordant)
)
fold_df <- fold_df %>% mutate(
  robust_computational_hit = as_flag(robust_computational_hit),
  fold_direction_consistent = as_flag(fold_direction_consistent),
  fold_main_direction_concordant = as_flag(fold_main_direction_concordant),
  leave_one_direction_stable = as_flag(leave_one_direction_stable)
)

dat <- summary_df %>%
  left_join(
    fold_df %>% select(gene_id, scenario, fold_alignment_mean, fold_alignment_min,
                       fold_alignment_max, leave_one_direction_stable),
    by = c("gene_id", "scenario")
  ) %>%
  mutate(
    robust_loo = robust_computational_hit & leave_one_direction_stable,
    intervention = paste(gene_id, scenario),
    label = paste0(gene_id, " (", tf_family, ") ", scenario),
    evidence_class = case_when(
      robust_loo ~ "Robust + leave-one stable",
      robust_computational_hit ~ "Robust, leave-one unstable",
      TRUE ~ "Other tested intervention"
    )
  )

priority_ids <- c(
  "Vitvi08g00761 KO", "Vitvi08g00761 KD50",
  "Vitvi13g00639 KO", "Vitvi13g00639 KD50", "Vitvi13g00639 OE1SD",
  "Vitvi13g00563 KO", "Vitvi13g00563 KD50", "Vitvi05g00053 OE1SD"
)
priority <- dat %>%
  filter(intervention %in% priority_ids, robust_loo) %>%
  mutate(label = factor(label, levels = rev(label[order(global_effect_rms_z)])))

# a: seed and leave-one ranges for the compact priority set
range_long <- bind_rows(
  priority %>% transmute(label, module_color, estimate = forest_alignment_mean,
                         lo = forest_alignment_min, hi = forest_alignment_max,
                         refit = "Three tree seeds"),
  priority %>% transmute(label, module_color, estimate = fold_alignment_mean,
                         lo = fold_alignment_min, hi = fold_alignment_max,
                         refit = "Four leave-one-variety refits")
) %>% mutate(refit = factor(refit, levels = c("Three tree seeds", "Four leave-one-variety refits")))
write_source(range_long, "figure_6a_priority_alignment_ranges.tsv")

p_a <- ggplot(range_long, aes(estimate, label, colour = module_color, shape = refit)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = COL$muted) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 position = position_dodge(width = 0.45), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8, stroke = 0.35) +
  scale_colour_manual(values = module_cols) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Trajectory-alignment score", y = NULL, colour = "TF module", shape = NULL) +
  guides(colour = "none") +
  theme_nature() + theme(legend.position = "top", legend.justification = "left")

# b: independent model comparison for all tested interventions
label_points <- dat %>%
  filter(intervention %in% c("Vitvi08g00761 KO", "Vitvi13g00639 KO",
                             "Vitvi13g00639 OE1SD", "Vitvi13g00563 KO", "Vitvi05g00053 OE1SD"))
rho_models <- cor(dat$forest_alignment_mean, dat$sparse_alignment, method = "spearman")
write_source(dat %>% select(gene_id, tf_family, module_color, scenario, forest_alignment_mean,
                            sparse_alignment, robust_computational_hit, leave_one_direction_stable,
                            global_effect_rms_z, permutation_fdr),
             "figure_6b_model_concordance.tsv")

p_b <- ggplot(dat, aes(forest_alignment_mean, sparse_alignment)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = COL$grid) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = COL$grid) +
  geom_point(aes(fill = evidence_class), shape = 21, size = 1.55, stroke = 0.25,
             colour = "white", alpha = 0.85) +
  ggrepel::geom_text_repel(data = label_points, aes(label = paste0(gene_id, " ", scenario)),
                           size = 2.15, max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.25, box.padding = 0.3) +
  scale_fill_manual(values = c(
    "Robust + leave-one stable" = COL$stable,
    "Robust, leave-one unstable" = COL$accent,
    "Other tested intervention" = COL$unstable
  ), labels = c(
    "Robust + leave-one stable" = "Stable",
    "Robust, leave-one unstable" = "Robust only",
    "Other tested intervention" = "Other"
  )) +
  labs(x = "Tree-ensemble alignment", y = "Sparse-model alignment", fill = NULL) +
  theme_nature() + theme(legend.position = "top", legend.justification = "left",
                         legend.text = element_text(size = 6.2))

# c: module displacement vectors for representative interventions
representative <- dat %>%
  filter(intervention %in% c("Vitvi08g00761 KO", "Vitvi13g00639 KO",
                             "Vitvi13g00639 OE1SD", "Vitvi13g00563 KO", "Vitvi05g00053 OE1SD")) %>%
  mutate(short = paste0(gene_id, " ", scenario))
write_source(representative %>% select(gene_id, tf_family, module_color, scenario,
                                       delta_blue_mean, delta_turquoise_mean,
                                       forest_alignment_mean, global_effect_rms_z),
             "figure_6c_module_displacement.tsv")

p_c <- ggplot(representative, aes(delta_blue_mean, delta_turquoise_mean, colour = module_color)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#E7F3EF", alpha = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = COL$muted) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL$muted) +
  geom_segment(aes(x = 0, y = 0, xend = delta_blue_mean, yend = delta_turquoise_mean),
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"), linewidth = 0.7) +
  ggrepel::geom_text_repel(aes(label = short), size = 2.1, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.25) +
  scale_colour_manual(values = module_cols) +
  labs(x = expression(Delta*" blue module score"),
       y = expression(Delta*" turquoise module score"), colour = "TF module") +
  theme_nature() + theme(legend.position = "none")

# d: cross-variety predictive adequacy and its limitation
validation_long <- validation %>%
  select(method, held_out_variety, delta_spearman, rmse_improvement_fraction) %>%
  pivot_longer(c(delta_spearman, rmse_improvement_fraction), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, delta_spearman = "Delta-expression Spearman rho",
                    rmse_improvement_fraction = "RMSE improvement vs persistence"),
    held_out_variety = factor(held_out_variety, levels = c("Vamu", "Vvcs", "Vvri", "Vrip"))
  )
write_source(validation_long, "figure_6d_leave_one_prediction.tsv")

p_d <- ggplot(validation_long, aes(held_out_variety, value, colour = method, group = method)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = COL$muted) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.8) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = method_cols, labels = c(forest = "Tree ensemble", sparse = "Sparse model")) +
  labs(x = "Held-out variety", y = NULL, colour = NULL) +
  theme_nature() +
  theme(legend.position = "top", legend.justification = "left", strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold"))

figure6 <- (panel_tag(p_a, "a") | panel_tag(p_b, "b")) /
  (wrap_elements(full = panel_tag(p_c, "c")) | panel_tag(p_d, "d")) +
  plot_layout(heights = c(1.05, 1), widths = c(1.1, 1))
save_figure(figure6, "figure_6_virtual_perturbation", 180, 175)

# Supplementary Figure S2: effect-size and robustness boundaries
write_source(dat %>% select(gene_id, tf_family, module_color, scenario, evidence_class,
                            global_effect_rms_z, affected_genes, forest_alignment_mean,
                            sparse_alignment, permutation_fdr),
             "figure_s2a_effect_size_distribution.tsv")
p_s2a <- ggplot(dat, aes(global_effect_rms_z, fill = evidence_class)) +
  geom_histogram(binwidth = 0.004, position = "identity", alpha = 0.72, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.04, linetype = "dashed", linewidth = 0.45, colour = COL$ink) +
  scale_fill_manual(values = c(
    "Robust + leave-one stable" = COL$stable,
    "Robust, leave-one unstable" = COL$accent,
    "Other tested intervention" = COL$unstable
  )) +
  labs(title = "Propagated effects are small", subtitle = "Directly clamped TF excluded from all effect sizes",
       x = "Downstream global RMS effect (z units)", y = "Interventions", fill = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_nature() + theme(
    legend.position = "top", legend.justification = "left",
    legend.box = "vertical", plot.margin = margin(5, 3, 5, 2)
  )

fold_priority <- fold_all %>%
  mutate(intervention = paste(gene_id, scenario)) %>%
  filter(intervention %in% priority_ids) %>%
  select(held_out_variety, gene_id, tf_family, module_color, scenario, trajectory_alignment) %>%
  mutate(label = paste0(gene_id, " ", scenario))
write_source(fold_priority, "figure_s2b_fold_priority_alignment.tsv")
p_s2b <- ggplot(fold_priority, aes(held_out_variety, trajectory_alignment, colour = module_color, group = label)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = COL$muted) +
  geom_line(linewidth = 0.45, alpha = 0.75) + geom_point(size = 1.35) +
  facet_wrap(~label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = module_cols) +
  labs(title = "Priority directions across leave-one refits", x = "Held-out variety",
       y = "Trajectory-alignment score", colour = NULL) +
  theme_nature(7) + theme(legend.position = "none", strip.background = element_blank(),
                           strip.text = element_text(size = 6.2, face = "bold"),
                           axis.text.x = element_text(angle = 35, hjust = 1))

scenario_stable <- dat %>% filter(robust_loo) %>%
  mutate(gene_label = paste0(gene_id, " (", tf_family, ")"))
write_source(scenario_stable %>% select(gene_id, tf_family, module_color, scenario,
                                        forest_alignment_mean, global_effect_rms_z),
             "figure_s2c_stable_scenarios.tsv")
p_s2c <- ggplot(scenario_stable, aes(forest_alignment_mean, reorder(gene_label, forest_alignment_mean),
                                     colour = scenario, size = global_effect_rms_z)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = COL$muted) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = c(KO = "#356D9A", KD50 = "#2A8C68", OE1SD = "#B05A47")) +
  scale_size_continuous(range = c(1.2, 3.2)) +
  labs(title = "Twenty interventions pass all gates", x = "Tree-ensemble alignment",
       y = NULL, colour = "Scenario", size = "RMS z") +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    size = guide_legend(order = 2, nrow = 1)
  ) +
  theme_nature() + theme(
    legend.position = "top", legend.justification = "left", legend.box = "vertical"
  )

gate_df <- data.frame(
  gate = factor(c("Tested", "Tree-seed direction", "Independent-model direction",
                  "Permutation FDR < 0.05", "Top-quartile RMS", "Leave-one direction"),
                levels = rev(c("Tested", "Tree-seed direction", "Independent-model direction",
                               "Permutation FDR < 0.05", "Top-quartile RMS", "Leave-one direction"))),
  n = c(
    nrow(dat),
    sum(dat$seed_direction_consistent),
    sum(dat$seed_direction_consistent & dat$model_direction_concordant),
    sum(dat$seed_direction_consistent & dat$model_direction_concordant & dat$permutation_fdr < 0.05),
    sum(dat$robust_computational_hit),
    sum(dat$robust_loo)
  )
)
write_source(gate_df, "figure_s2d_evidence_gates.tsv")
p_s2d <- ggplot(gate_df, aes(n, gate)) +
  geom_col(width = 0.68, fill = COL$stable) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Evidence gates narrow candidates", x = "Interventions retained", y = NULL) +
  theme_nature()

figure_s2 <- (wrap_elements(full = panel_tag(p_s2a, "a")) | panel_tag(p_s2b, "b")) /
  (panel_tag(p_s2c, "c") | panel_tag(p_s2d, "d")) +
  plot_layout(widths = c(1, 1.15), heights = c(1, 1))
save_figure(figure_s2, "figure_s2_virtual_perturbation_boundaries", 180, 185)

writeLines(c(
  sprintf("Forest/sparse Spearman rho: %.6f", rho_models),
  sprintf("Robust interventions: %d", sum(dat$robust_computational_hit)),
  sprintf("Robust and leave-one direction stable: %d", sum(dat$robust_loo)),
  sprintf("Priority interventions plotted: %d", nrow(priority)),
  "The directly clamped regulator was excluded from all outcome effect sizes."
), file.path(out_root, "FIGURE_SUMMARY.txt"))
