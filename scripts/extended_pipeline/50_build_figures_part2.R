# WARNING: The Figure 4 statistics exported by this historical plotting workflow
# are not final manuscript inference. Use frozen method_revision_v2 sources and
# scripts/method_revision_v2/plot_figure4_full_v2.R for final Figure 4.
# The retained Figure 3, Figure 5 and Supplementary Figure S1 sources remain documented.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
})

set.seed(20260731)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: 50_build_figures_part2.R <advanced_project> <discovery_project> [output_dir]")
}
project <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
legacy <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
project_results <- if (dir.exists(file.path(project, "results"))) file.path(project, "results") else file.path(project, "04_分析结果")
legacy_results <- file.path(legacy, "results_corrected")
package_source <- file.path(project, "03_图件源数据")
out_root <- if (length(args) == 3) args[[3]] else file.path(project, "figure_build")
fig_dir <- file.path(out_root, "figures")
src_dir <- file.path(out_root, "source_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

COL <- list(
  ink = "#202124", muted = "#6B7280", grid = "#D9DDE3", pale = "#F5F6F7",
  blue = "#356D9A", turquoise = "#159D98", brown = "#9A6A3A", grey = "#A7ADB4",
  up = "#C4473A", down = "#356D9A", accent = "#D49A2A"
)
module_cols <- c(blue = COL$blue, turquoise = COL$turquoise, brown = COL$brown, grey = COL$grey)
variety_cols <- c(Vamu = "#C4473A", Vvcs = "#356D9A", Vvri = "#2A8C68", Vrip = "#8A5AA3")
condition_cols <- c(Controlled_4C = "#356D9A", Natural = "#C4473A")
cultivar_cols <- c(Brianna = "#8A5AA3", Marquette = "#2A8C68")

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

save_figure <- function(plot, stem, width_mm = 180, height_mm = 180) {
  ggsave(file.path(fig_dir, paste0(stem, ".svg")), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, bg = "white")
  ragg::agg_png(file.path(fig_dir, paste0(stem, ".png")), width = width_mm, height = height_mm,
                units = "mm", res = 300, background = "white")
  print(plot); dev.off()
  ragg::agg_tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = width_mm, height = height_mm,
                 units = "mm", res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
}

read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "", ...)
}

write_source <- function(x, name) {
  write.table(x, file.path(src_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

hc_segments <- function(hc) {
  n <- length(hc$order)
  leaf_x <- numeric(n); leaf_x[hc$order] <- seq_len(n)
  node_x <- numeric(n - 1); node_h <- hc$height
  segs <- vector("list", 3 * (n - 1)); k <- 1L
  for (i in seq_len(n - 1)) {
    kids <- hc$merge[i, ]; child_x <- child_h <- numeric(2)
    for (j in 1:2) {
      if (kids[j] < 0) { child_x[j] <- leaf_x[-kids[j]]; child_h[j] <- 0 }
      else { child_x[j] <- node_x[kids[j]]; child_h[j] <- node_h[kids[j]] }
    }
    node_x[i] <- mean(child_x)
    segs[[k]] <- data.frame(x = child_x[1], xend = child_x[1], y = child_h[1], yend = node_h[i]); k <- k + 1L
    segs[[k]] <- data.frame(x = child_x[2], xend = child_x[2], y = child_h[2], yend = node_h[i]); k <- k + 1L
    segs[[k]] <- data.frame(x = child_x[1], xend = child_x[2], y = node_h[i], yend = node_h[i]); k <- k + 1L
  }
  bind_rows(segs)
}

tag_theme <- theme(
  plot.background = element_rect(fill = "white", colour = NA),
  plot.tag = element_text(family = "Arial", face = "bold", size = 8)
)

panel_tag <- function(p, tag) {
  p + labs(tag = tag) + theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8),
    plot.tag.position = c(0, 1)
  )
}

# Figure 3: resampling stability and biological trajectories
stab_dir <- file.path(project_results, "12_wgcna_stability")
boot <- read_tsv(file.path(stab_dir, "02_bootstrap_metrics.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise", "brown")) %>%
  mutate(abs_time_correlation = abs(time_correlation))
boot_long <- boot %>%
  select(replicate_id, module_color, mean_pairwise_correlation, median_abs_kme,
         pc1_variance_explained, abs_time_correlation) %>%
  pivot_longer(cols = -c(replicate_id, module_color), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    mean_pairwise_correlation = "Mean pairwise r",
    median_abs_kme = "Median |kME|",
    pc1_variance_explained = "PC1 variance",
    abs_time_correlation = "|Time correlation|"
  ))
boot_summary <- boot_long %>% group_by(module_color, metric) %>%
  summarise(median = median(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), .groups = "drop")
boot_summary$metric <- factor(boot_summary$metric, levels = c("Mean pairwise r", "Median |kME|", "PC1 variance", "|Time correlation|"))
boot_summary$module_color <- factor(boot_summary$module_color, levels = c("brown", "blue", "turquoise"))
write_source(boot_summary, "figure_3a_bootstrap_stability_heatmap.tsv")

p_boot <- ggplot(boot_summary, aes(metric, module_color, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", median)), size = 2.35, family = "Arial") +
  scale_fill_gradientn(colours = c("#F1F3F4", "#A9C4D4", "#356D9A"), limits = c(0.45, 0.9), oob = scales::squish) +
  labs(title = "Bootstrap module stability", x = NULL, y = NULL, fill = "Median") +
  theme_heatmap() + theme(axis.text.x = element_text(angle = 32, hjust = 1), legend.position = "right")

loo_map <- read_tsv(file.path(stab_dir, "06_leave_one_variety_module_mapping.tsv")) %>%
  filter(reference_module %in% c("blue", "turquoise", "brown")) %>%
  mutate(reference_module = factor(reference_module, levels = c("brown", "blue", "turquoise")),
         omitted_variety = factor(omitted_variety, levels = c("Vamu", "Vvcs", "Vvri", "Vrip")))
write_source(loo_map, "figure_3b_leave_one_variety_module_matching.tsv")

p_loo <- ggplot(loo_map, aes(omitted_variety, reference_module, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.4, family = "Arial") +
  scale_fill_gradientn(colours = c("#EFE5DA", "#F7F7F4", "#159D98"), values = scales::rescale(c(0.25, 0.5, 0.95)),
                       limits = c(0.25, 0.95), oob = scales::squish) +
  labs(title = "Leave-one-variety module recovery", subtitle = "Cell values are Jaccard indices", x = "Omitted variety", y = NULL, fill = "Jaccard") +
  theme_heatmap() + theme(legend.position = "right")

design <- read_tsv(file.path(legacy_results, "01_sample_design_and_qc", "sample_design.tsv")) %>%
  filter(tolower(as.character(qc_B)) %in% c("true", "1"))
eig <- read_tsv(file.path(legacy_results, "07_wgcna_fixed", "GSE124820", "module_eigengenes.txt"), row.names = 1)
eig$sample_id <- rownames(eig)
traj <- eig %>%
  select(sample_id, MEblue, MEturquoise, MEbrown) %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module_color", values_to = "eigengene") %>%
  mutate(module_color = sub("^ME", "", module_color)) %>%
  left_join(design[, c("sample_id", "variety", "time")], by = "sample_id") %>%
  group_by(module_color, variety, time) %>%
  summarise(mean = mean(eigengene), sem = sd(eigengene) / sqrt(n()), n = n(), .groups = "drop") %>%
  mutate(module_color = factor(module_color, levels = c("blue", "turquoise", "brown")))
write_source(traj, "figure_3c_module_eigengene_trajectories.tsv")

p_traj <- ggplot(traj, aes(time, mean, colour = variety, group = variety)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = COL$grid) +
  geom_ribbon(data = traj %>% filter(is.finite(sem)),
              aes(ymin = mean - sem, ymax = mean + sem, fill = variety),
              colour = NA, alpha = 0.10) +
  geom_line(linewidth = 0.75) + geom_point(size = 1.3, stroke = 0) +
  facet_wrap(~module_color, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = variety_cols) + scale_fill_manual(values = variety_cols) +
  labs(title = "Module eigengene trajectories across four grape materials", x = "Time (days)", y = "Module eigengene", colour = NULL, fill = NULL) +
  theme_nature() + theme(legend.position = "top", legend.justification = "left", panel.spacing = grid::unit(4, "mm"))

fig3 <- (panel_tag(p_boot, "a") | panel_tag(p_loo, "b")) /
  panel_tag(p_traj, "c") + plot_layout(heights = c(0.8, 1.2))
save_figure(fig3, "Figure_3_module_stability_and_trajectories_R", 180, 150)

# Figure 4: cross-dataset whole-module validation
val_dir <- file.path(project_results, "17_gse337039_independent_validation", "module_validation_four_datasets")
coh <- read_tsv(file.path(val_dir, "04_whole_module_coherence.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise")) %>%
  mutate(dataset = factor(dataset, levels = c("GSE273240", "GSE184114", "GSE277812", "GSE337039")),
         module_color = factor(module_color, levels = c("blue", "turquoise")))
write_source(coh, "figure_4a_whole_module_coherence.tsv")

p_coh <- ggplot(coh, aes(dataset, module_color, fill = observed_mean_pairwise_correlation)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("r = %.2f\n%.0f%% genes", observed_mean_pairwise_correlation, 100 * gene_coverage)),
            size = 2.2, lineheight = 0.9, family = "Arial") +
  scale_fill_gradientn(colours = c("#F1F3F4", "#A9C4D4", "#159D98"), limits = c(0, 0.6), oob = scales::squish) +
  labs(title = "Whole-module coherence across datasets", x = NULL, y = NULL, fill = "Mean r") +
  theme_heatmap() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

null <- read_tsv(file.path(val_dir, "05_coherence_null_permutations.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise"))
write_source(null, "figure_4b_coherence_null_distributions.tsv")
p_null <- ggplot(null, aes(dataset, null_mean_pairwise_correlation, fill = module_color)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.4, alpha = 0.5, position = position_dodge(width = 0.66)) +
  geom_point(data = coh, aes(x = dataset, y = observed_mean_pairwise_correlation, colour = module_color),
             shape = 21, fill = "white", size = 2.2, stroke = 0.75, position = position_dodge(width = 0.66), inherit.aes = FALSE) +
  scale_fill_manual(values = module_cols) + scale_colour_manual(values = module_cols) +
  labs(title = "Observed coherence exceeds matched null modules", x = NULL, y = "Mean pairwise correlation", fill = NULL, colour = NULL) +
  theme_nature() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top", legend.justification = "left")

scores <- read_tsv(file.path(val_dir, "03_independent_projected_scores.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise"))
score_order <- scores %>% filter(module_color == "blue") %>% group_by(dataset) %>% arrange(projected_module_score, .by_group = TRUE) %>%
  mutate(sample_index = row_number()) %>% ungroup() %>% select(dataset, sample_id, sample_index)
scores <- left_join(scores, score_order, by = c("dataset", "sample_id")) %>%
  mutate(module_color = factor(module_color, levels = c("turquoise", "blue")),
         dataset_short = factor(sub("^GSE", "", dataset), levels = c("273240", "184114", "277812", "337039")))
write_source(scores, "figure_4c_projected_module_score_heatmap.tsv")
p_scores <- ggplot(scores, aes(sample_index, module_color, fill = projected_module_score)) +
  geom_tile() + facet_grid(. ~ dataset_short, scales = "free_x") +
  scale_fill_gradient2(low = "#356D9A", mid = "#F7F7F4", high = "#C4473A", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(title = "Projected module scores in independent samples", x = "Samples ordered within dataset", y = NULL, fill = "Score") +
  theme_heatmap() + theme(axis.text.x = element_blank(), strip.text = element_text(size = 6.6), legend.position = "right")

effects <- read_tsv(file.path(val_dir, "06_projected_score_effects.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise"), !is.na(partial_eta_squared), fdr < 0.05) %>%
  group_by(dataset, module_color) %>% slice_max(partial_eta_squared, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(effect_short = case_when(
    grepl(":", effect) ~ "interaction",
    grepl("time|day", effect, ignore.case = TRUE) ~ "time",
    grepl("condition|treatment|phase", effect, ignore.case = TRUE) ~ "condition",
    TRUE ~ "design effect"
  ), module_color = factor(module_color, levels = c("blue", "turquoise")))
write_source(effects, "figure_4d_strongest_context_effects.tsv")
p_effects <- ggplot(effects, aes(dataset, module_color)) +
  geom_tile(fill = "#F3F4F5", colour = "white", linewidth = 0.5) +
  geom_point(aes(size = partial_eta_squared, fill = -log10(fdr)), shape = 21, colour = "white", stroke = 0.35) +
  geom_text(aes(label = effect_short), nudge_y = -0.27, size = 1.8, family = "Arial", colour = COL$ink) +
  scale_fill_gradient(low = "#E7ECEF", high = "#C4473A",
                      guide = guide_colorbar(barheight = grid::unit(12, "mm"),
                                             barwidth = grid::unit(2.4, "mm"),
                                             title.position = "top")) +
  scale_size_continuous(range = c(3, 7), limits = c(0, 1), breaks = c(0.5, 1.0)) +
  labs(title = "Strongest significant context effect", x = NULL, y = NULL, size = expression(partial~eta^2), fill = expression(-log[10](FDR))) +
  theme_heatmap() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right", legend.box = "vertical",
        legend.spacing.y = grid::unit(1.2, "mm"),
        legend.key.height = grid::unit(3.2, "mm"))

lodo <- read_tsv(file.path(val_dir, "07_leave_one_dataset_out_coherence.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise")) %>%
  mutate(omitted_dataset = ifelse(omitted_dataset == "NONE", "All datasets", paste("Without", omitted_dataset)),
         omitted_dataset = factor(omitted_dataset, levels = rev(c("All datasets", "Without GSE273240", "Without GSE184114", "Without GSE277812", "Without GSE337039"))),
         module_color = factor(module_color, levels = c("blue", "turquoise")))
write_source(lodo, "figure_4e_leave_one_dataset_out.tsv")
p_lodo <- ggplot(lodo, aes(weighted_stouffer_z, omitted_dataset, colour = module_color)) +
  geom_vline(xintercept = 1.645, linetype = 2, linewidth = 0.4, colour = COL$muted) +
  geom_point(size = 2.5) + facet_wrap(~module_color, nrow = 1) +
  scale_colour_manual(values = module_cols) +
  labs(title = "Leave-one-dataset-out coherence", x = "Weighted Stouffer Z", y = NULL, colour = NULL) +
  theme_nature() + theme(legend.position = "none")

fig4 <- (panel_tag(p_coh, "a") | panel_tag(p_null, "b")) /
  (panel_tag(p_scores, "c") | panel_tag(p_effects, "d")) /
  panel_tag(p_lodo, "e") + plot_layout(heights = c(0.86, 0.82, 0.72))
save_figure(fig4, "Figure_4_cross_dataset_module_validation_R", 180, 210)

# Figure 5: independent FASTQ validation in GSE337039
gse_dir <- file.path(project_results, "17_gse337039_independent_validation")
vst337 <- as.matrix(read_tsv(file.path(gse_dir, "04_vst_samples_by_vitvi_genes.tsv"), row.names = 1))
storage.mode(vst337) <- "double"
meta337 <- read_tsv(file.path(gse_dir, "03_pca_scores.tsv"))
meta337 <- meta337[match(rownames(vst337), meta337$run_accession), ]
stopifnot(nrow(vst337) == 60, all(rownames(vst337) == meta337$run_accession))

cor337 <- cor(t(vst337), use = "pairwise.complete.obs")
hc337 <- hclust(as.dist(1 - cor337), method = "average")
ord337 <- hc337$order
cor_o <- cor337[ord337, ord337]
cor_long <- as.data.frame(as.table(cor_o)); names(cor_long) <- c("sample_y", "sample_x", "correlation")
cor_long$x <- match(cor_long$sample_x, colnames(cor_o)); cor_long$y <- match(cor_long$sample_y, rownames(cor_o))
write_source(cor_long, "figure_5a_gse337039_sample_correlation.tsv")
p337_heat <- ggplot(cor_long, aes(x, y, fill = correlation)) +
  geom_raster() + coord_fixed(expand = FALSE) +
  scale_fill_gradientn(colours = c("#24476D", "#E7ECF0", "#F7F7F4", "#D58B73", "#8D2831"),
                       values = scales::rescale(c(0, 0.45, 0.65, 0.82, 1)), limits = c(0, 1), oob = scales::squish) +
  labs(x = NULL, y = NULL, fill = "Pearson r") + theme_heatmap() +
  theme(axis.text = element_blank(), plot.margin = margin(0, 0, 4, 3))
p337_tree <- ggplot(hc_segments(hc337)) +
  geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = 0.22, colour = COL$ink) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(title = "GSE337039 sample clustering") + theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8, face = "bold", margin = margin(l = 0)),
    plot.margin = margin(4, 14.5, 0, 14.5)
  )

# Reserve a dedicated legend column so the dendrogram and heatmap data panels
# share exactly the same x extent.  Without this, the heatmap legend shrinks
# only the lower panel and makes the tree appear wider than the heatmap.
p337_cluster <- p337_tree + p337_heat + guide_area() +
  plot_layout(
    design = "
AC
BC
",
    heights = c(0.75, 4), widths = c(1, 0.22), guides = "collect"
  ) & theme(legend.position = "right")

write_source(meta337, "figure_5b_gse337039_pca.tsv")
p337_pca <- ggplot(meta337, aes(PC1, PC2, colour = cultivar, shape = condition)) +
  geom_point(size = 2, alpha = 0.86, stroke = 0.35) +
  scale_colour_manual(values = cultivar_cols) + scale_shape_manual(values = c(Controlled_4C = 16, Natural = 17)) +
  labs(title = "Independent cohort PCA", x = "PC1", y = "PC2", colour = NULL, shape = NULL) +
  guides(shape = guide_legend(order = 1, nrow = 1), colour = guide_legend(order = 2, nrow = 1)) +
  theme_nature() + theme(legend.position = "top", legend.justification = "left", legend.box = "vertical")

three <- read_tsv(file.path(gse_dir, "01_three_way_interaction_LRT.tsv"))
bri <- read_tsv(file.path(gse_dir, "lrt_condition_by_time__Brianna.tsv"))
mar <- read_tsv(file.path(gse_dir, "lrt_condition_by_time__Marquette.tsv"))
top_genes <- three %>% filter(!is.na(padj), padj < 0.05, vitvi_id %in% colnames(vst337)) %>%
  arrange(padj) %>% slice_head(n = 60) %>% pull(vitvi_id)
z <- t(scale(t(t(vst337[, top_genes, drop = FALSE]))))
z[!is.finite(z)] <- 0
sample_order <- meta337 %>% arrange(cultivar, condition, time_index, replicate) %>% pull(run_accession)
z <- z[, sample_order, drop = FALSE]
hm <- as.data.frame(as.table(z)); names(hm) <- c("gene_id", "sample_id", "z")
hm <- hm %>% left_join(meta337[, c("run_accession", "cultivar", "condition", "time_index")], by = c("sample_id" = "run_accession")) %>%
  mutate(group = paste(cultivar, ifelse(condition == "Controlled_4C", "4C", "Nat.")))
write_source(hm, "figure_5c_three_way_interaction_heatmap.tsv")
p337_degheat <- ggplot(hm, aes(sample_id, gene_id, fill = z)) +
  geom_tile() + facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#356D9A", mid = "#F7F7F4", high = "#C4473A", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(title = "Top three-way interaction genes", x = "Samples ordered by time", y = "60 genes", fill = "Row z-score") +
  theme_heatmap() +
  theme(
    axis.text = element_blank(), strip.text = element_text(size = 5.8),
    legend.position = "right", axis.title.y = element_text(margin = margin(r = 2)),
    plot.title = element_text(margin = margin(l = 14)),
    plot.margin = margin(5, 3, 5, 1)
  )

deg_sets <- list(
  `Cultivar × condition × time` = three$vitvi_id[!is.na(three$padj) & three$padj < 0.05],
  `Brianna condition × time` = bri$vitvi_id[!is.na(bri$padj) & bri$padj < 0.05],
  `Marquette condition × time` = mar$vitvi_id[!is.na(mar$padj) & mar$padj < 0.05]
)
all_deg <- sort(unique(unlist(deg_sets)))
mem <- sapply(deg_sets, function(s) all_deg %in% s)
pat <- apply(mem, 1, function(v) paste(as.integer(v), collapse = ""))
pc <- sort(table(pat), decreasing = TRUE)
tp <- names(pc)[seq_len(min(7, length(pc)))]
bars <- data.frame(pattern = factor(tp, levels = tp), count = as.integer(pc[tp]))
matp <- bind_rows(lapply(seq_along(tp), function(i) data.frame(
  pattern = factor(tp[i], levels = tp), set = factor(names(deg_sets), levels = rev(names(deg_sets))),
  member = as.integer(strsplit(tp[i], "", fixed = TRUE)[[1]])
)))
write_source(cbind(gene_id = all_deg, as.data.frame(mem)), "figure_5d_deg_set_membership.tsv")
p5bar <- ggplot(bars, aes(pattern, count)) + geom_col(width = 0.62, fill = COL$ink) +
  geom_text(aes(label = scales::comma(count)), vjust = -0.25, size = 2.1, family = "Arial") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "DEG intersections", x = NULL, y = "Genes") + theme_nature() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 5, 0, 5))
p5mat <- ggplot(matp, aes(pattern, set)) + geom_point(colour = "#D7DADF", size = 1.7) +
  geom_line(data = matp %>% filter(member == 1) %>% group_by(pattern) %>% filter(n() > 1), aes(group = pattern), linewidth = 0.45, colour = COL$ink) +
  geom_point(data = matp %>% filter(member == 1), colour = COL$ink, size = 1.9) +
  labs(x = NULL, y = NULL) + theme_nature() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0, 5, 4, 5))
p337_upset <- p5bar / p5mat + plot_layout(heights = c(2.2, 1.15))

timecourse_path <- file.path(project_results, "18_submission_figures", "figure_8_gse337039", "source_data", "panel_f_module_timecourse.tsv")
if (!file.exists(timecourse_path)) timecourse_path <- file.path(package_source, "figure_5e_module_timecourse.tsv")
timecourse <- read_tsv(timecourse_path) %>%
  filter(module_color %in% c("blue", "turquoise")) %>%
  mutate(module_color = factor(module_color, levels = c("blue", "turquoise")))
write_source(timecourse, "figure_5e_module_timecourse.tsv")
p337_time <- ggplot(timecourse, aes(time_index, mean, colour = condition, group = condition)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = COL$grid) +
  geom_ribbon(data = timecourse %>% filter(is.finite(sem)),
              aes(ymin = mean - sem, ymax = mean + sem, fill = condition),
              alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.45) +
  facet_grid(module_color ~ cultivar, switch = "y") + scale_colour_manual(values = condition_cols) + scale_fill_manual(values = condition_cols) +
  labs(title = "Independent module trajectories", x = "Time index", y = "Projected module score", colour = NULL, fill = NULL) +
  theme_nature() + theme(legend.position = "top", legend.justification = "left", strip.placement = "outside",
                         strip.text.y.left = element_text(angle = 0, hjust = 1))

interact <- read_tsv(file.path(val_dir, "06_projected_score_effects.tsv")) %>%
  filter(dataset == "GSE337039", module_color %in% c("blue", "turquoise"), grepl(":", effect)) %>%
  mutate(cultivar = ifelse(grepl("Brianna", analysis_family), "Brianna", "Marquette"),
         label = sprintf("FDR %.1e", fdr),
         row = paste(cultivar, module_color, sep = " · "))
write_source(interact, "figure_5f_module_interaction_effects.tsv")
p337_interact <- ggplot(interact, aes(partial_eta_squared, reorder(row, partial_eta_squared), colour = module_color)) +
  geom_point(size = 2.8) + geom_text(aes(label = label), nudge_x = -0.015, hjust = 1, size = 2, family = "Arial", colour = COL$ink) +
  scale_colour_manual(values = module_cols) + coord_cartesian(xlim = c(0.72, 1.01)) +
  labs(title = "Condition × time effects", x = expression(partial~eta^2), y = NULL, colour = NULL) +
  theme_nature() + theme(legend.position = "none", panel.grid.major.x = element_line(colour = COL$grid, linewidth = 0.3))

fig5 <- (panel_tag(wrap_elements(full = p337_cluster), "a") | panel_tag(p337_pca, "b")) /
  (panel_tag(wrap_elements(full = p337_degheat), "c") |
     panel_tag(wrap_elements(full = p337_upset), "d")) /
  (panel_tag(p337_time, "e") | panel_tag(p337_interact, "f")) +
  plot_layout(heights = c(0.88, 0.92, 1.12), widths = c(0.95, 1.05))
save_figure(fig5, "Figure_5_GSE337039_independent_validation_R", 180, 235)

# Supplementary Figure S1: annotation scope and exploratory brown module
go_cov <- read_tsv(file.path(legacy_results, "08_wgcna_postprocessing", "01_go_mapping_audit", "module_go_coverage_audit.tsv"))
go_cov_long <- go_cov %>%
  select(module_color, module_genes, successfully_mapped_genes, go_annotated_genes) %>%
  pivot_longer(cols = -module_color, names_to = "stage", values_to = "genes") %>%
  mutate(stage = recode(stage, module_genes = "Module genes", successfully_mapped_genes = "Mapped IDs", go_annotated_genes = "GO annotated"),
         stage = factor(stage, levels = c("Module genes", "Mapped IDs", "GO annotated")))
write_source(go_cov_long, "figure_S1a_go_annotation_coverage.tsv")
p_go_cov <- ggplot(go_cov_long, aes(stage, genes, group = module_color, colour = module_color)) +
  geom_line(linewidth = 0.75) + geom_point(size = 2) +
  scale_colour_manual(values = module_cols) + scale_y_continuous(labels = scales::comma) +
  labs(title = "Annotation attrition by module", x = NULL, y = "Genes", colour = NULL) +
  theme_nature() + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

go <- read_tsv(file.path(legacy_results, "08_wgcna_postprocessing", "02_go_enrichment", "go_enrichment_all_tests.tsv"))
go <- go %>% filter(module_color %in% c("blue", "turquoise", "brown"), is.finite(p_value), p_value > 0) %>%
  group_by(module_color) %>% arrange(p_value, .by_group = TRUE) %>% mutate(rank = row_number(), expected = (rank - 0.5) / n()) %>% ungroup()
write_source(go, "figure_S1b_go_enrichment_tests.tsv")
p_go_qq <- ggplot(go, aes(-log10(expected), -log10(p_value), colour = module_color)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = 2, colour = COL$muted) +
  geom_point(size = 0.75, alpha = 0.55) + facet_wrap(~module_color, nrow = 1) +
  scale_colour_manual(values = module_cols, guide = "none") +
  labs(title = "GO enrichment P-value distributions", x = expression(Expected~-log[10](P)), y = expression(Observed~-log[10](P))) +
  theme_nature()

p_go_scatter <- ggplot(go, aes(fold_enrichment, -log10(p_value), colour = module_color, shape = ontology)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35, colour = COL$muted) +
  geom_point(size = 1, alpha = 0.58) +
  scale_colour_manual(values = module_cols) + coord_cartesian(xlim = c(0, quantile(go$fold_enrichment, 0.99, na.rm = TRUE))) +
  labs(title = "Nominal enrichment without FDR support", subtitle = "0 of 4,032 tests reached FDR < 0.05",
       x = "Fold enrichment", y = expression(-log[10](P)), colour = NULL, shape = "Ontology") +
  theme_nature() + theme(legend.position = "top", legend.justification = "left")

stab_sum <- read_tsv(file.path(stab_dir, "07_module_stability_summary.tsv")) %>%
  filter(module_color %in% c("blue", "turquoise", "brown"))
write_source(stab_sum, "figure_S1d_brown_module_stability.tsv")
p_brown <- ggplot(stab_sum, aes(median_jaccard, reorder(module_color, median_jaccard), colour = module_color)) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.4, colour = COL$muted) +
  geom_segment(aes(x = minimum_jaccard, xend = median_jaccard, yend = reorder(module_color, median_jaccard)), linewidth = 1.0) +
  geom_point(aes(x = minimum_jaccard), shape = 21, fill = "white", size = 2.1, stroke = 0.7) +
  geom_point(size = 2.5) + scale_colour_manual(values = module_cols, guide = "none") +
  scale_x_continuous(limits = c(0.2, 1), breaks = seq(0.2, 1, 0.2)) +
  labs(title = "Brown module fails the stability criterion", subtitle = "Open: minimum; filled: median Jaccard", x = "Leave-one-variety Jaccard", y = NULL) +
  theme_nature()

figs1 <- (panel_tag(p_go_cov, "a") | panel_tag(p_go_qq, "b")) /
  (panel_tag(p_go_scatter, "c") | panel_tag(p_brown, "d")) +
  plot_layout(heights = c(0.92, 1.08))
save_figure(figs1, "Supplementary_Figure_S1_annotation_and_brown_module_R", 180, 155)

writeLines(c(
  "R_NATURE_FIGURES_PART2_OK",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Figures: 3, 4, 5, Supplementary S1"
), file.path(out_root, "BUILD_STATUS_PART2.txt"))
