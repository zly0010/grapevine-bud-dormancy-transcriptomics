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

set.seed(20260731)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: 49_build_figures_part1.R <advanced_project> <discovery_project> [output_dir]")
}
project <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
legacy <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
project_results <- if (dir.exists(file.path(project, "results"))) file.path(project, "results") else file.path(project, "04_分析结果")
legacy_results <- file.path(legacy, "results_corrected")
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

module_cols <- c(
  black = "#282828", blue = COL$blue, brown = COL$brown, green = "#3F8C5A",
  grey = COL$grey, magenta = "#B24C8E", pink = "#DC7C9C", purple = "#7655A3",
  red = "#C4473A", turquoise = COL$turquoise, yellow = "#D2A72C"
)
variety_cols <- c(Vamu = "#C4473A", Vvcs = "#356D9A", Vvri = "#2A8C68", Vrip = "#8A5AA3")

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

theme_heatmap <- function(base_size = 7.5) {
  theme_nature(base_size) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      panel.background = element_blank()
    )
}

panel_tag <- function(p, tag) {
  p + labs(tag = tag) + theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8),
    plot.tag.position = c(0, 1)
  )
}

save_figure <- function(plot, stem, width_mm = 180, height_mm = 180) {
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

read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "", ...)
}

hc_segments <- function(hc) {
  n <- length(hc$order)
  leaf_x <- numeric(n)
  leaf_x[hc$order] <- seq_len(n)
  node_x <- numeric(n - 1)
  node_h <- hc$height
  segs <- vector("list", 3 * (n - 1))
  k <- 1L
  for (i in seq_len(n - 1)) {
    kids <- hc$merge[i, ]
    child_x <- numeric(2)
    child_h <- numeric(2)
    for (j in 1:2) {
      if (kids[j] < 0) {
        child_x[j] <- leaf_x[-kids[j]]
        child_h[j] <- 0
      } else {
        child_x[j] <- node_x[kids[j]]
        child_h[j] <- node_h[kids[j]]
      }
    }
    node_x[i] <- mean(child_x)
    segs[[k]] <- data.frame(x = child_x[1], xend = child_x[1], y = child_h[1], yend = node_h[i]); k <- k + 1L
    segs[[k]] <- data.frame(x = child_x[2], xend = child_x[2], y = child_h[2], yend = node_h[i]); k <- k + 1L
    segs[[k]] <- data.frame(x = child_x[1], xend = child_x[2], y = node_h[i], yend = node_h[i]); k <- k + 1L
  }
  bind_rows(segs)
}

write_source <- function(x, name) {
  write.table(x, file.path(src_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Figure 1: discovery cohort structure and co-expression architecture
design <- read_tsv(file.path(legacy_results, "01_sample_design_and_qc", "sample_design.tsv"))
design <- design %>% filter(tolower(as.character(qc_B)) %in% c("true", "1"))
vst_all <- as.matrix(read_tsv(file.path(legacy_results, "07_wgcna_fixed", "GSE124820_vst_fixed.txt"), row.names = 1))
storage.mode(vst_all) <- "double"
assignments <- read_tsv(file.path(legacy_results, "07_wgcna_fixed", "GSE124820", "module_assignments.txt"))
stopifnot(nrow(assignments) == 5000, all(assignments$gene_id %in% colnames(vst_all)))
vst <- vst_all[, assignments$gene_id, drop = FALSE]
design <- design[match(rownames(vst), design$sample_id), ]
stopifnot(nrow(vst) == 172, ncol(vst) == 5000, all(design$sample_id == rownames(vst)))

pca <- prcomp(vst, center = TRUE, scale. = FALSE)
pca_df <- cbind(design, as.data.frame(pca$x[, 1:3, drop = FALSE]))
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
write_source(pca_df, "figure_1b_discovery_pca.tsv")

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = variety, shape = variety)) +
  geom_point(size = 1.9, alpha = 0.82, stroke = 0.25) +
  scale_colour_manual(values = variety_cols) +
  scale_shape_manual(values = c(Vamu = 16, Vvcs = 17, Vvri = 15, Vrip = 18)) +
  labs(title = "Discovery cohort PCA", x = sprintf("PC1 (%.1f%%)", pve[1]),
       y = sprintf("PC2 (%.1f%%)", pve[2]), colour = NULL, shape = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.4))) + theme_nature() +
  theme(legend.position = "top", legend.justification = "left")

sample_cor <- cor(t(vst), use = "pairwise.complete.obs")
sample_hc <- hclust(as.dist(1 - sample_cor), method = "average")
sample_order <- sample_hc$order
sample_cor_o <- sample_cor[sample_order, sample_order]
sample_long <- as.data.frame(as.table(sample_cor_o))
names(sample_long) <- c("sample_y", "sample_x", "correlation")
sample_long$x <- match(sample_long$sample_x, colnames(sample_cor_o))
sample_long$y <- match(sample_long$sample_y, rownames(sample_cor_o))
write_source(sample_long[, c("sample_x", "sample_y", "correlation")], "figure_1a_sample_correlation.tsv")

p_corr_heat <- ggplot(sample_long, aes(x, y, fill = correlation)) +
  geom_raster() + coord_fixed(expand = FALSE) +
  scale_fill_gradientn(colours = c("#24476D", "#E7ECF0", "#F7F7F4", "#D58B73", "#8D2831"),
                       values = scales::rescale(c(-0.1, 0.35, 0.6, 0.8, 1)), limits = c(-0.1, 1), oob = scales::squish) +
  labs(x = NULL, y = NULL, fill = "Pearson r") + theme_heatmap() +
  theme(axis.text = element_blank(), legend.position = "right", plot.margin = margin(0, 3, 4, 3))

sample_seg <- hc_segments(sample_hc)
p_corr_tree <- ggplot(sample_seg) +
  geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = 0.18, colour = COL$ink, lineend = "square") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Sample correlation and clustering") + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold", colour = COL$ink,
                                  margin = margin(l = 14)),
        plot.margin = margin(4, 3, 0, 3))

p_corr <- p_corr_tree / p_corr_heat + plot_layout(heights = c(0.75, 4))

assignments <- assignments[match(colnames(vst), assignments$gene_id), ]
stopifnot(all(assignments$gene_id == colnames(vst)))

gene_hc <- hclust(as.dist(1 - cor(vst, use = "pairwise.complete.obs")), method = "average")
gene_seg <- hc_segments(gene_hc)
gene_bar <- data.frame(
  x = seq_along(gene_hc$order), y = 0,
  module_color = assignments$module_color[gene_hc$order],
  gene_id = assignments$gene_id[gene_hc$order]
)
write_source(gene_bar, "figure_1c_gene_dendrogram_order.tsv")

p_gene_tree <- ggplot(gene_seg) +
  geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = 0.08, colour = "#34373A", lineend = "square") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Co-expression gene dendrogram") + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold", colour = COL$ink,
                                  margin = margin(l = 14)),
        plot.margin = margin(4, 3, 0, 3))

p_gene_bar <- ggplot(gene_bar, aes(x, y, fill = module_color)) +
  geom_tile(width = 1.02, height = 1) +
  scale_fill_manual(values = module_cols, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) + theme_void() +
  theme(plot.margin = margin(0, 3, 4, 3))

p_gene <- p_gene_tree / p_gene_bar + plot_layout(heights = c(5.2, 0.38))

trait_df <- read_tsv(file.path(legacy_results, "07_wgcna_fixed", "GSE124820", "module_trait_correlation.txt"))
trait_names <- c("variety_Vamu", "variety_Vvcs", "variety_Vvri", "variety_Vrip", "time")
trait_long <- trait_df %>%
  select(module, all_of(trait_names), all_of(paste0("p_BH.", trait_names))) %>%
  pivot_longer(cols = all_of(trait_names), names_to = "trait", values_to = "correlation") %>%
  mutate(fdr = unlist(lapply(seq_len(n()), function(i) trait_df[match(module[i], trait_df$module), paste0("p_BH.", trait[i])])),
         module = sub("^ME", "", module),
         trait = recode(trait, variety_Vamu = "V. amurensis", variety_Vvcs = "V. vinifera CS",
                        variety_Vvri = "V. vinifera Riesling", variety_Vrip = "V. riparia", time = "Time"),
         label = sprintf("%.2f%s", correlation, ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "")))))
trait_long$module <- factor(trait_long$module, levels = rev(c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "grey")))
trait_long$trait <- factor(trait_long$trait, levels = c("V. amurensis", "V. vinifera CS", "V. vinifera Riesling", "V. riparia", "Time"))
write_source(trait_long, "figure_1d_module_trait_correlations.tsv")

p_trait <- ggplot(trait_long, aes(trait, module, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.05, family = "Arial", colour = COL$ink) +
  scale_fill_gradient2(low = "#356D9A", mid = "#F7F7F4", high = "#C4473A", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Module–trait associations", x = NULL, y = NULL, fill = "Pearson r") +
  theme_heatmap() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(face = "plain"), legend.position = "right")

fig1 <- (panel_tag(wrap_elements(full = p_corr), "a") | panel_tag(p_pca, "b")) /
  (panel_tag(wrap_elements(full = p_gene), "c") | panel_tag(p_trait, "d")) +
  plot_layout(widths = c(1.08, 0.92), heights = c(1.06, 0.94))

save_figure(fig1, "Figure_1_discovery_structure_and_WGCNA_R", 180, 177)

# Figure 2: cross-variety differential expression and meta-analysis
lrt_paths <- c(
  Vamu = file.path(legacy_results, "02_deseq2_gse124820", "Vamu_B_no_fail", "LRT_Vamu_B.txt"),
  Vvcs = file.path(legacy_results, "02_deseq2_gse124820", "Vvcs_B_no_fail", "LRT_Vvcs_B.txt"),
  Vvri = file.path(legacy_results, "02_deseq2_gse124820", "Vvri_B_no_fail", "LRT_Vvri_B.txt"),
  Vrip = file.path(legacy_results, "02_deseq2_gse124820", "Vrip_B_no_fail", "LRT_Vrip_B.txt")
)
deg_sets <- lapply(lrt_paths, function(p) {
  x <- read_tsv(p)
  x$gene_id[!is.na(x$padj) & x$padj < 0.05]
})
all_deg <- sort(unique(unlist(deg_sets)))
membership <- sapply(deg_sets, function(s) all_deg %in% s)
pattern <- apply(membership, 1, function(z) paste(as.integer(z), collapse = ""))
pattern_counts <- sort(table(pattern), decreasing = TRUE)
top_patterns <- names(pattern_counts)[seq_len(min(12, length(pattern_counts)))]
upset_bars <- data.frame(pattern = factor(top_patterns, levels = top_patterns), count = as.integer(pattern_counts[top_patterns])) %>%
  mutate(label = ifelse(row_number() <= 4, scales::comma(count), ""))
upset_matrix <- bind_rows(lapply(seq_along(top_patterns), function(i) {
  bits <- as.integer(strsplit(top_patterns[i], "", fixed = TRUE)[[1]])
  data.frame(pattern = factor(top_patterns[i], levels = top_patterns), variety = factor(names(deg_sets), levels = rev(names(deg_sets))), member = bits)
}))
write_source(cbind(gene_id = all_deg, as.data.frame(membership)), "figure_2a_deg_set_membership.tsv")

p_up_bar <- ggplot(upset_bars, aes(pattern, count)) +
  geom_col(width = 0.66, fill = COL$ink) +
  geom_text(aes(label = label), vjust = -0.25, size = 2.05, family = "Arial") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Time-responsive gene intersections", x = NULL, y = "Intersection size") +
  theme_nature() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 5, 0, 5))

p_up_mat <- ggplot(upset_matrix, aes(pattern, variety)) +
  geom_point(colour = "#D7DADF", size = 1.65) +
  geom_line(data = upset_matrix %>% filter(member == 1) %>% group_by(pattern) %>% filter(n() > 1),
            aes(group = pattern), linewidth = 0.45, colour = COL$ink) +
  geom_point(data = upset_matrix %>% filter(member == 1), colour = COL$ink, size = 1.8) +
  labs(x = NULL, y = NULL) + theme_nature() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0, 5, 4, 5))

p_upset <- p_up_bar / p_up_mat + plot_layout(heights = c(2.5, 1.25))

meta <- read_tsv(file.path(project_results, "10_cross_variety_meta", "02_gene_reml_meta.tsv"))
effects <- read_tsv(file.path(project_results, "10_cross_variety_meta", "01_per_variety_effects.tsv"))
meta <- meta %>% left_join(assignments[, c("gene_id", "module_color")], by = "gene_id") %>%
  mutate(class = case_when(
    robust_meta_gene & module_color == "blue" ~ "Blue module",
    robust_meta_gene & module_color == "turquoise" ~ "Turquoise module",
    robust_meta_gene ~ "Other robust",
    TRUE ~ "Not robust"
  ), neglog10_fdr = -log10(pmax(padj, 1e-300)))
write_source(meta, "figure_2c_meta_analysis_genes.tsv")

p_volcano <- ggplot(meta, aes(estimate, neglog10_fdr)) +
  geom_point(data = meta %>% filter(class == "Not robust"), colour = "#B9BDC2", size = 0.65, alpha = 0.38) +
  geom_point(data = meta %>% filter(class == "Other robust"), colour = "#555B63", size = 0.8, alpha = 0.65) +
  geom_point(data = meta %>% filter(class %in% c("Blue module", "Turquoise module")),
             aes(colour = class), size = 0.9, alpha = 0.78) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = COL$muted) +
  scale_colour_manual(values = c("Blue module" = COL$blue, "Turquoise module" = COL$turquoise)) +
  coord_cartesian(xlim = quantile(meta$estimate, c(0.005, 0.995), na.rm = TRUE)) +
  labs(title = "Cross-variety random-effects meta-analysis", x = "Pooled log2 fold change", y = expression(-log[10](FDR)), colour = NULL) +
  theme_nature() + theme(legend.position = "top", legend.justification = "left")

target <- meta %>% filter(robust_meta_gene, module_color %in% c("blue", "turquoise")) %>%
  group_by(module_color) %>% slice_max(abs(estimate), n = 25, with_ties = FALSE) %>% ungroup()
heat <- effects %>% filter(gene_id %in% target$gene_id) %>%
  select(gene_id, variety, log2_fold_change) %>%
  pivot_wider(names_from = variety, values_from = log2_fold_change) %>%
  left_join(target[, c("gene_id", "estimate", "module_color")], by = "gene_id")
heat_mat <- as.matrix(heat[, c("Vamu", "Vvcs", "Vvri", "Vrip", "estimate")])
rownames(heat_mat) <- heat$gene_id
heat_order <- hclust(dist(heat_mat))$order
heat <- heat[heat_order, ]
heat$y <- seq_len(nrow(heat))
heat_long <- heat %>%
  select(gene_id, y, module_color, Vamu, Vvcs, Vvri, Vrip, estimate) %>%
  pivot_longer(cols = c(Vamu, Vvcs, Vvri, Vrip, estimate), names_to = "source", values_to = "effect") %>%
  mutate(source = recode(source, estimate = "Meta"), source = factor(source, levels = c("Vamu", "Vvcs", "Vvri", "Vrip", "Meta")))
write_source(heat_long, "figure_2b_cross_variety_effect_heatmap.tsv")

heat_labels <- heat %>%
  group_by(module_color) %>%
  summarise(y = mean(y), .groups = "drop")

p_heat_strip <- ggplot(heat, aes(1, y, fill = module_color)) +
  geom_tile() +
  geom_text(data = heat_labels, aes(x = 1, y = y, label = module_color),
            inherit.aes = FALSE, angle = 90, colour = "white", family = "Arial",
            fontface = "bold", size = 1.75) +
  scale_fill_manual(values = module_cols, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) + theme_void()

p_effect_heat <- ggplot(heat_long, aes(source, y, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#356D9A", mid = "#F7F7F4", high = "#C4473A", midpoint = 0,
                       limits = c(-8, 8), oob = scales::squish) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Concordant effects across varieties", x = NULL, y = NULL, fill = "log2 FC") +
  theme_heatmap() + theme(axis.text.y = element_blank(), legend.position = "right")

p_heat <- p_heat_strip + p_effect_heat + plot_layout(widths = c(0.17, 1))

module_meta <- meta %>% filter(!is.na(module_color), module_color != "grey") %>%
  group_by(module_color) %>%
  summarise(n = n(), median_effect = median(estimate, na.rm = TRUE), robust_fraction = mean(robust_meta_gene, na.rm = TRUE), .groups = "drop")
boot_ci <- lapply(split(meta$estimate[!is.na(meta$module_color) & meta$module_color != "grey"],
                        meta$module_color[!is.na(meta$module_color) & meta$module_color != "grey"]), function(x) {
  q <- quantile(replicate(1000, median(sample(x, replace = TRUE), na.rm = TRUE)), c(0.025, 0.975), na.rm = TRUE)
  data.frame(ci_low = q[1], ci_high = q[2])
})
boot_ci <- bind_rows(boot_ci, .id = "module_color")
module_meta <- left_join(module_meta, boot_ci, by = "module_color") %>% arrange(median_effect)
module_meta$module_color <- factor(module_meta$module_color, levels = module_meta$module_color)
write_source(module_meta, "figure_2d_module_meta_effects.tsv")

p_module <- ggplot(module_meta, aes(median_effect, module_color)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = COL$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.55, colour = COL$ink) +
  geom_point(aes(size = robust_fraction, fill = module_color), shape = 21, colour = "white", stroke = 0.35) +
  scale_fill_manual(values = module_cols, guide = "none") +
  scale_size_continuous(range = c(2, 5), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Module-level pooled effects", x = "Median pooled log2 fold change", y = NULL, size = "Robust genes") +
  theme_nature() + theme(legend.position = "right")

i2_df <- meta %>% filter(!is.na(i2)) %>%
  mutate(group = case_when(module_color == "blue" ~ "Blue", module_color == "turquoise" ~ "Turquoise", TRUE ~ "All genes"))
write_source(i2_df[, c("gene_id", "i2", "group", "robust_meta_gene")], "figure_2e_heterogeneity.tsv")
p_i2 <- ggplot(i2_df, aes(i2, colour = group)) +
  stat_ecdf(linewidth = 0.75) +
  geom_vline(xintercept = 0.75, linetype = 2, linewidth = 0.35, colour = COL$muted) +
  scale_colour_manual(values = c("All genes" = "#777C82", Blue = COL$blue, Turquoise = COL$turquoise)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Between-variety heterogeneity", x = expression(I^2), y = "Cumulative fraction", colour = NULL) +
  theme_nature() + theme(legend.position = "top", legend.justification = "left")

fig2 <- (panel_tag(wrap_elements(full = p_upset), "a") | panel_tag(p_volcano, "b")) /
  (panel_tag(wrap_elements(full = p_heat), "c") |
     (panel_tag(p_module, "d") / panel_tag(p_i2, "e"))) +
  plot_layout(widths = c(0.92, 1.08), heights = c(0.86, 1.14))

save_figure(fig2, "Figure_2_cross_variety_meta_analysis_R", 180, 190)

writeLines(c(
  "R_NATURE_FIGURES_PART1_OK",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("R:", R.version.string),
  paste("ggplot2:", as.character(packageVersion("ggplot2"))),
  paste("WGCNA input:", nrow(vst), "samples x", ncol(vst), "genes")
), file.path(out_root, "BUILD_STATUS.txt"))
