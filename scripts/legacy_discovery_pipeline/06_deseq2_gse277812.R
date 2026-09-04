#!/usr/bin/env Rscript
# GSE277812 DESeq2 Analysis
# design = ~ stage * node
# Tests: stage effect, node effect, stage x node interaction

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(apeglm)
    library(ashr)
})

cat("========================================\n")
cat("GSE277812 DESeq2 Analysis (~ stage * node)\n")
cat("========================================\n\n")

args <- commandArgs(trailingOnly = TRUE)
base_dir <- normalizePath(
    if (length(args) >= 1L) args[[1]] else Sys.getenv("GRAPE_LEGACY_PROJECT", "."),
    winslash = "/", mustWork = TRUE
)
results_dir <- file.path(base_dir, "results_corrected/05_deseq2_gse277812")
figures_dir <- file.path(base_dir, "results_corrected/figures")
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(figures_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================
# 1. Load data and build design
# ============================================
cat("[1] Loading data...\n")

counts <- read.table(file.path(base_dir, "data/processed/GSE277812/counts_matrix.txt"),
                     header=TRUE, row.names=1, sep="\t", check.names=FALSE)
cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Parse sample names: T{stage}_{node}_{rep}
sample_names <- colnames(counts)
design <- do.call(rbind, lapply(sample_names, function(name) {
    parts <- strsplit(name, "_")[[1]]
    data.frame(
        sample_id = name,
        stage = parts[1],
        node = as.integer(parts[2]),
        replicate = as.integer(parts[3]),
        stringsAsFactors = FALSE
    )
}))

design$stage_factor <- factor(design$stage, levels=c("T1", "T2", "T3"))
design$node_factor <- factor(design$node, levels=c(2, 5, 10))
design$group <- paste0(design$stage, "_Node", design$node)

cat("  Stages:", paste(levels(design$stage_factor), collapse=", "), "\n")
cat("  Nodes:", paste(levels(design$node_factor), collapse=", "), "\n")
cat("  Samples per group:\n")
for (g in sort(unique(design$group))) {
    cat(sprintf("    %s: %d\n", g, sum(design$group == g)))
}

# Save design
write.table(design, file.path(results_dir, "sample_design_gse277812.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ============================================
# 2. Filter low-expression genes
# ============================================
cat("\n[2] Filtering...\n")

counts_mat <- data.matrix(counts)
min_samples <- round(ncol(counts_mat) * 0.1)
keep <- rowSums(counts_mat >= 10) >= min_samples
counts_f <- counts_mat[keep, ]
cat("  Genes after filtering:", nrow(counts_f), "/", nrow(counts_mat), "\n")

# ============================================
# 3. DESeq2 with interaction model
# ============================================
cat("\n[3] Running DESeq2 (~ stage_factor * node_factor)...\n")

dds <- DESeqDataSetFromMatrix(
    countData = counts_f,
    colData = design,
    design = ~ stage_factor * node_factor
)

# Full DESeq2
dds <- DESeq(dds, quiet=TRUE)
all_rn <- resultsNames(dds)
cat("  ResultsNames:\n")
for (rn in all_rn) cat("    ", rn, "\n")

# ============================================
# 4. LRT tests
# ============================================
cat("\n[4] LRT tests...\n")

# 4a. Overall model LRT
cat("  Overall model: full = ~ stage * node, reduced = ~ 1\n")
dds_overall <- DESeq(dds, test="LRT", reduced=~ 1, quiet=TRUE)
res_overall <- results(dds_overall)
n_overall <- sum(res_overall$padj < 0.05, na.rm=TRUE)
cat("  Overall LRT sig:", n_overall, "genes\n")

# Save overall LRT
overall_df <- as.data.frame(res_overall)
overall_df$gene_id <- rownames(overall_df)
overall_df <- overall_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(overall_df, file.path(results_dir, "LRT_overall_model.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# 4b. Stage main effect
cat("  Stage effect: full = ~ stage * node, reduced = ~ node\n")
dds_stage <- DESeq(dds, test="LRT", reduced=~ node_factor, quiet=TRUE)
res_stage <- results(dds_stage)
n_stage <- sum(res_stage$padj < 0.05, na.rm=TRUE)
cat("  Stage LRT sig:", n_stage, "genes\n")

stage_df <- as.data.frame(res_stage)
stage_df$gene_id <- rownames(stage_df)
stage_df <- stage_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(stage_df, file.path(results_dir, "LRT_stage_effect.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# 4c. Node main effect
cat("  Node effect: full = ~ stage * node, reduced = ~ stage\n")
dds_node <- DESeq(dds, test="LRT", reduced=~ stage_factor, quiet=TRUE)
res_node <- results(dds_node)
n_node <- sum(res_node$padj < 0.05, na.rm=TRUE)
cat("  Node LRT sig:", n_node, "genes\n")

node_df <- as.data.frame(res_node)
node_df$gene_id <- rownames(node_df)
node_df <- node_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(node_df, file.path(results_dir, "LRT_node_effect.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# 4d. Interaction effect
cat("  Interaction: full = ~ stage * node, reduced = ~ stage + node\n")
dds_int <- DESeq(dds, test="LRT", reduced=~ stage_factor + node_factor, quiet=TRUE)
res_int <- results(dds_int)
n_int <- sum(res_int$padj < 0.05, na.rm=TRUE)
cat("  Interaction LRT sig:", n_int, "genes\n")

int_df <- as.data.frame(res_int)
int_df$gene_id <- rownames(int_df)
int_df <- int_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(int_df, file.path(results_dir, "LRT_interaction_effect.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ============================================
# 5. Wald contrasts (using re-fit from full model)
# ============================================
cat("\n[5] Wald contrasts...\n")

# Re-fit full model
dds_full <- DESeq(dds, quiet=TRUE)

# 5a. Stage pairwise comparisons (T2 vs T1, T3 vs T1, T3 vs T2)
stage_comps <- list(
    c("T2", "T1"),
    c("T3", "T1"),
    c("T3", "T2")
)

for (comp in stage_comps) {
    label <- paste0(comp[1], "_vs_", comp[2])
    cat("  Stage:", label, "...")
    
    tryCatch({
        coef_name <- grep(paste0("stage_factor", comp[1], "_vs_", comp[2]), all_rn, value=TRUE)
        if (length(coef_name) > 0) {
            res_raw <- results(dds_full, name=coef_name[1])
            res_shrink <- lfcShrink(dds_full, coef=coef_name[1], type="apeglm", quiet=TRUE)
        } else {
            res_raw <- results(dds_full, contrast=c("stage_factor", comp[1], comp[2]))
            res_shrink <- lfcShrink(dds_full, contrast=c("stage_factor", comp[1], comp[2]),
                                    type="ashr", quiet=TRUE)
        }
        
        combined <- data.frame(
            gene_id = rownames(res_raw),
            baseMean = res_raw$baseMean,
            log2FC_raw = res_raw$log2FoldChange,
            lfcSE_raw = res_raw$lfcSE,
            stat = res_raw$stat,
            pvalue_raw = res_raw$pvalue,
            padj_raw = res_raw$padj,
            log2FC_shrink = res_shrink$log2FoldChange,
            lfcSE_shrink = res_shrink$lfcSE,
            padj_shrink = res_shrink$padj
        )
        combined <- combined[order(combined$padj_shrink), ]
        
        n_up <- sum(combined$padj_shrink < 0.05 & combined$log2FC_shrink > 1, na.rm=TRUE)
        n_down <- sum(combined$padj_shrink < 0.05 & combined$log2FC_shrink < -1, na.rm=TRUE)
        cat(" Up:", n_up, "Down:", n_down, "\n")
        
        write.table(combined,
                    file.path(results_dir, paste0("DEG_stage_", label, ".txt")),
                    sep="\t", row.names=FALSE, quote=FALSE)
    }, error=function(e) {
        cat(" ERROR:", e$message, "\n")
    })
}

# 5b. Node pairwise comparisons
node_comps <- list(
    c("5", "2"),
    c("10", "2"),
    c("10", "5")
)

for (comp in node_comps) {
    label <- paste0("Node", comp[1], "_vs_Node", comp[2])
    cat("  Node:", label, "...")
    
    tryCatch({
        coef_name <- grep(paste0("node_factor", comp[1], "_vs_", comp[2]), all_rn, value=TRUE)
        if (length(coef_name) > 0) {
            res_raw <- results(dds_full, name=coef_name[1])
            res_shrink <- lfcShrink(dds_full, coef=coef_name[1], type="apeglm", quiet=TRUE)
        } else {
            res_raw <- results(dds_full, contrast=c("node_factor", comp[1], comp[2]))
            res_shrink <- lfcShrink(dds_full, contrast=c("node_factor", comp[1], comp[2]),
                                    type="ashr", quiet=TRUE)
        }
        
        combined <- data.frame(
            gene_id = rownames(res_raw),
            baseMean = res_raw$baseMean,
            log2FC_raw = res_raw$log2FoldChange,
            lfcSE_raw = res_raw$lfcSE,
            stat = res_raw$stat,
            pvalue_raw = res_raw$pvalue,
            padj_raw = res_raw$padj,
            log2FC_shrink = res_shrink$log2FoldChange,
            lfcSE_shrink = res_shrink$lfcSE,
            padj_shrink = res_shrink$padj
        )
        combined <- combined[order(combined$padj_shrink), ]
        
        n_up <- sum(combined$padj_shrink < 0.05 & combined$log2FC_shrink > 1, na.rm=TRUE)
        n_down <- sum(combined$padj_shrink < 0.05 & combined$log2FC_shrink < -1, na.rm=TRUE)
        cat(" Up:", n_up, "Down:", n_down, "\n")
        
        write.table(combined,
                    file.path(results_dir, paste0("DEG_", label, ".txt")),
                    sep="\t", row.names=FALSE, quote=FALSE)
    }, error=function(e) {
        cat(" ERROR:", e$message, "\n")
    })
}

# ============================================
# 6. VST and PCA
# ============================================
cat("\n[6] VST and PCA...\n")

vsd <- vst(dds_full, blind=FALSE)

# Save VST matrix
vst_mat <- assay(vsd)
write.table(as.data.frame(vst_mat), file.path(results_dir, "VST_matrix.txt"),
            sep="\t", quote=FALSE)

# PCA colored by stage
pca_data <- plotPCA(vsd, intgroup=c("stage_factor", "node_factor"), returnData=TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

p1 <- ggplot(pca_data, aes(x=PC1, y=PC2, color=stage_factor, shape=node_factor)) +
    geom_point(size=3, alpha=0.8) +
    labs(title="GSE277812 - VST PCA (colored by Stage)",
         x=paste0("PC1 (", pct_var[1], "%)"),
         y=paste0("PC2 (", pct_var[2], "%)"),
         color="Stage", shape="Node") +
    theme_bw()

ggsave(file.path(results_dir, "PCA_by_stage.png"), p1, width=8, height=6, dpi=300)

# PCA colored by node
p2 <- ggplot(pca_data, aes(x=PC1, y=PC2, color=node_factor, shape=stage_factor)) +
    geom_point(size=3, alpha=0.8) +
    labs(title="GSE277812 - VST PCA (colored by Node)",
         x=paste0("PC1 (", pct_var[1], "%)"),
         y=paste0("PC2 (", pct_var[2], "%)"),
         color="Node", shape="Stage") +
    theme_bw()

ggsave(file.path(results_dir, "PCA_by_node.png"), p2, width=8, height=6, dpi=300)

# ============================================
# 7. Summary
# ============================================
cat("\n========================================\n")
cat("[7] Summary\n")
cat("========================================\n")

sink(file.path(results_dir, "SUMMARY_GSE277812.txt"))
cat("GSE277812 DESeq2 Analysis Summary\n")
cat("==================================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Design: ~ stage_factor * node_factor\n")
cat("Stages: T1, T2, T3\n")
cat("Nodes: 2, 5, 10\n")
cat("Samples: 27 (3 stages x 3 nodes x 3 reps)\n")
cat("Genes after filtering:", nrow(counts_f), "\n\n")
cat("LRT Results:\n")
cat(sprintf("  Overall model (stage * node): %d genes\n", n_overall))
cat(sprintf("  Stage effect: %d genes\n", n_stage))
cat(sprintf("  Node effect: %d genes\n", n_node))
cat(sprintf("  Stage x Node interaction: %d genes\n", n_int))
sink()

cat("  Overall LRT sig:", n_overall, "genes\n")
cat("  Stage LRT sig:", n_stage, "genes\n")
cat("  Node LRT sig:", n_node, "genes\n")
cat("  Interaction LRT sig:", n_int, "genes\n")

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(base_dir, "results_corrected/logs/sessionInfo_gse277812.txt"))

cat("\nDONE\n")
