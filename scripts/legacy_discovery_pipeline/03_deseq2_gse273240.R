#!/usr/bin/env Rscript
# GSE273240 Corrected DESeq2 Analysis (FIXED v2)
# Per deac phase: design = ~ day_factor * treatment
# LRT for time x treatment interaction
# Wald contrasts via list: main treatment + interaction per time point

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(apeglm)
    library(ashr)
})

cat("========================================\n")
cat("GSE273240 Corrected DESeq2 Analysis (v2)\n")
cat("========================================\n\n")

args <- commandArgs(trailingOnly = TRUE)
base_dir <- normalizePath(
    if (length(args) >= 1L) args[[1]] else Sys.getenv("GRAPE_LEGACY_PROJECT", "."),
    winslash = "/", mustWork = TRUE
)
results_dir <- file.path(base_dir, "results_corrected/03_deseq2_gse273240")
figures_dir <- file.path(base_dir, "results_corrected/figures")
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(figures_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================
# 1. Load data
# ============================================
cat("[1] Loading data...\n")

counts <- read.table(file.path(base_dir, "data/processed/GSE273240/counts_matrix.txt"),
                     header=TRUE, row.names=1, sep="\t", check.names=FALSE)
cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Parse sample info
sample_names <- colnames(counts)
parse_sample <- function(name) {
    parts <- strsplit(name, "_")[[1]]
    data.frame(
        sample_id = name,
        deac_phase = parts[1],
        day = as.numeric(gsub("d", "", parts[2])),
        treatment = parts[3],
        replicate = as.numeric(gsub("rep", "", parts[4])),
        stringsAsFactors = FALSE
    )
}

design <- do.call(rbind, lapply(sample_names, parse_sample))
design$day_factor <- paste0("Day", design$day)
design$treatment <- factor(design$treatment, levels=c("Control", "tetralone-ABA"))

cat("  Deac phases:", paste(unique(design$deac_phase), collapse=", "), "\n")
cat("  Days:", paste(sort(unique(design$day)), collapse=", "), "\n")
cat("  Treatments:", paste(unique(design$treatment), collapse=", "), "\n")
cat("  Samples per phase:\n")
for (p in unique(design$deac_phase)) {
    sub <- design[design$deac_phase == p, ]
    cat("    ", p, ":", nrow(sub), "samples\n")
}

# Save design
write.table(design, file.path(results_dir, "sample_design_gse273240.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ============================================
# 2. Analysis per deac phase
# ============================================
cat("\n[2] Running DESeq2 per deac phase...\n")

all_summaries <- list()

for (phase in unique(design$deac_phase)) {
    cat("\n========================================\n")
    cat("PHASE:", phase, "\n")
    cat("========================================\n")
    
    idx <- design$deac_phase == phase
    sub_counts <- counts[, idx]
    sub_design <- design[idx, ]
    
    cat("  Samples:", nrow(sub_design), "\n")
    cat("  Days:", paste(sort(unique(sub_design$day)), collapse=", "), "\n")
    
    # Filter
    counts_mat <- data.matrix(sub_counts)
    min_samples <- round(ncol(counts_mat) * 0.1)
    keep <- rowSums(counts_mat >= 10) >= min_samples
    counts_f <- counts_mat[keep, ]
    cat("  Genes after filtering:", nrow(counts_f), "/", nrow(counts_mat), "\n")
    
    # Factor setup
    days_sorted <- sort(unique(sub_design$day))
    sub_design$day_factor <- factor(sub_design$day_factor,
                                     levels=paste0("Day", days_sorted))
    sub_design$day_factor <- relevel(sub_design$day_factor, ref="Day0")
    
    outdir <- file.path(results_dir, phase)
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    
    # ---- Model 1: Interaction model ----
    cat("\n  --- MODEL 1: Interaction (~ day_factor * treatment) ---\n")
    
    dds_int <- DESeqDataSetFromMatrix(
        countData = counts_f,
        colData = sub_design,
        design = ~ day_factor * treatment
    )
    
    # LRT for interaction: does treatment effect vary by time?
    cat("  LRT: full = ~ day_factor * treatment, reduced = ~ day_factor + treatment\n")
    dds_int <- DESeq(dds_int, test="LRT",
                      reduced=~ day_factor + treatment, quiet=TRUE)
    res_lrt_int <- results(dds_int)
    res_lrt_int <- res_lrt_int[order(res_lrt_int$padj), ]
    n_int <- sum(res_lrt_int$padj < 0.05, na.rm=TRUE)
    cat("  Interaction LRT sig:", n_int, "genes\n")
    
    lrt_int_df <- as.data.frame(res_lrt_int)
    lrt_int_df$gene_id <- rownames(lrt_int_df)
    lrt_int_df <- lrt_int_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.table(lrt_int_df, file.path(outdir, paste0("LRT_interaction_", phase, ".txt")),
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # LRT for time main effect
    cat("  LRT: full = ~ day_factor * treatment, reduced = ~ treatment\n")
    dds_lrt_time <- DESeq(dds_int, test="LRT",
                           reduced=~ treatment, quiet=TRUE)
    res_lrt_time <- results(dds_lrt_time)
    res_lrt_time <- res_lrt_time[order(res_lrt_time$padj), ]
    n_time <- sum(res_lrt_time$padj < 0.05, na.rm=TRUE)
    cat("  Time LRT sig:", n_time, "genes\n")
    
    lrt_time_df <- as.data.frame(res_lrt_time)
    lrt_time_df$gene_id <- rownames(lrt_time_df)
    lrt_time_df <- lrt_time_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.table(lrt_time_df, file.path(outdir, paste0("LRT_time_", phase, ".txt")),
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Re-fit full model for Wald tests
    cat("  Fitting full model for Wald tests...\n")
    dds_full <- DESeq(dds_int, quiet=TRUE)
    all_rn <- resultsNames(dds_full)
    cat("  Coefficients:", paste(all_rn, collapse="\n    "), "\n")
    
    # Find coefficient names
    # R replaces - with . in factor names, so "tetralone-ABA" becomes "tetralone.ABA"
    trt_coef <- grep("^treatment", all_rn, value=TRUE)
    cat("  Treatment coef:", trt_coef, "\n")
    
    # Wald contrasts per time point
    cat("\n  --- Wald contrasts: tetralone-ABA vs Control per time point ---\n")
    
    for (d in days_sorted) {
        day_label <- paste0("Day", d)
        cat("    ", day_label, "...")
        
        tryCatch({
            if (d == 0) {
                # Reference day: just main treatment effect
                res_raw <- results(dds_full, name=trt_coef[1])
                res_shrink <- lfcShrink(dds_full, coef=trt_coef[1], type="apeglm", quiet=TRUE)
            } else {
                # Non-reference day: treatment main + interaction
                int_coef <- grep(paste0("day_factorDay", d, "\\.treatment"), all_rn, value=TRUE)
                if (length(int_coef) > 0) {
                    cat(" [coef:", trt_coef[1], "+", int_coef[1], "]")
                    res_raw <- results(dds_full, contrast=list(c(trt_coef[1], int_coef[1])))
                    res_shrink <- lfcShrink(dds_full, contrast=list(c(trt_coef[1], int_coef[1])),
                                            type="ashr", quiet=TRUE)
                } else {
                    cat(" [fallback]")
                    res_raw <- results(dds_full, name=trt_coef[1])
                    res_shrink <- lfcShrink(dds_full, coef=trt_coef[1], type="apeglm", quiet=TRUE)
                }
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
                        file.path(outdir, paste0("DEG_", phase, "_", day_label, "_ABA_vs_Control.txt")),
                        sep="\t", row.names=FALSE, quote=FALSE)
        }, error=function(e) {
            cat(" ERROR:", e$message, "\n")
        })
    }
    
    # Sensitivity: per-time-point ~treatment
    cat("\n  --- SENSITIVITY: Per-time-point ~treatment ---\n")
    for (d in days_sorted) {
        day_label <- paste0("Day", d)
        day_idx <- sub_design$day == d
        day_counts <- counts_f[, day_idx]
        day_info <- sub_design[day_idx, ]
        
        if (ncol(day_counts) < 4) next
        
        tryCatch({
            dds_sens <- DESeqDataSetFromMatrix(
                countData = day_counts,
                colData = day_info,
                design = ~ treatment
            )
            dds_sens <- DESeq(dds_sens, quiet=TRUE)
            res_sens <- results(dds_sens, contrast=c("treatment", "tetralone-ABA", "Control"))
            res_sens <- res_sens[order(res_sens$padj), ]
            
            n_up <- sum(res_sens$padj < 0.05 & res_sens$log2FoldChange > 1, na.rm=TRUE)
            n_down <- sum(res_sens$padj < 0.05 & res_sens$log2FoldChange < -1, na.rm=TRUE)
            cat("    ", day_label, ": Up", n_up, "Down", n_down, "\n")
            
            sens_df <- as.data.frame(res_sens)
            sens_df$gene_id <- rownames(sens_df)
            sens_df <- sens_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
            write.table(sens_df,
                        file.path(outdir, paste0("SENS_", phase, "_", day_label, "_ABA_vs_Control.txt")),
                        sep="\t", row.names=FALSE, quote=FALSE)
        }, error=function(e) {
            # Skip
        })
    }
    
    # PCA
    cat("  Generating VST-PCA...\n")
    vsd <- vst(dds_full, blind=FALSE)
    pca_data <- plotPCA(vsd, intgroup=c("day_factor", "treatment"), returnData=TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))
    
    p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=day_factor, shape=treatment)) +
        geom_point(size=2.5, alpha=0.8) +
        labs(title=paste("GSE273240", phase, "- Interaction Model VST-PCA"),
             x=paste0("PC1 (", pct_var[1], "%)"),
             y=paste0("PC2 (", pct_var[2], "%)"),
             color="Day", shape="Treatment") +
        scale_shape_manual(values=c("Control"=16, "tetralone-ABA"=17)) +
        theme_bw()
    
    ggsave(file.path(outdir, paste0("PCA_", phase, "_interaction.png")), p, width=8, height=6, dpi=300)
    
    # Save summary
    all_summaries[[phase]] <- list(
        n_samples = nrow(sub_design),
        n_genes_filtered = nrow(counts_f),
        n_interaction_lrt = n_int,
        n_time_lrt = n_time
    )
}

# ============================================
# 3. Summary
# ============================================
cat("\n========================================\n")
cat("[3] Summary\n")
cat("========================================\n")

sink(file.path(results_dir, "SUMMARY_GSE273240.txt"))
cat("GSE273240 Corrected Analysis Summary (v2)\n")
cat("==========================================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Model: ~ day_factor * treatment\n")
cat("LRT: interaction (reduced = ~ day_factor + treatment)\n")
cat("Wald: main treatment + interaction per time point\n")
cat("Sensitivity: per-time-point ~treatment\n\n")

for (phase in names(all_summaries)) {
    s <- all_summaries[[phase]]
    cat(sprintf("--- %s ---\n", phase))
    cat(sprintf("  Samples: %d\n", s$n_samples))
    cat(sprintf("  Genes after filtering: %d\n", s$n_genes_filtered))
    cat(sprintf("  Interaction LRT sig: %d\n", s$n_interaction_lrt))
    cat(sprintf("  Time main effect LRT sig: %d\n", s$n_time_lrt))
    cat("\n")
}
sink()

for (phase in names(all_summaries)) {
    s <- all_summaries[[phase]]
    cat(sprintf("  %s: %d samples, %d genes, %d interaction-sig, %d time-sig\n",
                phase, s$n_samples, s$n_genes_filtered, s$n_interaction_lrt, s$n_time_lrt))
}

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(base_dir, "results_corrected/logs/sessionInfo_gse273240.txt"))

cat("\nDONE\n")
