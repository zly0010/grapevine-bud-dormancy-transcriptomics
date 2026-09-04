#!/usr/bin/env Rscript
# GSE124820 DESeq2 - Single variety, all 3 QC schemes
# Usage: Rscript 02a_deseq2_one_variety.R <Vamu|Vvcs|Vvri|Vrip> [project_root]

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Usage: Rscript 02a_deseq2_one_variety.R <variety> [project_root]")
variety <- args[1]
cat("Analyzing variety:", variety, "\n")

suppressPackageStartupMessages({
    library(DESeq2)
    library(apeglm)
    library(ashr)
})

base_dir <- normalizePath(
    if (length(args) >= 2L) args[[2]] else Sys.getenv("GRAPE_LEGACY_PROJECT", "."),
    winslash = "/", mustWork = TRUE
)
results_base <- file.path(base_dir, "results_corrected/02_deseq2_gse124820")

# Load data
counts <- read.table(file.path(base_dir, "data/processed/GSE124820/counts_matrix.txt"),
                     header=TRUE, row.names=1, sep="\t", check.names=FALSE)
design <- read.table(file.path(base_dir, "results_corrected/01_sample_design_and_qc/sample_design.tsv"),
                     header=TRUE, sep="\t", check.names=FALSE)

# Fix Python True/False
for (col in c("qc_A", "qc_B", "qc_C")) {
    design[[col]] <- as.logical(toupper(design[[col]]))
}

counts <- counts[, design$sample_id]

# Subset to this variety
idx <- design$variety == variety
variety_counts <- counts[, idx]
variety_design <- design[idx, ]
cat("Samples:", nrow(variety_design), "\n")

# QC schemes
qc_schemes <- list(
    A = list(name="A_all", label="All", samples=design$sample_id[design$qc_A & idx]),
    B = list(name="B_no_fail", label="No fail", samples=design$sample_id[design$qc_B & idx]),
    C = list(name="C_pass_only", label="Pass only", samples=design$sample_id[design$qc_C & idx])
)

for (qc_name in names(qc_schemes)) {
    sample_ids <- qc_schemes[[qc_name]]$samples
    cat("\n=== QC", qc_name, ":", qc_schemes[[qc_name]]$label, "- n=", length(sample_ids), "===\n")
    
    sub_counts <- variety_counts[, sample_ids]
    sub_design <- variety_design[variety_design$sample_id %in% sample_ids, ]
    
    # Filter
    counts_mat <- data.matrix(sub_counts)
    min_samples <- round(ncol(counts_mat) * 0.1)
    keep <- rowSums(counts_mat >= 10) >= min_samples
    counts_f <- counts_mat[keep, ]
    cat("  Genes:", nrow(counts_f), "/", nrow(counts_mat), "\n")
    
    # Factor
    sub_design$time_factor <- as.character(sub_design$time_factor)
    all_days <- sort(unique(sub_design$time_factor))
    all_time_levels <- c("Day0", all_days[all_days != "Day0"])
    sub_design$time_factor <- factor(sub_design$time_factor, levels=all_time_levels)
    sub_design$time_factor <- relevel(sub_design$time_factor, ref = "Day0")
    
    dds <- DESeqDataSetFromMatrix(countData=counts_f, colData=sub_design, design=~time_factor)
    
    # LRT
    cat("  LRT...\n")
    dds_lrt <- DESeq(dds, test="LRT", reduced=~1, quiet=TRUE)
    res_lrt <- results(dds_lrt)
    res_lrt <- res_lrt[order(res_lrt$padj), ]
    cat("  LRT sig:", sum(res_lrt$padj < 0.05, na.rm=TRUE), "\n")
    
    lrt_df <- as.data.frame(res_lrt)
    lrt_df$gene_id <- rownames(lrt_df)
    lrt_df <- lrt_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    outdir <- file.path(results_base, paste0(variety, "_", qc_schemes[[qc_name]]$name))
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    write.table(lrt_df, file.path(outdir, paste0("LRT_", variety, "_", qc_name, ".txt")),
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Wald + lfcShrink
    dds_wald <- DESeq(dds, quiet=TRUE)
    
    days <- all_days[all_days != "Day0"]
    for (day in days) {
        cat("   ", day, "vs Day0...")
        res <- results(dds_wald, contrast=c("time_factor", day, "Day0"))
        
        coef_name <- paste0("time_factor_", day, "_vs_Day0")
        tryCatch({
            res_shrink <- lfcShrink(dds_wald, coef=coef_name, type="apeglm", quiet=TRUE)
        }, error=function(e) {
            res_shrink <<- lfcShrink(dds_wald, contrast=c("time_factor", day, "Day0"), type="ashr", quiet=TRUE)
        })
        
        combined <- data.frame(
            gene_id = rownames(res),
            baseMean = res$baseMean,
            log2FC_raw = res$log2FoldChange,
            lfcSE_raw = res$lfcSE,
            stat = res$stat,
            pvalue_raw = res$pvalue,
            padj_raw = res$padj,
            log2FC_shrink = res_shrink$log2FoldChange,
            lfcSE_shrink = res_shrink$lfcSE,
            padj_shrink = res_shrink$padj
        )
        
        n_up <- sum(combined$padj_shrink < 0.05 & combined$log2FC_shrink > 1, na.rm=TRUE)
        n_down <- sum(combined$padj_shrink < 0.05 & combined$log2FC_shrink < -1, na.rm=TRUE)
        cat(" Up:", n_up, "Down:", n_down, "\n")
        
        write.table(combined, file.path(outdir, paste0("DEG_", variety, "_", day, "_vs_Day0_", qc_name, ".txt")),
                    sep="\t", row.names=FALSE, quote=FALSE)
    }
}

cat("\nDONE:", variety, "\n")
