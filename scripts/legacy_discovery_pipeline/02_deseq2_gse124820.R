#!/usr/bin/env Rscript
# GSE124820 Corrected DESeq2 Analysis
# 4 varieties x 3 QC schemes x LRT + lfcShrink

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(apeglm)
    library(ashr)
})

cat("========================================\n")
cat("GSE124820 Corrected DESeq2 Analysis\n")
cat("========================================\n\n")

args <- commandArgs(trailingOnly = TRUE)
base_dir <- normalizePath(
    if (length(args) >= 1L) args[[1]] else Sys.getenv("GRAPE_LEGACY_PROJECT", "."),
    winslash = "/", mustWork = TRUE
)
results_base <- file.path(base_dir, "results_corrected/02_deseq2_gse124820")
figures_dir <- file.path(base_dir, "results_corrected/figures")
dir.create(results_base, showWarnings=FALSE, recursive=TRUE)
dir.create(figures_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================
# 1. Load data
# ============================================
cat("[1] Loading data...\n")

counts_file <- file.path(base_dir, "data/processed/GSE124820/counts_matrix.txt")
counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

design_file <- file.path(base_dir, "results_corrected/01_sample_design_and_qc/sample_design.tsv")
design <- read.table(design_file, header=TRUE, sep="\t", check.names=FALSE)

cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("  Design:", nrow(design), "samples\n")

# Verify column match
if (!all(design$sample_id %in% colnames(counts))) {
    stop("Sample IDs in design do not match count matrix columns!")
}

# Reorder counts to match design
counts <- counts[, design$sample_id]
cat("  Verified: columns match design\n\n")

# Fix Python True/False -> R TRUE/FALSE
for (col in c("qc_A", "qc_B", "qc_C")) {
    design[[col]] <- as.logical(toupper(design[[col]]))
}

# ============================================
# 2. Define QC schemes
# ============================================
cat("[2] QC schemes:\n")
qc_schemes <- list(
    A = list(name="A_all", label="All samples (189)", samples=design$sample_id),
    B = list(name="B_no_fail", label="No fail (172)", samples=design$sample_id[design$qc_B]),
    C = list(name="C_pass_only", label="Pass only (147)", samples=design$sample_id[design$qc_C])
)
for (nm in names(qc_schemes)) {
    cat("  Scheme", nm, ":", qc_schemes[[nm]]$label, "\n")
}
cat("\n")

# ============================================
# 3. Define varieties
# ============================================
varieties <- c("Vamu", "Vvcs", "Vvri", "Vrip")
variety_names <- c(
    Vamu = "Vitis amurensis PI588635",
    Vvcs = "Vitis vinifera Cabernet Sauvignon",
    Vvri = "Vitis vinifera Riesling",
    Vrip = "Vitis riparia PI588275"
)

# ============================================
# 4. Helper functions
# ============================================

run_deseq2_analysis <- function(counts, design_sub, variety, qc_name, output_dir) {
    cat("\n  ---", variety, "| QC", qc_name, "---\n")
    
    dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
    
    # Filter low-expression genes
    counts_mat <- data.matrix(counts)
    min_samples <- round(ncol(counts_mat) * 0.1)
    keep <- rowSums(counts_mat >= 10) >= min_samples
    counts_f <- counts_mat[keep, ]
    cat("    Genes after filtering:", nrow(counts_f), "/", nrow(counts_mat), "\n")
    cat("    Count matrix mode:", mode(counts_f), "class:", class(counts_f), "\n")
    
    # Ensure time_factor is factor with Day0 as reference
    design_sub$time_factor <- as.character(design_sub$time_factor)
    all_days <- sort(unique(design_sub$time_factor))
    all_time_levels <- c("Day0", all_days[all_days != "Day0"])
    design_sub$time_factor <- factor(design_sub$time_factor, levels=all_time_levels)
    if (!"Day0" %in% levels(design_sub$time_factor)) {
        stop(paste("Day0 not found in time_factor levels for", variety, qc_name))
    }
    design_sub$time_factor <- relevel(design_sub$time_factor, ref = "Day0")
    design_sub$variety <- factor(design_sub$variety)
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
        countData = counts_f,
        colData = design_sub,
        design = ~ time_factor
    )
    
    # --- LRT test ---
    cat("    Running LRT...\n")
    dds_lrt <- DESeq(dds, test="LRT", reduced=~1, quiet=TRUE)
    
    res_lrt <- results(dds_lrt)
    res_lrt <- res_lrt[order(res_lrt$padj), ]
    
    n_sig <- sum(res_lrt$padj < 0.05, na.rm=TRUE)
    cat("    LRT significant (padj<0.05):", n_sig, "genes\n")
    
    # Save LRT results
    lrt_df <- as.data.frame(res_lrt)
    lrt_df$gene_id <- rownames(lrt_df)
    lrt_df <- lrt_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.table(lrt_df, file.path(output_dir, paste0("LRT_", variety, "_", qc_name, ".txt")),
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # --- Wald tests: each time vs Day0 with lfcShrink ---
    cat("    Running Wald tests with lfcShrink...\n")
    
    # Re-run with standard DESeq (not LRT) for Wald tests
    dds_wald <- DESeq(dds, quiet=TRUE)
    
    days <- sort(unique(as.character(design_sub$time_factor)))
    days <- days[days != "Day0"]
    
    for (day in days) {
        cat("      ", day, "vs Day0...")
        
        res <- results(dds_wald, contrast=c("time_factor", day, "Day0"))
        
        # lfcShrink with apeglm
        coef_name <- paste0("time_factor_", day, "_vs_Day0")
        tryCatch({
            res_shrink <- lfcShrink(dds_wald, coef=coef_name, type="apeglm", quiet=TRUE)
        }, error=function(e) {
            cat(" apeglm failed, trying ashr...")
            res_shrink <<- lfcShrink(dds_wald, contrast=c("time_factor", day, "Day0"), type="ashr", quiet=TRUE)
        })
        
        # Combine results
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
        
        write.table(combined, 
                    file.path(output_dir, paste0("DEG_", variety, "_", day, "_vs_Day0_", qc_name, ".txt")),
                    sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    # --- PCA ---
    cat("    Generating PCA...\n")
    vsd <- vst(dds_wald, blind=FALSE)
    pca_data <- plotPCA(vsd, intgroup="time_factor", returnData=TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))
    
    p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=time_factor)) +
        geom_point(size=2.5, alpha=0.8) +
        labs(title=paste(variety_names[variety], "-", qc_name),
             subtitle=paste("VST, Top 500 variable genes"),
             x=paste0("PC1 (", pct_var[1], "%)"),
             y=paste0("PC2 (", pct_var[2], "%)"),
             color="Time") +
        theme_bw() +
        scale_color_brewer(palette="Set3")
    
    ggsave(file.path(output_dir, paste0("PCA_", variety, "_", qc_name, ".png")), p, width=8, height=6, dpi=300)
    
    # --- Sample correlation heatmap ---
    cor_mat <- cor(assay(vsd), method="pearson")
    anno_col <- data.frame(Time=design_sub$time_factor, row.names=design_sub$sample_id)
    
    png(file.path(output_dir, paste0("Correlation_", variety, "_", qc_name, ".png")),
        width=max(800, ncol(cor_mat)*8), height=max(800, ncol(cor_mat)*8), res=150)
    pheatmap(cor_mat,
             annotation_col=anno_col,
             color=colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             main=paste(variety_names[variety], "-", qc_name, "Sample Correlation"),
             show_rownames=FALSE, show_colnames=FALSE,
             fontsize=6)
    dev.off()
    
    return(list(dds=dds_wald, vsd=vsd, res_lrt=res_lrt, n_sig_lrt=n_sig))
}

# ============================================
# 5. Run analysis for each variety x QC scheme
# ============================================
cat("[3] Running DESeq2 analysis...\n")

all_results <- list()

for (variety in varieties) {
    cat("\n========================================\n")
    cat("VARIETY:", variety, "-", variety_names[variety], "\n")
    cat("========================================\n")
    
    idx <- design$variety == variety
    variety_counts <- counts[, idx]
    variety_design <- design[idx, ]
    
    for (qc_name in names(qc_schemes)) {
        scheme <- qc_schemes[[qc_name]]
        sample_ids <- intersect(scheme$samples, variety_design$sample_id)
        
        sub_counts <- variety_counts[, sample_ids]
        sub_design <- variety_design[variety_design$sample_id %in% sample_ids, ]
        
        out_dir <- file.path(results_base, paste0(variety, "_", qc_schemes[[qc_name]]$name))
        
        result <- run_deseq2_analysis(sub_counts, sub_design, variety, qc_schemes[[qc_name]]$name, out_dir)
        
        key <- paste0(variety, "_", qc_name)
        all_results[[key]] <- result
    }
}

# ============================================
# 6. Cross-variety consensus (QC scheme B)
# ============================================
cat("\n========================================\n")
cat("[4] Cross-variety consensus (QC scheme B)\n")
cat("========================================\n")

consensus_dir <- file.path(results_base, "cross_variety_consensus")
dir.create(consensus_dir, showWarnings=FALSE, recursive=TRUE)

# Read LRT results for scheme B
lrt_genes <- list()
for (variety in varieties) {
    lrt_file <- file.path(results_base, paste0(variety, "_B_no_fail"), paste0("LRT_", variety, "_B_no_fail.txt"))
    if (file.exists(lrt_file)) {
        lrt <- read.table(lrt_file, header=TRUE, sep="\t", check.names=FALSE)
        lrt_genes[[variety]] <- setNames(lrt$padj < 0.05, lrt$gene_id)
        cat("  ", variety, ": LRT significant =", sum(lrt_genes[[variety]], na.rm=TRUE), "\n")
    }
}

# Find genes significant in at least 3/4 varieties
all_genes <- unique(unlist(lapply(lrt_genes, names)))
consensus_table <- data.frame(gene_id=all_genes)

for (variety in varieties) {
    if (variety %in% names(lrt_genes)) {
        consensus_table[[variety]] <- lrt_genes[[variety]][consensus_table$gene_id]
    } else {
        consensus_table[[variety]] <- NA
    }
}
consensus_table[is.na(consensus_table)] <- FALSE

# Count support
consensus_table$n_support <- rowSums(consensus_table[, varieties])

cat("\n  LRT gene overlap:\n")
cat("    4/4 varieties:", sum(consensus_table$n_support == 4), "\n")
cat("    >=3/4 varieties:", sum(consensus_table$n_support >= 3), "\n")
cat("    >=2/4 varieties:", sum(consensus_table$n_support >= 2), "\n")
cat("    Any:", sum(consensus_table$n_support >= 1), "\n")

write.table(consensus_table, file.path(consensus_dir, "LRT_consensus_table_B.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ============================================
# 7. Time-response gene overlap at Day10
# ============================================
cat("\n[5] Day10 vs Day0 DEG overlap (scheme B, shrunken):\n")

day10_genes <- list()
for (variety in varieties) {
    deg_dir <- file.path(results_base, paste0(variety, "_B_no_fail"))
    deg_file <- file.path(deg_dir, paste0("DEG_", variety, "_Day10_vs_Day0_B_no_fail.txt"))
    if (file.exists(deg_file)) {
        deg <- read.table(deg_file, header=TRUE, sep="\t", check.names=FALSE)
        sig <- deg[deg$padj_shrink < 0.05 & abs(deg$log2FC_shrink) > 1, ]
        day10_genes[[variety]] <- setNames(sig$log2FC_shrink, sig$gene_id)
        cat("  ", variety, ": DEGs =", length(day10_genes[[variety]]), "\n")
    }
}

# Overlap
if (length(day10_genes) == 4) {
    common_all <- Reduce(intersect, lapply(day10_genes, names))
    cat("\n  Day10 DEGs common to all 4 varieties:", length(common_all), "\n")
    
    # Direction consistency
    if (length(common_all) > 0) {
        direction_df <- data.frame(gene_id=common_all)
        for (v in varieties) {
            direction_df[[v]] <- day10_genes[[v]][common_all]
        }
        
        all_up <- apply(direction_df[, varieties], 1, function(x) all(x > 0))
        all_down <- apply(direction_df[, varieties], 1, function(x) all(x < 0))
        mixed <- !(all_up | all_down)
        
        cat("    Consistent up:", sum(all_up), "\n")
        cat("    Consistent down:", sum(all_down), "\n")
        cat("    Mixed direction:", sum(mixed), "\n")
        
        direction_df$consensus <- ifelse(all_up, "up", ifelse(all_down, "down", "mixed"))
        direction_df$mean_abs_lfc <- rowMeans(abs(direction_df[, varieties]))
        direction_df <- direction_df[order(-direction_df$mean_abs_lfc), ]
        
        write.table(direction_df, file.path(consensus_dir, "Day10_consensus_DEGs_B.txt"),
                    sep="\t", row.names=FALSE, quote=FALSE)
        
        cat("\n  Top 20 consensus Day10 genes:\n")
        print(head(direction_df[, c("gene_id", varieties, "consensus", "mean_abs_lfc")], 20))
    }
}

# ============================================
# 8. QC scheme comparison
# ============================================
cat("\n========================================\n")
cat("[6] QC scheme comparison\n")
cat("========================================\n")

qc_summary <- data.frame()

for (variety in varieties) {
    for (qc_name in names(qc_schemes)) {
        key <- paste0(variety, "_", qc_name)
        if (key %in% names(all_results)) {
            r <- all_results[[key]]
            n_lrt <- r$n_sig_lrt
            n_samples <- ncol(r$vsd)
            
            # Count DEGs at Day10
            deg_dir <- file.path(results_base, paste0(variety, "_", qc_schemes[[qc_name]]$name))
            deg_file <- file.path(deg_dir, paste0("DEG_", variety, "_Day10_vs_Day0_", qc_schemes[[qc_name]]$name, ".txt"))
            n_deg10 <- 0
            if (file.exists(deg_file)) {
                deg <- read.table(deg_file, header=TRUE, sep="\t", check.names=FALSE)
                n_deg10 <- sum(deg$padj_shrink < 0.05 & abs(deg$log2FC_shrink) > 1, na.rm=TRUE)
            }
            
            qc_summary <- rbind(qc_summary, data.frame(
                variety=variety,
                qc_scheme=qc_name,
                n_samples=n_samples,
                n_lrt_sig=n_lrt,
                n_deg_day10=n_deg10
            ))
        }
    }
}

cat("\n  QC sensitivity comparison:\n")
print(qc_summary)

write.table(qc_summary, file.path(results_base, "qc_sensitivity_comparison.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ============================================
# 9. Session info
# ============================================
cat("\n[7] Saving session info...\n")
writeLines(capture.output(sessionInfo()),
           file.path(base_dir, "results_corrected/sessionInfo.txt"))

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
