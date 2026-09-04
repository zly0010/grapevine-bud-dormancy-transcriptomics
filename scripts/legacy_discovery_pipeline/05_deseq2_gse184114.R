#!/usr/bin/env Rscript
# GSE184114 Corrected DESeq2 Analysis
# Acclimation and Deacclimation analyzed separately
# Time units: hours only (never "Day")
# 0h only has Control baseline (no ABA_0h)

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(apeglm)
    library(ashr)
})

cat("========================================\n")
cat("GSE184114 Corrected DESeq2 Analysis\n")
cat("========================================\n\n")

args <- commandArgs(trailingOnly = TRUE)
base_dir <- normalizePath(
    if (length(args) >= 1L) args[[1]] else Sys.getenv("GRAPE_LEGACY_PROJECT", "."),
    winslash = "/", mustWork = TRUE
)
results_dir <- file.path(base_dir, "results_corrected/04_deseq2_gse184114")
figures_dir <- file.path(base_dir, "results_corrected/figures")
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(figures_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================
# 1. Load data and build design
# ============================================
cat("[1] Loading data...\n")

counts_raw <- read.table(file.path(base_dir, "data/processed/GSE184114/counts_matrix.txt"),
                         header=TRUE, row.names=1, sep="\t", check.names=FALSE)
cat("  Count matrix:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")

# Parse sample names
sample_names <- colnames(counts_raw)
parse_sample <- function(name) {
    m <- regmatches(name, regexec("^(Accl|Deaccl)_(Acclimation|Deacclimation)_(Control|ABA)_(\\d+)_(\\d+)(?:\\.(\\d+))?", name))[[1]]
    if (length(m) == 0) {
        return(data.frame(sample_id=name, phase=NA, treatment=NA, time_h=NA, rep=NA, dup=NA, stringsAsFactors=FALSE))
    }
    data.frame(
        sample_id = name,
        phase = m[3],
        treatment = m[4],
        time_h = as.integer(m[5]),
        rep = as.integer(m[6]),
        dup = ifelse(is.na(m[7]), 0L, as.integer(m[7])),
        stringsAsFactors = FALSE
    )
}

design_full <- do.call(rbind, lapply(sample_names, parse_sample))
cat("  Total samples:", nrow(design_full), "\n")
cat("  Parsed OK:", sum(!is.na(design_full$phase)), "/", nrow(design_full), "\n")

for (phase in c("Acclimation", "Deacclimation")) {
    sub <- design_full[design_full$phase == phase, ]
    cat(sprintf("  %s: %d samples\n", phase, nrow(sub)))
    cat(sprintf("    Control_0h: %d samples (baseline)\n",
                sum(sub$treatment == "Control" & sub$time_h == 0)))
    times_post <- sort(unique(sub$time_h[sub$time_h > 0]))
    for (t in times_post) {
        n_ctrl <- sum(sub$treatment == "Control" & sub$time_h == t)
        n_aba <- sum(sub$treatment == "ABA" & sub$time_h == t)
        cat(sprintf("    %dh: Control=%d ABA=%d\n", t, n_ctrl, n_aba))
    }
}

# ============================================
# Helper function for one phase
# ============================================
run_phase <- function(phase_name, interaction_times) {
    cat("\n========================================\n")
    cat("PHASE:", phase_name, "\n")
    cat("========================================\n")
    
    idx <- design_full$phase == phase_name
    sub_design <- design_full[idx, ]
    
    outdir <- file.path(results_dir, phase_name)
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    
    # ---- MODEL 1: Interaction (post-treatment only) ----
    cat("\n--- MODEL 1: Interaction ---\n")
    cat("Post-treatment times:", paste(interaction_times, collapse=", "), "h\n")
    cat("Model formula: ~ time_factor * treatment\n")
    cat("Reference: 2h Control (Acclimation) / 6h Control (Deacclimation)\n")
    
    post_idx <- sub_design$time_h %in% interaction_times
    post_design <- sub_design[post_idx, ]
    post_counts <- data.matrix(counts_raw[, idx][, post_idx])
    
    post_design$time_factor <- factor(post_design$time_h, levels=as.character(interaction_times))
    post_design$treatment <- factor(post_design$treatment, levels=c("Control", "ABA"))
    
    cat("  Samples:", nrow(post_design), "\n")
    cat("  Groups:\n")
    for (t in interaction_times) {
        n_ctrl <- sum(post_design$time_h == t & post_design$treatment == "Control")
        n_aba <- sum(post_design$time_h == t & post_design$treatment == "ABA")
        cat(sprintf("    %dh: Control=%d ABA=%d\n", t, n_ctrl, n_aba))
    }
    
    # Filter
    min_s <- round(ncol(post_counts) * 0.1)
    keep <- rowSums(post_counts >= 10) >= min_s
    counts_f <- post_counts[keep, ]
    cat("  Genes after filtering:", nrow(counts_f), "/", nrow(post_counts), "\n")
    
    # DESeq2 interaction model
    dds_int <- DESeqDataSetFromMatrix(countData = counts_f,
                                       colData = post_design,
                                       design = ~ time_factor * treatment)
    
    # LRT for interaction
    cat("  LRT: full = ~ time_factor * treatment, reduced = ~ time_factor + treatment\n")
    dds_lrt <- DESeq(dds_int, test="LRT",
                      reduced=~ time_factor + treatment, quiet=TRUE)
    res_lrt <- results(dds_lrt)
    res_lrt <- res_lrt[order(res_lrt$padj), ]
    n_int <- sum(res_lrt$padj < 0.05, na.rm=TRUE)
    cat("  Interaction LRT sig:", n_int, "genes\n")
    
    lrt_df <- as.data.frame(res_lrt)
    lrt_df$gene_id <- rownames(lrt_df)
    lrt_df <- lrt_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.table(lrt_df, file.path(outdir, paste0("LRT_interaction_", phase_name, ".txt")),
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # LRT for time main effect
    cat("  LRT: full = ~ time_factor * treatment, reduced = ~ treatment\n")
    dds_lrt_time <- DESeq(dds_int, test="LRT", reduced=~ treatment, quiet=TRUE)
    res_lrt_time <- results(dds_lrt_time)
    res_lrt_time <- res_lrt_time[order(res_lrt_time$padj), ]
    n_time <- sum(res_lrt_time$padj < 0.05, na.rm=TRUE)
    cat("  Time LRT sig:", n_time, "genes\n")
    
    lrt_time_df <- as.data.frame(res_lrt_time)
    lrt_time_df$gene_id <- rownames(lrt_time_df)
    lrt_time_df <- lrt_time_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.table(lrt_time_df, file.path(outdir, paste0("LRT_time_", phase_name, ".txt")),
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Full model for Wald
    cat("  Fitting full model for Wald tests...\n")
    dds_full <- DESeq(dds_int, quiet=TRUE)
    all_rn <- resultsNames(dds_full)
    cat("  Coefficients:", paste(all_rn, collapse="\n    "), "\n")
    
    # Find treatment coefficient
    trt_coef <- grep("^treatment_ABA_vs_Control$", all_rn, value=TRUE)
    cat("  Treatment coef:", trt_coef, "\n")
    
    # Wald contrasts: ABA vs Control at each time point
    cat("  Wald contrasts (ABA vs Control per time point):\n")
    for (t in interaction_times) {
        time_label <- paste0(t, "h")
        cat("    ", time_label, "...")
        
        tryCatch({
            if (t == interaction_times[1]) {
                # Reference time: main treatment effect
                res_raw <- results(dds_full, name=trt_coef[1])
                res_shrink <- lfcShrink(dds_full, coef=trt_coef[1], type="apeglm", quiet=TRUE)
            } else {
                # Non-reference time: main + interaction
                int_coef <- grep(paste0("time_factor", t, "\\.treatmentABA"), all_rn, value=TRUE)
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
                        file.path(outdir, paste0("DEG_", phase_name, "_", time_label, "_ABA_vs_Control.txt")),
                        sep="\t", row.names=FALSE, quote=FALSE)
        }, error=function(e) {
            cat(" ERROR:", e$message, "\n")
        })
    }
    
    # PCA on interaction model
    cat("  Generating VST-PCA (interaction model)...\n")
    vsd <- vst(dds_full, blind=FALSE)
    pca_data <- plotPCA(vsd, intgroup=c("time_factor", "treatment"), returnData=TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))
    
    p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=time_factor, shape=treatment)) +
        geom_point(size=2.5, alpha=0.8) +
        labs(title=paste("GSE184114", phase_name, "- Interaction Model VST-PCA"),
             subtitle=paste("Times:", paste(interaction_times, collapse=", "), "h"),
             x=paste0("PC1 (", pct_var[1], "%)"),
             y=paste0("PC2 (", pct_var[2], "%)"),
             color="Time (h)", shape="Treatment") +
        scale_shape_manual(values=c("Control"=16, "ABA"=17)) +
        theme_bw()
    
    ggsave(file.path(outdir, paste0("PCA_", phase_name, "_interaction.png")), p,
           width=8, height=6, dpi=300)
    
    # ---- MODEL 2: Group comparison vs 0h baseline ----
    cat("\n--- MODEL 2: Group comparison vs 0h baseline ---\n")
    cat("Model formula: ~ group_factor\n")
    cat("Reference group: 0h_Control\n")
    
    # All time points including 0h
    all_times <- sort(unique(sub_design$time_h))
    cat("All times:", paste(all_times, collapse=", "), "h\n")
    
    sub_design$group_label <- paste0(sub_design$time_h, "h_", sub_design$treatment)
    sub_design$group_factor <- factor(sub_design$group_label)
    sub_design$group_factor <- relevel(sub_design$group_factor, ref="0h_Control")
    
    cat("Groups:\n")
    for (g in levels(sub_design$group_factor)) {
        n <- sum(sub_design$group_factor == g)
        cat(sprintf("  %s: %d\n", g, n))
    }
    
    all_counts <- data.matrix(counts_raw[, idx])
    min_s <- round(ncol(all_counts) * 0.1)
    keep <- rowSums(all_counts >= 10) >= min_s
    counts_f2 <- all_counts[keep, ]
    cat("  Genes after filtering:", nrow(counts_f2), "/", nrow(all_counts), "\n")
    
    dds_group <- DESeqDataSetFromMatrix(countData = counts_f2,
                                         colData = sub_design,
                                         design = ~ group_factor)
    
    cat("  Running LRT (group)...\n")
    dds_group <- DESeq(dds_group, quiet=TRUE)
    res_group_lrt <- results(dds_group)
    n_group <- sum(res_group_lrt$padj < 0.05, na.rm=TRUE)
    cat("  LRT sig:", n_group, "genes\n")
    
    # Wald contrasts: each group vs 0h_Control
    cat("  Wald contrasts (each group vs 0h_Control):\n")
    all_rn_g <- resultsNames(dds_group)
    
    group_levels <- levels(sub_design$group_factor)
    non_ref_groups <- group_levels[group_levels != "0h_Control"]
    
    # Sort by time then treatment for systematic output
    non_ref_sorted <- non_ref_groups[order(as.integer(gsub("h_.*", "", non_ref_groups)),
                                            gsub("\\d+h_", "", non_ref_groups))]
    
    for (g in non_ref_sorted) {
        # Exact coefficient name with underscore after group_factor
        coef_name <- paste0("group_factor_", g, "_vs_0h_Control")
        # R replaces - with . in factor names
        coef_name_dot <- gsub("-", ".", coef_name)
        
        if (coef_name %in% all_rn_g) {
            actual_coef <- coef_name
        } else if (coef_name_dot %in% all_rn_g) {
            actual_coef <- coef_name_dot
        } else {
            # Exact match only - no grep substring matching
            cat("    ", g, ": coef not found, skipped\n")
            next
        }
        
        cat("    ", g, "vs 0h_Control...")
        tryCatch({
            # Check if apeglm supports this coef, fallback to ashr
            tryCatch({
                res_shrink <- lfcShrink(dds_group, coef=actual_coef, type="apeglm", quiet=TRUE)
            }, error=function(e2) {
                assign("res_shrink",
                       lfcShrink(dds_group, name=actual_coef, type="ashr", quiet=TRUE),
                       envir=parent.frame(2))
            })
            
            res_raw <- results(dds_group, name=actual_coef)
            
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
            
            safe_name <- gsub("-", "_", g)
            write.table(combined,
                        file.path(outdir, paste0("DEG_", phase_name, "_", safe_name, "_vs_0h_Control.txt")),
                        sep="\t", row.names=FALSE, quote=FALSE)
        }, error=function(e) {
            cat(" ERROR:", e$message, "\n")
        })
    }
    
    # PCA on group model
    cat("  Generating VST-PCA (group model)...\n")
    vsd_g <- vst(dds_group, blind=FALSE)
    pca_data_g <- plotPCA(vsd_g, intgroup=c("group_factor", "treatment"), returnData=TRUE)
    pct_var_g <- round(100 * attr(pca_data_g, "percentVar"))
    pca_data_g$time_h <- as.integer(gsub("h_.*", "", pca_data_g$group_factor))
    
    p2 <- ggplot(pca_data_g, aes(x=PC1, y=PC2, color=factor(time_h), shape=treatment)) +
        geom_point(size=2.5, alpha=0.8) +
        labs(title=paste("GSE184114", phase_name, "- Group Model VST-PCA"),
             subtitle="Reference: 0h Control",
             x=paste0("PC1 (", pct_var_g[1], "%)"),
             y=paste0("PC2 (", pct_var_g[2], "%)"),
             color="Time (h)", shape="Treatment") +
        scale_shape_manual(values=c("Control"=16, "ABA"=17)) +
        theme_bw()
    
    ggsave(file.path(outdir, paste0("PCA_", phase_name, "_group_model.png")), p2,
           width=8, height=6, dpi=300)
    
    return(list(
        phase = phase_name,
        interaction_times = interaction_times,
        n_total = nrow(sub_design),
        n_post_treatment = nrow(post_design),
        n_interaction_model_genes = nrow(counts_f),
        n_group_model_genes = nrow(counts_f2),
        n_interaction_lrt_sig = n_int,
        n_time_lrt_sig = n_time,
        n_group_lrt_sig = n_group
    ))
}

# ============================================
# Run both phases
# ============================================
cat("\n[2] Running DESeq2 analysis...\n")

summary_accl <- run_phase("Acclimation", c(2, 4, 24, 48))
summary_deaccl <- run_phase("Deacclimation", c(6, 12, 24, 48, 72))

# ============================================
# 3. Summary
# ============================================
cat("\n========================================\n")
cat("[3] Summary Report\n")
cat("========================================\n")

sink(file.path(results_dir, "SUMMARY_GSE184114.txt"))
cat("GSE184114 Corrected Analysis Summary\n")
cat("=====================================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

for (s in list(summary_accl, summary_deaccl)) {
    cat(sprintf("--- %s ---\n", s$phase))
    cat(sprintf("  Total samples: %d\n", s$n_total))
    cat(sprintf("  Post-treatment samples (interaction model): %d\n", s$n_post_treatment))
    cat(sprintf("  Interaction model genes: %d\n", s$n_interaction_model_genes))
    cat(sprintf("  Group model genes: %d\n", s$n_group_model_genes))
    cat(sprintf("  Interaction LRT sig: %d\n", s$n_interaction_lrt_sig))
    cat(sprintf("  Time main effect LRT sig: %d\n", s$n_time_lrt_sig))
    cat(sprintf("  Group LRT sig: %d\n", s$n_group_lrt_sig))
    cat(sprintf("  Interaction model: ~ time_factor * treatment\n"))
    cat(sprintf("  Group model: ~ group_factor (ref: 0h_Control)\n"))
    cat("\n")
}
sink()

for (s in list(summary_accl, summary_deaccl)) {
    cat(sprintf("  %s: %d total, %d post-tx, %d int-sig, %d time-sig, %d group-sig\n",
                s$phase, s$n_total, s$n_post_treatment,
                s$n_interaction_lrt_sig, s$n_time_lrt_sig, s$n_group_lrt_sig))
}

writeLines(capture.output(sessionInfo()),
           file.path(base_dir, "results_corrected/logs/sessionInfo_gse184114.txt"))

cat("\nDONE\n")
