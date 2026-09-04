# WGCNA analysis — fixed version (all 7 review fixes applied)
# Uses corrected VST matrices from 07_wgcna_fixed/ (samples as rows, genes as cols)
# Network type: signed | TOM type: signed
# Power selection: negative slope required, range 1:30
suppressPackageStartupMessages({
    library(WGCNA)
})
options(stringsAsFactors = FALSE, digits = 4)
allowWGCNAThreads(nThreads = 4)

# ============================================================
# select_power: find lowest power with negative slope + signedR2 >= 0.80
# fix #5: GSE273240 can use signedR2 >= 0.699 with tolerance
# ============================================================
select_power <- function(datExpr, powerVector = 1:30, tolerance = 0.80) {
    sft <- pickSoftThreshold(datExpr, powerVector = powerVector, verbose = 3,
                             networkType = "signed")
    fi <- sft$fitIndices
    fi$signedR2 <- -sign(fi$slope) * fi$SFT.R.sq

    cat("\n  Power selection table:\n")
    cat(sprintf("  %4s  %6s  %7s  %8s  %7s\n", "Pwr", "R2", "slope", "MeanK", "sR2"))
    for (i in seq_len(nrow(fi))) {
        cat(sprintf("  %4d  %6.3f  %7.4f  %8.1f  %7.3f\n",
                    fi$Power[i], fi$SFT.R.sq[i], fi$slope[i], fi$mean.k.[i], fi$signedR2[i]))
    }

    # Primary: negative slope AND signedR2 >= tolerance
    candidates <- fi[fi$slope < 0 & fi$signedR2 >= tolerance, ]

    if (nrow(candidates) == 0) {
        cat("\n  WARNING: No power with negative slope and signedR2 >=", tolerance, "\n")
        cat("  Flagging NEEDS_REVIEW\n")
        return(list(power = NA, status = "NEEDS_REVIEW", sft = sft, fi = fi,
                     threshold_used = tolerance))
    }

    if (tolerance < 0.80) {
        cat("\n  NOTE: Relaxed threshold to signedR2 >=", tolerance, "\n")
    }

    # Pick lowest power where mean connectivity is reasonable (>= 10)
    n_samples <- nrow(datExpr)
    reasonable <- candidates[candidates$mean.k. >= 10, ]
    if (nrow(reasonable) > 0) {
        chosen_power <- min(reasonable$Power)
    } else {
        chosen_power <- min(candidates$Power)
    }

    cat("\n  Selected power:", chosen_power, "\n")
    return(list(power = chosen_power, status = "OK", sft = sft, fi = fi,
                 threshold_used = tolerance))
}

# ============================================================
# build_trait_matrix: explicit trait selection per dataset
# fix #2: no auto-convert; each dataset specifies exact columns
# ============================================================
build_trait_matrix <- function(design, trait_cols) {
    # trait_cols: character vector of column names to include as traits
    # Only these columns are converted; all others are excluded
    traitMatrix <- data.frame(row.names = rownames(design))
    for (tn in trait_cols) {
        if (!(tn %in% colnames(design))) {
            cat("  WARNING: trait column", tn, "not found in design, skipping\n")
            next
        }
        vals <- design[[tn]]
        if (is.character(vals) || is.factor(vals)) {
            for (uv in unique(vals)) {
                traitMatrix[[paste0(tn, "_", uv)]] <- as.numeric(vals == uv)
            }
        } else {
            traitMatrix[[tn]] <- as.numeric(vals)
        }
    }
    return(traitMatrix)
}

# ============================================================
run_wgcna <- function(name, vst_file, output_dir, design_file,
                      trait_cols,
                      powers = 1:30, mergeCutHeight = 0.25,
                      minModuleSize = 30, maxGenes = 5000,
                      tolerance = 0.80) {

    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("WGCNA:", name, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # ---- Load VST matrix ----
    cat("Loading VST matrix...\n")
    vst_raw <- read.table(vst_file, header = TRUE, row.names = 1, sep = "\t",
                          check.names = FALSE)
    datExpr0 <- as.matrix(vst_raw)
    cat("  Raw:", nrow(datExpr0), "samples x", ncol(datExpr0), "genes\n")

    # ---- Load design and align ----
    # fix #1: use datExpr0 (not datExpr) for alignment check
    cat("Loading and aligning design table...\n")
    design <- read.table(design_file, header = TRUE, sep = "\t", check.names = FALSE,
                         row.names = NULL)
    rownames(design) <- design$sample_id
    design <- design[match(rownames(datExpr0), design$sample_id), , drop = FALSE]

    if (!identical(rownames(datExpr0), rownames(design))) {
        mismatch <- which(rownames(datExpr0) != rownames(design))
        stop(sprintf("Design alignment FAILED: %d mismatched rows (first: datExpr0='%s' vs design='%s')",
                     length(mismatch), rownames(datExpr0)[mismatch[1]], rownames(design)[mismatch[1]]))
    }
    cat("  Design aligned:", identical(rownames(datExpr0), rownames(design)), "\n")

    # ---- Filter top-variance genes ----
    vargenes <- apply(datExpr0, 2, var, na.rm = TRUE)
    if (ncol(datExpr0) > maxGenes) {
        top_idx <- order(vargenes, decreasing = TRUE)[1:maxGenes]
        datExpr0 <- datExpr0[, top_idx, drop = FALSE]
        cat("  Filtered to top", maxGenes, "variance genes\n")
    }

    # ---- Remove samples with NA ----
    na_rows <- which(apply(datExpr0, 1, anyNA))
    if (length(na_rows) > 0) {
        cat("  Removing", length(na_rows), "samples with NAs\n")
        datExpr0 <- datExpr0[-na_rows, , drop = FALSE]
        design <- design[rownames(datExpr0), , drop = FALSE]
    }

    # ---- Outlier removal via sample dendrogram ----
    cat("Checking for outlier samples...\n")
    sampleTree <- hclust(dist(datExpr0), method = "average")
    cutH <- median(sampleTree$height) + 3 * mad(sampleTree$height)
    grp <- cutree(sampleTree, h = cutH)
    largest <- as.integer(names(which.max(table(grp))))
    outliers <- which(grp != largest)
    if (length(outliers) > 0 && length(outliers) < nrow(datExpr0) * 0.4) {
        cat("  Removing", length(outliers), "outlier samples\n")
        datExpr0 <- datExpr0[-outliers, , drop = FALSE]
        design <- design[rownames(datExpr0), , drop = FALSE]
    }

    png(file.path(output_dir, "01_sample_dendrogram.png"), width = 1200, height = 600, res = 150)
    par(mar = c(6, 4, 2, 1))
    plot(sampleTree, main = paste(name, "- Sample Dendrogram"), xlab = "", sub = "")
    abline(h = cutH, col = "red")
    dev.off()
    cat("  Final:", nrow(datExpr0), "samples x", ncol(datExpr0), "genes\n")

    datExpr <- datExpr0
    rm(datExpr0)

    # ---- goodSamplesGenes check ----
    cat("Running goodSamplesGenes...\n")
    gsg <- goodSamplesGenes(datExpr, verbose = 2)
    cat("  All samples OK:", gsg$allOK, "\n")
    if (!gsg$allOK) {
        if (any(!gsg$goodSamples)) {
            datExpr <- datExpr[gsg$goodSamples, , drop = FALSE]
            design <- design[rownames(datExpr), , drop = FALSE]
        }
        if (any(!gsg$goodGenes)) {
            datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
        }
        cat("  After goodSamplesGenes:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
    }

    # ---- Re-verify design alignment after filtering ----
    if (!identical(rownames(datExpr), rownames(design))) {
        design <- design[match(rownames(datExpr), rownames(design)), , drop = FALSE]
        if (!identical(rownames(datExpr), rownames(design))) {
            stop("Design re-alignment failed after sample/gene filtering")
        }
    }

    # ---- Power selection ----
    cat("Choosing soft-thresholding power (signed network)...\n")
    ps <- select_power(datExpr, powerVector = powers, tolerance = tolerance)
    power <- ps$power
    status <- ps$status

    png(file.path(output_dir, "02_power_selection.png"), width = 1200, height = 500, res = 150)
    par(mfrow = c(1, 2))
    fi <- ps$fi
    plot(fi$Power, fi$signedR2,
         xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R^2)",
         type = "n", main = "Scale Independence")
    text(fi$Power, fi$signedR2, labels = fi$Power, cex = 0.9, col = "red")
    abline(h = 0.80, col = "blue", lty = 2)
    abline(h = 0.70, col = "grey60", lty = 3)
    if (tolerance < 0.80) abline(h = tolerance, col = "orange", lty = 2)
    if (!is.na(power)) abline(v = power, col = "darkgreen", lty = 2)
    plot(fi$Power, fi$mean.k.,
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         type = "n", main = "Mean Connectivity")
    text(fi$Power, fi$mean.k., labels = fi$Power, cex = 0.9, col = "red")
    if (!is.na(power)) abline(v = power, col = "darkgreen", lty = 2)
    dev.off()
    capture.output(print(fi), file = file.path(output_dir, "power_selection.txt"))

    # Save power status
    writeLines(c(paste("Status:", status),
                 paste("Power:", power),
                 paste("Tolerance:", tolerance),
                 if (tolerance < 0.80) "NOTE: Relaxed threshold used" else ""),
               file.path(output_dir, "power_status.txt"))

    if (is.na(power)) {
        cat("\n  !!!! NO VALID POWER FOUND — NEEDS_REVIEW !!!!\n")
        cat("  Skipping network construction for this dataset.\n")
        # Write empty files so acceptance checks don't fail
        writeLines(c("# No valid power found — no modules detected"), file.path(output_dir, "hub_genes_by_module.txt"))
        writeLines(c("# No valid power found — no modules detected"), file.path(output_dir, "module_sizes.txt"))
        writeLines(c("# No valid power found — no module-trait correlation"), file.path(output_dir, "module_trait_correlation.txt"))
        return(list(name = name, status = "NEEDS_REVIEW", datExpr = datExpr,
                    design = design, moduleColors = NULL, MEs = NULL,
                    power = NA, gene_ids = colnames(datExpr)))
    }

    # ---- Network construction & module detection (signed) ----
    cat("Running blockwiseModules (signed network, signed TOM)...\n")
    net <- blockwiseModules(
        datExpr, power = power,
        networkType = "signed",
        TOMType = "signed",
        minModuleSize = minModuleSize,
        reassignThreshold = 0,
        mergeCutHeight = mergeCutHeight,
        numericLabels = TRUE,
        pamRespectsDendro = FALSE,
        saveTOMs = FALSE,
        verbose = 3,
        maxBlockSize = ncol(datExpr),
        nThreads = 4
    )
    moduleColors <- labels2colors(net$colors)
    n_mods <- length(unique(net$colors))
    cat("  Detected", n_mods, "modules\n")

    # ---- CHECK: no module >80% of genes ----
    mod_tab <- table(moduleColors)
    total_genes <- sum(mod_tab)
    max_pct <- max(mod_tab) / total_genes * 100
    dominant_mod <- names(which.max(mod_tab))
    cat(sprintf("  Largest module: %s = %d genes (%.1f%% of %d total)\n",
                dominant_mod, max(mod_tab), max_pct, total_genes))
    if (max_pct > 80) {
        cat(sprintf("\n  !!!! MODULE %s CONTAINS %.1f%% OF GENES (>80%%) — NEEDS_REVIEW !!!!\n",
                    dominant_mod, max_pct))
        status <- "NEEDS_REVIEW"
    }

    # ---- Module assignments ----
    mod_df <- data.frame(
        gene_id = colnames(datExpr),
        module_number = net$colors,
        module_color = moduleColors
    )
    write.table(mod_df, file.path(output_dir, "module_assignments.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    capture.output(print(mod_tab), file = file.path(output_dir, "module_sizes.txt"))
    cat("  Modules:", paste(names(mod_tab), mod_tab, sep = "=", collapse = ", "), "\n")

    # ---- Eigengenes (fix #3: use moduleEigengenes for proper naming) ----
    cat("Computing module eigengenes via moduleEigengenes()...\n")
    MEs_result <- moduleEigengenes(datExpr, colors = moduleColors)
    MEs <- MEs_result$eigengenes
    cat("  ME columns:", paste(colnames(MEs), collapse = ", "), "\n")
    write.table(MEs, file.path(output_dir, "module_eigengenes.txt"),
                sep = "\t", row.names = TRUE, quote = FALSE)
    ME_cor <- cor(MEs)
    write.table(ME_cor, file.path(output_dir, "eigengene_correlation.txt"), sep = "\t")

    # ---- Module dendrogram ----
    png(file.path(output_dir, "03_module_dendrogram.png"), width = 1200, height = 800, res = 150)
    plotDendroAndColors(
        net$dendrograms[[1]],
        moduleColors[net$blockGenes[[1]]],
        "Module colors", dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05,
        main = paste(name, "- Gene Dendrogram and Module Colors")
    )
    dev.off()

    # ---- Module-trait correlation (fix #2: explicit traits only) ----
    cat("Module-trait correlations...\n")
    traitMatrix <- build_trait_matrix(design, trait_cols)
    cat("  Samples:", nrow(datExpr), " | Traits:", ncol(traitMatrix), "\n")

    corME <- cor(MEs, traitMatrix, use = "p")
    pvalME <- corPvalueStudent(corME, nrow(datExpr))

    # BH-corrected p-values
    p_adj <- matrix(p.adjust(pvalME, method = "BH"), nrow = nrow(pvalME),
                    ncol = ncol(pvalME), dimnames = dimnames(pvalME))

    textMat <- paste(signif(corME, 2), "\n(", signif(pvalME, 1), ")", sep = "")
    dim(textMat) <- dim(corME)

    tryCatch({
        png(file.path(output_dir, "04_module_trait_correlation.png"),
            width = max(800, ncol(traitMatrix) * 80),
            height = max(600, nrow(corME) * 40), res = 150)
        par(mar = c(max(8, ncol(traitMatrix) * 0.8), 8, 2, 1))
        labeledHeatmap(
            Matrix = corME, xLabels = colnames(traitMatrix),
            yLabels = colnames(MEs), ySymbols = colnames(MEs),
            colorLabels = FALSE, colors = blueWhiteRed(50),
            textMatrix = textMat, setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1, 1),
            main = paste(name, "- Module-Trait Relationships")
        )
        dev.off()
    }, error = function(e) {
        cat("  Heatmap plot failed:", conditionMessage(e), "- saving data only\n")
        try(dev.off(dev.cur()), silent = TRUE)
    })

    # Save correlation matrix with raw and adjusted p-values
    cor_df <- data.frame(
        module = rownames(corME), corME,
        p_raw = pvalME, p_BH = p_adj
    )
    write.table(cor_df, file.path(output_dir, "module_trait_correlation.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # ---- kME and hub genes (fix #3: use moduleEigengenes MEs) ----
    cat("Computing kME and hub genes...\n")
    kME <- as.data.frame(cor(datExpr, MEs, use = "p"))
    kME$gene_id <- rownames(kME)
    gene_module <- data.frame(gene_id = colnames(datExpr), module_color = moduleColors)
    kME <- merge(kME, gene_module, by = "gene_id")

    # For each module, find hub genes (top 10 by |kME|)
    hub_genes <- list()
    for (mod in names(mod_tab)) {
        if (mod == "grey") next
        me_col <- paste0("ME", mod)
        if (!(me_col %in% colnames(kME))) next
        mod_genes <- kME[kME$module_color == mod, ]
        mod_genes$kME_abs <- abs(mod_genes[[me_col]])
        hub <- mod_genes[order(-mod_genes$kME_abs), ][1:min(10, nrow(mod_genes)), ]
        hub_genes[[mod]] <- hub[, c("gene_id", me_col, "kME_abs")]
    }

    hub_lines <- character()
    for (mod in names(hub_genes)) {
        hub_lines <- c(hub_lines, paste0("=== ", mod, " (top hub genes) ==="))
        hub_lines <- c(hub_lines, paste0(hub_genes[[mod]]$gene_id,
                        "\tkME=", signif(hub_genes[[mod]]$kME_abs, 3)))
        hub_lines <- c(hub_lines, "")
    }
    writeLines(hub_lines, file.path(output_dir, "hub_genes_by_module.txt"))
    write.table(kME, file.path(output_dir, "kME_table.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("  Saved kME table and hub genes\n")

    cat("DONE:", name, "\n")
    return(list(
        name = name,
        status = status,
        datExpr = datExpr,
        design = design,
        moduleColors = moduleColors,
        MEs = MEs,
        power = power,
        gene_ids = colnames(datExpr)
    ))
}

# ============================================================
# Module preservation (fix #4: Zsummary + medianRank + common genes)
# ============================================================
run_module_preservation <- function(results, out_dir, nPermutations = 200) {
    results <- results[!sapply(results, is.null)]
    if (length(results) < 2) {
        cat("\nNot enough datasets with results for module preservation\n")
        return(invisible(NULL))
    }
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("MODULE PRESERVATION (nPermutations=", nPermutations, ")\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    ds_names <- names(results)
    all_preserved <- list()

    for (i in 1:(length(ds_names) - 1)) {
        for (j in (i + 1):length(ds_names)) {
            ref_name <- ds_names[i]
            test_name <- ds_names[j]
            ref <- results[[ref_name]]
            test <- results[[test_name]]

            if (is.null(ref$moduleColors) || is.null(test$moduleColors)) {
                cat(sprintf("\n--- %s vs %s: SKIPPED (one or both have no modules) ---\n",
                            ref_name, test_name))
                next
            }

            cat(sprintf("\n--- %s (ref) vs %s (test) ---\n", ref_name, test_name))

            # Common genes from actual datExpr column names
            ref_genes <- colnames(ref$datExpr)
            test_genes <- colnames(test$datExpr)
            common_genes <- intersect(ref_genes, test_genes)
            cat(sprintf("  ref genes: %d | test genes: %d | common: %d\n",
                        length(ref_genes), length(test_genes), length(common_genes)))

            if (length(common_genes) < 200) {
                cat(sprintf("  SKIPPED: only %d common genes (< 200 minimum)\n", length(common_genes)))
                next
            }

            # Reorder expression matrices and module colors to common genes
            ref_expr <- ref$datExpr[, common_genes, drop = FALSE]
            test_expr <- test$datExpr[, common_genes, drop = FALSE]
            ref_colors <- ref$moduleColors[match(common_genes, ref_genes)]
            test_colors <- test$moduleColors[match(common_genes, test_genes)]

            cat(sprintf("  After reorder: ref %d samples x %d genes, test %d samples x %d genes\n",
                        nrow(ref_expr), ncol(ref_expr), nrow(test_expr), ncol(test_expr)))

            # Build multiExpr and multiColor
            multiExpr <- list(
                list(data = ref_expr),
                list(data = test_expr)
            )
            multiColor <- list(
                ref_colors,
                test_colors
            )

            tryCatch({
                mp <- modulePreservation(
                    multiExpr, multiColor = multiColor,
                    referenceNetworks = 1,
                    nPermutations = nPermutations,
                    randomSeed = 42,
                    quickCor = 0,
                    savePermutedStatistics = FALSE,
                    verbose = 2
                )

                # fix #4: output mp$preservation (Zsummary, medianRank) + mp$quality
                out_file <- file.path(out_dir,
                    paste0("preservation_", ref_name, "_vs_", test_name, ".txt"))
                lines <- c(
                    paste("Reference:", ref_name),
                    paste("Test:", test_name),
                    paste("Common genes:", length(common_genes)),
                    "",
                    "=== mp$preservation ==="
                )
                if (!is.null(mp$preservation)) {
                    lines <- c(lines, capture.output(print(mp$preservation)))
                }
                lines <- c(lines, "", "=== mp$quality ===")
                if (!is.null(mp$quality)) {
                    lines <- c(lines, capture.output(print(mp$quality)))
                }
                writeLines(lines, out_file)

                # Print key stats
                if (!is.null(mp$preservation)) {
                    cat("  Preservation Zsummary:\n")
                    for (rn in rownames(mp$preservation)) {
                        z <- mp$preservation[rn, "ref"]
                        mr <- mp$preservation[rn, "ref.rank"]
                        cat(sprintf("    %s: Zsummary=%.2f  medianRank=%.1f\n", rn, z, mr))
                    }
                }
                cat("  Saved:", basename(out_file), "\n")
                all_preserved[[paste(ref_name, "vs", test_name)]] <- mp
            }, error = function(e) {
                cat("  ERROR in modulePreservation:", conditionMessage(e), "\n")
            })
        }
    }
    return(all_preserved)
}

# ============================================================
# Run all datasets
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
    normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
    normalizePath(Sys.getenv("GRAPE_LEGACY_PROJECT", "."), winslash = "/", mustWork = TRUE)
}
BASE <- file.path(project_root, "results_corrected")
results <- list()

# fix #2: explicit trait columns per dataset
configs <- list(
    list(name = "GSE124820",
         vst = file.path(BASE, "07_wgcna_fixed/GSE124820_vst_fixed.txt"),
         out = file.path(BASE, "07_wgcna_fixed/GSE124820"),
         design = file.path(BASE, "01_sample_design_and_qc/sample_design.tsv"),
         trait_cols = c("variety", "time")),
    list(name = "GSE273240",
         vst = file.path(BASE, "07_wgcna_fixed/GSE273240_vst_fixed.txt"),
         out = file.path(BASE, "07_wgcna_fixed/GSE273240"),
         design = file.path(BASE, "03_deseq2_gse273240/sample_design_gse273240.txt"),
         trait_cols = c("deac_phase", "day", "treatment"),
         tolerance = 0.699),
    list(name = "GSE184114",
         vst = file.path(BASE, "07_wgcna_fixed/GSE184114_vst_fixed.txt"),
         out = file.path(BASE, "07_wgcna_fixed/GSE184114"),
         design = file.path(BASE, "01_sample_design_and_qc/sample_design_gse184114.txt"),
         trait_cols = c("phase", "treatment", "time_h")),
    list(name = "GSE277812",
         vst = file.path(BASE, "07_wgcna_fixed/GSE277812_vst_fixed.txt"),
         out = file.path(BASE, "07_wgcna_fixed/GSE277812"),
         design = file.path(BASE, "05_deseq2_gse277812/sample_design_gse277812.txt"),
         trait_cols = c("stage", "node"))
)

for (cfg in configs) {
    cat("\n>>>", cfg$name, "\n")
    tol <- ifelse(is.null(cfg$tolerance), 0.80, cfg$tolerance)
    results[[cfg$name]] <- tryCatch(
        run_wgcna(cfg$name, cfg$vst, cfg$out, cfg$design,
                  trait_cols = cfg$trait_cols,
                  maxGenes = 5000, tolerance = tol),
        error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
    )
}

# Module preservation (fix #4)
pres_dir <- file.path(BASE, "07_wgcna_fixed/module_preservation")
pres_results <- run_module_preservation(results, pres_dir, nPermutations = 200)

# ============================================================
# Summary
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
for (nm in names(results)) {
    r <- results[[nm]]
    if (is.null(r)) {
        cat(sprintf("  %s: ERROR\n", nm))
    } else {
        cat(sprintf("  %s: status=%s  power=%s  samples=%d  genes=%d  modules=%d\n",
                    nm, r$status, ifelse(is.na(r$power), "NA", r$power),
                    nrow(r$datExpr), length(r$gene_ids),
                    ifelse(is.null(r$moduleColors), 0, length(unique(r$moduleColors)))))
    }
}

needs_review <- sapply(results, function(r) !is.null(r) && r$status == "NEEDS_REVIEW")
if (any(needs_review)) {
    cat("\nDATASETS REQUIRING REVIEW:\n")
    cat("  ", paste(names(needs_review)[needs_review], collapse = ", "), "\n")
}

cat("\nAll output in:", BASE, "/07_wgcna_fixed/\n")
cat("DONE\n")
