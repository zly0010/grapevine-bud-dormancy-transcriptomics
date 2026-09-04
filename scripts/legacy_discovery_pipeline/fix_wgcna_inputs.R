# WGCNA Input Audit and Fix
# Reads existing VST matrices + raw counts for GSE184114,
# applies make.unique where needed, filters qc_B for GSE124820,
# transposes all to samples-as-rows, validates strictly.
# Output: results_corrected/07_wgcna_fixed/ with corrected matrices + INPUT_AUDIT.txt
# Does NOT run full DESeq2 pipeline or WGCNA — only vst() for GSE184114.

suppressPackageStartupMessages({
    library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
PROJECT_ROOT <- normalizePath(
    if (length(args) >= 1L) args[[1]] else Sys.getenv("GRAPE_LEGACY_PROJECT", "."),
    winslash = "/", mustWork = TRUE
)
OUT_DIR <- file.path(PROJECT_ROOT, "results_corrected/07_wgcna_fixed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

audit_lines <- character()
audit <- function(...) {
    line <- paste0(...)
    cat(line, "\n")
    audit_lines <<- c(audit_lines, line)
}

audit("=== WGCNA INPUT AUDIT ===")
audit("Generated:", Sys.time())
audit("")

PASS <- TRUE
check <- function(label, condition, detail = "") {
    status <- if (condition) "PASS" else "FAIL"
    if (!condition) PASS <<- FALSE
    msg <- sprintf("[%s] %s %s", status, label, detail)
    audit(msg)
    return(invisible(condition))
}

# ============================================================
# 1. GSE124820 — filter to qc_B=True (172 samples), transpose
# ============================================================
audit("--- GSE124820 ---")

vst124 <- read.table(file.path(PROJECT_ROOT, "results_corrected/02_deseq2_gse124820/vst_matrix.txt"),
                      header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("  Raw VST:", nrow(vst124), "genes x", ncol(vst124), "samples\n")

design124 <- read.table(file.path(PROJECT_ROOT, "results_corrected/01_sample_design_and_qc/sample_design.tsv"),
                         header = TRUE, sep = "\t", check.names = FALSE, row.names = NULL)

qc_b_samples <- design124$sample_id[design124$qc_B == "True"]
cat("  qc_B=True samples:", length(qc_b_samples), "\n")

# Check design sample_ids are unique
design_dup_ids <- any(duplicated(design124$sample_id))
cat("  Design duplicate IDs:", design_dup_ids, "\n")

# Check VST sample_ids (column names) — they use the original GEO IDs
# VST colnames = design$sample_id? Check
vst_colnames <- colnames(vst124)
cat("  VST colnames (first 5):", paste(head(vst_colnames, 5), collapse = ", "), "\n")
cat("  Design sample_id (first 5):", paste(head(design124$sample_id, 5), collapse = ", "), "\n")

# Filter VST to qc_B=True samples
common_124 <- intersect(qc_b_samples, vst_colnames)
cat("  Common qc_B samples in VST:", length(common_124), "\n")

vst124_filtered <- vst124[, common_124, drop = FALSE]

# Transpose: samples as rows, genes as columns
vst124_t <- as.data.frame(t(as.matrix(vst124_filtered)))

# Check for missing values
n_nas_124 <- sum(is.na(vst124_t))
check("GSE124820: samples == 172", nrow(vst124_t) == 172,
      sprintf("(got %d)", nrow(vst124_t)))
check("GSE124820: duplicate sample IDs == 0", !any(duplicated(rownames(vst124_t))),
      sprintf("(dupes: %d)", sum(duplicated(rownames(vst124_t)))))
check("GSE124820: missing values == 0", n_nas_124 == 0,
      sprintf("(got %d)", n_nas_124))

# Check unmatched (samples in design but not in VST)
unmatched_124 <- setdiff(qc_b_samples, vst_colnames)
check("GSE124820: unmatched samples == 0", length(unmatched_124) == 0,
      sprintf("(unmatched: %d)", length(unmatched_124)))

audit(sprintf("  Output: %d samples x %d genes", nrow(vst124_t), ncol(vst124_t)))

# Save
write.table(vst124_t, file.path(OUT_DIR, "GSE124820_vst_fixed.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# ============================================================
# 2. GSE184114 — make.unique on raw counts, rebuild VST (74 samples)
# ============================================================
audit("")
audit("--- GSE184114 ---")

# Read raw counts
counts184_raw <- read.table(file.path(PROJECT_ROOT, "data/processed/GSE184114/counts_matrix.txt"),
                             header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("  Raw counts:", nrow(counts184_raw), "genes x", ncol(counts184_raw), "samples\n")
cat("  Raw colnames (first 8):", paste(head(colnames(counts184_raw), 8), collapse = ", "), "\n")

# Check duplicates before make.unique
raw_dup <- any(duplicated(colnames(counts184_raw)))
cat("  Raw duplicates:", raw_dup, "\n")
if (raw_dup) {
    dup_names <- names(which(table(colnames(counts184_raw)) > 1))
    cat("  Duplicated:", paste(dup_names, collapse = ", "), "\n")
}

# Apply make.unique
colnames(counts184_raw) <- make.unique(colnames(counts184_raw), sep = ".")
cat("  After make.unique:", ncol(counts184_raw), "samples, all unique:",
    !any(duplicated(colnames(counts184_raw))), "\n")

# Load design
design184 <- read.table(file.path(PROJECT_ROOT, "results_corrected/01_sample_design_and_qc/sample_design_gse184114.txt"),
                         header = TRUE, sep = "\t", check.names = FALSE, row.names = NULL)
cat("  Design:", nrow(design184), "samples\n")

# Verify all design sample_ids now match counts colnames
common184 <- intersect(design184$sample_id, colnames(counts184_raw))
cat("  Common samples:", length(common184), "\n")
unmatched184 <- setdiff(design184$sample_id, colnames(counts184_raw))
if (length(unmatched184) > 0) {
    cat("  UNMATCHED in counts:", paste(unmatched184, collapse = ", "), "\n")
}
extra184 <- setdiff(colnames(counts184_raw), design184$sample_id)
if (length(extra184) > 0) {
    cat("  EXTRA in counts:", paste(extra184, collapse = ", "), "\n")
}

# Subset to common samples and match design order
counts184 <- counts184_raw[, common184, drop = FALSE]
design184_sub <- design184[match(common184, design184$sample_id), ]

cat("  Matched counts:", nrow(counts184), "genes x", ncol(counts184), "samples\n")

# Filter low-count genes (same as original: >=10 counts in >=50% samples)
min_count <- 10
keep <- rowSums(counts184 >= min_count) >= ncol(counts184) * 0.5
counts184_f <- counts184[keep, ]
cat("  After filtering:", nrow(counts184_f), "genes\n")

# Build DESeqDataSet and compute VST
colData184 <- data.frame(design184_sub, row.names = design184_sub$sample_id)
dds184 <- DESeqDataSetFromMatrix(
    countData = round(counts184_f),
    colData = colData184,
    design = ~ 1
)
vsd184 <- vst(dds184, blind = FALSE)
vst184_mat <- assay(vsd184)
cat("  VST computed:", nrow(vst184_mat), "genes x", ncol(vst184_mat), "samples\n")

# Transpose: samples as rows, genes as columns
vst184_t <- as.data.frame(t(vst184_mat))

n_nas_184 <- sum(is.na(vst184_t))
check("GSE184114: samples == 74", nrow(vst184_t) == 74,
      sprintf("(got %d)", nrow(vst184_t)))
check("GSE184114: duplicate sample IDs == 0", !any(duplicated(rownames(vst184_t))),
      sprintf("(dupes: %d)", sum(duplicated(rownames(vst184_t)))))
check("GSE184114: missing values == 0", n_nas_184 == 0,
      sprintf("(got %d)", n_nas_184))
check("GSE184114: unmatched samples == 0", length(unmatched184) == 0,
      sprintf("(unmatched: %d)", length(unmatched184)))

audit(sprintf("  Output: %d samples x %d genes", nrow(vst184_t), ncol(vst184_t)))

# Save
write.table(vst184_t, file.path(OUT_DIR, "GSE184114_vst_fixed.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# ============================================================
# 3. GSE273240 — read, transpose, validate (90 samples)
# ============================================================
audit("")
audit("--- GSE273240 ---")

vst273 <- read.table(file.path(PROJECT_ROOT, "results_corrected/03_deseq2_gse273240/vst_matrix.txt"),
                      header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("  Raw VST:", nrow(vst273), "genes x", ncol(vst273), "samples\n")

# Verify sample IDs match design
design273 <- read.table(file.path(PROJECT_ROOT, "results_corrected/03_deseq2_gse273240/sample_design_gse273240.txt"),
                         header = TRUE, sep = "\t", check.names = FALSE)
cat("  Design:", nrow(design273), "samples\n")
common273 <- intersect(design273$sample_id, colnames(vst273))
cat("  Common:", length(common273), "\n")
unmatched273 <- setdiff(design273$sample_id, colnames(vst273))
if (length(unmatched273) > 0) cat("  UNMATCHED:", paste(unmatched273, collapse = ", "), "\n")

vst273_sub <- vst273[, common273, drop = FALSE]
vst273_t <- as.data.frame(t(as.matrix(vst273_sub)))

n_nas_273 <- sum(is.na(vst273_t))
check("GSE273240: samples == 90", nrow(vst273_t) == 90,
      sprintf("(got %d)", nrow(vst273_t)))
check("GSE273240: duplicate sample IDs == 0", !any(duplicated(rownames(vst273_t))),
      sprintf("(dupes: %d)", sum(duplicated(rownames(vst273_t)))))
check("GSE273240: missing values == 0", n_nas_273 == 0,
      sprintf("(got %d)", n_nas_273))
check("GSE273240: unmatched samples == 0", length(unmatched273) == 0,
      sprintf("(unmatched: %d)", length(unmatched273)))

audit(sprintf("  Output: %d samples x %d genes", nrow(vst273_t), ncol(vst273_t)))

write.table(vst273_t, file.path(OUT_DIR, "GSE273240_vst_fixed.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# ============================================================
# 4. GSE277812 — read, transpose, validate (27 samples)
# ============================================================
audit("")
audit("--- GSE277812 ---")

vst277 <- read.table(file.path(PROJECT_ROOT, "results_corrected/05_deseq2_gse277812/VST_matrix.txt"),
                      header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("  Raw VST:", nrow(vst277), "genes x", ncol(vst277), "samples\n")

# Verify sample IDs match design
design277 <- read.table(file.path(PROJECT_ROOT, "results_corrected/05_deseq2_gse277812/sample_design_gse277812.txt"),
                         header = TRUE, sep = "\t", check.names = FALSE)
cat("  Design:", nrow(design277), "samples\n")
common277 <- intersect(design277$sample_id, colnames(vst277))
cat("  Common:", length(common277), "\n")
unmatched277 <- setdiff(design277$sample_id, colnames(vst277))
if (length(unmatched277) > 0) cat("  UNMATCHED:", paste(unmatched277, collapse = ", "), "\n")

vst277_sub <- vst277[, common277, drop = FALSE]
vst277_t <- as.data.frame(t(as.matrix(vst277_sub)))

n_nas_277 <- sum(is.na(vst277_t))
check("GSE277812: samples == 27", nrow(vst277_t) == 27,
      sprintf("(got %d)", nrow(vst277_t)))
check("GSE277812: duplicate sample IDs == 0", !any(duplicated(rownames(vst277_t))),
      sprintf("(dupes: %d)", sum(duplicated(rownames(vst277_t)))))
check("GSE277812: missing values == 0", n_nas_277 == 0,
      sprintf("(got %d)", n_nas_277))
check("GSE277812: unmatched samples == 0", length(unmatched277) == 0,
      sprintf("(unmatched: %d)", length(unmatched277)))

audit(sprintf("  Output: %d samples x %d genes", nrow(vst277_t), ncol(vst277_t)))

write.table(vst277_t, file.path(OUT_DIR, "GSE277812_vst_fixed.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# ============================================================
# 5. Cross-dataset: check gene ID uniqueness per dataset
# ============================================================
audit("")
audit("--- Gene ID Uniqueness ---")
for (ds in c("GSE124820", "GSE184114", "GSE273240", "GSE277812")) {
    f <- file.path(OUT_DIR, paste0(ds, "_vst_fixed.txt"))
    m <- read.table(f, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, nrows = 2)
    n_gene_dup <- sum(duplicated(colnames(m)))
    check(sprintf("%s: gene IDs unique (dupes=%d)", ds, n_gene_dup), n_gene_dup == 0)
}

# ============================================================
# Summary
# ============================================================
audit("")
audit("=== FINAL VERDICT ===")
if (PASS) {
    audit("ALL CHECKS PASSED — WGCNA inputs are valid")
} else {
    audit("FAILED — fix issues above before proceeding")
}
audit("")
audit("Expected: GSE124820=172, GSE184114=74, GSE273240=90, GSE277812=27")
audit("Matrix format: rows=samples, columns=Vitvi gene IDs, all unique")

# Write audit file
writeLines(audit_lines, file.path(OUT_DIR, "INPUT_AUDIT.txt"))
cat("\nAudit saved to:", file.path(OUT_DIR, "INPUT_AUDIT.txt"), "\n")

if (!PASS) {
    stop("AUDIT FAILED — see INPUT_AUDIT.txt for details")
}
