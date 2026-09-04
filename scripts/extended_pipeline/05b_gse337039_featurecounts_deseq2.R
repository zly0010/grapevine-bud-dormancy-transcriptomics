#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 6)
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) {
  stop(paste(
    "Usage: 05b_gse337039_featurecounts_deseq2.R",
    "<work_dir> <manifest.tsv> <refseq_to_vitvi.tsv> <output_dir>",
    "[auto|unstranded|forward]"
  ))
}

work_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
mapping_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
output_dir <- args[[4]]
strand_request <- if (length(args) == 5) args[[5]] else "auto"
if (!strand_request %in% c("auto", "unstranded", "forward")) {
  stop("Strandedness must be auto, unstranded, or forward")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_path, check.names = FALSE)
required_meta <- c("run_accession", "cultivar", "condition", "time_index", "collection_date", "replicate")
if (!all(required_meta %in% names(manifest))) stop("Manifest is missing required columns")
if (nrow(manifest) != 60 || anyDuplicated(manifest$run_accession)) {
  stop("Manifest must contain 60 unique runs")
}

count_dir <- file.path(work_dir, "star_counts")
read_fc <- function(run, strand) {
  path <- file.path(count_dir, paste0(run, ".", strand, ".txt"))
  if (!file.exists(path)) stop("Missing featureCounts file: ", path)
  tab <- read.delim(path, comment.char = "#", check.names = FALSE)
  if (!all(c("Geneid", "Length") %in% names(tab))) stop("Invalid featureCounts file: ", path)
  count_col <- names(tab)[ncol(tab)]
  setNames(as.numeric(tab[[count_col]]), tab$Geneid)
}

read_assigned_rate <- function(run, strand) {
  path <- file.path(count_dir, paste0(run, ".", strand, ".txt.summary"))
  if (!file.exists(path)) return(NA_real_)
  tab <- read.delim(path, check.names = FALSE)
  value_col <- names(tab)[ncol(tab)]
  values <- as.numeric(tab[[value_col]])
  assigned <- values[tab$Status == "Assigned"]
  if (length(assigned) != 1 || sum(values) == 0) return(NA_real_)
  assigned / sum(values)
}

qc <- do.call(rbind, lapply(manifest$run_accession, function(run) {
  data.frame(
    run_accession = run,
    unstranded_assigned_rate = read_assigned_rate(run, "unstranded"),
    forward_assigned_rate = read_assigned_rate(run, "forward")
  )
}))

if (strand_request == "auto") {
  med_unstranded <- median(qc$unstranded_assigned_rate, na.rm = TRUE)
  med_forward <- median(qc$forward_assigned_rate, na.rm = TRUE)
  if (!is.finite(med_unstranded) || !is.finite(med_forward)) stop("Cannot infer strandedness")
  # In a forward-stranded library, forward counts approach the unstranded upper bound.
  # In an unstranded library, forward counts are typically about half of that bound.
  selected_strand <- if (med_forward / med_unstranded >= 0.80) "forward" else "unstranded"
} else {
  selected_strand <- strand_request
}
qc$selected_strand <- selected_strand
qc$selected_assigned_rate <- qc[[paste0(selected_strand, "_assigned_rate")]]

message("Merging ", selected_strand, " featureCounts files")
count_list <- lapply(manifest$run_accession, read_fc, strand = selected_strand)
all_genes <- Reduce(union, lapply(count_list, names))
counts <- vapply(count_list, function(x) {
  out <- numeric(length(all_genes))
  names(out) <- all_genes
  out[names(x)] <- x
  out
}, numeric(length(all_genes)))
colnames(counts) <- manifest$run_accession
storage.mode(counts) <- "integer"

mapping <- read.delim(mapping_path, check.names = FALSE)
mapping <- mapping[
  mapping$one_to_one %in% c(TRUE, "True", "TRUE", 1) &
    mapping$present_in_ncbi_gtf %in% c(TRUE, "True", "TRUE", 1),
  c("refseq_gene_id", "vitvi_id")
]
mapping <- mapping[!duplicated(mapping$refseq_gene_id) & !duplicated(mapping$vitvi_id), ]
refseq_to_vitvi <- setNames(mapping$vitvi_id, mapping$refseq_gene_id)

coldata <- manifest
rownames(coldata) <- coldata$run_accession
coldata <- coldata[colnames(counts), , drop = FALSE]
coldata$cultivar <- relevel(factor(coldata$cultivar), "Marquette")
coldata$condition <- relevel(factor(coldata$condition), "Natural")
coldata$time_factor <- relevel(factor(coldata$time_index), "0")

full_design <- ~ cultivar + condition + time_factor + cultivar:condition +
  cultivar:time_factor + condition:time_factor + cultivar:condition:time_factor
reduced_design <- ~ cultivar + condition + time_factor + cultivar:condition +
  cultivar:time_factor + condition:time_factor

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = full_design)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message("Retained ", sum(keep), " genes after count filter")

message("Fitting global three-way interaction LRT")
dds_lrt <- DESeq(dds, test = "LRT", reduced = reduced_design, parallel = FALSE)
three_way <- as.data.frame(results(dds_lrt))
three_way$refseq_gene_id <- rownames(three_way)
three_way$vitvi_id <- unname(refseq_to_vitvi[three_way$refseq_gene_id])
three_way <- three_way[, c("refseq_gene_id", "vitvi_id", setdiff(names(three_way), c("refseq_gene_id", "vitvi_id")))]

write_result <- function(result, label) {
  out <- as.data.frame(result)
  out$refseq_gene_id <- rownames(out)
  out$vitvi_id <- unname(refseq_to_vitvi[out$refseq_gene_id])
  out$contrast <- label
  out <- out[, c("refseq_gene_id", "vitvi_id", "contrast", setdiff(names(out), c("refseq_gene_id", "vitvi_id", "contrast")))]
  write.table(out, file.path(output_dir, paste0(label, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

contrast_index <- list()
cursor <- 1L
message("Fitting chilling-response time-course LRT within each cultivar")
for (cultivar_value in levels(coldata$cultivar)) {
  selected <- coldata$cultivar == cultivar_value
  dds_sub <- dds[, selected]
  dds_sub$condition <- droplevels(dds_sub$condition)
  dds_sub$time_factor <- droplevels(dds_sub$time_factor)
  design(dds_sub) <- ~ condition + time_factor + condition:time_factor
  fit <- DESeq(dds_sub, test = "LRT", reduced = ~ condition + time_factor, parallel = FALSE)
  label <- paste("lrt_condition_by_time", cultivar_value, sep = "__")
  write_result(results(fit), label)
  contrast_index[[cursor]] <- data.frame(contrast = label, family = "condition_by_time_LRT", cultivar = cultivar_value, condition = "Controlled_4C_vs_Natural", time = "all")
  cursor <- cursor + 1L
}

message("Fitting controlled-vs-natural contrasts at matched cultivar/time")
for (cultivar_value in levels(coldata$cultivar)) {
  for (time_value in levels(coldata$time_factor)) {
    selected <- coldata$cultivar == cultivar_value & coldata$time_factor == time_value
    dds_sub <- dds[, selected]
    dds_sub$condition <- droplevels(dds_sub$condition)
    design(dds_sub) <- ~ condition
    fit <- DESeq(dds_sub, test = "Wald", parallel = FALSE)
    label <- paste("condition", cultivar_value, paste0("T", time_value), "Controlled_vs_Natural", sep = "__")
    write_result(results(fit, contrast = c("condition", "Controlled_4C", "Natural")), label)
    contrast_index[[cursor]] <- data.frame(contrast = label, family = "condition_at_time", cultivar = cultivar_value, condition = "Controlled_4C_vs_Natural", time = time_value)
    cursor <- cursor + 1L
  }
}

message("Generating VST and PCA")
vsd <- vst(dds_lrt, blind = FALSE)
vst_refseq <- t(assay(vsd))
mapped <- unname(refseq_to_vitvi[colnames(vst_refseq)])
mapped_keep <- !is.na(mapped) & nzchar(mapped)
vst_vitvi <- vst_refseq[, mapped_keep, drop = FALSE]
colnames(vst_vitvi) <- mapped[mapped_keep]
pca <- prcomp(vst_refseq, center = TRUE, scale. = FALSE)
pca_table <- data.frame(
  run_accession = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  coldata[rownames(pca$x), c("cultivar", "condition", "time_index", "collection_date", "replicate")],
  row.names = NULL
)

mapped_counts <- counts[rownames(counts) %in% names(refseq_to_vitvi), , drop = FALSE]
rownames(mapped_counts) <- unname(refseq_to_vitvi[rownames(mapped_counts)])
write.table(three_way, file.path(output_dir, "01_three_way_interaction_LRT.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(qc, file.path(output_dir, "02_featurecounts_sample_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(pca_table, file.path(output_dir, "03_pca_scores.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(sample_id = rownames(vst_vitvi), vst_vitvi, check.names = FALSE), file.path(output_dir, "04_vst_samples_by_vitvi_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(refseq_gene_id = rownames(counts), counts, check.names = FALSE), file.path(output_dir, "05_refseq_gene_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(vitvi_id = rownames(mapped_counts), mapped_counts, check.names = FALSE), file.path(output_dir, "05b_vitvi_gene_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, contrast_index), file.path(output_dir, "06_contrast_index.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(selected_strand, file.path(output_dir, "SELECTED_STRANDEDNESS.txt"))
saveRDS(dds_lrt, file.path(output_dir, "gse337039_dds_lrt.rds"))
saveRDS(vsd, file.path(output_dir, "gse337039_vsd.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

if (any(qc$selected_assigned_rate < 0.50, na.rm = TRUE)) {
  warning(sum(qc$selected_assigned_rate < 0.50, na.rm = TRUE), " samples have featureCounts assignment below 50%")
}
message("GSE337039_FEATURECOUNTS_DESEQ2_OK")
