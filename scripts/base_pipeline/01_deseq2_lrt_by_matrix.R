suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript R/01_deseq2_lrt_by_matrix.R <count_matrix.tsv.gz> <design.tsv> <output_dir>")
}

count_path <- args[[1]]
design_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readr::read_tsv(count_path, show_col_types = FALSE)
gene_id <- counts[[1]]
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- gene_id
storage.mode(count_matrix) <- "integer"

design <- readr::read_tsv(design_path, show_col_types = FALSE)
if (!"qc_status" %in% names(design)) {
  design$qc_status <- "pass"
}
design <- design %>%
  filter(is.na(qc_status) | qc_status != "fail")

keep_samples <- intersect(colnames(count_matrix), design$sample_accession)
count_matrix <- count_matrix[, keep_samples, drop = FALSE]
design <- design %>%
  filter(sample_accession %in% keep_samples) %>%
  distinct(sample_accession, .keep_all = TRUE)
design <- design[match(colnames(count_matrix), design$sample_accession), ]
stopifnot(all(colnames(count_matrix) == design$sample_accession))

design <- design %>%
  mutate(
    genotype = factor(genotype),
    treatment = factor(treatment),
    experiment = factor(experiment),
    time_factor = factor(time, levels = unique(time)),
    replicate = factor(replicate)
  )

candidate_terms <- c("genotype", "treatment", "experiment")
available_terms <- candidate_terms[vapply(design[candidate_terms], function(x) nlevels(factor(x)) > 1, logical(1))]
full_terms <- c(available_terms, "time_factor")
reduced_terms <- available_terms
full_formula <- as.formula(paste("~", paste(full_terms, collapse = " + ")))
reduced_formula <- if (length(reduced_terms) > 0) {
  as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
} else {
  ~ 1
}

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = as.data.frame(design),
  design = full_formula
)

keep_genes <- rowSums(counts(dds) >= 10) >= max(3, floor(ncol(dds) * 0.10))
dds <- dds[keep_genes, ]

dds_lrt <- DESeq(dds, test = "LRT", reduced = reduced_formula)
lrt <- results(dds_lrt)
lrt_tbl <- as.data.frame(lrt) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)
readr::write_tsv(lrt_tbl, file.path(output_dir, "deseq2_lrt_time.tsv"))

vsd <- vst(dds, blind = FALSE)
vst_tbl <- as.data.frame(assay(vsd)) %>%
  rownames_to_column("gene_id")
readr::write_tsv(vst_tbl, file.path(output_dir, "vst_expression.tsv.gz"))

pca <- plotPCA(vsd, intgroup = c("genotype", "treatment", "time_factor"), returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = genotype, shape = treatment, label = time_factor)) +
  geom_point(size = 2.5) +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw()
ggsave(file.path(output_dir, "pca_vst.png"), p, width = 7, height = 5, dpi = 300)

saveRDS(dds_lrt, file.path(output_dir, "dds_lrt.rds"))
