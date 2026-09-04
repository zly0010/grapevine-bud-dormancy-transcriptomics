#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 08b_go_enrichment.R <project_root> <output_dir>")
}

project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

mapping_file <- file.path(project_root, "results_corrected/06_gene_annotation/vitvi_to_ensembl_mapping.tsv")
go_file <- file.path(project_root, "results_corrected/06_gene_annotation/grape_go_annotation.tsv")
module_file <- file.path(project_root, "results_corrected/07_wgcna_fixed/GSE124820/module_assignments.txt")
audit_file <- file.path(project_root, "results_corrected/08_wgcna_postprocessing/01_go_mapping_audit/module_go_coverage_audit.tsv")
required_files <- c(mapping_file, go_file, module_file, audit_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Required files are missing: ", paste(missing_files, collapse = "; "))
}

message("Checking prerequisite mapping audit")
audit <- read.delim(audit_file, check.names = FALSE, quote = "", fileEncoding = "UTF-8")
if (!"audit_pass" %in% names(audit) || nrow(audit) != 3L || !all(as.logical(audit$audit_pass))) {
  stop("GO mapping audit prerequisite is absent or did not pass")
}

message("Reading formal modules and authoritative annotations")
mapping <- read.delim(mapping_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
go <- read.delim(go_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
modules <- read.delim(module_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")

required_mapping_columns <- c("vitvi_id", "ensembl_id")
required_go_columns <- c("ensembl_gene_id", "go_id", "go_term", "go_type")
required_module_columns <- c("gene_id", "module_color")
if (!all(required_mapping_columns %in% names(mapping))) stop("Invalid mapping columns")
if (!all(required_go_columns %in% names(go))) stop("Invalid GO annotation columns")
if (!all(required_module_columns %in% names(modules))) stop("Invalid module assignment columns")

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}
for (column in required_mapping_columns) mapping[[column]] <- clean_id(mapping[[column]])
for (column in required_go_columns) go[[column]] <- clean_id(go[[column]])
for (column in required_module_columns) modules[[column]] <- clean_id(modules[[column]])

if (nrow(modules) != 5000L || anyNA(modules$gene_id) || anyDuplicated(modules$gene_id)) {
  stop("Formal GSE124820 WGCNA background must contain exactly 5,000 unique non-missing genes")
}

target_modules <- c("blue", "turquoise", "brown")
ontology_map <- c(P = "BP", F = "MF", C = "CC")
min_overlap <- 3L
fdr_threshold <- 0.05

module_sizes <- table(modules$module_color)
expected_sizes <- c(blue = 1045L, turquoise = 2333L, brown = 680L)
if (!all(as.integer(module_sizes[names(expected_sizes)]) == expected_sizes)) {
  stop("Target module sizes differ from the audited formal network")
}

mapping_valid <- unique(mapping[!is.na(mapping$vitvi_id) & !is.na(mapping$ensembl_id), required_mapping_columns, drop = FALSE])
go_valid <- unique(go[!is.na(go$ensembl_gene_id) & !is.na(go$go_id) & !is.na(go$go_type), required_go_columns, drop = FALSE])
if (!all(go_valid$go_type %in% names(ontology_map))) {
  stop("GO annotation contains unsupported ontology codes")
}
go_valid$ontology <- unname(ontology_map[go_valid$go_type])

type_count <- aggregate(go_type ~ go_id, data = unique(go_valid[, c("go_id", "go_type")]), FUN = length)
if (any(type_count$go_type != 1L)) {
  stop("At least one GO ID is assigned to multiple ontology types")
}

term_records <- unique(go_valid[, c("go_id", "go_term", "ontology")])
term_name_count <- aggregate(go_term ~ go_id, data = term_records, FUN = function(x) length(unique(x[!is.na(x)])))
conflicting_term_ids <- term_name_count$go_id[term_name_count$go_term > 1L]
if (length(conflicting_term_ids) > 0L) {
  stop("Conflicting term names detected for GO IDs: ", paste(head(conflicting_term_ids, 20L), collapse = ", "))
}
term_info <- term_records[!duplicated(term_records$go_id), c("go_id", "go_term", "ontology")]

gene_map <- merge(modules[, c("gene_id", "module_color")], mapping_valid,
                  by.x = "gene_id", by.y = "vitvi_id", all = FALSE, sort = FALSE)
gene_go <- merge(gene_map, go_valid,
                 by.x = "ensembl_id", by.y = "ensembl_gene_id", all = FALSE, sort = FALSE)
gene_go <- unique(gene_go[, c("gene_id", "module_color", "ensembl_id", "go_id", "go_term", "ontology")])
gene_term <- unique(gene_go[, c("gene_id", "module_color", "go_id", "go_term", "ontology")])

if (length(unique(gene_term$gene_id)) != 1037L) {
  stop("GO-annotated background gene count differs from passed audit: observed ", length(unique(gene_term$gene_id)))
}

background_size <- nrow(modules)
all_results <- list()
result_index <- 0L

message("Running one-sided hypergeometric over-representation tests")
for (module_name in target_modules) {
  module_size <- sum(modules$module_color == module_name)
  for (ontology_name in c("BP", "MF", "CC")) {
    ontology_terms <- term_info[term_info$ontology == ontology_name, , drop = FALSE]
    ontology_gene_term <- gene_term[gene_term$ontology == ontology_name, , drop = FALSE]
    background_counts <- table(ontology_gene_term$go_id)
    ontology_terms <- ontology_terms[ontology_terms$go_id %in% names(background_counts), , drop = FALSE]
    ontology_terms <- ontology_terms[order(ontology_terms$go_id), , drop = FALSE]

    module_gene_term <- ontology_gene_term[ontology_gene_term$module_color == module_name, , drop = FALSE]
    module_counts <- table(module_gene_term$go_id)
    term_ids <- ontology_terms$go_id
    K <- as.integer(background_counts[term_ids])
    k <- as.integer(module_counts[term_ids])
    k[is.na(k)] <- 0L

    p_value <- phyper(k - 1L, K, background_size - K, module_size, lower.tail = FALSE)
    fdr <- p.adjust(p_value, method = "BH")
    fold_enrichment <- (k / module_size) / (K / background_size)
    fold_enrichment[k == 0L] <- 0

    overlap_split <- split(module_gene_term$gene_id, module_gene_term$go_id)
    overlap_genes <- vapply(term_ids, function(term_id) {
      values <- sort(unique(overlap_split[[term_id]]))
      if (length(values) == 0L) "" else paste(values, collapse = ";")
    }, character(1))

    result <- data.frame(
      module_color = module_name,
      ontology = ontology_name,
      go_id = term_ids,
      go_term = ontology_terms$go_term,
      background_size = background_size,
      module_size = module_size,
      background_term_genes = K,
      module_term_genes = k,
      background_ratio = K / background_size,
      module_ratio = k / module_size,
      fold_enrichment = fold_enrichment,
      p_value = p_value,
      fdr_bh = fdr,
      min_overlap_pass = k >= min_overlap,
      significant = (k >= min_overlap) & (fdr < fdr_threshold),
      overlap_gene_ids = overlap_genes,
      stringsAsFactors = FALSE
    )
    result <- result[order(result$fdr_bh, result$p_value, -result$fold_enrichment, result$go_id), ]

    write.table(result, file.path(output_dir, paste0(module_name, "_", ontology_name, "_all_terms.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    write.table(result[result$significant, , drop = FALSE],
                file.path(output_dir, paste0(module_name, "_", ontology_name, "_significant.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "")

    result_index <- result_index + 1L
    all_results[[result_index]] <- result
  }
}

combined <- do.call(rbind, all_results)
combined <- combined[order(match(combined$module_color, target_modules), match(combined$ontology, c("BP", "MF", "CC")), combined$fdr_bh, combined$p_value), ]
significant <- combined[combined$significant, , drop = FALSE]

write.table(combined, file.path(output_dir, "go_enrichment_all_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(significant, file.path(output_dir, "go_enrichment_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(gene_term, file.path(output_dir, "audited_background_gene_to_go.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

summary_grid <- expand.grid(module_color = target_modules, ontology = c("BP", "MF", "CC"), stringsAsFactors = FALSE)
summary_grid <- summary_grid[order(match(summary_grid$module_color, target_modules), match(summary_grid$ontology, c("BP", "MF", "CC"))), ]
summary_grid$terms_tested <- mapply(function(m, o) sum(combined$module_color == m & combined$ontology == o), summary_grid$module_color, summary_grid$ontology)
summary_grid$terms_overlap_ge_3 <- mapply(function(m, o) sum(combined$module_color == m & combined$ontology == o & combined$min_overlap_pass), summary_grid$module_color, summary_grid$ontology)
summary_grid$significant_terms <- mapply(function(m, o) sum(combined$module_color == m & combined$ontology == o & combined$significant), summary_grid$module_color, summary_grid$ontology)
write.table(summary_grid, file.path(output_dir, "enrichment_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameters <- data.frame(
  parameter = c("mapping_chain", "target_modules", "background", "background_size", "test", "alternative", "multiple_testing", "correction_scope", "fdr_threshold", "minimum_overlap_genes", "significance_rule"),
  value = c(
    "gene_id -> vitvi_id -> ensembl_id -> ensembl_gene_id -> GO term",
    paste(target_modules, collapse = ","),
    "All genes in formal GSE124820 module_assignments.txt", as.character(background_size),
    "Hypergeometric over-representation test", "greater", "Benjamini-Hochberg",
    "Separately within each module x ontology family", as.character(fdr_threshold), as.character(min_overlap),
    "fdr_bh < 0.05 and module_term_genes >= 3"
  ),
  stringsAsFactors = FALSE
)
write.table(parameters, file.path(output_dir, "parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

source_manifest <- data.frame(
  role = c("Vitvi_to_Ensembl_mapping", "GO_annotation", "formal_GSE124820_module_assignment", "passed_GO_mapping_audit"),
  path = required_files,
  bytes = file.info(required_files)$size,
  modified_time = format(file.info(required_files)$mtime, "%Y-%m-%d %H:%M:%S %z"),
  md5 = unname(tools::md5sum(required_files)),
  stringsAsFactors = FALSE
)
write.table(source_manifest, file.path(output_dir, "source_file_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

session_connection <- file(file.path(output_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_connection)
close(session_connection)

escape_markdown <- function(x) {
  x <- ifelse(is.na(x), "", x)
  gsub("\\|", "\\\\|", x)
}

report <- c(
  "# GO enrichment report",
  "",
  "- Status: **COMPLETE**",
  "- Analysis: one-sided hypergeometric over-representation analysis.",
  "- Background: all 5,000 genes in the formal GSE124820 WGCNA network.",
  "- Mapping sources: only the audited authoritative Vitvi-to-Ensembl and grape GO tables.",
  "- Multiple testing: Benjamini-Hochberg within each module-by-ontology family.",
  "- Significance: FDR < 0.05 and at least 3 overlapping module genes.",
  "",
  "## Significant term counts",
  "",
  "| Module | Ontology | Terms tested | Terms with overlap >=3 | Significant terms |",
  "|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(summary_grid)), function(i) {
    sprintf("| %s | %s | %d | %d | %d |", summary_grid$module_color[i], summary_grid$ontology[i],
            summary_grid$terms_tested[i], summary_grid$terms_overlap_ge_3[i], summary_grid$significant_terms[i])
  }, character(1)),
  "",
  "## Top significant terms",
  ""
)

for (module_name in target_modules) {
  report <- c(report, paste0("### ", module_name), "")
  module_significant <- significant[significant$module_color == module_name, , drop = FALSE]
  if (nrow(module_significant) == 0L) {
    report <- c(report, "No terms met both significance criteria.", "")
  } else {
    module_significant <- module_significant[order(module_significant$fdr_bh, -module_significant$fold_enrichment), ]
    top <- head(module_significant, 15L)
    report <- c(report,
      "| Ontology | GO ID | Term | Overlap | Fold enrichment | FDR |",
      "|---|---|---|---:|---:|---:|",
      vapply(seq_len(nrow(top)), function(i) {
        sprintf("| %s | %s | %s | %d | %.3f | %.3g |",
                top$ontology[i], top$go_id[i], escape_markdown(top$go_term[i]), top$module_term_genes[i],
                top$fold_enrichment[i], top$fdr_bh[i])
      }, character(1)), "")
  }
}

report <- c(report,
  "## Quality checks",
  "",
  paste0("- Audited GO-annotated background genes: ", length(unique(gene_term$gene_id)), "."),
  paste0("- Unique tested background GO terms: ", length(unique(combined$go_id)), "."),
  paste0("- Total module-by-ontology tests: ", nrow(combined), "."),
  paste0("- Total significant module-term associations: ", nrow(significant), "."),
  "- All gene-term associations were deduplicated at the `gene_id + GO ID` level before testing.",
  "- Unannotated genes remain in the 5,000-gene universe and module-size denominators, as pre-specified.",
  "",
  "## Output interpretation",
  "",
  "`go_enrichment_all_tests.tsv` preserves the complete tested universe. `go_enrichment_significant.tsv` contains only associations passing both the FDR and overlap thresholds. Per-module/per-ontology files are provided for direct figure generation and manuscript auditing."
)

report_connection <- file(file.path(output_dir, "GO_ENRICHMENT_REPORT.md"), open = "wt", encoding = "UTF-8")
writeLines(report, report_connection, useBytes = TRUE)
close(report_connection)

message("GO enrichment completed")
print(summary_grid, row.names = FALSE)
message("Total significant associations: ", nrow(significant))

