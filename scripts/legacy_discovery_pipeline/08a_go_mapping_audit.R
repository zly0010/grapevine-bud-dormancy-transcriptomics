#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 08a_go_mapping_audit.R <project_root> <output_dir>")
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

required_files <- c(mapping_file, go_file, module_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Required files are missing: ", paste(missing_files, collapse = "; "))
}

message("Reading authoritative mapping and annotation files")
mapping <- read.delim(mapping_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
go <- read.delim(go_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
modules <- read.delim(module_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")

required_mapping_columns <- c("vitvi_id", "ensembl_id")
required_go_columns <- c("ensembl_gene_id", "go_id", "go_term", "go_type")
required_module_columns <- c("gene_id", "module_color")
if (!all(required_mapping_columns %in% names(mapping))) {
  stop("Mapping file lacks required columns: ", paste(setdiff(required_mapping_columns, names(mapping)), collapse = ", "))
}
if (!all(required_go_columns %in% names(go))) {
  stop("GO file lacks required columns: ", paste(setdiff(required_go_columns, names(go)), collapse = ", "))
}
if (!all(required_module_columns %in% names(modules))) {
  stop("Module file lacks required columns: ", paste(setdiff(required_module_columns, names(modules)), collapse = ", "))
}

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

mapping$vitvi_id <- clean_id(mapping$vitvi_id)
mapping$ensembl_id <- clean_id(mapping$ensembl_id)
go$ensembl_gene_id <- clean_id(go$ensembl_gene_id)
go$go_id <- clean_id(go$go_id)
go$go_term <- clean_id(go$go_term)
go$go_type <- clean_id(go$go_type)
modules$gene_id <- clean_id(modules$gene_id)
modules$module_color <- clean_id(modules$module_color)

if (anyNA(modules$gene_id) || anyNA(modules$module_color)) {
  stop("Formal module assignment contains missing gene_id or module_color")
}
if (anyDuplicated(modules$gene_id)) {
  stop("Formal module assignment contains duplicated gene_id values")
}
if (nrow(modules) != 5000L) {
  stop("Expected 5,000 formal WGCNA genes but found ", nrow(modules))
}

mapping_valid <- mapping[!is.na(mapping$vitvi_id), c("vitvi_id", "ensembl_id"), drop = FALSE]
go_valid <- go[!is.na(go$ensembl_gene_id) & !is.na(go$go_id) & !is.na(go$go_type), required_go_columns, drop = FALSE]

valid_go_types <- c("P", "F", "C")
unexpected_go_types <- setdiff(unique(go_valid$go_type), valid_go_types)
if (length(unexpected_go_types) > 0L) {
  stop("Unexpected GO type values: ", paste(unexpected_go_types, collapse = ", "))
}

mapping_split <- split(mapping_valid$ensembl_id, mapping_valid$vitvi_id)
mapping_split <- lapply(mapping_split, function(x) sort(unique(x[!is.na(x)])))
go_gene_set <- unique(go_valid$ensembl_gene_id)

collapse_values <- function(x) {
  if (length(x) == 0L) return(NA_character_)
  paste(x, collapse = ";")
}

audit <- modules[, c("gene_id", "module_color"), drop = FALSE]
mapped_ids <- lapply(audit$gene_id, function(g) mapping_split[[g]])
mapped_ids <- lapply(mapped_ids, function(x) if (is.null(x)) character(0) else x)
go_ids <- lapply(mapped_ids, function(x) intersect(x, go_gene_set))

audit$mapping_row_present <- audit$gene_id %in% mapping_valid$vitvi_id
audit$successfully_mapped <- lengths(mapped_ids) > 0L
audit$ensembl_id <- vapply(mapped_ids, collapse_values, character(1))
audit$go_annotated <- lengths(go_ids) > 0L
audit$go_annotated_ensembl_id <- vapply(go_ids, collapse_values, character(1))

target_modules <- c("blue", "turquoise", "brown")
expected <- data.frame(
  module_color = target_modules,
  expected_module_genes = c(1045L, 2333L, 680L),
  expected_mapped_genes = c(912L, 2123L, 558L),
  expected_go_annotated_genes = c(250L, 518L, 121L),
  expected_coverage_percent_1dp = c(23.9, 22.2, 17.8),
  stringsAsFactors = FALSE
)

summarize_module <- function(module_name) {
  x <- audit[audit$module_color == module_name, , drop = FALSE]
  data.frame(
    module_color = module_name,
    module_genes = nrow(x),
    mapping_row_present_genes = sum(x$mapping_row_present),
    successfully_mapped_genes = sum(x$successfully_mapped),
    go_annotated_genes = sum(x$go_annotated),
    mapping_rate_percent = 100 * mean(x$successfully_mapped),
    go_coverage_percent = 100 * mean(x$go_annotated),
    stringsAsFactors = FALSE
  )
}

coverage <- do.call(rbind, lapply(target_modules, summarize_module))
coverage <- merge(coverage, expected, by = "module_color", sort = FALSE)
coverage <- coverage[match(target_modules, coverage$module_color), ]
coverage$module_size_match <- coverage$module_genes == coverage$expected_module_genes
coverage$mapped_count_match <- coverage$successfully_mapped_genes == coverage$expected_mapped_genes
coverage$go_count_match <- coverage$go_annotated_genes == coverage$expected_go_annotated_genes
coverage$coverage_1dp_match <- round(coverage$go_coverage_percent, 1) == coverage$expected_coverage_percent_1dp
coverage$audit_pass <- coverage$module_size_match & coverage$mapped_count_match & coverage$go_count_match & coverage$coverage_1dp_match

background <- data.frame(
  set = "GSE124820_formal_WGCNA_background",
  genes = nrow(audit),
  mapping_row_present_genes = sum(audit$mapping_row_present),
  successfully_mapped_genes = sum(audit$successfully_mapped),
  go_annotated_genes = sum(audit$go_annotated),
  mapping_rate_percent = 100 * mean(audit$successfully_mapped),
  go_coverage_percent = 100 * mean(audit$go_annotated),
  stringsAsFactors = FALSE
)

source_manifest <- data.frame(
  role = c("Vitvi_to_Ensembl_mapping", "GO_annotation", "formal_GSE124820_module_assignment"),
  path = required_files,
  bytes = file.info(required_files)$size,
  modified_time = format(file.info(required_files)$mtime, "%Y-%m-%d %H:%M:%S %z"),
  md5 = unname(tools::md5sum(required_files)),
  stringsAsFactors = FALSE
)

parameters <- data.frame(
  parameter = c(
    "mapping_chain", "target_modules", "background_definition", "coverage_denominator",
    "expected_blue_coverage_percent_1dp", "expected_turquoise_coverage_percent_1dp",
    "expected_brown_coverage_percent_1dp", "pass_rule"
  ),
  value = c(
    "gene_id -> vitvi_id -> ensembl_id -> ensembl_gene_id -> GO term",
    paste(target_modules, collapse = ","),
    "All 5,000 genes in formal GSE124820 module_assignments.txt",
    "All genes assigned to the named module, including genes without mapping/GO",
    "23.9", "22.2", "17.8",
    "Exact module/mapped/GO counts and one-decimal coverage must match the project summary"
  ),
  stringsAsFactors = FALSE
)

write.table(audit, file.path(output_dir, "gene_level_mapping_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(coverage, file.path(output_dir, "module_go_coverage_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(background, file.path(output_dir, "background_go_coverage_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(source_manifest, file.path(output_dir, "source_file_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(parameters, file.path(output_dir, "parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

overall_pass <- all(coverage$audit_pass)
report <- c(
  "# GO mapping coverage audit",
  "",
  paste0("- Status: **", if (overall_pass) "PASS" else "FAIL", "**"),
  "- Scope: formal GSE124820 WGCNA only; no DESeq2 or WGCNA was rerun.",
  "- Mapping chain: `gene_id -> vitvi_id -> ensembl_id -> ensembl_gene_id -> GO term`.",
  "- Sources: only the two authoritative files in `results_corrected/06_gene_annotation` plus formal module assignments in `results_corrected/07_wgcna_fixed`.",
  "- Background for later enrichment: all 5,000 genes in the formal GSE124820 network.",
  "",
  "## Coverage results",
  "",
  "| Module | Genes | Successfully mapped | GO-annotated | GO coverage | Audit |",
  "|---|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(coverage)), function(i) {
    sprintf("| %s | %d | %d | %d | %.3f%% | %s |",
            coverage$module_color[i], coverage$module_genes[i], coverage$successfully_mapped_genes[i],
            coverage$go_annotated_genes[i], coverage$go_coverage_percent[i],
            if (coverage$audit_pass[i]) "PASS" else "FAIL")
  }, character(1)),
  "",
  "## Background coverage",
  "",
  sprintf("The 5,000-gene background contains %d successfully mapped genes and %d GO-annotated genes (%.3f%% coverage).",
          background$successfully_mapped_genes, background$go_annotated_genes, background$go_coverage_percent),
  "",
  "## Integrity checks",
  "",
  paste0("- Formal module assignment rows: ", nrow(modules), "; unique gene IDs: ", length(unique(modules$gene_id)), "."),
  paste0("- Mapping rows: ", nrow(mapping), "; non-empty Vitvi IDs: ", sum(!is.na(mapping$vitvi_id)), "; non-empty Ensembl mappings: ", sum(!is.na(mapping$ensembl_id)), "."),
  paste0("- GO rows: ", nrow(go), "; unique annotated Ensembl genes: ", length(unique(go_valid$ensembl_gene_id)), "; unique GO terms: ", length(unique(go_valid$go_id)), "."),
  paste0("- Observed GO types: ", paste(sort(unique(go_valid$go_type)), collapse = ", "), "."),
  "- File sizes, modification times and MD5 checksums are recorded in `source_file_manifest.tsv`.",
  "",
  "## Decision",
  "",
  if (overall_pass) {
    "The mapping audit passed. BP/MF/CC over-representation analysis may proceed using the audited 5,000-gene background."
  } else {
    "The mapping audit failed. Enrichment must not proceed until the discrepancy is resolved."
  }
)

report_connection <- file(file.path(output_dir, "GO_MAPPING_AUDIT_REPORT.md"), open = "wt", encoding = "UTF-8")
writeLines(report, report_connection, useBytes = TRUE)
close(report_connection)

message("Audit status: ", if (overall_pass) "PASS" else "FAIL")
print(coverage[, c("module_color", "module_genes", "successfully_mapped_genes", "go_annotated_genes", "go_coverage_percent", "audit_pass")], row.names = FALSE)
print(background, row.names = FALSE)

if (!overall_pass) {
  quit(save = "no", status = 2L)
}

