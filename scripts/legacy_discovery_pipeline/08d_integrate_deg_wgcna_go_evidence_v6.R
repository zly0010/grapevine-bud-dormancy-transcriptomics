#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08d_integrate_deg_wgcna_go_evidence.R <project_root> <output_dir>")
project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

files <- c(
  kme = file.path(project_root, "results_corrected/07_wgcna_fixed/GSE124820/kME_table.txt"),
  module_trait = file.path(project_root, "results_corrected/07_wgcna_fixed/GSE124820/module_trait_correlation.txt"),
  preservation = file.path(project_root, "results_corrected/07_wgcna_fixed/module_preservation/preservation_GSE124820_vs_GSE273240.txt"),
  lrt_consensus = file.path(project_root, "results_corrected/02_deseq2_gse124820/cross_variety_consensus/LRT_consensus_table_B.txt"),
  day10_consensus = file.path(project_root, "results_corrected/02_deseq2_gse124820/cross_variety_consensus/Day10_consensus_DEGs_B.txt"),
  mapping = file.path(project_root, "results_corrected/06_gene_annotation/vitvi_to_ensembl_mapping.tsv"),
  go = file.path(project_root, "results_corrected/06_gene_annotation/grape_go_annotation.tsv"),
  go_audit = file.path(project_root, "results_corrected/08_wgcna_postprocessing/01_go_mapping_audit/module_go_coverage_audit.tsv"),
  go_enrichment = file.path(project_root, "results_corrected/08_wgcna_postprocessing/02_go_enrichment/go_enrichment_significant.tsv")
)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0L) stop("Missing required files: ", paste(missing_files, collapse = "; "))

message("Reading corrected DEG, formal WGCNA, preservation and audited GO evidence")
kme <- read.delim(files[["kme"]], check.names = FALSE, quote = "", fileEncoding = "UTF-8")
module_trait <- read.delim(files[["module_trait"]], check.names = FALSE, quote = "", fileEncoding = "UTF-8")
lrt <- read.delim(files[["lrt_consensus"]], check.names = FALSE, quote = "", fileEncoding = "UTF-8")
day10 <- read.delim(files[["day10_consensus"]], check.names = FALSE, quote = "", fileEncoding = "UTF-8")
mapping <- read.delim(files[["mapping"]], check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
go <- read.delim(files[["go"]], check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
go_audit <- read.delim(files[["go_audit"]], check.names = FALSE, quote = "", fileEncoding = "UTF-8")
go_enrichment <- read.delim(files[["go_enrichment"]], check.names = FALSE, quote = "", fileEncoding = "UTF-8")

target_modules <- c("blue", "turquoise", "brown")
if (nrow(kme) != 5000L || anyDuplicated(kme$gene_id) || !all(c("gene_id", "module_color") %in% names(kme))) {
  stop("Formal kME table is not the expected 5,000-gene unique table")
}
expected_sizes <- c(blue = 1045L, turquoise = 2333L, brown = 680L)
observed_sizes <- table(kme$module_color)
if (!all(as.integer(observed_sizes[names(expected_sizes)]) == expected_sizes)) stop("Target module sizes do not match formal results")
if (nrow(lrt) != 22147L || sum(lrt$n_support == 4L) != 13345L) stop("Corrected LRT consensus table fails expected count checks")
if (nrow(day10) != 4210L || anyDuplicated(day10$gene_id)) stop("Corrected Day10 consensus table fails expected count checks")
if (!all(as.logical(go_audit$audit_pass))) stop("GO mapping audit did not pass")
if (nrow(go_enrichment) != 0L) message("Significant module-level GO terms present: ", nrow(go_enrichment))

target <- kme[kme$module_color %in% target_modules, , drop = FALSE]
module_kme_columns <- paste0("ME", target$module_color)
column_index <- match(module_kme_columns, names(target))
if (anyNA(column_index)) stop("Module-specific kME columns are missing")
target$module_kme <- as.numeric(as.matrix(target)[cbind(seq_len(nrow(target)), column_index)])
target$abs_module_kme <- abs(target$module_kme)
target$high_kme_0_8 <- target$abs_module_kme >= 0.8

trait_rows <- module_trait[module_trait$module %in% paste0("ME", target_modules), c("module", "time", "p_BH.time")]
trait_evidence <- data.frame(
  module_color = sub("^ME", "", trait_rows$module),
  module_time_correlation = trait_rows$time,
  module_time_fdr = trait_rows$p_BH.time,
  stringsAsFactors = FALSE
)
if (nrow(trait_evidence) != 3L || anyDuplicated(trait_evidence$module_color)) stop("Module-time evidence is incomplete")

preservation_lines <- readLines(files[["preservation"]], warn = FALSE, encoding = "UTF-8")
candidate_markers <- grep("$Z$ref.Set_1$inColumnsAlsoPresentIn.Set_2", trimws(preservation_lines), fixed = TRUE)
marker <- candidate_markers[vapply(candidate_markers, function(i) {
  i < length(preservation_lines) && grepl("Zsummary.pres", preservation_lines[i + 1L], fixed = TRUE)
}, logical(1))]
if (length(marker) != 1L) stop("Could not uniquely locate the Zsummary.pres section")
preservation_tail <- trimws(preservation_lines[(marker + 2L):length(preservation_lines)])
next_table_header <- grep("^Z\\.propVarExplained\\.pres", preservation_tail)
if (length(next_table_header) < 1L) stop("Could not locate the end of the primary Zsummary.pres table")
preservation_section <- preservation_tail[seq_len(next_table_header[1] - 1L)]
parse_preservation <- function(module_name) {
  line <- grep(paste0("^", module_name, "[[:space:]]+"), preservation_section, value = TRUE)
  if (length(line) != 1L) stop("Could not uniquely parse preservation row for ", module_name)
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1]]
  data.frame(module_color = fields[1], common_gene_module_size = as.integer(fields[2]), preservation_zsummary = as.numeric(fields[3]), stringsAsFactors = FALSE)
}
preservation_evidence <- do.call(rbind, lapply(target_modules, parse_preservation))
expected_z <- c(blue = 27.8399, turquoise = 34.1040, brown = 14.4094)
if (any(abs(preservation_evidence$preservation_zsummary - expected_z[preservation_evidence$module_color]) > 1e-4)) {
  stop("Parsed preservation Zsummary values differ from formal reported values")
}
preservation_evidence$strong_preservation_z_gt_10 <- preservation_evidence$preservation_zsummary > 10

target <- merge(target, trait_evidence, by = "module_color", all.x = TRUE, sort = FALSE)
target <- merge(target, preservation_evidence, by = "module_color", all.x = TRUE, sort = FALSE)
if (anyNA(target[, c("module_time_correlation", "module_time_fdr", "preservation_zsummary")])) stop("Failed to attach module-level evidence")

lrt_evidence <- lrt[, c("gene_id", "n_support")]
names(lrt_evidence)[2] <- "lrt_variety_support"
target <- merge(target, lrt_evidence, by = "gene_id", all.x = TRUE, sort = FALSE)
target$lrt_variety_support[is.na(target$lrt_variety_support)] <- 0L
target$common_time_response_4of4 <- target$lrt_variety_support == 4L

day10_columns <- c("gene_id", "Vamu", "Vvcs", "Vvri", "Vrip", "consensus", "mean_abs_lfc")
if (!all(day10_columns %in% names(day10))) stop("Day10 consensus table lacks expected columns")
day10_evidence <- day10[, day10_columns]
names(day10_evidence)[2:5] <- paste0("day10_log2fc_", names(day10_evidence)[2:5])
names(day10_evidence)[6:7] <- c("day10_consensus_direction", "day10_mean_abs_lfc")
target <- merge(target, day10_evidence, by = "gene_id", all.x = TRUE, sort = FALSE)
target$day10_common_deg <- !is.na(target$day10_consensus_direction)
target$day10_direction_consistent <- target$day10_common_deg & target$day10_consensus_direction %in% c("up", "down")
target$deg_consensus_evidence <- target$common_time_response_4of4 | target$day10_common_deg

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

mapping_valid <- unique(mapping[!is.na(mapping$vitvi_id) & !is.na(mapping$ensembl_id), c("vitvi_id", "ensembl_id")])
go_valid <- unique(go[!is.na(go$ensembl_gene_id) & !is.na(go$go_id), c("ensembl_gene_id", "go_id", "go_term", "go_type")])
gene_go <- merge(mapping_valid, go_valid, by.x = "ensembl_id", by.y = "ensembl_gene_id", all = FALSE, sort = FALSE)
gene_go <- gene_go[gene_go$vitvi_id %in% target$gene_id, ]

collapse_unique <- function(x) {
  values <- sort(unique(x[!is.na(x) & x != ""]))
  if (length(values) == 0L) NA_character_ else paste(values, collapse = ";")
}
if (nrow(gene_go) > 0L) {
  go_by_gene <- do.call(rbind, lapply(split(gene_go, gene_go$vitvi_id), function(x) {
    data.frame(
      gene_id = x$vitvi_id[1], ensembl_id = collapse_unique(x$ensembl_id), go_ids = collapse_unique(x$go_id),
      go_terms = collapse_unique(x$go_term), go_types = collapse_unique(x$go_type), stringsAsFactors = FALSE
    )
  }))
} else {
  go_by_gene <- data.frame(gene_id = character(), ensembl_id = character(), go_ids = character(), go_terms = character(), go_types = character())
}
target <- merge(target, go_by_gene, by = "gene_id", all.x = TRUE, sort = FALSE)
target$go_annotated <- !is.na(target$go_ids)

go_text <- tolower(ifelse(is.na(target$go_terms), "", target$go_terms))
theme_patterns <- c(
  dormancy_bud = "dorman|bud",
  aba = "abscisic|\\baba\\b",
  cold_temperature = "cold|chilling|freez|temperature",
  redox_ros = "redox|oxid|reactive oxygen|peroxid|glutathione",
  cell_wall = "cell wall|pectin|cellulose|xyloglucan|lignin",
  stress = "stress"
)
for (theme_name in names(theme_patterns)) {
  target[[paste0("go_theme_", theme_name)]] <- grepl(theme_patterns[[theme_name]], go_text, perl = TRUE)
}
theme_columns <- paste0("go_theme_", names(theme_patterns))
target$go_priority_theme_count <- rowSums(target[, theme_columns, drop = FALSE])
target$go_priority_theme_any <- target$go_priority_theme_count > 0L

module_sig_counts <- table(factor(go_enrichment$module_color, levels = target_modules))
module_go <- data.frame(module_color = target_modules, module_significant_go_terms = as.integer(module_sig_counts), stringsAsFactors = FALSE)
target <- merge(target, module_go, by = "module_color", all.x = TRUE, sort = FALSE)
target$module_significant_go_enrichment <- target$module_significant_go_terms > 0L

target$core_evidence_eligible <- target$high_kme_0_8 & target$strong_preservation_z_gt_10 &
  (target$module_time_fdr < 0.05) & target$deg_consensus_evidence & target$go_annotated

target <- target[order(match(target$module_color, target_modules), -target$abs_module_kme, target$gene_id), ]
core <- target[target$core_evidence_eligible, , drop = FALSE]

write.table(target, file.path(output_dir, "integrated_gene_evidence_all_target_modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(core, file.path(output_dir, "core_evidence_eligible_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(trait_evidence, file.path(output_dir, "module_time_evidence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(preservation_evidence, file.path(output_dir, "module_preservation_evidence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_rows <- do.call(rbind, lapply(target_modules, function(module_name) {
  x <- target[target$module_color == module_name, ]
  data.frame(
    module_color = module_name, genes = nrow(x), high_kme = sum(x$high_kme_0_8),
    common_time_response_4of4 = sum(x$common_time_response_4of4), day10_common_deg = sum(x$day10_common_deg),
    go_annotated = sum(x$go_annotated), go_priority_theme = sum(x$go_priority_theme_any),
    core_evidence_eligible = sum(x$core_evidence_eligible), stringsAsFactors = FALSE
  )
}))
write.table(summary_rows, file.path(output_dir, "evidence_integration_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameters <- data.frame(
  parameter = c("target_modules", "high_kme_threshold", "strong_preservation_threshold", "module_time_threshold", "deg_evidence", "go_sources", "go_priority_theme_patterns", "core_eligibility"),
  value = c(
    paste(target_modules, collapse = ","), "abs(module-specific kME) >= 0.8", "Zsummary > 10", "BH-adjusted module-time P < 0.05",
    "4/4 LRT time-response support OR membership in corrected 4,210-gene Day10 consensus DEG set",
    "Only audited vitvi_to_ensembl_mapping.tsv and grape_go_annotation.tsv",
    paste(paste(names(theme_patterns), theme_patterns, sep = "="), collapse = ";"),
    "target module AND high kME AND strong preservation AND significant module-time association AND DEG consensus evidence AND GO annotation"
  ), stringsAsFactors = FALSE
)
write.table(parameters, file.path(output_dir, "parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

source_manifest <- data.frame(role = names(files), path = unname(files), bytes = file.info(unname(files))$size,
                              modified_time = format(file.info(unname(files))$mtime, "%Y-%m-%d %H:%M:%S %z"),
                              md5 = unname(tools::md5sum(unname(files))), stringsAsFactors = FALSE)
write.table(source_manifest, file.path(output_dir, "source_file_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
session_connection <- file(file.path(output_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_connection)
close(session_connection)

report <- c(
  "# DEG-WGCNA-GO evidence integration",
  "",
  "- Status: **COMPLETE**",
  "- Scope: blue, turquoise and brown genes from the formal GSE124820 WGCNA network.",
  "- No DESeq2 or WGCNA computation was rerun; only existing corrected/formal outputs were joined.",
  "- GO joins used only the two audited authoritative annotation files.",
  "- Module-level GO enrichment contributed no positive term because no association passed FDR < 0.05; gene-level GO annotations remain available as evidence.",
  "",
  "## Evidence counts",
  "",
  "| Module | Genes | abs(kME)>=0.8 | 4/4 time response | Day10 common DEG | GO annotated | Priority GO theme | Core eligible |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary_rows)), function(i) {
    sprintf("| %s | %d | %d | %d | %d | %d | %d | %d |",
            summary_rows$module_color[i], summary_rows$genes[i], summary_rows$high_kme[i],
            summary_rows$common_time_response_4of4[i], summary_rows$day10_common_deg[i],
            summary_rows$go_annotated[i], summary_rows$go_priority_theme[i], summary_rows$core_evidence_eligible[i])
  }, character(1)),
  "",
  "## Module-level evidence",
  "",
  "| Module | Time correlation | Time FDR | Preservation Zsummary | Strong preservation |",
  "|---|---:|---:|---:|---|",
  vapply(target_modules, function(module_name) {
    trow <- trait_evidence[trait_evidence$module_color == module_name, ]
    prow <- preservation_evidence[preservation_evidence$module_color == module_name, ]
    sprintf("| %s | %.4f | %.3g | %.4f | %s |", module_name, trow$module_time_correlation, trow$module_time_fdr,
            prow$preservation_zsummary, if (prow$strong_preservation_z_gt_10) "yes" else "no")
  }, character(1)),
  "",
  "## Core eligibility rule",
  "",
  "A gene is marked core-eligible only if it belongs to a target module, has `abs(kME) >= 0.8`, is in a strongly preserved and time-associated module, has corrected cross-variety DEG evidence, and has an audited GO annotation. Ranking is intentionally deferred to the next stage."
)
report_connection <- file(file.path(output_dir, "EVIDENCE_INTEGRATION_REPORT.md"), open = "wt", encoding = "UTF-8")
writeLines(report, report_connection, useBytes = TRUE)
close(report_connection)

message("Evidence integration completed")
print(summary_rows, row.names = FALSE)
message("Core evidence-eligible genes: ", nrow(core))
