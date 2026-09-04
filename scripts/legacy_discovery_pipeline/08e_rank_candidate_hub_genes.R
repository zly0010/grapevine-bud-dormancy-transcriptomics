#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08e_rank_candidate_hub_genes.R <project_root> <output_dir>")
project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

integration_file <- file.path(project_root, "results_corrected/08_wgcna_postprocessing/04_evidence_integration_v6/integrated_gene_evidence_all_target_modules.tsv")
integration_report <- file.path(project_root, "results_corrected/08_wgcna_postprocessing/04_evidence_integration_v6/EVIDENCE_INTEGRATION_REPORT.md")
required_files <- c(integration_file, integration_report)
if (!all(file.exists(required_files))) stop("Required evidence integration outputs are missing")

message("Reading passed integrated evidence table")
x <- read.delim(integration_file, check.names = FALSE, quote = "", na.strings = c("", "NA"), fileEncoding = "UTF-8")
required_columns <- c(
  "gene_id", "module_color", "abs_module_kme", "high_kme_0_8", "strong_preservation_z_gt_10",
  "module_time_correlation", "module_time_fdr", "common_time_response_4of4", "day10_common_deg",
  "day10_direction_consistent", "go_annotated", "go_priority_theme_any", "core_evidence_eligible",
  "module_significant_go_enrichment", "go_terms"
)
if (!all(required_columns %in% names(x))) stop("Integrated evidence table lacks required columns")
if (nrow(x) != 4058L || anyDuplicated(x$gene_id)) stop("Integrated evidence table fails row/uniqueness checks")

core <- x[as.logical(x$core_evidence_eligible), , drop = FALSE]
if (nrow(core) != 464L) stop("Expected 464 core evidence-eligible genes but found ", nrow(core))
logical_columns <- c(
  "high_kme_0_8", "strong_preservation_z_gt_10", "common_time_response_4of4", "day10_common_deg",
  "day10_direction_consistent", "go_annotated", "go_priority_theme_any", "module_significant_go_enrichment",
  "core_evidence_eligible"
)
for (column in logical_columns) core[[column]] <- as.logical(core[[column]])

weights <- c(
  abs_kme = 35,
  common_time_response_4of4 = 15,
  day10_common_deg = 10,
  day10_direction_consistent = 5,
  strong_module_preservation = 10,
  absolute_module_time_correlation = 10,
  go_annotation = 7.5,
  priority_go_theme = 7.5
)

core$score_abs_kme <- weights[["abs_kme"]] * pmin(core$abs_module_kme, 1)
core$score_common_time_response <- weights[["common_time_response_4of4"]] * core$common_time_response_4of4
core$score_day10_common_deg <- weights[["day10_common_deg"]] * core$day10_common_deg
core$score_day10_direction_consistent <- weights[["day10_direction_consistent"]] * core$day10_direction_consistent
core$score_strong_preservation <- weights[["strong_module_preservation"]] * core$strong_preservation_z_gt_10
core$score_module_time_correlation <- weights[["absolute_module_time_correlation"]] * pmin(abs(core$module_time_correlation), 1)
core$score_go_annotation <- weights[["go_annotation"]] * core$go_annotated
core$score_priority_go_theme <- weights[["priority_go_theme"]] * core$go_priority_theme_any

score_columns <- grep("^score_", names(core), value = TRUE)
core$composite_score_0_100 <- rowSums(core[, score_columns, drop = FALSE])
if (anyNA(core$composite_score_0_100) || any(core$composite_score_0_100 < 0 | core$composite_score_0_100 > 100)) {
  stop("Composite scores are missing or outside 0-100")
}

theme_columns <- grep("^go_theme_", names(core), value = TRUE)
theme_labels <- sub("^go_theme_", "", theme_columns)
core$priority_go_theme_labels <- vapply(seq_len(nrow(core)), function(i) {
  hits <- theme_labels[as.logical(core[i, theme_columns, drop = TRUE])]
  if (length(hits) == 0L) "" else paste(hits, collapse = ";")
}, character(1))

core <- core[order(-core$composite_score_0_100, -core$abs_module_kme, core$gene_id), ]
core$overall_rank <- seq_len(nrow(core))
core$within_module_rank <- ave(core$composite_score_0_100, core$module_color,
                               FUN = function(v) rank(-v, ties.method = "first"))
core$within_module_rank <- as.integer(core$within_module_rank)

front_columns <- c(
  "overall_rank", "within_module_rank", "gene_id", "module_color", "composite_score_0_100",
  "abs_module_kme", "module_kme", "preservation_zsummary", "module_time_correlation", "module_time_fdr",
  "common_time_response_4of4", "day10_common_deg", "day10_consensus_direction", "day10_mean_abs_lfc",
  "go_annotated", "go_priority_theme_any", "priority_go_theme_labels", "ensembl_id", "go_ids", "go_terms"
)
front_columns <- front_columns[front_columns %in% names(core)]
remaining_columns <- setdiff(names(core), front_columns)
core <- core[, c(front_columns, remaining_columns)]

top_50 <- head(core, 50L)
top_by_module <- do.call(rbind, lapply(c("blue", "turquoise", "brown"), function(module_name) {
  head(core[core$module_color == module_name, , drop = FALSE], 20L)
}))

write.table(core, file.path(output_dir, "candidate_hub_gene_ranking_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(top_50, file.path(output_dir, "top_50_candidate_hub_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(top_by_module, file.path(output_dir, "top_20_candidate_hubs_by_module.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

weight_table <- data.frame(component = names(weights), maximum_points = as.numeric(weights), stringsAsFactors = FALSE)
write.table(weight_table, file.path(output_dir, "ranking_weights.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_rows <- do.call(rbind, lapply(c("blue", "turquoise", "brown"), function(module_name) {
  z <- core[core$module_color == module_name, ]
  data.frame(
    module_color = module_name, ranked_candidates = nrow(z),
    median_score = median(z$composite_score_0_100), maximum_score = max(z$composite_score_0_100),
    priority_theme_candidates = sum(z$go_priority_theme_any), both_deg_evidence = sum(z$common_time_response_4of4 & z$day10_common_deg),
    stringsAsFactors = FALSE
  )
}))
write.table(summary_rows, file.path(output_dir, "ranking_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameters <- data.frame(
  parameter = c("input_population", "eligibility", "ranking_direction", "tie_breakers", "module_go_enrichment", "score_range", "interpretation"),
  value = c(
    "464 core evidence-eligible genes from evidence integration v6",
    "abs(kME)>=0.8; strongly preserved and time-associated target module; corrected consensus DEG evidence; audited GO annotation",
    "Descending composite score", "Descending abs(kME), then ascending gene_id",
    "No positive points because zero module-term associations passed FDR < 0.05",
    "0 to 100", "Prioritization score, not a posterior probability or causal effect estimate"
  ), stringsAsFactors = FALSE
)
write.table(parameters, file.path(output_dir, "parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

source_manifest <- data.frame(role = c("integrated_evidence_table", "integration_report"), path = required_files,
                              bytes = file.info(required_files)$size,
                              modified_time = format(file.info(required_files)$mtime, "%Y-%m-%d %H:%M:%S %z"),
                              md5 = unname(tools::md5sum(required_files)), stringsAsFactors = FALSE)
write.table(source_manifest, file.path(output_dir, "source_file_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
session_connection <- file(file.path(output_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_connection)
close(session_connection)

escape_md <- function(x) gsub("\\|", "\\\\|", ifelse(is.na(x), "", x))
report <- c(
  "# Candidate hub gene ranking",
  "",
  "- Status: **COMPLETE**",
  "- Ranked population: 464 genes passing all pre-specified core eligibility criteria.",
  "- Score: transparent 0-100 prioritization index; it is not a probability and does not establish causality.",
  "- Module-level GO enrichment supplied no positive score because no term passed FDR < 0.05.",
  "",
  "## Score weights",
  "",
  "| Component | Maximum points |",
  "|---|---:|",
  vapply(seq_len(nrow(weight_table)), function(i) sprintf("| %s | %.1f |", weight_table$component[i], weight_table$maximum_points[i]), character(1)),
  "",
  "## Ranking summary",
  "",
  "| Module | Candidates | Median score | Maximum score | Priority GO theme | Both DEG evidence types |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary_rows)), function(i) {
    sprintf("| %s | %d | %.3f | %.3f | %d | %d |", summary_rows$module_color[i], summary_rows$ranked_candidates[i],
            summary_rows$median_score[i], summary_rows$maximum_score[i], summary_rows$priority_theme_candidates[i], summary_rows$both_deg_evidence[i])
  }, character(1)),
  ""
)

for (module_name in c("blue", "turquoise", "brown")) {
  report <- c(report, paste0("## Top 10: ", module_name), "")
  top <- head(core[core$module_color == module_name, , drop = FALSE], 10L)
  report <- c(report,
    "| Overall rank | Gene | Score | abs(kME) | 4/4 time | Day10 | Direction | GO priority theme | GO terms |",
    "|---:|---|---:|---:|---|---|---|---|---|",
    vapply(seq_len(nrow(top)), function(i) {
      terms <- top$go_terms[i]
      if (!is.na(terms) && nchar(terms) > 120L) terms <- paste0(substr(terms, 1L, 117L), "...")
      sprintf("| %d | %s | %.3f | %.3f | %s | %s | %s | %s | %s |",
              top$overall_rank[i], top$gene_id[i], top$composite_score_0_100[i], top$abs_module_kme[i],
              if (top$common_time_response_4of4[i]) "yes" else "no", if (top$day10_common_deg[i]) "yes" else "no",
              ifelse(is.na(top$day10_consensus_direction[i]), "", top$day10_consensus_direction[i]),
              escape_md(top$priority_go_theme_labels[i]), escape_md(terms))
    }, character(1)), "")
}

report <- c(report,
  "## Guardrails",
  "",
  "- Rankings use only corrected DEG outputs, formal WGCNA/preservation outputs, and audited GO mappings.",
  "- GO keyword themes are transparent text tags derived from authoritative GO term names; they are not de novo annotations.",
  "- Candidates without GO annotation were excluded by the pre-specified core rule but remain in the full integration table for audit.",
  "- Experimental validation is required before describing any ranked gene as causal."
)
report_connection <- file(file.path(output_dir, "CANDIDATE_HUB_GENE_RANKING_REPORT.md"), open = "wt", encoding = "UTF-8")
writeLines(report, report_connection, useBytes = TRUE)
close(report_connection)

message("Candidate hub ranking completed")
print(summary_rows, row.names = FALSE)
print(top_50[, c("overall_rank", "gene_id", "module_color", "composite_score_0_100", "abs_module_kme", "priority_go_theme_labels")], row.names = FALSE)

