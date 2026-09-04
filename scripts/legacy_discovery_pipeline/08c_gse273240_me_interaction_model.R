#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 08c_gse273240_me_interaction_model.R <project_root> <output_dir>")
}

project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

me_file <- file.path(project_root, "results_corrected/07_wgcna_fixed/GSE273240/module_eigengenes.txt")
design_file <- file.path(project_root, "results_corrected/03_deseq2_gse273240/sample_design_gse273240.txt")
required_files <- c(me_file, design_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) stop("Missing required files: ", paste(missing_files, collapse = "; "))
if (!requireNamespace("car", quietly = TRUE)) stop("The installed R environment lacks the required 'car' package")

message("Reading formal GSE273240 module eigengenes and sample design")
me <- read.delim(me_file, row.names = 1, check.names = FALSE, quote = "", fileEncoding = "UTF-8")
design <- read.delim(design_file, check.names = FALSE, quote = "", fileEncoding = "UTF-8")

required_design_columns <- c("sample_id", "deac_phase", "day", "treatment", "replicate", "day_factor")
if (!all(required_design_columns %in% names(design))) {
  stop("Design file lacks columns: ", paste(setdiff(required_design_columns, names(design)), collapse = ", "))
}
if (nrow(me) != 90L || nrow(design) != 90L) stop("Expected 90 samples in both eigengene and design tables")
if (anyDuplicated(rownames(me)) || anyDuplicated(design$sample_id)) stop("Duplicated sample IDs detected")
if (!setequal(rownames(me), design$sample_id)) {
  stop("Module eigengene and design sample IDs do not match exactly")
}
me <- me[design$sample_id, , drop = FALSE]
if (anyNA(me) || any(!is.finite(as.matrix(me)))) stop("Module eigengene table contains missing or non-finite values")

formal_me_columns <- setdiff(names(me), "MEgrey")
expected_me_columns <- paste0("ME", c("black", "blue", "brown", "green", "magenta", "pink", "red", "turquoise", "yellow"))
if (!setequal(formal_me_columns, expected_me_columns)) {
  stop("Unexpected formal module eigengene columns: ", paste(formal_me_columns, collapse = ", "))
}
formal_me_columns <- expected_me_columns

phase_levels <- c("Deac1", "Deac2", "Deac3")
day_table <- unique(design[, c("day", "day_factor")])
day_table <- day_table[order(as.numeric(day_table$day)), ]
day_levels <- day_table$day_factor
treatment_levels <- c("Control", "tetralone-ABA")

design$deac_phase <- factor(design$deac_phase, levels = phase_levels)
design$day_factor <- factor(design$day_factor, levels = day_levels)
design$treatment <- factor(design$treatment, levels = treatment_levels)
if (anyNA(design[, c("deac_phase", "day_factor", "treatment")])) stop("Unexpected factor levels in design")

cell_counts <- as.data.frame(xtabs(~ deac_phase + day_factor + treatment, data = design))
names(cell_counts) <- c("deac_phase", "day_factor", "treatment", "n")
if (nrow(cell_counts) != 30L || any(cell_counts$n != 3L)) {
  stop("The 3 x 5 x 2 design is not complete and balanced with three replicates per cell")
}
write.table(cell_counts, file.path(output_dir, "design_balance_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

old_contrasts <- options("contrasts")
on.exit(options(old_contrasts), add = TRUE)
options(contrasts = c("contr.sum", "contr.poly"))

effect_results <- list()
coefficient_results <- list()
sample_results <- list()
cell_results <- list()
contrast_results <- list()
type_comparisons <- list()

effect_index <- 0L
coefficient_index <- 0L
sample_index <- 0L
cell_index <- 0L
contrast_index <- 0L
comparison_index <- 0L

effect_order <- c(
  "deac_phase", "day_factor", "treatment", "deac_phase:day_factor",
  "deac_phase:treatment", "day_factor:treatment", "deac_phase:day_factor:treatment"
)

message("Fitting full-factorial Type III models with sum-to-zero contrasts")
for (me_column in formal_me_columns) {
  module_color <- sub("^ME", "", me_column)
  model_data <- design
  model_data$ME <- me[[me_column]]
  model <- lm(ME ~ deac_phase * day_factor * treatment, data = model_data)
  type3 <- as.data.frame(car::Anova(model, type = 3, test.statistic = "F"))
  type3$effect <- rownames(type3)
  rownames(type3) <- NULL
  names(type3) <- sub("Sum Sq", "sum_sq", names(type3), fixed = TRUE)
  names(type3) <- sub("F value", "f_value", names(type3), fixed = TRUE)
  names(type3) <- sub("Pr(>F)", "p_value", names(type3), fixed = TRUE)
  type3 <- type3[type3$effect %in% effect_order, , drop = FALSE]
  type3 <- type3[match(effect_order, type3$effect), , drop = FALSE]
  residual_sum_sq <- sum(residuals(model)^2)
  residual_df <- df.residual(model)
  type3$module_color <- module_color
  type3$module_eigengene <- me_column
  type3$residual_sum_sq <- residual_sum_sq
  type3$residual_df <- residual_df
  type3$partial_eta_squared <- type3$sum_sq / (type3$sum_sq + residual_sum_sq)
  effect_index <- effect_index + 1L
  effect_results[[effect_index]] <- type3[, c(
    "module_color", "module_eigengene", "effect", "sum_sq", "Df", "f_value", "p_value",
    "partial_eta_squared", "residual_sum_sq", "residual_df"
  )]

  type1 <- as.data.frame(anova(model))
  type1$effect <- rownames(type1)
  rownames(type1) <- NULL
  names(type1) <- sub("Pr(>F)", "type1_p_value", names(type1), fixed = TRUE)
  type1 <- type1[type1$effect %in% effect_order, c("effect", "type1_p_value"), drop = FALSE]
  comparison <- merge(type3[, c("effect", "p_value")], type1, by = "effect", all = TRUE, sort = FALSE)
  comparison$module_color <- module_color
  comparison$absolute_p_difference <- abs(comparison$p_value - comparison$type1_p_value)
  comparison_index <- comparison_index + 1L
  type_comparisons[[comparison_index]] <- comparison[, c("module_color", "effect", "p_value", "type1_p_value", "absolute_p_difference")]

  coefficients <- as.data.frame(summary(model)$coefficients)
  coefficients$coefficient <- rownames(coefficients)
  rownames(coefficients) <- NULL
  names(coefficients) <- c("estimate", "std_error", "t_value", "p_value", "coefficient")
  coefficients$module_color <- module_color
  coefficient_index <- coefficient_index + 1L
  coefficient_results[[coefficient_index]] <- coefficients[, c("module_color", "coefficient", "estimate", "std_error", "t_value", "p_value")]

  sample_output <- model_data[, c("sample_id", "deac_phase", "day", "day_factor", "treatment", "replicate")]
  sample_output$module_color <- module_color
  sample_output$module_eigengene <- me_column
  sample_output$observed_me <- model_data$ME
  sample_output$fitted_me <- fitted(model)
  sample_output$residual <- residuals(model)
  sample_index <- sample_index + 1L
  sample_results[[sample_index]] <- sample_output

  group_key <- interaction(model_data$deac_phase, model_data$day_factor, model_data$treatment, drop = TRUE)
  split_values <- split(model_data$ME, group_key)
  split_rows <- split(seq_len(nrow(model_data)), group_key)
  cell_output <- do.call(rbind, lapply(names(split_values), function(key) {
    rows <- split_rows[[key]]
    values <- split_values[[key]]
    data.frame(
      module_color = module_color,
      deac_phase = as.character(model_data$deac_phase[rows[1]]),
      day = model_data$day[rows[1]],
      day_factor = as.character(model_data$day_factor[rows[1]]),
      treatment = as.character(model_data$treatment[rows[1]]),
      n = length(values), mean_me = mean(values), sd_me = sd(values), se_me = sd(values) / sqrt(length(values)),
      stringsAsFactors = FALSE
    )
  }))
  cell_index <- cell_index + 1L
  cell_results[[cell_index]] <- cell_output

  mse <- deviance(model) / df.residual(model)
  phases <- levels(model_data$deac_phase)
  days <- levels(model_data$day_factor)
  contrast_output <- do.call(rbind, lapply(phases, function(phase) {
    do.call(rbind, lapply(days, function(day_name) {
      control_values <- model_data$ME[model_data$deac_phase == phase & model_data$day_factor == day_name & model_data$treatment == "Control"]
      treated_values <- model_data$ME[model_data$deac_phase == phase & model_data$day_factor == day_name & model_data$treatment == "tetralone-ABA"]
      difference <- mean(treated_values) - mean(control_values)
      standard_error <- sqrt(mse * (1 / length(control_values) + 1 / length(treated_values)))
      t_value <- difference / standard_error
      p_value <- 2 * pt(abs(t_value), df = df.residual(model), lower.tail = FALSE)
      data.frame(
        module_color = module_color, deac_phase = phase, day_factor = day_name,
        contrast = "tetralone-ABA_minus_Control", estimate = difference,
        std_error = standard_error, t_value = t_value, df = df.residual(model), p_value = p_value,
        stringsAsFactors = FALSE
      )
    }))
  }))
  contrast_index <- contrast_index + 1L
  contrast_results[[contrast_index]] <- contrast_output
}

effects <- do.call(rbind, effect_results)
effects$fdr_bh_global <- p.adjust(effects$p_value, method = "BH")
effects$significant_fdr_0_05 <- effects$fdr_bh_global < 0.05
effects <- effects[order(match(effects$effect, effect_order), effects$fdr_bh_global, effects$module_color), ]

coefficients <- do.call(rbind, coefficient_results)
samples <- do.call(rbind, sample_results)
cells <- do.call(rbind, cell_results)
contrasts <- do.call(rbind, contrast_results)
contrasts$fdr_bh_global <- p.adjust(contrasts$p_value, method = "BH")
contrasts$significant_fdr_0_05 <- contrasts$fdr_bh_global < 0.05
comparisons <- do.call(rbind, type_comparisons)

max_type_p_difference <- max(comparisons$absolute_p_difference, na.rm = TRUE)
orthogonality_check_pass <- is.finite(max_type_p_difference) && max_type_p_difference < 1e-10
if (!orthogonality_check_pass) {
  stop("Type I and Type III p-values differ in the supposedly balanced design; max difference = ", max_type_p_difference)
}

primary_effects <- c("day_factor:treatment", "deac_phase:day_factor:treatment")
primary <- effects[effects$effect %in% primary_effects, , drop = FALSE]

write.table(effects, file.path(output_dir, "me_type3_effect_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(primary, file.path(output_dir, "primary_interaction_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coefficients, file.path(output_dir, "model_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples, file.path(output_dir, "sample_fitted_values_and_residuals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cells, file.path(output_dir, "module_eigengene_cell_means.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contrasts, file.path(output_dir, "phase_day_treatment_contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comparisons, file.path(output_dir, "type1_type3_orthogonality_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameters <- data.frame(
  parameter = c(
    "response", "model", "day_encoding", "modules_tested", "grey_module", "anova_type", "factor_contrasts",
    "effect_multiple_testing", "effect_correction_family_size", "primary_effects", "simple_contrast_multiple_testing",
    "significance_threshold", "balance_requirement"
  ),
  value = c(
    "Formal GSE273240 module eigengenes", "ME ~ deac_phase * day_factor * treatment",
    paste0("Categorical, chronological levels: ", paste(day_levels, collapse = ",")),
    paste(sub("^ME", "", formal_me_columns), collapse = ","), "Excluded from biological hypothesis tests",
    "Type III F tests via car::Anova", "contr.sum for unordered factors",
    "Benjamini-Hochberg across all formal modules and all seven model effects", as.character(nrow(effects)),
    paste(primary_effects, collapse = ","), "Benjamini-Hochberg across all module x phase x day treatment contrasts",
    "FDR < 0.05", "Complete 3 x 5 x 2 design with exactly 3 replicates per cell"
  ),
  stringsAsFactors = FALSE
)
write.table(parameters, file.path(output_dir, "parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

source_manifest <- data.frame(
  role = c("formal_GSE273240_module_eigengenes", "formal_GSE273240_sample_design"),
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

significant_effects <- effects[effects$significant_fdr_0_05, , drop = FALSE]
significant_primary <- primary[primary$significant_fdr_0_05, , drop = FALSE]

report <- c(
  "# GSE273240 module eigengene interaction model",
  "",
  "- Status: **COMPLETE**",
  "- Input: formal module eigengenes from `results_corrected/07_wgcna_fixed/GSE273240`; the WGCNA network was not rerun.",
  "- Model: `ME ~ deac_phase * day_factor * treatment`.",
  "- Day treatment: categorical to retain the non-linear five-time-point design.",
  "- Test: Type III F tests with sum-to-zero contrasts in the complete balanced design.",
  "- Multiplicity: BH correction across all 9 formal modules and all 7 effects (63 tests).",
  "- Grey was excluded from biological tests.",
  "",
  "## Design audit",
  "",
  "- Samples: 90; complete cells: 30; replicates per cell: 3.",
  paste0("- Phase levels: ", paste(phase_levels, collapse = ", "), "."),
  paste0("- Day levels: ", paste(day_levels, collapse = ", "), "."),
  paste0("- Treatment levels: ", paste(treatment_levels, collapse = ", "), "."),
  sprintf("- Type I/Type III orthogonality check: PASS (maximum absolute P-value difference %.3g).", max_type_p_difference),
  "",
  "## Significant omnibus effects",
  ""
)

if (nrow(significant_effects) == 0L) {
  report <- c(report, "No omnibus effect passed global FDR < 0.05.", "")
} else {
  report <- c(report,
    "| Module | Effect | F | Partial eta-squared | P | Global FDR |",
    "|---|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(significant_effects)), function(i) {
      sprintf("| %s | %s | %.4f | %.4f | %.3g | %.3g |",
              significant_effects$module_color[i], significant_effects$effect[i], significant_effects$f_value[i],
              significant_effects$partial_eta_squared[i], significant_effects$p_value[i], significant_effects$fdr_bh_global[i])
    }, character(1)), "")
}

report <- c(report, "## Pre-specified interaction focus", "")
if (nrow(significant_primary) == 0L) {
  report <- c(report, "Neither the day-by-treatment nor the phase-by-day-by-treatment effect passed the global FDR threshold in any module.", "")
} else {
  report <- c(report,
    "| Module | Effect | F | Partial eta-squared | P | Global FDR |",
    "|---|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(significant_primary)), function(i) {
      sprintf("| %s | %s | %.4f | %.4f | %.3g | %.3g |",
              significant_primary$module_color[i], significant_primary$effect[i], significant_primary$f_value[i],
              significant_primary$partial_eta_squared[i], significant_primary$p_value[i], significant_primary$fdr_bh_global[i])
    }, character(1)), "")
}

report <- c(report,
  "## Interpretation guardrails",
  "",
  "- A treatment main effect is not interpreted in isolation when a treatment interaction is significant.",
  "- Cell means and phase-by-day treatment contrasts are saved for directionality; their contrast FDR family is separate from the 63 omnibus tests.",
  "- Module eigengene sign is mathematically arbitrary, so biological direction is interpreted jointly with module gene-expression patterns, not sign alone."
)

report_connection <- file(file.path(output_dir, "GSE273240_ME_INTERACTION_REPORT.md"), open = "wt", encoding = "UTF-8")
writeLines(report, report_connection, useBytes = TRUE)
close(report_connection)

message("GSE273240 module eigengene interaction models completed")
message("Significant omnibus effects: ", nrow(significant_effects), " / ", nrow(effects))
message("Significant primary interactions: ", nrow(significant_primary), " / ", nrow(primary))
print(primary[, c("module_color", "effect", "f_value", "p_value", "fdr_bh_global", "significant_fdr_0_05")], row.names = FALSE)

