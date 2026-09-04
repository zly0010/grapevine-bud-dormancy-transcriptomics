#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 5)
suppressPackageStartupMessages(library(WGCNA))
allowWGCNAThreads(nThreads = 4)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: 12_wgcna_resampling_stability.R <source_project> <output_dir> [bootstrap_replicates]")
}

source_project <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
n_boot <- if (length(args) == 3) as.integer(args[[3]]) else 100L
if (!is.finite(n_boot) || n_boot < 10) stop("bootstrap_replicates must be >= 10")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260718)

vst_path <- file.path(source_project, "results_corrected", "07_wgcna_fixed", "GSE124820_vst_fixed.txt")
design_path <- file.path(source_project, "results_corrected", "01_sample_design_and_qc", "sample_design.tsv")
module_path <- file.path(source_project, "results_corrected", "07_wgcna_fixed", "GSE124820", "module_assignments.txt")
for (p in c(vst_path, design_path, module_path)) if (!file.exists(p)) stop("Missing input: ", p)

message("Loading fixed GSE124820 VST matrix and reference modules")
vst <- read.delim(vst_path, row.names = 1, check.names = FALSE)
design <- read.delim(design_path, check.names = FALSE)
modules <- read.delim(module_path, check.names = FALSE)

if (nrow(modules) != 5000 || anyDuplicated(modules$gene_id)) {
  stop("Reference module assignment is not the expected unique 5,000-gene set")
}
if (!all(modules$gene_id %in% colnames(vst))) stop("Reference genes are missing from VST")

design <- design[design$sample_id %in% rownames(vst), , drop = FALSE]
design <- design[match(rownames(vst), design$sample_id), , drop = FALSE]
if (!identical(rownames(vst), design$sample_id)) stop("VST/design sample alignment failed")
qc_b <- tolower(as.character(design$qc_B)) %in% c("true", "t", "1")
if (!all(qc_b)) stop("Fixed VST contains samples outside formal QC-B set")

dat <- as.matrix(vst[, modules$gene_id, drop = FALSE])
storage.mode(dat) <- "double"
target_modules <- c("blue", "turquoise", "brown")
module_genes <- setNames(lapply(target_modules, function(m) {
  modules$gene_id[modules$module_color == m]
}), target_modules)

module_metrics <- function(x, time_value, module_name, replicate_type, replicate_id) {
  z <- scale(x)
  if (any(!is.finite(z))) stop("Non-finite scaled expression in ", replicate_type, " ", replicate_id)
  p <- ncol(z)
  n <- nrow(z)

  # Sum of the correlation matrix without explicitly allocating p x p memory.
  total_cor <- sum(rowSums(z)^2) / (n - 1)
  mean_pairwise_cor <- (total_cor - p) / (p * (p - 1))

  pc <- irlba::irlba(z, nv = 1, nu = 1)
  eigengene <- pc$u[, 1] * pc$d[[1]]
  score <- rowMeans(z)
  if (cor(eigengene, score) < 0) eigengene <- -eigengene
  kme <- as.numeric(cor(z, eigengene))
  pc1_variance <- pc$d[[1]]^2 / sum(z^2)

  data.frame(
    replicate_type = replicate_type,
    replicate_id = as.character(replicate_id),
    module_color = module_name,
    samples = n,
    genes = p,
    mean_pairwise_correlation = mean_pairwise_cor,
    median_abs_kme = median(abs(kme)),
    pc1_variance_explained = pc1_variance,
    time_correlation = cor(eigengene, time_value, method = "pearson")
  )
}

reference_metrics <- do.call(rbind, lapply(target_modules, function(m) {
  module_metrics(dat[, module_genes[[m]], drop = FALSE], design$time, m, "reference", "all")
}))

message("Running ", n_boot, " stratified bootstrap module-coherence replicates")
strata <- interaction(design$variety, design$time, drop = TRUE)
strata_indices <- split(seq_len(nrow(design)), strata)
bootstrap_rows <- vector("list", n_boot * length(target_modules))
cursor <- 1L
for (b in seq_len(n_boot)) {
  idx <- unlist(lapply(strata_indices, function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
  for (m in target_modules) {
    bootstrap_rows[[cursor]] <- module_metrics(
      dat[idx, module_genes[[m]], drop = FALSE],
      design$time[idx], m, "stratified_bootstrap", b
    )
    cursor <- cursor + 1L
  }
  if (b %% 10 == 0) message("  bootstrap ", b, "/", n_boot)
}
bootstrap_metrics <- do.call(rbind, bootstrap_rows)

reference_time_sign <- setNames(sign(reference_metrics$time_correlation), reference_metrics$module_color)
bootstrap_summary <- do.call(rbind, lapply(split(bootstrap_metrics, bootstrap_metrics$module_color), function(x) {
  m <- x$module_color[[1]]
  summarize <- function(v) c(
    median = median(v),
    ci_low = unname(quantile(v, 0.025)),
    ci_high = unname(quantile(v, 0.975))
  )
  mp <- summarize(x$mean_pairwise_correlation)
  km <- summarize(x$median_abs_kme)
  pv <- summarize(x$pc1_variance_explained)
  tc <- summarize(x$time_correlation)
  data.frame(
    module_color = m,
    bootstrap_replicates = nrow(x),
    mean_pairwise_cor_median = mp[[1]],
    mean_pairwise_cor_ci_low = mp[[2]],
    mean_pairwise_cor_ci_high = mp[[3]],
    median_abs_kme_median = km[[1]],
    median_abs_kme_ci_low = km[[2]],
    median_abs_kme_ci_high = km[[3]],
    pc1_variance_median = pv[[1]],
    pc1_variance_ci_low = pv[[2]],
    pc1_variance_ci_high = pv[[3]],
    time_cor_median = tc[[1]],
    time_cor_ci_low = tc[[2]],
    time_cor_ci_high = tc[[3]],
    time_direction_retention = mean(sign(x$time_correlation) == reference_time_sign[[m]])
  )
}))

message("Running fixed-membership leave-one-variety-out coherence checks")
fixed_loo <- do.call(rbind, lapply(unique(design$variety), function(omitted) {
  idx <- which(design$variety != omitted)
  do.call(rbind, lapply(target_modules, function(m) {
    z <- module_metrics(dat[idx, module_genes[[m]], drop = FALSE], design$time[idx], m,
                        "leave_one_variety_fixed", omitted)
    z$omitted_variety <- omitted
    z
  }))
}))

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
message("Rebuilding four leave-one-variety-out signed networks")
recluster_maps <- list()
recluster_assignments <- list()
for (omitted in unique(design$variety)) {
  message("  recluster without ", omitted)
  idx <- which(design$variety != omitted)
  dat_sub <- dat[idx, , drop = FALSE]
  net <- blockwiseModules(
    dat_sub,
    power = 12,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    maxBlockSize = ncol(dat_sub),
    nThreads = 4,
    verbose = 2
  )
  new_colors <- labels2colors(net$colors)
  assignment <- data.frame(
    omitted_variety = omitted,
    gene_id = colnames(dat_sub),
    resampled_module = new_colors
  )
  recluster_assignments[[omitted]] <- assignment

  new_sets <- split(assignment$gene_id, assignment$resampled_module)
  map_rows <- lapply(target_modules, function(m) {
    ref <- module_genes[[m]]
    stats <- do.call(rbind, lapply(names(new_sets), function(new_m) {
      new <- new_sets[[new_m]]
      overlap <- length(intersect(ref, new))
      union_n <- length(union(ref, new))
      universe <- ncol(dat_sub)
      mat <- matrix(c(
        overlap,
        length(ref) - overlap,
        length(new) - overlap,
        universe - length(union(ref, new))
      ), nrow = 2)
      ft <- fisher.test(mat, alternative = "greater")
      data.frame(
        omitted_variety = omitted,
        reference_module = m,
        resampled_module = new_m,
        reference_size = length(ref),
        resampled_size = length(new),
        overlap = overlap,
        jaccard = overlap / union_n,
        recovery_fraction = overlap / length(ref),
        odds_ratio = unname(ft$estimate),
        fisher_p = ft$p.value
      )
    }))
    stats[which.max(stats$jaccard), , drop = FALSE]
  })
  recluster_maps[[omitted]] <- do.call(rbind, map_rows)
  rm(net)
  gc()
}

recluster_map <- do.call(rbind, recluster_maps)
recluster_map$fisher_fdr <- p.adjust(recluster_map$fisher_p, "BH")
recluster_map$passes_jaccard_0_5 <- recluster_map$jaccard >= 0.50
recluster_assignment <- do.call(rbind, recluster_assignments)

stability_summary <- do.call(rbind, lapply(split(recluster_map, recluster_map$reference_module), function(x) {
  data.frame(
    module_color = x$reference_module[[1]],
    leave_one_variety_runs = nrow(x),
    minimum_jaccard = min(x$jaccard),
    median_jaccard = median(x$jaccard),
    minimum_recovery_fraction = min(x$recovery_fraction),
    all_runs_jaccard_ge_0_5 = all(x$jaccard >= 0.50),
    all_runs_fdr_lt_0_05 = all(x$fisher_fdr < 0.05)
  )
}))

write.table(reference_metrics, file.path(output_dir, "01_reference_module_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bootstrap_metrics, file.path(output_dir, "02_bootstrap_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bootstrap_summary, file.path(output_dir, "03_bootstrap_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fixed_loo, file.path(output_dir, "04_leave_one_variety_fixed_membership.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(recluster_assignment, file.path(output_dir, "05_leave_one_variety_recluster_assignments.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(recluster_map, file.path(output_dir, "06_leave_one_variety_module_mapping.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(stability_summary, file.path(output_dir, "07_module_stability_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

summary_lines <- c(
  "# WGCNA resampling stability",
  "",
  paste0("- Stratified bootstrap replicates: ", n_boot, "."),
  "- Bootstrap strata: variety x time.",
  "- Leave-one-variety-out networks: 4, rebuilt with the original 5,000 genes and fixed formal parameters.",
  "- Stability rule: every omission has maximum mapped-module Jaccard >= 0.50 and enrichment FDR < 0.05.",
  "",
  paste(apply(stability_summary, 1, function(z) {
    paste0("- ", z[["module_color"]], ": minimum Jaccard ",
           sprintf("%.3f", as.numeric(z[["minimum_jaccard"]])),
           "; all-run criterion = ", z[["all_runs_jaccard_ge_0_5"]], ".")
  }), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "SUMMARY.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Completed WGCNA stability analysis")
