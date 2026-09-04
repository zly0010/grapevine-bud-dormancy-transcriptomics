#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 6)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 10_cross_variety_meta.R <source_project> <output_dir>")
}

source_project <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

varieties <- c("Vamu", "Vvcs", "Vvri", "Vrip")

input_file <- function(variety) {
  file.path(
    source_project,
    "results_corrected", "02_deseq2_gse124820",
    paste0(variety, "_B_no_fail"),
    paste0("DEG_", variety, "_Day10_vs_Day0_B.txt")
  )
}

files <- setNames(vapply(varieties, input_file, character(1)), varieties)
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing inputs: ", paste(missing, collapse = "; "))
}

read_effect <- function(path, variety) {
  x <- read.delim(path, check.names = FALSE)
  needed <- c("gene_id", "log2FC_raw", "lfcSE_raw")
  if (!all(needed %in% names(x))) {
    stop("Missing required columns in ", path)
  }
  x <- x[, needed]
  names(x)[2:3] <- paste0(c("yi_", "sei_"), variety)
  x
}

effects <- Reduce(
  function(x, y) merge(x, y, by = "gene_id", all = FALSE, sort = FALSE),
  Map(read_effect, files, names(files))
)

yi_cols <- paste0("yi_", varieties)
se_cols <- paste0("sei_", varieties)
valid <- apply(effects[, c(yi_cols, se_cols)], 1, function(z) {
  all(is.finite(z)) && all(z[(length(varieties) + 1):(2 * length(varieties))] > 0)
})
effects <- effects[valid, , drop = FALSE]

reml_meta <- function(y, se) {
  vi <- se^2
  k <- length(y)
  objective <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * y) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (y - mu)^2))
  }
  upper <- max(c(stats::var(y), vi, 0.01), na.rm = TRUE) * 100
  tau2 <- stats::optimize(objective, interval = c(0, upper))$minimum
  if (objective(0) <= objective(tau2) + 1e-10) tau2 <- 0

  w <- 1 / (vi + tau2)
  mu <- sum(w * y) / sum(w)
  q_hk <- sum(w * (y - mu)^2) / (k - 1)
  # Modified Hartung-Knapp avoids spuriously narrow intervals when q_hk < 1.
  se_hk <- sqrt(max(1, q_hk) / sum(w))
  t_value <- mu / se_hk
  p_value <- 2 * stats::pt(abs(t_value), df = k - 1, lower.tail = FALSE)
  crit <- stats::qt(0.975, df = k - 1)

  w_fixed <- 1 / vi
  mu_fixed <- sum(w_fixed * y) / sum(w_fixed)
  q <- sum(w_fixed * (y - mu_fixed)^2)
  i2 <- if (q > 0) max(0, (q - (k - 1)) / q) else 0

  c(
    k = k,
    estimate = mu,
    se_hk = se_hk,
    ci_low = mu - crit * se_hk,
    ci_high = mu + crit * se_hk,
    p_value = p_value,
    tau2 = tau2,
    i2 = i2,
    q = q,
    q_hk = q_hk
  )
}

message("Running four-variety REML meta-analysis for ", nrow(effects), " common genes")
full_matrix <- t(vapply(seq_len(nrow(effects)), function(i) {
  reml_meta(
    as.numeric(effects[i, yi_cols]),
    as.numeric(effects[i, se_cols])
  )
}, numeric(10)))

full <- data.frame(gene_id = effects$gene_id, full_matrix, check.names = FALSE)
full$padj <- p.adjust(full$p_value, method = "BH")
full$heterogeneity_high <- full$i2 >= 0.75

loo_frames <- vector("list", length(varieties))
for (j in seq_along(varieties)) {
  omitted <- varieties[[j]]
  keep <- setdiff(seq_along(varieties), j)
  message("Running leave-one-variety-out: omit ", omitted)
  loo_matrix <- t(vapply(seq_len(nrow(effects)), function(i) {
    reml_meta(
      as.numeric(effects[i, yi_cols[keep]]),
      as.numeric(effects[i, se_cols[keep]])
    )
  }, numeric(10)))
  loo <- data.frame(gene_id = effects$gene_id, omitted_variety = omitted, loo_matrix)
  loo$padj <- p.adjust(loo$p_value, method = "BH")
  loo_frames[[j]] <- loo
}
loo <- do.call(rbind, loo_frames)

per_variety_sign <- apply(effects[, yi_cols, drop = FALSE], 1, function(z) {
  s <- sign(as.numeric(z))
  length(unique(s[s != 0])) == 1 && !any(s == 0)
})
loo_signs <- tapply(sign(loo$estimate), loo$gene_id, function(z) {
  length(z) == length(varieties) && length(unique(z[z != 0])) == 1 && !any(z == 0)
})
loo_signs <- loo_signs[match(full$gene_id, names(loo_signs))]

full$all_varieties_same_sign <- per_variety_sign
full$all_loo_same_sign <- as.logical(loo_signs)
full$robust_meta_gene <- full$padj < 0.05 &
  full$all_varieties_same_sign &
  full$all_loo_same_sign

per_variety_long <- do.call(rbind, lapply(varieties, function(v) {
  data.frame(
    gene_id = effects$gene_id,
    variety = v,
    log2_fold_change = effects[[paste0("yi_", v)]],
    standard_error = effects[[paste0("sei_", v)]]
  )
}))

modules_path <- file.path(
  source_project, "results_corrected", "07_wgcna_fixed",
  "GSE124820", "module_assignments.txt"
)
module_summary <- NULL
if (file.exists(modules_path)) {
  modules <- read.delim(modules_path, check.names = FALSE)
  annotated <- merge(full, modules[, c("gene_id", "module_color")], by = "gene_id")
  total_n <- nrow(annotated)
  total_robust <- sum(annotated$robust_meta_gene)
  module_summary <- do.call(rbind, lapply(split(annotated, annotated$module_color), function(x) {
    robust_n <- sum(x$robust_meta_gene)
    module_n <- nrow(x)
    p_enrich <- if (total_robust > 0) {
      stats::phyper(robust_n - 1, module_n, total_n - module_n, total_robust, lower.tail = FALSE)
    } else {
      NA_real_
    }
    data.frame(
      module_color = x$module_color[[1]],
      genes_with_meta_effect = module_n,
      median_meta_log2fc = stats::median(x$estimate),
      mean_meta_log2fc = mean(x$estimate),
      same_sign_fraction = mean(x$all_varieties_same_sign),
      robust_meta_genes = robust_n,
      robust_fraction = robust_n / module_n,
      robust_enrichment_p = p_enrich
    )
  }))
  module_summary$robust_enrichment_fdr <- p.adjust(module_summary$robust_enrichment_p, "BH")
  module_summary <- module_summary[order(module_summary$robust_enrichment_fdr), ]
}

write.table(per_variety_long, file.path(output_dir, "01_per_variety_effects.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(full, file.path(output_dir, "02_gene_reml_meta.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(loo, file.path(output_dir, "03_leave_one_variety_out.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(full[full$robust_meta_gene, ], file.path(output_dir, "04_robust_meta_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
if (!is.null(module_summary)) {
  write.table(module_summary, file.path(output_dir, "05_module_meta_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

summary_lines <- c(
  "# Cross-variety effect-size synthesis",
  "",
  paste0("- Contrast: Day10 vs Day0 in GSE124820 QC-B samples."),
  paste0("- Common genes analyzed: ", nrow(full), "."),
  paste0("- REML + modified Hartung-Knapp FDR < 0.05: ", sum(full$padj < 0.05), "."),
  paste0("- Four-variety same-direction genes: ", sum(full$all_varieties_same_sign), "."),
  paste0("- Robust meta genes after all leave-one-variety-out direction checks: ", sum(full$robust_meta_gene), "."),
  paste0("- High heterogeneity (I2 >= 75%): ", sum(full$heterogeneity_high), "."),
  "",
  "Robust does not mean universally causal. It means the estimated Day10 response is statistically supported, directionally concordant across the four varieties, and does not reverse when any one variety is omitted."
)
writeLines(summary_lines, file.path(output_dir, "SUMMARY.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Completed: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
