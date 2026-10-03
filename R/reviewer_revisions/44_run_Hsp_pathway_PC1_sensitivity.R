############################################################
## 44_run_Hsp_pathway_PC1_sensitivity.R
##
## Reviewer sensitivity:
## Compare the canonical Hsp60/10 RNA mean-z pathway score
## with a PC1 / eigengene summary of the exact same
## standardized client-expression matrix.
############################################################

options(stringsAsFactors = FALSE)

## =========================================================
## 1. Load canonical production objects
## =========================================================

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R",
  "R/04_build_pathway_sets_all_clients.R"
)) {
  if (!file.exists(f)) {
    stop("Missing required production script: ", f)
  }
  source(f)
}

require_objects(
  c(
    "rna_mat",
    "rna_meta_adj",
    "rna_scores",
    "all_hsp60_10_clients"
  ),
  context = "44_run_Hsp_pathway_PC1_sensitivity.R"
)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "Hsp_pathway_PC1_sensitivity"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## =========================================================
## 2. Canonical RNA/client-set validation
## =========================================================

if (ncol(rna_mat) != 577) {
  stop(
    "Expected canonical RNA cohort n=577; found ",
    ncol(rna_mat)
  )
}

if (
  !identical(
    colnames(rna_mat),
    as.character(rna_meta_adj$sample_id)
  )
) {
  stop(
    "rna_mat columns are not exactly aligned to ",
    "rna_meta_adj$sample_id."
  )
}

sample_ids <- colnames(rna_mat)

genes_use <- intersect(
  clean_gene_symbols(all_hsp60_10_clients),
  rownames(rna_mat)
)

if (length(genes_use) != 297) {
  stop(
    "Expected 297 detected Hsp60/10 RNA clients; found ",
    length(genes_use)
  )
}

if (anyDuplicated(genes_use)) {
  stop("Duplicate genes in detected Hsp60/10 RNA client set.")
}

cat(
  "Canonical RNA samples:",
  length(sample_ids),
  "\n"
)

cat(
  "Detected Hsp60/10 RNA clients:",
  length(genes_use),
  "\n"
)


## =========================================================
## 3. Build exact standardized client matrix
## =========================================================
##
## This is the same transformation used by:
## score_pathway_mean_z()
##
## rows    = genes
## columns = participants

client_z <- zscore_rows(
  rna_mat[
    genes_use,
    sample_ids,
    drop = FALSE
  ]
)

if (
  nrow(client_z) != 297 ||
  ncol(client_z) != 577
) {
  stop(
    "Unexpected standardized client-matrix dimensions: ",
    nrow(client_z),
    " x ",
    ncol(client_z)
  )
}

nonfinite_by_gene <- rowSums(
  !is.finite(client_z)
)

if (any(nonfinite_by_gene > 0)) {
  bad_genes <- names(
    nonfinite_by_gene[
      nonfinite_by_gene > 0
    ]
  )

  stop(
    "Non-finite standardized values prevent exact PCA. Genes: ",
    paste(bad_genes, collapse = ", ")
  )
}


## =========================================================
## 4. Reproduce production mean-z score exactly
## =========================================================

mean_z_recomputed <- colMeans(
  client_z
)

production_index <- match(
  sample_ids,
  rna_scores$sample_id
)

if (anyNA(production_index)) {
  stop(
    "Could not align all RNA samples to production rna_scores."
  )
}

production_mean_z <- rna_scores$Hsp60_10_all_clients[
  production_index
]

if (any(!is.finite(production_mean_z))) {
  stop(
    "Production Hsp60/10 mean-z score contains ",
    "non-finite values."
  )
}

mean_z_max_abs_diff <- max(
  abs(
    mean_z_recomputed -
      production_mean_z
  )
)

cat(
  "Maximum absolute difference vs production mean-z:",
  format(
    mean_z_max_abs_diff,
    scientific = TRUE
  ),
  "\n"
)

if (mean_z_max_abs_diff > 1e-12) {
  stop(
    "Recomputed mean-z score does not exactly reproduce ",
    "the production pathway score."
  )
}


## =========================================================
## 5. PCA / eigengene
## =========================================================
##
## Participants are observations.
## Standardized genes are variables.
##
## Because genes are already centered and standardized,
## prcomp() must not center/scale them a second time.

pca_fit <- prcomp(
  t(client_z),
  center = FALSE,
  scale. = FALSE,
  retx = TRUE
)

if (
  nrow(pca_fit$x) != 577 ||
  nrow(pca_fit$rotation) != 297
) {
  stop("Unexpected PCA output dimensions.")
}

pc1_raw <- as.numeric(
  pca_fit$x[, 1]
)

pc1_loading_raw <- as.numeric(
  pca_fit$rotation[, 1]
)

names(pc1_loading_raw) <- rownames(
  pca_fit$rotation
)

initial_cor <- cor(
  mean_z_recomputed,
  pc1_raw,
  method = "pearson"
)

## PCA sign is arbitrary.
## Orient PC1 so higher PC1 means higher mean Hsp60/10
## client abundance.

orientation_multiplier <- ifelse(
  initial_cor < 0,
  -1,
  1
)

pc1_oriented_raw <- (
  pc1_raw *
    orientation_multiplier
)

pc1_loading <- (
  pc1_loading_raw *
    orientation_multiplier
)

## Standardize participant PC1 scores for interpretable
## participant-level comparison with the pathway score.
pc1_eigengene <- as.numeric(
  scale(pc1_oriented_raw)
)

names(pc1_eigengene) <- sample_ids


## =========================================================
## 6. PCA variance explained and concordance
## =========================================================

variance_explained <- (
  pca_fit$sdev^2 /
    sum(pca_fit$sdev^2)
)

pca_variance_tbl <- data.frame(
  PC = paste0(
    "PC",
    seq_along(variance_explained)
  ),
  variance_explained = variance_explained,
  cumulative_variance_explained = cumsum(
    variance_explained
  ),
  stringsAsFactors = FALSE
)

pearson_r <- cor(
  mean_z_recomputed,
  pc1_eigengene,
  method = "pearson"
)

spearman_rho <- cor(
  mean_z_recomputed,
  pc1_eigengene,
  method = "spearman"
)

concordance_tbl <- data.frame(
  metric = c(
    "pearson_r_mean_z_vs_PC1",
    "spearman_rho_mean_z_vs_PC1",
    "PC1_variance_explained",
    "PC1_orientation_multiplier",
    "mean_z_max_abs_diff_vs_production"
  ),
  value = c(
    pearson_r,
    spearman_rho,
    variance_explained[1],
    orientation_multiplier,
    mean_z_max_abs_diff
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("MEAN-Z / PC1 CONCORDANCE\n")
cat("============================================================\n")

cat(
  "Pearson r:",
  pearson_r,
  "\n"
)

cat(
  "Spearman rho:",
  spearman_rho,
  "\n"
)

cat(
  "PC1 variance explained:",
  variance_explained[1],
  "\n"
)


## =========================================================
## 7. Gene loadings
## =========================================================

loading_tbl <- data.frame(
  gene = names(pc1_loading),
  loading = as.numeric(pc1_loading),
  abs_loading = abs(
    as.numeric(pc1_loading)
  ),
  loading_sign = ifelse(
    pc1_loading > 0,
    "positive",
    ifelse(
      pc1_loading < 0,
      "negative",
      "zero"
    )
  ),
  stringsAsFactors = FALSE
)

loading_tbl <- loading_tbl[
  order(
    -loading_tbl$abs_loading
  ),
  ,
  drop = FALSE
]

loading_summary_tbl <- data.frame(
  n_genes = nrow(loading_tbl),
  n_positive = sum(
    loading_tbl$loading > 0
  ),
  n_negative = sum(
    loading_tbl$loading < 0
  ),
  n_zero = sum(
    loading_tbl$loading == 0
  ),
  fraction_positive = mean(
    loading_tbl$loading > 0
  ),
  mean_loading = mean(
    loading_tbl$loading
  ),
  median_loading = median(
    loading_tbl$loading
  ),
  mean_abs_loading = mean(
    loading_tbl$abs_loading
  ),
  stringsAsFactors = FALSE
)


## =========================================================
## 8. Participant-level score table
## =========================================================

stage_index <- match(
  sample_ids,
  rna_scores$sample_id
)

clinical_stage <- as.character(
  rna_scores$clinical_stage[
    stage_index
  ]
)

if (anyNA(clinical_stage)) {
  stop(
    "Missing canonical clinical stage in RNA cohort."
  )
}

expected_stage_counts <- c(
  NCI = 200L,
  MCI = 158L,
  AD = 219L
)

observed_stage_counts <- table(
  factor(
    clinical_stage,
    levels = names(
      expected_stage_counts
    )
  )
)

if (
  !identical(
    as.integer(observed_stage_counts),
    as.integer(expected_stage_counts)
  )
) {
  stop(
    "Unexpected RNA stage counts. Observed: ",
    paste(
      names(observed_stage_counts),
      observed_stage_counts,
      collapse = "; "
    )
  )
}

score_tbl <- data.frame(
  sample_id = sample_ids,
  clinical_stage = factor(
    clinical_stage,
    levels = c(
      "NCI",
      "MCI",
      "AD"
    )
  ),
  mean_z = mean_z_recomputed,
  PC1_eigengene = pc1_eigengene,
  PC1_raw_oriented = pc1_oriented_raw,
  stringsAsFactors = FALSE
)


## =========================================================
## 9. Stage summaries
## =========================================================

scores_long <- dplyr::bind_rows(
  score_tbl |>
    dplyr::transmute(
      sample_id,
      clinical_stage,
      score_type = "mean_z",
      score = mean_z
    ),
  score_tbl |>
    dplyr::transmute(
      sample_id,
      clinical_stage,
      score_type = "PC1_eigengene",
      score = PC1_eigengene
    )
)

stage_summary_tbl <- scores_long |>
  dplyr::group_by(
    score_type,
    clinical_stage
  ) |>
  dplyr::summarise(
    n = sum(
      is.finite(score)
    ),
    mean_score = mean(
      score,
      na.rm = TRUE
    ),
    sd_score = sd(
      score,
      na.rm = TRUE
    ),
    sem = sd_score / sqrt(n),
    median_score = median(
      score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## =========================================================
## 10. Exact Figure-1 stage contrasts
## =========================================================

comparison_tbl <- data.frame(
  comparison = c(
    "NCI_vs_MCI",
    "MCI_vs_AD"
  ),
  group1 = c(
    "NCI",
    "MCI"
  ),
  group2 = c(
    "MCI",
    "AD"
  ),
  stringsAsFactors = FALSE
)

run_stage_contrasts <- function(
  dat,
  score_name
) {
  out <- lapply(
    seq_len(
      nrow(comparison_tbl)
    ),
    function(i) {
      g1 <- comparison_tbl$group1[i]
      g2 <- comparison_tbl$group2[i]

      x <- dat$score[
        as.character(
          dat$clinical_stage
        ) == g1
      ]

      y <- dat$score[
        as.character(
          dat$clinical_stage
        ) == g2
      ]

      p <- suppressWarnings(
        wilcox.test(
          x,
          y,
          exact = FALSE
        )$p.value
      )

      data.frame(
        score_type = score_name,
        comparison = comparison_tbl$comparison[i],
        group1 = g1,
        group2 = g2,
        n1 = sum(
          is.finite(x)
        ),
        n2 = sum(
          is.finite(y)
        ),
        mean1 = mean(
          x,
          na.rm = TRUE
        ),
        mean2 = mean(
          y,
          na.rm = TRUE
        ),
        median1 = median(
          x,
          na.rm = TRUE
        ),
        median2 = median(
          y,
          na.rm = TRUE
        ),
        delta_mean = (
          mean(
            y,
            na.rm = TRUE
          ) -
            mean(
              x,
              na.rm = TRUE
            )
        ),
        delta_median = (
          median(
            y,
            na.rm = TRUE
          ) -
            median(
              x,
              na.rm = TRUE
            )
        ),
        p_value = p,
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(out) |>
    dplyr::mutate(
      p_fdr_bh_within_score = p.adjust(
        p_value,
        method = "BH"
      )
    )
}

mean_z_stage_tests <- run_stage_contrasts(
  scores_long |>
    dplyr::filter(
      score_type == "mean_z"
    ),
  "mean_z"
)

pc1_stage_tests <- run_stage_contrasts(
  scores_long |>
    dplyr::filter(
      score_type == "PC1_eigengene"
    ),
  "PC1_eigengene"
)

stage_test_tbl <- dplyr::bind_rows(
  mean_z_stage_tests,
  pc1_stage_tests
)


## =========================================================
## 11. Direction-concordance audit
## =========================================================

direction_tbl <- stage_test_tbl |>
  dplyr::mutate(
    direction = dplyr::case_when(
      delta_mean > 0 ~ "increase",
      delta_mean < 0 ~ "decrease",
      TRUE ~ "no_change"
    )
  ) |>
  dplyr::select(
    score_type,
    comparison,
    delta_mean,
    delta_median,
    direction,
    p_value,
    p_fdr_bh_within_score
  )

direction_wide <- direction_tbl |>
  dplyr::select(
    score_type,
    comparison,
    direction
  ) |>
  tidyr::pivot_wider(
    names_from = score_type,
    values_from = direction
  )

if (
  any(
    direction_wide$mean_z !=
      direction_wide$PC1_eigengene
  )
) {
  warning(
    "Mean-z and PC1 stage-effect directions are not fully concordant."
  )
}


## =========================================================
## 12. Write outputs
## =========================================================

write.csv(
  data.frame(
    gene = genes_use,
    stringsAsFactors = FALSE
  ),
  file.path(
    OUTDIR,
    "Hsp60_10_RNA_client_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  score_tbl,
  file.path(
    OUTDIR,
    "participant_mean_z_and_PC1_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  concordance_tbl,
  file.path(
    OUTDIR,
    "mean_z_vs_PC1_concordance.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_variance_tbl,
  file.path(
    OUTDIR,
    "PCA_variance_explained.csv"
  ),
  row.names = FALSE
)

write.csv(
  loading_tbl,
  file.path(
    OUTDIR,
    "PC1_gene_loadings.csv"
  ),
  row.names = FALSE
)

write.csv(
  loading_summary_tbl,
  file.path(
    OUTDIR,
    "PC1_loading_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_summary_tbl,
  file.path(
    OUTDIR,
    "stage_score_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_test_tbl,
  file.path(
    OUTDIR,
    "stage_Wilcoxon_mean_z_vs_PC1.csv"
  ),
  row.names = FALSE
)

write.csv(
  direction_tbl,
  file.path(
    OUTDIR,
    "stage_direction_concordance.csv"
  ),
  row.names = FALSE
)


## =========================================================
## 13. Provenance
## =========================================================

git_head <- tryCatch(
  system2(
    "git",
    c(
      "rev-parse",
      "HEAD"
    ),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) NA_character_
)

provenance <- data.frame(
  item = c(
    "RNA_expression_object",
    "RNA_samples",
    "RNA_detected_Hsp60_10_clients",
    "pathway_score_method_primary",
    "sensitivity_method",
    "PCA_input",
    "PCA_center",
    "PCA_scale",
    "stage_levels",
    "stage_counts",
    "stage_test",
    "stage_test_primary_p",
    "stage_test_supplemental_adjustment",
    "git_HEAD"
  ),
  value = c(
    "rna_mat",
    ncol(rna_mat),
    length(genes_use),
    "mean of gene-wise z-scores",
    "PC1/eigengene",
    "same gene-wise z-scored client matrix",
    "FALSE",
    "FALSE",
    "NCI;MCI;AD",
    paste(
      names(expected_stage_counts),
      expected_stage_counts,
      sep = "=",
      collapse = ";"
    ),
    "Wilcoxon rank-sum; exact=FALSE",
    "raw p-value to reproduce Figure 1",
    "BH across 2 contrasts within score type",
    paste(
      git_head,
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(
    OUTDIR,
    "analysis_provenance.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)


## =========================================================
## 14. Final hard validation
## =========================================================

if (
  nrow(score_tbl) != 577 ||
  anyDuplicated(score_tbl$sample_id)
) {
  stop(
    "Participant score-table validation failed."
  )
}

if (
  any(
    !is.finite(
      score_tbl$mean_z
    )
  ) ||
  any(
    !is.finite(
      score_tbl$PC1_eigengene
    )
  )
) {
  stop(
    "Non-finite participant scores remain."
  )
}

if (
  abs(
    mean(
      score_tbl$PC1_eigengene
    )
  ) > 1e-10
) {
  stop(
    "Standardized PC1 eigengene is not centered at zero."
  )
}

if (
  abs(
    sd(
      score_tbl$PC1_eigengene
    ) -
      1
  ) > 1e-10
) {
  stop(
    "Standardized PC1 eigengene does not have SD=1."
  )
}

if (pearson_r <= 0) {
  stop(
    "PC1 orientation failed: correlation with mean-z is non-positive."
  )
}

if (nrow(loading_tbl) != 297) {
  stop(
    "Expected 297 PC1 gene loadings; found ",
    nrow(loading_tbl)
  )
}

if (nrow(stage_test_tbl) != 4) {
  stop(
    "Expected 4 total stage tests ",
    "(2 contrasts x 2 score types)."
  )
}

cat("\n============================================================\n")
cat("Hsp60/10 PATHWAY PC1 SENSITIVITY COMPLETED AND VALIDATED\n")
cat("============================================================\n")

cat(
  "RNA samples:",
  nrow(score_tbl),
  "\n"
)

cat(
  "Detected clients:",
  length(genes_use),
  "\n"
)

cat(
  "Mean-z reproduction max abs diff:",
  mean_z_max_abs_diff,
  "\n"
)

cat(
  "Mean-z vs PC1 Pearson r:",
  pearson_r,
  "\n"
)

cat(
  "Mean-z vs PC1 Spearman rho:",
  spearman_rho,
  "\n"
)

cat(
  "PC1 variance explained:",
  variance_explained[1],
  "\n\n"
)

print(
  stage_summary_tbl
)

cat("\nStage contrasts:\n")

print(
  stage_test_tbl
)

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)
