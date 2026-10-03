options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  source(f)
}

INDIR <- "outputs/reviewer_revisions/variancePartition_source_of_variation"
OUTDIR <- "outputs/reviewer_revisions/variancePartition_postmortem"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("======================================================================\n")
cat("VARIANCEPARTITION POSTMORTEM\n")
cat("======================================================================\n")


# =============================================================================
# 1. Confirm exactly what happened to the 1,998 protein failures
# =============================================================================

protein_manifest <- read.csv(
  file.path(
    INDIR,
    "protein_variancePartition_feature_manifest.csv"
  ),
  check.names = FALSE
)

protein_vp <- readRDS(
  file.path(
    INDIR,
    "protein_variancePartition_results.rds"
  )
)

returned_ids <- rownames(
  protein_vp
)

protein_manifest$returned_by_variancePartition <- (
  protein_manifest$feature_id %in%
    returned_ids
)

cat("\n======================================================================\n")
cat("PROTEIN FAILURE STRUCTURE\n")
cat("======================================================================\n")

cat(
  "Intended >=80% proteins:",
  sum(
    protein_manifest$included_primary_80pct
  ),
  "\n"
)

cat(
  "Returned protein models:",
  length(returned_ids),
  "\n"
)

cat(
  "Failed / omitted:",
  sum(
    protein_manifest$included_primary_80pct &
      !protein_manifest$returned_by_variancePartition
  ),
  "\n"
)

cat(
  "\nComplete all 398 AND returned:\n"
)

print(
  table(
    complete_all_398 =
      protein_manifest$complete_all_398,
    returned =
      protein_manifest$returned_by_variancePartition
  )
)

failed <- protein_manifest[
  protein_manifest$included_primary_80pct &
    !protein_manifest$returned_by_variancePartition,
  ,
  drop = FALSE
]

cat(
  "\nFailed proteins observed-sample distribution:\n"
)

print(
  summary(
    failed$n_observed
  )
)

cat(
  "\nReturned proteins observed-sample distribution:\n"
)

print(
  summary(
    protein_manifest$n_observed[
      protein_manifest$returned_by_variancePartition
    ]
  )
)

write.csv(
  protein_manifest,
  file.path(
    OUTDIR,
    "protein_returned_vs_missingness_manifest.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 2. Inspect variancePartition error attributes
# =============================================================================

cat("\n======================================================================\n")
cat("PROTEIN ERROR ATTRIBUTE\n")
cat("======================================================================\n")

errors <- attr(
  protein_vp,
  "errors"
)

if (is.null(errors)) {

  cat("No 'errors' attribute present on saved object.\n")

} else {

  print(
    head(
      errors,
      30
    )
  )

  capture.output(
    print(errors),
    file = file.path(
      OUTDIR,
      "protein_errors_full.txt"
    )
  )
}


# =============================================================================
# 3. Fix/check sample-name issue explicitly
# =============================================================================

cat("\n======================================================================\n")
cat("SAMPLE NAME / ORDER AUDIT\n")
cat("======================================================================\n")

rna_ids <- as.character(
  rna_meta_adj$sample_id
)

prot_ids <- as.character(
  prot_meta_adj$sample_id
)

cat(
  "RNA column order exact:",
  identical(
    colnames(rna_mat_raw),
    rna_ids
  ),
  "\n"
)

cat(
  "Protein column order exact:",
  identical(
    colnames(prot_mat_raw),
    prot_ids
  ),
  "\n"
)

cat(
  "RNA metadata rownames currently equal sample IDs:",
  identical(
    rownames(rna_meta_adj),
    rna_ids
  ),
  "\n"
)

cat(
  "Protein metadata rownames currently equal sample IDs:",
  identical(
    rownames(prot_meta_adj),
    prot_ids
  ),
  "\n"
)

cat(
  "\nThe variancePartition warning is therefore ",
  "a labeling issue if column order is exact, but the final script ",
  "should explicitly set metadata rownames to sample IDs.\n",
  sep = ""
)


# =============================================================================
# 4. RNA batch x RIN structure
# =============================================================================

cat("\n======================================================================\n")
cat("RNA BATCH x RIN STRUCTURE\n")
cat("======================================================================\n")

rna_d <- data.frame(
  sample_id =
    as.character(
      rna_meta_adj$sample_id
    ),

  batch =
    factor(
      as.character(
        rna_meta$library_batch
      )
    ),

  rin =
    as.numeric(
      rna_meta_adj$rin_num
    ),

  age =
    as.numeric(
      rna_meta_adj$age_num
    ),

  pmi =
    as.numeric(
      rna_meta_adj$pmi_num
    ),

  sex =
    as.character(
      rna_meta_adj$sex_label
    ),

  stage =
    as.character(
      rna_meta_adj$clinical_stage
    ),

  braak =
    as.numeric(
      rna_meta_adj$braak_num
    ),

  cerad =
    as.numeric(
      rna_meta_adj$cerad_num
    ),

  stringsAsFactors = FALSE
)


batch_summary <- do.call(
  rbind,
  lapply(
    split(
      rna_d,
      rna_d$batch
    ),
    function(x) {

      data.frame(
        batch =
          as.character(
            x$batch[1]
          ),

        n =
          nrow(x),

        rin_mean =
          mean(
            x$rin,
            na.rm = TRUE
          ),

        rin_sd =
          sd(
            x$rin,
            na.rm = TRUE
          ),

        rin_min =
          min(
            x$rin,
            na.rm = TRUE
          ),

        rin_max =
          max(
            x$rin,
            na.rm = TRUE
          ),

        age_mean =
          mean(
            x$age,
            na.rm = TRUE
          ),

        pmi_mean =
          mean(
            x$pmi,
            na.rm = TRUE
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(batch_summary) <- NULL

cat("\nRIN by RNA batch:\n")

print(
  batch_summary,
  row.names = FALSE
)

write.csv(
  batch_summary,
  file.path(
    OUTDIR,
    "RNA_batch_covariate_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 5. RNA batch x disease stage
# =============================================================================

cat("\n======================================================================\n")
cat("RNA BATCH x CLINICAL STAGE\n")
cat("======================================================================\n")

stage_tab <- table(
  batch = rna_d$batch,
  stage = rna_d$stage,
  useNA = "ifany"
)

print(
  stage_tab
)

stage_prop <- prop.table(
  stage_tab,
  margin = 1
)

cat("\nRow proportions:\n")

print(
  round(
    stage_prop,
    3
  )
)

write.csv(
  as.data.frame.matrix(
    stage_tab
  ),
  file.path(
    OUTDIR,
    "RNA_batch_by_stage_counts.csv"
  )
)

write.csv(
  as.data.frame.matrix(
    stage_prop
  ),
  file.path(
    OUTDIR,
    "RNA_batch_by_stage_row_proportions.csv"
  )
)


# =============================================================================
# 6. Quantify how strongly batch predicts RIN
# =============================================================================

cat("\n======================================================================\n")
cat("HOW MUCH RIN VARIANCE IS EXPLAINED BY BATCH?\n")
cat("======================================================================\n")

m_rin <- lm(
  rin ~ batch,
  data = rna_d
)

rin_batch_r2 <- summary(
  m_rin
)$r.squared

rin_batch_adj_r2 <- summary(
  m_rin
)$adj.r.squared

cat(
  "RIN ~ batch R2:",
  rin_batch_r2,
  "\n"
)

cat(
  "RIN ~ batch adjusted R2:",
  rin_batch_adj_r2,
  "\n"
)


# =============================================================================
# 7. Does batch remain in the CURRENT adjusted RNA matrix?
# =============================================================================

cat("\n======================================================================\n")
cat("BATCH SIGNAL IN CURRENT ADJUSTED RNA MATRIX\n")
cat("======================================================================\n")

# Current rna_mat_adj was adjusted using production covariates.
# Test per-gene batch R2 descriptively using ordinary ANOVA.
#
# This is NOT the final biological analysis; it asks whether the current
# residualization left substantial batch-associated variation.

batch_r2_one_gene <- function(y, batch) {

  fit <- lm(
    y ~ batch
  )

  summary(fit)$r.squared
}

set.seed(1300)

audit_genes <- unique(
  round(
    seq(
      1,
      nrow(rna_mat_adj),
      length.out = min(
        2000,
        nrow(rna_mat_adj)
      )
    )
  )
)

batch_r2_raw <- vapply(
  audit_genes,
  function(i) {
    batch_r2_one_gene(
      as.numeric(
        rna_mat_raw[i, ]
      ),
      rna_d$batch
    )
  },
  numeric(1)
)

batch_r2_adjusted <- vapply(
  audit_genes,
  function(i) {
    batch_r2_one_gene(
      as.numeric(
        rna_mat_adj[i, ]
      ),
      rna_d$batch
    )
  },
  numeric(1)
)

batch_residual_audit <- data.frame(
  feature_id =
    rownames(rna_mat_raw)[
      audit_genes
    ],

  raw_batch_R2 =
    batch_r2_raw,

  adjusted_batch_R2 =
    batch_r2_adjusted,

  stringsAsFactors = FALSE
)

write.csv(
  batch_residual_audit,
  file.path(
    OUTDIR,
    "RNA_batch_R2_raw_vs_current_adjusted_2000gene_audit.csv"
  ),
  row.names = FALSE
)

cat(
  "\nRaw RNA batch R2:\n"
)

print(
  summary(
    batch_r2_raw
  )
)

cat(
  "\nCurrent adjusted RNA batch R2:\n"
)

print(
  summary(
    batch_r2_adjusted
  )
)

cat(
  "\nSpearman raw vs adjusted batch R2:",
  cor(
    batch_r2_raw,
    batch_r2_adjusted,
    method = "spearman"
  ),
  "\n"
)


# =============================================================================
# 8. Current production adjustment reminder
# =============================================================================

cat("\n======================================================================\n")
cat("CURRENT PRODUCTION RNA ADJUSTMENT\n")
cat("======================================================================\n")

cat(
  "rna_covars: ",
  paste(
    rna_covars,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "\nRNA batch is NOT currently listed in rna_covars.\n"
)


# =============================================================================
# 9. Save compact decision table
# =============================================================================

decision <- data.frame(
  metric = c(
    "protein_intended_80pct",
    "protein_returned",
    "protein_failed_or_omitted",
    "protein_complete_all_398",
    "RNA_RIN_explained_by_batch_R2",
    "RNA_raw_batch_R2_median_2000genes",
    "RNA_adjusted_batch_R2_median_2000genes"
  ),

  value = c(
    sum(
      protein_manifest$included_primary_80pct
    ),

    length(returned_ids),

    sum(
      protein_manifest$included_primary_80pct &
        !protein_manifest$returned_by_variancePartition
    ),

    sum(
      protein_manifest$complete_all_398
    ),

    rin_batch_r2,

    median(
      batch_r2_raw,
      na.rm = TRUE
    ),

    median(
      batch_r2_adjusted,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  decision,
  file.path(
    OUTDIR,
    "variancePartition_postmortem_decision_metrics.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("DECISION METRICS\n")
cat("======================================================================\n")

print(
  decision,
  row.names = FALSE
)

cat(
  "\nPOSTMORTEM COMPLETED.\n"
)
