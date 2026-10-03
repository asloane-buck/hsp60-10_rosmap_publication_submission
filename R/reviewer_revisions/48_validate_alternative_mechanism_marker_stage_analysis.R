############################################################
## 48_validate_alternative_mechanism_marker_stage_analysis.R
##
## Hard validation of Reviewer 3.1 alternative-mechanism
## marker-panel analysis.
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "alternative_mechanism_marker_panel"
)

required_files <- c(
  "prespecified_marker_detectability_manifest.csv",
  "marker_stage_limma_results.csv",
  "marker_stage_observation_counts.csv",
  "marker_category_stage_summary.csv",
  "undetected_prespecified_markers.csv",
  "marker_model_estimability.csv",
  "marker_model_estimability_summary.csv",
  "non_estimable_marker_contrasts.csv",
  "marker_model_warning_audit.csv",
  "marker_partial_na_coefficient_details.csv",
  "marker_partial_na_gene_audit.csv",
  "marker_stage_analysis_provenance.csv",
  "marker_stage_sessionInfo.txt"
)

missing_files <- required_files[
  !file.exists(file.path(OUTDIR, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing required outputs: ",
    paste(missing_files, collapse = ", ")
  )
}

manifest <- read.csv(
  file.path(
    OUTDIR,
    "prespecified_marker_detectability_manifest.csv"
  ),
  stringsAsFactors = FALSE
)

res <- read.csv(
  file.path(
    OUTDIR,
    "marker_stage_limma_results.csv"
  ),
  stringsAsFactors = FALSE
)

counts <- read.csv(
  file.path(
    OUTDIR,
    "marker_stage_observation_counts.csv"
  ),
  stringsAsFactors = FALSE
)

undetected <- read.csv(
  file.path(
    OUTDIR,
    "undetected_prespecified_markers.csv"
  ),
  stringsAsFactors = FALSE
)

estimability <- read.csv(
  file.path(
    OUTDIR,
    "marker_model_estimability.csv"
  ),
  stringsAsFactors = FALSE
)

nonestimable <- read.csv(
  file.path(
    OUTDIR,
    "non_estimable_marker_contrasts.csv"
  ),
  stringsAsFactors = FALSE
)

warning_audit <- read.csv(
  file.path(
    OUTDIR,
    "marker_model_warning_audit.csv"
  ),
  stringsAsFactors = FALSE
)

partial_coef <- read.csv(
  file.path(
    OUTDIR,
    "marker_partial_na_coefficient_details.csv"
  ),
  stringsAsFactors = FALSE
)

partial_gene <- read.csv(
  file.path(
    OUTDIR,
    "marker_partial_na_gene_audit.csv"
  ),
  stringsAsFactors = FALSE
)

provenance <- read.csv(
  file.path(
    OUTDIR,
    "marker_stage_analysis_provenance.csv"
  ),
  stringsAsFactors = FALSE
)

checks <- data.frame(
  check = character(),
  passed = logical(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_check <- function(name, passed, detail) {
  checks <<- rbind(
    checks,
    data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  )
}

## ---------------------------------------------------------
## Marker universe
## ---------------------------------------------------------

add_check(
  "prespecified_markers_44",
  nrow(manifest) == 44,
  paste("observed =", nrow(manifest))
)

detected <- manifest[
  manifest$detected_raw_TMT,
  ,
  drop = FALSE
]

add_check(
  "detected_markers_39",
  nrow(detected) == 39,
  paste("observed =", nrow(detected))
)

add_check(
  "undetected_markers_5",
  nrow(undetected) == 5,
  paste("observed =", nrow(undetected))
)

expected_undetected <- sort(
  c(
    "ESRRA",
    "PPARGC1A",
    "PPARGC1B",
    "PINK1",
    "PRKN"
  )
)

add_check(
  "undetected_marker_identity",
  identical(
    sort(undetected$gene),
    expected_undetected
  ),
  paste(
    sort(undetected$gene),
    collapse = ";"
  )
)

expected_detected_categories <- c(
  autophagy_lysosome = 8L,
  mitochondrial_biogenesis = 4L,
  mitochondrial_mass_import = 7L,
  mitophagy = 6L,
  proteasome_20S_core = 14L
)

observed_category_counts <- table(
  detected$category
)

add_check(
  "detected_category_counts",
  identical(
    as.integer(
      observed_category_counts[
        names(expected_detected_categories)
      ]
    ),
    as.integer(expected_detected_categories)
  ),
  paste(
    names(observed_category_counts),
    observed_category_counts,
    collapse = "; "
  )
)

## ---------------------------------------------------------
## Result dimensions / contrasts
## ---------------------------------------------------------

add_check(
  "result_rows_117",
  nrow(res) == 39 * 3,
  paste("observed =", nrow(res))
)

expected_contrasts <- c(
  "MCI_vs_NCI",
  "AD_vs_NCI",
  "AD_vs_MCI"
)

add_check(
  "three_expected_contrasts",
  identical(
    sort(unique(res$contrast)),
    sort(expected_contrasts)
  ),
  paste(
    sort(unique(res$contrast)),
    collapse = ";"
  )
)

results_per_contrast <- table(res$contrast)

add_check(
  "39_results_per_contrast",
  all(
    results_per_contrast[
      expected_contrasts
    ] == 39
  ),
  paste(
    names(results_per_contrast),
    results_per_contrast,
    collapse = "; "
  )
)

## ---------------------------------------------------------
## Finite statistical results
## ---------------------------------------------------------

numeric_result_cols <- c(
  "effect",
  "std_error",
  "conf_low",
  "conf_high",
  "statistic",
  "p_value",
  "marker_panel_fdr"
)

finite_ok <- all(
  vapply(
    res[numeric_result_cols],
    function(x) all(is.finite(x)),
    logical(1)
  )
)

add_check(
  "all_statistics_finite",
  finite_ok,
  "effect/SE/CI/t/P/FDR"
)

add_check(
  "valid_p_values",
  all(
    res$p_value >= 0 &
      res$p_value <= 1
  ),
  paste(
    "range =",
    paste(range(res$p_value), collapse = " to ")
  )
)

add_check(
  "valid_fdr_values",
  all(
    res$marker_panel_fdr >= 0 &
      res$marker_panel_fdr <= 1
  ),
  paste(
    "range =",
    paste(
      range(res$marker_panel_fdr),
      collapse = " to "
    )
  )
)

## ---------------------------------------------------------
## Estimability / structured missingness
## ---------------------------------------------------------

add_check(
  "no_nonestimable_stage_coefficients",
  !any(estimability$any_na_stage_coefficient),
  paste(
    "n =",
    sum(estimability$any_na_stage_coefficient)
  )
)

add_check(
  "no_nonestimable_contrasts",
  nrow(nonestimable) == 0,
  paste("n =", nrow(nonestimable))
)

affected_genes <- sort(
  unique(partial_coef$gene)
)

expected_partial_genes <- sort(
  c(
    "TIMM23",
    "GABPA",
    "GABPB1",
    "NRF1",
    "BNIP3L",
    "CALCOCO2",
    "FUNDC1",
    "ATG5",
    "MAP1LC3B"
  )
)

add_check(
  "nine_partial_na_markers",
  identical(
    affected_genes,
    expected_partial_genes
  ),
  paste(
    affected_genes,
    collapse = ";"
  )
)

add_check(
  "partial_na_batch_only",
  all(partial_coef$is_batch_coefficient),
  paste(
    "non-batch =",
    sum(!partial_coef$is_batch_coefficient)
  )
)

add_check(
  "partial_na_counts_match_zero_batches",
  all(partial_gene$counts_match),
  paste(
    "failures =",
    sum(!partial_gene$counts_match)
  )
)

warning_ok <- (
  nrow(warning_audit) == 1 &&
  warning_audit$n_warning_messages[1] == 1 &&
  identical(
    warning_audit$warning_messages[1],
    "Partial NA coefficients for 9 probe(s)"
  )
)

add_check(
  "expected_limma_warning_only",
  warning_ok,
  paste(
    warning_audit$n_warning_messages[1],
    warning_audit$warning_messages[1]
  )
)

## ---------------------------------------------------------
## Observation-count sanity
## ---------------------------------------------------------

add_check(
  "observation_count_rows_39",
  nrow(counts) == 39,
  paste("observed =", nrow(counts))
)

add_check(
  "stage_counts_do_not_exceed_cohort",
  all(
    counts$n_NCI <= 167 &
      counts$n_MCI <= 96 &
      counts$n_AD <= 109
  ),
  "NCI<=167; MCI<=96; AD<=109"
)

add_check(
  "total_n_matches_stage_sum",
  all(
    counts$n_total ==
      counts$n_NCI +
      counts$n_MCI +
      counts$n_AD
  ),
  paste(
    "mismatches =",
    sum(
      counts$n_total !=
        counts$n_NCI +
        counts$n_MCI +
        counts$n_AD
    )
  )
)

## ---------------------------------------------------------
## FDR-significant results
## ---------------------------------------------------------

sig <- res[
  res$significant_marker_panel_fdr05,
  ,
  drop = FALSE
]

add_check(
  "nine_fdr_significant_results",
  nrow(sig) == 9,
  paste("observed =", nrow(sig))
)

add_check(
  "all_significant_are_AD_vs_NCI",
  all(sig$contrast == "AD_vs_NCI"),
  paste(
    unique(sig$contrast),
    collapse = ";"
  )
)

expected_sig_genes <- sort(
  c(
    "TOMM22",
    "TOMM40",
    "BNIP3",
    "BNIP3L",
    "ATG7",
    "PSMA2",
    "PSMB1",
    "PSMB5",
    "PSMB6"
  )
)

add_check(
  "significant_gene_identity",
  identical(
    sort(sig$gene),
    expected_sig_genes
  ),
  paste(
    sort(sig$gene),
    collapse = ";"
  )
)

expected_decreased <- sort(
  c(
    "TOMM22",
    "TOMM40",
    "BNIP3",
    "BNIP3L",
    "ATG7"
  )
)

observed_decreased <- sort(
  sig$gene[sig$effect < 0]
)

add_check(
  "expected_decreased_markers",
  identical(
    observed_decreased,
    expected_decreased
  ),
  paste(
    observed_decreased,
    collapse = ";"
  )
)

expected_increased <- sort(
  c(
    "PSMA2",
    "PSMB1",
    "PSMB5",
    "PSMB6"
  )
)

observed_increased <- sort(
  sig$gene[sig$effect > 0]
)

add_check(
  "expected_increased_markers",
  identical(
    observed_increased,
    expected_increased
  ),
  paste(
    observed_increased,
    collapse = ";"
  )
)

## ---------------------------------------------------------
## Provenance
## ---------------------------------------------------------

prov <- setNames(
  as.character(provenance$value),
  provenance$item
)

add_check(
  "provenance_model_samples_372",
  identical(
    unname(prov["model_samples"]),
    "372"
  ),
  paste("value =", prov["model_samples"])
)

add_check(
  "provenance_covariates",
  identical(
    unname(prov["protein_covariates"]),
    "age_num;sex_factor;pmi_num;batch_factor"
  ),
  paste(
    "value =",
    prov["protein_covariates"]
  )
)

add_check(
  "provenance_BH",
  identical(
    unname(prov["multiple_testing_method"]),
    "Benjamini-Hochberg"
  ),
  paste(
    "value =",
    prov["multiple_testing_method"]
  )
)

## ---------------------------------------------------------
## Report
## ---------------------------------------------------------

print(
  checks,
  row.names = FALSE
)

n_pass <- sum(checks$passed)
n_total <- nrow(checks)

cat(
  "\nValidation:",
  n_pass,
  "/",
  n_total,
  "checks passed\n"
)

write.csv(
  checks,
  file.path(
    OUTDIR,
    "marker_stage_validation_checks.csv"
  ),
  row.names = FALSE
)

if (!all(checks$passed)) {
  failed <- checks$check[
    !checks$passed
  ]

  stop(
    "R3.1 marker-panel validation failed: ",
    paste(failed, collapse = ", ")
  )
}

cat(
  "\n============================================================\n"
)
cat(
  "R3.1 ALTERNATIVE-MECHANISM MARKER PANEL VALIDATED\n"
)
cat(
  "============================================================\n"
)
