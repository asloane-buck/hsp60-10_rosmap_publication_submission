############################################################
## 51_validate_APOE_sensitivity.R
##
## Validation for Reviewer 3.6 APOE-e4 sensitivity analysis.
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "APOE_sensitivity"
)

required_files <- c(
  "APOE_participant_lookup.csv",
  "APOE_cohort_summary.csv",
  "APOE_e4_carrier_by_stage.csv",
  "APOE_sensitivity_cohort_audit.csv",
  "APOE_sensitivity_late_stage_clients.csv",
  "APOE_sensitivity_Braak_clients.csv",
  "APOE_sensitivity_CERAD_clients.csv",
  "APOE_sensitivity_production_model_validation.csv",
  "APOE_sensitivity_same_subset_paired_effects.csv",
  "APOE_sensitivity_effect_summary.csv",
  "APOE_sensitivity_vulnerability_results.csv",
  "APOE_sensitivity_vulnerability_paired.csv",
  "APOE_sensitivity_vulnerability_summary.csv",
  "APOE_sensitivity_top_quartile_overlap.csv",
  "APOE_sensitivity_full_vs_subset_summary.csv",
  "APOE_sensitivity_provenance.csv",
  "APOE_sensitivity_sessionInfo.txt"
)

missing_files <- required_files[
  !file.exists(file.path(OUTDIR, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing APOE sensitivity output(s): ",
    paste(missing_files, collapse = ", ")
  )
}

read_out <- function(name) {
  utils::read.csv(
    file.path(OUTDIR, name),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

cohort <- read_out("APOE_cohort_summary.csv")
stage <- read_out("APOE_e4_carrier_by_stage.csv")
audit <- read_out("APOE_sensitivity_cohort_audit.csv")
late <- read_out("APOE_sensitivity_late_stage_clients.csv")
braak <- read_out("APOE_sensitivity_Braak_clients.csv")
cerad <- read_out("APOE_sensitivity_CERAD_clients.csv")
prod <- read_out("APOE_sensitivity_production_model_validation.csv")
paired <- read_out("APOE_sensitivity_same_subset_paired_effects.csv")
effect_summary <- read_out("APOE_sensitivity_effect_summary.csv")
vuln <- read_out("APOE_sensitivity_vulnerability_results.csv")
vuln_paired <- read_out("APOE_sensitivity_vulnerability_paired.csv")
vuln_summary <- read_out("APOE_sensitivity_vulnerability_summary.csv")
topq <- read_out("APOE_sensitivity_top_quartile_overlap.csv")
full_subset <- read_out("APOE_sensitivity_full_vs_subset_summary.csv")
prov <- read_out("APOE_sensitivity_provenance.csv")

checks <- list()

add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

get_audit <- function(metric) {
  hit <- audit$value[audit$metric == metric]
  if (length(hit) != 1L) return(NA_character_)
  hit
}

get_prov <- function(item) {
  hit <- prov$value[prov$item == item]
  if (length(hit) != 1L) return(NA_character_)
  hit
}

models_expected <- c(
  "full_original",
  "APOE_complete_original_covariates",
  "APOE_complete_plus_e4"
)

## ---------------------------------------------------------
## Cohort bookkeeping
## ---------------------------------------------------------

add_check(
  "full_protein_n_400",
  identical(get_audit("full_protein_samples"), "400"),
  paste0("observed = ", get_audit("full_protein_samples"))
)

add_check(
  "APOE_complete_n_338",
  identical(get_audit("APOE_complete_protein_samples"), "338"),
  paste0("observed = ", get_audit("APOE_complete_protein_samples"))
)

add_check(
  "APOE_carrier_counts_271_67",
  identical(get_audit("APOE_e4_noncarriers"), "271") &&
    identical(get_audit("APOE_e4_carriers"), "67"),
  paste0(
    "noncarrier=", get_audit("APOE_e4_noncarriers"),
    "; carrier=", get_audit("APOE_e4_carriers")
  )
)

add_check(
  "APOE_stage_nuisance_complete_317",
  identical(
    get_audit("APOE_complete_stage_nuisance_complete"),
    "NCI=148;MCI=86;AD=83"
  ),
  paste0(
    "observed = ",
    get_audit("APOE_complete_stage_nuisance_complete")
  )
)

add_check(
  "Hsp_client_n_306",
  identical(get_audit("Hsp60_10_client_proteins"), "306"),
  paste0("observed = ", get_audit("Hsp60_10_client_proteins"))
)

protein_cohort_row <- cohort[
  cohort$cohort == "Protein_full_400",
  ,
  drop = FALSE
]
rna_cohort_row <- cohort[
  cohort$cohort == "RNA_corrected_577",
  ,
  drop = FALSE
]

add_check(
  "protein_APOE_availability_338_of_400",
  nrow(protein_cohort_row) == 1L &&
    protein_cohort_row$n_total == 400L &&
    protein_cohort_row$n_apoe_available == 338L &&
    protein_cohort_row$n_apoe_missing == 62L,
  if (nrow(protein_cohort_row) == 1L) {
    paste0(
      protein_cohort_row$n_apoe_available,
      "/", protein_cohort_row$n_total,
      "; missing=", protein_cohort_row$n_apoe_missing
    )
  } else {
    "Protein_full_400 row missing/duplicated"
  }
)

add_check(
  "RNA_APOE_availability_576_of_577",
  nrow(rna_cohort_row) == 1L &&
    rna_cohort_row$n_total == 577L &&
    rna_cohort_row$n_apoe_available == 576L &&
    rna_cohort_row$n_apoe_missing == 1L,
  if (nrow(rna_cohort_row) == 1L) {
    paste0(
      rna_cohort_row$n_apoe_available,
      "/", rna_cohort_row$n_total,
      "; missing=", rna_cohort_row$n_apoe_missing
    )
  } else {
    "RNA_corrected_577 row missing/duplicated"
  }
)

## ---------------------------------------------------------
## Result dimensions / model identities
## ---------------------------------------------------------

for (nm in c("late", "braak", "cerad")) {
  tbl <- get(nm)
  add_check(
    paste0(nm, "_rows_918"),
    nrow(tbl) == 918L,
    paste0("observed = ", nrow(tbl))
  )

  model_counts <- table(tbl$model)
  ok <- setequal(names(model_counts), models_expected) &&
    all(model_counts[models_expected] == 306L)

  add_check(
    paste0(nm, "_306_per_model"),
    ok,
    paste(
      paste(names(model_counts), as.integer(model_counts), sep = "="),
      collapse = "; "
    )
  )

  stat_cols <- c("effect", "se", "ci_low", "ci_high", "p_value", "fdr")
  finite_ok <- all(vapply(
    stat_cols,
    function(x) all(is.finite(tbl[[x]])),
    logical(1)
  ))

  add_check(
    paste0(nm, "_all_statistics_finite"),
    finite_ok,
    paste(stat_cols, collapse = "/")
  )

  add_check(
    paste0(nm, "_valid_p_and_fdr"),
    all(tbl$p_value >= 0 & tbl$p_value <= 1) &&
      all(tbl$fdr >= 0 & tbl$fdr <= 1),
    paste0(
      "p=", signif(min(tbl$p_value), 4), "..",
      signif(max(tbl$p_value), 4),
      "; FDR=", signif(min(tbl$fdr), 4), "..",
      signif(max(tbl$fdr), 4)
    )
  )
}

## ---------------------------------------------------------
## Exact production helper equivalence
## ---------------------------------------------------------

add_check(
  "production_validation_six_rows",
  nrow(prod) == 6L,
  paste0("observed = ", nrow(prod))
)

add_check(
  "production_validation_zero_effect_difference",
  nrow(prod) == 6L &&
    all(prod$max_effect_difference == 0),
  paste0(
    "max = ",
    if (nrow(prod) > 0) max(prod$max_effect_difference) else NA
  )
)

add_check(
  "production_validation_zero_p_difference",
  nrow(prod) == 6L &&
    all(prod$max_p_difference == 0),
  paste0(
    "max = ",
    if (nrow(prod) > 0) max(prod$max_p_difference) else NA
  )
)

add_check(
  "production_validation_zero_n_mismatch",
  nrow(prod) == 6L &&
    all(prod$n_mismatches == 0L),
  paste0(
    "total mismatches = ",
    if (nrow(prod) > 0) sum(prod$n_mismatches) else NA
  )
)

## ---------------------------------------------------------
## Same-subset comparison integrity
## ---------------------------------------------------------

add_check(
  "effect_summary_three_endpoints",
  nrow(effect_summary) == 3L &&
    setequal(
      effect_summary$endpoint,
      c("Late_AD_minus_MCI", "Braak", "CERAD")
    ),
  paste(effect_summary$endpoint, collapse = ";")
)

add_check(
  "paired_effect_rows_918",
  nrow(paired) == 918L,
  paste0("observed = ", nrow(paired))
)

add_check(
  "same_subset_all_306_shared",
  all(effect_summary$n_shared == 306L),
  paste(
    paste(effect_summary$endpoint, effect_summary$n_shared, sep = "="),
    collapse = "; "
  )
)

late_sum <- effect_summary[
  effect_summary$endpoint == "Late_AD_minus_MCI",
  ,
  drop = FALSE
]
braak_sum <- effect_summary[
  effect_summary$endpoint == "Braak",
  ,
  drop = FALSE
]
cerad_sum <- effect_summary[
  effect_summary$endpoint == "CERAD",
  ,
  drop = FALSE
]

add_check(
  "late_APOE_effect_high_concordance",
  nrow(late_sum) == 1L &&
    late_sum$pearson_effect > 0.999 &&
    late_sum$spearman_effect > 0.995 &&
    late_sum$direction_concordance == 1,
  if (nrow(late_sum) == 1L) {
    paste0(
      "Pearson=", signif(late_sum$pearson_effect, 6),
      "; Spearman=", signif(late_sum$spearman_effect, 6),
      "; direction=", signif(late_sum$direction_concordance, 6)
    )
  } else {
    "late summary missing"
  }
)

add_check(
  "Braak_APOE_effect_high_concordance",
  nrow(braak_sum) == 1L &&
    braak_sum$pearson_effect > 0.98 &&
    braak_sum$spearman_effect > 0.98 &&
    braak_sum$direction_concordance > 0.95,
  if (nrow(braak_sum) == 1L) {
    paste0(
      "Pearson=", signif(braak_sum$pearson_effect, 6),
      "; Spearman=", signif(braak_sum$spearman_effect, 6),
      "; direction=", signif(braak_sum$direction_concordance, 6)
    )
  } else {
    "Braak summary missing"
  }
)

add_check(
  "CERAD_APOE_effect_high_concordance",
  nrow(cerad_sum) == 1L &&
    cerad_sum$pearson_effect > 0.97 &&
    cerad_sum$spearman_effect > 0.97 &&
    cerad_sum$direction_concordance > 0.95,
  if (nrow(cerad_sum) == 1L) {
    paste0(
      "Pearson=", signif(cerad_sum$pearson_effect, 6),
      "; Spearman=", signif(cerad_sum$spearman_effect, 6),
      "; direction=", signif(cerad_sum$direction_concordance, 6)
    )
  } else {
    "CERAD summary missing"
  }
)

## Observed canonical FDR counts from the validated input run.
add_check(
  "late_FDR_counts_21_to_18",
  nrow(late_sum) == 1L &&
    late_sum$n_FDR_without_APOE == 21L &&
    late_sum$n_FDR_with_APOE == 18L,
  if (nrow(late_sum) == 1L) {
    paste0(
      late_sum$n_FDR_without_APOE,
      " -> ",
      late_sum$n_FDR_with_APOE
    )
  } else {
    "late summary missing"
  }
)

add_check(
  "Braak_FDR_counts_14_to_7",
  nrow(braak_sum) == 1L &&
    braak_sum$n_FDR_without_APOE == 14L &&
    braak_sum$n_FDR_with_APOE == 7L,
  if (nrow(braak_sum) == 1L) {
    paste0(
      braak_sum$n_FDR_without_APOE,
      " -> ",
      braak_sum$n_FDR_with_APOE
    )
  } else {
    "Braak summary missing"
  }
)

add_check(
  "CERAD_FDR_counts_2_to_0",
  nrow(cerad_sum) == 1L &&
    cerad_sum$n_FDR_without_APOE == 2L &&
    cerad_sum$n_FDR_with_APOE == 0L,
  if (nrow(cerad_sum) == 1L) {
    paste0(
      cerad_sum$n_FDR_without_APOE,
      " -> ",
      cerad_sum$n_FDR_with_APOE
    )
  } else {
    "CERAD summary missing"
  }
)

## ---------------------------------------------------------
## Vulnerability-score stability
## ---------------------------------------------------------

add_check(
  "vulnerability_rows_918",
  nrow(vuln) == 918L,
  paste0("observed = ", nrow(vuln))
)

add_check(
  "vulnerability_paired_rows_306",
  nrow(vuln_paired) == 306L,
  paste0("observed = ", nrow(vuln_paired))
)

vuln_score_row <- vuln_summary[
  vuln_summary$metric == "pathology_vulnerability_score",
  ,
  drop = FALSE
]

add_check(
  "vulnerability_score_high_concordance",
  nrow(vuln_score_row) == 1L &&
    vuln_score_row$spearman > 0.98 &&
    vuln_score_row$pearson > 0.97,
  if (nrow(vuln_score_row) == 1L) {
    paste0(
      "Spearman=", signif(vuln_score_row$spearman, 6),
      "; Pearson=", signif(vuln_score_row$pearson, 6)
    )
  } else {
    "vulnerability summary row missing"
  }
)

late_topq <- topq[topq$axis == "late_decline", , drop = FALSE]
braak_topq <- topq[topq$axis == "inverse_braak", , drop = FALSE]

add_check(
  "late_top_quartile_retention_76_of_77",
  nrow(late_topq) == 1L &&
    late_topq$n_topq_without_APOE == 77L &&
    late_topq$n_topq_with_APOE == 77L &&
    late_topq$n_overlap == 76L,
  if (nrow(late_topq) == 1L) {
    paste0(
      "overlap=", late_topq$n_overlap,
      "/", late_topq$n_topq_without_APOE,
      " (", signif(late_topq$pct_original_topq_retained, 4), "%)"
    )
  } else {
    "late top-quartile row missing"
  }
)

add_check(
  "Braak_top_quartile_retention_71_of_77",
  nrow(braak_topq) == 1L &&
    braak_topq$n_topq_without_APOE == 77L &&
    braak_topq$n_topq_with_APOE == 77L &&
    braak_topq$n_overlap == 71L,
  if (nrow(braak_topq) == 1L) {
    paste0(
      "overlap=", braak_topq$n_overlap,
      "/", braak_topq$n_topq_without_APOE,
      " (", signif(braak_topq$pct_original_topq_retained, 4), "%)"
    )
  } else {
    "Braak top-quartile row missing"
  }
)

## ---------------------------------------------------------
## Cohort restriction is explicitly separated from APOE
## ---------------------------------------------------------

add_check(
  "full_vs_subset_three_endpoints",
  nrow(full_subset) == 3L &&
    setequal(
      full_subset$endpoint,
      c("Late_AD_minus_MCI", "Braak", "CERAD")
    ),
  paste(full_subset$endpoint, collapse = ";")
)

add_check(
  "full_vs_subset_all_306_shared",
  all(full_subset$n_shared == 306L),
  paste(
    paste(full_subset$endpoint, full_subset$n_shared, sep = "="),
    collapse = "; "
  )
)

## ---------------------------------------------------------
## Provenance
## ---------------------------------------------------------

add_check(
  "provenance_APOE_definition",
  identical(get_prov("APOE_carrier_definition"), "24,34,44"),
  paste0("value = ", get_prov("APOE_carrier_definition"))
)

add_check(
  "provenance_stage_covariates",
  identical(
    get_prov("stage_primary_covariates"),
    "age_num;sex_factor;pmi_num;batch_factor"
  ),
  paste0("value = ", get_prov("stage_primary_covariates"))
)

add_check(
  "provenance_APOE_factor_covariate",
  grepl(
    "apoe4_carrier_factor",
    get_prov("stage_APOE_covariates"),
    fixed = TRUE
  ),
  paste0("value = ", get_prov("stage_APOE_covariates"))
)

add_check(
  "provenance_direct_Braak",
  identical(get_prov("Braak_adjusted_for_CERAD"), "FALSE"),
  paste0("value = ", get_prov("Braak_adjusted_for_CERAD"))
)

## ---------------------------------------------------------
## Final report
## ---------------------------------------------------------

validation <- do.call(rbind, checks)
print(validation, row.names = FALSE)

n_passed <- sum(validation$passed)
n_total <- nrow(validation)

cat(
  "\nValidation: ",
  n_passed,
  " / ",
  n_total,
  " checks passed\n",
  sep = ""
)

failed <- validation[!validation$passed, , drop = FALSE]

if (nrow(failed) > 0L) {
  cat("\nFAILED CHECKS:\n")
  print(failed, row.names = FALSE)
  stop(
    "APOE sensitivity validation failed: ",
    nrow(failed),
    " check(s).",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("R3.6 APOE SENSITIVITY VALIDATED\n")
cat("============================================================\n")
