############################################################
## 59_validate_remaining_analysis_gates.R
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "remaining_analysis_gates"
)

required <- c(
  "01_pathology_numeric_category_counts.csv",
  "02_pathology_categorical_vs_continuous_gene_results.csv",
  "03_pathology_categorical_vs_continuous_summary.csv",
  "04_Braak_three_bin_summary.csv",
  "05_Fig1B_numerical_audit_rows.csv",
  "06_Fig1B_scale_audit.csv",
  "07_Fig1B_scale_summary.csv",
  "08_matched_cohort_summary.csv",
  "09_matched_cohort_stage_counts.csv",
  "10_matched_cohort_stage_column_audit.csv",
  "11_FDR_numerical_consistency_audit.csv",
  "12_subtraction_code_inventory.csv",
  "13_subtraction_output_inventory.csv",
  "14_subtraction_policy.csv",
  "15_remaining_gate_status.csv",
  "16_provenance.csv",
  "17_sessionInfo.txt"
)

missing <- required[
  !file.exists(file.path(OUTDIR, required))
]
if (length(missing)) {
  stop("Missing remaining-gate output(s): ", paste(missing, collapse = ", "))
}

read_out <- function(x) read.csv(
  file.path(OUTDIR, x),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts <- read_out("01_pathology_numeric_category_counts.csv")
gene <- read_out("02_pathology_categorical_vs_continuous_gene_results.csv")
sumry <- read_out("03_pathology_categorical_vs_continuous_summary.csv")
bin3 <- read_out("04_Braak_three_bin_summary.csv")
fig1 <- read_out("05_Fig1B_numerical_audit_rows.csv")
fig1checks <- read_out("06_Fig1B_scale_audit.csv")
fig1sum <- read_out("07_Fig1B_scale_summary.csv")
matched <- read_out("08_matched_cohort_summary.csv")
matched_counts <- read_out("09_matched_cohort_stage_counts.csv")
fdr <- read_out("11_FDR_numerical_consistency_audit.csv")
subout <- read_out("13_subtraction_output_inventory.csv")
subpolicy <- read_out("14_subtraction_policy.csv")
status <- read_out("15_remaining_gate_status.csv")
prov <- read_out("16_provenance.csv")

tol <- 1e-10
checks <- list()

add <- function(name, ok, detail) {
  checks[[length(checks)+1L]] <<- data.frame(
    check = name,
    passed = isTRUE(ok),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

getp <- function(item) {
  z <- prov$value[prov$item == item]
  if (length(z) != 1L) return(NA_character_)
  as.character(z)
}

add(
  "pathology_612_gene_endpoint_rows",
  nrow(gene) == 612L &&
    sum(gene$endpoint == "Braak") == 306L &&
    sum(gene$endpoint == "CERAD") == 306L,
  paste0(
    "rows=", nrow(gene),
    "; Braak=", sum(gene$endpoint == "Braak"),
    "; CERAD=", sum(gene$endpoint == "CERAD")
  )
)

add(
  "pathology_two_endpoint_summaries",
  nrow(sumry) == 2L &&
    setequal(sumry$endpoint, c("Braak", "CERAD")),
  paste(sumry$endpoint, collapse = ";")
)

add(
  "continuous_models_reproduce_production",
  all(sumry$max_abs_beta_difference_vs_production < tol) &&
    all(sumry$max_abs_p_difference_vs_production < tol) &&
    all(sumry$max_abs_fdr_difference_vs_production < tol),
  paste0(
    "max beta=", max(sumry$max_abs_beta_difference_vs_production),
    "; max P=", max(sumry$max_abs_p_difference_vs_production),
    "; max FDR=", max(sumry$max_abs_fdr_difference_vs_production)
  )
)

add(
  "categorical_statistics_finite",
  all(is.finite(gene$categorical_omnibus_p)) &&
    all(is.finite(gene$categorical_omnibus_fdr_bh)) &&
    all(is.finite(gene$high_vs_low_numeric_effect)) &&
    all(is.finite(gene$nonlinearity_p)) &&
    all(is.finite(gene$nonlinearity_fdr_bh)),
  "all categorical sensitivity statistics finite"
)

add(
  "categorical_P_FDR_valid",
  all(gene$categorical_omnibus_p >= 0 & gene$categorical_omnibus_p <= 1) &&
    all(gene$categorical_omnibus_fdr_bh >= 0 & gene$categorical_omnibus_fdr_bh <= 1) &&
    all(gene$nonlinearity_p >= 0 & gene$nonlinearity_p <= 1) &&
    all(gene$nonlinearity_fdr_bh >= 0 & gene$nonlinearity_fdr_bh <= 1),
  "all P/FDR values within [0,1]"
)

add(
  "pathology_category_counts_present",
  nrow(counts) >= 6L &&
    all(c("Braak", "CERAD") %in% counts$endpoint),
  paste0("rows=", nrow(counts))
)

add(
  "Fig1B_20_rows",
  nrow(fig1) == 20L,
  paste0("rows=", nrow(fig1))
)

add(
  "Fig1B_all_checks_pass",
  all(fig1checks$passed),
  paste(
    fig1checks$check[!fig1checks$passed],
    collapse = ";"
  )
)

add(
  "Fig1B_zero_arithmetic_error",
  fig1sum$max_abs_point_arithmetic_error < tol,
  paste0("max=", fig1sum$max_abs_point_arithmetic_error)
)

add(
  "strict_matched_n_198",
  nrow(matched) == 1L &&
    matched$n_rows == 198L &&
    matched$n_unique_participants == 198L &&
    matched$duplicate_participant_rows == 0L,
  if (nrow(matched) == 1L) {
    paste0(
      "rows=", matched$n_rows,
      "; unique=", matched$n_unique_participants,
      "; duplicates=", matched$duplicate_participant_rows
    )
  } else {
    "matched summary row missing"
  }
)

add(
  "strict_matched_stage_counts_94_56_48",
  nrow(matched) == 1L &&
    matched$n_NCI == 94L &&
    matched$n_MCI == 56L &&
    matched$n_AD == 48L &&
    matched$counts_match_expected,
  if (nrow(matched) == 1L) {
    paste0(
      "NCI=", matched$n_NCI,
      "; MCI=", matched$n_MCI,
      "; AD=", matched$n_AD
    )
  } else {
    "matched summary row missing"
  }
)

add(
  "FDR_families_all_pass",
  nrow(fdr) >= 8L &&
    all(fdr$passed) &&
    all(fdr$max_abs_fdr_difference < tol),
  paste0(
    "families=", nrow(fdr),
    "; maxdiff=", max(fdr$max_abs_fdr_difference)
  )
)

add(
  "canonical_no_subtraction_outputs_exist",
  sum(
    subout$classification ==
      "preferred_side_by_side_no_subtraction"
  ) >= 3L,
  paste0(
    "preferred files=",
    sum(
      subout$classification ==
        "preferred_side_by_side_no_subtraction"
    )
  )
)

add(
  "Fig1B_subtraction_explicitly_deemphasized",
  any(
    subout$classification ==
      "canonical_Fig1B_descriptive_subtraction_to_deemphasize"
  ) &&
    any(grepl(
      "de-emphasize|remove",
      subpolicy$manuscript_action_later,
      ignore.case = TRUE
    )),
  "descriptive Fig1B subtraction flagged for later manuscript/figure edit"
)

add(
  "four_gate_rows",
  nrow(status) == 4L,
  paste0("rows=", nrow(status))
)

add(
  "no_failed_analysis_gates",
  !any(status$analytical_status == "FAILED"),
  paste(status$analytical_status, collapse = ";")
)

add(
  "three_fully_revalidated_one_edit_pending",
  sum(
    status$analytical_status ==
      "ANALYTICALLY_REVALIDATED"
  ) == 3L &&
    sum(
      status$analytical_status ==
        "ANALYTICALLY_RESOLVED_MANUSCRIPT_EDIT_PENDING"
    ) == 1L,
  paste(
    status$gate,
    status$analytical_status,
    sep = "=",
    collapse = "; "
  )
)

add(
  "Braak_not_CERAD_adjusted",
  identical(getp("Braak_adjusted_for_CERAD"), "FALSE"),
  paste0("value=", getp("Braak_adjusted_for_CERAD"))
)

add(
  "protein_covariates_exact",
  identical(
    getp("protein_covariates"),
    "age_num;sex_factor;pmi_num;batch_factor"
  ),
  paste0("value=", getp("protein_covariates"))
)

validation <- do.call(rbind, checks)
print(validation, row.names = FALSE)

npass <- sum(validation$passed)
ntotal <- nrow(validation)

cat(
  "\nValidation: ",
  npass,
  " / ",
  ntotal,
  " checks passed\n",
  sep = ""
)

failed <- validation[
  !validation$passed,
  ,
  drop = FALSE
]

if (nrow(failed)) {
  cat("\nFAILED CHECKS:\n")
  print(failed, row.names = FALSE)
  stop(
    "Remaining-analysis-gates validation failed: ",
    nrow(failed),
    " check(s).",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("REMAINING ANALYSIS GATES REVALIDATED / RESOLVED\n")
cat("============================================================\n")
cat(
  "Three gates are analytically revalidated. The subtraction gate is\n",
  "analytically resolved, with manuscript/figure de-emphasis still pending.\n",
  sep = ""
)
