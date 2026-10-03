############################################################
## 57_validate_reporting_requirements_resolution.R
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "reporting_requirements_resolution"
)

required <- c(
  "01_gap_classification_and_resolution.csv",
  "02_pathology_Braak_CERAD_effect_SE_CI_n_P_FDR.csv",
  "03_pathology_exact_model_validation.csv",
  "04_Fig1B_reporting_complete.csv",
  "05_PC1_reporting_complete.csv",
  "06_network_cognition_reporting_complete.csv",
  "07_regional_Hsp_reporting_complete.csv",
  "08_cognition_null_reporting_complete.csv",
  "09_multiplicity_policy.csv",
  "10_final_reporting_status.csv",
  "11_unresolved_reporting_items.csv",
  "12_provenance.csv",
  "13_sessionInfo.txt"
)

missing <- required[!file.exists(file.path(OUTDIR, required))]
if (length(missing)) {
  stop("Missing reporting-resolution output(s): ", paste(missing, collapse=", "))
}

read_out <- function(x) read.csv(
  file.path(OUTDIR, x),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gaps <- read_out("01_gap_classification_and_resolution.csv")
path <- read_out("02_pathology_Braak_CERAD_effect_SE_CI_n_P_FDR.csv")
val <- read_out("03_pathology_exact_model_validation.csv")
fig1 <- read_out("04_Fig1B_reporting_complete.csv")
pc1 <- read_out("05_PC1_reporting_complete.csv")
network <- read_out("06_network_cognition_reporting_complete.csv")
regional <- read_out("07_regional_Hsp_reporting_complete.csv")
null <- read_out("08_cognition_null_reporting_complete.csv")
policy <- read_out("09_multiplicity_policy.csv")
status <- read_out("10_final_reporting_status.csv")
unresolved <- read_out("11_unresolved_reporting_items.csv")
prov <- read_out("12_provenance.csv")

checks <- list()

add <- function(name, ok, detail) {
  checks[[length(checks)+1L]] <<- data.frame(
    check=name,
    passed=isTRUE(ok),
    detail=as.character(detail),
    stringsAsFactors=FALSE
  )
}

getp <- function(item) {
  x <- prov$value[prov$item == item]
  if (length(x) != 1L) return(NA_character_)
  as.character(x)
}

tol <- 1e-10

add(
  "seven_current_gaps_accounted_for",
  nrow(gaps) == 7L,
  paste0("observed=", nrow(gaps))
)

add(
  "three_false_negatives_four_genuine",
  sum(gaps$gap_type == "audit_detector_false_negative") == 3L &&
    sum(gaps$gap_type == "genuine_reporting_gap") == 4L,
  paste0(
    "false_negative=",
    sum(gaps$gap_type == "audit_detector_false_negative"),
    "; genuine=",
    sum(gaps$gap_type == "genuine_reporting_gap")
  )
)

add(
  "pathology_612_rows",
  nrow(path) == 612L,
  paste0("observed=", nrow(path))
)

add(
  "pathology_306_each_endpoint",
  sum(path$endpoint == "Braak") == 306L &&
    sum(path$endpoint == "CERAD") == 306L,
  paste0(
    "Braak=", sum(path$endpoint == "Braak"),
    "; CERAD=", sum(path$endpoint == "CERAD")
  )
)

stat_cols <- c(
  "effect","std_error","conf_low_95","conf_high_95",
  "n","p_value","fdr_bh"
)

add(
  "pathology_statistics_finite",
  all(vapply(stat_cols, function(x) {
    all(is.finite(path[[x]]))
  }, logical(1))),
  paste(stat_cols, collapse=";")
)

add(
  "pathology_positive_SE",
  all(path$std_error > 0),
  paste0("min_SE=", min(path$std_error))
)

add(
  "pathology_CI_contains_effect",
  all(
    path$conf_low_95 <= path$effect &
      path$effect <= path$conf_high_95
  ),
  "all 612 effects inside 95% CI"
)

add(
  "pathology_valid_P_FDR",
  all(path$p_value >= 0 & path$p_value <= 1) &&
    all(path$fdr_bh >= 0 & path$fdr_bh <= 1),
  "P/FDR within [0,1]"
)

add(
  "pathology_exact_model_zero_difference",
  nrow(val) == 2L &&
    all(val$max_effect_difference <= tol) &&
    all(val$max_p_difference <= tol) &&
    all(val$max_fdr_difference <= tol) &&
    all(val$n_mismatches == 0L),
  paste0(
    "max effect=", max(val$max_effect_difference),
    "; max P=", max(val$max_p_difference),
    "; max FDR=", max(val$max_fdr_difference),
    "; n mismatch=", sum(val$n_mismatches)
  )
)

add(
  "Fig1B_40_rows",
  nrow(fig1) == 40L,
  paste0("observed=", nrow(fig1))
)

add(
  "Fig1B_current_schema",
  all(
    c(
      "pathway", "pathway_label", "pathway_short",
      "modality", "shift", "effect",
      "effect_label", "inferential_status",
      "primary_inference_note"
    ) %in% names(fig1)
  ),
  paste(names(fig1), collapse = ";")
)

fig1_keys <- paste(
  fig1$pathway,
  fig1$modality,
  fig1$shift,
  sep = "||"
)

add(
  "Fig1B_balanced_10x2x2",
  length(unique(fig1$pathway)) == 10L &&
    setequal(unique(fig1$modality), c("RNA", "Protein")) &&
    setequal(unique(fig1$shift), c("NCI to MCI", "MCI to AD")) &&
    length(unique(fig1_keys)) == 40L,
  paste0(
    "pathways=", length(unique(fig1$pathway)),
    "; modalities=", length(unique(fig1$modality)),
    "; shifts=", length(unique(fig1$shift)),
    "; unique keys=", length(unique(fig1_keys))
  )
)

add(
  "Fig1B_effects_finite",
  all(is.finite(fig1$effect)),
  "all current modality-specific stage changes finite"
)

add(
  "Fig1B_no_cross_modal_inference",
  all(fig1$inferential_status == "descriptive_cross_modal_display") &&
    !any(
      c(
        "conf_low_95", "conf_high_95", "p_value",
        "fdr_bh_all_20_fig1b_tests"
      ) %in% names(fig1)
    ),
  "no RNA-protein subtraction CI/P/FDR fields"
)

add(
  "PC1_four_rows",
  nrow(pc1) == 4L,
  paste0("observed=", nrow(pc1))
)

add(
  "PC1_effect_is_delta_mean",
  max(abs(pc1$effect - pc1$delta_mean)) <= tol,
  paste0("maxdiff=", max(abs(pc1$effect-pc1$delta_mean)))
)

add(
  "PC1_effect_matches_mean_difference",
  max(abs(pc1$effect - (pc1$mean2-pc1$mean1))) <= tol,
  "effect = mean2 - mean1"
)

add(
  "PC1_FDR_present",
  all(is.finite(pc1$p_fdr_bh_within_score)),
  "existing within-score BH retained"
)

add(
  "network_three_outcomes",
  nrow(network) == 3L,
  paste0("observed=", nrow(network))
)

add(
  "network_complete_effect_SE_CI_n_P",
  all(is.finite(network$n)) &&
    all(is.finite(network$estimate)) &&
    all(is.finite(network$std.error)) &&
    all(is.finite(network$conf.low)) &&
    all(is.finite(network$conf.high)) &&
    all(is.finite(network$p.value)),
  "all 3 outcomes complete"
)

add(
  "network_FDR_valid",
  all(
    network$fdr_bh_three_network_outcomes >= 0 &
      network$fdr_bh_three_network_outcomes <= 1
  ),
  "BH across three outcomes"
)

add(
  "regional_Hsp_two_primary_tests",
  nrow(regional) == 2L,
  paste0("observed=", nrow(regional))
)

add(
  "regional_Hsp_existing_FDR_present",
  all(is.finite(regional$BH_FDR_two_primary_tests)),
  "existing 2-test BH retained"
)

add(
  "regional_Hsp_FDR_reproduces",
  max(abs(
    regional$BH_FDR_two_primary_tests -
      p.adjust(regional$p_region_x_predictor, method="BH")
  )) <= tol,
  "BH exact"
)

add(
  "null_FDR_valid",
  all(
    null$fdr_bh_across_null_metrics >= 0 &
      null$fdr_bh_across_null_metrics <= 1
  ),
  paste0("rows=", nrow(null))
)

add(
  "multiplicity_policy_10_families",
  nrow(policy) == 10L &&
    anyDuplicated(policy$analysis_id) == 0L,
  paste0("rows=", nrow(policy))
)

add(
  "final_status_17_analyses",
  nrow(status) == 17L &&
    anyDuplicated(status$analysis_id) == 0L,
  paste0("rows=", nrow(status))
)

allowed <- c("resolved", "not_applicable")
status_values <- unlist(status[, -1, drop=FALSE])

add(
  "final_status_no_unresolved",
  all(status_values %in% allowed),
  paste(sort(unique(status_values)), collapse=";")
)

add(
  "unresolved_table_empty",
  nrow(unresolved) == 0L,
  paste0("rows=", nrow(unresolved))
)

add(
  "provenance_zero_unresolved",
  identical(getp("unresolved_after_resolution"), "0"),
  paste0("value=", getp("unresolved_after_resolution"))
)

add(
  "provenance_covariates_exact",
  identical(
    getp("protein_covariates"),
    "age_num;sex_factor;pmi_num;batch_factor"
  ),
  paste0("value=", getp("protein_covariates"))
)

validation <- do.call(rbind, checks)
print(validation, row.names=FALSE)

npass <- sum(validation$passed)
ntotal <- nrow(validation)

cat(
  "\nValidation: ", npass, " / ", ntotal,
  " checks passed\n",
  sep=""
)

failed <- validation[!validation$passed,,drop=FALSE]

if (nrow(failed)) {
  cat("\nFAILED CHECKS:\n")
  print(failed, row.names=FALSE)
  stop(
    "Reporting-requirements resolution validation failed: ",
    nrow(failed),
    " check(s).",
    call.=FALSE
  )
}

cat("\n============================================================\n")
cat("REPORTING REQUIREMENTS ANALYTICALLY REVALIDATED\n")
cat("============================================================\n")
