#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

source("R/00_config.R")

out_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential_final"
)

required_files <- c(
  "04_RNA_DESeq2_all_contrasts.csv",
  "05_protein_feature_universe_manifest.csv",
  "06_protein_feature_universe_counts.csv",
  "07_protein_PRIMARY_complete400_limma.csv",
  "08_protein_SENSITIVITY_missing_batches_le2_limma.csv",
  "09_protein_SENSITIVITY_missingness_le20_limma.csv",
  "10_protein_limma_warning_audit.csv",
  "11_protein_contrast_estimability_audit.csv",
  "12_PRIMARY_reviewer_contrast_summary.csv",
  "13_protein_sensitivity_contrast_summary.csv",
  "14_PRIMARY_Hsp60_10_client_differential_results.csv",
  "15_PRIMARY_HSPD1_HSPE1_results.csv",
  "17_PRIMARY_cross_modal_Hsp60_10_side_by_side_NO_SUBTRACTION.csv",
  "18_PRIMARY_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv",
  "19_protein_sensitivity_stability_vs_primary.csv",
  "25_Hsp60_10_primary_universe_manifest.csv",
  "26_Hsp60_10_primary_universe_counts.csv",
  "20_sample_stage_audit.csv",
  "21_analysis_policy_audit.csv",
  "22_model_specification_audit.csv",
  "23_REVIEWER_REPORT_SUMMARY.txt",
  "24_validation_checks.csv",
  "sessionInfo.txt",
  "run_provenance.txt"
)

missing <- required_files[!file.exists(file.path(out_dir, required_files))]
if (length(missing) > 0) stop("Missing final differential-analysis outputs: ", paste(missing, collapse = ", "), call. = FALSE)

read_final <- function(filename) {
  readr::read_csv(file.path(out_dir, filename), show_col_types = FALSE)
}

rna <- read_final("04_RNA_DESeq2_all_contrasts.csv")
primary <- read_final("07_protein_PRIMARY_complete400_limma.csv")
batch2 <- read_final("08_protein_SENSITIVITY_missing_batches_le2_limma.csv")
missing20 <- read_final("09_protein_SENSITIVITY_missingness_le20_limma.csv")
universe_counts <- read_final("06_protein_feature_universe_counts.csv")
warnings <- read_final("10_protein_limma_warning_audit.csv")
estimability <- read_final("11_protein_contrast_estimability_audit.csv")
summary_tbl <- read_final("12_PRIMARY_reviewer_contrast_summary.csv")
sample_audit <- read_final("20_sample_stage_audit.csv")
policy <- read_final("21_analysis_policy_audit.csv")
core <- read_final("15_PRIMARY_HSPD1_HSPE1_results.csv")
cross_hsp <- read_final("17_PRIMARY_cross_modal_Hsp60_10_side_by_side_NO_SUBTRACTION.csv")
hsp_manifest <- read_final("25_Hsp60_10_primary_universe_manifest.csv")
hsp_counts <- read_final("26_Hsp60_10_primary_universe_counts.csv")
internal <- read_final("24_validation_checks.csv")

expected_contrasts <- sort(c("MCI_vs_NCI", "AD_vs_NCI", "AD_vs_MCI"))

get_universe_n <- function(name) {
  x <- universe_counts$n_features[universe_counts$universe == name]
  if (length(x) != 1) return(NA_integer_)
  as.integer(x)
}

get_stage_n <- function(modality, stage) {
  x <- sample_audit$n[sample_audit$modality == modality & sample_audit$clinical_stage == stage]
  if (length(x) != 1) return(NA_integer_)
  as.integer(x)
}

primary_warning <- warnings |>
  dplyr::filter(.data$universe == "primary_complete400")

primary_estimability <- estimability |>
  dplyr::filter(.data$universe == "primary_complete400")

get_hsp_category_n <- function(category_name) {
  x <- hsp_counts$n[
    hsp_counts$category == category_name
  ]
  if (length(x) != 1) return(NA_integer_)
  as.integer(x)
}

rna_hsp_total <- sum(
  hsp_manifest$in_rna_primary,
  na.rm = TRUE
)

protein_hsp_total <- sum(
  hsp_manifest$in_protein_primary,
  na.rm = TRUE
)

checks <- tibble::tibble(
  check = c(
    "RNA has 3 contrasts",
    "Primary protein has 3 contrasts",
    "Batch<=2 sensitivity has 3 contrasts",
    "<=20% missing sensitivity has 3 contrasts",
    "Primary protein universe n",
    "Batch<=2 protein universe n",
    "<=20% missing protein universe n",
    "RNA NCI n",
    "RNA MCI n",
    "RNA AD n",
    "Protein NCI n",
    "Protein MCI n",
    "Protein AD n",
    "Primary protein has zero NA stage coefficients",
    "Primary protein contrasts all FDR-testable",
    "Primary summary has 6 modality/contrast rows",
    "HSPD1/HSPE1 primary table has 12 rows",
    "Cross-modal Hsp table has no subtraction column",
    "Canonical Hsp60/10 manifest total",
    "RNA primary Hsp60/10 total",
    "Protein primary Hsp60/10 total",
    "Shared primary Hsp60/10 total",
    "RNA-only primary Hsp60/10 total",
    "Protein-only primary Hsp60/10 total",
    "Neither-primary Hsp60/10 total",
    "Internal validation passed",
    "Policy labels primary protein complete400",
    "Reviewer report is nonempty"
  ),
  observed = c(
    paste(sort(unique(rna$contrast)), collapse = ";"),
    paste(sort(unique(primary$contrast)), collapse = ";"),
    paste(sort(unique(batch2$contrast)), collapse = ";"),
    paste(sort(unique(missing20$contrast)), collapse = ";"),
    get_universe_n("primary_complete400"),
    get_universe_n("sensitivity_missing_batches_le2"),
    get_universe_n("sensitivity_missingness_le20"),
    get_stage_n("RNA", "NCI"),
    get_stage_n("RNA", "MCI"),
    get_stage_n("RNA", "AD"),
    get_stage_n("Protein", "NCI"),
    get_stage_n("Protein", "MCI"),
    get_stage_n("Protein", "AD"),
    primary_warning$n_features_with_na_stage_coefficient,
    all(primary_estimability$n_fdr_finite == primary_estimability$n_features),
    nrow(summary_tbl),
    nrow(core),
    !any(grepl("minus|difference", colnames(cross_hsp), ignore.case = TRUE)),
    nrow(hsp_manifest),
    rna_hsp_total,
    protein_hsp_total,
    get_hsp_category_n("Both"),
    get_hsp_category_n("RNA only"),
    get_hsp_category_n("Protein only"),
    get_hsp_category_n("Neither"),
    all(internal$passed),
    any(policy$role == "primary" & policy$universe == "primary_complete400" & policy$n_features == 5088),
    file.info(file.path(out_dir, "23_REVIEWER_REPORT_SUMMARY.txt"))$size > 0
  ),
  expected = c(
    paste(expected_contrasts, collapse = ";"),
    paste(expected_contrasts, collapse = ";"),
    paste(expected_contrasts, collapse = ";"),
    paste(expected_contrasts, collapse = ";"),
    "5088",
    "5899",
    "7159",
    "200",
    "158",
    "219",
    "167",
    "96",
    "109",
    "0",
    "TRUE",
    "6",
    "12",
    "TRUE",
    "321",
    "297",
    "274",
    "257",
    "40",
    "17",
    "7",
    "TRUE",
    "TRUE",
    "TRUE"
  )
) |>
  dplyr::mutate(passed = as.character(.data$observed) == as.character(.data$expected))

print(checks, n = Inf)

validation_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential_final_validation"
)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(checks, file.path(validation_dir, "validation_checks.csv"))

if (any(!checks$passed)) stop("FINAL conventional stage differential-analysis validation FAILED.", call. = FALSE)

cat("\nPRIMARY REVIEWER SUMMARY\n")
print(tibble::as_tibble(summary_tbl), n = Inf)

cat("\nPROTEIN WARNING AUDIT\n")
print(tibble::as_tibble(warnings), n = Inf)

cat("\nHSPD1 / HSPE1 PRIMARY RESULTS\n")
print(
  tibble::as_tibble(
    core |>
      dplyr::select("modality", "contrast", "gene_symbol", "effect", "conf_low", "conf_high", "p_value", "fdr")
  ),
  n = Inf
)

cat("\nALL FINAL CONVENTIONAL STAGE DIFFERENTIAL-ANALYSIS CHECKS PASSED.\n")
