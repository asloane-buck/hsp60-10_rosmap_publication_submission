#!/usr/bin/env Rscript

############################################################
## 22_validate_conventional_stage_differential_analysis.R
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

source("R/00_config.R")

out_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential"
)

required_files <- c(
  "02_RNA_DESeq2_all_contrasts.csv",
  "03_protein_limma_all_contrasts.csv",
  "04_reviewer_contrast_summary.csv",
  "05_HSPD1_HSPE1_differential_results.csv",
  "06_Hsp60_10_client_differential_results.csv",
  "08_cross_modal_shared_features_side_by_side_NO_SUBTRACTION.csv",
  "09_cross_modal_Hsp60_10_clients_side_by_side_NO_SUBTRACTION.csv",
  "10_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv",
  "11_sample_stage_audit.csv",
  "12_model_specification_audit.csv",
  "14_REVIEWER_REPORT_SUMMARY.txt",
  "15_validation_checks.csv",
  "sessionInfo.txt",
  "run_provenance.txt"
)

missing <- required_files[
  !file.exists(file.path(out_dir, required_files))
]

if (length(missing) > 0) {
  stop(
    "Missing conventional differential-analysis output(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

rna <- readr::read_csv(
  file.path(out_dir, "02_RNA_DESeq2_all_contrasts.csv"),
  show_col_types = FALSE
)
protein <- readr::read_csv(
  file.path(out_dir, "03_protein_limma_all_contrasts.csv"),
  show_col_types = FALSE
)
summary_tbl <- readr::read_csv(
  file.path(out_dir, "04_reviewer_contrast_summary.csv"),
  show_col_types = FALSE
)
core <- readr::read_csv(
  file.path(out_dir, "05_HSPD1_HSPE1_differential_results.csv"),
  show_col_types = FALSE
)
cross_hsp <- readr::read_csv(
  file.path(out_dir, "09_cross_modal_Hsp60_10_clients_side_by_side_NO_SUBTRACTION.csv"),
  show_col_types = FALSE
)
sample_audit <- readr::read_csv(
  file.path(out_dir, "11_sample_stage_audit.csv"),
  show_col_types = FALSE
)
model_audit <- readr::read_csv(
  file.path(out_dir, "12_model_specification_audit.csv"),
  show_col_types = FALSE
)
internal_validation <- readr::read_csv(
  file.path(out_dir, "15_validation_checks.csv"),
  show_col_types = FALSE
)

expected_contrasts <- c(
  "MCI_vs_NCI",
  "AD_vs_NCI",
  "AD_vs_MCI"
)

get_stage_n <- function(modality, stage) {
  x <- sample_audit |>
    dplyr::filter(
      .data$modality == modality,
      .data$clinical_stage == stage
    ) |>
    dplyr::pull(.data$n)

  if (length(x) != 1) return(NA_integer_)
  as.integer(x)
}

checks <- tibble::tibble(
  check = c(
    "RNA contains all three contrasts",
    "Protein contains all three contrasts",
    "Contrast summary has 6 modality-contrast rows",
    "HSPD1/HSPE1 has 12 expected rows",
    "RNA NCI sample count",
    "RNA MCI sample count",
    "RNA AD sample count",
    "Protein NCI sample count",
    "Protein MCI sample count",
    "Protein AD sample count",
    "RNA method is DESeq2",
    "Protein method is limma",
    "RNA raw-count input documented",
    "Protein raw processed TMT input documented",
    "Cross-modal Hsp client table has no subtraction column",
    "All internal validation checks passed",
    "Reviewer report exists and is nonempty"
  ),
  observed = c(
    paste(sort(unique(rna$contrast)), collapse = ";"),
    paste(sort(unique(protein$contrast)), collapse = ";"),
    nrow(summary_tbl),
    nrow(core),
    get_stage_n("RNA", "NCI"),
    get_stage_n("RNA", "MCI"),
    get_stage_n("RNA", "AD"),
    get_stage_n("Protein", "NCI"),
    get_stage_n("Protein", "MCI"),
    get_stage_n("Protein", "AD"),
    model_audit$method[model_audit$modality == "RNA"],
    model_audit$method[model_audit$modality == "Protein"],
    model_audit$abundance_input[model_audit$modality == "RNA"],
    model_audit$abundance_input[model_audit$modality == "Protein"],
    !any(
      grepl(
        "minus|difference",
        colnames(cross_hsp),
        ignore.case = TRUE
      )
    ),
    all(internal_validation$passed),
    file.info(
      file.path(out_dir, "14_REVIEWER_REPORT_SUMMARY.txt")
    )$size > 0
  ),
  expected = c(
    paste(sort(expected_contrasts), collapse = ";"),
    paste(sort(expected_contrasts), collapse = ";"),
    "6",
    "12",
    "200",
    "158",
    "220",
    "167",
    "96",
    "109",
    "DESeq2",
    "limma",
    "raw integer RNA-seq counts",
    "processed non-residualized log2 TMT abundance",
    "TRUE",
    "TRUE",
    "TRUE"
  )
) |>
  dplyr::mutate(
    passed = as.character(.data$observed) == as.character(.data$expected)
  )

print(checks, n = Inf)

validation_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential_validation"
)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  checks,
  file.path(validation_dir, "validation_checks.csv")
)

if (any(!checks$passed)) {
  stop(
    "Conventional stage differential-analysis validation FAILED.",
    call. = FALSE
  )
}

cat("\nCONTRAST SUMMARY\n")
print(summary_tbl, n = Inf)

cat("\nHSPD1 / HSPE1\n")
print(
  core |>
    dplyr::select(
      .data$modality,
      .data$contrast,
      .data$gene_symbol,
      .data$effect,
      .data$conf_low,
      .data$conf_high,
      .data$p_value,
      .data$fdr
    ),
  n = Inf
)

cat("\nALL CONVENTIONAL STAGE DIFFERENTIAL-ANALYSIS CHECKS PASSED.\n")
