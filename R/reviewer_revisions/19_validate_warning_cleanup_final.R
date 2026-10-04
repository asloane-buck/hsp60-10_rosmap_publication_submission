#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")

read_script <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

runner <- read_script("R/supplemental/90_run_supplemental_pipeline.R")
s2 <- read_script("R/supplemental/11_make_supplementary_figure_2_matched_individual_sensitivity.R")
s3 <- read_script("R/supplemental/12_make_supplementary_figure_3_mitochondrial_specificity.R")
s4 <- read_script("R/supplemental/13_make_supplementary_figure_4_pathology_model_robustness.R")
s5 <- read_script("R/supplemental/14_make_supplementary_figure_5_msbb_cross_cohort_validation.R")

checks <- tibble::tribble(
  ~check, ~passed,
  "Runner cognition select uses character names",
  grepl(
    'dplyr::select\\("gene", "cognition_priority_score"\\)',
    runner
  ),
  "Supp Fig 2 pivot_wider tidyselect cleaned",
  grepl(
    'pivot_wider\\(names_from = "Modality", values_from = "Score"\\)',
    s2
  ),
  "Supp Fig 3 null_mean rename cleaned",
  grepl(
    'rename\\(null_perm_mean = "null_mean"\\)',
    s3
  ),
  "Supp Fig 4 model 1 select cleaned",
  grepl(
    'select\\("y_z", "braak_z", "age_z", "pmi_z", "sex_model", "batch_model"\\)',
    s4
  ),
  "Supp Fig 4 model 2 select cleaned",
  grepl(
    'select\\("y_z", "braak_z", "cerad_z", "age_z", "pmi_z", "sex_model", "batch_model"\\)',
    s4
  ),
  "Supp Fig 5 scalar case_when removed",
  !grepl(
    '"collapse_magnitude_positive" %in% colnames\\(msbb_collapse\\) ~',
    s5
  ),
  "Supp Fig 5 group select cleaned",
  grepl(
    'select\\("group", value = dplyr::all_of\\(value_col\\)\\)',
    s5
  )
)

print(checks, n = Inf)

if (any(!checks$passed)) {
  stop("Final warning-cleanup static validation FAILED.", call. = FALSE)
}

message("FINAL WARNING-CLEANUP STATIC CHECKS PASSED.")
