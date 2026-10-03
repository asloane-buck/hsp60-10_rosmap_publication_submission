############################################################
## 49_build_APOE_cohort_demographics.R
##
## Reviewer 3.6 / Reviewer 2.1a
## APOE genotype and APOE-e4 carrier bookkeeping for the
## corrected RNA and protein cohorts.
############################################################

options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  if (!file.exists(f)) {
    stop("Missing required production script: ", f)
  }
  source(f)
}

require_objects(
  c(
    "analysis_meta",
    "rna_meta_adj",
    "prot_meta_adj",
    "protein_covars"
  ),
  context = "49_build_APOE_cohort_demographics.R"
)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "APOE_sensitivity"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

## =========================================================
## 1. Build unique participant-level APOE table
## =========================================================

require_columns(
  analysis_meta,
  c("individual_id", "apoe"),
  "analysis_meta"
)

apoe_source <- analysis_meta |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    apoe_raw = suppressWarnings(
      as.integer(.data$apoe)
    )
  ) |>
  dplyr::filter(
    !is.na(.data$individual_id),
    nzchar(.data$individual_id)
  )

## Check for conflicting APOE calls for the same participant.
apoe_conflicts <- apoe_source |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_distinct_nonmissing = dplyr::n_distinct(
      .data$apoe_raw[
        !is.na(.data$apoe_raw)
      ]
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    .data$n_distinct_nonmissing > 1
  )

if (nrow(apoe_conflicts) > 0) {
  stop(
    "Conflicting APOE genotype calls detected for ",
    nrow(apoe_conflicts),
    " participant(s)."
  )
}

apoe_unique <- apoe_source |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    apoe_genotype = {
      vals <- unique(
        .data$apoe_raw[
          !is.na(.data$apoe_raw)
        ]
      )

      if (length(vals) == 0) {
        NA_integer_
      } else {
        vals[1]
      }
    },
    .groups = "drop"
  )

allowed_genotypes <- c(
  22L, 23L, 24L,
  33L, 34L, 44L
)

unexpected_genotypes <- sort(
  unique(
    apoe_unique$apoe_genotype[
      !is.na(apoe_unique$apoe_genotype) &
        !apoe_unique$apoe_genotype %in%
          allowed_genotypes
    ]
  )
)

if (length(unexpected_genotypes) > 0) {
  stop(
    "Unexpected APOE genotype coding detected: ",
    paste(
      unexpected_genotypes,
      collapse = ", "
    )
  )
}

apoe_unique <- apoe_unique |>
  dplyr::mutate(
    apoe_genotype = as.character(
      .data$apoe_genotype
    ),

    apoe4_carrier = dplyr::case_when(
      is.na(.data$apoe_genotype) ~ NA,
      .data$apoe_genotype %in%
        c("24", "34", "44") ~ TRUE,
      .data$apoe_genotype %in%
        c("22", "23", "33") ~ FALSE,
      TRUE ~ NA
    ),

    apoe4_allele_count = dplyr::case_when(
      is.na(.data$apoe_genotype) ~ NA_integer_,
      .data$apoe_genotype %in%
        c("22", "23", "33") ~ 0L,
      .data$apoe_genotype %in%
        c("24", "34") ~ 1L,
      .data$apoe_genotype == "44" ~ 2L,
      TRUE ~ NA_integer_
    )
  )

if (anyDuplicated(apoe_unique$individual_id)) {
  stop(
    "APOE participant table is not unique by individual_id."
  )
}

## =========================================================
## 2. Exact corrected RNA cohort
## =========================================================

require_columns(
  rna_meta_adj,
  c(
    "sample_id",
    "individual_id",
    "clinical_stage"
  ),
  "rna_meta_adj"
)

rna_apoe <- rna_meta_adj |>
  dplyr::transmute(
    sample_id = as.character(.data$sample_id),
    individual_id = as.character(.data$individual_id),
    clinical_stage = as.character(.data$clinical_stage)
  ) |>
  dplyr::left_join(
    apoe_unique,
    by = "individual_id"
  )

if (nrow(rna_apoe) != 577L) {
  stop(
    "Expected corrected RNA cohort n=577; found ",
    nrow(rna_apoe)
  )
}

expected_rna_stage <- c(
  NCI = 200L,
  MCI = 158L,
  AD = 219L
)

observed_rna_stage <- table(
  factor(
    rna_apoe$clinical_stage,
    levels = names(expected_rna_stage)
  )
)

if (
  !identical(
    as.integer(observed_rna_stage),
    as.integer(expected_rna_stage)
  )
) {
  stop(
    "Corrected RNA stage counts changed."
  )
}

## =========================================================
## 3. Full corrected protein cohort
## =========================================================

require_columns(
  prot_meta_adj,
  c(
    "SampleID",
    "individual_id",
    "clinical_stage"
  ),
  "prot_meta_adj"
)

protein_full_apoe <- prot_meta_adj |>
  dplyr::transmute(
    sample_id = as.character(.data$SampleID),
    individual_id = as.character(.data$individual_id),
    clinical_stage = as.character(.data$clinical_stage)
  ) |>
  dplyr::left_join(
    apoe_unique,
    by = "individual_id"
  )

if (nrow(protein_full_apoe) != 400L) {
  stop(
    "Expected full protein cohort n=400; found ",
    nrow(protein_full_apoe)
  )
}

if (anyDuplicated(protein_full_apoe$individual_id)) {
  stop(
    "Full protein cohort is not unique by participant."
  )
}

## =========================================================
## 4. Nuisance-complete protein stage cohort
## =========================================================

stage_levels <- c(
  "NCI",
  "MCI",
  "AD"
)

protein_stage_meta <- prot_meta_adj |>
  dplyr::filter(
    .data$clinical_stage %in%
      stage_levels
  )

protein_complete_covars <- complete.cases(
  protein_stage_meta[
    ,
    c(
      protein_covars,
      "clinical_stage"
    ),
    drop = FALSE
  ]
)

protein_stage_meta <- protein_stage_meta[
  protein_complete_covars,
  ,
  drop = FALSE
]

protein_stage_apoe <- protein_stage_meta |>
  dplyr::transmute(
    sample_id = as.character(.data$SampleID),
    individual_id = as.character(.data$individual_id),
    clinical_stage = as.character(.data$clinical_stage)
  ) |>
  dplyr::left_join(
    apoe_unique,
    by = "individual_id"
  )

if (nrow(protein_stage_apoe) != 372L) {
  stop(
    "Expected nuisance-complete protein stage cohort n=372; found ",
    nrow(protein_stage_apoe)
  )
}

expected_protein_stage <- c(
  NCI = 167L,
  MCI = 96L,
  AD = 109L
)

observed_protein_stage <- table(
  factor(
    protein_stage_apoe$clinical_stage,
    levels = names(expected_protein_stage)
  )
)

if (
  !identical(
    as.integer(observed_protein_stage),
    as.integer(expected_protein_stage)
  )
) {
  stop(
    "Protein nuisance-complete stage counts changed."
  )
}

## =========================================================
## 5. Combine cohorts for summaries
## =========================================================

cohort_long <- dplyr::bind_rows(
  rna_apoe |>
    dplyr::mutate(
      cohort = "RNA_corrected_577"
    ),

  protein_full_apoe |>
    dplyr::mutate(
      cohort = "Protein_full_400"
    ),

  protein_stage_apoe |>
    dplyr::mutate(
      cohort = "Protein_stage_nuisance_complete_372"
    )
)

## =========================================================
## 6. Overall APOE availability / carrier summary
## =========================================================

cohort_summary <- cohort_long |>
  dplyr::group_by(.data$cohort) |>
  dplyr::summarise(
    n_total = dplyr::n(),

    n_apoe_available = sum(
      !is.na(.data$apoe_genotype)
    ),

    n_apoe_missing = sum(
      is.na(.data$apoe_genotype)
    ),

    pct_apoe_available =
      100 *
      .data$n_apoe_available /
      .data$n_total,

    n_e4_noncarrier = sum(
      .data$apoe4_carrier == FALSE,
      na.rm = TRUE
    ),

    n_e4_carrier = sum(
      .data$apoe4_carrier == TRUE,
      na.rm = TRUE
    ),

    pct_e4_carrier_among_available =
      100 *
      .data$n_e4_carrier /
      .data$n_apoe_available,

    n_e4_allele0 = sum(
      .data$apoe4_allele_count == 0,
      na.rm = TRUE
    ),

    n_e4_allele1 = sum(
      .data$apoe4_allele_count == 1,
      na.rm = TRUE
    ),

    n_e4_allele2 = sum(
      .data$apoe4_allele_count == 2,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

## =========================================================
## 7. Genotype counts
## =========================================================

genotype_levels <- c(
  "22", "23", "24",
  "33", "34", "44"
)

genotype_counts <- cohort_long |>
  dplyr::mutate(
    apoe_genotype_display = dplyr::if_else(
      is.na(.data$apoe_genotype),
      "Missing",
      .data$apoe_genotype
    )
  ) |>
  dplyr::count(
    .data$cohort,
    .data$apoe_genotype_display,
    name = "n"
  ) |>
  dplyr::group_by(.data$cohort) |>
  dplyr::mutate(
    pct_of_cohort =
      100 * .data$n / sum(.data$n)
  ) |>
  dplyr::ungroup()

## =========================================================
## 8. Stage-stratified e4 frequencies
## =========================================================

stage_carrier_summary <- cohort_long |>
  dplyr::filter(
    .data$clinical_stage %in%
      stage_levels
  ) |>
  dplyr::group_by(
    .data$cohort,
    .data$clinical_stage
  ) |>
  dplyr::summarise(
    n_total = dplyr::n(),

    n_apoe_available = sum(
      !is.na(.data$apoe4_carrier)
    ),

    n_missing = sum(
      is.na(.data$apoe4_carrier)
    ),

    n_e4_noncarrier = sum(
      .data$apoe4_carrier == FALSE,
      na.rm = TRUE
    ),

    n_e4_carrier = sum(
      .data$apoe4_carrier == TRUE,
      na.rm = TRUE
    ),

    pct_e4_carrier_among_available =
      100 *
      .data$n_e4_carrier /
      .data$n_apoe_available,

    .groups = "drop"
  )

## =========================================================
## 9. Write outputs
## =========================================================

write.csv(
  apoe_unique,
  file.path(
    OUTDIR,
    "APOE_participant_lookup.csv"
  ),
  row.names = FALSE
)

write.csv(
  rna_apoe,
  file.path(
    OUTDIR,
    "APOE_RNA_corrected_cohort_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  protein_full_apoe,
  file.path(
    OUTDIR,
    "APOE_protein_full_cohort_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  protein_stage_apoe,
  file.path(
    OUTDIR,
    "APOE_protein_stage_cohort_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  cohort_summary,
  file.path(
    OUTDIR,
    "APOE_cohort_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  genotype_counts,
  file.path(
    OUTDIR,
    "APOE_genotype_counts_by_cohort.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_carrier_summary,
  file.path(
    OUTDIR,
    "APOE_e4_carrier_by_stage.csv"
  ),
  row.names = FALSE
)

## =========================================================
## 10. Provenance
## =========================================================

git_head <- tryCatch(
  system2(
    "git",
    c("rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) NA_character_
)

provenance <- data.frame(
  item = c(
    "APOE_source_object",
    "APOE_source_column",
    "APOE_genotype_encoding",
    "APOE_e4_carrier_definition",
    "APOE_e4_allele_count_definition",
    "RNA_cohort_n",
    "protein_full_cohort_n",
    "protein_stage_cohort_n",
    "protein_stage_covariates",
    "git_HEAD"
  ),

  value = c(
    "analysis_meta",
    "apoe",
    "22;23;24;33;34;44",
    "24,34,44",
    "0=22/23/33;1=24/34;2=44",
    "577",
    "400",
    "372",
    paste(
      protein_covars,
      collapse = ";"
    ),
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
    "APOE_demographics_provenance.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "APOE_demographics_sessionInfo.txt"
  )
)

## =========================================================
## 11. Console report
## =========================================================

cat(
  "\n============================================================\n"
)

cat(
  "APOE DEMOGRAPHICS — CORRECTED COHORTS\n"
)

cat(
  "============================================================\n\n"
)

cat("Cohort summary:\n")
print(
  cohort_summary,
  n = Inf,
  width = Inf
)

cat("\nStage-stratified e4 carrier summary:\n")
print(
  stage_carrier_summary,
  n = Inf,
  width = Inf
)

cat("\nGenotype counts:\n")
print(
  genotype_counts,
  n = Inf,
  width = Inf
)

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)
