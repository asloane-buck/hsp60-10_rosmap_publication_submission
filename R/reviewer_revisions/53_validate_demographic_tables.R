############################################################
## 53_validate_demographic_tables.R
##
## Validation for canonical demographic tables.
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "demographic_tables"
)

required_files <- c(
  "demographics_cohort_counts.csv",
  "demographics_stage_group_counts.csv",
  "demographics_continuous_summary.csv",
  "demographics_sex_summary.csv",
  "demographics_APOE_summary.csv",
  "demographics_regional_region_counts.csv",
  "demographics_variable_availability.csv",
  "demographics_manuscript_ready_table.csv",
  "demographics_RNA_manifest.csv",
  "demographics_protein_full_manifest.csv",
  "demographics_protein_primary_manifest.csv",
  "demographics_protein_stage_complete_manifest.csv",
  "demographics_regional_paired_manifest.csv",
  "demographics_provenance.csv",
  "demographics_sessionInfo.txt"
)

missing_files <- required_files[
  !file.exists(file.path(OUTDIR, required_files))
]

if (length(missing_files) > 0L) {
  stop(
    "Missing demographic output(s): ",
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

counts <- read_out("demographics_cohort_counts.csv")
stages <- read_out("demographics_stage_group_counts.csv")
cont <- read_out("demographics_continuous_summary.csv")
sex <- read_out("demographics_sex_summary.csv")
apoe <- read_out("demographics_APOE_summary.csv")
regions <- read_out("demographics_regional_region_counts.csv")
availability <- read_out("demographics_variable_availability.csv")
display <- read_out("demographics_manuscript_ready_table.csv")
rna_manifest <- read_out("demographics_RNA_manifest.csv")
prot_full_manifest <- read_out("demographics_protein_full_manifest.csv")
prot_primary_manifest <- read_out("demographics_protein_primary_manifest.csv")
prot_complete_manifest <- read_out(
  "demographics_protein_stage_complete_manifest.csv"
)
regional_manifest <- read_out("demographics_regional_paired_manifest.csv")
prov <- read_out("demographics_provenance.csv")

checks <- list()

add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

get_count <- function(cohort, group) {
  hit <- counts$n[
    counts$cohort == cohort &
      counts$group == group
  ]
  if (length(hit) != 1L) return(NA_integer_)
  as.integer(hit)
}

get_stage <- function(cohort, group) {
  hit <- stages$n[
    stages$cohort == cohort &
      stages$group == group
  ]
  if (length(hit) != 1L) return(NA_integer_)
  as.integer(hit)
}

get_apoe <- function(cohort, group = "Overall") {
  apoe[
    apoe$cohort == cohort &
      apoe$group == group,
    ,
    drop = FALSE
  ]
}

get_prov <- function(item) {
  hit <- prov$value[prov$item == item]
  if (length(hit) != 1L) return(NA_character_)
  hit
}

## ---------------------------------------------------------
## Core cohort sizes
## ---------------------------------------------------------

add_check(
  "RNA_n_577",
  identical(get_count("RNA_corrected_577", "Overall"), 577L),
  paste0("observed = ", get_count("RNA_corrected_577", "Overall"))
)

add_check(
  "RNA_stage_200_158_219",
  identical(get_stage("RNA_corrected_577", "NCI"), 200L) &&
    identical(get_stage("RNA_corrected_577", "MCI"), 158L) &&
    identical(get_stage("RNA_corrected_577", "AD"), 219L),
  paste0(
    "NCI=", get_stage("RNA_corrected_577", "NCI"),
    "; MCI=", get_stage("RNA_corrected_577", "MCI"),
    "; AD=", get_stage("RNA_corrected_577", "AD")
  )
)

add_check(
  "protein_full_n_400",
  identical(get_count("Protein_full_400", "Overall"), 400L),
  paste0("observed = ", get_count("Protein_full_400", "Overall"))
)

add_check(
  "protein_primary_n_374",
  identical(get_count("Protein_primary_374", "Overall"), 374L),
  paste0("observed = ", get_count("Protein_primary_374", "Overall"))
)

add_check(
  "protein_primary_stage_168_97_109",
  identical(get_stage("Protein_primary_374", "NCI"), 168L) &&
    identical(get_stage("Protein_primary_374", "MCI"), 97L) &&
    identical(get_stage("Protein_primary_374", "AD"), 109L),
  paste0(
    "NCI=", get_stage("Protein_primary_374", "NCI"),
    "; MCI=", get_stage("Protein_primary_374", "MCI"),
    "; AD=", get_stage("Protein_primary_374", "AD")
  )
)

add_check(
  "protein_complete_n_372",
  identical(
    get_count("Protein_stage_nuisance_complete_372", "Overall"),
    372L
  ),
  paste0(
    "observed = ",
    get_count("Protein_stage_nuisance_complete_372", "Overall")
  )
)

add_check(
  "protein_complete_stage_167_96_109",
  identical(
    get_stage("Protein_stage_nuisance_complete_372", "NCI"),
    167L
  ) &&
    identical(
      get_stage("Protein_stage_nuisance_complete_372", "MCI"),
      96L
    ) &&
    identical(
      get_stage("Protein_stage_nuisance_complete_372", "AD"),
      109L
    ),
  paste0(
    "NCI=", get_stage("Protein_stage_nuisance_complete_372", "NCI"),
    "; MCI=", get_stage("Protein_stage_nuisance_complete_372", "MCI"),
    "; AD=", get_stage("Protein_stage_nuisance_complete_372", "AD")
  )
)

add_check(
  "regional_paired_n_215",
  identical(
    get_count("Regional_paired_DLPFC_STG_215", "Overall"),
    215L
  ),
  paste0(
    "observed = ",
    get_count("Regional_paired_DLPFC_STG_215", "Overall")
  )
)

## ---------------------------------------------------------
## Manifest dimensions / uniqueness
## ---------------------------------------------------------

manifest_specs <- list(
  RNA = list(tbl = rna_manifest, n = 577L),
  protein_full = list(tbl = prot_full_manifest, n = 400L),
  protein_primary = list(tbl = prot_primary_manifest, n = 374L),
  protein_complete = list(tbl = prot_complete_manifest, n = 372L),
  regional = list(tbl = regional_manifest, n = 215L)
)

for (nm in names(manifest_specs)) {
  x <- manifest_specs[[nm]]$tbl
  n_expected <- manifest_specs[[nm]]$n

  add_check(
    paste0(nm, "_manifest_n"),
    nrow(x) == n_expected,
    paste0("observed = ", nrow(x))
  )

  add_check(
    paste0(nm, "_manifest_unique_individual"),
    "individual_id" %in% names(x) &&
      anyDuplicated(x$individual_id) == 0L,
    if ("individual_id" %in% names(x)) {
      paste0("duplicates = ", anyDuplicated(x$individual_id))
    } else {
      "individual_id missing"
    }
  )
}

## ---------------------------------------------------------
## APOE counts
## ---------------------------------------------------------

rna_apoe <- get_apoe("RNA_corrected_577")
prot_full_apoe <- get_apoe("Protein_full_400")
prot_primary_apoe <- get_apoe("Protein_primary_374")
prot_complete_apoe <- get_apoe(
  "Protein_stage_nuisance_complete_372"
)

add_check(
  "RNA_APOE_576_available_1_missing",
  nrow(rna_apoe) == 1L &&
    rna_apoe$n_apoe_available == 576L &&
    rna_apoe$n_apoe_missing == 1L,
  if (nrow(rna_apoe) == 1L) {
    paste0(
      rna_apoe$n_apoe_available,
      " available; ",
      rna_apoe$n_apoe_missing,
      " missing"
    )
  } else {
    "row missing"
  }
)

add_check(
  "protein_full_APOE_338_available_62_missing",
  nrow(prot_full_apoe) == 1L &&
    prot_full_apoe$n_apoe_available == 338L &&
    prot_full_apoe$n_apoe_missing == 62L,
  if (nrow(prot_full_apoe) == 1L) {
    paste0(
      prot_full_apoe$n_apoe_available,
      " available; ",
      prot_full_apoe$n_apoe_missing,
      " missing"
    )
  } else {
    "row missing"
  }
)

add_check(
  "protein_primary_APOE_318_available",
  nrow(prot_primary_apoe) == 1L &&
    prot_primary_apoe$n_apoe_available == 318L,
  if (nrow(prot_primary_apoe) == 1L) {
    paste0("available = ", prot_primary_apoe$n_apoe_available)
  } else {
    "row missing"
  }
)

add_check(
  "protein_complete_APOE_317_available",
  nrow(prot_complete_apoe) == 1L &&
    prot_complete_apoe$n_apoe_available == 317L,
  if (nrow(prot_complete_apoe) == 1L) {
    paste0("available = ", prot_complete_apoe$n_apoe_available)
  } else {
    "row missing"
  }
)

regional_apoe <- get_apoe("Regional_paired_DLPFC_STG_215")

add_check(
  "regional_APOE_not_linkable_not_missing",
  nrow(regional_apoe) == 1L &&
    regional_apoe$n_total == 215L &&
    is.na(regional_apoe$n_apoe_available) &&
    is.na(regional_apoe$n_apoe_missing) &&
    is.na(regional_apoe$n_e4_noncarrier) &&
    is.na(regional_apoe$n_e4_carrier),
  if (nrow(regional_apoe) == 1L) {
    paste0(
      "n=", regional_apoe$n_total,
      "; available=", regional_apoe$n_apoe_available,
      "; missing=", regional_apoe$n_apoe_missing
    )
  } else {
    "regional APOE row missing"
  }
)

## ---------------------------------------------------------
## Sex summaries internally consistent
## ---------------------------------------------------------

sex_consistent <- all(
  sex$n_sex_available + sex$n_sex_missing == sex$n_total
) &&
  all(
    sex$n_female + sex$n_male <= sex$n_sex_available
  )

add_check(
  "sex_summary_internal_consistency",
  sex_consistent,
  paste0("rows = ", nrow(sex))
)

## ---------------------------------------------------------
## Continuous summaries internally consistent
## ---------------------------------------------------------

cont_counts_ok <- all(
  cont$n_nonmissing + cont$n_missing == cont$n_total
)

finite_rows <- cont$n_nonmissing > 0L
cont_finite_ok <- all(
  is.finite(cont$mean[finite_rows]) &
    is.finite(cont$median[finite_rows]) &
    is.finite(cont$q1[finite_rows]) &
    is.finite(cont$q3[finite_rows])
)

add_check(
  "continuous_summary_count_consistency",
  cont_counts_ok,
  paste0("rows = ", nrow(cont))
)

add_check(
  "continuous_summary_finite_when_observed",
  cont_finite_ok,
  paste0(
    "observed-variable rows = ",
    sum(finite_rows)
  )
)

## ---------------------------------------------------------
## Regional coverage
## ---------------------------------------------------------

add_check(
  "regional_two_regions",
  nrow(regions) == 2L &&
    setequal(regions$region, c("DLPFC", "STG")),
  paste(regions$region, collapse = ";")
)

add_check(
  "regional_215_each_region",
  nrow(regions) == 2L &&
    all(regions$n_participants == 215L),
  paste(
    paste(regions$region, regions$n_participants, sep = "="),
    collapse = "; "
  )
)

add_check(
  "regional_paired_design_true",
  nrow(regions) == 2L &&
    all(regions$paired_design %in% TRUE),
  paste0(
    "paired = ",
    paste(regions$paired_design, collapse = ";")
  )
)

## ---------------------------------------------------------
## Variable-policy checks
## ---------------------------------------------------------

education_rows <- availability[
  availability$variable == "Education",
  ,
  drop = FALSE
]

add_check(
  "education_not_silently_inferred",
  nrow(education_rows) >= 3L &&
    all(!(education_rows$available %in% TRUE)),
  paste0(
    "rows = ", nrow(education_rows),
    "; available_true = ",
    sum(education_rows$available %in% TRUE)
  )
)

regional_cerad <- availability[
  availability$cohort == "Regional_paired_DLPFC_STG_215" &
    availability$variable == "CERAD",
  ,
  drop = FALSE
]

add_check(
  "regional_bScore_not_relabelled_CERAD",
  nrow(regional_cerad) == 1L &&
    !(regional_cerad$available %in% TRUE) &&
    grepl(
      "not inferred",
      regional_cerad$note,
      ignore.case = TRUE
    ),
  if (nrow(regional_cerad) == 1L) {
    regional_cerad$note
  } else {
    "regional CERAD policy row missing"
  }
)

## ---------------------------------------------------------
## Display table
## ---------------------------------------------------------

expected_display_rows <- c(
  "N",
  "Age at death, mean (SD)",
  "Female, n (%)",
  "PMI, median [IQR]",
  "RIN, mean (SD)",
  "Braak, median [IQR]",
  "CERAD, median [IQR]",
  "APOE e4 carrier, n/N (%)",
  "APOE missing, n"
)

add_check(
  "display_table_expected_rows",
  identical(display$row, expected_display_rows),
  paste(display$row, collapse = ";")
)

add_check(
  "display_table_has_core_columns",
  all(
    c(
      "RNA_corrected_577__Overall",
      "RNA_corrected_577__NCI",
      "RNA_corrected_577__MCI",
      "RNA_corrected_577__AD",
      "Protein_full_400__Overall",
      "Protein_primary_374__Overall",
      "Protein_stage_nuisance_complete_372__Overall",
      "Regional_paired_DLPFC_STG_215__Overall"
    ) %in% names(display)
  ),
  paste(names(display), collapse = ";")
)

## ---------------------------------------------------------
## Provenance
## ---------------------------------------------------------

add_check(
  "provenance_RNA_counts",
  identical(
    get_prov("RNA_stage_counts"),
    "NCI=200;MCI=158;AD=219"
  ),
  paste0("value = ", get_prov("RNA_stage_counts"))
)

add_check(
  "provenance_protein_primary_counts",
  identical(
    get_prov("protein_primary_stage_counts"),
    "NCI=168;MCI=97;AD=109"
  ),
  paste0("value = ", get_prov("protein_primary_stage_counts"))
)

add_check(
  "provenance_protein_complete_counts",
  identical(
    get_prov("protein_stage_complete_counts"),
    "NCI=167;MCI=96;AD=109"
  ),
  paste0("value = ", get_prov("protein_stage_complete_counts"))
)

add_check(
  "provenance_regional_215_paired",
  identical(
    get_prov("regional_regions"),
    "DLPFC=215;STG=215;paired=TRUE"
  ),
  paste0("value = ", get_prov("regional_regions"))
)

add_check(
  "provenance_protein_covariates",
  identical(
    get_prov("protein_nuisance_covariates"),
    "age_num;sex_factor;pmi_num;batch_factor"
  ),
  paste0("value = ", get_prov("protein_nuisance_covariates"))
)

## ---------------------------------------------------------
## Final report
## ---------------------------------------------------------

validation <- do.call(rbind, checks)

print(
  validation,
  row.names = FALSE
)

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

failed <- validation[
  !validation$passed,
  ,
  drop = FALSE
]

if (nrow(failed) > 0L) {
  cat("\nFAILED CHECKS:\n")
  print(failed, row.names = FALSE)

  stop(
    "Demographic-table validation failed: ",
    nrow(failed),
    " check(s).",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("CANONICAL DEMOGRAPHIC TABLES VALIDATED\n")
cat("============================================================\n")
