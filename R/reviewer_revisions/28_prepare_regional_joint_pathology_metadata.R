#!/usr/bin/env Rscript

############################################################
## Prepare regional participant metadata with joint
## Braak + CERAD pathology.
##
## IMPORTANT:
##   - Regional cohort = 1,669 individuals.
##   - Authoritative pathology source is the local
##     AMP-AD Diverse Cohorts individual metadata.
##   - CERAD is stored there as amyCerad using C0-C3:
##       C0 = None/No AD
##       C1 = Sparse/Possible
##       C2 = Moderate/Probable
##       C3 = Frequent/Definite
##
##   Existing ROSMAP/AMP-AD protein-analysis convention:
##       cerad_num = 1 = worst pathology
##       cerad_num = 4 = least pathology
##
##   Therefore:
##       cerad_num = 4 - authoritative C code
##
## No production file is written unless all structural
## validation gates pass.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source("R/00_config.R")
source("R/01_utils.R")

cat(
  "\n============================================================\n",
  "REGIONAL JOINT-PATHOLOGY METADATA PREPARATION\n",
  "============================================================\n",
  sep = ""
)

############################################################
## Paths
############################################################

regional_path <- paste0(
  "outputs/reviewer_revisions/regional_braak_fix/",
  "regional_participant_metadata_corrected_braak.csv"
)

authoritative_path <- paste0(
  "../Brain Region Specificity/Metadata/",
  "AMP-AD_DiverseCohorts_individual_metadata.csv"
)

output_dir <- (
  "outputs/reviewer_revisions/regional_joint_pathology_fix"
)

output_path <- file.path(
  output_dir,
  "regional_participant_metadata_joint_pathology.csv"
)

audit_path <- file.path(
  output_dir,
  "regional_joint_pathology_join_audit.csv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## Load regional Braak metadata
############################################################

if (!file.exists(regional_path)) {
  stop(
    "Regional Braak metadata not found:\n",
    regional_path,
    call. = FALSE
  )
}

regional <- readr::read_csv(
  regional_path,
  show_col_types = FALSE
)

cat("\n===== REGIONAL BRAAK METADATA =====\n")
cat("Rows:", nrow(regional), "\n")
cat(
  "Unique individuals:",
  dplyr::n_distinct(regional$individual_id),
  "\n"
)

if (!"individual_id" %in% names(regional)) {
  stop(
    "SAFETY GATE FAILED: regional metadata lacks individual_id.",
    call. = FALSE
  )
}

if (anyDuplicated(regional$individual_id) > 0) {
  stop(
    "SAFETY GATE FAILED: regional metadata contains duplicate ",
    "individual_id values.",
    call. = FALSE
  )
}

############################################################
## Load authoritative AMP-AD individual metadata
############################################################

if (!file.exists(authoritative_path)) {
  stop(
    "Authoritative AMP-AD metadata not found:\n",
    authoritative_path,
    call. = FALSE
  )
}

cat(
  "\nAuthoritative AMP-AD metadata source:\n",
  authoritative_path,
  "\n",
  sep = ""
)

ampad <- readr::read_csv(
  authoritative_path,
  show_col_types = FALSE,
  na = c("", "NA", "NaN"),
  col_types = readr::cols(
    individualID_AMPAD_1.0 = readr::col_character(),
    .default = readr::col_guess()
  )
)

cat("\n===== AMP-AD METADATA =====\n")
cat("Rows:", nrow(ampad), "\n")
cat("Columns:", ncol(ampad), "\n")

required_ampad <- c(
  "individualID",
  "amyCerad",
  "Braak"
)

missing_required <- setdiff(
  required_ampad,
  names(ampad)
)

if (length(missing_required) > 0) {
  stop(
    "SAFETY GATE FAILED: authoritative AMP-AD metadata is missing: ",
    paste(missing_required, collapse = ", "),
    call. = FALSE
  )
}

############################################################
## Authoritative AMP-AD ID audit
############################################################

cat("\n===== AMP-AD ID AUDIT =====\n")

cat(
  "Unique AMP-AD IDs:",
  dplyr::n_distinct(ampad$individualID),
  "\n"
)

if (anyDuplicated(ampad$individualID) > 0) {
  stop(
    "SAFETY GATE FAILED: authoritative AMP-AD metadata contains ",
    "duplicate individualID values.",
    call. = FALSE
  )
}

regional_ids <- unique(regional$individual_id)
ampad_ids <- unique(ampad$individualID)

id_overlap <- intersect(
  regional_ids,
  ampad_ids
)

cat(
  "Regional individuals:",
  length(regional_ids),
  "\n"
)

cat(
  "AMP-AD individuals:",
  length(ampad_ids),
  "\n"
)

cat(
  "ID overlap:",
  length(id_overlap),
  "\n"
)

cat(
  "Regional ID overlap:",
  round(
    100 * length(id_overlap) / length(regional_ids),
    2
  ),
    "%\n"
)

if (length(id_overlap) != length(regional_ids)) {

  missing_ids <- setdiff(
    regional_ids,
    ampad_ids
  )

  cat(
    "\nFirst missing regional IDs:\n"
  )

  print(
    head(missing_ids, 25)
  )

  stop(
    "SAFETY GATE FAILED: not all regional individuals are ",
    "present in authoritative AMP-AD metadata.",
    call. = FALSE
  )
}

############################################################
## Authoritative CERAD coding audit
############################################################

cat("\n===== AUTHORITATIVE CERAD CODING AUDIT =====\n")

cerad_levels <- c(
  "None/No AD/C0",
  "Sparse/Possible/C1",
  "Moderate/Probable/C2",
  "Frequent/Definite/C3"
)

print(
  table(
    ampad$amyCerad,
    useNA = "always"
  )
)

unexpected_cerad <- setdiff(
  unique(
    ampad$amyCerad[
      !is.na(ampad$amyCerad)
    ]
  ),
  c(
    cerad_levels,
    "missing or unknown"
  )
)

if (length(unexpected_cerad) > 0) {
  stop(
    "SAFETY GATE FAILED: unexpected amyCerad values found:\n",
    paste(unexpected_cerad, collapse = "\n"),
    call. = FALSE
  )
}

############################################################
## Build authoritative CERAD variables
############################################################

ampad_pathology <- ampad |>
  dplyr::transmute(
    individual_id = .data$individualID,

    ## Preserve the original AMP-AD pathology label.
    cerad = .data$amyCerad,

    ## Authoritative source coding:
    ## C0 = none/no AD
    ## C1 = sparse/possible
    ## C2 = moderate/probable
    ## C3 = frequent/definite
    cerad_source_c0_c3 = dplyr::case_when(
      .data$amyCerad == "None/No AD/C0" ~ 0,
      .data$amyCerad == "Sparse/Possible/C1" ~ 1,
      .data$amyCerad == "Moderate/Probable/C2" ~ 2,
      .data$amyCerad == "Frequent/Definite/C3" ~ 3,
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::mutate(

    ## Existing analysis convention:
    ##
    ##   1 = Frequent/Definite
    ##   2 = Moderate/Probable
    ##   3 = Sparse/Possible
    ##   4 = None/No AD
    ##
    ## Higher values therefore indicate LESS pathology.
    cerad_num = dplyr::if_else(
      is.finite(.data$cerad_source_c0_c3),
      4 - .data$cerad_source_c0_c3,
      NA_real_
    )
  )

############################################################
## Join
############################################################

regional_joint <- regional |>
  dplyr::left_join(
    ampad_pathology,
    by = "individual_id"
  )

############################################################
## Join validation
############################################################

cat("\n===== POST-JOIN AUDIT =====\n")

joined_rows <- nrow(regional_joint)

cerad_available <- sum(
  is.finite(regional_joint$cerad_num)
)

cerad_missing <- sum(
  !is.finite(regional_joint$cerad_num)
)

cerad_coverage <- (
  cerad_available / joined_rows
)

cat(
  "Regional rows:",
  nrow(regional),
  "\n"
)

cat(
  "Joined rows:",
  joined_rows,
  "\n"
)

cat(
  "CERAD available:",
  cerad_available,
  "\n"
)

cat(
  "CERAD missing/unknown:",
  cerad_missing,
  "\n"
)

cat(
  "CERAD coverage:",
  round(100 * cerad_coverage, 2),
  "%\n"
)

############################################################
## Hard structural gates
############################################################

if (joined_rows != nrow(regional)) {
  stop(
    "SAFETY GATE FAILED: row count changed during metadata join.",
    call. = FALSE
  )
}

if (
  dplyr::n_distinct(regional_joint$individual_id) !=
    dplyr::n_distinct(regional$individual_id)
) {
  stop(
    "SAFETY GATE FAILED: unique individual count changed.",
    call. = FALSE
  )
}

if (cerad_available != 1333L) {
  stop(
    "SAFETY GATE FAILED: expected 1,333 valid CERAD values in ",
    "the authoritative 1,669-person metadata, but found ",
    cerad_available,
    ".",
    call. = FALSE
  )
}

if (cerad_missing != 336L) {
  stop(
    "SAFETY GATE FAILED: expected 336 missing/unknown CERAD values, ",
    "but found ",
    cerad_missing,
    ".",
    call. = FALSE
  )
}

############################################################
## CERAD crosswalk audit
############################################################

cat("\n===== CERAD CROSSWALK =====\n")

print(
  regional_joint |>
    dplyr::filter(
      !is.na(.data$cerad)
    ) |>
    dplyr::count(
      .data$cerad,
      .data$cerad_source_c0_c3,
      .data$cerad_num,
      sort = FALSE
    )
)

############################################################
## Braak/CERAD availability together
############################################################

braak_available <- sum(
  is.finite(regional_joint$braak_num)
)

joint_available <- sum(
  is.finite(regional_joint$braak_num) &
    is.finite(regional_joint$cerad_num)
)

cat("\n===== JOINT PATHOLOGY AVAILABILITY =====\n")

cat(
  "Braak available:",
  braak_available,
  "\n"
)

cat(
  "CERAD available:",
  cerad_available,
  "\n"
)

cat(
  "Both Braak + CERAD available:",
  joint_available,
  "\n"
)

############################################################
## Audit table
############################################################

audit <- tibble::tibble(
  metric = c(
    "regional_rows",
    "regional_unique_individuals",
    "ampad_rows",
    "ampad_unique_individuals",
    "id_overlap",
    "cerad_available",
    "cerad_missing_unknown",
    "cerad_coverage_percent",
    "braak_available",
    "joint_braak_cerad_available"
  ),
  value = c(
    nrow(regional),
    dplyr::n_distinct(regional$individual_id),
    nrow(ampad),
    dplyr::n_distinct(ampad$individualID),
    length(id_overlap),
    cerad_available,
    cerad_missing,
    100 * cerad_coverage,
    braak_available,
    joint_available
  )
)

############################################################
## Write production outputs
############################################################

readr::write_csv(
  regional_joint,
  output_path,
  na = "NA"
)

readr::write_csv(
  audit,
  audit_path,
  na = "NA"
)

cat(
  "\n============================================================\n",
  "PASS: REGIONAL JOINT-PATHOLOGY METADATA PREPARED\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nProduction metadata:\n",
  output_path,
  "\n",
  sep = ""
)

cat(
  "\nJoin audit:\n",
  audit_path,
  "\n",
  sep = ""
)

cat(
  "\nKey results:\n",
  "  Regional individuals: 1,669\n",
  "  AMP-AD ID overlap: 1,669 / 1,669\n",
  "  CERAD available: 1,333\n",
  "  CERAD missing/unknown: 336\n",
  "  CERAD coverage: ",
  round(100 * cerad_coverage, 2),
  "%\n",
  "  Both Braak + CERAD: ",
  joint_available,
  "\n",
  sep = ""
)

cat(
  "\nCERAD convention:\n",
  "  cerad_source_c0_c3: C0 = least pathology -> C3 = greatest pathology\n",
  "  cerad_num: 1 = greatest pathology -> 4 = least pathology\n",
  sep = ""
)

