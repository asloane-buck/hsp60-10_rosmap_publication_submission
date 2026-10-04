############################################################
## Reviewer revision audit 01
## Matched RNA/protein participant and stage harmonization
##
## Purpose:
## 1. Identify the exact participants represented in both ROSMAP RNA-seq
##    and TMT proteomics matrices.
## 2. Assign ONE common clinical stage to each participant using the
##    analysis_meta/cogdx-derived diagnosis_stage variable.
## 3. Compare that common stage with the legacy proteomics
##    Control/AsymAD/AD label used in the submitted figures.
## 4. Audit duplicates, missing diagnoses, and the source of unequal
##    RNA/protein stage counts in the previous matched-individual figure.
##
## Run from repository root:
##   Rscript R/reviewer_revisions/01_audit_matched_rna_protein_stage_harmonization.R
############################################################

options(stringsAsFactors = FALSE)

if (!file.exists(file.path("R", "00_config.R"))) {
  stop(
    "Run this script from the repository root. Expected to find R/00_config.R.",
    call. = FALSE
  )
}

source(file.path("R", "00_config.R"))
source(file.path("R", "01_utils.R"))
source(file.path("R", "02_load_data.R"))
source(file.path("R", "03_build_adjusted_core_objects.R"))

revision_dir <- file.path(cfg$output_root, "reviewer_revisions", "stage_harmonization")
dir.create(revision_dir, recursive = TRUE, showWarnings = FALSE)

write_revision <- function(x, filename) {
  out <- file.path(revision_dir, filename)
  readr::write_csv(x, out)
  message("Wrote: ", out)
  invisible(x)
}

collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = ";")
}

harmonize_common_stage <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("Control", "NCI") ~ "NCI",
    x %in% c("Early_AD", "MCI", "Intermediate") ~ "MCI",
    x == "AD" ~ "AD",
    TRUE ~ NA_character_
  )
}

harmonize_legacy_protein_stage <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "Control" ~ "NCI",
    x == "AsymAD" ~ "MCI",
    x == "AD" ~ "AD",
    TRUE ~ NA_character_
  )
}

# ------------------------------------------------------------
# 1. One common cogdx-derived diagnosis per participant
# ------------------------------------------------------------

common_stage_raw <- analysis_meta |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    diagnosis_stage_raw = as.character(.data$diagnosis_stage),
    common_stage = harmonize_common_stage(.data$diagnosis_stage)
  ) |>
  dplyr::filter(!is.na(.data$individual_id), nzchar(.data$individual_id))

common_stage_conflicts <- common_stage_raw |>
  dplyr::filter(!is.na(.data$common_stage)) |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_distinct_common_stages = dplyr::n_distinct(.data$common_stage),
    common_stage_values = collapse_unique(.data$common_stage),
    .groups = "drop"
  ) |>
  dplyr::filter(.data$n_distinct_common_stages > 1)

write_revision(common_stage_conflicts, "01_common_stage_conflicts.csv")

if (nrow(common_stage_conflicts) > 0) {
  stop(
    "At least one participant has multiple cogdx-derived common stages. See 01_common_stage_conflicts.csv before proceeding.",
    call. = FALSE
  )
}

common_stage_by_id <- common_stage_raw |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    common_stage = {
      vals <- unique(.data$common_stage[!is.na(.data$common_stage)])
      if (length(vals) == 0) NA_character_ else vals[[1]]
    },
    diagnosis_stage_raw = collapse_unique(.data$diagnosis_stage_raw),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 2. Matrix-backed RNA and protein sample maps
# ------------------------------------------------------------

rna_map <- rna_meta_adj |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    rna_sample_id = as.character(.data$sample_id),
    rna_stage_from_metadata = harmonize_common_stage(.data$diagnosis_stage)
  ) |>
  dplyr::filter(
    .data$rna_sample_id %in% colnames(rna_mat_adj),
    !is.na(.data$individual_id),
    nzchar(.data$individual_id)
  )

protein_map <- prot_meta_adj |>
  dplyr::transmute(
    individual_id = as.character(.data$IndividualID),
    protein_sample_id = as.character(.data$SampleID),
    protein_legacy_dx = as.character(.data$EmoryStrictDx.2019),
    protein_legacy_stage = harmonize_legacy_protein_stage(.data$EmoryStrictDx.2019)
  ) |>
  dplyr::filter(
    .data$protein_sample_id %in% colnames(prot_mat_adj),
    !is.na(.data$individual_id),
    nzchar(.data$individual_id)
  )

rna_by_id <- rna_map |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_rna_samples = dplyr::n_distinct(.data$rna_sample_id),
    rna_sample_ids = collapse_unique(.data$rna_sample_id),
    rna_stage_from_metadata = collapse_unique(.data$rna_stage_from_metadata),
    .groups = "drop"
  )

protein_by_id <- protein_map |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_protein_samples = dplyr::n_distinct(.data$protein_sample_id),
    protein_sample_ids = collapse_unique(.data$protein_sample_id),
    protein_legacy_dx = collapse_unique(.data$protein_legacy_dx),
    protein_legacy_stage = collapse_unique(.data$protein_legacy_stage),
    .groups = "drop"
  )

rna_duplicates <- rna_by_id |>
  dplyr::filter(.data$n_rna_samples != 1)
protein_duplicates <- protein_by_id |>
  dplyr::filter(.data$n_protein_samples != 1)

write_revision(rna_duplicates, "02_rna_duplicate_or_nonunique_participants.csv")
write_revision(protein_duplicates, "03_protein_duplicate_or_nonunique_participants.csv")

# ------------------------------------------------------------
# 3. Shared participant audit
# ------------------------------------------------------------

shared_ids <- intersect(rna_by_id$individual_id, protein_by_id$individual_id)

pair_audit <- tibble::tibble(individual_id = sort(shared_ids)) |>
  dplyr::left_join(common_stage_by_id, by = "individual_id") |>
  dplyr::left_join(rna_by_id, by = "individual_id") |>
  dplyr::left_join(protein_by_id, by = "individual_id") |>
  dplyr::mutate(
    has_common_stage = !is.na(.data$common_stage),
    one_rna_sample = .data$n_rna_samples == 1,
    one_protein_sample = .data$n_protein_samples == 1,
    strict_pair = .data$has_common_stage & .data$one_rna_sample & .data$one_protein_sample,
    rna_agrees_with_common = is.na(.data$rna_stage_from_metadata) |
      .data$rna_stage_from_metadata == .data$common_stage,
    legacy_protein_label_matches_common = is.na(.data$protein_legacy_stage) |
      .data$protein_legacy_stage == .data$common_stage
  ) |>
  dplyr::arrange(
    factor(.data$common_stage, levels = c("NCI", "MCI", "AD")),
    .data$individual_id
  )

write_revision(pair_audit, "04_matched_participant_stage_audit.csv")

# ------------------------------------------------------------
# 4. Reproduce the OLD matched-figure counting logic
#    to show exactly why RNA and protein group Ns differed.
# ------------------------------------------------------------

old_rna_counts <- rna_map |>
  dplyr::filter(.data$individual_id %in% shared_ids, !is.na(.data$rna_stage_from_metadata)) |>
  dplyr::count(stage = .data$rna_stage_from_metadata, name = "n") |>
  dplyr::mutate(count_definition = "OLD RNA stage assignment")

old_protein_counts <- protein_map |>
  dplyr::filter(.data$individual_id %in% shared_ids, !is.na(.data$protein_legacy_stage)) |>
  dplyr::count(stage = .data$protein_legacy_stage, name = "n") |>
  dplyr::mutate(count_definition = "OLD protein Control/AsymAD/AD assignment")

strict_common_counts <- pair_audit |>
  dplyr::filter(.data$strict_pair) |>
  dplyr::count(stage = .data$common_stage, name = "n") |>
  dplyr::mutate(count_definition = "STRICT shared participants; one common cogdx-derived stage")

count_comparison <- dplyr::bind_rows(
  old_rna_counts,
  old_protein_counts,
  strict_common_counts
) |>
  dplyr::mutate(stage = factor(.data$stage, levels = c("NCI", "MCI", "AD"))) |>
  dplyr::arrange(.data$count_definition, .data$stage)

write_revision(count_comparison, "05_old_vs_harmonized_stage_counts.csv")

# ------------------------------------------------------------
# 5. Cross-tab common clinical stage vs legacy protein label
# ------------------------------------------------------------

stage_crosswalk <- pair_audit |>
  dplyr::filter(.data$strict_pair) |>
  dplyr::count(
    common_stage = .data$common_stage,
    protein_legacy_stage = .data$protein_legacy_stage,
    protein_legacy_dx = .data$protein_legacy_dx,
    name = "n"
  ) |>
  dplyr::arrange(
    factor(.data$common_stage, levels = c("NCI", "MCI", "AD")),
    .data$protein_legacy_dx
  )

write_revision(stage_crosswalk, "06_common_stage_vs_legacy_protein_stage_crosswalk.csv")

# ------------------------------------------------------------
# 6. Strict-pair sample manifest for downstream revised Supp Fig 2
# ------------------------------------------------------------

strict_manifest <- pair_audit |>
  dplyr::filter(.data$strict_pair) |>
  dplyr::transmute(
    individual_id = .data$individual_id,
    common_stage = factor(.data$common_stage, levels = c("NCI", "MCI", "AD")),
    rna_sample_id = .data$rna_sample_ids,
    protein_sample_id = .data$protein_sample_ids,
    protein_legacy_dx = .data$protein_legacy_dx,
    protein_legacy_stage = .data$protein_legacy_stage,
    legacy_protein_label_matches_common = .data$legacy_protein_label_matches_common
  ) |>
  dplyr::arrange(.data$common_stage, .data$individual_id)

write_revision(strict_manifest, "07_strict_matched_rna_protein_manifest.csv")

# ------------------------------------------------------------
# 7. Compact summary for terminal / reviewer-planning audit
# ------------------------------------------------------------

summary_tbl <- tibble::tibble(
  metric = c(
    "RNA matrix-backed unique participants",
    "Protein matrix-backed unique participants",
    "Participants present in BOTH modalities",
    "Shared participants missing common cogdx stage",
    "Shared participants with exactly one RNA + one protein sample and common stage",
    "Strict pairs whose legacy protein stage differs from common cogdx stage"
  ),
  value = c(
    nrow(rna_by_id),
    nrow(protein_by_id),
    length(shared_ids),
    sum(!pair_audit$has_common_stage),
    sum(pair_audit$strict_pair),
    sum(pair_audit$strict_pair & !pair_audit$legacy_protein_label_matches_common, na.rm = TRUE)
  )
)

write_revision(summary_tbl, "00_stage_harmonization_summary.csv")

message("\n============================================================")
message("MATCHED RNA/PROTEIN STAGE HARMONIZATION AUDIT")
message("============================================================")
print(summary_tbl, n = Inf)

message("\nOLD VS HARMONIZED STAGE COUNTS")
print(count_comparison, n = Inf)

message("\nCOMMON cogdx STAGE VS LEGACY PROTEIN LABEL")
print(stage_crosswalk, n = Inf)

message("\nPrimary file to inspect next:")
message(file.path(revision_dir, "04_matched_participant_stage_audit.csv"))
message("\nStrict manifest for the corrected matched-individual analysis:")
message(file.path(revision_dir, "07_strict_matched_rna_protein_manifest.csv"))
message("\nAudit complete.")
