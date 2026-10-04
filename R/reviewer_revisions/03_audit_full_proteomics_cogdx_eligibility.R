############################################################
## Reviewer revision audit 03
## Full ROSMAP proteomics cohort eligibility before legacy-stage filtering
##
## WHY THIS AUDIT EXISTS
## The submitted pipeline filters TMT proteomics samples to legacy
## Control/AsymAD/AD labels inside R/02_load_data.R before downstream
## clinical-stage analyses are performed. If the revised primary stage
## definition is cogdx-derived NCI/MCI/AD, cohort eligibility must be
## evaluated BEFORE that legacy filter.
##
## This script is audit-only. It does not modify the manuscript pipeline.
##
## It answers:
## 1. How many matrix-backed proteomics samples exist before the legacy
##    Control/AsymAD/AD filter?
## 2. How many have cogdx = 1/2/4 in the canonical ROSMAP clinical file?
## 3. Which samples would be gained/lost when eligibility is defined by
##    cogdx rather than the legacy proteomics label?
## 4. Are there duplicate participants or region/tissue complications?
## 5. Does cogdx in the clinical file agree with the value represented in
##    Analysis_Meta_Merged for overlapping participants?
##
## Reproducibility:
## - Run from repository root.
## - Uses repository configuration and canonical utility functions.
## - Uses no absolute hard-coded input paths.
## - Writes manifests, source-path audit, Git provenance, and sessionInfo().
## - No stochastic procedures are used.
##
## Run:
##   Rscript R/reviewer_revisions/03_audit_full_proteomics_cogdx_eligibility.R
############################################################

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Repository checks
# ------------------------------------------------------------

required_scripts <- file.path("R", c("00_config.R", "01_utils.R"))
missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "Run this script from the repository root. Missing required script(s):\n",
    paste(missing_scripts, collapse = "\n"),
    call. = FALSE
  )
}

for (script_i in required_scripts) {
  source(script_i)
}

require_objects(c("cfg"), context = "reviewer audit 03")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "full_proteomics_cogdx_eligibility"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_review_csv <- function(x, filename) {
  out <- file.path(out_dir, filename)
  readr::write_csv(x, out)
  message("Wrote: ", out)
  invisible(x)
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- sort(unique(x))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

clinical_stage_from_cogdx <- function(x) {
  x <- safe_num(x)
  dplyr::case_when(
    x == 1 ~ "NCI",
    x == 2 ~ "MCI",
    x == 4 ~ "AD",
    TRUE ~ NA_character_
  )
}

# ------------------------------------------------------------
# 1. Resolve the exact input files using the same candidates
#    as the manuscript pipeline.
# ------------------------------------------------------------

analysis_meta_file <- first_existing(
  c(
    file.path(cfg$metadata_dir, "Analysis_Meta_Merged.csv"),
    file.path(cfg$metadata_dir, "Analysis_Meta_Merged_copy.csv"),
    file.path(cfg$derived_dir, "Analysis_Meta_Merged_copy.csv")
  ),
  "Analysis_Meta_Merged metadata"
)

protein_matrix_file <- first_existing(
  c(
    file.path(
      cfg$derived_dir,
      "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_copy.csv"
    ),
    file.path(
      cfg$proteomics_input_dir,
      "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"
    ),
    file.path(
      cfg$proteomics_input_dir,
      "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_copy.csv"
    )
  ),
  "full ROSMAP TMT protein matrix"
)

protein_metadata_file <- first_existing(
  c(
    file.path(cfg$derived_dir, "matched_metadata_copy.csv"),
    file.path(cfg$proteomics_input_dir, "matched_metadata.csv"),
    file.path(cfg$proteomics_input_dir, "matched_metadata_copy.csv")
  ),
  "ROSMAP TMT matched metadata"
)

clinical_file <- cfg$rosmap_clinical_file
if (!file.exists(clinical_file)) {
  stop(
    "ROSMAP clinical file not found:\n",
    clinical_file,
    "\nSet ROSMAP_CLINICAL_FILE or place the documented controlled-access file in the configured metadata directory.",
    call. = FALSE
  )
}

source_paths <- tibble::tibble(
  input = c(
    "Analysis_Meta_Merged",
    "Full TMT matrix",
    "TMT matched metadata",
    "ROSMAP clinical"
  ),
  path = c(
    analysis_meta_file,
    protein_matrix_file,
    protein_metadata_file,
    clinical_file
  )
)

write_review_csv(source_paths, "01_input_source_paths.csv")

# ------------------------------------------------------------
# 2. Canonical clinical cogdx table
# ------------------------------------------------------------

clinical_raw <- readr::read_csv(
  clinical_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

clinical_id_col <- require_alias_col(
  clinical_raw,
  c("projid", "individual_id", "individualid", "individualID", "IndividualID"),
  "clinical participant ID",
  context = "ROSMAP clinical metadata"
)

clinical_cogdx_col <- require_alias_col(
  clinical_raw,
  c("cogdx", "COGDX", "diagnosis"),
  "cogdx",
  context = "ROSMAP clinical metadata"
)

clinical_long <- tibble::tibble(
  individual_id = as.character(clinical_raw[[clinical_id_col]]),
  clinical_cogdx = safe_num(clinical_raw[[clinical_cogdx_col]])
) |>
  dplyr::filter(
    !is.na(.data$individual_id),
    nzchar(.data$individual_id)
  )

clinical_conflicts <- clinical_long |>
  dplyr::filter(is.finite(.data$clinical_cogdx)) |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_distinct_cogdx = dplyr::n_distinct(.data$clinical_cogdx),
    cogdx_values = collapse_unique(.data$clinical_cogdx),
    .groups = "drop"
  ) |>
  dplyr::filter(.data$n_distinct_cogdx > 1)

write_review_csv(clinical_conflicts, "02_clinical_cogdx_conflicts.csv")

if (nrow(clinical_conflicts) > 0) {
  stop(
    "The ROSMAP clinical file contains participants with multiple cogdx values. ",
    "Inspect 02_clinical_cogdx_conflicts.csv before defining a canonical stage.",
    call. = FALSE
  )
}

clinical_by_id <- clinical_long |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    clinical_cogdx = {
      vals <- unique(.data$clinical_cogdx[is.finite(.data$clinical_cogdx)])
      if (length(vals) == 0) NA_real_ else vals[[1]]
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    clinical_stage = clinical_stage_from_cogdx(.data$clinical_cogdx),
    clinical_stage_eligible = !is.na(.data$clinical_stage),
    clinical_cogdx_status = dplyr::case_when(
      is.na(.data$clinical_cogdx) ~ "Missing cogdx",
      .data$clinical_cogdx %in% c(1, 2, 4) ~ "Primary cogdx 1/2/4",
      TRUE ~ "Non-primary cogdx"
    )
  )

# ------------------------------------------------------------
# 3. Independent check against Analysis_Meta_Merged
# ------------------------------------------------------------

analysis_meta <- readr::read_csv(
  analysis_meta_file,
  show_col_types = FALSE
) |>
  janitor::clean_names() |>
  standardize_common_aliases(
    context = "Analysis_Meta_Merged metadata",
    require_ids = TRUE,
    require_diagnosis = TRUE
  ) |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    analysis_meta_cogdx = safe_num(.data$diagnosis)
  ) |>
  dplyr::filter(
    !is.na(.data$individual_id),
    nzchar(.data$individual_id)
  )

analysis_meta_conflicts <- analysis_meta |>
  dplyr::filter(is.finite(.data$analysis_meta_cogdx)) |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_distinct_cogdx = dplyr::n_distinct(.data$analysis_meta_cogdx),
    cogdx_values = collapse_unique(.data$analysis_meta_cogdx),
    .groups = "drop"
  ) |>
  dplyr::filter(.data$n_distinct_cogdx > 1)

write_review_csv(
  analysis_meta_conflicts,
  "03_analysis_meta_cogdx_conflicts.csv"
)

analysis_meta_by_id <- analysis_meta |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    analysis_meta_cogdx = {
      vals <- unique(.data$analysis_meta_cogdx[is.finite(.data$analysis_meta_cogdx)])
      if (length(vals) == 0) NA_real_ else vals[[1]]
    },
    .groups = "drop"
  )

cogdx_source_agreement <- clinical_by_id |>
  dplyr::inner_join(analysis_meta_by_id, by = "individual_id") |>
  dplyr::mutate(
    comparable = is.finite(.data$clinical_cogdx) &
      is.finite(.data$analysis_meta_cogdx),
    agrees = dplyr::if_else(
      .data$comparable,
      .data$clinical_cogdx == .data$analysis_meta_cogdx,
      NA
    )
  )

write_review_csv(
  cogdx_source_agreement,
  "04_clinical_vs_analysis_meta_cogdx_agreement.csv"
)

# ------------------------------------------------------------
# 4. Load ALL matrix-backed proteomics samples WITHOUT applying
#    the legacy Control/AsymAD/AD filter.
# ------------------------------------------------------------

protein_mat_full <- read_gene_matrix_csv(
  protein_matrix_file,
  gene_col = 1
)

protein_meta_full <- readr::read_csv(
  protein_metadata_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  janitor::clean_names()

names(protein_meta_full) <- stringr::str_replace_all(
  names(protein_meta_full),
  "\\.",
  "_"
)

protein_meta_full <- protein_meta_full |>
  add_alias_column(
    "batch_channel",
    c("batch_channel", "batch.channel", "batch", "channel", "tmt_channel"),
    context = "full protein matched metadata"
  ) |>
  add_alias_column(
    "sample_id",
    alias_sets$sample_id,
    context = "full protein matched metadata"
  ) |>
  add_alias_column(
    "individual_id",
    alias_sets$individual_id,
    context = "full protein matched metadata"
  ) |>
  add_alias_column(
    "emory_strict_dx_2019",
    c(
      "emory_strict_dx_2019",
      "emorystrictdx_2019",
      "EmoryStrictDx.2019",
      "diagnosis",
      "cogdx"
    ),
    context = "full protein matched metadata"
  )

protein_meta_full <- protein_meta_full |>
  dplyr::mutate(
    batch_channel = as.character(.data$batch_channel),
    sample_id = as.character(.data$sample_id),
    individual_id = as.character(.data$individual_id),
    protein_legacy_dx = as.character(.data$emory_strict_dx_2019)
  ) |>
  dplyr::distinct(.data$batch_channel, .keep_all = TRUE)

matrix_channels <- colnames(protein_mat_full)
metadata_channels <- protein_meta_full$batch_channel
common_channels <- intersect(matrix_channels, metadata_channels)

if (length(common_channels) < 10) {
  stop(
    "Too few full TMT matrix columns matched the protein metadata batch/channel identifier.",
    call. = FALSE
  )
}

unmatched_matrix_channels <- setdiff(matrix_channels, metadata_channels)
unmatched_metadata_channels <- setdiff(metadata_channels, matrix_channels)

channel_alignment_audit <- tibble::tibble(
  metric = c(
    "Full protein matrix columns",
    "Unique metadata batch/channel IDs",
    "Matrix columns matched to metadata",
    "Matrix columns without metadata match",
    "Metadata channels without matrix match"
  ),
  value = c(
    length(matrix_channels),
    dplyr::n_distinct(metadata_channels),
    length(common_channels),
    length(unmatched_matrix_channels),
    length(unmatched_metadata_channels)
  )
)

write_review_csv(
  channel_alignment_audit,
  "05_full_protein_channel_alignment_summary.csv"
)

write_review_csv(
  tibble::tibble(batch_channel = unmatched_matrix_channels),
  "06_matrix_channels_without_metadata.csv"
)

write_review_csv(
  tibble::tibble(batch_channel = unmatched_metadata_channels),
  "07_metadata_channels_without_matrix.csv"
)

protein_meta_matrix_backed <- protein_meta_full[
  match(common_channels, protein_meta_full$batch_channel),
  ,
  drop = FALSE
]

if (!identical(
  as.character(protein_meta_matrix_backed$batch_channel),
  common_channels
)) {
  stop("Full protein metadata could not be aligned to matrix channels.", call. = FALSE)
}

# Detect region/tissue column if present. Do not infer or recode it.
region_candidates <- c(
  "region",
  "brain_region",
  "tissue",
  "tissue_region",
  "structure",
  "brainregion"
)
region_col <- pick_existing(protein_meta_matrix_backed, region_candidates)

protein_manifest <- tibble::as_tibble(protein_meta_matrix_backed) |>
  dplyr::transmute(
    batch_channel = as.character(.data$batch_channel),
    protein_sample_id = as.character(.data$sample_id),
    individual_id = as.character(.data$individual_id),
    protein_legacy_dx = as.character(.data$protein_legacy_dx),
    region_or_tissue = if (!is.na(region_col)) {
      as.character(.data[[region_col]])
    } else {
      NA_character_
    }
  ) |>
  dplyr::left_join(clinical_by_id, by = "individual_id") |>
  dplyr::left_join(analysis_meta_by_id, by = "individual_id") |>
  dplyr::mutate(
    legacy_pipeline_eligible = .data$protein_legacy_dx %in%
      c("Control", "AsymAD", "AD"),
    proposed_clinical_stage_eligible = .data$clinical_stage_eligible,
    cohort_transition = dplyr::case_when(
      .data$legacy_pipeline_eligible &
        .data$proposed_clinical_stage_eligible ~ "Retained by both definitions",
      !.data$legacy_pipeline_eligible &
        .data$proposed_clinical_stage_eligible ~ "Gained under cogdx definition",
      .data$legacy_pipeline_eligible &
        !.data$proposed_clinical_stage_eligible ~ "Lost under cogdx definition",
      TRUE ~ "Excluded by both definitions"
    )
  )

write_review_csv(
  protein_manifest,
  "08_full_matrix_backed_protein_manifest.csv"
)

# ------------------------------------------------------------
# 5. Duplicate participant audit
# ------------------------------------------------------------

participant_sample_counts <- protein_manifest |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_matrix_samples = dplyr::n(),
    protein_sample_ids = collapse_unique(.data$protein_sample_id),
    batch_channels = collapse_unique(.data$batch_channel),
    regions_or_tissues = collapse_unique(.data$region_or_tissue),
    clinical_cogdx = collapse_unique(.data$clinical_cogdx),
    clinical_stage = collapse_unique(.data$clinical_stage),
    legacy_dx_values = collapse_unique(.data$protein_legacy_dx),
    any_legacy_eligible = any(.data$legacy_pipeline_eligible),
    any_clinical_stage_eligible = any(.data$proposed_clinical_stage_eligible),
    .groups = "drop"
  )

duplicate_participants <- participant_sample_counts |>
  dplyr::filter(.data$n_matrix_samples > 1)

write_review_csv(
  duplicate_participants,
  "09_duplicate_participants_in_full_protein_matrix.csv"
)

# ------------------------------------------------------------
# 6. Cohort-definition summaries
# ------------------------------------------------------------

cohort_transition_counts <- protein_manifest |>
  dplyr::count(.data$cohort_transition, name = "n_samples") |>
  dplyr::arrange(dplyr::desc(.data$n_samples))

write_review_csv(
  cohort_transition_counts,
  "10_legacy_vs_cogdx_cohort_transition_counts.csv"
)

clinical_stage_counts <- protein_manifest |>
  dplyr::mutate(
    clinical_stage_display = dplyr::case_when(
      !is.na(.data$clinical_stage) ~ .data$clinical_stage,
      .data$clinical_cogdx_status == "Non-primary cogdx" ~ "Excluded non-primary cogdx",
      TRUE ~ "Missing cogdx"
    )
  ) |>
  dplyr::count(
    .data$clinical_stage_display,
    name = "n_samples"
  )

write_review_csv(
  clinical_stage_counts,
  "11_full_protein_cogdx_stage_counts.csv"
)

legacy_dx_counts <- protein_manifest |>
  dplyr::mutate(
    protein_legacy_dx = dplyr::if_else(
      is.na(.data$protein_legacy_dx) | !nzchar(.data$protein_legacy_dx),
      "Missing legacy diagnosis",
      .data$protein_legacy_dx
    )
  ) |>
  dplyr::count(
    .data$protein_legacy_dx,
    name = "n_samples",
    sort = TRUE
  )

write_review_csv(
  legacy_dx_counts,
  "12_full_protein_legacy_diagnosis_counts.csv"
)

full_crosswalk <- protein_manifest |>
  dplyr::mutate(
    clinical_stage_display = dplyr::case_when(
      !is.na(.data$clinical_stage) ~ .data$clinical_stage,
      .data$clinical_cogdx_status == "Non-primary cogdx" ~ "Excluded non-primary cogdx",
      TRUE ~ "Missing cogdx"
    ),
    legacy_dx_display = dplyr::if_else(
      is.na(.data$protein_legacy_dx) | !nzchar(.data$protein_legacy_dx),
      "Missing legacy diagnosis",
      .data$protein_legacy_dx
    )
  ) |>
  dplyr::count(
    .data$clinical_stage_display,
    .data$legacy_dx_display,
    name = "n_samples"
  ) |>
  dplyr::arrange(
    .data$clinical_stage_display,
    dplyr::desc(.data$n_samples)
  )

write_review_csv(
  full_crosswalk,
  "13_full_cogdx_vs_legacy_diagnosis_crosswalk.csv"
)

region_counts <- if (!is.na(region_col)) {
  protein_manifest |>
    dplyr::mutate(
      region_or_tissue = dplyr::if_else(
        is.na(.data$region_or_tissue) | !nzchar(.data$region_or_tissue),
        "Missing",
        .data$region_or_tissue
      )
    ) |>
    dplyr::count(
      .data$region_or_tissue,
      .data$proposed_clinical_stage_eligible,
      name = "n_samples"
    )
} else {
  tibble::tibble(
    region_or_tissue = "No region/tissue column detected in matched metadata",
    proposed_clinical_stage_eligible = NA,
    n_samples = nrow(protein_manifest)
  )
}

write_review_csv(region_counts, "14_region_tissue_audit.csv")

# ------------------------------------------------------------
# 7. Source agreement specifically among matrix-backed samples
# ------------------------------------------------------------

matrix_source_agreement <- protein_manifest |>
  dplyr::filter(
    is.finite(.data$clinical_cogdx),
    is.finite(.data$analysis_meta_cogdx)
  ) |>
  dplyr::mutate(
    cogdx_agrees = .data$clinical_cogdx == .data$analysis_meta_cogdx
  )

source_agreement_summary <- tibble::tibble(
  metric = c(
    "Matrix-backed samples with cogdx in both clinical and Analysis_Meta_Merged",
    "Agreement count",
    "Disagreement count",
    "Agreement fraction"
  ),
  value = c(
    nrow(matrix_source_agreement),
    sum(matrix_source_agreement$cogdx_agrees),
    sum(!matrix_source_agreement$cogdx_agrees),
    if (nrow(matrix_source_agreement) > 0) {
      mean(matrix_source_agreement$cogdx_agrees)
    } else {
      NA_real_
    }
  )
)

write_review_csv(
  source_agreement_summary,
  "15_matrix_backed_cogdx_source_agreement_summary.csv"
)

write_review_csv(
  matrix_source_agreement |>
    dplyr::filter(!.data$cogdx_agrees),
  "16_matrix_backed_cogdx_source_disagreements.csv"
)

# ------------------------------------------------------------
# 8. Compact audit summary
# ------------------------------------------------------------

current_legacy_n <- sum(protein_manifest$legacy_pipeline_eligible)
proposed_cogdx_n <- sum(
  protein_manifest$proposed_clinical_stage_eligible,
  na.rm = TRUE
)

summary_tbl <- tibble::tibble(
  metric = c(
    "Full protein matrix columns",
    "Matrix-backed samples with matched protein metadata",
    "Unique participants represented in full matrix-backed metadata",
    "Current pipeline legacy Control/AsymAD/AD eligible samples",
    "Proposed cogdx 1/2/4 eligible samples",
    "Samples retained by both definitions",
    "Samples gained under cogdx definition",
    "Samples lost under cogdx definition",
    "Samples excluded by both definitions",
    "Matrix-backed samples missing clinical cogdx",
    "Matrix-backed samples with non-primary cogdx",
    "Participants with >1 matrix-backed protein sample",
    "Detected region/tissue column"
  ),
  value = c(
    as.character(length(matrix_channels)),
    as.character(length(common_channels)),
    as.character(dplyr::n_distinct(protein_manifest$individual_id)),
    as.character(current_legacy_n),
    as.character(proposed_cogdx_n),
    as.character(sum(
      protein_manifest$cohort_transition == "Retained by both definitions"
    )),
    as.character(sum(
      protein_manifest$cohort_transition == "Gained under cogdx definition"
    )),
    as.character(sum(
      protein_manifest$cohort_transition == "Lost under cogdx definition"
    )),
    as.character(sum(
      protein_manifest$cohort_transition == "Excluded by both definitions"
    )),
    as.character(sum(
      protein_manifest$clinical_cogdx_status == "Missing cogdx"
    )),
    as.character(sum(
      protein_manifest$clinical_cogdx_status == "Non-primary cogdx"
    )),
    as.character(nrow(duplicate_participants)),
    if (is.na(region_col)) "NONE DETECTED" else region_col
  )
)

write_review_csv(
  summary_tbl,
  "00_full_proteomics_cogdx_eligibility_summary.csv"
)

# ------------------------------------------------------------
# 9. Provenance
# ------------------------------------------------------------

git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

git_branch <- tryCatch(
  system2("git", c("branch", "--show-current"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

git_status <- tryCatch(
  system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
  error = function(e) "git status unavailable"
)

provenance_lines <- c(
  "Reviewer revision audit 03 provenance",
  paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Repository root: ", normalizePath(cfg$project_dir, mustWork = FALSE)),
  paste0("Git branch: ", paste(git_branch, collapse = " ")),
  paste0("Git HEAD: ", paste(git_head, collapse = " ")),
  paste0("Working tree dirty: ", length(git_status) > 0),
  "",
  "Canonical repository scripts sourced:",
  paste0("  ", required_scripts),
  "",
  "Input paths resolved through repository configuration:",
  paste0("  ", source_paths$input, ": ", source_paths$path),
  "",
  "Important: this script intentionally bypasses the legacy",
  "Control/AsymAD/AD eligibility filter in R/02_load_data.R.",
  "It does not alter the manuscript pipeline.",
  "",
  "Git status --porcelain:",
  if (length(git_status) == 0) "  <clean>" else paste0("  ", git_status)
)

writeLines(
  provenance_lines,
  con = file.path(out_dir, "17_run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "18_sessionInfo.txt")
)

# ------------------------------------------------------------
# 10. Terminal output
# ------------------------------------------------------------

message("\n============================================================")
message("FULL PROTEOMICS COGDX ELIGIBILITY AUDIT")
message("============================================================\n")

message("SUMMARY")
print(summary_tbl, n = Inf)

message("\nCOHORT TRANSITION COUNTS")
print(cohort_transition_counts, n = Inf)

message("\nCOGDX STAGE COUNTS ACROSS ALL MATRIX-BACKED PROTEOMICS SAMPLES")
print(clinical_stage_counts, n = Inf)

message("\nLEGACY PROTEOMICS DIAGNOSIS COUNTS")
print(legacy_dx_counts, n = Inf)

message("\nCOGDX SOURCE AGREEMENT")
print(source_agreement_summary, n = Inf)

message("\nPrimary files to inspect next:")
message(file.path(out_dir, "00_full_proteomics_cogdx_eligibility_summary.csv"))
message(file.path(out_dir, "08_full_matrix_backed_protein_manifest.csv"))
message(file.path(out_dir, "10_legacy_vs_cogdx_cohort_transition_counts.csv"))
message(file.path(out_dir, "11_full_protein_cogdx_stage_counts.csv"))
message(file.path(out_dir, "13_full_cogdx_vs_legacy_diagnosis_crosswalk.csv"))
message(file.path(out_dir, "14_region_tissue_audit.csv"))
message("\nAudit complete.")
