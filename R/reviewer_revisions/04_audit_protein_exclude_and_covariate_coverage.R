############################################################
## Reviewer revision audit 04
## Clinical-stage cohort eligibility, legacy Exclude, and covariate coverage
##
## Purpose
## -------
## Audit the final remaining gate before changing the production pipeline:
##
## 1. Characterize participants labeled "Exclude" by the legacy
##    EmoryStrictDx.2019 proteomics grouping.
## 2. Determine whether cogdx-eligible participants currently labeled
##    "Exclude" have adequate protein data and nuisance covariates for
##    inclusion in the revised clinical-stage analysis.
## 3. Quantify covariate/pathology availability for the proposed
##    cogdx 1/2/4 proteomics cohort.
## 4. Surface any explicit QC/exclusion/reason fields in matched_metadata
##    rather than assuming what "Exclude" means.
##
## This is AUDIT-ONLY and does not modify the production pipeline.
##
## Run from repository root:
##   Rscript R/reviewer_revisions/04_audit_protein_exclude_and_covariate_coverage.R
############################################################

options(stringsAsFactors = FALSE)

required_scripts <- file.path("R", c("00_config.R", "01_utils.R"))
missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "Run this script from the repository root. Missing:\n",
    paste(missing_scripts, collapse = "\n"),
    call. = FALSE
  )
}

for (script_i in required_scripts) source(script_i)

require_objects("cfg", context = "reviewer audit 04")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "protein_exclude_covariate_audit"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_review_csv <- function(x, filename) {
  out <- file.path(out_dir, filename)
  readr::write_csv(x, out)
  message("Wrote: ", out)
  invisible(x)
}

canonical_batch_channel <- function(x) {
  x <- stringr::str_trim(as.character(x))
  hit <- stringr::str_match(
    x,
    stringr::regex("^b0*([0-9]+)\\.(.+)$", ignore_case = TRUE)
  )
  out <- x
  ok <- !is.na(hit[, 1])
  out[ok] <- paste0(
    "b",
    as.integer(hit[ok, 2]),
    ".",
    hit[ok, 3]
  )
  out
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

is_present_chr <- function(x) {
  !is.na(x) & nzchar(stringr::str_trim(as.character(x)))
}

# ------------------------------------------------------------
# 1. Resolve inputs
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

for (p in c(analysis_meta_file, protein_matrix_file, protein_metadata_file, clinical_file)) {
  if (!file.exists(p)) stop("Required input file not found: ", p, call. = FALSE)
}

write_review_csv(
  tibble::tibble(
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
  ),
  "01_input_source_paths.csv"
)

# ------------------------------------------------------------
# 2. Load full TMT matrix and matched metadata
# ------------------------------------------------------------

protein_mat <- read_gene_matrix_csv(protein_matrix_file, gene_col = 1)

meta_raw <- readr::read_csv(
  protein_metadata_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

meta <- meta_raw |>
  janitor::clean_names()

names(meta) <- stringr::str_replace_all(names(meta), "\\.", "_")

meta <- meta |>
  add_alias_column(
    "batch_channel",
    c("batch_channel", "batch.channel", "batch", "channel", "tmt_channel"),
    context = "protein matched metadata"
  ) |>
  add_alias_column(
    "individual_id",
    alias_sets$individual_id,
    context = "protein matched metadata"
  ) |>
  add_alias_column(
    "legacy_dx",
    c(
      "emory_strict_dx_2019",
      "emorystrictdx_2019",
      "EmoryStrictDx.2019",
      "diagnosis"
    ),
    context = "protein matched metadata"
  ) |>
  dplyr::mutate(
    batch_channel_raw = as.character(.data$batch_channel),
    batch_channel_canonical = canonical_batch_channel(.data$batch_channel),
    individual_id = as.character(.data$individual_id),
    legacy_dx = as.character(.data$legacy_dx)
  )

if (anyDuplicated(meta$batch_channel_canonical) > 0) {
  stop("Duplicate canonical batch/channel IDs in matched metadata.", call. = FALSE)
}

matrix_channels <- tibble::tibble(
  matrix_col_index = seq_along(colnames(protein_mat)),
  matrix_channel_raw = colnames(protein_mat),
  batch_channel_canonical = canonical_batch_channel(colnames(protein_mat))
)

if (anyDuplicated(matrix_channels$batch_channel_canonical) > 0) {
  stop("Duplicate canonical batch/channel IDs in protein matrix.", call. = FALSE)
}

channel_map <- matrix_channels |>
  dplyr::left_join(
    meta |>
      dplyr::select(
        "batch_channel_canonical",
        "batch_channel_raw",
        "individual_id",
        "legacy_dx",
        dplyr::everything()
      ),
    by = "batch_channel_canonical",
    suffix = c("_matrix", "_meta")
  )

if (any(is.na(channel_map$individual_id))) {
  stop("Not all 400 TMT matrix columns map to metadata after canonicalization.", call. = FALSE)
}

if (nrow(channel_map) != ncol(protein_mat)) {
  stop("Channel map is not one row per protein matrix column.", call. = FALSE)
}

# ------------------------------------------------------------
# 3. Canonical clinical data
# ------------------------------------------------------------

clinical_raw <- readr::read_csv(
  clinical_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

protein_ids <- unique(channel_map$individual_id)

clinical_id_candidates <- intersect(
  c(
    "individual_id", "individualid", "individualID", "IndividualID",
    "projid", "proj_id", "subject_id", "SubjectID"
  ),
  colnames(clinical_raw)
)

if (length(clinical_id_candidates) == 0) {
  stop("No recognized participant ID column in ROSMAP clinical file.", call. = FALSE)
}

clinical_id_overlap <- purrr::map_dfr(
  clinical_id_candidates,
  function(id_col) {
    ids <- as.character(clinical_raw[[id_col]])
    tibble::tibble(
      id_column = id_col,
      n_protein_ids_matched = length(
        intersect(protein_ids, unique(ids[is_present_chr(ids)]))
      )
    )
  }
) |>
  dplyr::arrange(dplyr::desc(.data$n_protein_ids_matched), .data$id_column)

write_review_csv(clinical_id_overlap, "02_clinical_id_overlap.csv")

clinical_id_col <- clinical_id_overlap$id_column[[1]]

if (clinical_id_overlap$n_protein_ids_matched[[1]] != length(protein_ids)) {
  stop(
    "Best clinical ID column does not match all protein participant IDs.",
    call. = FALSE
  )
}

cogdx_col <- require_alias_col(
  clinical_raw,
  c("cogdx", "COGDX", "diagnosis"),
  "cogdx",
  context = "ROSMAP clinical metadata"
)

pick_optional <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

clinical_age_col <- pick_optional(
  clinical_raw,
  c("age_death", "age", "age_at_death")
)
clinical_sex_col <- pick_optional(
  clinical_raw,
  c("msex", "sex", "gender")
)
clinical_pmi_col <- pick_optional(
  clinical_raw,
  c("pmi", "pmihours", "postmortem_interval")
)
clinical_braak_col <- pick_optional(
  clinical_raw,
  c("braaksc", "braak", "Braak")
)
clinical_cerad_col <- pick_optional(
  clinical_raw,
  c("ceradsc", "cerad", "CERAD")
)

clinical_tbl <- tibble::tibble(
  individual_id = as.character(clinical_raw[[clinical_id_col]]),
  cogdx = safe_num(clinical_raw[[cogdx_col]]),
  clinical_stage = clinical_stage_from_cogdx(clinical_raw[[cogdx_col]]),
  clinical_age = if (!is.na(clinical_age_col)) {
    safe_num(clinical_raw[[clinical_age_col]])
  } else {
    NA_real_
  },
  clinical_sex = if (!is.na(clinical_sex_col)) {
    as.character(clinical_raw[[clinical_sex_col]])
  } else {
    NA_character_
  },
  clinical_pmi = if (!is.na(clinical_pmi_col)) {
    safe_num(clinical_raw[[clinical_pmi_col]])
  } else {
    NA_real_
  },
  clinical_braak = if (!is.na(clinical_braak_col)) {
    safe_num(clinical_raw[[clinical_braak_col]])
  } else {
    NA_real_
  },
  clinical_cerad = if (!is.na(clinical_cerad_col)) {
    safe_num(clinical_raw[[clinical_cerad_col]])
  } else {
    NA_real_
  }
) |>
  dplyr::filter(is_present_chr(.data$individual_id))

if (anyDuplicated(clinical_tbl$individual_id) > 0) {
  stop(
    "Selected clinical participant ID is not unique. ",
    "Resolve before proceeding.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 4. Analysis_Meta_Merged covariates
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
  )

analysis_meta <- analysis_meta |>
  dplyr::mutate(
    individual_id = as.character(.data$individual_id),
    analysis_age = safe_num(
      coalesce_alias_chr(.data, alias_sets$age_death)
    ),
    analysis_pmi = safe_num(
      coalesce_alias_chr(.data, alias_sets$pmi)
    ),
    analysis_sex = coalesce_alias_chr(
      .data,
      c("sex_label", alias_sets$sex)
    ),
    analysis_braak = safe_num(
      coalesce_alias_chr(.data, alias_sets$braak)
    ),
    analysis_cerad = safe_num(
      coalesce_alias_chr(.data, alias_sets$cerad)
    )
  ) |>
  dplyr::select(
    "individual_id",
    "analysis_age",
    "analysis_sex",
    "analysis_pmi",
    "analysis_braak",
    "analysis_cerad"
  ) |>
  dplyr::distinct(.data$individual_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 5. Protein data completeness per matrix column
# ------------------------------------------------------------

finite_n <- colSums(is.finite(protein_mat))
finite_fraction <- finite_n / nrow(protein_mat)

protein_coverage <- tibble::tibble(
  matrix_col_index = seq_along(colnames(protein_mat)),
  matrix_channel_raw = colnames(protein_mat),
  n_finite_proteins = as.integer(finite_n),
  fraction_finite_proteins = as.numeric(finite_fraction)
)

# ------------------------------------------------------------
# 6. Unified 400-sample audit manifest
# ------------------------------------------------------------

manifest <- channel_map |>
  dplyr::select(
    "matrix_col_index",
    "matrix_channel_raw",
    "batch_channel_canonical",
    "individual_id",
    "legacy_dx"
  ) |>
  dplyr::left_join(clinical_tbl, by = "individual_id") |>
  dplyr::left_join(analysis_meta, by = "individual_id") |>
  dplyr::left_join(
    protein_coverage,
    by = c("matrix_col_index", "matrix_channel_raw")
  ) |>
  dplyr::mutate(
    cogdx_eligible = .data$cogdx %in% c(1, 2, 4),
    legacy_eligible = .data$legacy_dx %in% c("Control", "AsymAD", "AD"),
    cohort_transition = dplyr::case_when(
      .data$cogdx_eligible & .data$legacy_eligible ~
        "Retained by both",
      .data$cogdx_eligible & !.data$legacy_eligible ~
        "Gained under cogdx",
      !.data$cogdx_eligible & .data$legacy_eligible ~
        "Lost under cogdx",
      TRUE ~
        "Excluded by both"
    ),
    nuisance_complete_analysis_meta =
      is.finite(.data$analysis_age) &
      is_present_chr(.data$analysis_sex) &
      is.finite(.data$analysis_pmi),
    nuisance_complete_clinical =
      is.finite(.data$clinical_age) &
      is_present_chr(.data$clinical_sex) &
      is.finite(.data$clinical_pmi),
    pathology_complete_analysis_meta =
      is.finite(.data$analysis_braak) &
      is.finite(.data$analysis_cerad),
    pathology_complete_clinical =
      is.finite(.data$clinical_braak) &
      is.finite(.data$clinical_cerad)
  )

write_review_csv(manifest, "03_full_sample_covariate_manifest.csv")

# ------------------------------------------------------------
# 7. Characterize legacy "Exclude"
# ------------------------------------------------------------

exclude_manifest <- manifest |>
  dplyr::filter(.data$legacy_dx == "Exclude")

write_review_csv(exclude_manifest, "04_legacy_exclude_sample_manifest.csv")

exclude_summary <- exclude_manifest |>
  dplyr::count(
    .data$clinical_stage,
    .data$cogdx,
    .data$cogdx_eligible,
    name = "n_samples"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$cogdx_eligible),
    .data$clinical_stage,
    .data$cogdx
  )

write_review_csv(exclude_summary, "05_legacy_exclude_by_cogdx.csv")

# ------------------------------------------------------------
# 8. Look explicitly for QC / reason / exclusion metadata fields
# ------------------------------------------------------------

candidate_flag_regex <- stringr::regex(
  "exclude|reason|qc|quality|flag|outlier|contam|diagnos|dx|braak|cerad|cog|path|region|tissue",
  ignore_case = TRUE
)

flag_columns <- names(meta)[
  stringr::str_detect(names(meta), candidate_flag_regex)
]

flag_schema <- tibble::tibble(
  column = flag_columns,
  class = vapply(
    meta[flag_columns],
    function(x) paste(class(x), collapse = ";"),
    character(1)
  ),
  n_nonmissing = vapply(
    meta[flag_columns],
    function(x) sum(!is.na(x)),
    integer(1)
  ),
  n_unique = vapply(
    meta[flag_columns],
    function(x) dplyr::n_distinct(x, na.rm = TRUE),
    integer(1)
  )
)

write_review_csv(flag_schema, "06_candidate_qc_exclusion_columns.csv")

low_cardinality_flags <- flag_columns[
  vapply(
    meta[flag_columns],
    function(x) {
      n <- dplyr::n_distinct(x, na.rm = TRUE)
      is.finite(n) && n >= 1 && n <= 25
    },
    logical(1)
  )
]

flag_value_counts <- if (length(low_cardinality_flags) > 0) {
  purrr::map_dfr(
    low_cardinality_flags,
    function(col_i) {
      meta |>
        dplyr::transmute(
          legacy_dx = .data$legacy_dx,
          field = col_i,
          value = as.character(.data[[col_i]])
        ) |>
        dplyr::mutate(
          value = dplyr::if_else(
            is.na(.data$value) | !nzchar(.data$value),
            "<MISSING>",
            .data$value
          )
        ) |>
        dplyr::count(
          .data$legacy_dx,
          .data$field,
          .data$value,
          name = "n"
        )
    }
  )
} else {
  tibble::tibble(
    legacy_dx = character(),
    field = character(),
    value = character(),
    n = integer()
  )
}

write_review_csv(
  flag_value_counts,
  "07_candidate_qc_exclusion_value_counts.csv"
)

# ------------------------------------------------------------
# 9. Coverage and covariate comparisons
# ------------------------------------------------------------

coverage_by_transition <- manifest |>
  dplyr::group_by(.data$cohort_transition) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_finite_proteins = stats::median(.data$n_finite_proteins),
    min_finite_proteins = min(.data$n_finite_proteins),
    max_finite_proteins = max(.data$n_finite_proteins),
    median_fraction_finite = stats::median(.data$fraction_finite_proteins),
    nuisance_complete_analysis_meta =
      sum(.data$nuisance_complete_analysis_meta),
    nuisance_complete_clinical =
      sum(.data$nuisance_complete_clinical),
    braak_available_analysis_meta =
      sum(is.finite(.data$analysis_braak)),
    braak_available_clinical =
      sum(is.finite(.data$clinical_braak)),
    cerad_available_analysis_meta =
      sum(is.finite(.data$analysis_cerad)),
    cerad_available_clinical =
      sum(is.finite(.data$clinical_cerad)),
    .groups = "drop"
  )

write_review_csv(
  coverage_by_transition,
  "08_coverage_covariates_by_cohort_transition.csv"
)

coverage_by_stage <- manifest |>
  dplyr::filter(.data$cogdx_eligible) |>
  dplyr::group_by(.data$clinical_stage) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_finite_proteins = stats::median(.data$n_finite_proteins),
    median_fraction_finite = stats::median(.data$fraction_finite_proteins),
    nuisance_complete_analysis_meta =
      sum(.data$nuisance_complete_analysis_meta),
    nuisance_complete_clinical =
      sum(.data$nuisance_complete_clinical),
    braak_available_analysis_meta =
      sum(is.finite(.data$analysis_braak)),
    braak_available_clinical =
      sum(is.finite(.data$clinical_braak)),
    cerad_available_analysis_meta =
      sum(is.finite(.data$analysis_cerad)),
    cerad_available_clinical =
      sum(is.finite(.data$clinical_cerad)),
    .groups = "drop"
  )

write_review_csv(
  coverage_by_stage,
  "09_covariate_coverage_by_cogdx_stage.csv"
)

gained <- manifest |>
  dplyr::filter(.data$cohort_transition == "Gained under cogdx")

write_review_csv(
  gained,
  "10_gained_under_cogdx_manifest.csv"
)

gained_summary <- tibble::tibble(
  metric = c(
    "Gained samples",
    "Gained NCI",
    "Gained MCI",
    "Gained AD",
    "Median detected proteins in gained samples",
    "Minimum detected proteins in gained samples",
    "Gained with complete age/sex/PMI in Analysis_Meta_Merged",
    "Gained with complete age/sex/PMI in ROSMAP clinical",
    "Gained with Braak in Analysis_Meta_Merged",
    "Gained with Braak in ROSMAP clinical",
    "Gained with CERAD in Analysis_Meta_Merged",
    "Gained with CERAD in ROSMAP clinical"
  ),
  value = c(
    nrow(gained),
    sum(gained$clinical_stage == "NCI", na.rm = TRUE),
    sum(gained$clinical_stage == "MCI", na.rm = TRUE),
    sum(gained$clinical_stage == "AD", na.rm = TRUE),
    stats::median(gained$n_finite_proteins),
    min(gained$n_finite_proteins),
    sum(gained$nuisance_complete_analysis_meta),
    sum(gained$nuisance_complete_clinical),
    sum(is.finite(gained$analysis_braak)),
    sum(is.finite(gained$clinical_braak)),
    sum(is.finite(gained$analysis_cerad)),
    sum(is.finite(gained$clinical_cerad))
  )
)

write_review_csv(
  gained_summary,
  "00_gained_cogdx_and_covariate_summary.csv"
)

# ------------------------------------------------------------
# 10. Provenance
# ------------------------------------------------------------

git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

git_branch <- tryCatch(
  system2("git", c("branch", "--show-current"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

provenance <- c(
  "Reviewer revision audit 04 provenance",
  paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Repository root: ", normalizePath(cfg$project_dir, mustWork = FALSE)),
  paste0("Git branch: ", paste(git_branch, collapse = " ")),
  paste0("Git HEAD: ", paste(git_head, collapse = " ")),
  "",
  "Canonical scripts sourced:",
  paste0("  ", required_scripts),
  "",
  paste0("Selected clinical participant ID: ", clinical_id_col),
  paste0("Selected cogdx column: ", cogdx_col),
  "",
  "This script is audit-only and does not modify production objects."
)

writeLines(
  provenance,
  file.path(out_dir, "11_run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "12_sessionInfo.txt")
)

# ------------------------------------------------------------
# 11. Terminal summary
# ------------------------------------------------------------

message("\n============================================================")
message("PROTEIN EXCLUDE + COVARIATE COVERAGE AUDIT")
message("============================================================\n")

message("LEGACY EXCLUDE BY COGDX")
print(exclude_summary, n = Inf)

message("\nGAINED COGDX-ELIGIBLE SAMPLES")
print(gained_summary, n = Inf)

message("\nCOVERAGE / COVARIATES BY COHORT TRANSITION")
print(coverage_by_transition, n = Inf)

message("\nCANDIDATE QC/EXCLUSION FIELDS IN MATCHED METADATA")
print(flag_schema, n = Inf)

message("\nPrimary files to inspect:")
message(file.path(out_dir, "00_gained_cogdx_and_covariate_summary.csv"))
message(file.path(out_dir, "04_legacy_exclude_sample_manifest.csv"))
message(file.path(out_dir, "06_candidate_qc_exclusion_columns.csv"))
message(file.path(out_dir, "07_candidate_qc_exclusion_value_counts.csv"))
message(file.path(out_dir, "08_coverage_covariates_by_cohort_transition.csv"))
message(file.path(out_dir, "09_covariate_coverage_by_cogdx_stage.csv"))
message(file.path(out_dir, "10_gained_under_cogdx_manifest.csv"))
message("\nAudit complete.")
