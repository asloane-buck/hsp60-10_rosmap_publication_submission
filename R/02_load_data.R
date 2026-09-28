############################################################
## 02_load_data.R
## Load RNA, proteomics, and metadata.
##
## Revised clinical-stage / proteomics ingestion architecture:
## - RNA primary stage: ROSMAP cogdx 1/2/4 -> NCI/MCI/AD.
## - Protein source: all matrix-backed TMT samples, with no legacy
##   Control/AsymAD/AD eligibility filter at ingestion.
## - Protein clinical metadata: ROSMAP_clinical.csv joined by individualID.
## - TMT batch/channel IDs are canonicalized before matrix/metadata matching.
## - EmoryStrictDx.2019 is retained only as protein_legacy_dx.
############################################################

require_objects("cfg", context = "02_load_data.R")

############################################################
## 1. Shared / RNA metadata
############################################################

load_analysis_meta <- function() {
  meta <- readr::read_csv(
    first_existing(
      c(
        file.path(cfg$metadata_dir, "Analysis_Meta_Merged.csv"),
        file.path(cfg$metadata_dir, "Analysis_Meta_Merged_copy.csv"),
        file.path(cfg$derived_dir, "Analysis_Meta_Merged_copy.csv")
      ),
      "Analysis_Meta_Merged metadata"
    ),
    show_col_types = FALSE
  ) |>
    janitor::clean_names()

  meta <- standardize_common_aliases(
    meta,
    context = "Analysis_Meta_Merged metadata",
    require_ids = TRUE,
    require_diagnosis = TRUE
  )

  meta |>
    dplyr::mutate(
      sample_id = as.character(.data$sample_id),
      individual_id = as.character(.data$individual_id),
      cogdx_num = safe_num(.data$diagnosis),
      clinical_stage = clinical_stage_primary(.data$diagnosis),
      clinical_stage_broad = clinical_stage_broad(.data$diagnosis),
      clinical_stage = factor(
        .data$clinical_stage,
        levels = c("NCI", "MCI", "AD")
      ),
      clinical_stage_broad = factor(
        .data$clinical_stage_broad,
        levels = c("NCI", "MCI", "AD")
      ),
      ## Transitional alias. This now carries the SAME NCI/MCI/AD labels
      ## as clinical_stage; it is not the old Control/Early_AD/AD coding.
      diagnosis_stage = .data$clinical_stage,
      braak_num = safe_num(.data$braak),
      cerad_num = safe_num(.data$cerad),
      sex_label = dplyr::case_when(
        .data$sex %in% c(0, "0", "Female", "female", "F", "f") ~ "Female",
        .data$sex %in% c(1, "1", "Male", "male", "M", "m") ~ "Male",
        TRUE ~ as.character(.data$sex)
      )
    )
}

load_rna_matrix <- function() {
  rna_file <- first_existing(
    c(
      file.path(cfg$derived_dir, "ROSMAP_vst_gene_symbol_matrix_STAGE.csv"),
      file.path(
        cfg$project_dir,
        "ROSMAP_pathway_scores_robustness",
        "ROSMAP_vst_gene_symbol_matrix_STAGE.csv"
      )
    ),
    "ROSMAP VST gene-symbol RNA matrix"
  )

  read_gene_matrix_csv(rna_file, gene_col = "gene_symbol")
}

############################################################
## 2. Canonical ROSMAP clinical metadata for proteomics
############################################################

load_rosmap_clinical_core <- function() {
  if (!file.exists(cfg$rosmap_clinical_file)) {
    stop(
      "Missing ROSMAP clinical file: ",
      cfg$rosmap_clinical_file,
      call. = FALSE
    )
  }

  clinical_raw <- readr::read_csv(
    cfg$rosmap_clinical_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  ## The proteomics IndividualID values were empirically shown
  ## to match ROSMAP_clinical$individualID 400/400. Do not substitute projid.
  if (!"individualID" %in% colnames(clinical_raw)) {
    stop(
      "ROSMAP_clinical.csv is missing required proteomics join key ",
      "'individualID'. Available columns: ",
      paste(colnames(clinical_raw), collapse = ", "),
      call. = FALSE
    )
  }

  cogdx_col <- require_alias_col(
    clinical_raw,
    c("cogdx", "COGDX"),
    "cogdx",
    context = "ROSMAP clinical metadata"
  )
  age_col <- require_alias_col(
    clinical_raw,
    alias_sets$age_death,
    "age at death",
    context = "ROSMAP clinical metadata"
  )
  sex_col <- require_alias_col(
    clinical_raw,
    alias_sets$sex,
    "sex",
    context = "ROSMAP clinical metadata"
  )
  pmi_col <- require_alias_col(
    clinical_raw,
    alias_sets$pmi,
    "PMI",
    context = "ROSMAP clinical metadata"
  )
  braak_col <- require_alias_col(
    clinical_raw,
    alias_sets$braak,
    "Braak",
    context = "ROSMAP clinical metadata"
  )
  cerad_col <- require_alias_col(
    clinical_raw,
    alias_sets$cerad,
    "CERAD",
    context = "ROSMAP clinical metadata"
  )

  sex_raw <- clinical_raw[[sex_col]]

  out <- tibble::tibble(
    IndividualID = as.character(clinical_raw$individualID),
    cogdx_num = safe_num(clinical_raw[[cogdx_col]]),
    clinical_stage = clinical_stage_primary(clinical_raw[[cogdx_col]]),
    clinical_stage_broad = clinical_stage_broad(clinical_raw[[cogdx_col]]),
    age_death = safe_num(clinical_raw[[age_col]]),
    sex_label = dplyr::case_when(
      sex_raw %in% c(0, "0", "Female", "female", "F", "f") ~ "Female",
      sex_raw %in% c(1, "1", "Male", "male", "M", "m") ~ "Male",
      TRUE ~ as.character(sex_raw)
    ),
    pmi = safe_num(clinical_raw[[pmi_col]]),
    braak = safe_num(clinical_raw[[braak_col]]),
    cerad = safe_num(clinical_raw[[cerad_col]])
  ) |>
    dplyr::filter(
      !is.na(.data$IndividualID),
      nzchar(.data$IndividualID)
    ) |>
    dplyr::mutate(
      clinical_stage = factor(
        .data$clinical_stage,
        levels = c("NCI", "MCI", "AD")
      ),
      clinical_stage_broad = factor(
        .data$clinical_stage_broad,
        levels = c("NCI", "MCI", "AD")
      ),
      braak_num = safe_num(.data$braak),
      cerad_num = safe_num(.data$cerad)
    )

  if (anyDuplicated(out$IndividualID) > 0) {
    dup <- unique(out$IndividualID[duplicated(out$IndividualID)])
    stop(
      "ROSMAP clinical individualID is not unique. Duplicate IDs: ",
      paste(head(dup, 10), collapse = ", "),
      call. = FALSE
    )
  }

  out
}

############################################################
## 3. Full proteomics ingestion
############################################################

load_protein <- function(rosmap_clinical_core) {
  prot_file <- first_existing(
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
    "ROSMAP TMT protein matrix"
  )

  meta_file <- first_existing(
    c(
      file.path(cfg$derived_dir, "matched_metadata_copy.csv"),
      file.path(cfg$proteomics_input_dir, "matched_metadata.csv"),
      file.path(cfg$proteomics_input_dir, "matched_metadata_copy.csv")
    ),
    "ROSMAP TMT matched metadata"
  )

  mat <- read_gene_matrix_csv(prot_file, gene_col = 1)

  expected_n <- suppressWarnings(
    as.integer(Sys.getenv("ROSMAP_TMT_EXPECTED_N", unset = "400"))
  )
  if (is.finite(expected_n) && ncol(mat) != expected_n) {
    stop(
      "Unexpected number of TMT matrix columns. Expected ",
      expected_n,
      " for the manuscript dataset but found ",
      ncol(mat),
      ". Set ROSMAP_TMT_EXPECTED_N only if intentionally using a ",
      "different data release.",
      call. = FALSE
    )
  }

  meta <- readr::read_csv(
    meta_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  ) |>
    janitor::clean_names()

  names(meta) <- stringr::str_replace_all(names(meta), "\\.", "_")

  meta <- add_alias_column(
    meta,
    "batch_channel",
    c("batch_channel", "batch.channel", "channel", "tmt_channel"),
    context = "protein matched metadata"
  )
  meta <- add_alias_column(
    meta,
    "sample_id",
    alias_sets$sample_id,
    context = "protein matched metadata"
  )
  meta <- add_alias_column(
    meta,
    "individual_id",
    alias_sets$individual_id,
    context = "protein matched metadata"
  )
  meta <- add_alias_column(
    meta,
    "emory_strict_dx_2019",
    c(
      "emory_strict_dx_2019",
      "emorystrictdx_2019",
      "EmoryStrictDx.2019"
    ),
    context = "protein matched metadata"
  )

  require_columns(
    meta,
    c(
      "batch_channel",
      "sample_id",
      "individual_id",
      "emory_strict_dx_2019"
    ),
    "protein matched metadata"
  )

  meta <- meta |>
    dplyr::mutate(
      batch_channel_raw = as.character(.data$batch_channel),
      batch_channel_canonical = canonical_batch_channel(.data$batch_channel),
      SampleID = as.character(.data$sample_id),
      IndividualID = as.character(.data$individual_id),
      protein_legacy_dx = as.character(.data$emory_strict_dx_2019),
      EmoryStrictDx.2019 = .data$protein_legacy_dx,
      tmt_batch = protein_batch_from_channel(.data$batch_channel_canonical),
      ## prep_covariates() preferentially looks for "batch".
      ## This is the TMT plex/batch, NOT the unique batch-channel identifier.
      batch = .data$tmt_batch
    )

  if (anyDuplicated(meta$batch_channel_canonical) > 0) {
    stop(
      "Canonical protein metadata batch/channel IDs are not unique.",
      call. = FALSE
    )
  }

  matrix_channel_raw <- colnames(mat)
  matrix_channel_canonical <- canonical_batch_channel(matrix_channel_raw)

  if (anyDuplicated(matrix_channel_canonical) > 0) {
    stop(
      "Canonical protein matrix batch/channel IDs are not unique.",
      call. = FALSE
    )
  }

  idx <- match(
    matrix_channel_canonical,
    meta$batch_channel_canonical
  )

  if (any(is.na(idx))) {
    missing_channels <- matrix_channel_raw[is.na(idx)]
    stop(
      "Protein matrix columns remain unmatched after canonical ",
      "batch/channel normalization: ",
      paste(head(missing_channels, 20), collapse = ", "),
      call. = FALSE
    )
  }

  meta <- meta[idx, , drop = FALSE]

  if (nrow(meta) != ncol(mat)) {
    stop(
      "Protein metadata does not map one-to-one to all matrix columns.",
      call. = FALSE
    )
  }

  if (anyDuplicated(meta$SampleID) > 0) {
    stop("Protein SampleID values are not unique.", call. = FALSE)
  }

  if (anyDuplicated(meta$IndividualID) > 0) {
    stop(
      "More than one matrix-backed protein sample exists for at least ",
      "one participant. Resolve before continuing.",
      call. = FALSE
    )
  }

  meta$batch_channel_matrix_raw <- matrix_channel_raw

  ## Rename matrix columns only AFTER verified one-to-one channel mapping.
  colnames(mat) <- meta$SampleID

  meta <- checked_left_join(
    meta,
    rosmap_clinical_core,
    by = "IndividualID",
    label = "protein metadata to ROSMAP clinical"
  )

  if (any(is.na(meta$cogdx_num))) {
    stop(
      "At least one matrix-backed protein participant is missing cogdx ",
      "after the individualID clinical join.",
      call. = FALSE
    )
  }

  if (!all(meta$SampleID == colnames(mat))) {
    stop(
      "Protein metadata/sample ordering changed during clinical join.",
      call. = FALSE
    )
  }

  list(
    mat = mat,
    meta = meta,
    matrix_channel_raw = matrix_channel_raw,
    matrix_channel_canonical = matrix_channel_canonical
  )
}

############################################################
## 4. Build production objects
############################################################

analysis_meta <- load_analysis_meta()
rosmap_clinical_core <- load_rosmap_clinical_core()

rna_mat_raw <- load_rna_matrix()

protein <- load_protein(rosmap_clinical_core)
prot_mat_raw <- protein$mat
prot_meta_path <- protein$meta

## RNA primary clinical-stage cohort remains the narrow cogdx 1/2/4
##
## RNA technical-QC exclusion:
## 492_120515 has the ambiguous composite sequencing-batch annotation
## "0, 6, 7" and cannot be assigned to one technical batch.
rna_excluded_ambiguous_batch_samples <- c("492_120515")

if (!all(rna_excluded_ambiguous_batch_samples %in% analysis_meta$sample_id)) {
  stop(
    "Expected ambiguous-batch RNA sample was not present in analysis_meta.",
    call. = FALSE
  )
}

## cohort used in the submitted RNA analysis.
rna_meta <- analysis_meta |>
  dplyr::filter(
    .data$sample_id %in% colnames(rna_mat_raw),
    !is.na(.data$clinical_stage),
    !.data$sample_id %in% rna_excluded_ambiguous_batch_samples
  ) |>
  dplyr::distinct(.data$sample_id, .keep_all = TRUE)

rna_mat_raw <- rna_mat_raw[, rna_meta$sample_id, drop = FALSE]
rna_meta <- rna_meta[
  match(colnames(rna_mat_raw), rna_meta$sample_id),
  ,
  drop = FALSE
]

if (!all(colnames(rna_mat_raw) == rna_meta$sample_id)) {
  stop(
    "RNA metadata could not be aligned to RNA matrix sample_id.",
    call. = FALSE
  )
}

## Hard validation of canonical RNA technical-QC cohort.
if ("492_120515" %in% rna_meta$sample_id ||
    "492_120515" %in% colnames(rna_mat_raw)) {
  stop(
    "Ambiguous sequencing-batch sample 492_120515 remains in RNA cohort.",
    call. = FALSE
  )
}

if (!"sequencing_batch" %in% colnames(rna_meta)) {
  stop(
    "RNA metadata is missing sequencing_batch.",
    call. = FALSE
  )
}

if (anyNA(rna_meta$sequencing_batch)) {
  stop(
    "Canonical RNA cohort contains missing sequencing_batch values.",
    call. = FALSE
  )
}

if (any(as.character(rna_meta$sequencing_batch) == "0, 6, 7")) {
  stop(
    'Canonical RNA cohort still contains ambiguous batch "0, 6, 7".',
    call. = FALSE
  )
}

if (ncol(rna_mat_raw) != 577L || nrow(rna_meta) != 577L) {
  stop(
    "Expected canonical RNA cohort n=577; found matrix n=",
    ncol(rna_mat_raw),
    " and metadata n=",
    nrow(rna_meta),
    ".",
    call. = FALSE
  )
}

expected_rna_stage_counts <- c(
  NCI = 200L,
  MCI = 158L,
  AD = 219L
)

observed_rna_stage_counts <- vapply(
  names(expected_rna_stage_counts),
  function(stage_name) {
    sum(as.character(rna_meta$clinical_stage) == stage_name)
  },
  integer(1)
)

if (!identical(
  unname(observed_rna_stage_counts),
  unname(expected_rna_stage_counts)
)) {
  stop(
    "Unexpected canonical RNA stage counts: ",
    paste(
      names(observed_rna_stage_counts),
      observed_rna_stage_counts,
      sep = "=",
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (dplyr::n_distinct(rna_meta$sequencing_batch) != 9L) {
  stop(
    "Expected 9 valid RNA sequencing batches after QC exclusion; found ",
    dplyr::n_distinct(rna_meta$sequencing_batch),
    ".",
    call. = FALSE
  )
}

## Backward-compatible matrix names. No stage filtering is applied to protein.
rna_mat <- rna_mat_raw
prot_mat <- prot_mat_raw

############################################################
## 5. Aggregate ingestion audits
############################################################

protein_primary_counts <- prot_meta_path |>
  dplyr::filter(!is.na(.data$clinical_stage)) |>
  dplyr::count(.data$clinical_stage, name = "n") |>
  dplyr::rename(stage = clinical_stage)

protein_broad_counts <- prot_meta_path |>
  dplyr::filter(!is.na(.data$clinical_stage_broad)) |>
  dplyr::count(.data$clinical_stage_broad, name = "n") |>
  dplyr::rename(stage = clinical_stage_broad)

rna_primary_counts <- rna_meta |>
  dplyr::count(.data$clinical_stage, name = "n") |>
  dplyr::rename(stage = clinical_stage)

protein_batch_counts <- prot_meta_path |>
  dplyr::count(.data$tmt_batch, name = "n_samples") |>
  dplyr::arrange(.data$tmt_batch)

ingestion_audit <- tibble::tibble(
  metric = c(
    "RNA primary-stage samples",
    "Protein TMT matrix columns",
    "Protein metadata rows mapped to matrix",
    "Unique protein participants",
    "Protein participants matched to ROSMAP clinical",
    "Protein primary cogdx 1/2/4 eligible",
    "Protein broad cogdx 1/(2+3)/(4+5) eligible",
    "Protein legacy Control/AsymAD/AD",
    "Protein legacy Exclude",
    "TMT batches represented",
    "TMT columns requiring batch-number zero-padding normalization"
  ),
  value = c(
    ncol(rna_mat_raw),
    ncol(prot_mat_raw),
    nrow(prot_meta_path),
    dplyr::n_distinct(prot_meta_path$IndividualID),
    sum(!is.na(prot_meta_path$cogdx_num)),
    sum(!is.na(prot_meta_path$clinical_stage)),
    sum(!is.na(prot_meta_path$clinical_stage_broad)),
    sum(
      prot_meta_path$protein_legacy_dx %in% c("Control", "AsymAD", "AD"),
      na.rm = TRUE
    ),
    sum(prot_meta_path$protein_legacy_dx == "Exclude", na.rm = TRUE),
    dplyr::n_distinct(prot_meta_path$tmt_batch),
    sum(
      prot_meta_path$batch_channel_matrix_raw !=
        prot_meta_path$batch_channel_canonical
    )
  )
)

write_tbl(ingestion_audit, "protein_ingestion_audit")
write_tbl(protein_primary_counts, "protein_primary_cogdx_stage_counts")
write_tbl(protein_broad_counts, "protein_broad_cogdx_stage_counts")
write_tbl(rna_primary_counts, "rna_primary_cogdx_stage_counts")
write_tbl(protein_batch_counts, "protein_tmt_batch_counts")

write_tbl(
  tibble::tibble(
    object = c(
      "rna_mat_raw",
      "prot_mat_raw",
      "rna_meta",
      "prot_meta_path"
    ),
    nrow = c(
      nrow(rna_mat_raw),
      nrow(prot_mat_raw),
      nrow(rna_meta),
      nrow(prot_meta_path)
    ),
    ncol = c(
      ncol(rna_mat_raw),
      ncol(prot_mat_raw),
      ncol(rna_meta),
      ncol(prot_meta_path)
    )
  ),
  "load_data_dimensions"
)

save_obj(analysis_meta, "analysis_meta")
save_obj(rosmap_clinical_core, "rosmap_clinical_core")
save_obj(rna_mat_raw, "rna_mat_raw")
save_obj(prot_mat_raw, "prot_mat_raw")
save_obj(rna_meta, "rna_meta")
save_obj(prot_meta_path, "prot_meta_path")

message("Loaded 02_load_data.R")
