############################################################
## 03_build_adjusted_core_objects.R
## Covariate-adjust RNA/protein matrices and build core metadata.
##
## Protein metadata now come directly from the correctly mapped full
## TMT cohort joined to ROSMAP_clinical by individualID.
## Analysis_Meta_Merged remains the RNA metadata/covariate source.
############################################################

require_objects(
  c(
    "rna_mat_raw",
    "prot_mat_raw",
    "rna_meta",
    "prot_meta_path",
    "analysis_meta"
  ),
  context = "03_build_adjusted_core_objects.R"
)

############################################################
## 1. RNA covariates
############################################################

rna_meta_adj <- rna_meta |>
  add_analysis_meta_covars(
    id_col = "individual_id",
    analysis_meta = analysis_meta
  ) |>
  prep_covariates()

require_columns(
  rna_meta_adj,
  c("sample_id", "clinical_stage", "sequencing_batch"),
  "RNA adjusted metadata"
)

if (anyNA(rna_meta_adj$sequencing_batch)) {
  stop(
    "RNA adjusted metadata contains missing sequencing_batch values.",
    call. = FALSE
  )
}

if (any(as.character(rna_meta_adj$sequencing_batch) == "0, 6, 7")) {
  stop(
    'RNA adjusted metadata contains excluded ambiguous batch "0, 6, 7".',
    call. = FALSE
  )
}

rna_meta_adj$sequencing_batch_factor <- factor(
  as.character(rna_meta_adj$sequencing_batch)
)

if (nlevels(rna_meta_adj$sequencing_batch_factor) != 9L) {
  stop(
    "Expected 9 RNA sequencing-batch levels; found ",
    nlevels(rna_meta_adj$sequencing_batch_factor),
    ".",
    call. = FALSE
  )
}

rna_meta_adj <- rna_meta_adj[
  match(colnames(rna_mat_raw), rna_meta_adj$sample_id),
  ,
  drop = FALSE
]

if (!all(colnames(rna_mat_raw) == rna_meta_adj$sample_id)) {
  stop(
    "RNA metadata could not be aligned to RNA matrix sample_id.",
    call. = FALSE
  )
}

############################################################
## 2. Protein covariates
############################################################

## Do NOT rejoin Analysis_Meta_Merged here. Protein clinical/pathology
## metadata were joined directly from ROSMAP_clinical in 02_load_data.R.
prot_meta_adj <- prot_meta_path |>
  prep_covariates()

require_columns(
  prot_meta_adj,
  c(
    "SampleID",
    "IndividualID",
    "cogdx_num",
    "clinical_stage",
    "clinical_stage_broad",
    "protein_legacy_dx",
    "tmt_batch",
    "batch"
  ),
  "protein adjusted metadata"
)

prot_meta_adj <- prot_meta_adj[
  match(colnames(prot_mat_raw), prot_meta_adj$SampleID),
  ,
  drop = FALSE
]

if (!all(colnames(prot_mat_raw) == prot_meta_adj$SampleID)) {
  stop(
    "Protein metadata could not be aligned to protein matrix SampleID.",
    call. = FALSE
  )
}

############################################################
## 3. Select nuisance covariates
############################################################

rna_covars <- select_usable_covars(
  rna_meta_adj,
  c(
    "age_num",
    "sex_factor",
    "pmi_num",
    "rin_num",
    "sequencing_batch_factor"
  )
)

required_rna_covars <- c(
  "age_num",
  "sex_factor",
  "pmi_num",
  "rin_num",
  "sequencing_batch_factor"
)

missing_required_rna_covars <- setdiff(
  required_rna_covars,
  rna_covars
)

if (length(missing_required_rna_covars) > 0L) {
  stop(
    "Required RNA nuisance covariates were not usable: ",
    paste(missing_required_rna_covars, collapse = ", "),
    call. = FALSE
  )
}

protein_covars <- select_usable_covars(
  prot_meta_adj,
  c(
    "age_num",
    "sex_factor",
    "pmi_num",
    "batch_factor"
  )
)

required_protein_covars <- c(
  "age_num",
  "sex_factor",
  "pmi_num",
  "batch_factor"
)

missing_protein_covars <- setdiff(
  required_protein_covars,
  protein_covars
)

if (length(missing_protein_covars) > 0) {
  stop(
    "Required protein nuisance covariates were not usable: ",
    paste(missing_protein_covars, collapse = ", "),
    ". Selected protein covariates: ",
    paste(protein_covars, collapse = ", "),
    call. = FALSE
  )
}

message(
  "RNA covariates used: ",
  paste(rna_covars, collapse = ", ")
)
message(
  "Protein covariates used: ",
  paste(protein_covars, collapse = ", ")
)

############################################################
## 4. Covariate adjustment
############################################################

if (isTRUE(cfg$adjust_main_figures)) {
  rna_mat_adj <- residualize_matrix(
    mat = rna_mat_raw,
    meta_df = rna_meta_adj,
    sample_col = "sample_id",
    covars = rna_covars
  )

  prot_mat_adj <- residualize_matrix(
    mat = prot_mat_raw,
    meta_df = prot_meta_adj,
    sample_col = "SampleID",
    covars = protein_covars
  )
} else {
  rna_mat_adj <- rna_mat_raw
  prot_mat_adj <- prot_mat_raw
}

## From this point onward, main figures use adjusted matrices by default.
rna_mat <- rna_mat_adj
prot_mat <- prot_mat_adj

############################################################
## 5. Explicit sample eligibility / completeness audits
############################################################

protein_nuisance_complete <- stats::complete.cases(
  prot_meta_adj[, protein_covars, drop = FALSE]
)

protein_has_adjusted_data <- colSums(is.finite(prot_mat_adj)) > 0

protein_primary_eligible <- !is.na(prot_meta_adj$clinical_stage)
protein_broad_eligible <- !is.na(prot_meta_adj$clinical_stage_broad)

protein_primary_analyzable <-
  protein_primary_eligible &
  protein_nuisance_complete &
  protein_has_adjusted_data

protein_broad_analyzable <-
  protein_broad_eligible &
  protein_nuisance_complete &
  protein_has_adjusted_data

protein_primary_analysis_counts <- tibble::tibble(
  stage = factor(
    as.character(prot_meta_adj$clinical_stage),
    levels = c("NCI", "MCI", "AD")
  ),
  phenotype_eligible = protein_primary_eligible,
  nuisance_complete = protein_nuisance_complete,
  adjusted_data_available = protein_has_adjusted_data,
  analyzable = protein_primary_analyzable
) |>
  dplyr::filter(.data$phenotype_eligible) |>
  dplyr::group_by(.data$stage, .drop = FALSE) |>
  dplyr::summarise(
    n_phenotype_eligible = dplyr::n(),
    n_nuisance_complete = sum(.data$nuisance_complete),
    n_adjusted_data_available = sum(.data$adjusted_data_available),
    n_analyzable = sum(.data$analyzable),
    .groups = "drop"
  )

protein_broad_analysis_counts <- tibble::tibble(
  stage = factor(
    as.character(prot_meta_adj$clinical_stage_broad),
    levels = c("NCI", "MCI", "AD")
  ),
  phenotype_eligible = protein_broad_eligible,
  nuisance_complete = protein_nuisance_complete,
  adjusted_data_available = protein_has_adjusted_data,
  analyzable = protein_broad_analyzable
) |>
  dplyr::filter(.data$phenotype_eligible) |>
  dplyr::group_by(.data$stage, .drop = FALSE) |>
  dplyr::summarise(
    n_phenotype_eligible = dplyr::n(),
    n_nuisance_complete = sum(.data$nuisance_complete),
    n_adjusted_data_available = sum(.data$adjusted_data_available),
    n_analyzable = sum(.data$analyzable),
    .groups = "drop"
  )

adjustment_audit <- tibble::tibble(
  modality = c("RNA", "Protein"),
  original_matrix = c("rna_mat_raw", "prot_mat_raw"),
  adjusted_matrix = c("rna_mat_adj", "prot_mat_adj"),
  n_input_samples = c(
    ncol(rna_mat_raw),
    ncol(prot_mat_raw)
  ),
  n_samples_with_adjusted_data = c(
    sum(colSums(is.finite(rna_mat_adj)) > 0),
    sum(protein_has_adjusted_data)
  ),
  covariates_used = c(
    paste(rna_covars, collapse = ";"),
    paste(protein_covars, collapse = ";")
  ),
  biological_variables_not_adjusted_here = c(
    "clinical_stage;braak_num;cerad_num",
    paste(
      c(
        "clinical_stage",
        "clinical_stage_broad",
        "cogdx_num",
        "braak_num",
        "cerad_num"
      ),
      collapse = ";"
    )
  ),
  note = c(
    "Clinical stage/pathology variables are tested downstream, not regressed out during nuisance residualization.",
    "Protein adjustment uses ROSMAP_clinical age/sex/PMI plus TMT batch derived from canonical batch-channel IDs."
  )
)

write_tbl(
  adjustment_audit,
  "covariate_adjustment_audit"
)
write_tbl(
  protein_primary_analysis_counts,
  "protein_primary_stage_analysis_counts"
)
write_tbl(
  protein_broad_analysis_counts,
  "protein_broad_stage_analysis_counts"
)

core_adjustment_checkpoint <- tibble::tibble(
  metric = c(
    "Protein input samples",
    "Protein samples with adjusted abundance data",
    "Protein primary-stage phenotype eligible",
    "Protein primary-stage analyzable",
    "Protein broad-stage phenotype eligible",
    "Protein broad-stage analyzable",
    "Protein TMT batch levels"
  ),
  value = c(
    ncol(prot_mat_raw),
    sum(protein_has_adjusted_data),
    sum(protein_primary_eligible),
    sum(protein_primary_analyzable),
    sum(protein_broad_eligible),
    sum(protein_broad_analyzable),
    dplyr::n_distinct(prot_meta_adj$tmt_batch)
  )
)

write_tbl(
  core_adjustment_checkpoint,
  "core_adjustment_checkpoint"
)

save_obj(rna_meta_adj, "rna_meta_adj")
save_obj(prot_meta_adj, "prot_meta_adj")
save_obj(rna_covars, "rna_covars")
save_obj(protein_covars, "protein_covars")
save_obj(rna_mat_adj, "rna_mat_adj")
save_obj(prot_mat_adj, "prot_mat_adj")

message("Loaded 03_build_adjusted_core_objects.R")
