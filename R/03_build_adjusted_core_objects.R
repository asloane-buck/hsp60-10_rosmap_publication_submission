############################################################
## 03_build_adjusted_core_objects.R
## Covariate-adjust RNA/protein matrices and build core metadata.
############################################################

require_objects(
  c("rna_mat_raw", "prot_mat_raw", "rna_meta", "prot_meta_path", "analysis_meta"),
  context = "03_build_adjusted_core_objects.R"
)

rna_meta_adj <- rna_meta |>
  add_analysis_meta_covars(id_col = "individual_id", analysis_meta = analysis_meta) |>
  prep_covariates()

require_columns(rna_meta_adj, "sample_id", "RNA adjusted metadata")
rna_meta_adj <- rna_meta_adj[match(colnames(rna_mat_raw), rna_meta_adj$sample_id), , drop = FALSE]
if (!all(colnames(rna_mat_raw) == rna_meta_adj$sample_id)) {
  stop(
    "RNA metadata could not be aligned to RNA matrix sample_id. Available columns: ",
    available_cols_msg(rna_meta_adj),
    call. = FALSE
  )
}

prot_meta_adj <- prot_meta_path |>
  add_analysis_meta_covars(id_col = "IndividualID", analysis_meta = analysis_meta) |>
  prep_covariates()

require_columns(prot_meta_adj, "SampleID", "protein adjusted metadata")
prot_meta_adj <- prot_meta_adj[match(colnames(prot_mat_raw), prot_meta_adj$SampleID), , drop = FALSE]
if (!all(colnames(prot_mat_raw) == prot_meta_adj$SampleID)) {
  stop(
    "Protein metadata could not be aligned to protein matrix SampleID. Available columns: ",
    available_cols_msg(prot_meta_adj),
    call. = FALSE
  )
}

rna_covars <- select_usable_covars(
  rna_meta_adj,
  c("age_num", "sex_factor", "pmi_num", "rin_num", "batch_factor")
)

protein_covars <- select_usable_covars(
  prot_meta_adj,
  c("age_num", "sex_factor", "pmi_num", "batch_factor")
)

message("RNA covariates used: ", paste(rna_covars, collapse = ", "))
message("Protein covariates used: ", paste(protein_covars, collapse = ", "))

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

adjustment_audit <- tibble(
  modality = c("RNA", "Protein"),
  original_matrix = c("rna_mat_raw", "prot_mat_raw"),
  adjusted_matrix = c("rna_mat_adj", "prot_mat_adj"),
  covariates_used = c(paste(rna_covars, collapse = ";"), paste(protein_covars, collapse = ";")),
  biological_variables_not_adjusted_here = c("diagnosis_stage;braak_num;cerad_num", "EmoryStrictDx.2019;braak_num;cerad_num"),
  note = "Stage/pathology variables are tested downstream, not regressed out during nuisance residualization."
)

write_tbl(adjustment_audit, "covariate_adjustment_audit")
save_obj(rna_meta_adj, "rna_meta_adj")
save_obj(prot_meta_adj, "prot_meta_adj")
save_obj(rna_covars, "rna_covars")
save_obj(protein_covars, "protein_covars")
save_obj(rna_mat_adj, "rna_mat_adj")
save_obj(prot_mat_adj, "prot_mat_adj")

message("Loaded 03_build_adjusted_core_objects.R")
