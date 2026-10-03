############################################################
## 06_build_cognition_objects.R
## Build the cognition-model source table from the corrected full
## proteomics cohort.
##
## Key invariants:
## - Hsp60/10 abundance score is recomputed from processed non-residualized prot_mat_raw.
## - Protein sample/participant mapping comes from the verified 400-sample
##   TMT ingestion.
## - Cognition outcomes/covariates come directly from ROSMAP_clinical.
## - No EmoryStrictDx/AsymAD classification is used here.
############################################################

require_objects(
  c(
    "cfg",
    "prot_mat_raw",
    "prot_meta_adj",
    "prot_scores",
    "prot_meta_path",
    "all_hsp60_10_client_tbl"
  ),
  context = "06_build_cognition_objects.R"
)

require_columns(
  prot_scores,
  c(
    "SampleID",
    "IndividualID",
    "Hsp60_10_all_clients",
    "cogdx_num",
    "clinical_stage",
    "clinical_stage_broad"
  ),
  "protein pathway-score table"
)

if (!file.exists(cfg$rosmap_clinical_file)) {
  stop(
    "Missing ROSMAP clinical file required for cognition analysis: ",
    cfg$rosmap_clinical_file,
    call. = FALSE
  )
}

############################################################
## 1. Verified Hsp60/10 client network score
############################################################

n_total_clients_reference <- nrow(all_hsp60_10_client_tbl)
n_detected_clients_used <- sum(
  all_hsp60_10_client_tbl$detected_in_protein,
  na.rm = TRUE
)

protein_mat_raw_cog <- as.matrix(prot_mat_raw)
mode(protein_mat_raw_cog) <- "numeric"

canonical_protein_rows <- canonical_gene_symbol(
  rownames(protein_mat_raw_cog)
)

if (anyDuplicated(canonical_protein_rows) > 0) {
  duplicated_genes <- unique(
    canonical_protein_rows[duplicated(canonical_protein_rows)]
  )
  stop(
    "Canonical gene-symbol conversion created duplicate protein rows in ",
    "cognition builder: ",
    paste(head(duplicated_genes, 20), collapse = ", "),
    call. = FALSE
  )
}

rownames(protein_mat_raw_cog) <- canonical_protein_rows

client_genes_for_cognition <- all_hsp60_10_client_tbl |>
  dplyr::transmute(
    gene = canonical_gene_symbol(.data$gene)
  ) |>
  dplyr::filter(
    !is.na(.data$gene),
    nzchar(.data$gene)
  ) |>
  dplyr::distinct(.data$gene) |>
  dplyr::pull(.data$gene)

detected_client_genes_for_cognition <- intersect(
  client_genes_for_cognition,
  rownames(protein_mat_raw_cog)
)

if (
  length(detected_client_genes_for_cognition) !=
    n_detected_clients_used
) {
  stop(
    "Raw-protein cognition score detected ",
    length(detected_client_genes_for_cognition),
    " Hsp60/10 clients, but canonical production detection expects ",
    n_detected_clients_used,
    ".",
    call. = FALSE
  )
}

client_mat_raw_cog <- protein_mat_raw_cog[
  detected_client_genes_for_cognition,
  ,
  drop = FALSE
]

client_mat_raw_cog_z <- t(scale(t(client_mat_raw_cog)))
client_mat_raw_cog_z <- client_mat_raw_cog_z[
  rowSums(is.finite(client_mat_raw_cog_z)) > 0,
  ,
  drop = FALSE
]

if (nrow(client_mat_raw_cog_z) != n_detected_clients_used) {
  stop(
    "One or more detected Hsp60/10 clients had no usable variance ",
    "for the raw-protein cognition score.",
    call. = FALSE
  )
}

raw_hsp60_network_score <- colMeans(
  client_mat_raw_cog_z,
  na.rm = TRUE
)

raw_score_tbl <- tibble::tibble(
  SampleID = names(raw_hsp60_network_score),
  Hsp60_client_score = as.numeric(raw_hsp60_network_score)
)

if (
  nrow(raw_score_tbl) != 400 ||
  anyDuplicated(raw_score_tbl$SampleID) > 0
) {
  stop(
    "Raw Hsp60/10 cognition score is not one unique row per 400 protein samples.",
    call. = FALSE
  )
}

score_tbl <- prot_scores |>
  dplyr::transmute(
    SampleID = as.character(.data$SampleID),
    individualID = as.character(.data$IndividualID),
    cogdx_from_protein_metadata = safe_num(.data$cogdx_num),
    clinical_stage = .data$clinical_stage,
    clinical_stage_broad = .data$clinical_stage_broad,
    n_detected_clients_used = n_detected_clients_used,
    n_total_clients_reference = n_total_clients_reference
  ) |>
  checked_left_join(
    raw_score_tbl,
    by = "SampleID",
    label = "raw Hsp60/10 cognition score to protein sample metadata"
  )

if (any(!is.finite(score_tbl$Hsp60_client_score))) {
  stop(
    "At least one matrix-backed protein sample lacks a finite raw ",
    "Hsp60/10 cognition network score.",
    call. = FALSE
  )
}

if (nrow(score_tbl) != nrow(prot_meta_path)) {
  stop(
    "Cognition score table is not one row per matrix-backed protein sample: ",
    nrow(score_tbl), " versus ", nrow(prot_meta_path),
    call. = FALSE
  )
}

if (anyDuplicated(score_tbl$SampleID) > 0) {
  stop("Duplicate protein SampleID in cognition score table.", call. = FALSE)
}

if (anyDuplicated(score_tbl$individualID) > 0) {
  stop(
    "Duplicate protein participant ID in cognition score table.",
    call. = FALSE
  )
}

score_mean <- mean(score_tbl$Hsp60_client_score, na.rm = TRUE)
score_sd <- stats::sd(score_tbl$Hsp60_client_score, na.rm = TRUE)

if (!is.finite(score_sd) || score_sd <= 0) {
  stop("Hsp60/10 client score has no usable variance.", call. = FALSE)
}

score_tbl <- score_tbl |>
  dplyr::mutate(
    Hsp60_client_score_z =
      (.data$Hsp60_client_score - score_mean) / score_sd
  )

############################################################
## 2. ROSMAP clinical cognition variables
############################################################

clinical_raw <- readr::read_csv(
  cfg$rosmap_clinical_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (!"individualID" %in% colnames(clinical_raw)) {
  stop(
    "ROSMAP_clinical.csv must contain individualID for proteomics.",
    call. = FALSE
  )
}

clinical_cogdx_col <- require_alias_col(
  clinical_raw,
  c("cogdx", "COGDX"),
  "cogdx",
  context = "ROSMAP clinical metadata"
)

clinical_dcfdx_col <- require_alias_col(
  clinical_raw,
  c("dcfdx_lv"),
  "dcfdx_lv",
  context = "ROSMAP clinical metadata"
)

clinical_mmse_col <- require_alias_col(
  clinical_raw,
  c("cts_mmse30_lv", "mmse_last_valid"),
  "last-valid MMSE",
  context = "ROSMAP clinical metadata"
)

clinical_mmse_ad_col <- require_alias_col(
  clinical_raw,
  c("cts_mmse30_first_ad_dx", "mmse_first_ad_dx"),
  "first-AD-dx MMSE",
  context = "ROSMAP clinical metadata"
)

clinical_age_col <- require_alias_col(
  clinical_raw,
  alias_sets$age_death,
  "age at death",
  context = "ROSMAP clinical metadata"
)

clinical_educ_col <- require_alias_col(
  clinical_raw,
  c("educ", "education"),
  "education",
  context = "ROSMAP clinical metadata"
)

clinical_sex_col <- require_alias_col(
  clinical_raw,
  alias_sets$sex,
  "sex",
  context = "ROSMAP clinical metadata"
)

clinical_pmi_col <- require_alias_col(
  clinical_raw,
  alias_sets$pmi,
  "PMI",
  context = "ROSMAP clinical metadata"
)

clinical_braak_col <- require_alias_col(
  clinical_raw,
  alias_sets$braak,
  "Braak",
  context = "ROSMAP clinical metadata"
)

clinical_cerad_col <- require_alias_col(
  clinical_raw,
  alias_sets$cerad,
  "CERAD",
  context = "ROSMAP clinical metadata"
)

clinical_proj_col <- pick_existing(
  clinical_raw,
  c("projid", "proj_id")
)

clinical_apoe_col <- pick_existing(
  clinical_raw,
  c("apoe_genotype", "apoe")
)

clinical_sex_raw <- clinical_raw[[clinical_sex_col]]

clinical_tbl <- tibble::tibble(
  individualID = as.character(clinical_raw$individualID),
  projid = if (!is.na(clinical_proj_col)) {
    as.character(clinical_raw[[clinical_proj_col]])
  } else {
    NA_character_
  },
  cogdx = safe_num(clinical_raw[[clinical_cogdx_col]]),
  dcfdx_lv = safe_num(clinical_raw[[clinical_dcfdx_col]]),
  mmse_last_valid = safe_num(clinical_raw[[clinical_mmse_col]]),
  mmse_first_ad_dx = safe_num(clinical_raw[[clinical_mmse_ad_col]]),
  age_death = safe_num(clinical_raw[[clinical_age_col]]),
  educ = safe_num(clinical_raw[[clinical_educ_col]]),
  sex = factor(
    dplyr::case_when(
      clinical_sex_raw %in% c(1, "1", "Male", "male", "M", "m") ~ "Male",
      clinical_sex_raw %in% c(0, "0", "Female", "female", "F", "f") ~ "Female",
      TRUE ~ as.character(clinical_sex_raw)
    )
  ),
  pmi = safe_num(clinical_raw[[clinical_pmi_col]]),
  Braak = safe_num(clinical_raw[[clinical_braak_col]]),
  CERAD = safe_num(clinical_raw[[clinical_cerad_col]]),
  apoe_genotype = if (!is.na(clinical_apoe_col)) {
    safe_num(clinical_raw[[clinical_apoe_col]])
  } else {
    NA_real_
  }
) |>
  dplyr::filter(
    !is.na(.data$individualID),
    nzchar(.data$individualID)
  )

if (anyDuplicated(clinical_tbl$individualID) > 0) {
  stop(
    "ROSMAP clinical individualID is not unique in cognition builder.",
    call. = FALSE
  )
}

############################################################
## 3. Final cognition model table
############################################################

require_columns(
  prot_meta_adj,
  c("SampleID", "batch_factor"),
  "protein adjusted metadata for cognition batch"
)

batch_tbl <- prot_meta_adj |>
  dplyr::transmute(
    SampleID = as.character(.data$SampleID),
    batch_factor = as.factor(.data$batch_factor)
  ) |>
  dplyr::distinct(.data$SampleID, .keep_all = TRUE)

if (
  nrow(batch_tbl) != 400 ||
  sum(!is.na(batch_tbl$batch_factor)) != 400
) {
  stop(
    "Expected non-missing TMT batch_factor for all 400 protein samples ",
    "in cognition builder.",
    call. = FALSE
  )
}

cognition_model_df <- score_tbl |>
  checked_left_join(
    batch_tbl,
    by = "SampleID",
    label = "corrected raw Hsp60 score to TMT batch"
  ) |>
  checked_left_join(
    clinical_tbl,
    by = "individualID",
    label = "corrected raw Hsp60 score to ROSMAP cognition metadata"
  )

if (nrow(cognition_model_df) != nrow(score_tbl)) {
  stop("Cognition join changed protein row count.", call. = FALSE)
}

if (any(is.na(cognition_model_df$cogdx))) {
  stop(
    "At least one matrix-backed protein participant is missing cogdx ",
    "after cognition clinical join.",
    call. = FALSE
  )
}

cogdx_agreement <- cognition_model_df |>
  dplyr::filter(
    is.finite(.data$cogdx),
    is.finite(.data$cogdx_from_protein_metadata)
  ) |>
  dplyr::summarise(
    n_compared = dplyr::n(),
    n_agree = sum(.data$cogdx == .data$cogdx_from_protein_metadata),
    n_disagree = sum(.data$cogdx != .data$cogdx_from_protein_metadata)
  )

if (cogdx_agreement$n_disagree != 0) {
  stop(
    "Cognition clinical cogdx disagrees with protein-ingestion cogdx.",
    call. = FALSE
  )
}

cognition_availability <- tibble::tibble(
  variable = c(
    "Hsp60_client_score",
    "cogdx",
    "dcfdx_lv",
    "mmse_last_valid",
    "mmse_first_ad_dx",
    "age_death",
    "educ",
    "sex",
    "pmi",
    "Braak",
    "CERAD",
    "batch_factor"
  ),
  n_available = c(
    sum(is.finite(cognition_model_df$Hsp60_client_score)),
    sum(is.finite(cognition_model_df$cogdx)),
    sum(is.finite(cognition_model_df$dcfdx_lv)),
    sum(is.finite(cognition_model_df$mmse_last_valid)),
    sum(is.finite(cognition_model_df$mmse_first_ad_dx)),
    sum(is.finite(cognition_model_df$age_death)),
    sum(is.finite(cognition_model_df$educ)),
    sum(!is.na(cognition_model_df$sex)),
    sum(is.finite(cognition_model_df$pmi)),
    sum(is.finite(cognition_model_df$Braak)),
    sum(is.finite(cognition_model_df$CERAD)),
    sum(!is.na(cognition_model_df$batch_factor))
  ),
  n_total = nrow(cognition_model_df)
)

write_tbl(
  cognition_model_df,
  "cognition_model_df_COVARIATE_ADJUSTED_CLIENT_SCORE"
)

write_tbl(
  cognition_availability,
  "cognition_model_df_variable_availability"
)

write_tbl(
  cogdx_agreement,
  "cognition_model_df_cogdx_join_audit"
)

cognition_model_specification <- tibble::tibble(
  protein_abundance_source = "prot_mat_raw",
  protein_pre_residualized = FALSE,
  direct_model_covariates =
    "age_death;sex;educ;pmi;Braak;CERAD;batch_factor",
  note = paste(
    "Cognition inference uses processed non-residualized protein",
    "abundance with nuisance/pathology covariates and TMT batch",
    "included directly in each regression."
  )
)

write_tbl(
  cognition_model_specification,
  "cognition_model_specification"
)

save_obj(cognition_model_df, "cognition_model_df")

message("Loaded 06_build_cognition_objects.R")
