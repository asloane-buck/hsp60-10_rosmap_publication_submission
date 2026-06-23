
############################################################
## 06_build_cognition_objects.R
## Literal migration of the upstream Figure 4 cognition builder.
## Uses adjusted prot_mat from 03/05 for Hsp60/10 client score.
############################################################

require_objects(
  c("cfg", "prot_mat", "prot_scores", "all_hsp60_10_client_tbl", "prot_meta_path", "analysis_meta"),
  context = "06_build_cognition_objects.R"
)

if (!file.exists(cfg$rosmap_clinical_file)) {
  stop("Missing ROSMAP clinical file required for Figure 4: ", cfg$rosmap_clinical_file,
       "\nUpdate cfg$rosmap_clinical_file in 00_config.R.", call. = FALSE)
}

rosmap_clinical <- readr::read_csv(cfg$rosmap_clinical_file, show_col_types = FALSE) |>
  janitor::clean_names()

## Backward-compatible aliases used by the original cognition builder.
protein_meta_raw <- prot_meta_path |>
  standardize_common_aliases(context = "protein metadata for cognition", require_ids = FALSE) |>
  dplyr::mutate(
    sample_id = coalesce_alias_chr(cur_data(), c("SampleID", "sample_id")),
    individual_id = coalesce_alias_chr(cur_data(), c("IndividualID", "individual_id", "individualID", "projid"))
  )

individual_meta_raw <- analysis_meta |>
  standardize_common_aliases(context = "individual metadata for cognition", require_ids = TRUE) |>
  dplyr::mutate(individual_id = as.character(individual_id)) |>
  dplyr::distinct(individual_id, .keep_all = TRUE)

############################################################
## Exact upstream builder from ROSMAP cognition.Rmd
############################################################

needed_objects <- c(
  "prot_mat",
  "prot_scores",
  "protein_meta_raw",
  "individual_meta_raw",
  "all_hsp60_10_client_tbl"
)

missing_objects <- needed_objects[!needed_objects %in% ls(envir = .GlobalEnv)]

if (length(missing_objects) > 0) {
  stop(
    "Cognition builder missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    ". Available objects: ", available_objects_msg(),
    call. = FALSE
  )
}

all_hsp_clients <- all_hsp60_10_client_tbl %>%
  transmute(
    gene = toupper(as.character(gene))
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)

cat("\nCurated Hsp60/10 clients:", nrow(all_hsp_clients), "\n")

protein_mat <- prot_mat
rownames(protein_mat) <- toupper(rownames(protein_mat))

detected_clients <- intersect(
  rownames(protein_mat),
  all_hsp_clients$gene
)

missing_clients <- setdiff(
  all_hsp_clients$gene,
  detected_clients
)

cat("Detected Hsp60/10 clients in prot_mat:", length(detected_clients), "\n")
cat("Missing Hsp60/10 clients:", length(missing_clients), "\n")

missing_clients_tbl <- tibble(
  missing_client = missing_clients
)

print(missing_clients_tbl, n = Inf)

client_mat <- protein_mat[detected_clients, , drop = FALSE]

client_mat <- apply(client_mat, 2, as.numeric)
rownames(client_mat) <- detected_clients

client_mat_z <- t(scale(t(client_mat)))

client_mat_z <- client_mat_z[
  rowSums(!is.na(client_mat_z)) > 0,
  ,
  drop = FALSE
]

manual_all_client_score <- colMeans(client_mat_z, na.rm = TRUE)

manual_score_tbl <- tibble(
  SampleID = names(manual_all_client_score),
  verified_all_detected_Hsp60_10_client_score = as.numeric(manual_all_client_score),
  n_detected_clients_used = nrow(client_mat_z),
  n_total_clients_reference = nrow(all_hsp_clients)
)

coalesce_existing_chr <- function(data, candidates) {
  present <- candidates[candidates %in% colnames(data)]
  if (length(present) == 0) return(rep(NA_character_, nrow(data)))
  out <- rep(NA_character_, nrow(data))
  for (cc in present) out <- coalesce(out, as.character(data[[cc]]))
  out
}

coalesce_existing_num <- function(data, candidates) {
  readr::parse_number(coalesce_existing_chr(data, candidates))
}

df <- manual_score_tbl %>%
  checked_left_join(
    prot_scores %>%
      mutate(SampleID = as.character(SampleID)),
    by = "SampleID",
    label = "manual Hsp60 score to pathway scores"
  ) %>%
  checked_left_join(
    protein_meta_raw,
    by = c("SampleID" = "sample_id"),
    label = "manual Hsp60 score to protein metadata",
    suffix = c("_score", "_meta")
  )

df <- df %>%
  mutate(
    join_individual_id = coalesce_existing_chr(
      cur_data(),
      c(
        "IndividualID_score",
        "IndividualID_meta",
        "IndividualID",
        "individual_id_score",
        "individual_id_meta",
        "individual_id"
      )
    )
  )

individual_join_meta <- tibble::tibble(
  join_individual_id = as.character(individual_meta_raw$individual_id),
  age_death_individual = coalesce_alias_chr(individual_meta_raw, alias_sets$age_death),
  pmi_individual = coalesce_alias_chr(individual_meta_raw, alias_sets$pmi),
  sex_individual = coalesce_alias_chr(individual_meta_raw, c("sex_label", alias_sets$sex))
)

df2 <- df %>%
  checked_left_join(
    individual_join_meta,
    by = "join_individual_id",
    label = "sample cognition metadata to individual metadata"
  )

model_df <- df2 %>%
  mutate(
    Hsp60_client_score = as.numeric(verified_all_detected_Hsp60_10_client_score),
    Hsp60_client_score_z = as.numeric(scale(Hsp60_client_score)),
    age_death = coalesce_existing_num(cur_data(), c("age_death", "age_death_individual", "age", "age_meta", "age_death_meta")),
    pmi = coalesce_existing_num(cur_data(), c("pmi_individual", "pmi", "pmi_meta", "pmi_score")),
    sex = as.factor(coalesce_existing_chr(cur_data(), c("sex_label", "sex_label_meta", "sex_individual", "sex_meta", "sex"))),
    Braak = coalesce_existing_num(cur_data(), c("braak_num_score", "braak_num_meta", "braak_num", "braak", "braak_meta")),
    CERAD = coalesce_existing_num(cur_data(), c("cerad_num_score", "cerad_num_meta", "cerad_num", "cerad", "cerad_meta")),
    diagnosis = as.factor(coalesce_existing_chr(cur_data(), c("EmoryStrictDx.2019_score", "emory_strict_dx_2019", "EmoryStrictDx.2019_meta", "diagnosis", "emory_coded_dx_2017")))
  ) %>%
  select(
    SampleID,
    join_individual_id,
    Hsp60_client_score,
    Hsp60_client_score_z,
    age_death,
    sex,
    pmi,
    Braak,
    CERAD,
    diagnosis,
    n_detected_clients_used,
    n_total_clients_reference
  )

model_df2 <- model_df %>%
  mutate(
    diagnosis_numeric = case_when(
      diagnosis == "Control" ~ 0,
      diagnosis == "AsymAD" ~ 1,
      diagnosis == "Intermediate" ~ 1,
      diagnosis == "AD" ~ 2,
      TRUE ~ NA_real_
    )
  )

score_tbl <- model_df2 %>%
  transmute(
    SampleID,
    individualID = as.character(join_individual_id),
    Hsp60_client_score = as.numeric(Hsp60_client_score),
    Hsp60_client_score_z = as.numeric(scale(Hsp60_client_score)),
    n_detected_clients_used,
    n_total_clients_reference
  ) %>%
  distinct()

clinical_id_col <- require_alias_col(
  rosmap_clinical,
  c("individual_id", "individualid", "individualID", "projid"),
  "clinical individual ID",
  context = "ROSMAP clinical metadata"
)
clinical_proj_col <- pick_existing(rosmap_clinical, c("projid", "individual_id", "individualid", "individualID"))
clinical_cogdx_col <- require_alias_col(rosmap_clinical, c("cogdx", "diagnosis"), "cogdx", "ROSMAP clinical metadata")
clinical_dcfdx_col <- require_alias_col(rosmap_clinical, c("dcfdx_lv"), "dcfdx_lv", "ROSMAP clinical metadata")
clinical_mmse_col <- require_alias_col(rosmap_clinical, c("cts_mmse30_lv", "mmse_last_valid"), "last-valid MMSE", "ROSMAP clinical metadata")
clinical_mmse_ad_col <- require_alias_col(rosmap_clinical, c("cts_mmse30_first_ad_dx", "mmse_first_ad_dx"), "first-AD-dx MMSE", "ROSMAP clinical metadata")
clinical_age_col <- require_alias_col(rosmap_clinical, alias_sets$age_death, "age at death", "ROSMAP clinical metadata")
clinical_educ_col <- require_alias_col(rosmap_clinical, c("educ", "education"), "education", "ROSMAP clinical metadata")
clinical_sex_col <- require_alias_col(rosmap_clinical, alias_sets$sex, "sex", "ROSMAP clinical metadata")
clinical_pmi_col <- require_alias_col(rosmap_clinical, alias_sets$pmi, "PMI", "ROSMAP clinical metadata")
clinical_braak_col <- require_alias_col(rosmap_clinical, alias_sets$braak, "Braak", "ROSMAP clinical metadata")
clinical_cerad_col <- require_alias_col(rosmap_clinical, alias_sets$cerad, "CERAD", "ROSMAP clinical metadata")
clinical_apoe_col <- pick_existing(rosmap_clinical, c("apoe_genotype", "apoe"))

clinical_sex_raw <- rosmap_clinical[[clinical_sex_col]]

clinical_tbl <- tibble::tibble(
  individualID = as.character(rosmap_clinical[[clinical_id_col]]),
  projid = if (!is.na(clinical_proj_col)) as.character(rosmap_clinical[[clinical_proj_col]]) else NA_character_,
  cogdx = safe_num(rosmap_clinical[[clinical_cogdx_col]]),
  dcfdx_lv = safe_num(rosmap_clinical[[clinical_dcfdx_col]]),
  mmse_last_valid = safe_num(rosmap_clinical[[clinical_mmse_col]]),
  mmse_first_ad_dx = safe_num(rosmap_clinical[[clinical_mmse_ad_col]]),
  age_death = safe_num(rosmap_clinical[[clinical_age_col]]),
  educ = safe_num(rosmap_clinical[[clinical_educ_col]]),
  sex = factor(case_when(clinical_sex_raw == 1 ~ "Male", clinical_sex_raw == 0 ~ "Female", TRUE ~ as.character(clinical_sex_raw))),
  pmi = safe_num(rosmap_clinical[[clinical_pmi_col]]),
  Braak = safe_num(rosmap_clinical[[clinical_braak_col]]),
  CERAD = safe_num(rosmap_clinical[[clinical_cerad_col]]),
  apoe_genotype = if (!is.na(clinical_apoe_col)) safe_num(rosmap_clinical[[clinical_apoe_col]]) else NA_real_
)

cognition_model_df <- score_tbl %>%
  checked_left_join(clinical_tbl, by = "individualID", label = "client score to ROSMAP clinical")

write_tbl(cognition_model_df, "cognition_model_df_COVARIATE_ADJUSTED_CLIENT_SCORE")
save_obj(cognition_model_df, "cognition_model_df")
message("Loaded 06_build_cognition_objects.R")
