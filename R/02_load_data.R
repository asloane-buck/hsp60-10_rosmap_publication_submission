############################################################
## 02_load_data.R
## Load RNA, proteomics, and metadata.
############################################################

require_objects("cfg", context = "02_load_data.R")

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
    mutate(
      sample_id = as.character(sample_id),
      individual_id = as.character(individual_id),
      diagnosis_stage = case_when(
        diagnosis %in% c(1, "1", "Control") ~ "Control",
        diagnosis %in% c(2, "2", "Early_AD", "MCI", "Intermediate") ~ "Early_AD",
        diagnosis %in% c(4, "4", "AD") ~ "AD",
        TRUE ~ NA_character_
      ),
      diagnosis_stage = factor(diagnosis_stage, levels = c("Control", "Early_AD", "AD")),
      braak_num = safe_num(braak),
      cerad_num = safe_num(cerad),
      sex_label = case_when(
        sex %in% c(0, "0", "Female", "female", "F", "f") ~ "Female",
        sex %in% c(1, "1", "Male", "male", "M", "m") ~ "Male",
        TRUE ~ as.character(sex)
      )
    )
}

load_rna_matrix <- function() {
  rna_file <- first_existing(
    c(
      file.path(cfg$derived_dir, "ROSMAP_vst_gene_symbol_matrix_STAGE.csv"),
      file.path(cfg$project_dir, "ROSMAP_pathway_scores_robustness", "ROSMAP_vst_gene_symbol_matrix_STAGE.csv")
    ),
    "ROSMAP VST gene-symbol RNA matrix"
  )
  read_gene_matrix_csv(rna_file, gene_col = "gene_symbol")
}

load_protein <- function() {
  prot_file <- first_existing(
    c(
      file.path(cfg$derived_dir, "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_copy.csv"),
      file.path(cfg$proteomics_input_dir, "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"),
      file.path(cfg$proteomics_input_dir, "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_copy.csv")
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

  meta <- readr::read_csv(meta_file, show_col_types = FALSE, name_repair = "minimal") |>
    janitor::clean_names()

  names(meta) <- stringr::str_replace_all(names(meta), "\\.", "_")

  meta <- add_alias_column(
    meta,
    "batch_channel",
    c("batch_channel", "batch.channel", "batch", "channel", "tmt_channel"),
    context = "protein matched metadata"
  )
  meta <- add_alias_column(meta, "sample_id", alias_sets$sample_id, context = "protein matched metadata")
  meta <- add_alias_column(meta, "individual_id", alias_sets$individual_id, context = "protein matched metadata")
  meta <- add_alias_column(
    meta,
    "emory_strict_dx_2019",
    c("emory_strict_dx_2019", "emorystrictdx_2019", "EmoryStrictDx.2019", "diagnosis", "cogdx"),
    context = "protein matched metadata"
  )

  required <- c("batch_channel", "sample_id", "individual_id", "emory_strict_dx_2019")
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop(
      "Protein metadata missing: ", paste(missing, collapse = ", "),
      ". Available columns: ", available_cols_msg(meta),
      call. = FALSE
    )
  }

  meta <- meta |>
    mutate(
      batch_channel = as.character(batch_channel),
      SampleID = as.character(sample_id),
      IndividualID = as.character(individual_id),
      EmoryStrictDx.2019 = as.character(emory_strict_dx_2019)
    ) |>
    filter(EmoryStrictDx.2019 %in% c("Control", "AsymAD", "AD")) |>
    distinct(batch_channel, .keep_all = TRUE)

  common_channels <- intersect(colnames(mat), meta$batch_channel)
  if (length(common_channels) < 10) stop("Too few protein matrix columns matched metadata batch_channel.", call. = FALSE)

  mat <- mat[, common_channels, drop = FALSE]
  meta <- meta[match(common_channels, meta$batch_channel), , drop = FALSE]
  colnames(mat) <- meta$SampleID

  meta <- meta |>
    mutate(
      SampleID = colnames(mat),
      EmoryStrictDx.2019 = factor(EmoryStrictDx.2019, levels = c("Control", "AsymAD", "AD"))
    )

  list(mat = mat, meta = meta)
}

analysis_meta <- load_analysis_meta()

rna_mat_raw <- load_rna_matrix()
protein <- load_protein()
prot_mat_raw <- protein$mat
prot_meta <- protein$meta

rna_meta <- analysis_meta |>
  filter(sample_id %in% colnames(rna_mat_raw), !is.na(diagnosis_stage)) |>
  distinct(sample_id, .keep_all = TRUE)

rna_mat_raw <- rna_mat_raw[, rna_meta$sample_id, drop = FALSE]
rna_meta <- rna_meta[match(colnames(rna_mat_raw), rna_meta$sample_id), , drop = FALSE]
if (!all(colnames(rna_mat_raw) == rna_meta$sample_id)) {
  stop(
    "RNA metadata could not be aligned to RNA matrix sample_id. Available columns: ",
    available_cols_msg(rna_meta),
    call. = FALSE
  )
}

prot_meta_path <- prot_meta |>
  checked_left_join(
    analysis_meta |>
      select(individual_id, braak_num, cerad_num, diagnosis, sex_label, any_of(c("region", "tissue", "age", "age_death", "pmi", "rin"))) |>
      distinct(individual_id, .keep_all = TRUE),
    by = c("IndividualID" = "individual_id"),
    label = "protein metadata to analysis metadata"
  )

## Backward-compatible names after adjustment scripts.
rna_mat <- rna_mat_raw
prot_mat <- prot_mat_raw

write_tbl(tibble(object = c("rna_mat_raw", "prot_mat_raw", "rna_meta", "prot_meta_path"),
                 nrow = c(nrow(rna_mat_raw), nrow(prot_mat_raw), nrow(rna_meta), nrow(prot_meta_path)),
                 ncol = c(ncol(rna_mat_raw), ncol(prot_mat_raw), ncol(rna_meta), ncol(prot_meta_path))),
          "load_data_dimensions")

save_obj(analysis_meta, "analysis_meta")
save_obj(rna_mat_raw, "rna_mat_raw")
save_obj(prot_mat_raw, "prot_mat_raw")
save_obj(rna_meta, "rna_meta")
save_obj(prot_meta_path, "prot_meta_path")

message("Loaded 02_load_data.R")
