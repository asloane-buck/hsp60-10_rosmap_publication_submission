# ============================================================
# Load supplemental figure inputs
# ============================================================

if (!exists("candidate_csv_paths")) {
  stop("Run 00_supplemental_config.R before 01_supplemental_load_inputs.R.", call. = FALSE)
}

inputs <- list()
input_source_audit <- tibble::tibble(
  input_name = character(),
  source_type = character(),
  path = character(),
  status = character(),
  details = character()
)

first_or_na <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA_character_
  } else {
    x[1]
  }
}

record_input_source <- function(input_name, source_type, path, status, details = "") {
  input_source_audit <<- dplyr::bind_rows(
    input_source_audit,
    tibble::tibble(
      input_name = input_name,
      source_type = source_type,
      path = path,
      status = status,
      details = details
    )
  )
}

store_input <- function(input_name, object, source_type, path) {
  if (input_name %in% names(inputs) && !isTRUE(overwrite_inputs)) {
    message(
      "Input already loaded and overwrite_inputs is FALSE, keeping existing: ",
      input_name
    )
    record_input_source(
      input_name,
      source_type,
      path,
      "skipped_duplicate",
      "Existing object retained"
    )
    return(invisible(FALSE))
  }

  inputs[[input_name]] <<- object
  object_dim <- paste(dim(object), collapse = " x ")
  if (identical(object_dim, "")) {
    object_dim <- paste("length", length(object))
  }

  message("Loaded input: ", input_name, " [", object_dim, "]")
  record_input_source(input_name, source_type, path, "loaded", object_dim)
  invisible(TRUE)
}

first_existing_path <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  first_or_na(paths[file.exists(paths)])
}

read_table_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
}

clean_gene_for_loading <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_to_upper(x)
  x <- stringr::str_replace(x, "\\.\\d+$", "")
  x <- dplyr::case_when(
    x == "ATP5B" ~ "ATP5F1B",
    x == "ATP5A1" ~ "ATP5F1A",
    x == "ATP5C1" ~ "ATP5F1C",
    TRUE ~ x
  )
  x[x == ""] <- NA_character_
  x
}

pick_present_col <- function(df, candidates) {
  selected <- candidates[candidates %in% colnames(df)]
  if (length(selected) == 0) NA_character_ else selected[1]
}

# ============================================================
# Load CSV/TSV candidates
# ============================================================

for (input_name in names(candidate_csv_paths)) {
  candidate_paths <- candidate_csv_paths[[input_name]]
  if (length(candidate_paths) == 0) {
    message("No CSV candidates configured for: ", input_name)
    record_input_source(input_name, "csv", "", "not_configured", "")
    next
  }

  selected_path <- first_existing_path(candidate_paths)
  if (is.na(selected_path)) {
    message("No existing CSV/TSV found for: ", input_name)
    message("Candidates tried: ", paste(candidate_paths, collapse = " | "))
    record_input_source(
      input_name,
      "csv",
      paste(candidate_paths, collapse = " | "),
      "missing",
      "No candidate file exists"
    )
    next
  }

  message("Reading table for ", input_name, ": ", selected_path)
  object <- read_table_file(selected_path)
  message("  Dimensions: ", paste(dim(object), collapse = " x "))
  message("  Columns: ", paste(colnames(object), collapse = ", "))
  store_input(input_name, object, "csv", selected_path)
}

# ============================================================
# Load RDS candidates
# ============================================================

for (input_name in names(candidate_rds_paths)) {
  path <- candidate_rds_paths[[input_name]]
  if (is.na(path) || !nzchar(path)) {
    message("No RDS path configured for: ", input_name)
    record_input_source(input_name, "rds", "", "not_configured", "")
    next
  }

  if (!file.exists(path)) {
    message("RDS file not found for ", input_name, ": ", path)
    record_input_source(input_name, "rds", path, "missing", "")
    next
  }

  message("Reading RDS for ", input_name, ": ", path)
  object <- readRDS(path)
  store_input(input_name, object, "rds", path)

  if (is.list(object) && !is.data.frame(object) && !is.null(names(object))) {
    retained_names <- intersect(names(object), candidate_object_names)
    for (element_name in retained_names) {
      store_input(element_name, object[[element_name]], "rds_list_element", path)
    }
  }
}

# ============================================================
# Load selected objects from RData candidates
# ============================================================

for (workspace_name in names(candidate_rdata_paths)) {
  path <- candidate_rdata_paths[[workspace_name]]
  if (is.na(path) || !nzchar(path)) {
    message("No RData path configured for: ", workspace_name)
    record_input_source(workspace_name, "rdata", "", "not_configured", "")
    next
  }

  if (!file.exists(path)) {
    message("RData file not found for ", workspace_name, ": ", path)
    record_input_source(workspace_name, "rdata", path, "missing", "")
    next
  }

  message("Loading RData workspace: ", path)
  workspace_env <- new.env(parent = baseenv())
  loaded_names <- load(path, envir = workspace_env)
  relevant_names <- intersect(loaded_names, candidate_object_names)

  if (length(relevant_names) == 0) {
    message("  No relevant objects retained from workspace.")
    record_input_source(workspace_name, "rdata", path, "loaded_none_retained", "")
    next
  }

  message("  Retaining relevant objects: ", paste(relevant_names, collapse = ", "))
  for (object_name in relevant_names) {
    store_input(object_name, get(object_name, envir = workspace_env), "rdata", path)
  }
}

# ============================================================
# Canonicalize HSP60/HSP10 interactor inventory input
# ============================================================

interactor_inventory_object_names <- c(
  "hsp60_hsp10_interactor_inventory",
  "all_hsp60_10_client_tbl",
  "hsp60_client_tbl",
  "hsp_client_tbl",
  "hsp_clients_all",
  "hsp60_clients"
)

canonicalize_interactor_inventory <- function() {
  selected_name <- interactor_inventory_object_names[interactor_inventory_object_names %in% names(inputs)][1]
  if (is.na(selected_name)) {
    message("No HSP60/HSP10 interactor inventory input detected yet.")
    return(invisible(FALSE))
  }

  if (!"hsp60_hsp10_interactor_inventory" %in% names(inputs)) {
    inputs$hsp60_hsp10_interactor_inventory <<- inputs[[selected_name]]
    record_input_source(
      "hsp60_hsp10_interactor_inventory",
      "alias",
      selected_name,
      "loaded",
      "Canonical alias for Bie et al. 2020 HSP60/HSP10 interactor inventory"
    )
    message("Canonical HSP60/HSP10 interactor inventory input: hsp60_hsp10_interactor_inventory")
    message("  Source object/file input name: ", selected_name)
  } else {
    message("Canonical HSP60/HSP10 interactor inventory input detected: hsp60_hsp10_interactor_inventory")
  }

  invisible(TRUE)
}

canonicalize_interactor_inventory()

# ============================================================
# Align protein metadata to protein matrix samples
# ============================================================

coerce_matrix_for_loading <- function(object) {
  if (is.null(object)) {
    return(NULL)
  }
  if (is.matrix(object)) {
    return(object)
  }
  if (is.data.frame(object)) {
    df <- as.data.frame(object)
    numeric_cols <- vapply(df, is.numeric, logical(1))
    if (!any(numeric_cols)) {
      return(NULL)
    }
    mat <- as.matrix(df[, numeric_cols, drop = FALSE])
    rownames(mat) <- rownames(df)
    return(mat)
  }
  NULL
}

get_interactor_genes_for_loading <- function() {
  selected_name <- interactor_inventory_object_names[interactor_inventory_object_names %in% names(inputs)][1]
  if (is.na(selected_name)) {
    return(character())
  }

  object <- inputs[[selected_name]]
  if (is.atomic(object) && !is.matrix(object)) {
    return(unique(stats::na.omit(clean_gene_for_loading(object))))
  }

  df <- as.data.frame(object)
  gene_col <- pick_present_col(
    df,
    c("gene", "gene_symbol", "symbol", "hgnc_symbol", "curated_gene", "display_gene")
  )
  if (is.na(gene_col)) {
    return(character())
  }
  unique(stats::na.omit(clean_gene_for_loading(df[[gene_col]])))
}

infer_protein_sample_ids_for_loading <- function() {
  matrix_names <- c(
    "prot_mat",
    "protein_matrix",
    "protein_mat",
    "tmt_mat",
    "prot_expr_mat",
    "prot_mat_in",
    "prot_mat_raw"
  )
  selected_matrix_name <- matrix_names[matrix_names %in% names(inputs)][1]
  if (is.na(selected_matrix_name)) {
    message("No protein matrix available for protein metadata alignment.")
    return(NULL)
  }

  protein_mat <- coerce_matrix_for_loading(inputs[[selected_matrix_name]])
  if (is.null(protein_mat)) {
    message("Protein matrix candidate could not be coerced for metadata alignment: ", selected_matrix_name)
    return(NULL)
  }

  interactor_genes <- get_interactor_genes_for_loading()
  row_gene_hits <- if (length(interactor_genes) > 0) {
    sum(clean_gene_for_loading(rownames(protein_mat)) %in% interactor_genes)
  } else {
    0
  }
  column_gene_hits <- if (length(interactor_genes) > 0) {
    sum(clean_gene_for_loading(colnames(protein_mat)) %in% interactor_genes)
  } else {
    0
  }

  if (row_gene_hits >= column_gene_hits) {
    orientation <- "genes_rows_samples_columns"
    sample_ids <- colnames(protein_mat)
  } else {
    orientation <- "genes_columns_samples_rows"
    sample_ids <- rownames(protein_mat)
  }

  message("Protein matrix selected for metadata alignment: ", selected_matrix_name)
  message("Protein matrix dimensions: ", paste(dim(protein_mat), collapse = " x "))
  message("Protein matrix orientation: ", orientation)
  message("Protein matrix row gene hits: ", row_gene_hits)
  message("Protein matrix column gene hits: ", column_gene_hits)
  message("First 10 protein matrix sample IDs: ", paste(head(sample_ids, 10), collapse = ", "))

  list(
    matrix_name = selected_matrix_name,
    matrix = protein_mat,
    sample_ids = as.character(sample_ids),
    orientation = orientation
  )
}

protein_metadata_id_candidates <- c(
  "SampleID",
  "sample_id",
  "sample",
  "Sample",
  "IndividualID",
  "individualID",
  "projid",
  "SpecimenID",
  "specimenID",
  "specimen_id",
  "batch.channel",
  "batch_channel",
  "MulticonsensusStudyFileID",
  "protein_id",
  "TMT_channel",
  "channel",
  "id"
)

protein_stage_candidates <- c(
  "clinical_stage",
  "clinical_stage_broad",
  "cogdx_num",
  "cogdx",
  "diagnosis_stage",
  "diagnosis"
)

protein_covariate_candidates <- list(
  age = c("age", "age_at_death", "age_death", "age_c", "Age", "ageAtDeath"),
  sex = c("sex", "msex", "sex_label", "Sex", "gender"),
  PMI = c("pmi", "PMI", "pmi_c", "postmortem_interval", "postmortem.interval"),
  RIN = c("rin", "RIN", "rin_c", "RINvalue", "rna_integrity_number"),
  Braak = c("braak", "Braak", "braak_num", "braaksc", "braak_stage", "braak_stage_num", "Braak.Stage"),
  CERAD = c("cerad", "CERAD", "cerad_num", "ceradsc", "cerad_score", "CERAD.Score")
)

score_metadata_covariates <- function(df) {
  detected_covariates <- vapply(
    protein_covariate_candidates,
    function(candidates) !is.na(pick_present_col(df, candidates)),
    logical(1)
  )
  sum(detected_covariates)
}

test_protein_metadata_overlap <- function(df, object_name, id_col, protein_sample_ids) {
  ids <- as.character(df[[id_col]])
  ids <- ids[!is.na(ids) & nzchar(ids)]
  overlap_ids <- intersect(protein_sample_ids, ids)
  non_overlap <- setdiff(protein_sample_ids, ids)
  duplicate_ids <- ids[duplicated(ids)]
  duplicate_overlap_count <- sum(unique(duplicate_ids) %in% protein_sample_ids)

  tibble::tibble(
    metadata_object = object_name,
    id_column = id_col,
    n_metadata_rows = nrow(df),
    n_unique_ids = dplyr::n_distinct(ids),
    overlap_count = length(overlap_ids),
    percent_overlap = length(overlap_ids) / length(protein_sample_ids),
    duplicate_overlap_ids = duplicate_overlap_count,
    exact_one_row_per_protein_sample = nrow(df) == length(protein_sample_ids) &&
      dplyr::n_distinct(ids) == nrow(df),
    stage_available = !is.na(pick_present_col(df, protein_stage_candidates)),
    n_detected_covariate_columns = score_metadata_covariates(df),
    id_priority = match(id_col, protein_metadata_id_candidates),
    example_overlapping_ids = paste(head(overlap_ids, 5), collapse = "; "),
    example_non_overlapping_protein_sample_ids = paste(head(non_overlap, 5), collapse = "; ")
  )
}

build_protein_metadata_candidates <- function() {
  candidate_names <- names(inputs)[vapply(inputs, is.data.frame, logical(1))]
  candidate_names <- candidate_names[
    grepl(
      "prot|protein|tmt|sample_meta|clinical|biospecimen",
      candidate_names,
      ignore.case = TRUE
    ) |
      vapply(inputs[candidate_names], function(df) {
        any(protein_metadata_id_candidates %in% colnames(df))
      }, logical(1))
  ]

  candidate_names
}

align_protein_metadata <- function() {
  protein_info <- infer_protein_sample_ids_for_loading()
  if (is.null(protein_info)) {
    return(invisible(FALSE))
  }

  protein_sample_ids <- protein_info$sample_ids
  candidate_names <- build_protein_metadata_candidates()

  if (length(candidate_names) == 0) {
    message("No protein metadata candidates found for alignment.")
    return(invisible(FALSE))
  }

  message("Protein metadata candidate objects: ", paste(candidate_names, collapse = ", "))
  for (candidate_name in candidate_names) {
    df <- as.data.frame(inputs[[candidate_name]])
    message("Candidate protein metadata object: ", candidate_name)
    message("  Dimensions: ", paste(dim(df), collapse = " x "))
    message("  Columns: ", paste(colnames(df), collapse = ", "))
  }

  overlap_tbl <- purrr::map_dfr(candidate_names, function(candidate_name) {
    df <- as.data.frame(inputs[[candidate_name]])
    ids_to_test <- intersect(protein_metadata_id_candidates, colnames(df))
    if (!is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
      df <- tibble::rownames_to_column(df, ".rownames")
      ids_to_test <- c(ids_to_test, ".rownames")
    }

    purrr::map_dfr(
      ids_to_test,
      ~ test_protein_metadata_overlap(df, candidate_name, .x, protein_sample_ids)
    )
  })

  if (nrow(overlap_tbl) == 0) {
    message("No candidate protein metadata ID columns were testable.")
    return(invisible(FALSE))
  }

  message("Protein metadata overlap diagnostics:")
  print(
    overlap_tbl |>
      dplyr::arrange(dplyr::desc(.data$overlap_count), .data$metadata_object, .data$id_column),
    n = Inf
  )

  best <- overlap_tbl |>
    dplyr::filter(.data$duplicate_overlap_ids == 0) |>
    dplyr::arrange(
      dplyr::desc(.data$overlap_count),
      dplyr::desc(.data$percent_overlap),
      dplyr::desc(.data$exact_one_row_per_protein_sample),
      dplyr::desc(.data$n_detected_covariate_columns),
      dplyr::desc(.data$stage_available),
      .data$id_priority,
      .data$metadata_object
    ) |>
    dplyr::slice(1)

  if (nrow(best) == 0 || best$percent_overlap < 0.90) {
    message("Best protein metadata overlap was poor; not guessing.")
    print(overlap_tbl |> dplyr::arrange(dplyr::desc(.data$overlap_count)) |> head(20))
    stop("No strong protein metadata overlap found.", call. = FALSE)
  }

  selected_name <- best$metadata_object
  selected_id_col <- best$id_column
  selected_df <- as.data.frame(inputs[[selected_name]])
  if (identical(selected_id_col, ".rownames")) {
    selected_df <- tibble::rownames_to_column(selected_df, ".rownames")
  }
  selected_df[[selected_id_col]] <- as.character(selected_df[[selected_id_col]])

  aligned_meta <- selected_df[match(protein_sample_ids, selected_df[[selected_id_col]]), , drop = FALSE]
  rownames(aligned_meta) <- protein_sample_ids

  message("Selected protein metadata object name: ", selected_name)
  message("Selected protein sample ID column: ", selected_id_col)
  message("Protein metadata overlap count: ", best$overlap_count)
  message("Protein metadata dimensions before alignment: ", paste(dim(selected_df), collapse = " x "))
  message("Protein metadata dimensions after alignment: ", paste(dim(aligned_meta), collapse = " x "))

  if (!identical(as.character(aligned_meta[[selected_id_col]]), protein_sample_ids)) {
    stop("Protein metadata alignment failed: selected ID column is not in protein matrix sample order.", call. = FALSE)
  }

  message(
    "Confirmed protein metadata alignment: nrow(inputs$protein_meta) equals protein sample count and ",
    selected_id_col,
    " matches protein matrix sample IDs in order."
  )

  clinical_join_details <- "not_attempted"
  if ("protein_clinical_meta" %in% names(inputs)) {
    clinical_df <- as.data.frame(inputs$protein_clinical_meta)
    if ("IndividualID" %in% colnames(aligned_meta) && "individualID" %in% colnames(clinical_df)) {
      join_overlap <- sum(unique(as.character(aligned_meta$IndividualID)) %in% unique(as.character(clinical_df$individualID)))
      clinical_join_details <- paste0("IndividualID overlap ", join_overlap, "/", nrow(aligned_meta))
      message("Protein clinical metadata join check: ", clinical_join_details)
      if (join_overlap / nrow(aligned_meta) >= 0.90) {
        clinical_join_df <- clinical_df |>
          dplyr::distinct(.data$individualID, .keep_all = TRUE)
        aligned_meta <- aligned_meta |>
          dplyr::left_join(
            clinical_join_df,
            by = c("IndividualID" = "individualID"),
            suffix = c("", "_clinical")
          )
        rownames(aligned_meta) <- protein_sample_ids
        clinical_join_details <- paste0(clinical_join_details, "; joined")
      }
    } else if ("projid" %in% colnames(aligned_meta) && "projid" %in% colnames(clinical_df)) {
      join_overlap <- sum(unique(as.character(aligned_meta$projid)) %in% unique(as.character(clinical_df$projid)))
      clinical_join_details <- paste0("projid overlap ", join_overlap, "/", nrow(aligned_meta))
      message("Protein clinical metadata join check: ", clinical_join_details)
      if (join_overlap / nrow(aligned_meta) >= 0.90) {
        clinical_join_df <- clinical_df |>
          dplyr::distinct(.data$projid, .keep_all = TRUE)
        aligned_meta <- aligned_meta |>
          dplyr::left_join(clinical_join_df, by = "projid", suffix = c("", "_clinical"))
        rownames(aligned_meta) <- protein_sample_ids
        clinical_join_details <- paste0(clinical_join_details, "; joined")
      }
    } else {
      message("Protein clinical metadata loaded, but no validated join key was available.")
    }
  }

  final_stage_col <- pick_present_col(aligned_meta, protein_stage_candidates)
  final_covariate_tbl <- purrr::imap_dfr(protein_covariate_candidates, function(candidates, covariate_name) {
    detected_col <- pick_present_col(aligned_meta, candidates)
    tibble::tibble(
      covariate = covariate_name,
      detected_column = detected_col,
      n_non_missing = if (!is.na(detected_col)) sum(!is.na(aligned_meta[[detected_col]])) else NA_integer_
    )
  })

  message("Detected protein stage column: ", final_stage_col)
  message("Detected protein covariate columns and non-missing counts:")
  print(final_covariate_tbl)

  inputs$protein_meta <<- aligned_meta
  inputs$prot_meta_aligned <<- aligned_meta
  inputs$protein_meta_alignment_audit <<- overlap_tbl
  inputs$protein_meta_selected_alignment <<- tibble::tibble(
    selected_metadata_object = selected_name,
    selected_id_column = selected_id_col,
    overlap_count = best$overlap_count,
    percent_overlap = best$percent_overlap,
    protein_matrix_name = protein_info$matrix_name,
    protein_matrix_orientation = protein_info$orientation,
    n_protein_samples = length(protein_sample_ids),
    n_metadata_rows_before_alignment = nrow(selected_df),
    n_metadata_rows_after_alignment = nrow(aligned_meta),
    clinical_join_details = clinical_join_details
  )

  record_input_source(
    "protein_meta",
    "aligned_metadata",
    selected_name,
    "loaded",
    paste0("Aligned by ", selected_id_col, "; overlap ", best$overlap_count, "/", length(protein_sample_ids))
  )

  invisible(TRUE)
}

align_protein_metadata()

# ============================================================
# Final diagnostics
# ============================================================

message("Loaded input names:")
if (length(inputs) == 0) {
  message("  <none>")
} else {
  message("  ", paste(names(inputs), collapse = ", "))
}

message("Input source audit dimensions: ", paste(dim(input_source_audit), collapse = " x "))
