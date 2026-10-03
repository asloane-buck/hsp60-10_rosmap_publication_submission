# ============================================================
# Supplementary Figure 2: matched-individual RNA/protein sensitivity
# ============================================================

# This script replaces the previous robustness-dashboard design.
# Purpose: show RNA and protein remodeling trajectories in the same
# matched-individual subset while retaining modality-specific inference.
# Direct RNA-vs-protein significance testing is intentionally avoided because
# the modality-specific standardized scores are not interpreted as a direct
# cross-modal effect difference.

if (!exists("inputs")) {
  stop("Run 01_supplemental_load_inputs.R before 11_make_supplementary_figure_2_matched_individual_sensitivity.R.", call. = FALSE)
}
if (!exists("pick_col")) {
  stop("Run 02_supplemental_helper_functions.R before 11_make_supplementary_figure_2_matched_individual_sensitivity.R.", call. = FALSE)
}

# ============================================================
# Main Figure 1 style settings
# ============================================================

fig1_pathway_colors <- c(
  Hsp60_10_all_clients = "#B23A48",
  Broad_MitoCarta_non_Hsp60_10 = "#4C78A8",
  Mito_translation_non_Hsp60_10 = "#7A3E9D",
  OXPHOS_ETC_non_Hsp60_10 = "#2E86DE",
  TCA_pyruvate_metabolism_non_Hsp60_10 = "#2CA02C",
  FAO_metabolism_non_Hsp60_10 = "#E67E22",
  Mito_protein_quality_control_non_Hsp60_10 = "#A65E2E",
  Proteasome_core_non_Hsp60_10 = "#7D3C98",
  Lysosome_core_non_Hsp60_10 = "#1E8449",
  UPRmt_core_non_Hsp60_10 = "#F39C12"
)

modality_colors <- c(
  "Protein" = "#8E1B1B",
  "RNA" = "#1B4F9C"
)

if (!requireNamespace("ggh4x", quietly = TRUE)) {
  message("Package ggh4x is not installed. Facet strip text colors will fall back to default black. Run install.packages('ggh4x') if you want pathway-colored facet labels.")
}

# ============================================================
# Local helpers
# ============================================================

sf2_get_input <- function(candidate_names, required = FALSE, label = "input") {
  selected_name <- candidate_names[candidate_names %in% names(inputs)][1]

  if (!is.na(selected_name)) {
    message("Detected ", label, ": ", selected_name)
    return(list(name = selected_name, object = inputs[[selected_name]]))
  }

  if (isTRUE(required)) {
    stop(
      "Missing required input for ", label, ". Candidate object names tried: ",
      paste(candidate_names, collapse = ", "),
      ". Available input names: ",
      paste(names(inputs), collapse = ", "),
      call. = FALSE
    )
  }

  message("Optional input unavailable for ", label, ". Tried: ", paste(candidate_names, collapse = ", "))
  list(name = NA_character_, object = NULL)
}

sf2_as_gene_table <- function(object, label) {
  if (is.null(object)) {
    return(tibble::tibble(gene = character()))
  }

  if (is.atomic(object) && !is.matrix(object)) {
    out <- tibble::tibble(gene = clean_gene(object))
  } else {
    df <- as.data.frame(object)
    gene_col <- pick_col(
      df,
      c(
        "gene", "gene_symbol", "symbol", "hgnc_symbol",
        "curated_gene", "display_gene", "Gene", "SYMBOL", "GENE"
      ),
      required = TRUE,
      label = paste(label, "gene symbol")
    )
    out <- tibble::tibble(gene = clean_gene(df[[gene_col]]))
  }

  out |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::distinct(.data$gene)
}

sf2_metadata_sample_ids <- function(meta, label) {
  if (is.null(meta)) {
    return(character())
  }

  df <- as.data.frame(meta)
  sample_col <- pick_col(
    df,
    c(
      "sample_id", "SampleID", "sampleid", "sample", "sample_name",
      "specimenID", "specimen_id", "SpecimenID", "projid", "individualID",
      "IndividualID", "individual_id", "rna_id", "protein_id", "batch.channel"
    ),
    required = FALSE,
    label = paste(label, "sample identifier")
  )

  if (!is.na(sample_col)) {
    message("Detected sample column for ", label, ": ", sample_col)
    return(as.character(df[[sample_col]]))
  }

  if (!is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
    message("Using row names as sample identifiers for ", label)
    return(rownames(df))
  }

  message("No sample identifier detected for ", label)
  character()
}

sf2_prepare_matrix_input <- function(matrix_candidates, meta_candidates, interactor_genes, label) {
  matrix_input <- sf2_get_input(matrix_candidates, required = FALSE, label = paste(label, "matrix"))
  meta_input <- sf2_get_input(meta_candidates, required = FALSE, label = paste(label, "metadata"))
  sample_ids <- sf2_metadata_sample_ids(meta_input$object, label)

  oriented <- tryCatch(
    infer_and_orient_matrix(
      object = matrix_input$object,
      gene_symbols = interactor_genes,
      sample_ids = sample_ids,
      label = label,
      required = FALSE
    ),
    error = function(err) {
      message("Could not prepare ", label, ": ", conditionMessage(err))
      list(
        matrix = NULL,
        diagnostics = tibble::tibble(
          matrix_label = label,
          status = "unavailable",
          reason = conditionMessage(err)
        )
      )
    }
  )

  list(
    matrix_name = matrix_input$name,
    meta_name = meta_input$name,
    matrix = oriented$matrix,
    orientation_diagnostics = oriented$diagnostics %||% tibble::tibble(),
    metadata = meta_input$object,
    sample_ids = sample_ids
  )
}

sf2_first_present_col <- function(df, candidates) {
  selected <- candidates[candidates %in% colnames(df)]
  if (length(selected) == 0) NA_character_ else selected[1]
}

sf2_stage_candidates <- c(
  "clinical_stage",
  "cogdx_num",
  "cogdx",
  "diagnosis_stage",
  "diagnosis"
)

sf2_individual_id_candidates <- c(
  "IndividualID", "individualID", "individual_id", "projid", "proj_id",
  "projid.x", "projid.y", "subject_id", "SubjectID", "id"
)

sf2_sample_id_candidates <- c(
  "SampleID", "sample_id", "sampleid", "sample", "sample_name",
  "SpecimenID", "specimenID", "specimen_id", "batch.channel",
  "rna_id", "protein_id", "id"
)

sf2_harmonize_stage <- function(stage) {
  stage_chr <- stringr::str_squish(as.character(stage))
  stage_clean <- stringr::str_to_lower(stage_chr)

  dplyr::case_when(
    stage_clean %in% c("1", "nci") ~ "NCI",
    stage_clean %in% c("2", "mci") ~ "MCI",
    stage_clean %in% c("4", "ad") ~ "AD",
    TRUE ~ NA_character_
  )
}

sf2_build_sample_meta <- function(meta, matrix_sample_ids, modality) {
  if (is.null(meta)) {
    stop(modality, " metadata is unavailable; cannot build matched-individual sensitivity figure.", call. = FALSE)
  }

  df <- as.data.frame(meta)

  sample_col <- sf2_first_present_col(df, sf2_sample_id_candidates)
  id_col <- sf2_first_present_col(df, sf2_individual_id_candidates)
  stage_col <- sf2_first_present_col(df, sf2_stage_candidates)

  # If sample IDs are stored as rownames, expose them as a usable column.
  if (is.na(sample_col) && !is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
    df$.sample_id_from_rownames <- rownames(df)
    sample_col <- ".sample_id_from_rownames"
  }

  if (is.na(id_col)) {
    stop(
      modality, " metadata is missing an individual ID column. Tried: ",
      paste(sf2_individual_id_candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(stage_col)) {
    stop(
      modality, " metadata is missing a diagnosis/stage column. Tried: ",
      paste(sf2_stage_candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(sample_col)) {
    stop(
      modality, " metadata is missing a sample ID column/rownames. Tried: ",
      paste(sf2_sample_id_candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(df) |>
    dplyr::transmute(
      sample_id = as.character(.data[[sample_col]]),
      individual_id = as.character(.data[[id_col]]),
      original_stage = as.character(.data[[stage_col]]),
      Stage = sf2_harmonize_stage(.data[[stage_col]]),
      Modality = modality
    ) |>
    dplyr::filter(.data$sample_id %in% matrix_sample_ids) |>
    dplyr::mutate(
      sample_id = factor(.data$sample_id, levels = matrix_sample_ids)
    ) |>
    dplyr::arrange(.data$sample_id) |>
    dplyr::mutate(sample_id = as.character(.data$sample_id)) |>
    dplyr::filter(
      !is.na(.data$individual_id),
      nzchar(.data$individual_id),
      !is.na(.data$Stage)
    )

  if (anyDuplicated(out$sample_id) > 0) {
    stop(modality, " metadata contain duplicate matrix sample IDs.", call. = FALSE)
  }

  if (anyDuplicated(out$individual_id) > 0) {
    stop(
      modality,
      " metadata contain more than one primary-stage sample per participant.",
      call. = FALSE
    )
  }

  message(modality, " metadata columns:")
  message("  sample column: ", sample_col)
  message("  individual column: ", id_col)
  message("  stage column: ", stage_col)
  message("  rows matched to matrix samples: ", nrow(out), "/", length(matrix_sample_ids))

  attr(out, "sample_col") <- sample_col
  attr(out, "id_col") <- id_col
  attr(out, "stage_col") <- stage_col
  out
}

sf2_score_pathway_mean_z <- function(expr_mat, genes) {
  genes_use <- intersect(clean_gene(genes), rownames(expr_mat))
  if (length(genes_use) < 3) {
    return(rep(NA_real_, ncol(expr_mat)))
  }

  mat <- expr_mat[genes_use, , drop = FALSE]
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- NA_real_
  score <- colMeans(z, na.rm = TRUE)
  score[is.nan(score)] <- NA_real_
  score
}

sf2_sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

sf2_p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

sf2_safe_unpaired_wilcox <- function(df, value_col, group_col, g1, g2) {
  sub <- df |>
    dplyr::filter(.data[[group_col]] %in% c(g1, g2)) |>
    dplyr::select(dplyr::all_of(c(value_col, group_col))) |>
    tidyr::drop_na()
  x <- sub |> dplyr::filter(.data[[group_col]] == g1) |> dplyr::pull(.data[[value_col]])
  y <- sub |> dplyr::filter(.data[[group_col]] == g2) |> dplyr::pull(.data[[value_col]])

  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}


sf2_make_pathway_sets <- function(interactor_tbl, mito_tbl) {
  interactor_genes <- interactor_tbl$gene
  
  mito_genes <- if (!is.null(mito_tbl) && nrow(mito_tbl) > 0) {
    setdiff(mito_tbl$gene, interactor_genes)
  } else {
    character()
  }
  
  tca_input <- sf2_get_input(
    c(
      "tca_pyruvate_tbl",
      "TCA_pyruvate_metabolism_non_Hsp60_10",
      "tca_pyruvate_metabolism_non_hsp60_10"
    ),
    required = TRUE,
    label = "TCA/pyruvate Figure 1 pathway table"
  )
  
  tca_tbl <- sf2_as_gene_table(
    tca_input$object,
    "TCA/pyruvate Figure 1 pathway table"
  )
  
  tca_pyruvate_genes <- setdiff(tca_tbl$gene, interactor_genes)
  
  tibble::tibble(
    Pathway = c(
      "Hsp60/10 client network",
      "Broad mitochondrial background",
      "TCA/pyruvate metabolism"
    ),
    facet_label = c(
      "Hsp60/10\nclients",
      "Mito\nbackground",
      "TCA/pyruvate\nnon-client"
    ),
    pathway_color_key = c(
      "Hsp60_10_all_clients",
      "Broad_MitoCarta_non_Hsp60_10",
      "TCA_pyruvate_metabolism_non_Hsp60_10"
    ),
    genes = list(
      interactor_genes,
      mito_genes,
      tca_pyruvate_genes
    )
  ) |>
    dplyr::mutate(n_genes_defined = purrr::map_int(.data$genes, length))
}

sf2_make_paper_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey20", linewidth = 0.55),
      axis.line = ggplot2::element_line(color = "grey20", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "grey20", linewidth = 0.35),
      axis.text = ggplot2::element_text(size = base_size * 0.85, color = "grey15"),
      axis.title = ggplot2::element_text(size = base_size * 0.95, color = "grey10", face = "bold"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold", size = base_size * 0.95),
      legend.text = ggplot2::element_text(size = base_size * 0.9),
      strip.background = ggplot2::element_rect(fill = "grey94", color = "grey25", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "plain", size = base_size * 0.82),
      panel.spacing = grid::unit(0.55, "lines")
    )
}

# ============================================================
# Main builder
# ============================================================

make_supfig2_matched_individual_sensitivity <- function(inputs) {
  message("Starting Supplementary Figure 2: matched-individual RNA/protein sensitivity.")

  interactor_input <- sf2_get_input(
    c(
      "hsp60_hsp10_interactor_inventory",
      "hsp60_hsp10_client_inventory",
      "all_hsp60_10_client_tbl",
      "hsp60_client_tbl",
      "hsp_client_tbl",
      "hsp_clients_all",
      "hsp60_clients"
    ),
    required = TRUE,
    label = "Hsp60/10 client inventory table"
  )
  interactor_tbl <- sf2_as_gene_table(interactor_input$object, "Hsp60/10 client inventory table")

  mito_input <- sf2_get_input(
    c("mito_background_tbl", "mitochondrial_background_tbl"),
    required = FALSE,
    label = "mitochondrial background table"
  )
  mito_tbl <- if (!is.null(mito_input$object)) {
    sf2_as_gene_table(mito_input$object, "mitochondrial background table")
  } else {
    NULL
  }

  pathway_tbl <- sf2_make_pathway_sets(interactor_tbl, mito_tbl)

  rna_info <- sf2_prepare_matrix_input(
    matrix_candidates = c("rna_mat", "vst_mat", "rna_vst", "rna_matrix", "vst_symbol_mat", "rna_mat_stage"),
    meta_candidates = c("rna_meta", "meta_rna_sub", "rna_metadata", "rna_meta_in", "meta_rna"),
    interactor_genes = interactor_tbl$gene,
    label = "RNA matrix"
  )

  protein_info <- sf2_prepare_matrix_input(
    matrix_candidates = c("prot_mat", "protein_matrix", "protein_mat", "tmt_mat", "prot_expr_mat", "prot_mat_in", "prot_mat_raw"),
    meta_candidates = c("protein_meta", "prot_meta_aligned", "prot_meta_path", "prot_meta", "meta_protein_sub", "prot_meta_in"),
    interactor_genes = interactor_tbl$gene,
    label = "Protein matrix"
  )

  if (is.null(rna_info$matrix)) {
    stop("RNA matrix is unavailable; cannot make Supplementary Figure 2.", call. = FALSE)
  }
  if (is.null(protein_info$matrix)) {
    stop("Protein matrix is unavailable; cannot make Supplementary Figure 2.", call. = FALSE)
  }

  rna_meta <- sf2_build_sample_meta(rna_info$metadata, colnames(rna_info$matrix), "RNA")
  protein_meta <- sf2_build_sample_meta(protein_info$metadata, colnames(protein_info$matrix), "Protein")

  matched_stage_crosswalk <- protein_meta |>
    dplyr::select(
      individual_id,
      protein_sample_id = sample_id,
      protein_stage = Stage
    ) |>
    dplyr::inner_join(
      rna_meta |>
        dplyr::select(
          individual_id,
          rna_sample_id = sample_id,
          rna_stage = Stage
        ),
      by = "individual_id"
    ) |>
    dplyr::mutate(
      stage_agrees = as.character(.data$protein_stage) ==
        as.character(.data$rna_stage)
    )

  if (any(!matched_stage_crosswalk$stage_agrees)) {
    bad <- matched_stage_crosswalk |>
      dplyr::filter(!.data$stage_agrees)

    readr::write_csv(
      bad,
      file.path(audits_dir, "SuppFig2_STAGE_DISAGREEMENTS_ERROR.csv")
    )

    stop(
      "RNA and protein canonical clinical stages disagree for ",
      nrow(bad),
      " matched participants. See SuppFig2_STAGE_DISAGREEMENTS_ERROR.csv.",
      call. = FALSE
    )
  }

  overlap_ids <- matched_stage_crosswalk$individual_id

  if (length(overlap_ids) < 10) {
    stop(
      "Too few primary-stage matched individuals found across RNA and protein: ",
      length(overlap_ids),
      ". Check individual ID/stage columns.",
      call. = FALSE
    )
  }

  message(
    "Matched individuals with RNA and protein and identical canonical stage: ",
    length(overlap_ids)
  )

  readr::write_csv(
    matched_stage_crosswalk,
    file.path(audits_dir, "SuppFig2_matched_stage_crosswalk.csv")
  )

  readr::write_csv(
    dplyr::bind_rows(rna_meta, protein_meta) |>
      dplyr::mutate(in_matched_subset = .data$individual_id %in% overlap_ids),
    file.path(audits_dir, "SuppFig2_matched_sample_availability_audit.csv")
  )

  make_modality_long <- function(expr_mat, sample_meta, modality) {
    purrr::pmap_dfr(
      pathway_tbl,
      function(Pathway, facet_label, pathway_color_key, genes, n_genes_defined) {
        genes_found <- intersect(clean_gene(genes), rownames(expr_mat))
        score <- sf2_score_pathway_mean_z(expr_mat, genes_found)

        tibble::tibble(
          sample_id = colnames(expr_mat),
          Score = as.numeric(score),
          Pathway = Pathway,
          facet_label = facet_label,
          pathway_color_key = pathway_color_key,
          n_genes_defined = n_genes_defined,
          n_genes_detected = length(genes_found),
          Modality = modality
        ) |>
          dplyr::left_join(sample_meta, by = c("sample_id", "Modality"))
      }
    )
  }

  matched_long <- dplyr::bind_rows(
    make_modality_long(rna_info$matrix, rna_meta, "RNA"),
    make_modality_long(protein_info$matrix, protein_meta, "Protein")
  ) |>
    dplyr::filter(
      .data$individual_id %in% overlap_ids,
      !is.na(.data$Stage),
      !is.na(.data$Score)
    ) |>
    dplyr::mutate(
      Stage = factor(.data$Stage, levels = c("NCI", "MCI", "AD")),
      Modality = factor(.data$Modality, levels = c("Protein", "RNA")),
      facet_label = factor(.data$facet_label, levels = pathway_tbl$facet_label),
      x = as.numeric(.data$Stage)
    )

  readr::write_csv(matched_long, file.path(audits_dir, "SuppFig2_matched_individual_long_scores.csv"))

  matched_summary <- matched_long |>
    dplyr::group_by(.data$Pathway, .data$facet_label, .data$pathway_color_key, .data$Modality, .data$Stage) |>
    dplyr::summarise(
      mean_score = mean(.data$Score, na.rm = TRUE),
      sem = sf2_sem(.data$Score),
      n = sum(!is.na(.data$Score)),
      n_genes_detected = dplyr::first(.data$n_genes_detected),
      .groups = "drop"
    ) |>
    dplyr::mutate(x = as.numeric(.data$Stage))

  readr::write_csv(matched_summary, file.path(audits_dir, "SuppFig2_matched_individual_summary.csv"))

  within_sig <- matched_long |>
    dplyr::group_by(.data$Pathway, .data$facet_label, .data$Modality) |>
    dplyr::group_modify(~{
      d <- .x
      key <- .y
      summ <- matched_summary |>
        dplyr::filter(.data$Pathway == key$Pathway[[1]], .data$Modality == key$Modality[[1]]) |>
        dplyr::arrange(.data$x)

      all_top <- max(summ$mean_score + summ$sem, na.rm = TRUE)
      all_bot <- min(summ$mean_score - summ$sem, na.rm = TRUE)
      rng <- all_top - all_bot
      if (!is.finite(rng) || rng == 0) rng <- 0.12

      y1_top <- max(summ$mean_score[summ$x %in% c(1, 2)] + summ$sem[summ$x %in% c(1, 2)], na.rm = TRUE)
      y2_top <- max(summ$mean_score[summ$x %in% c(2, 3)] + summ$sem[summ$x %in% c(2, 3)], na.rm = TRUE)
      x_shift <- if (as.character(key$Modality[[1]]) == "Protein") -0.10 else 0.10

      p_early <- sf2_safe_unpaired_wilcox(d, "Score", "Stage", "NCI", "MCI")
      p_late <- sf2_safe_unpaired_wilcox(d, "Score", "Stage", "MCI", "AD")

      tibble::tibble(
        comparison = c("NCI_vs_MCI", "MCI_vs_AD"),
        x1 = c(1, 2) + x_shift,
        x2 = c(2, 3) + x_shift,
        y = if (as.character(key$Modality[[1]]) == "Protein") {
          c(y1_top + 0.08 * rng, y2_top + 0.20 * rng)
        } else {
          c(y1_top + 0.18 * rng, y2_top + 0.10 * rng)
        },
        p_value_raw = c(p_early, p_late),
        label = sf2_p_to_stars(c(p_early, p_late))
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$label != "")

  readr::write_csv(
    within_sig,
    file.path(
      audits_dir,
      "SuppFig2_within_modality_significance.csv"
    )
  )

  plot_range <- range(
    c(
      matched_summary$mean_score - matched_summary$sem,
      matched_summary$mean_score + matched_summary$sem,
      within_sig$y
    ),
    na.rm = TRUE
  )
  plot_pad <- diff(plot_range) * 0.08
  if (!is.finite(plot_pad) || plot_pad == 0) plot_pad <- 0.08
  
  facet_color_tbl <- matched_summary |>
    dplyr::distinct(.data$facet_label, .data$pathway_color_key) |>
    dplyr::arrange(factor(.data$facet_label, levels = pathway_tbl$facet_label))
  
  strip_fill_values <- fig1_pathway_colors[facet_color_tbl$pathway_color_key]
  strip_fill_values[is.na(strip_fill_values)] <- "grey90"
  
  p_matched <- ggplot2::ggplot(
    matched_summary,
    ggplot2::aes(
      x = .data$x,
      y = .data$mean_score,
      color = .data$Modality,
      group = .data$Modality
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$mean_score - .data$sem, ymax = .data$mean_score + .data$sem),
      width = 0.06,
      linewidth = 0.8
    ) +
    ggplot2::geom_segment(
      data = within_sig,
      ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y, color = .data$Modality),
      inherit.aes = FALSE,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = within_sig,
      ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y - 0.012, yend = .data$y, color = .data$Modality),
      inherit.aes = FALSE,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = within_sig,
      ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y - 0.012, yend = .data$y, color = .data$Modality),
      inherit.aes = FALSE,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = within_sig,
      ggplot2::aes(x = (.data$x1 + .data$x2) / 2, y = .data$y + 0.010, label = .data$label, color = .data$Modality),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.9,
      show.legend = FALSE
    ) +
    {
      if (requireNamespace("ggh4x", quietly = TRUE)) {
        ggh4x::facet_wrap2(
          ~ facet_label,
          ncol = 3,
          strip = ggh4x::strip_themed(
            background_x = lapply(
              strip_fill_values,
              function(col_i) ggplot2::element_rect(
                fill = "grey94",
                color = "grey80",
                linewidth = 0.45
              )
            ),
            text_x = lapply(
              strip_fill_values,
              function(col_i) ggplot2::element_text(
                face = "plain",
                color = col_i,
                size = 7.2
              )
            )
          )
        )
      } else {
        ggplot2::facet_wrap(~ facet_label, ncol = 3)
      }
    } +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2, 3),
      labels = c("NCI", "MCI", "AD")
    ) +
    ggplot2::scale_color_manual(
      values = modality_colors[c("Protein", "RNA")],
      breaks = c("Protein", "RNA"),
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(plot_range[1] - plot_pad, plot_range[2] + plot_pad),
      clip = "off"
    ) +
    ggplot2::labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = "Mean pathway score",
      color = "Modality"
    ) +
    sf2_make_paper_theme(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin = ggplot2::margin(8, 16, 8, 8)
    )

  save_plot_set(
    p_matched,
    "Supplementary_Figure_2_matched_individual_sensitivity",
    width = 13.2,
    height = 6.6,
    output_dir = figures_dir
  )

  obsolete_paths <- file.path(
    figures_dir,
    paste0("Supplementary_Figure_2_protein_remodeling_robustness", c(".pdf", ".png"))
  )
  obsolete_paths <- obsolete_paths[file.exists(obsolete_paths)]
  if (length(obsolete_paths) > 0) {
    unlink(obsolete_paths)
    message("Removed obsolete duplicate Supplementary Figure 2 files: ", paste(obsolete_paths, collapse = ", "))
  }

  save_panel_set(
    p_matched,
    "SuppFig2_matched_individual_sensitivity",
    width = 13.2,
    height = 6.6,
    output_dir = panels_dir
  )

  run_summary <- c(
    "Supplementary Figure 2 run summary",
    paste0("generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "figure: matched-individual RNA/protein sensitivity analysis",
    paste0("matched primary-stage individuals with identical RNA/protein stage: ", length(overlap_ids)),
    paste0("RNA matrix: ", nrow(rna_info$matrix), " genes x ", ncol(rna_info$matrix), " samples"),
    paste0("Protein matrix: ", nrow(protein_info$matrix), " genes x ", ncol(protein_info$matrix), " samples"),
    paste0("RNA metadata input: ", rna_info$meta_name),
    paste0("Protein metadata input: ", protein_info$meta_name),
    paste0("Client inventory size after harmonization: ", nrow(interactor_tbl)),
    paste0("Mitochondrial background size after harmonization: ", ifelse(is.null(mito_tbl), 0, nrow(mito_tbl))),
    "pathways plotted:",
    paste0("  - ", pathway_tbl$Pathway, " (defined n=", pathway_tbl$n_genes_defined, ")", collapse = "\n")
  )
  writeLines(run_summary, file.path(audits_dir, "SuppFig2_run_summary.txt"))

  message("Supplementary Figure 2 complete.")
  invisible(list(
    plot = p_matched,
    matched_long = matched_long,
    matched_summary = matched_summary,
    within_sig = within_sig,
    pathway_tbl = pathway_tbl,
    rna_meta = rna_meta,
    protein_meta = protein_meta,
    matched_stage_crosswalk = matched_stage_crosswalk,
    overlap_ids = overlap_ids
  ))
}

supfig2_outputs <- make_supfig2_matched_individual_sensitivity(inputs)
