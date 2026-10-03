# ============================================================
# Supplementary Figure 1: cohort, detection, and client inventory
# ============================================================

if (!exists("inputs")) {
  stop("Run 01_supplemental_load_inputs.R before 10_make_supplementary_figure_1_cohort_detection.R.", call. = FALSE)
}

if (!exists("pick_col")) {
  stop("Run 02_supplemental_helper_functions.R before 10_make_supplementary_figure_1_cohort_detection.R.", call. = FALSE)
}


# ============================================================
# Main-figure style settings
# ============================================================

# Match the primary manuscript convention:
# Protein = red, RNA = blue.
SUP_PROTEIN_COL <- "#B2182B"
SUP_RNA_COL <- "#1B4F9C"
SUP_PROTEIN_LIGHT <- "#F6D7D7"
SUP_RNA_LIGHT <- "#DCEAF7"
SUP_NEUTRAL_DARK <- "#2F3A45"
SUP_NEUTRAL_MID <- "#8A97A3"
SUP_NEUTRAL_LIGHT <- "#EEF0F2"

theme_mainfigure_panel <- function(base_size = 8.5) {
  theme_supplement(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", color = "#111827"),
      axis.text = ggplot2::element_text(color = "#1F2937"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.25)
    )
}


# ============================================================
# Input selection helpers
# ============================================================

get_input <- function(candidate_names, required = FALSE, label = "input") {
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

as_gene_table <- function(object, label) {
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

metadata_sample_ids <- function(meta, label) {
  if (is.null(meta)) {
    return(character())
  }

  df <- as.data.frame(meta)
  sample_col <- pick_col(
    df,
    c(
      "sample_id", "SampleID", "sampleid", "sample", "sample_name",
      "specimenID", "specimen_id", "projid", "individualID",
      "individual_id", "rna_id", "protein_id"
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

prepare_matrix_input <- function(matrix_candidates, meta_candidates, interactor_genes, label) {
  matrix_input <- get_input(matrix_candidates, required = FALSE, label = paste(label, "matrix"))
  meta_input <- get_input(meta_candidates, required = FALSE, label = paste(label, "metadata"))
  sample_ids <- metadata_sample_ids(meta_input$object, label)

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

unavailable_plot <- function(title, detail) {
  ggplot2::ggplot(tibble::tibble(x = 0, y = 0), ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_label(
      ggplot2::aes(label = detail),
      size = 4,
      label.size = 0.25,
      fill = "white",
      color = "grey20"
    ) +
    ggplot2::labs(title = title) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    theme_supplement(base_size = 8.5) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ============================================================
# Panel A: workflow schematic
# ============================================================

make_workflow_panel <- function(inputs, rna_info, protein_info, n_interactors) {
  nodes <- tibble::tibble(
    node = c("inventory", "rna", "protein", "matched", "pathology", "validation"),
    label = c(
      "Client\ninventory",
      "ROSMAP\nRNA-seq",
      "ROSMAP\nproteomics",
      "Matched\nRNA/protein",
      "Pathology\nanalyses",
      "External/regional\nvalidation"
    ),
    x = c(1.05, 3.05, 5.05, 7.05, 9.00, 11.60),
    y = 1,
    fill_group = c("inventory", "rna", "protein", "matched", "analysis", "validation")
  ) |>
    dplyr::mutate(
      node_width = dplyr::if_else(.data$node == "validation", 1.00, 0.80),
      node_height = 0.34
    )

  edges <- tibble::tibble(
    x = nodes$x[-nrow(nodes)] + nodes$node_width[-nrow(nodes)] + 0.16,
    xend = nodes$x[-1] - nodes$node_width[-1] - 0.16,
    y = 1,
    yend = 1
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      arrow = grid::arrow(length = grid::unit(0.095, "inches"), type = "closed"),
      linewidth = 0.45,
      color = "grey35"
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(
        xmin = .data$x - .data$node_width,
        xmax = .data$x + .data$node_width,
        ymin = .data$y - .data$node_height,
        ymax = .data$y + .data$node_height,
        fill = .data$fill_group
      ),
      color = "#374151",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.45,
      lineheight = 0.90,
      color = "#111827"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        inventory = "#F3F4F6",
        rna = SUP_RNA_LIGHT,
        protein = SUP_PROTEIN_LIGHT,
        matched = "#F5F5F5",
        analysis = "#F7F1E8",
        validation = "#F3F4F6"
      ),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 13.05), ylim = c(0.34, 1.66), clip = "off") +
    ggplot2::labs(title = "Analysis workflow") +
    theme_mainfigure_panel(base_size = 8.5) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 20, 6, 20)
    )
}

# ============================================================
# Panel B: client detection counts
# ============================================================

make_detection_panel <- function(interactor_tbl, rna_info, protein_info) {
  rna_available <- !is.null(rna_info$matrix)
  protein_available <- !is.null(protein_info$matrix)

  rna_genes <- if (rna_available) rownames(rna_info$matrix) else character()
  protein_genes <- if (protein_available) rownames(protein_info$matrix) else character()

  detection_audit <- interactor_tbl |>
    dplyr::mutate(
      source_inventory = "Bie et al. 2020 Hsp60/10 interactor inventory",
      original_reported_inventory = 323L,
      analyzed_after_gene_symbol_harmonization = nrow(interactor_tbl),
      rna_input_available = rna_available,
      protein_input_available = protein_available,
      detected_rna = if (rna_available) .data$gene %in% rna_genes else NA,
      detected_protein = if (protein_available) .data$gene %in% protein_genes else NA,
      detection_group = dplyr::case_when(
        rna_available & protein_available & .data$detected_rna & .data$detected_protein ~ "Both RNA and protein",
        rna_available & protein_available & .data$detected_rna & !.data$detected_protein ~ "RNA only",
        rna_available & protein_available & !.data$detected_rna & .data$detected_protein ~ "Protein only",
        rna_available & protein_available & !.data$detected_rna & !.data$detected_protein ~ "Neither",
        rna_available & .data$detected_rna ~ "Detected in RNA",
        rna_available & !.data$detected_rna ~ "Not detected in RNA",
        protein_available & .data$detected_protein ~ "Detected in protein",
        protein_available & !.data$detected_protein ~ "Not detected in protein",
        TRUE ~ "No expression matrix loaded"
      )
    )

  detection_audit_paths <- file.path(audits_dir, "SuppFig1_interactor_detection_audit.csv")
  purrr::walk(detection_audit_paths, ~ readr::write_csv(detection_audit, .x))

  detection_counts <- detection_audit |>
    dplyr::count(.data$detection_group, name = "n_interactors") |>
    dplyr::mutate(
      display_group = dplyr::recode(
        as.character(.data$detection_group),
        "Both RNA and protein" = "RNA + protein",
        "Neither" = "Not detected",
        .default = as.character(.data$detection_group)
      ),
      detection_group = factor(
        .data$detection_group,
        levels = c(
          "RNA only",
          "Protein only",
          "Both RNA and protein",
          "Neither",
          "Detected in RNA",
          "Not detected in RNA",
          "Detected in protein",
          "Not detected in protein",
          "No expression matrix loaded"
        )
      )
    ) |>
    dplyr::arrange(.data$detection_group)

  readr::write_csv(
    detection_counts,
    file.path(tables_dir, "SuppFig1_detection_counts.csv")
  )

  message("Source inventory: Bie et al. 2020 Hsp60/10 interactor inventory")
  message("Original reported inventory: 323 interactors")
  message("Analyzed after gene-symbol harmonization: ", nrow(interactor_tbl), " interactors")
  if (rna_available) {
    message("Hsp60/10 clients detected in RNA: ", sum(detection_audit$detected_rna, na.rm = TRUE))
  }
  if (protein_available) {
    message("Hsp60/10 clients detected in protein: ", sum(detection_audit$detected_protein, na.rm = TRUE))
  }

  plot_counts <- detection_counts |>
    dplyr::mutate(
      plot_group = dplyr::recode(
        as.character(.data$display_group),
        "RNA only" = "RNA",
        "Protein only" = "Protein",
        "RNA + protein" = "Both",
        "Not detected" = "Neither"
      ),
      plot_group = factor(
        .data$plot_group,
        levels = c("RNA", "Protein", "Both", "Neither")
      )
    )

  p <- ggplot2::ggplot(
    plot_counts,
    ggplot2::aes(.data$plot_group, .data$n_interactors, fill = .data$plot_group)
  ) +
    ggplot2::geom_col(width = 0.62, color = "grey25", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$n_interactors),
      vjust = -0.35,
      size = 3.4,
      color = "#111827"
    ) +
    ggplot2::labs(
      title = "Client detection by modality",
      subtitle = paste0("Published Hsp60/10 client inventory, n = ", nrow(interactor_tbl), "."),
      x = NULL,
      y = "Number of clients"
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(
      values = c(
        "RNA" = SUP_RNA_COL,
        "Protein" = SUP_PROTEIN_COL,
        "Both" = "#6B4C7A",
        "Neither" = SUP_NEUTRAL_MID
      ),
      guide = "none",
      drop = FALSE
    ) +
    theme_mainfigure_panel(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = ggplot2::element_blank()
    )

  attr(p, "detection_counts") <- detection_counts
  attr(p, "detection_audit") <- detection_audit
  p
}

# ============================================================
# Panel C: metadata summaries
# ============================================================

stage_column_candidates <- c(
  "clinical_stage",
  "cogdx_num",
  "cogdx",
  "diagnosis_stage",
  "diagnosis"
)

covariate_candidates <- list(
  age = c("age", "age_at_death", "age_death", "age_c", "Age", "ageAtDeath"),
  sex = c("sex", "msex", "sex_label", "Sex", "gender"),
  RIN = c("rin", "RIN", "rin_c", "RINvalue", "rna_integrity_number"),
  PMI = c("pmi", "PMI", "pmi_c", "postmortem_interval", "postmortem.interval"),
  Braak = c("braak", "Braak", "braak_num", "braaksc", "braak_stage", "braak_stage_num", "Braak.Stage"),
  CERAD = c("cerad", "CERAD", "cerad_num", "ceradsc", "cerad_score", "CERAD.Score")
)

harmonize_stage_for_display <- function(stage) {
  stage_chr <- stringr::str_squish(as.character(stage))
  stage_clean <- stringr::str_to_lower(stage_chr)

  dplyr::case_when(
    is.na(stage_chr) | stage_chr == "" | stage_clean == "missing" ~ NA_character_,
    stage_clean %in% c("1", "nci") ~ "NCI",
    stage_clean %in% c("2", "mci") ~ "MCI",
    stage_clean %in% c("4", "ad") ~ "AD",
    TRUE ~ NA_character_
  )
}

summarize_stage <- function(meta, modality) {
  if (is.null(meta)) {
    return(tibble::tibble(
      modality = modality,
      stage_column = NA_character_,
      stage = "Metadata not loaded",
      n_samples = 0
    ))
  }

  df <- as.data.frame(meta)
  stage_col <- pick_col(df, stage_column_candidates, required = FALSE, label = paste(modality, "stage"))

  if (is.na(stage_col)) {
    return(tibble::tibble(
      modality = modality,
      stage_column = NA_character_,
      stage = "Stage column not found",
      n_samples = nrow(df)
    ))
  }

  df |>
    tibble::as_tibble() |>
    dplyr::mutate(stage = as.character(.data[[stage_col]])) |>
    dplyr::mutate(stage = dplyr::if_else(is.na(.data$stage) | !nzchar(.data$stage), "Missing", .data$stage)) |>
    dplyr::count(.data$stage, name = "n_samples") |>
    dplyr::mutate(modality = modality, stage_column = stage_col, .before = 1)
}

summarize_covariates <- function(meta, modality) {
  if (is.null(meta)) {
    return(tibble::tibble(
      modality = modality,
      covariate = names(covariate_candidates),
      detected_column = NA_character_,
      column_available = FALSE,
      n_samples = 0,
      n_non_missing = NA_integer_,
      pct_non_missing = NA_real_
    ))
  }

  df <- as.data.frame(meta)
  purrr::imap_dfr(covariate_candidates, function(candidates, covariate_name) {
    detected_col <- pick_col(
      df,
      candidates,
      required = FALSE,
      label = paste(modality, covariate_name)
    )

    if (is.na(detected_col)) {
      return(tibble::tibble(
        modality = modality,
        covariate = covariate_name,
        detected_column = NA_character_,
        column_available = FALSE,
        n_samples = nrow(df),
        n_non_missing = NA_integer_,
        pct_non_missing = NA_real_
      ))
    }

    n_non_missing <- sum(!is.na(df[[detected_col]]))
    tibble::tibble(
      modality = modality,
      covariate = covariate_name,
      detected_column = detected_col,
      column_available = TRUE,
      n_samples = nrow(df),
      n_non_missing = n_non_missing,
      pct_non_missing = n_non_missing / nrow(df)
    )
  })
}

make_metadata_panels <- function(rna_meta, protein_meta) {
  stage_audit <- dplyr::bind_rows(
    summarize_stage(rna_meta, "RNA"),
    summarize_stage(protein_meta, "Protein")
  )

  covariate_audit <- dplyr::bind_rows(
    summarize_covariates(rna_meta, "RNA"),
    summarize_covariates(protein_meta, "Protein")
  )

  readr::write_csv(stage_audit, file.path(audits_dir, "SuppFig1_stage_sample_size_audit.csv"))
  readr::write_csv(covariate_audit, file.path(audits_dir, "SuppFig1_covariate_availability_audit.csv"))

  message("Stage sample size audit:")
  print(stage_audit)
  message("Covariate availability audit:")
  print(covariate_audit)

  # Metadata availability and analyzable primary sample counts are distinct.
  # Keep stage_audit above unchanged. The manuscript panel uses the same
  # frozen client-score sample counts as final Main Figure 1.
  if (isTRUE(get0("SUPP_USE_FROZEN_MODELS", ifnotfound = FALSE))) {
    primary_source <- file.path(main_outputs_dir, "tables", "main_fig1_all_clients",
      "Fig1A_all_clients_overlay_summary.csv")
    primary <- readr::read_csv(primary_source, show_col_types = FALSE)
    required <- c("Pathway", "Modality", "Stage", "n")
    if (!all(required %in% names(primary))) {
      stop("S1 final Figure 1 source lacks Pathway/Modality/Stage/n columns.", call. = FALSE)
    }
    plotted_stage_audit <- primary |>
      dplyr::filter(.data$Pathway == "Hsp60_10_all_clients",
        .data$Modality %in% c("RNA", "Protein"), .data$Stage %in% c("NCI", "MCI", "AD")) |>
      dplyr::transmute(modality = as.character(.data$Modality),
        stage = as.character(.data$Stage), n_samples = as.numeric(.data$n),
        source_file = primary_source,
        population = "Analyzable primary Hsp60/10 client-score samples")
    expected <- c(RNA_NCI = 200, RNA_MCI = 158, RNA_AD = 219,
      Protein_NCI = 167, Protein_MCI = 96, Protein_AD = 109)
    key <- paste(plotted_stage_audit$modality, plotted_stage_audit$stage, sep = "_")
    actual <- setNames(plotted_stage_audit$n_samples, key)[names(expected)]
    if (nrow(plotted_stage_audit) != 6 || anyDuplicated(key) ||
        anyNA(actual) || any(actual != expected)) {
      stop("S1 final Figure 1 source differs from verified primary analyzable counts.", call. = FALSE)
    }
    readr::write_csv(plotted_stage_audit,
      file.path(audits_dir, "SuppFig1_plotted_primary_stage_counts.csv"))
    message("S1 manuscript panel: final analyzable client-score counts; metadata availability is audited separately.")
    print(plotted_stage_audit |> dplyr::select(-"source_file"))
    usable_stage <- plotted_stage_audit |>
      dplyr::mutate(display_stage = factor(.data$stage, levels = rev(c("NCI", "MCI", "AD"))),
        modality = factor(.data$modality, levels = c("RNA", "Protein")))
  } else {
  usable_stage <- stage_audit |>
    dplyr::filter(!.data$stage %in% c("Metadata not loaded", "Stage column not found")) |>
    dplyr::mutate(
      display_stage = harmonize_stage_for_display(.data$stage),
      display_stage = factor(.data$display_stage, levels = rev(c("NCI", "MCI", "AD"))),
      modality = factor(.data$modality, levels = c("RNA", "Protein"))
    ) |>
    dplyr::filter(!is.na(.data$display_stage)) |>
    dplyr::group_by(.data$modality, .data$display_stage) |>
    dplyr::summarise(n_samples = sum(.data$n_samples), .groups = "drop")

  }

  if (nrow(usable_stage) > 0) {
    x_max <- max(usable_stage$n_samples, na.rm = TRUE)
    sample_plot <- ggplot2::ggplot(
      usable_stage,
      ggplot2::aes(.data$n_samples, .data$display_stage, fill = .data$modality)
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.68),
        width = 0.56
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$n_samples),
        position = ggplot2::position_dodge(width = 0.68),
        hjust = -0.18,
        size = 3.2,
        color = "#111827"
      ) +
      ggplot2::scale_fill_manual(values = c(RNA = SUP_RNA_COL, Protein = SUP_PROTEIN_COL)) +
      ggplot2::scale_x_continuous(
        limits = c(0, x_max * 1.16),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        title = "Analyzable samples by clinical stage",
        x = "Samples",
        y = NULL,
        fill = NULL
      ) +
      theme_mainfigure_panel(base_size = 8.5) +
      ggplot2::theme(
        legend.position = "top",
        legend.justification = "left",
        panel.grid.major.y = ggplot2::element_blank()
      )
  } else {
    sample_plot <- unavailable_plot(
      "Samples by diagnosis/stage",
      "No stage column detected in loaded metadata"
    )
  }

  tile_df <- covariate_audit |>
    dplyr::mutate(
      label = dplyr::case_when(
        .data$column_available ~ scales::percent(.data$pct_non_missing, accuracy = 1),
        TRUE ~ "not available"
      ),
      modality = factor(.data$modality, levels = c("Protein", "RNA")),
      covariate = factor(.data$covariate, levels = rev(c("age", "sex", "PMI", "RIN", "Braak", "CERAD")))
    )

  covariate_plot <- ggplot2::ggplot(
    tile_df,
    ggplot2::aes(.data$modality, .data$covariate, fill = .data$pct_non_missing)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::scale_fill_gradient(
      low = "#F3F4F6",
      high = "#6B7280",
      limits = c(0, 1),
      labels = scales::percent,
      na.value = "#E5E7EB"
    ) +
    ggplot2::labs(
      title = "Covariate availability",
      x = NULL,
      y = NULL
    ) +
    theme_mainfigure_panel(base_size = 8.5) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  list(
    sample_plot = sample_plot,
    covariate_plot = covariate_plot,
    stage_audit = stage_audit,
    covariate_audit = covariate_audit
  )
}

# ============================================================
# Protein abundance audit
# ============================================================

make_abundance_panel <- function(interactor_tbl, mito_tbl, protein_info) {
  audit_path <- file.path(audits_dir, "SuppFig1_abundance_distribution_audit.csv")

  if (is.null(protein_info$matrix)) {
    write_skip_audit(audit_path, "SuppFig1D", "Protein matrix was not loaded or could not be oriented.")
    return(NULL)
  }

  if (is.null(mito_tbl) || nrow(mito_tbl) == 0) {
    write_skip_audit(audit_path, "SuppFig1D", "Mitochondrial background input was not loaded.")
    return(NULL)
  }

  protein_mat <- protein_info$matrix
  interactor_genes <- interactor_tbl$gene
  mito_genes <- mito_tbl$gene
  non_inventory_mito <- setdiff(mito_genes, interactor_genes)

  abundance_tbl <- tibble::tibble(
    gene = rownames(protein_mat),
    mean_protein_abundance = rowMeans(protein_mat, na.rm = TRUE),
    n_non_missing = rowSums(!is.na(protein_mat)),
    n_samples = ncol(protein_mat)
  ) |>
    dplyr::mutate(
      mean_protein_abundance = dplyr::if_else(
        is.nan(.data$mean_protein_abundance),
        NA_real_,
        .data$mean_protein_abundance
      )
    )

  group_tbl <- tibble::tibble(
    gene = c(interactor_genes, non_inventory_mito),
    group = c(
        rep("Hsp60/10 client", length(interactor_genes)),
      rep("Mitochondrial protein not in inventory", length(non_inventory_mito))
    )
  ) |>
    dplyr::distinct(.data$gene, .keep_all = TRUE)

  plot_tbl <- group_tbl |>
    dplyr::inner_join(abundance_tbl, by = "gene") |>
    dplyr::filter(!is.na(.data$mean_protein_abundance))

  readr::write_csv(plot_tbl, audit_path)

  if (nrow(plot_tbl) == 0 || dplyr::n_distinct(plot_tbl$group) < 1) {
    write_skip_audit(audit_path, "SuppFig1D", "No abundance values available for requested gene groups.")
    return(NULL)
  }

  message("Abundance distribution audit rows: ", nrow(plot_tbl))
  print(plot_tbl |> dplyr::count(.data$group, name = "n_genes"))

  invisible(plot_tbl)
}

# ============================================================
# Panel E: proteomics coverage
# ============================================================

make_detection_rate_panel <- function(interactor_tbl, mito_tbl, protein_info) {
  audit_path <- file.path(audits_dir, "SuppFig1_detection_rate_audit.csv")
  gene_audit_path <- file.path(audits_dir, "SuppFig1_detection_rate_gene_audit.csv")

  if (is.null(protein_info$matrix)) {
    write_skip_audit(audit_path, "SuppFig1E", "Protein matrix was not loaded or could not be oriented.")
    return(NULL)
  }

  protein_mat <- protein_info$matrix
  interactor_genes <- interactor_tbl$gene

  if (!is.null(mito_tbl) && nrow(mito_tbl) > 0) {
    non_inventory_mito <- setdiff(mito_tbl$gene, interactor_genes)
    group_tbl <- tibble::tibble(
      gene = c(interactor_genes, non_inventory_mito),
      group = c(
        rep("Hsp60/10 clients", length(interactor_genes)),
        rep("Non-client mitochondrial proteins", length(non_inventory_mito))
      )
    ) |>
      dplyr::distinct(.data$gene, .keep_all = TRUE)
    background_subtitle <- "Mitochondrial background: non-inventory MitoCarta proteins"
  } else {
    group_tbl <- tibble::tibble(
      gene = interactor_genes,
      group = "Hsp60/10 clients"
    )
    background_subtitle <- "Mitochondrial background unavailable; Hsp60/10 clients shown"
  }

  present_genes <- intersect(group_tbl$gene, rownames(protein_mat))
  rate_tbl_present <- tibble::tibble(
    gene = present_genes,
    n_non_missing = rowSums(!is.na(protein_mat[present_genes, , drop = FALSE])),
    n_samples = ncol(protein_mat)
  ) |>
    dplyr::mutate(detection_rate = .data$n_non_missing / .data$n_samples)

  gene_rate_tbl <- group_tbl |>
    dplyr::left_join(rate_tbl_present, by = "gene") |>
    dplyr::mutate(
      in_protein_matrix = .data$gene %in% present_genes,
      n_non_missing = dplyr::coalesce(.data$n_non_missing, 0),
      n_samples = dplyr::coalesce(.data$n_samples, ncol(protein_mat)),
      detection_rate = dplyr::coalesce(.data$detection_rate, 0),
      detection_rate_category = dplyr::case_when(
        .data$detection_rate == 0 ~ "0%",
        .data$detection_rate > 0 & .data$detection_rate <= 0.50 ~ ">0-50%",
        .data$detection_rate > 0.50 & .data$detection_rate < 0.95 ~ "50-95%",
        .data$detection_rate >= 0.95 ~ ">=95%",
        TRUE ~ NA_character_
      )
    )

  category_levels <- c("0%", ">0-50%", "50-95%", ">=95%")
  group_levels <- c("Non-client mitochondrial proteins", "Hsp60/10 clients")

  coverage_counts <- gene_rate_tbl |>
    dplyr::mutate(
      group = factor(.data$group, levels = group_levels),
      detection_rate_category = factor(.data$detection_rate_category, levels = category_levels)
    ) |>
    dplyr::count(.data$group, .data$detection_rate_category, name = "n_genes") |>
    tidyr::complete(
      group,
      detection_rate_category,
      fill = list(n_genes = 0)
    ) |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(
      group_n = sum(.data$n_genes),
      fraction = dplyr::if_else(.data$group_n > 0, .data$n_genes / .data$group_n, NA_real_),
      group_label = paste0(as.character(.data$group), "\n(n = ", .data$group_n, ")")
    ) |>
    dplyr::ungroup()

  readr::write_csv(coverage_counts, audit_path)
  readr::write_csv(gene_rate_tbl, gene_audit_path)

  message("Detection-rate category audit rows: ", nrow(coverage_counts))
  print(coverage_counts |> dplyr::select("group", "detection_rate_category", "n_genes", "group_n"))

  group_label_levels <- coverage_counts |>
    dplyr::distinct(.data$group, .data$group_label) |>
    dplyr::arrange(.data$group) |>
    dplyr::pull(.data$group_label)

  p <- ggplot2::ggplot(
    coverage_counts,
    ggplot2::aes(
      x = .data$fraction,
      y = factor(.data$group_label, levels = group_label_levels),
      fill = .data$detection_rate_category
    )
  ) +
    ggplot2::geom_col(width = 0.52, color = "white", linewidth = 0.35) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = c(
        "0%" = "#ECECEC",
        ">0-50%" = "#D2D6DB",
        "50-95%" = "#A7AFB7",
        ">=95%" = "#6B7280"
      ),
      drop = FALSE,
      labels = c("0", "0-50", "50-95", "95+"),
      name = "Rate",
      guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = "Protein detection-rate categories",
      subtitle = NULL,
      x = "Share of genes",
      y = NULL
    ) +
    theme_mainfigure_panel(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(0.75, "lines"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  attr(p, "coverage_counts") <- coverage_counts
  attr(p, "gene_rate_audit") <- gene_rate_tbl
  p
}

# ============================================================
# Run summary audit
# ============================================================

first_non_missing <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) NA_character_ else as.character(x[1])
}

format_covariate_columns <- function(covariate_audit, modality_label) {
  covariate_audit |>
    dplyr::filter(.data$modality == .env$modality_label) |>
    dplyr::mutate(
      detected_column = dplyr::if_else(
        is.na(.data$detected_column),
        "not available",
        .data$detected_column
      ),
      summary = paste0(.data$covariate, "=", .data$detected_column)
    ) |>
    dplyr::pull(.data$summary) |>
    paste(collapse = ", ")
}

input_source_description <- function(input_name_value, fallback = NA_character_) {
  if (
    exists("input_source_audit") &&
      is.data.frame(input_source_audit) &&
      input_name_value %in% input_source_audit$input_name
  ) {
    source_path <- input_source_audit |>
      dplyr::filter(.data$input_name == .env$input_name_value, .data$status == "loaded") |>
      dplyr::pull(.data$path) |>
      utils::tail(1)

    if (length(source_path) == 1 && !is.na(source_path) && nzchar(source_path)) {
      return(source_path)
    }
  }

  fallback
}

write_supfig1_run_summary <- function(
  interactor_tbl,
  detection_counts,
  rna_info,
  protein_info,
  stage_audit,
  covariate_audit,
  mito_tbl,
  mito_input_name,
  panels_generated,
  panels_skipped
) {
  count_for <- function(display_group) {
    out <- detection_counts |>
      dplyr::filter(.data$display_group == .env$display_group) |>
      dplyr::summarise(n = sum(.data$n_interactors), .groups = "drop") |>
      dplyr::pull(.data$n)
    if (length(out) == 0 || is.na(out)) 0 else out
  }

  rna_stage_column <- stage_audit |>
    dplyr::filter(.data$modality == "RNA") |>
    dplyr::pull(.data$stage_column) |>
    first_non_missing()

  protein_stage_column <- stage_audit |>
    dplyr::filter(.data$modality == "Protein") |>
    dplyr::pull(.data$stage_column) |>
    first_non_missing()

  mito_source <- if (!is.null(mito_tbl) && nrow(mito_tbl) > 0) {
    paste0(
      input_source_description(mito_input_name, fallback = mito_input_name),
      " (",
      nrow(mito_tbl),
      " genes after gene-symbol harmonization)"
    )
  } else {
    "not available"
  }

  summary_lines <- c(
    "Supplementary Figure 1 run summary",
    paste0("generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "source inventory: Bie et al. Hsp60/10 client inventory",
    "original reported inventory size: 323",
    paste0("analyzed size after gene-symbol harmonization: ", nrow(interactor_tbl)),
    paste0("RNA-only count: ", count_for("RNA only")),
    paste0("protein-only count: ", count_for("Protein only")),
    paste0("RNA + protein count: ", count_for("RNA + protein")),
    paste0("not-detected count: ", count_for("Not detected")),
    paste0("RNA sample count: ", ifelse(is.null(rna_info$matrix), 0, ncol(rna_info$matrix))),
    paste0("protein sample count: ", ifelse(is.null(protein_info$matrix), 0, ncol(protein_info$matrix))),
    paste0("stage column used for RNA: ", rna_stage_column %||% "not available"),
    paste0("stage column used for protein: ", protein_stage_column %||% "not available"),
    paste0("covariate columns detected for RNA: ", format_covariate_columns(covariate_audit, "RNA")),
    paste0("covariate columns detected for protein: ", format_covariate_columns(covariate_audit, "Protein")),
    paste0("mitochondrial background source: ", mito_source),
    paste0("panels generated: ", paste(panels_generated, collapse = "; ")),
    paste0(
      "panels skipped: ",
      ifelse(length(panels_skipped) == 0, "none", paste(panels_skipped, collapse = "; "))
    )
  )

  run_summary_path <- file.path(audits_dir, "SuppFig1_run_summary.txt")
  writeLines(summary_lines, run_summary_path)
  message("Saved run summary: ", run_summary_path)
  invisible(run_summary_path)
}

# ============================================================
# Build and save Supplementary Figure 1
# ============================================================

make_supfig1_cohort_detection <- function(inputs) {
  message("Starting Supplementary Figure 1.")

  interactor_input <- get_input(
    c(
      "hsp60_hsp10_interactor_inventory",
      "all_hsp60_10_client_tbl",
      "hsp60_client_tbl",
      "hsp_client_tbl",
      "hsp_clients_all",
      "hsp60_clients"
    ),
    required = TRUE,
    label = "Hsp60/10 interactor inventory table"
  )
  interactor_tbl <- as_gene_table(interactor_input$object, "Hsp60/10 interactor inventory table")

  if (nrow(interactor_tbl) == 0) {
    stop("Hsp60/10 interactor inventory table loaded but contained zero usable gene symbols.", call. = FALSE)
  }

  mito_input <- get_input(
    c("mito_background_tbl", "mitochondrial_background_tbl"),
    required = FALSE,
    label = "mitochondrial background table"
  )
  mito_tbl <- if (!is.null(mito_input$object)) {
    as_gene_table(mito_input$object, "mitochondrial background table")
  } else {
    NULL
  }

  message("Source inventory: Bie et al. 2020 Hsp60/10 interactor inventory")
  message("Original reported inventory: 323 interactors")
  message("Analyzed after gene-symbol harmonization: ", nrow(interactor_tbl), " interactors")
  message("Interactor inventory dimensions after gene-symbol harmonization: ", paste(dim(interactor_tbl), collapse = " x "))
  if (!is.null(mito_tbl)) {
    message("Mitochondrial background dimensions after gene cleaning: ", paste(dim(mito_tbl), collapse = " x "))
  }

  rna_info <- prepare_matrix_input(
    matrix_candidates = c("rna_mat", "vst_mat", "rna_vst", "rna_matrix", "vst_symbol_mat", "rna_mat_stage"),
    meta_candidates = c("rna_meta", "meta_rna_sub", "rna_metadata", "rna_meta_in", "meta_rna"),
    interactor_genes = interactor_tbl$gene,
    label = "RNA matrix"
  )

  protein_info <- prepare_matrix_input(
    matrix_candidates = c("prot_mat", "protein_matrix", "protein_mat", "tmt_mat", "prot_expr_mat", "prot_mat_in", "prot_mat_raw"),
    meta_candidates = c("protein_meta", "prot_meta_aligned", "prot_meta_path", "prot_meta", "meta_protein_sub", "prot_meta_in"),
    interactor_genes = interactor_tbl$gene,
    label = "Protein matrix"
  )

  if (!is.null(rna_info$matrix)) {
    message("RNA matrix genes/samples: ", nrow(rna_info$matrix), " genes x ", ncol(rna_info$matrix), " samples")
  }
  if (!is.null(protein_info$matrix)) {
    message("Protein matrix genes/samples: ", nrow(protein_info$matrix), " genes x ", ncol(protein_info$matrix), " samples")
  }

  panel_a <- make_workflow_panel(inputs, rna_info, protein_info, n_interactors = nrow(interactor_tbl))
  panel_b <- make_detection_panel(interactor_tbl, rna_info, protein_info)
  metadata_panels <- make_metadata_panels(rna_info$metadata, protein_info$metadata)
  panel_c <- metadata_panels$sample_plot
  panel_d <- metadata_panels$covariate_plot
  abundance_audit <- make_abundance_panel(interactor_tbl, mito_tbl, protein_info)
  panel_e <- make_detection_rate_panel(interactor_tbl, mito_tbl, protein_info)

  panel_a_labeled <- panel_label(panel_a, "A")
  panel_b_labeled <- panel_label(panel_b, "B")
  panel_c_labeled <- panel_label(panel_c, "C")
  panel_d_labeled <- panel_label(panel_d, "D")

  panels_generated <- c(
    "A Analysis workflow",
    "B Client detection across RNA and proteomics",
    "C Samples by diagnosis/stage",
    "D Covariate availability"
  )
  panels_skipped <- character()

  if (is.null(panel_e)) {
    panels_skipped <- c(panels_skipped, "E Protein detection-rate categories")
    panel_e <- unavailable_plot(
      "Protein detection-rate categories",
      "Protein matrix was not loaded or could not be oriented"
    )
  } else {
    panels_generated <- c(panels_generated, "E Protein detection-rate categories")
  }

  panel_e_labeled <- panel_label(panel_e, "E")

  panel_paths <- c(
    save_panel_set(panel_a_labeled, "SuppFig1A_workflow", width = 12.6, height = 2.35, output_dir = panels_dir),
    save_panel_set(panel_b_labeled, "SuppFig1B_detection_counts", width = 5.9, height = 3.8, output_dir = panels_dir),
    save_panel_set(panel_c_labeled, "SuppFig1C_sample_counts", width = 5.9, height = 3.8, output_dir = panels_dir),
    save_panel_set(panel_d_labeled, "SuppFig1D_covariate_availability", width = 5.9, height = 3.8, output_dir = panels_dir),
    save_panel_set(panel_e_labeled, "SuppFig1E_proteomics_coverage", width = 5.9, height = 3.8, output_dir = panels_dir)
  )

  composite_body <- panel_a_labeled /
    (panel_b_labeled | panel_c_labeled) /
    (panel_d_labeled | panel_e_labeled) +
    patchwork::plot_layout(heights = c(0.56, 1, 1))

  composite <- composite_body

  save_plot_set(
    composite,
    "Supplementary_Figure_1_cohort_detection",
    width = 12.6,
    height = 10.2,
    output_dir = figures_dir
  )

  run_summary_path <- write_supfig1_run_summary(
    interactor_tbl = interactor_tbl,
    detection_counts = attr(panel_b, "detection_counts"),
    rna_info = rna_info,
    protein_info = protein_info,
    stage_audit = metadata_panels$stage_audit,
    covariate_audit = metadata_panels$covariate_audit,
    mito_tbl = mito_tbl,
    mito_input_name = mito_input$name,
    panels_generated = panels_generated,
    panels_skipped = panels_skipped
  )

  output_summary <- tibble::tibble(
    output_type = c("composite_pdf", "composite_png", "panels_dir", "audits_dir", "tables_dir", "run_summary"),
    path = c(
      file.path(figures_dir, "Supplementary_Figure_1_cohort_detection.pdf"),
      file.path(figures_dir, "Supplementary_Figure_1_cohort_detection.png"),
      panels_dir,
      audits_dir,
      tables_dir,
      run_summary_path
    )
  )

  message("Supplementary Figure 1 complete.")
  print(output_summary)

  invisible(list(
    interactor_tbl = interactor_tbl,
    mito_tbl = mito_tbl,
    rna_info = rna_info,
    protein_info = protein_info,
    panels = list(panel_a_labeled, panel_b_labeled, panel_c_labeled, panel_d_labeled, panel_e_labeled),
    panel_paths = panel_paths,
    abundance_audit = abundance_audit,
    output_summary = output_summary
  ))
}

supfig1_outputs <- make_supfig1_cohort_detection(inputs)
