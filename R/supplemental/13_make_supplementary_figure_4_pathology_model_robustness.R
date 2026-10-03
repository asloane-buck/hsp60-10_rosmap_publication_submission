# ============================================================
# Supplementary Figure 4: pathology-model robustness
# ============================================================

# Purpose:
# Test whether Hsp60/10 pathology coupling persists beyond simple
# diagnostic/stage grouping and is robust to covariate/CERAD adjustment.

if (!exists("inputs")) {
  stop("Run 01_supplemental_load_inputs.R before 13_make_supplementary_figure_4_pathology_model_robustness.R.", call. = FALSE)
}
if (!exists("pick_col")) {
  stop("Run 02_supplemental_helper_functions.R before 13_make_supplementary_figure_4_pathology_model_robustness.R.", call. = FALSE)
}

# ============================================================
# Main Figure style settings
# ============================================================

sf4_hsp_col <- "#B23A48"
sf4_background_col <- "#4C78A8"
sf4_hsp_fill <- "#F4D9DE"
sf4_background_fill <- "#DCE7F4"
sf4_text <- "#111827"

sf4_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey25", linewidth = 0.45),
      axis.line = ggplot2::element_line(color = "grey25", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "grey25", linewidth = 0.35),
      axis.text = ggplot2::element_text(color = "grey15"),
      axis.title = ggplot2::element_text(color = "grey10", face = "bold"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

sf4_panel_label <- function(plot, label, tag_position = c(-0.08, 1.07)) {
  plot +
    ggplot2::labs(tag = label) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 9.5),
      plot.tag.position = tag_position,
      plot.margin = ggplot2::margin(13, 8, 6, 18)
    )
}

sf4_group_colors <- c(
  "Non-client mitochondrial proteins" = sf4_background_col,
  "Hsp60/10 clients" = sf4_hsp_col
)

sf4_group_fills <- c(
  "Non-client mitochondrial proteins" = sf4_background_fill,
  "Hsp60/10 clients" = sf4_hsp_fill
)

# ============================================================
# Input helpers
# ============================================================

sf4_get_input <- function(candidate_names, required = TRUE, label = "input") {
  selected_name <- candidate_names[candidate_names %in% names(inputs)][1]

  if (!is.na(selected_name)) {
    message("Detected ", label, ": ", selected_name)
    return(list(name = selected_name, object = inputs[[selected_name]]))
  }

  selected_global <- candidate_names[vapply(candidate_names, exists, logical(1), envir = .GlobalEnv)][1]
  if (!is.na(selected_global)) {
    message("Detected ", label, " in global environment: ", selected_global)
    return(list(name = selected_global, object = get(selected_global, envir = .GlobalEnv)))
  }

  if (isTRUE(required)) {
    stop(
      "Missing required input for ", label, ". Tried: ",
      paste(candidate_names, collapse = ", "),
      ". Available input names: ",
      paste(names(inputs), collapse = ", "),
      call. = FALSE
    )
  }

  list(name = NA_character_, object = NULL)
}

sf4_validate_null_tbl <- function(tbl, label) {
  required_cols <- c(
    "gene",
    "protein_late_decline_magnitude",
    "inverse_braak_magnitude"
  )
  
  missing_cols <- setdiff(required_cols, colnames(tbl))
  if (length(missing_cols) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(tbl), collapse = ", "),
      call. = FALSE
    )
  }
  
  out <- tibble::as_tibble(tbl) |>
    dplyr::mutate(gene = clean_gene(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene))
  
  if (!"braak_rho" %in% colnames(out)) {
    out <- out |>
      dplyr::mutate(braak_rho = NA_real_)
  }
  
  out
}

sf4_as_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

sf4_z <- function(x) {
  x <- sf4_as_numeric(x)
  if (sum(is.finite(x)) < 3 || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

sf4_pick_col <- function(df, candidates, required = TRUE, label = "column") {
  pick_col(df, candidates, required = required, label = label)
}

sf4_metadata_sample_ids <- function(meta, label) {
  if (is.null(meta)) return(character())

  df <- as.data.frame(meta)
  sample_col <- sf4_pick_col(
    df,
    c(
      "sample_id", "SampleID", "sampleid", "sample", "sample_name",
      "SpecimenID", "specimenID", "specimen_id", "batch.channel", ".sample_id",
      "rna_id", "protein_id", "id"
    ),
    required = FALSE,
    label = paste(label, "sample identifier")
  )

  if (!is.na(sample_col)) {
    return(as.character(df[[sample_col]]))
  }

  if (!is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
    return(rownames(df))
  }

  character()
}

sf4_prepare_protein_matrix <- function(hsp_genes) {
  matrix_input <- sf4_get_input(
    c("prot_mat_raw"),
    required = TRUE,
    label = "protein matrix"
  )

  meta_input <- sf4_get_input(
    c("protein_meta", "prot_meta_aligned", "prot_meta", "meta_protein_sub", "prot_meta_in", "protein_metadata"),
    required = TRUE,
    label = "protein metadata"
  )

  sample_ids <- sf4_metadata_sample_ids(meta_input$object, "Protein matrix")

  oriented <- infer_and_orient_matrix(
    object = matrix_input$object,
    gene_symbols = hsp_genes,
    sample_ids = sample_ids,
    label = "Protein matrix",
    required = TRUE
  )

  if (is.null(oriented$matrix)) {
    stop("Protein matrix could not be oriented.", call. = FALSE)
  }

  list(
    matrix_name = matrix_input$name,
    meta_name = meta_input$name,
    matrix = oriented$matrix,
    metadata = meta_input$object,
    diagnostics = oriented$diagnostics
  )
}

sf4_prepare_model_metadata <- function(meta, matrix_sample_ids) {
  df <- as.data.frame(meta)

  sample_col <- sf4_pick_col(
    df,
    c(
      "sample_id", "SampleID", "sampleid", "sample", "sample_name",
      "SpecimenID", "specimenID", "specimen_id", "batch.channel", ".sample_id",
      "rna_id", "protein_id", "id"
    ),
    required = FALSE,
    label = "protein metadata sample ID"
  )

  if (is.na(sample_col) && !is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
    df$.sample_id_from_rownames <- rownames(df)
    sample_col <- ".sample_id_from_rownames"
  }

  if (is.na(sample_col)) {
    stop(
      "Could not find a protein metadata sample ID column. Available columns: ",
      paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  braak_col <- sf4_pick_col(df, c("braak_num", "Braak", "braak", "braaksc", "braak_stage", "Braak.Stage"), TRUE, "Braak")
  cerad_col <- sf4_pick_col(df, c("cerad_num", "CERAD", "cerad", "ceradsc", "cerad_score", "CERAD.Score"), TRUE, "CERAD")
  age_col <- sf4_pick_col(df, c("age_num", "age_death", "age", "age_at_death", "age_c", "Age"), TRUE, "age")
  sex_col <- sf4_pick_col(df, c("sex_factor", "sex_label", "sex", "msex", "Sex", "gender"), TRUE, "sex")
  pmi_col <- sf4_pick_col(df, c("pmi_num", "pmi", "PMI", "pmi_c", "postmortem_interval", "postmortem.interval"), TRUE, "PMI")
  batch_col <- sf4_pick_col(df, c("batch_factor", "tmt_batch", "batch", "Batch"), TRUE, "TMT batch")

  meta_out <- tibble::as_tibble(df) |>
    dplyr::mutate(sample_id = as.character(.data[[sample_col]])) |>
    dplyr::filter(.data$sample_id %in% matrix_sample_ids) |>
    dplyr::distinct(.data$sample_id, .keep_all = TRUE) |>
    dplyr::mutate(sample_id = factor(.data$sample_id, levels = matrix_sample_ids)) |>
    dplyr::arrange(.data$sample_id) |>
    dplyr::mutate(
      sample_id = as.character(.data$sample_id),
      braak_num_model = sf4_as_numeric(.data[[braak_col]]),
      cerad_num_model = sf4_as_numeric(.data[[cerad_col]]),
      age_model = sf4_as_numeric(.data[[age_col]]),
      pmi_model = sf4_as_numeric(.data[[pmi_col]]),
      sex_model = as.factor(.data[[sex_col]]),
      batch_model = as.factor(.data[[batch_col]])
    )

  readr::write_csv(
    tibble::tibble(
      conceptual_variable = c("sample_id", "braak", "cerad", "age", "sex", "pmi", "tmt_batch"),
      detected_column = c(sample_col, braak_col, cerad_col, age_col, sex_col, pmi_col, batch_col),
      n_matrix_samples = length(matrix_sample_ids),
      n_metadata_rows_matched = nrow(meta_out)
    ),
    file.path(audits_dir, "SuppFig4_metadata_column_audit.csv")
  )

  message("Supplementary Figure 4 metadata alignment:")
  message("  matrix samples: ", length(matrix_sample_ids))
  message("  metadata rows matched: ", nrow(meta_out))
  message("  Braak column: ", braak_col)
  message("  CERAD column: ", cerad_col)
  message("  age column: ", age_col)
  message("  sex column: ", sex_col)
  message("  PMI column: ", pmi_col)
  message("  TMT batch column: ", batch_col)

  meta_out
}

# ============================================================
# Model fitting
# ============================================================

sf4_extract_coef <- function(fit, coef_name) {
  sm <- summary(fit)$coefficients
  if (!coef_name %in% rownames(sm)) {
    return(tibble::tibble(beta = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  tibble::tibble(
    beta = unname(sm[coef_name, "Estimate"]),
    se = unname(sm[coef_name, "Std. Error"]),
    p_value = unname(sm[coef_name, "Pr(>|t|)"])
  )
}

sf4_fit_gene <- function(gene, protein_mat, meta_model) {
  y <- as.numeric(protein_mat[gene, meta_model$sample_id])

  model_df <- meta_model |>
    dplyr::mutate(
      y_raw = y,
      y_z = sf4_z(.data$y_raw),
      braak_z = sf4_z(.data$braak_num_model),
      cerad_z = sf4_z(.data$cerad_num_model),
      age_z = sf4_z(.data$age_model),
      pmi_z = sf4_z(.data$pmi_model)
    )

  # Model 1: Braak + core covariates.
  m1_df <- model_df |>
    dplyr::select("y_z", "braak_z", "age_z", "pmi_z", "sex_model", "batch_model") |>
    tidyr::drop_na("y_z", "braak_z")

  m1_terms <- c("braak_z")
  if (sum(is.finite(m1_df$age_z)) >= 10 && stats::sd(m1_df$age_z, na.rm = TRUE) > 0) m1_terms <- c(m1_terms, "age_z")
  if (sum(is.finite(m1_df$pmi_z)) >= 10 && stats::sd(m1_df$pmi_z, na.rm = TRUE) > 0) m1_terms <- c(m1_terms, "pmi_z")
  if (dplyr::n_distinct(stats::na.omit(m1_df$sex_model)) > 1) m1_terms <- c(m1_terms, "sex_model")
  if (dplyr::n_distinct(stats::na.omit(m1_df$batch_model)) > 1) m1_terms <- c(m1_terms, "batch_model")

  m1_result <- tryCatch({
    if (nrow(m1_df) < 30 || stats::sd(m1_df$y_z, na.rm = TRUE) == 0) {
      tibble::tibble(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = nrow(m1_df))
    } else {
      fit <- stats::lm(stats::as.formula(paste("y_z ~", paste(m1_terms, collapse = " + "))), data = m1_df)
      sf4_extract_coef(fit, "braak_z") |>
        dplyr::mutate(n = stats::nobs(fit))
    }
  }, error = function(e) {
    tibble::tibble(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = nrow(m1_df))
  })

  # Model 2: Braak + CERAD + core covariates.
  m2_df <- model_df |>
    dplyr::select("y_z", "braak_z", "cerad_z", "age_z", "pmi_z", "sex_model", "batch_model") |>
    tidyr::drop_na("y_z", "braak_z", "cerad_z")

  m2_terms <- c("braak_z", "cerad_z")
  if (sum(is.finite(m2_df$age_z)) >= 10 && stats::sd(m2_df$age_z, na.rm = TRUE) > 0) m2_terms <- c(m2_terms, "age_z")
  if (sum(is.finite(m2_df$pmi_z)) >= 10 && stats::sd(m2_df$pmi_z, na.rm = TRUE) > 0) m2_terms <- c(m2_terms, "pmi_z")
  if (dplyr::n_distinct(stats::na.omit(m2_df$sex_model)) > 1) m2_terms <- c(m2_terms, "sex_model")
  if (dplyr::n_distinct(stats::na.omit(m2_df$batch_model)) > 1) m2_terms <- c(m2_terms, "batch_model")

  m2_result <- tryCatch({
    if (nrow(m2_df) < 30 || stats::sd(m2_df$y_z, na.rm = TRUE) == 0) {
      tibble::tibble(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = nrow(m2_df))
    } else {
      fit <- stats::lm(stats::as.formula(paste("y_z ~", paste(m2_terms, collapse = " + "))), data = m2_df)
      sf4_extract_coef(fit, "braak_z") |>
        dplyr::mutate(n = stats::nobs(fit))
    }
  }, error = function(e) {
    tibble::tibble(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = nrow(m2_df))
  })

  tibble::tibble(gene = gene) |>
    dplyr::bind_cols(
      m1_result |>
        dplyr::rename(
          braak_beta_adjusted = "beta",
          braak_se_adjusted = "se",
          braak_p_adjusted = "p_value",
          n_adjusted = "n"
        ),
      m2_result |>
        dplyr::rename(
          braak_beta_cerad_adjusted = "beta",
          braak_se_cerad_adjusted = "se",
          braak_p_cerad_adjusted = "p_value",
          n_cerad_adjusted = "n"
        )
    ) |>
    dplyr::mutate(
      inverse_braak_beta_adjusted = dplyr::case_when(
        !is.finite(.data$braak_beta_adjusted) ~ NA_real_,
        .data$braak_beta_adjusted < 0 ~ abs(.data$braak_beta_adjusted),
        TRUE ~ 0
      ),
      inverse_braak_beta_cerad_adjusted = dplyr::case_when(
        !is.finite(.data$braak_beta_cerad_adjusted) ~ NA_real_,
        .data$braak_beta_cerad_adjusted < 0 ~ abs(.data$braak_beta_cerad_adjusted),
        TRUE ~ 0
      )
    )
}

sf4_fit_all_models <- function(protein_mat, meta_model, genes) {
  genes_use <- intersect(clean_gene(genes), rownames(protein_mat))
  message("Fitting pathology robustness models for ", length(genes_use), " proteins.")

  purrr::map_dfr(genes_use, ~ sf4_fit_gene(.x, protein_mat, meta_model))
}

# ============================================================
# Plot helpers
# ============================================================

sf4_p_label <- function(p) {
  if (is.na(p)) return("Wilcoxon P = NA")
  if (p < 0.001) return(paste0("Wilcoxon P = ", formatC(p, format = "e", digits = 2)))
  paste0("Wilcoxon P = ", signif(p, 3))
}

sf4_wilcox_group_p <- function(df, value_col) {
  df2 <- df |>
    dplyr::filter(.data$group %in% names(sf4_group_colors)) |>
    dplyr::select("group", value = dplyr::all_of(value_col)) |>
    tidyr::drop_na()

  if (dplyr::n_distinct(df2$group) < 2) return(NA_real_)

  tryCatch(
    stats::wilcox.test(value ~ group, data = df2, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

sf4_make_distribution_plot <- function(plot_tbl, value_col, title, subtitle, ylab) {
  p_val <- sf4_wilcox_group_p(plot_tbl, value_col)
  y_max <- max(plot_tbl[[value_col]], na.rm = TRUE)
  y_min <- min(plot_tbl[[value_col]], na.rm = TRUE)
  y_rng <- y_max - y_min
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- 0.1

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$group, y = .data[[value_col]], fill = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_violin(width = 0.78, alpha = 0.75, color = NA, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95, color = "grey25", linewidth = 0.4) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data$group),
      width = 0.08,
      height = 0,
      size = 0.75,
      alpha = 0.28,
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = y_max + 0.08 * y_rng,
      label = sf4_p_label(p_val),
      size = 3.2,
      fontface = "bold",
      color = "grey20"
    ) +
    ggplot2::scale_fill_manual(values = sf4_group_fills, guide = "none") +
    ggplot2::scale_color_manual(values = sf4_group_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = c(
      "Non-client mitochondrial proteins" = "Non-client\nmitochondrial proteins",
      "Hsp60/10 clients" = "Hsp60/10\nclients"
    )) +
    ggplot2::coord_cartesian(ylim = c(y_min - 0.03 * y_rng, y_max + 0.16 * y_rng), clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = ylab
    ) +
    sf4_theme(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.25)
    )
}

sf4_make_scatter_plot <- function(plot_tbl) {
  scatter_tbl <- plot_tbl |>
    dplyr::filter(
      is.finite(.data$protein_late_decline_magnitude),
      is.finite(.data$inverse_braak_beta_cerad_adjusted)
    )

  rho <- tryCatch(
    stats::cor(
      scatter_tbl$protein_late_decline_magnitude,
      scatter_tbl$inverse_braak_beta_cerad_adjusted,
      method = "spearman",
      use = "complete.obs"
    ),
    error = function(e) NA_real_
  )

  ggplot2::ggplot(
    scatter_tbl,
    ggplot2::aes(
      x = .data$protein_late_decline_magnitude,
      y = .data$inverse_braak_beta_cerad_adjusted
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_point(
      data = scatter_tbl |> dplyr::filter(.data$group == "Non-client mitochondrial proteins"),
      color = sf4_background_col,
      alpha = 0.35,
      size = 1.4
    ) +
    ggplot2::geom_point(
      data = scatter_tbl |> dplyr::filter(.data$group == "Hsp60/10 clients"),
      color = sf4_hsp_col,
      alpha = 0.78,
      size = 1.8
    ) +
    ggplot2::annotate(
      "label",
      x = Inf,
      y = Inf,
      hjust = 1.03,
      vjust = 1.1,
      label = paste0("n = ", nrow(scatter_tbl), "; Spearman rho = ", signif(rho, 3),
        "\nP = ", formatC(stats::cor.test(scatter_tbl$protein_late_decline_magnitude,
          scatter_tbl$inverse_braak_beta_cerad_adjusted, method = "spearman", exact = FALSE)$p.value,
          format = "e", digits = 2)),
      size = 3.2,
      fill = "white",
      color = "grey20"
    ) +
    ggplot2::labs(
      title = "Late-stage decline vs adjusted Braak",
      subtitle = NULL,
      x = "Late-stage protein decline magnitude",
      y = "CERAD-adjusted inverse Braak"
    ) +
    sf4_theme(base_size = 8.5)
}

# ============================================================
# Main builder
# ============================================================

make_supfig4_pathology_model_robustness <- function(inputs) {
  message("Starting Supplementary Figure 4: pathology-model robustness.")

  hsp_input <- sf4_get_input(
    c("hsp_null_tbl", "hsp60_10_null_tbl", "hsp_specificity_tbl"),
    required = TRUE,
    label = "Hsp60/10 pathology table"
  )
  background_input <- sf4_get_input(
    c("background_null_pool", "mitochondrial_background_null_pool", "background_specificity_pool"),
    required = TRUE,
    label = "non-Hsp60/10 mitochondrial background table"
  )

  hsp_tbl <- sf4_validate_null_tbl(hsp_input$object, "hsp_null_tbl") |>
    dplyr::mutate(group = "Hsp60/10 clients")

  background_tbl <- sf4_validate_null_tbl(background_input$object, "background_null_pool") |>
    dplyr::filter(!.data$gene %in% hsp_tbl$gene) |>
    dplyr::mutate(group = "Non-client mitochondrial proteins")

  base_tbl <- dplyr::bind_rows(hsp_tbl, background_tbl) |>
    dplyr::mutate(
      group = factor(
        .data$group,
        levels = c("Non-client mitochondrial proteins", "Hsp60/10 clients")
      )
    )

  frozen <- isTRUE(get0("SUPP_USE_FROZEN_MODELS", ifnotfound = FALSE))
  if (frozen) {
    frozen_path <- file.path(audits_dir, "SuppFig4_plotted_gene_level_values.csv")
    plot_tbl <- readr::read_csv(frozen_path, show_col_types = FALSE)
    required <- c("gene", "group", "protein_late_decline_magnitude",
      "inverse_braak_beta_adjusted", "inverse_braak_beta_cerad_adjusted")
    if (!all(required %in% names(plot_tbl)) || nrow(plot_tbl) != 915 ||
        sum(plot_tbl$group == "Hsp60/10 clients") != 306 ||
        sum(plot_tbl$group == "Non-client mitochondrial proteins") != 609 ||
        anyDuplicated(plot_tbl$gene) ||
        any(!is.finite(plot_tbl$inverse_braak_beta_adjusted)) ||
        any(!is.finite(plot_tbl$inverse_braak_beta_cerad_adjusted))) {
      stop("S4 frozen source does not contain the verified 306-client/609-background models.", call. = FALSE)
    }
    plot_tbl$group <- factor(plot_tbl$group,
      levels = c("Non-client mitochondrial proteins", "Hsp60/10 clients"))
    adjusted_tbl <- plot_tbl
    message("S4: rendering frozen non-residualized abundance regression coefficients.")
  } else {
  protein_info <- sf4_prepare_protein_matrix(hsp_tbl$gene)
  protein_mat <- protein_info$matrix
  meta_model <- sf4_prepare_model_metadata(protein_info$metadata, colnames(protein_mat))

  if (nrow(meta_model) < 30) {
    stop("Too few protein metadata rows matched to protein matrix samples for modeling.", call. = FALSE)
  }

  protein_mat <- protein_mat[, meta_model$sample_id, drop = FALSE]

  genes_for_models <- intersect(base_tbl$gene, rownames(protein_mat))

  if (length(genes_for_models) == 0) {
    readr::write_csv(
      protein_info$diagnostics %||% tibble::tibble(),
      file.path(audits_dir, "SuppFig4_protein_matrix_orientation_failed_audit.csv")
    )
    stop(
      "Supp Fig 4 found zero overlap between base_tbl genes and protein matrix rownames after orientation. ",
      "Inspect the matrix orientation audit and the cleaned gene identifiers in prot_mat_raw. ",
      "Do not substitute a covariate-residualized matrix into a directly adjusted model.",
      call. = FALSE
    )
  }

  adjusted_tbl <- sf4_fit_all_models(protein_mat, meta_model, genes_for_models)

  # Do not let prior session patches or previously joined columns create .x/.y suffixes.
  # Supp Fig 4 should always use freshly fitted values from adjusted_tbl.
  base_tbl <- base_tbl |>
    dplyr::select(
      -dplyr::any_of(c(
        "inverse_braak_beta_adjusted",
        "braak_beta_adjusted",
        "braak_se_adjusted",
        "braak_p_adjusted",
        "n_adjusted",
        "braak_beta_cerad_adjusted",
        "braak_se_cerad_adjusted",
        "braak_p_cerad_adjusted",
        "n_cerad_adjusted",
        "inverse_braak_beta_cerad_adjusted",
        "inverse_braak_beta_adjusted.x",
        "inverse_braak_beta_adjusted.y",
        "braak_beta_adjusted.x",
        "braak_beta_adjusted.y",
        "braak_p_adjusted.x",
        "braak_p_adjusted.y",
        "n_adjusted.x",
        "n_adjusted.y",
        "braak_metric_source",
        "braak_metric_model"
      ))
    )

  required_adjusted_cols <- c(
    "gene",
    "braak_beta_adjusted",
    "braak_p_adjusted",
    "n_adjusted",
    "inverse_braak_beta_adjusted",
    "inverse_braak_beta_cerad_adjusted"
  )

  missing_adjusted_cols <- setdiff(required_adjusted_cols, colnames(adjusted_tbl))
  if (length(missing_adjusted_cols) > 0) {
    stop(
      "Adjusted model table is missing expected columns: ",
      paste(missing_adjusted_cols, collapse = ", "),
      "\nAvailable adjusted_tbl columns: ",
      paste(colnames(adjusted_tbl), collapse = ", "),
      call. = FALSE
    )
  }

  plot_tbl <- base_tbl |>
    dplyr::left_join(adjusted_tbl, by = "gene") |>
    dplyr::filter(.data$gene %in% genes_for_models)

  if (!"inverse_braak_beta_adjusted" %in% colnames(plot_tbl)) {
    stop(
      "plot_tbl lacks inverse_braak_beta_adjusted after joining adjusted_tbl. ",
      "This should not happen after stale adjusted columns were removed from base_tbl.",
      call. = FALSE
    )
  }

  readr::write_csv(base_tbl, file.path(audits_dir, "SuppFig4_input_pathology_table_audit.csv"))
  readr::write_csv(adjusted_tbl, file.path(audits_dir, "SuppFig4_adjusted_braak_model_results.csv"))
  readr::write_csv(plot_tbl, file.path(audits_dir, "SuppFig4_plotted_gene_level_values.csv"))

  }

  summary_tbl <- plot_tbl |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_inverse_braak = stats::median(.data$inverse_braak_beta_adjusted, na.rm = TRUE),
      median_adjusted_inverse_beta = stats::median(.data$inverse_braak_beta_adjusted, na.rm = TRUE),
      median_cerad_adjusted_inverse_beta = stats::median(.data$inverse_braak_beta_cerad_adjusted, na.rm = TRUE),
      median_late_decline = stats::median(.data$protein_late_decline_magnitude, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(summary_tbl, file.path(audits_dir, "SuppFig4_group_summary.csv"))

  # Three panels: nuisance-adjusted Braak, CERAD-adjusted Braak, and all-protein scatter.
  panel_a <- sf4_panel_label(sf4_make_distribution_plot(
    plot_tbl, "inverse_braak_beta_adjusted", "Adjusted Braak", NULL,
    "Adjusted inverse Braak magnitude"), "A")
  panel_b <- sf4_panel_label(sf4_make_distribution_plot(
    plot_tbl, "inverse_braak_beta_cerad_adjusted", "CERAD-adjusted Braak", NULL,
    "CERAD-adjusted inverse Braak magnitude"), "B")
  panel_c <- sf4_panel_label(sf4_make_scatter_plot(plot_tbl), "C")
  panel_paths <- c(
    save_panel_set(panel_a, "SuppFig4A_covariate_adjusted_braak_distribution", output_dir = panels_dir),
    save_panel_set(panel_b, "SuppFig4B_cerad_adjusted_braak_distribution", output_dir = panels_dir),
    save_panel_set(panel_c, "SuppFig4C_late_decline_vs_cerad_adjusted_braak", output_dir = panels_dir)
  )
  sf4_spearman_result <- stats::cor.test(plot_tbl$protein_late_decline_magnitude,
    plot_tbl$inverse_braak_beta_cerad_adjusted, method = "spearman", exact = FALSE)
  # Compute outside tibble's column mask: the output column `test` must not
  # shadow the correlation-result object during evaluation.
  sf4_check_values <- c(
    sf4_wilcox_group_p(plot_tbl, "inverse_braak_beta_adjusted"),
    sf4_wilcox_group_p(plot_tbl, "inverse_braak_beta_cerad_adjusted"),
    unname(sf4_spearman_result$estimate), sf4_spearman_result$p.value)
  checks <- tibble::tibble(
    test = c("A_two_sided_Wilcoxon", "B_two_sided_Wilcoxon", "C_all_915_Spearman_rho", "C_all_915_Spearman_P"),
    value = sf4_check_values
  )
  if (frozen && any(abs(checks$value - c(2.0563661173078336e-4,
      4.822979248565039e-4, .32095895898382526, 2.278288414378184e-23)) > 1e-8)) {
    stop("S4 frozen source differs from the verified caption statistics.", call. = FALSE)
  }
  readr::write_csv(checks, file.path(audits_dir, "SuppFig4_manuscript_panel_statistics.csv"))
  composite <- (panel_a | panel_b) / panel_c +
    patchwork::plot_layout(heights = c(1, 1))

  save_plot_set(
    composite,
    "Supplementary_Figure_4_pathology_model_robustness",
    width = 12.4,
    height = 8.6,
    output_dir = figures_dir
  )

  run_summary <- c(
    "Supplementary Figure 4 run summary",
    paste0("generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "figure: pathology-model robustness of Hsp60/10 client vulnerability",
    "primary adjusted model: protein_z ~ braak_z + age_z + sex + pmi_z + TMT batch",
    "CERAD sensitivity model: protein_z ~ braak_z + cerad_z + age_z + sex + pmi_z + TMT batch",
    paste0("model source: ", if (frozen) "frozen SuppFig4_plotted_gene_level_values.csv" else protein_info$matrix_name),
    paste0("Hsp60/10 genes plotted: ", sum(plot_tbl$group == "Hsp60/10 clients", na.rm = TRUE)),
    paste0("background genes plotted: ", sum(plot_tbl$group == "Non-client mitochondrial proteins", na.rm = TRUE)),
    "group summary:",
    paste0(
      "  - ", summary_tbl$group,
      ": n=", summary_tbl$n,
      ", median inverse Braak=", signif(summary_tbl$median_inverse_braak, 4),
      ", median adjusted inverse beta=", signif(summary_tbl$median_adjusted_inverse_beta, 4),
      ", median CERAD-adjusted inverse beta=", signif(summary_tbl$median_cerad_adjusted_inverse_beta, 4),
      collapse = "\n"
    )
  )

  writeLines(run_summary, file.path(audits_dir, "SuppFig4_run_summary.txt"))

  message("Supplementary Figure 4 complete.")
  print(summary_tbl)

  invisible(list(
    plot = composite,
    plot_tbl = plot_tbl,
    adjusted_tbl = adjusted_tbl,
    summary_tbl = summary_tbl,
    panels = list(A = panel_a, B = panel_b, C = panel_c),
    panel_paths = panel_paths
  ))
}

supfig4_outputs <- make_supfig4_pathology_model_robustness(inputs)
