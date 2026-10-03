# ============================================================
# Supplementary Figure 10:
# Alternative-mechanism marker abundance analysis
# ============================================================

if (!exists("project_dir")) {
  stop("Run 00_supplemental_config.R before this script.", call. = FALSE)
}

if (!exists("save_plot_set")) {
  stop("Run 02_supplemental_helper_functions.R before this script.", call. = FALSE)
}

sf10_input_dir <- file.path(
  project_dir,
  "outputs",
  "reviewer_revisions",
  "alternative_mechanism_marker_panel"
)

sf10_results_file <- file.path(
  sf10_input_dir,
  "marker_stage_limma_results.csv"
)

sf10_summary_file <- file.path(
  sf10_input_dir,
  "marker_category_stage_summary.csv"
)

sf10_undetected_file <- file.path(
  sf10_input_dir,
  "undetected_prespecified_markers.csv"
)

required_files <- c(
  sf10_results_file,
  sf10_summary_file,
  sf10_undetected_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required S10 input(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

sf10_raw <- readr::read_csv(
  sf10_results_file,
  show_col_types = FALSE
)

sf10_category_summary <- readr::read_csv(
  sf10_summary_file,
  show_col_types = FALSE
)

sf10_undetected <- readr::read_csv(
  sf10_undetected_file,
  show_col_types = FALSE
)

required_cols <- c(
  "category",
  "gene",
  "contrast",
  "effect",
  "marker_panel_fdr",
  "significant_marker_panel_fdr05"
)

missing_cols <- setdiff(
  required_cols,
  colnames(sf10_raw)
)

if (length(missing_cols) > 0) {
  stop(
    "S10 source table missing columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

category_levels <- c(
  "mitochondrial_mass_import",
  "mitochondrial_biogenesis",
  "mitophagy",
  "autophagy_lysosome",
  "proteasome_20S_core"
)

category_labels <- c(
  "mitochondrial_mass_import" = "Mass / import",
  "mitochondrial_biogenesis" = "Biogenesis",
  "mitophagy" = "Mitophagy",
  "autophagy_lysosome" = "Autophagy / lysosome",
  "proteasome_20S_core" = "20S proteasome"
)

contrast_levels <- c(
  "MCI_vs_NCI",
  "AD_vs_MCI",
  "AD_vs_NCI"
)

contrast_labels <- c(
  "MCI_vs_NCI" = "MCI vs NCI",
  "AD_vs_MCI" = "AD vs MCI",
  "AD_vs_NCI" = "AD vs NCI"
)

sf10 <- sf10_raw |>
  dplyr::filter(
    .data$category %in% category_levels,
    .data$contrast %in% contrast_levels
  ) |>
  dplyr::mutate(
    category = factor(
      .data$category,
      levels = category_levels
    ),
    category_label = category_labels[
      as.character(.data$category)
    ],
    category_label = factor(
      .data$category_label,
      levels = unname(
        category_labels[category_levels]
      )
    ),
    contrast_label = contrast_labels[
      .data$contrast
    ],
    contrast_label = factor(
      .data$contrast_label,
      levels = unname(
        contrast_labels[contrast_levels]
      )
    )
  )

gene_order_tbl <- sf10 |>
  dplyr::distinct(
    .data$category,
    .data$category_label,
    .data$gene
  ) |>
  dplyr::arrange(
    .data$category,
    .data$gene
  ) |>
  dplyr::mutate(
    gene_order = dplyr::row_number()
  )

gene_levels <- rev(gene_order_tbl$gene)

sf10 <- sf10 |>
  dplyr::mutate(
    gene = factor(
      .data$gene,
      levels = gene_levels
    )
  )

n_detected_genes <- dplyr::n_distinct(sf10$gene)
n_tests <- nrow(sf10)
n_fdr <- sum(
  sf10$significant_marker_panel_fdr05,
  na.rm = TRUE
)

if (n_detected_genes != 39) {
  stop(
    "Expected 39 detected marker genes; found ",
    n_detected_genes,
    ".",
    call. = FALSE
  )
}

if (n_tests != 117) {
  stop(
    "Expected 117 marker-by-contrast tests; found ",
    n_tests,
    ".",
    call. = FALSE
  )
}

if (n_fdr != 9) {
  stop(
    "Expected 9 FDR-significant marker tests; found ",
    n_fdr,
    ".",
    call. = FALSE
  )
}

sig_by_contrast <- sf10 |>
  dplyr::filter(
    .data$significant_marker_panel_fdr05
  ) |>
  dplyr::count(
    .data$contrast,
    name = "n_fdr"
  )

if (
  nrow(sig_by_contrast) != 1 ||
  as.character(sig_by_contrast$contrast[1]) != "AD_vs_NCI" ||
  sig_by_contrast$n_fdr[1] != 9
) {
  stop(
    "Expected all 9 FDR-significant marker tests in AD vs NCI.",
    call. = FALSE
  )
}

sf10_table_dir <- file.path(
  tables_dir,
  "supplementary_figure_10_marker_sensitivity"
)

ensure_dir(sf10_table_dir)

readr::write_csv(
  sf10 |>
    dplyr::arrange(
      .data$category,
      .data$gene,
      .data$contrast_label
    ),
  file.path(
    sf10_table_dir,
    "Supplementary_Figure_10_marker_results.csv"
  )
)

readr::write_csv(
  sf10_category_summary,
  file.path(
    sf10_table_dir,
    "Supplementary_Figure_10_category_summary.csv"
  )
)

readr::write_csv(
  sf10_undetected,
  file.path(
    sf10_table_dir,
    "Supplementary_Figure_10_undetected_markers.csv"
  )
)

max_abs_effect <- max(
  abs(sf10$effect),
  na.rm = TRUE
)

fill_limit <- ceiling(
  max_abs_effect * 100
) / 100

fill_limit <- max(
  fill_limit,
  0.12
)

sf10_theme <- function(base_size = 8.3) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        color = "grey15",
        face = "bold",
        size = 8.2,
        margin = ggplot2::margin(t = 4)
      ),
      axis.text.y = ggplot2::element_text(
        color = "grey15",
        size = 7.5
      ),
      axis.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "grey94",
        color = "grey75",
        linewidth = 0.35
      ),
      strip.text.y = ggplot2::element_text(
        face = "bold",
        color = "grey15",
        size = 7.5
      ),
      strip.text.y.left = ggplot2::element_text(
        face = "bold",
        color = "grey15",
        angle = 90,
        size = 7.5
      ),
      strip.placement = "outside",
      panel.spacing.y = grid::unit(
        0.05,
        "lines"
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 8
      ),
      legend.text = ggplot2::element_text(
        size = 7.5
      ),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(
        8, 10, 8, 10
      )
    )
}

sf10_plot <- ggplot2::ggplot(
  sf10,
  ggplot2::aes(
    x = .data$contrast_label,
    y = .data$gene,
    fill = .data$effect
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.5,
    width = 0.92,
    height = 0.92
  ) +
  ggplot2::geom_text(
    data = sf10 |>
      dplyr::filter(
        .data$significant_marker_panel_fdr05
      ),
    ggplot2::aes(label = "*"),
    color = "black",
    size = 3.2,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#B23A48",
    midpoint = 0,
    limits = c(
      -fill_limit,
      fill_limit
    ),
    oob = scales::squish,
    name = "limma effect\n(log2 abundance)"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(
      category_label
    ),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  sf10_theme() +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(
        0.32,
        "cm"
      ),
      barwidth = grid::unit(
        5.0,
        "cm"
      )
    )
  )

save_panel_set(
  sf10_plot,
  "Supplementary_Figure_10_marker_heatmap",
  width = 6.7,
  height = 7.0
)

save_plot_set(
  sf10_plot,
  "Supplementary_Figure_10_marker_sensitivity",
  width = 6.7,
  height = 7.0,
  output_dir = figures_dir
)

sf10_audit <- tibble::tibble(
  check = c(
    "results_file_exists",
    "summary_file_exists",
    "undetected_file_exists",
    "39_detected_markers",
    "117_marker_tests",
    "9_FDR_significant_tests",
    "all_FDR_significant_in_AD_vs_NCI",
    "5_undetected_prespecified_markers"
  ),
  pass = c(
    file.exists(sf10_results_file),
    file.exists(sf10_summary_file),
    file.exists(sf10_undetected_file),
    n_detected_genes == 39,
    n_tests == 117,
    n_fdr == 9,
    nrow(sig_by_contrast) == 1 &&
      as.character(
        sig_by_contrast$contrast[1]
      ) == "AD_vs_NCI" &&
      sig_by_contrast$n_fdr[1] == 9,
    nrow(sf10_undetected) == 5
  )
)

readr::write_csv(
  sf10_audit,
  file.path(
    audits_dir,
    "Supplementary_Figure_10_marker_sensitivity_audit.csv"
  )
)

if (!all(sf10_audit$pass)) {
  print(sf10_audit, n = Inf)
  stop(
    "Supplementary Figure 10 audit failed.",
    call. = FALSE
  )
}

supfig10_outputs <- list(
  marker_results = sf10,
  category_summary = sf10_category_summary,
  undetected = sf10_undetected,
  plot = sf10_plot,
  audit = sf10_audit
)

message("Supplementary Figure 10 complete.")
