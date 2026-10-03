# ============================================================
# Supplementary Figure 7:
# Conventional RNA differential expression and protein
# differential abundance across diagnostic-stage contrasts
# ============================================================

if (!exists("project_dir")) {
  stop(
    "Run 00_supplemental_config.R before this script.",
    call. = FALSE
  )
}

if (!exists("save_plot_set")) {
  stop(
    "Run 02_supplemental_helper_functions.R before this script.",
    call. = FALSE
  )
}

# ============================================================
# Paths
# ============================================================

sf7_input_dir <- file.path(
  project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential_final"
)

sf7_summary_file <- file.path(
  sf7_input_dir,
  "18_PRIMARY_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv"
)

sf7_side_by_side_file <- file.path(
  sf7_input_dir,
  "17_PRIMARY_cross_modal_Hsp60_10_side_by_side_NO_SUBTRACTION.csv"
)

required_files <- c(
  sf7_summary_file,
  sf7_side_by_side_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required S7 input(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ============================================================
# Style
# ============================================================

sf7_rna_col <- "#4C78A8"
sf7_protein_col <- "#B23A48"
sf7_both_col <- "#7A5195"
sf7_neither_col <- "#D9DDE3"
sf7_text <- "#111827"

sf7_class_colors <- c(
  "RNA only" = sf7_rna_col,
  "Protein only" = sf7_protein_col,
  "Both" = sf7_both_col,
  "Neither" = sf7_neither_col
)

sf7_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "grey25",
        linewidth = 0.45
      ),
      axis.line = ggplot2::element_line(
        color = "grey25",
        linewidth = 0.35
      ),
      axis.ticks = ggplot2::element_line(
        color = "grey25",
        linewidth = 0.35
      ),
      axis.text = ggplot2::element_text(
        color = "grey15"
      ),
      axis.title = ggplot2::element_text(
        color = "grey10",
        face = "bold"
      ),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.25
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

sf7_panel_label <- function(plot, label) {
  plot +
    ggplot2::labs(tag = label) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 9.5
      ),
      plot.tag.position = c(0, 1)
    )
}

# ============================================================
# Load canonical validated outputs
# ============================================================

sf7_summary_raw <- readr::read_csv(
  sf7_summary_file,
  show_col_types = FALSE
)

sf7_side_by_side_raw <- readr::read_csv(
  sf7_side_by_side_file,
  show_col_types = FALSE
)

# ============================================================
# Validate source schema
# ============================================================

summary_required <- c(
  "contrast",
  "n_shared_hsp_clients",
  "n_rna_only_fdr",
  "n_protein_only_fdr",
  "n_both_fdr",
  "n_neither_fdr",
  "n_direction_concordant",
  "pct_direction_concordant"
)

side_required <- c(
  "contrast",
  "gene_symbol",
  "rna_effect",
  "rna_std_error",
  "rna_conf_low",
  "rna_conf_high",
  "rna_p_value",
  "rna_fdr",
  "rna_significant_fdr05",
  "protein_effect",
  "protein_std_error",
  "protein_conf_low",
  "protein_conf_high",
  "protein_p_value",
  "protein_fdr",
  "protein_significant_fdr05",
  "significance_class",
  "direction_concordant"
)

missing_summary_cols <- setdiff(
  summary_required,
  colnames(sf7_summary_raw)
)

missing_side_cols <- setdiff(
  side_required,
  colnames(sf7_side_by_side_raw)
)

if (length(missing_summary_cols) > 0) {
  stop(
    "S7 summary missing columns: ",
    paste(missing_summary_cols, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_side_cols) > 0) {
  stop(
    "S7 side-by-side table missing columns: ",
    paste(missing_side_cols, collapse = ", "),
    call. = FALSE
  )
}

# ============================================================
# Display order
# ============================================================

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

sf7_summary <- sf7_summary_raw |>
  dplyr::filter(.data$contrast %in% contrast_levels) |>
  dplyr::mutate(
    contrast = factor(
      .data$contrast,
      levels = contrast_levels
    )
  ) |>
  dplyr::arrange(.data$contrast)

if (nrow(sf7_summary) != 3) {
  stop("Expected exactly three primary contrasts.", call. = FALSE)
}

if (!all(sf7_summary$n_shared_hsp_clients == 257)) {
  stop(
    "Expected 257 shared Hsp60/10 clients per contrast.",
    call. = FALSE
  )
}

sf7_count_check <- sf7_summary |>
  dplyr::mutate(
    category_sum =
      .data$n_rna_only_fdr +
      .data$n_protein_only_fdr +
      .data$n_both_fdr +
      .data$n_neither_fdr
  )

if (!all(
  sf7_count_check$category_sum ==
    sf7_count_check$n_shared_hsp_clients
)) {
  stop(
    "Significance counts do not sum to shared-client universe.",
    call. = FALSE
  )
}

sf7_side_count_check <- sf7_side_by_side_raw |>
  dplyr::filter(.data$contrast %in% contrast_levels) |>
  dplyr::count(.data$contrast, name = "n_rows")

if (!all(sf7_side_count_check$n_rows == 257)) {
  stop(
    "Expected 257 rows per contrast in side-by-side table.",
    call. = FALSE
  )
}

# ============================================================
# Panel A source
# ============================================================

sf7_panel_a_source <- sf7_summary |>
  dplyr::select(
    "contrast",
    "n_shared_hsp_clients",
    "n_rna_only_fdr",
    "n_protein_only_fdr",
    "n_both_fdr",
    "n_neither_fdr"
  ) |>
  tidyr::pivot_longer(
    cols = c(
      "n_rna_only_fdr",
      "n_protein_only_fdr",
      "n_both_fdr",
      "n_neither_fdr"
    ),
    names_to = "significance_class",
    values_to = "n_clients"
  ) |>
  dplyr::mutate(
    significance_class = dplyr::recode(
      .data$significance_class,
      n_rna_only_fdr = "RNA only",
      n_protein_only_fdr = "Protein only",
      n_both_fdr = "Both",
      n_neither_fdr = "Neither"
    ),
    significance_class = factor(
      .data$significance_class,
      levels = c(
        "Neither",
        "Protein only",
        "RNA only",
        "Both"
      )
    ),
    contrast_label = contrast_labels[
      as.character(.data$contrast)
    ],
    contrast_label = factor(
      .data$contrast_label,
      levels = unname(
        contrast_labels[contrast_levels]
      )
    )
  )

# ============================================================
# Panel B source
# ============================================================

sf7_panel_b_source <- sf7_summary |>
  dplyr::transmute(
    contrast = .data$contrast,
    contrast_label = contrast_labels[
      as.character(.data$contrast)
    ],
    n_shared_hsp_clients = .data$n_shared_hsp_clients,
    n_direction_concordant =
      .data$n_direction_concordant,
    pct_direction_concordant =
      .data$pct_direction_concordant,
    count_label = paste0(
      .data$n_direction_concordant,
      "/",
      .data$n_shared_hsp_clients
    )
  ) |>
  dplyr::mutate(
    contrast_label = factor(
      .data$contrast_label,
      levels = unname(
        contrast_labels[contrast_levels]
      )
    )
  )

# ============================================================
# Write source tables
# ============================================================

sf7_table_dir <- file.path(
  tables_dir,
  "supplementary_figure_7_conventional_differential"
)

ensure_dir(sf7_table_dir)

readr::write_csv(
  sf7_panel_a_source,
  file.path(
    sf7_table_dir,
    "Supplementary_Figure_7A_significance_categories.csv"
  )
)

readr::write_csv(
  sf7_panel_b_source,
  file.path(
    sf7_table_dir,
    "Supplementary_Figure_7B_direction_concordance.csv"
  )
)

readr::write_csv(
  sf7_side_by_side_raw |>
    dplyr::filter(.data$contrast %in% contrast_levels) |>
    dplyr::mutate(
      contrast = factor(
        .data$contrast,
        levels = contrast_levels
      )
    ) |>
    dplyr::arrange(
      .data$contrast,
      .data$gene_symbol
    ),
  file.path(
    sf7_table_dir,
    "Supplementary_Table_conventional_RNA_protein_clients.csv"
  )
)

# ============================================================
# Panel A
# ============================================================

sf7_panel_a_label_source <- sf7_panel_a_source |>
  dplyr::filter(.data$n_clients >= 20)

sf7_panel_a_callouts <- tibble::tribble(
  ~contrast_label, ~y, ~label,
  "AD vs MCI", 50, "RNA only 7; both 1",
  "AD vs NCI", 124, "RNA only 6; both 10"
) |>
  dplyr::mutate(
    contrast_label = factor(
      .data$contrast_label,
      levels = unname(
        contrast_labels[contrast_levels]
      )
    )
  )

sf7_panel_a <- ggplot2::ggplot(
  sf7_panel_a_source,
  ggplot2::aes(
    x = .data$contrast_label,
    y = .data$n_clients,
    fill = .data$significance_class
  )
) +
  ggplot2::geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    data = sf7_panel_a_label_source,
    ggplot2::aes(label = .data$n_clients),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 2.5,
    color = sf7_text
  ) +
  ggplot2::geom_text(
    data = sf7_panel_a_callouts,
    ggplot2::aes(
      x = .data$contrast_label,
      y = .data$y,
      label = .data$label
    ),
    inherit.aes = FALSE,
    size = 1.95,
    fontface = "plain",
    color = "grey25"
  ) +
  ggplot2::scale_fill_manual(
    values = sf7_class_colors,
    breaks = c(
      "RNA only",
      "Protein only",
      "Both",
      "Neither"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 270),
    breaks = c(0, 50, 100, 150, 200, 250),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Shared Hsp60/10 clients",
    fill = NULL
  ) +
  sf7_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.key.width = grid::unit(0.55, "cm"),
    legend.key.height = grid::unit(0.40, "cm"),
    legend.spacing.x = grid::unit(0.12, "cm"),
    legend.margin = ggplot2::margin(0, 0, 6, 0),
    plot.margin = ggplot2::margin(12, 8, 6, 12)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

sf7_panel_a <- sf7_panel_label(
  sf7_panel_a,
  "A"
)


# ============================================================
# Panel B
# ============================================================

sf7_panel_b <- ggplot2::ggplot(
  sf7_panel_b_source,
  ggplot2::aes(
    x = .data$contrast_label,
    y = .data$pct_direction_concordant
  )
) +
  ggplot2::geom_col(
    width = 0.62,
    fill = "#6B7280"
  ) +
  ggplot2::geom_text(
    data = sf7_panel_b_source |>
      dplyr::filter(.data$contrast != "MCI_vs_NCI"),
    ggplot2::aes(
      label = sprintf(
        "%.1f%%\n%s",
        .data$pct_direction_concordant,
        .data$count_label
      )
    ),
    vjust = -0.28,
    size = 2.3,
    color = sf7_text,
    lineheight = 0.95
  ) +
  ggplot2::geom_text(
    data = sf7_panel_b_source |>
      dplyr::filter(.data$contrast == "MCI_vs_NCI"),
    ggplot2::aes(
      y = .data$pct_direction_concordant - 1.8,
      label = sprintf(
        "%.1f%%\n%s",
        .data$pct_direction_concordant,
        .data$count_label
      )
    ),
    vjust = 1,
    size = 2.3,
    color = sf7_text,
    lineheight = 0.95
  ) +
  ggplot2::geom_hline(
    yintercept = 50,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey60"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 90),
    breaks = c(0, 20, 40, 60, 80),
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "RNA/protein effect-direction concordance"
  ) +
  sf7_theme() +
  ggplot2::theme(
    legend.position = "none"
  )

sf7_panel_b <- sf7_panel_label(
  sf7_panel_b,
  "B"
)

# ============================================================
# Save panels
# ============================================================

save_panel_set(
  sf7_panel_a,
  "Supplementary_Figure_7A_significance_categories",
  width = 4,
  height = 3.4
)

save_panel_set(
  sf7_panel_b,
  "Supplementary_Figure_7B_direction_concordance",
  width = 4,
  height = 3.4
)

# ============================================================
# Assemble final
# ============================================================

sf7_final <- (
  sf7_panel_a |
    sf7_panel_b
) +
  patchwork::plot_layout(
    widths = c(1.08, 0.92)
  ) &
  ggplot2::theme(
    plot.margin = ggplot2::margin(12, 10, 8, 10)
  )

save_plot_set(
  sf7_final,
  "Supplementary_Figure_7_conventional_differential",
  width = 6.7,
  height = 3.55,
  output_dir = figures_dir
)

# ============================================================
# Audit
# ============================================================

sf7_audit <- tibble::tibble(
  check = c(
    "summary_file_exists",
    "side_by_side_file_exists",
    "three_primary_contrasts",
    "shared_clients_257_each",
    "category_counts_sum_to_257",
    "side_by_side_rows_257_each"
  ),
  pass = c(
    file.exists(sf7_summary_file),
    file.exists(sf7_side_by_side_file),
    nrow(sf7_summary) == 3,
    all(sf7_summary$n_shared_hsp_clients == 257),
    all(
      sf7_count_check$category_sum ==
        sf7_count_check$n_shared_hsp_clients
    ),
    all(sf7_side_count_check$n_rows == 257)
  )
)

readr::write_csv(
  sf7_audit,
  file.path(
    audits_dir,
    "Supplementary_Figure_7_conventional_differential_audit.csv"
  )
)

if (!all(sf7_audit$pass)) {
  print(sf7_audit, n = Inf)
  stop("Supplementary Figure 7 audit failed.", call. = FALSE)
}

supfig7_outputs <- list(
  summary = sf7_summary,
  panel_a_source = sf7_panel_a_source,
  panel_b_source = sf7_panel_b_source,
  side_by_side = sf7_side_by_side_raw,
  panel_a = sf7_panel_a,
  panel_b = sf7_panel_b,
  final = sf7_final,
  audit = sf7_audit
)

message("Supplementary Figure 7 complete.")
