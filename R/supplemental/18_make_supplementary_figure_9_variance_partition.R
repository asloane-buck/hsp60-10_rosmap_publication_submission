# ============================================================
# Supplementary Figure 9:
# Source-of-variation analysis
# ============================================================

if (!exists("project_dir")) {
  stop("Run 00_supplemental_config.R before this script.", call. = FALSE)
}

if (!exists("save_plot_set")) {
  stop("Run 02_supplemental_helper_functions.R before this script.", call. = FALSE)
}

sf9_input_dir <- file.path(
  project_dir,
  "outputs",
  "reviewer_revisions",
  "variancePartition_source_of_variation"
)

sf9_summary_file <- file.path(
  sf9_input_dir,
  "MASTER_variance_fraction_summary.csv"
)

if (!file.exists(sf9_summary_file)) {
  stop("Missing S9 variance summary.", call. = FALSE)
}

sf9_raw <- readr::read_csv(
  sf9_summary_file,
  show_col_types = FALSE
)

required_cols <- c(
  "modality",
  "source",
  "n_features",
  "median_fraction",
  "q25_fraction",
  "q75_fraction"
)

missing_cols <- setdiff(required_cols, colnames(sf9_raw))

if (length(missing_cols) > 0) {
  stop(
    "S9 source table missing columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

sf9_source_label <- function(x) {
  dplyr::recode(
    as.character(x),
    sequencing_batch = "Sequencing batch",
    rin_z = "RIN",
    age_z = "Age",
    sex_binary = "Sex",
    pmi_z = "PMI",
    tmt_batch = "TMT batch",
    Residuals = "Residual",
    .default = as.character(x)
  )
}

sf9 <- sf9_raw |>
  dplyr::mutate(
    modality = dplyr::case_when(
      stringr::str_to_lower(.data$modality) == "rna" ~ "RNA",
      stringr::str_to_lower(.data$modality) == "protein" ~ "Protein",
      TRUE ~ as.character(.data$modality)
    ),
    source_label = sf9_source_label(.data$source),
    is_residual = stringr::str_detect(
      stringr::str_to_lower(.data$source),
      "^resid"
    )
  ) |>
  dplyr::filter(.data$modality %in% c("RNA", "Protein"))

sf9_residual <- sf9 |>
  dplyr::filter(.data$is_residual)

sf9_nonresidual <- sf9 |>
  dplyr::filter(!.data$is_residual)

rna_residual <- sf9_residual |>
  dplyr::filter(.data$modality == "RNA") |>
  dplyr::pull(.data$median_fraction)

protein_residual <- sf9_residual |>
  dplyr::filter(.data$modality == "Protein") |>
  dplyr::pull(.data$median_fraction)

if (length(rna_residual) != 1 || length(protein_residual) != 1) {
  stop("Could not uniquely identify residual fractions.", call. = FALSE)
}

rna_order <- c(
  "Sequencing batch",
  "RIN",
  "PMI",
  "Sex",
  "Age"
)

protein_order <- c(
  "PMI",
  "Age",
  "Sex",
  "TMT batch"
)

sf9_rna <- sf9_nonresidual |>
  dplyr::filter(.data$modality == "RNA") |>
  dplyr::mutate(
    source_label = factor(
      .data$source_label,
      levels = rev(rna_order)
    ),
    median_pct = 100 * .data$median_fraction,
    q25_pct = 100 * .data$q25_fraction,
    q75_pct = 100 * .data$q75_fraction
  )

sf9_protein <- sf9_nonresidual |>
  dplyr::filter(.data$modality == "Protein") |>
  dplyr::mutate(
    source_label = factor(
      .data$source_label,
      levels = rev(protein_order)
    ),
    median_pct = 100 * .data$median_fraction,
    q25_pct = 100 * .data$q25_fraction,
    q75_pct = 100 * .data$q75_fraction
  )

sf9_table_dir <- file.path(
  tables_dir,
  "supplementary_figure_9_variance_partition"
)

ensure_dir(sf9_table_dir)

readr::write_csv(
  sf9,
  file.path(
    sf9_table_dir,
    "Supplementary_Figure_9_variance_fraction_summary.csv"
  )
)

readr::write_csv(
  sf9_rna,
  file.path(
    sf9_table_dir,
    "Supplementary_Figure_9A_RNA_sources.csv"
  )
)

readr::write_csv(
  sf9_protein,
  file.path(
    sf9_table_dir,
    "Supplementary_Figure_9B_protein_sources.csv"
  )
)

sf9_theme <- function(base_size = 8.5) {
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
        color = "grey30",
        linewidth = 0.4
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(
        color = "grey35",
        linewidth = 0.3
      ),
      axis.text = ggplot2::element_text(
        color = "grey15"
      ),
      axis.title = ggplot2::element_text(
        color = "grey10",
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 9.5,
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        size = 8,
        color = "grey30",
        hjust = 0,
        margin = ggplot2::margin(b = 5)
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.25
      ),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
}

sf9_panel_label <- function(plot, label) {
  plot +
    ggplot2::labs(tag = label) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 9.5
      ),
      plot.tag.position = c(-0.05, 1.07)
    )
}

sf9_rna_label_x <- 43
sf9_protein_label_x <- 1.37

sf9_panel_a <- ggplot2::ggplot(
  sf9_rna,
  ggplot2::aes(
    y = .data$source_label
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$q25_pct,
      xend = .data$q75_pct,
      yend = .data$source_label
    ),
    linewidth = 1.0,
    color = "#4C78A8"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = .data$median_pct),
    size = 3.0,
    color = "#4C78A8"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = sf9_rna_label_x,
      label = sprintf("%.1f%%", .data$median_pct)
    ),
    hjust = 1,
    size = 2.45,
    color = "#111827"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-1.5, 45),
    breaks = c(0, 10, 20, 30, 40),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "RNA",
    subtitle = paste0(
      "Residual median: ",
      sprintf("%.1f%%", 100 * rna_residual),
      "   |   n = 21,433 features"
    ),
    x = "Variance fraction",
    y = NULL
  ) +
  sf9_theme() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      12, 14, 8, 18
    )
  )

sf9_panel_a <- sf9_panel_label(
  sf9_panel_a,
  "A"
)

sf9_panel_b <- ggplot2::ggplot(
  sf9_protein,
  ggplot2::aes(
    y = .data$source_label
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$q25_pct,
      xend = .data$q75_pct,
      yend = .data$source_label
    ),
    linewidth = 1.0,
    color = "#B23A48"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = .data$median_pct),
    size = 3.0,
    color = "#B23A48"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = sf9_protein_label_x,
      label = sprintf("%.1f%%", .data$median_pct)
    ),
    hjust = 1,
    size = 2.45,
    color = "#111827"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-0.06, 1.45),
    breaks = c(0, 0.4, 0.8, 1.2),
    labels = function(x) paste0(
      formatC(x, format = "f", digits = 1),
      "%"
    ),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Protein",
    subtitle = paste0(
      "Residual median: ",
      sprintf("%.1f%%", 100 * protein_residual),
      "   |   n = 5,088 proteins"
    ),
    x = "Variance fraction",
    y = NULL
  ) +
  sf9_theme() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      12, 14, 8, 18
    )
  )

sf9_panel_b <- sf9_panel_label(
  sf9_panel_b,
  "B"
)

save_panel_set(
  sf9_panel_a,
  "Supplementary_Figure_9A_RNA_variance_sources",
  width = 4.0,
  height = 2.8
)

save_panel_set(
  sf9_panel_b,
  "Supplementary_Figure_9B_protein_variance_sources",
  width = 4.0,
  height = 2.8
)

sf9_final <- (
  sf9_panel_a |
    sf9_panel_b
) +
  patchwork::plot_layout(
    widths = c(1, 1)
  )

save_plot_set(
  sf9_final,
  "Supplementary_Figure_9_variance_partition",
  width = 6.7,
  height = 2.8,
  output_dir = figures_dir
)

sf9_audit <- tibble::tibble(
  check = c(
    "summary_file_exists",
    "RNA_sources_present",
    "protein_sources_present",
    "RNA_residual_unique",
    "protein_residual_unique"
  ),
  pass = c(
    file.exists(sf9_summary_file),
    nrow(sf9_rna) == 5,
    nrow(sf9_protein) == 4,
    length(rna_residual) == 1,
    length(protein_residual) == 1
  )
)

readr::write_csv(
  sf9_audit,
  file.path(
    audits_dir,
    "Supplementary_Figure_9_variance_partition_audit.csv"
  )
)

if (!all(sf9_audit$pass)) {
  print(sf9_audit, n = Inf)
  stop("Supplementary Figure 9 audit failed.", call. = FALSE)
}

supfig9_outputs <- list(
  RNA = sf9_rna,
  protein = sf9_protein,
  panel_a = sf9_panel_a,
  panel_b = sf9_panel_b,
  final = sf9_final,
  audit = sf9_audit
)

message("Supplementary Figure 9 complete.")
