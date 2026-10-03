# ============================================================
# Supplementary Figure 8:
# PC1 sensitivity of the Hsp60/10 RNA client score
# ============================================================
#
# Purpose:
# Compare the production mean client z-score with an oriented
# PC1/eigengene sensitivity score computed from the same RNA
# Hsp60/10 client matrix.
#
# This script does NOT recompute PCA or stage tests.
# It reads validated reviewer-revision outputs produced by:
#   R/reviewer_revisions/44_run_Hsp_pathway_PC1_sensitivity.R
#
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

sf8_input_dir <- file.path(
  project_dir,
  "outputs",
  "reviewer_revisions",
  "Hsp_pathway_PC1_sensitivity"
)

sf8_participant_file <- file.path(
  sf8_input_dir,
  "participant_mean_z_and_PC1_scores.csv"
)

sf8_stage_summary_file <- file.path(
  sf8_input_dir,
  "stage_score_summary.csv"
)

sf8_stage_test_file <- file.path(
  sf8_input_dir,
  "stage_Wilcoxon_mean_z_vs_PC1.csv"
)

sf8_concordance_file <- file.path(
  sf8_input_dir,
  "mean_z_vs_PC1_concordance.csv"
)

sf8_required_files <- c(
  sf8_participant_file,
  sf8_stage_summary_file,
  sf8_stage_test_file,
  sf8_concordance_file
)

sf8_missing_files <- sf8_required_files[
  !file.exists(sf8_required_files)
]

if (length(sf8_missing_files) > 0) {
  stop(
    "Missing required S8 input(s):\n",
    paste(sf8_missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ============================================================
# Style
# ============================================================

sf8_mean_col <- "#4C78A8"
sf8_pc1_col <- "#7A5195"
sf8_point_col <- "#6B7280"
sf8_text <- "#111827"

sf8_theme <- function(base_size = 8.5) {
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
      strip.background = ggplot2::element_rect(
        fill = "grey96",
        color = "grey65",
        linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "grey15"
      ),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.25
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

sf8_panel_label <- function(plot, label) {
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
# Load data
# ============================================================

sf8_participant <- readr::read_csv(
  sf8_participant_file,
  show_col_types = FALSE
)

sf8_stage_summary <- readr::read_csv(
  sf8_stage_summary_file,
  show_col_types = FALSE
)

sf8_stage_tests <- readr::read_csv(
  sf8_stage_test_file,
  show_col_types = FALSE
)

sf8_concordance <- readr::read_csv(
  sf8_concordance_file,
  show_col_types = FALSE
)

# ============================================================
# Validate schemas
# ============================================================

participant_required <- c(
  "sample_id",
  "clinical_stage",
  "mean_z",
  "PC1_eigengene"
)

summary_required <- c(
  "score_type",
  "clinical_stage",
  "n",
  "mean_score",
  "sd_score",
  "sem",
  "median_score"
)

test_required <- c(
  "score_type",
  "comparison",
  "group1",
  "group2",
  "n1",
  "n2",
  "delta_mean",
  "p_value",
  "p_fdr_bh_within_score"
)

if (length(setdiff(
  participant_required,
  colnames(sf8_participant)
)) > 0) {
  stop("Participant PC1 source table schema mismatch.", call. = FALSE)
}

if (length(setdiff(
  summary_required,
  colnames(sf8_stage_summary)
)) > 0) {
  stop("Stage summary schema mismatch.", call. = FALSE)
}

if (length(setdiff(
  test_required,
  colnames(sf8_stage_tests)
)) > 0) {
  stop("Stage test schema mismatch.", call. = FALSE)
}

# ============================================================
# Pull frozen concordance metrics
# ============================================================

sf8_metric <- function(name) {
  x <- sf8_concordance |>
    dplyr::filter(.data$metric == name) |>
    dplyr::pull(.data$value)

  if (length(x) != 1) {
    stop(
      "Could not uniquely identify concordance metric: ",
      name,
      call. = FALSE
    )
  }

  as.numeric(x)
}

sf8_pearson <- sf8_metric(
  "pearson_r_mean_z_vs_PC1"
)

sf8_spearman <- sf8_metric(
  "spearman_rho_mean_z_vs_PC1"
)

sf8_pc1_var <- sf8_metric(
  "PC1_variance_explained"
)

# ============================================================
# Validate participant counts
# ============================================================

sf8_stage_levels <- c(
  "NCI",
  "MCI",
  "AD"
)

sf8_participant <- sf8_participant |>
  dplyr::filter(
    is.finite(.data$mean_z),
    is.finite(.data$PC1_eigengene),
    .data$clinical_stage %in% sf8_stage_levels
  ) |>
  dplyr::mutate(
    clinical_stage = factor(
      .data$clinical_stage,
      levels = sf8_stage_levels
    )
  )

sf8_n <- nrow(sf8_participant)

if (sf8_n != 577) {
  stop(
    "Expected 577 participants in PC1 sensitivity table; found ",
    sf8_n,
    ".",
    call. = FALSE
  )
}

# ============================================================
# Source tables
# ============================================================

sf8_stage_plot <- sf8_stage_summary |>
  dplyr::filter(
    .data$score_type %in% c(
      "mean_z",
      "PC1_eigengene"
    ),
    .data$clinical_stage %in% sf8_stage_levels
  ) |>
  dplyr::mutate(
    clinical_stage = factor(
      .data$clinical_stage,
      levels = sf8_stage_levels
    ),
    score_label = dplyr::recode(
      .data$score_type,
      mean_z = "Mean client z-score",
      PC1_eigengene = "Oriented PC1"
    ),
    score_label = factor(
      .data$score_label,
      levels = c(
        "Mean client z-score",
        "Oriented PC1"
      )
    )
  )

sf8_mci_ad_tests <- sf8_stage_tests |>
  dplyr::filter(
    .data$comparison == "MCI_vs_AD",
    .data$score_type %in% c(
      "mean_z",
      "PC1_eigengene"
    )
  ) |>
  dplyr::mutate(
    score_label = dplyr::recode(
      .data$score_type,
      mean_z = "Mean client z-score",
      PC1_eigengene = "Oriented PC1"
    ),
    annotation = dplyr::case_when(
      .data$p_fdr_bh_within_score < 0.001 ~
        "MCI vs AD\nBH-adjusted P < 0.001",
      TRUE ~ paste0(
        "MCI vs AD\nBH-adjusted P = ",
        formatC(
          .data$p_fdr_bh_within_score,
          format = "f",
          digits = 3
        )
      )
    )
  )

# ============================================================
# Write frozen source tables
# ============================================================

sf8_table_dir <- file.path(
  tables_dir,
  "supplementary_figure_8_PC1_sensitivity"
)

ensure_dir(sf8_table_dir)

readr::write_csv(
  sf8_participant,
  file.path(
    sf8_table_dir,
    "Supplementary_Figure_8A_participant_scores.csv"
  )
)

readr::write_csv(
  sf8_stage_plot,
  file.path(
    sf8_table_dir,
    "Supplementary_Figure_8B_stage_summary.csv"
  )
)

readr::write_csv(
  sf8_mci_ad_tests,
  file.path(
    sf8_table_dir,
    "Supplementary_Figure_8B_MCI_AD_tests.csv"
  )
)

# ============================================================
# Panel A:
# mean-z versus oriented PC1
# ============================================================

sf8_annotation <- paste0(
  "Pearson r = ",
  sprintf("%.3f", sf8_pearson),
  "\nSpearman rho = ",
  sprintf("%.3f", sf8_spearman),
  "\nPC1 variance = ",
  sprintf("%.1f%%", 100 * sf8_pc1_var)
)

sf8_panel_a <- ggplot2::ggplot(
  sf8_participant,
  ggplot2::aes(
    x = .data$mean_z,
    y = .data$PC1_eigengene
  )
) +
  ggplot2::geom_point(
    size = 1.35,
    alpha = 0.42,
    color = sf8_point_col
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.65,
    color = sf8_mean_col,
    fill = "#DCE7F4",
    alpha = 0.45
  ) +
  ggplot2::annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = sf8_annotation,
    hjust = -0.08,
    vjust = 1.15,
    size = 2.45,
    color = sf8_text
  ) +
  ggplot2::labs(
    x = "Mean Hsp60/10 client z-score",
    y = "Oriented PC1 eigengene"
  ) +
  sf8_theme() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      12, 8, 8, 12
    )
  )

sf8_panel_a <- sf8_panel_label(
  sf8_panel_a,
  "A"
)

# ============================================================
# Panel B:
# stage sensitivity
# ============================================================

sf8_annotation_y <- sf8_stage_plot |>
  dplyr::group_by(.data$score_label) |>
  dplyr::summarise(
    ymax = max(
      .data$mean_score + .data$sem,
      na.rm = TRUE
    ),
    ymin = min(
      .data$mean_score - .data$sem,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    y = .data$ymax +
      0.16 * pmax(
        .data$ymax - .data$ymin,
        0.2
      )
  ) |>
  dplyr::left_join(
    sf8_mci_ad_tests |>
      dplyr::select(
    "score_label",
    "annotation"
  ),
    by = "score_label"
  )

sf8_panel_b <- ggplot2::ggplot(
  sf8_stage_plot,
  ggplot2::aes(
    x = .data$clinical_stage,
    y = .data$mean_score,
    group = 1
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "grey70",
    linewidth = 0.35
  ) +
  ggplot2::geom_line(
    linewidth = 0.7,
    color = sf8_pc1_col
  ) +
  ggplot2::geom_point(
    size = 2.2,
    color = sf8_pc1_col
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = .data$mean_score - .data$sem,
      ymax = .data$mean_score + .data$sem
    ),
    width = 0.10,
    linewidth = 0.55,
    color = sf8_pc1_col
  ) +
  ggplot2::geom_text(
    data = sf8_annotation_y,
    ggplot2::aes(
      x = "MCI",
      y = .data$y,
      label = .data$annotation
    ),
    inherit.aes = FALSE,
    size = 2.15,
    lineheight = 0.95,
    color = sf8_text
  ) +
  ggplot2::facet_wrap(
    ~score_label,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Score, mean \u00b1 SEM"
  ) +
  sf8_theme() +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(
      12, 8, 8, 10
    )
  )

sf8_panel_b <- sf8_panel_label(
  sf8_panel_b,
  "B"
)

# ============================================================
# Save panels
# ============================================================

save_panel_set(
  sf8_panel_a,
  "Supplementary_Figure_8A_PC1_concordance",
  width = 3.8,
  height = 3.5
)

save_panel_set(
  sf8_panel_b,
  "Supplementary_Figure_8B_stage_sensitivity",
  width = 4.4,
  height = 3.5
)

# ============================================================
# Assemble final
# ============================================================

sf8_final <- (
  sf8_panel_a |
    sf8_panel_b
) +
  patchwork::plot_layout(
    widths = c(0.95, 1.05)
  ) &
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      12, 10, 8, 10
    )
  )

save_plot_set(
  sf8_final,
  "Supplementary_Figure_8_PC1_sensitivity",
  width = 6.7,
  height = 3.65,
  output_dir = figures_dir
)

# ============================================================
# Audit
# ============================================================

sf8_check_pearson <- stats::cor(
  sf8_participant$mean_z,
  sf8_participant$PC1_eigengene,
  method = "pearson"
)

sf8_check_spearman <- stats::cor(
  sf8_participant$mean_z,
  sf8_participant$PC1_eigengene,
  method = "spearman"
)

sf8_audit <- tibble::tibble(
  check = c(
    "participant_file_exists",
    "stage_summary_file_exists",
    "stage_test_file_exists",
    "concordance_file_exists",
    "participant_n_577",
    "pearson_matches_frozen_summary",
    "spearman_matches_frozen_summary",
    "PC1_variance_21_8_percent"
  ),
  pass = c(
    file.exists(sf8_participant_file),
    file.exists(sf8_stage_summary_file),
    file.exists(sf8_stage_test_file),
    file.exists(sf8_concordance_file),
    sf8_n == 577,
    abs(sf8_check_pearson - sf8_pearson) < 1e-10,
    abs(sf8_check_spearman - sf8_spearman) < 1e-10,
    abs(sf8_pc1_var - 0.218252067945806) < 1e-12
  )
)

readr::write_csv(
  sf8_audit,
  file.path(
    audits_dir,
    "Supplementary_Figure_8_PC1_sensitivity_audit.csv"
  )
)

if (!all(sf8_audit$pass)) {
  print(sf8_audit, n = Inf)
  stop(
    "Supplementary Figure 8 audit failed.",
    call. = FALSE
  )
}

supfig8_outputs <- list(
  participant = sf8_participant,
  stage_summary = sf8_stage_plot,
  stage_tests = sf8_stage_tests,
  panel_a = sf8_panel_a,
  panel_b = sf8_panel_b,
  final = sf8_final,
  audit = sf8_audit
)

message("Supplementary Figure 8 complete.")
