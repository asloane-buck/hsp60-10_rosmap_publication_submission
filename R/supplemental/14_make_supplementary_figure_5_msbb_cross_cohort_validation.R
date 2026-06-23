# ============================================================
# Supplementary Figure 5: MSBB cross-cohort validation
# ============================================================

# Purpose:
# Test whether Hsp60/10 client pathology coupling observed in ROSMAP
# is directionally reproduced in the independent MSBB proteomics cohort.

if (!exists("inputs")) {
  stop("Run 01_supplemental_load_inputs.R before 14_make_supplementary_figure_5_msbb_cross_cohort_validation.R.", call. = FALSE)
}
if (!exists("pick_col")) {
  stop("Run 02_supplemental_helper_functions.R before 14_make_supplementary_figure_5_msbb_cross_cohort_validation.R.", call. = FALSE)
}

# ============================================================
# Main Figure style settings
# ============================================================

sf5_hsp_col <- "#B23A48"
sf5_background_col <- "#4C78A8"
sf5_hsp_fill <- "#F4D9DE"
sf5_background_fill <- "#DCE7F4"
sf5_text <- "#111827"

sf5_group_colors <- c(
  "Non-client mitochondrial proteins" = sf5_background_col,
  "Hsp60/10 clients" = sf5_hsp_col
)

sf5_group_fills <- c(
  "Non-client mitochondrial proteins" = sf5_background_fill,
  "Hsp60/10 clients" = sf5_hsp_fill
)

sf5_concordance_colors <- c(
  "Concordant inverse\nBraak effect" = sf5_hsp_col,
  "Not concordant" = "grey70"
)

sf5_theme <- function(base_size = 8.5) {
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
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

sf5_panel_label <- function(plot, label, tag_position = c(0, 1), plot_margin = NULL) {
  tag_theme <- ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 9.5),
    plot.tag.position = tag_position
  )
  if (!is.null(plot_margin)) {
    tag_theme <- tag_theme + ggplot2::theme(plot.margin = plot_margin)
  }

  plot +
    ggplot2::labs(tag = label) +
    tag_theme
}

# ============================================================
# Input helpers
# ============================================================

sf5_get_input <- function(candidate_names, required = TRUE, label = "input") {
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

sf5_clean_gene <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "\\|.*$", "")
  clean_gene(x)
}

sf5_validate_tbl <- function(tbl, required_cols, label) {
  missing_cols <- setdiff(required_cols, colnames(tbl))
  if (length(missing_cols) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(tbl), collapse = ", "),
      call. = FALSE
    )
  }
  tibble::as_tibble(tbl)
}

sf5_fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("P = NA")
  if (p < 0.001) return(paste0("P = ", formatC(p, format = "e", digits = 2)))
  paste0("P = ", signif(p, 3))
}

sf5_wilcox_p <- function(df, value_col) {
  df2 <- df |>
    dplyr::filter(.data$group %in% names(sf5_group_colors)) |>
    dplyr::select(.data$group, value = dplyr::all_of(value_col)) |>
    tidyr::drop_na()

  if (dplyr::n_distinct(df2$group) < 2) return(NA_real_)

  tryCatch(
    stats::wilcox.test(value ~ group, data = df2, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

sf5_wilcox_directional_p <- function(df, value_col) {
  df2 <- df |>
    dplyr::filter(.data$group %in% names(sf5_group_colors)) |>
    dplyr::transmute(
      group = factor(
        as.character(.data$group),
        levels = c(
          "Hsp60/10 clients",
          "Non-client mitochondrial proteins"
        )
      ),
      value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) |>
    tidyr::drop_na()

  if (dplyr::n_distinct(df2$group) < 2) return(NA_real_)

  tryCatch(
    stats::wilcox.test(
      value ~ group,
      data = df2,
      alternative = "greater",
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

sf5_spearman_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3) {
    return(tibble::tibble(n = sum(ok), rho = NA_real_, p_value = NA_real_))
  }

  out <- tryCatch(
    stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(tibble::tibble(n = sum(ok), rho = NA_real_, p_value = NA_real_))
  }

  tibble::tibble(
    n = sum(ok),
    rho = unname(out$estimate),
    p_value = out$p.value
  )
}

# ============================================================
# MSBB input loading / standardization
# ============================================================

sf5_find_first_existing <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

sf5_load_msbb_combined <- function() {
  candidate_paths <- if (exists("msbb_validation_paths")) {
    msbb_validation_paths
  } else {
    c(
      file.path(supplemental_input_dir, "MSBB_covariate_adjusted_braak_and_collapse.csv"),
      file.path(supplemental_input_dir, "MSBB_standardized_braak_collapse_percentiles.csv"),
      file.path(supplemental_input_dir, "MSBB_covariate_adjusted_braak_collapse_percentiles.csv"),
      file.path(supplemental_input_dir, "COMBINED_ROSMAP_MSBB_validation_results.csv")
    )
  }
  
  selected_path <- sf5_find_first_existing(candidate_paths)
  
  if (is.na(selected_path)) {
    stop(
      "Could not find a saved MSBB validation results file.\n",
      "Tried:\n",
      paste(candidate_paths, collapse = "\n"),
      "\nPlace a MSBB validation table in data/supplemental_inputs/ or set MSBB_VALIDATION_FILE.",
      call. = FALSE
    )
  }
  
  message("Loading MSBB combined validation file: ", selected_path)
  
  out <- readr::read_csv(selected_path, show_col_types = FALSE)
  
  # If this is a combined ROSMAP/MSBB table, keep only MSBB rows.
  if ("cohort" %in% colnames(out)) {
    out <- out |>
      dplyr::filter(toupper(as.character(.data$cohort)) == "MSBB")
  }
  
  out
}

sf5_standardize_msbb_inputs <- function(msbb_combined) {
  msbb_combined <- tibble::as_tibble(msbb_combined)
  
  gene_col <- dplyr::case_when(
    "gene" %in% colnames(msbb_combined) ~ "gene",
    "gene_symbol" %in% colnames(msbb_combined) ~ "gene_symbol",
    "Gene" %in% colnames(msbb_combined) ~ "Gene",
    TRUE ~ NA_character_
  )
  
  if (is.na(gene_col)) {
    stop(
      "MSBB combined table lacks a gene column. Available columns: ",
      paste(colnames(msbb_combined), collapse = ", "),
      call. = FALSE
    )
  }
  
  braak_beta_col <- dplyr::case_when(
    "braak_beta" %in% colnames(msbb_combined) ~ "braak_beta",
    "beta_braak" %in% colnames(msbb_combined) ~ "beta_braak",
    TRUE ~ NA_character_
  )
  
  braak_p_col <- dplyr::case_when(
    "braak_p" %in% colnames(msbb_combined) ~ "braak_p",
    "p_braak" %in% colnames(msbb_combined) ~ "p_braak",
    TRUE ~ NA_character_
  )
  
  braak_n_col <- dplyr::case_when(
    "braak_n" %in% colnames(msbb_combined) ~ "braak_n",
    "n_braak" %in% colnames(msbb_combined) ~ "n_braak",
    "n_model" %in% colnames(msbb_combined) ~ "n_model",
    TRUE ~ NA_character_
  )
  
  braak_rho_col <- dplyr::case_when(
    "braak_rho_raw" %in% colnames(msbb_combined) ~ "braak_rho_raw",
    "rho_braak" %in% colnames(msbb_combined) ~ "rho_braak",
    "braak_rho" %in% colnames(msbb_combined) ~ "braak_rho",
    TRUE ~ NA_character_
  )
  
  high_minus_low_col <- dplyr::case_when(
    "high_minus_low" %in% colnames(msbb_combined) ~ "high_minus_low",
    "late_shift_high_minus_low" %in% colnames(msbb_combined) ~ "late_shift_high_minus_low",
    TRUE ~ NA_character_
  )
  
  collapse_col <- dplyr::case_when(
    "collapse_magnitude" %in% colnames(msbb_combined) ~ "collapse_magnitude",
    "collapse_magnitude_pos" %in% colnames(msbb_combined) ~ "collapse_magnitude_pos",
    "collapse_magnitude_positive" %in% colnames(msbb_combined) ~ "collapse_magnitude_positive",
    TRUE ~ NA_character_
  )
  
  missing_core <- c(
    if (is.na(braak_beta_col)) "braak_beta/beta_braak",
    if (is.na(braak_p_col)) "braak_p/p_braak",
    if (is.na(braak_n_col)) "braak_n/n_braak/n_model",
    if (is.na(collapse_col)) "collapse_magnitude/collapse_magnitude_pos"
  )
  
  if (length(missing_core) > 0) {
    stop(
      "MSBB combined table is missing required conceptual columns: ",
      paste(missing_core, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(msbb_combined), collapse = ", "),
      call. = FALSE
    )
  }
  
  msbb_braak_effects <- msbb_combined |>
    dplyr::transmute(
      gene = sf5_clean_gene(.data[[gene_col]]),
      braak_beta = suppressWarnings(as.numeric(.data[[braak_beta_col]])),
      braak_p = suppressWarnings(as.numeric(.data[[braak_p_col]])),
      braak_n = suppressWarnings(as.numeric(.data[[braak_n_col]])),
      braak_rho_raw = if (!is.na(braak_rho_col)) {
        suppressWarnings(as.numeric(.data[[braak_rho_col]]))
      } else {
        NA_real_
      }
    ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  
  msbb_collapse <- msbb_combined |>
    dplyr::transmute(
      gene = sf5_clean_gene(.data[[gene_col]]),
      high_minus_low = if (!is.na(high_minus_low_col)) {
        suppressWarnings(as.numeric(.data[[high_minus_low_col]]))
      } else {
        -suppressWarnings(as.numeric(.data[[collapse_col]]))
      },
      collapse_magnitude = suppressWarnings(as.numeric(.data[[collapse_col]])),
      collapse_magnitude_positive = suppressWarnings(as.numeric(.data[[collapse_col]])),
      collapse_n_low = if ("collapse_n_low" %in% colnames(msbb_combined)) {
        suppressWarnings(as.numeric(.data$collapse_n_low))
      } else {
        NA_real_
      },
      collapse_n_high = if ("collapse_n_high" %in% colnames(msbb_combined)) {
        suppressWarnings(as.numeric(.data$collapse_n_high))
      } else {
        NA_real_
      }
    ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  
  list(
    msbb_braak_effects = msbb_braak_effects,
    msbb_collapse = msbb_collapse
  )
}

# ============================================================
# Data preparation
# ============================================================

sf5_prepare_rosmap_tbl <- function(hsp_tbl, background_tbl) {
  hsp_tbl <- hsp_tbl |>
    sf5_validate_tbl(
      required_cols = c("gene", "inverse_braak_magnitude", "protein_collapse_magnitude"),
      label = "hsp_null_tbl"
    ) |>
    dplyr::mutate(
      gene = sf5_clean_gene(.data$gene),
      group = "Hsp60/10 clients"
    )

  background_tbl <- background_tbl |>
    sf5_validate_tbl(
      required_cols = c("gene", "inverse_braak_magnitude", "protein_collapse_magnitude"),
      label = "background_null_pool"
    ) |>
    dplyr::mutate(
      gene = sf5_clean_gene(.data$gene),
      group = "Non-client mitochondrial proteins"
    ) |>
    dplyr::filter(!.data$gene %in% hsp_tbl$gene)

  dplyr::bind_rows(hsp_tbl, background_tbl) |>
    dplyr::transmute(
      gene = .data$gene,
      group = factor(.data$group, levels = names(sf5_group_colors)),
      rosmap_inverse_braak = suppressWarnings(as.numeric(.data$inverse_braak_magnitude)),
      rosmap_collapse = suppressWarnings(as.numeric(.data$protein_collapse_magnitude))
    ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene))
}

sf5_prepare_msbb_tbl <- function(msbb_braak_effects, msbb_collapse, hsp_genes, background_genes) {
  braak_tbl <- msbb_braak_effects |>
    sf5_validate_tbl(
      required_cols = c("gene", "braak_beta", "braak_p", "braak_n", "braak_rho_raw"),
      label = "msbb_braak_effects"
    ) |>
    dplyr::transmute(
      gene = sf5_clean_gene(.data$gene),
      msbb_braak_beta = suppressWarnings(as.numeric(.data$braak_beta)),
      msbb_inverse_braak_beta = -.data$msbb_braak_beta,
      msbb_braak_p = suppressWarnings(as.numeric(.data$braak_p)),
      msbb_braak_n = suppressWarnings(as.numeric(.data$braak_n)),
      msbb_braak_rho_raw = suppressWarnings(as.numeric(.data$braak_rho_raw)),
      msbb_inverse_braak_rho = -.data$msbb_braak_rho_raw
    )

  collapse_tbl <- msbb_collapse |>
    sf5_validate_tbl(
      required_cols = c("gene", "high_minus_low", "collapse_magnitude"),
      label = "msbb_collapse"
    ) |>
    dplyr::mutate(gene = sf5_clean_gene(.data$gene)) |>
    dplyr::transmute(
      gene = .data$gene,
      msbb_high_minus_low = suppressWarnings(as.numeric(.data$high_minus_low)),
      msbb_collapse = dplyr::case_when(
        "collapse_magnitude_positive" %in% colnames(msbb_collapse) ~ suppressWarnings(as.numeric(.data$collapse_magnitude_positive)),
        TRUE ~ -suppressWarnings(as.numeric(.data$high_minus_low))
      ),
      msbb_collapse_n_low = if ("collapse_n_low" %in% colnames(msbb_collapse)) suppressWarnings(as.numeric(.data$collapse_n_low)) else NA_real_,
      msbb_collapse_n_high = if ("collapse_n_high" %in% colnames(msbb_collapse)) suppressWarnings(as.numeric(.data$collapse_n_high)) else NA_real_
    )

  braak_tbl |>
    dplyr::left_join(collapse_tbl, by = "gene") |>
    dplyr::mutate(
      group = dplyr::case_when(
        .data$gene %in% hsp_genes ~ "Hsp60/10 clients",
        .data$gene %in% background_genes ~ "Non-client mitochondrial proteins",
        TRUE ~ NA_character_
      ),
      group = factor(.data$group, levels = names(sf5_group_colors))
    ) |>
    dplyr::filter(!is.na(.data$group))
}

sf5_prepare_cross_tbl <- function(rosmap_tbl, msbb_tbl) {
  rosmap_tbl |>
    dplyr::select(.data$gene, .data$group, .data$rosmap_inverse_braak, .data$rosmap_collapse) |>
    dplyr::inner_join(
      msbb_tbl |>
        dplyr::select(
          .data$gene,
          .data$msbb_inverse_braak_beta,
          .data$msbb_inverse_braak_rho,
          .data$msbb_braak_p,
          .data$msbb_braak_n,
          .data$msbb_collapse
        ),
      by = "gene"
    ) |>
    dplyr::mutate(
      concordant_inverse_braak = is.finite(.data$rosmap_inverse_braak) &
        is.finite(.data$msbb_inverse_braak_beta) &
        .data$rosmap_inverse_braak > 0 &
        .data$msbb_inverse_braak_beta > 0,
      concordance_class = dplyr::if_else(
        .data$concordant_inverse_braak,
        "Concordant inverse\nBraak effect",
        "Not concordant"
      ),
      concordance_class = factor(.data$concordance_class, levels = names(sf5_concordance_colors))
    )
}

# ============================================================
# Plot helpers
# ============================================================

sf5_distribution_plot <- function(plot_tbl, value_col, title, subtitle, ylab) {
  plot_tbl <- plot_tbl |>
    dplyr::filter(is.finite(.data[[value_col]]))

  p_two <- sf5_wilcox_p(plot_tbl, value_col)

  p_directional <- if (value_col == "msbb_inverse_braak_beta") {
    sf5_wilcox_directional_p(plot_tbl, value_col)
  } else {
    NA_real_
  }

  label_txt <- if (value_col == "msbb_inverse_braak_beta") {
    paste0(
      "Directional Wilcoxon ", sf5_fmt_p(p_directional),
      "\nTwo-sided ", sf5_fmt_p(p_two)
    )
  } else {
    paste0("Wilcoxon ", sf5_fmt_p(p_two))
  }

  y_max <- max(plot_tbl[[value_col]], na.rm = TRUE)
  y_min <- min(plot_tbl[[value_col]], na.rm = TRUE)
  y_rng <- y_max - y_min
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- 0.1

  annotation_y <- y_max + if (value_col == "msbb_inverse_braak_beta") 0.13 * y_rng else 0.09 * y_rng
  upper_pad <- if (value_col == "msbb_inverse_braak_beta") 0.24 * y_rng else 0.18 * y_rng

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$group, y = .data[[value_col]], fill = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_violin(width = 0.78, alpha = 0.72, color = NA, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, color = "grey25", linewidth = 0.38, alpha = 0.95) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data$group),
      width = 0.09,
      height = 0,
      size = 0.68,
      alpha = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = annotation_y,
      label = label_txt,
      size = 3.15,
      fontface = "bold",
      color = "grey20",
      lineheight = 0.95
    ) +
    ggplot2::scale_fill_manual(values = sf5_group_fills, guide = "none") +
    ggplot2::scale_color_manual(values = sf5_group_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = c(
      "Non-client mitochondrial proteins" = "Non-client\nmitochondrial proteins",
      "Hsp60/10 clients" = "Hsp60/10\nclients"
    )) +
    ggplot2::coord_cartesian(ylim = c(y_min - 0.04 * y_rng, y_max + upper_pad), clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = ylab
    ) +
    sf5_theme(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold"),
      legend.position = "none",
      plot.margin = ggplot2::margin(5, 8, 5, 18)
    )
}

sf5_cross_cohort_scatter <- function(cross_tbl) {
  plot_tbl <- cross_tbl |>
    dplyr::filter(.data$group == "Hsp60/10 clients") |>
    dplyr::filter(is.finite(.data$rosmap_inverse_braak), is.finite(.data$msbb_inverse_braak_beta))

  cor_tbl <- sf5_spearman_test(plot_tbl$rosmap_inverse_braak, plot_tbl$msbb_inverse_braak_beta)
  rho <- cor_tbl$rho[1]
  p_cor <- cor_tbl$p_value[1]

  priority_label_genes <- c("NDUFA10", "PDHA1", "NFS1", "TIMM13", "SDHB", "IDH2", "NIPSNAP1", "PDHB")
  label_tbl <- plot_tbl |>
    dplyr::filter(.data$concordant_inverse_braak, .data$gene %in% priority_label_genes) |>
    dplyr::mutate(
      gene = factor(.data$gene, levels = priority_label_genes),
      nudge_x = dplyr::case_when(
        .data$gene == "NDUFA10" ~ -0.060,
        .data$gene == "NIPSNAP1" ~ -0.050,
        .data$gene == "TIMM13" ~ -0.030,
        .data$gene %in% c("IDH2", "PDHB", "SDHB") ~ 0.050,
        TRUE ~ 0.018
      ),
      nudge_y = dplyr::case_when(
        .data$gene == "NDUFA10" ~ 0.050,
        .data$gene == "SDHB" ~ 0.040,
        .data$gene == "TIMM13" ~ -0.024,
        .data$gene == "IDH2" ~ -0.006,
        .data$gene == "PDHB" ~ 0.050,
        .data$gene == "PDHA1" ~ 0.036,
        .data$gene == "NFS1" ~ -0.002,
        .data$gene == "NIPSNAP1" ~ -0.046,
        TRUE ~ 0
      )
    ) |>
    dplyr::arrange(.data$gene)

  ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = .data$rosmap_inverse_braak, y = .data$msbb_inverse_braak_beta)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$concordance_class),
      size = 1.8,
      alpha = 0.80
    ) +
    ggrepel::geom_text_repel(
      data = label_tbl,
      ggplot2::aes(label = .data$gene),
      size = 2.55,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.28,
      point.padding = 0.18,
      force = 3.0,
      force_pull = 0.25,
      seed = 730,
      nudge_x = label_tbl$nudge_x,
      nudge_y = label_tbl$nudge_y,
      segment.size = 0.22,
      color = "grey15"
    ) +
    ggplot2::scale_color_manual(values = sf5_concordance_colors, name = NULL, drop = FALSE) +
    ggplot2::annotate(
      "label",
      x = Inf,
      y = -Inf,
      hjust = 1.03,
      vjust = -0.35,
      label = paste0("Spearman rho = ", signif(rho, 3), "\n", sf5_fmt_p(p_cor)),
      size = 2.85,
      label.size = 0.25,
      fill = "white",
      color = "grey20"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.12))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.16))) +
    ggplot2::labs(
      title = "Cross-cohort comparison",
      subtitle = NULL,
      x = "ROSMAP inverse Braak magnitude",
      y = "MSBB inverse Braak beta"
    ) +
    sf5_theme(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "top",
      plot.margin = ggplot2::margin(5, 12, 5, 10)
    )
}

sf5_concordance_summary_plot <- function(cross_tbl) {
  hsp_tbl <- cross_tbl |>
    dplyr::filter(.data$group == "Hsp60/10 clients") |>
    dplyr::filter(is.finite(.data$rosmap_inverse_braak), is.finite(.data$msbb_inverse_braak_beta)) |>
    dplyr::mutate(
      category = dplyr::case_when(
        .data$rosmap_inverse_braak > 0 & .data$msbb_inverse_braak_beta > 0 ~ "Inverse in both cohorts",
        .data$rosmap_inverse_braak > 0 & .data$msbb_inverse_braak_beta <= 0 ~ "ROSMAP only",
        .data$rosmap_inverse_braak <= 0 & .data$msbb_inverse_braak_beta > 0 ~ "MSBB only",
        TRUE ~ "Neither / opposite"
      ),
      category = factor(
        .data$category,
        levels = c("Inverse in both cohorts", "ROSMAP only", "MSBB only", "Neither / opposite")
      )
    )

  count_tbl <- hsp_tbl |>
    dplyr::count(.data$category, name = "n") |>
    tidyr::complete(.data$category, fill = list(n = 0)) |>
    dplyr::mutate(
      fraction = .data$n / sum(.data$n),
      label = paste0(.data$n, " / ", sum(.data$n), "\n", scales::percent(.data$fraction, accuracy = 1))
    )

  ggplot2::ggplot(count_tbl, ggplot2::aes(x = .data$fraction, y = .data$category)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$category == "Inverse in both cohorts"),
      width = 0.58,
      color = "white",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      hjust = -0.04,
      size = 3.0,
      color = "grey15"
    ) +
    ggplot2::scale_fill_manual(values = c("TRUE" = sf5_hsp_col, "FALSE" = "grey72"), guide = "none") +
    ggplot2::scale_x_continuous(labels = scales::percent, expand = ggplot2::expansion(mult = c(0, 0.52))) +
    ggplot2::labs(
      title = "Gene-level concordance",
      subtitle = NULL,
      x = "Fraction of Hsp60/10 clients",
      y = NULL
    ) +
    sf5_theme(base_size = 8.5)
}

# ============================================================
# Main builder
# ============================================================

make_supfig5_msbb_cross_cohort_validation <- function(inputs) {
  message("Starting Supplementary Figure 5: MSBB cross-cohort validation.")

  hsp_input <- sf5_get_input(
    c("hsp_null_tbl", "hsp60_10_null_tbl", "hsp_specificity_tbl"),
    required = TRUE,
    label = "ROSMAP Hsp60/10 table"
  )
  background_input <- sf5_get_input(
    c("background_null_pool", "mitochondrial_background_null_pool", "background_specificity_pool"),
    required = TRUE,
    label = "ROSMAP non-Hsp60/10 mitochondrial background table"
  )
  msbb_braak_input <- sf5_get_input(
    c("msbb_braak_effects", "msbb_braak_results", "msbb_gene_braak_results"),
    required = FALSE,
    label = "MSBB Braak effects table"
  )
  
  msbb_collapse_input <- sf5_get_input(
    c("msbb_collapse", "msbb_collapse_results", "msbb_gene_collapse_results"),
    required = FALSE,
    label = "MSBB collapse table"
  )
  
  if (is.null(msbb_braak_input$object) || is.null(msbb_collapse_input$object)) {
    message("MSBB split objects not found in inputs/global environment; loading saved MSBB combined validation file.")
    msbb_combined <- sf5_load_msbb_combined()
    msbb_standardized <- sf5_standardize_msbb_inputs(msbb_combined)
    
    msbb_braak_input <- list(
      name = "MSBB combined file -> standardized Braak effects",
      object = msbb_standardized$msbb_braak_effects
    )
    
    msbb_collapse_input <- list(
      name = "MSBB combined file -> standardized collapse",
      object = msbb_standardized$msbb_collapse
    )
    
    inputs$msbb_braak_effects <- msbb_braak_input$object
    inputs$msbb_collapse <- msbb_collapse_input$object
  }


  rosmap_tbl <- sf5_prepare_rosmap_tbl(
    tibble::as_tibble(hsp_input$object),
    tibble::as_tibble(background_input$object)
  )

  hsp_genes <- rosmap_tbl |>
    dplyr::filter(.data$group == "Hsp60/10 clients") |>
    dplyr::pull(.data$gene) |>
    unique()

  background_genes <- rosmap_tbl |>
    dplyr::filter(.data$group == "Non-client mitochondrial proteins") |>
    dplyr::pull(.data$gene) |>
    unique()

  msbb_tbl <- sf5_prepare_msbb_tbl(
    tibble::as_tibble(msbb_braak_input$object),
    tibble::as_tibble(msbb_collapse_input$object),
    hsp_genes = hsp_genes,
    background_genes = background_genes
  )

  cross_tbl <- sf5_prepare_cross_tbl(rosmap_tbl, msbb_tbl)

  if (nrow(msbb_tbl) == 0) {
    stop("MSBB validation table is empty after matching to Hsp60/10/background genes.", call. = FALSE)
  }
  if (sum(cross_tbl$group == "Hsp60/10 clients", na.rm = TRUE) < 10) {
    stop("Too few Hsp60/10 clients overlap between ROSMAP and MSBB for cross-cohort comparison.", call. = FALSE)
  }

  readr::write_csv(rosmap_tbl, file.path(audits_dir, "SuppFig5_rosmap_discovery_values.csv"))
  readr::write_csv(msbb_tbl, file.path(audits_dir, "SuppFig5_msbb_validation_values.csv"))
  readr::write_csv(cross_tbl, file.path(audits_dir, "SuppFig5_cross_cohort_gene_values.csv"))

  rosmap_summary <- rosmap_tbl |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_inverse_braak = stats::median(.data$rosmap_inverse_braak, na.rm = TRUE),
      median_collapse = stats::median(.data$rosmap_collapse, na.rm = TRUE),
      .groups = "drop"
    )

  msbb_summary <- msbb_tbl |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_inverse_braak_beta = stats::median(.data$msbb_inverse_braak_beta, na.rm = TRUE),
      median_inverse_braak_rho = stats::median(.data$msbb_inverse_braak_rho, na.rm = TRUE),
      median_collapse = stats::median(.data$msbb_collapse, na.rm = TRUE),
      .groups = "drop"
    )

  hsp_cross_tbl <- cross_tbl |>
    dplyr::filter(.data$group == "Hsp60/10 clients") |>
    dplyr::filter(is.finite(.data$rosmap_inverse_braak), is.finite(.data$msbb_inverse_braak_beta))

  hsp_cross_cor <- sf5_spearman_test(
    hsp_cross_tbl$rosmap_inverse_braak,
    hsp_cross_tbl$msbb_inverse_braak_beta
  )

  concordance_summary <- hsp_cross_tbl |>
    dplyr::summarise(
      n_hsp_overlap = dplyr::n(),
      n_rosmap_inverse = sum(.data$rosmap_inverse_braak > 0, na.rm = TRUE),
      n_msbb_inverse_beta = sum(.data$msbb_inverse_braak_beta > 0, na.rm = TRUE),
      n_concordant_inverse_braak = sum(.data$concordant_inverse_braak, na.rm = TRUE),
      pct_concordant_inverse_braak = 100 * mean(.data$concordant_inverse_braak, na.rm = TRUE),
      expected_concordant_by_margins = (
        mean(.data$rosmap_inverse_braak > 0, na.rm = TRUE) *
          mean(.data$msbb_inverse_braak_beta > 0, na.rm = TRUE) *
          dplyr::n()
      ),
      spearman_rho = hsp_cross_cor$rho[1],
      spearman_p = hsp_cross_cor$p_value[1],
      .groups = "drop"
    )

  validation_test_summary <- tibble::tibble(
    test = c(
      "MSBB Hsp60/10 vs mitochondrial background, inverse Braak beta, two-sided Wilcoxon",
      "MSBB Hsp60/10 greater than mitochondrial background, inverse Braak beta, directional Wilcoxon",
      "ROSMAP-MSBB Hsp60/10 gene-level Spearman correlation, inverse Braak beta"
    ),
    statistic = c(
      NA_real_,
      NA_real_,
      concordance_summary$spearman_rho[1]
    ),
    p_value = c(
      sf5_wilcox_p(msbb_tbl, "msbb_inverse_braak_beta"),
      sf5_wilcox_directional_p(msbb_tbl, "msbb_inverse_braak_beta"),
      concordance_summary$spearman_p[1]
    ),
    interpretation = c(
      "Borderline two-sided group shift",
      "Pre-specified directional validation test",
      "Modest gene-level cross-cohort concordance"
    )
  )

  readr::write_csv(rosmap_summary, file.path(audits_dir, "SuppFig5_rosmap_group_summary.csv"))
  readr::write_csv(msbb_summary, file.path(audits_dir, "SuppFig5_msbb_group_summary.csv"))
  readr::write_csv(concordance_summary, file.path(audits_dir, "SuppFig5_concordance_summary.csv"))
  readr::write_csv(validation_test_summary, file.path(audits_dir, "SuppFig5_validation_test_summary.csv"))

  panel_a <- sf5_panel_label(
    sf5_distribution_plot(
      rosmap_tbl,
      value_col = "rosmap_inverse_braak",
      title = "ROSMAP discovery distribution",
      subtitle = NULL,
      ylab = "ROSMAP inverse Braak\nmagnitude"
    ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(size = 7.8, margin = ggplot2::margin(r = 4)),
        plot.margin = ggplot2::margin(15, 8, 5, 30)
      ),
    "A"
  )

  panel_b <- sf5_panel_label(
    sf5_distribution_plot(
      msbb_tbl,
      value_col = "msbb_inverse_braak_beta",
      title = "MSBB validation distribution",
      subtitle = NULL,
      ylab = "MSBB inverse Braak beta"
    ),
    "B",
    tag_position = c(-0.08, 1.06),
    plot_margin = ggplot2::margin(15, 8, 5, 28)
  )

  panel_c <- sf5_panel_label(sf5_cross_cohort_scatter(cross_tbl), "C")
  panel_d <- sf5_panel_label(sf5_concordance_summary_plot(cross_tbl), "D")

  panel_paths <- c(
    save_panel_set(panel_a, "SuppFig5A_rosmap_discovery_distribution", width = 5.9, height = 4.0, output_dir = panels_dir),
    save_panel_set(panel_b, "SuppFig5B_msbb_validation_distribution", width = 5.9, height = 4.0, output_dir = panels_dir),
    save_panel_set(panel_c, "SuppFig5C_cross_cohort_effect_comparison", width = 5.9, height = 4.0, output_dir = panels_dir),
    save_panel_set(panel_d, "SuppFig5D_gene_level_concordance", width = 5.9, height = 4.0, output_dir = panels_dir)
  )

  composite <- (panel_a | panel_b) /
    (panel_c | panel_d) +
    patchwork::plot_layout(heights = c(1, 1))

  save_plot_set(
    composite,
    "Supplementary_Figure_5_msbb_cross_cohort_validation",
    width = 12.4,
    height = 8.6,
    output_dir = figures_dir
  )

  run_summary <- c(
    "Supplementary Figure 5 run summary",
    paste0("generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "figure: external MSBB cross-cohort validation",
    paste0("ROSMAP Hsp60/10 input: ", hsp_input$name),
    paste0("ROSMAP background input: ", background_input$name),
    paste0("MSBB Braak input: ", msbb_braak_input$name),
    paste0("MSBB collapse input: ", msbb_collapse_input$name),
    paste0("ROSMAP genes: ", nrow(rosmap_tbl)),
    paste0("MSBB matched genes: ", nrow(msbb_tbl)),
    paste0("Cross-cohort genes: ", nrow(cross_tbl)),
    "ROSMAP summary:",
    paste0(
      "  - ", rosmap_summary$group,
      ": n=", rosmap_summary$n,
      ", median inverse Braak=", signif(rosmap_summary$median_inverse_braak, 4),
      ", median collapse=", signif(rosmap_summary$median_collapse, 4),
      collapse = "
"
    ),
    "MSBB summary:",
    paste0(
      "  - ", msbb_summary$group,
      ": n=", msbb_summary$n,
      ", median inverse Braak beta=", signif(msbb_summary$median_inverse_braak_beta, 4),
      ", median inverse Braak rho=", signif(msbb_summary$median_inverse_braak_rho, 4),
      ", median collapse=", signif(msbb_summary$median_collapse, 4),
      collapse = "
"
    ),
    paste0(
      "Hsp60/10 concordance: ",
      concordance_summary$n_concordant_inverse_braak,
      " / ",
      concordance_summary$n_hsp_overlap,
      " (",
      signif(concordance_summary$pct_concordant_inverse_braak, 3),
      "%) concordant inverse Braak; expected by marginal positivity=",
      signif(concordance_summary$expected_concordant_by_margins, 4),
      "; Spearman rho=",
      signif(concordance_summary$spearman_rho, 4),
      "; Spearman P=",
      signif(concordance_summary$spearman_p, 4)
    ),
    "Validation test summary:",
    paste0(
      "  - ", validation_test_summary$test,
      ": P=", signif(validation_test_summary$p_value, 4),
      collapse = "
"
    )
  )

  writeLines(run_summary, file.path(audits_dir, "SuppFig5_run_summary.txt"))

  message("Supplementary Figure 5 complete.")
  print(rosmap_summary)
  print(msbb_summary)
  print(concordance_summary)
  print(validation_test_summary)

  invisible(list(
    plot = composite,
    rosmap_tbl = rosmap_tbl,
    msbb_tbl = msbb_tbl,
    cross_tbl = cross_tbl,
    rosmap_summary = rosmap_summary,
    msbb_summary = msbb_summary,
    concordance_summary = concordance_summary,
    validation_test_summary = validation_test_summary,
    panels = list(A = panel_a, B = panel_b, C = panel_c, D = panel_d),
    panel_paths = panel_paths
  ))
}

supfig5_outputs <- make_supfig5_msbb_cross_cohort_validation(inputs)
