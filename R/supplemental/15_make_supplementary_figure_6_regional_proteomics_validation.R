# ============================================================
# Supplementary Figure 6: regional proteomics validation
# ============================================================

# Purpose:
# Summarize AMP-AD regional proteomics evidence for Hsp60/10 client
# vulnerability across DLPFC and STG.

if (!exists("inputs")) {
  stop("Run 01_supplemental_load_inputs.R before 15_make_supplementary_figure_6_regional_proteomics_validation.R.", call. = FALSE)
}
if (!exists("pick_col")) {
  stop("Run 02_supplemental_helper_functions.R before 15_make_supplementary_figure_6_regional_proteomics_validation.R.", call. = FALSE)
}

# ============================================================
# Main Figure style settings
# ============================================================

sf6_hsp_col <- "#B23A48"
sf6_background_col <- "#4C78A8"
sf6_hsp_fill <- "#F4D9DE"
sf6_background_fill <- "#DCE7F4"
sf6_stg_col <- "#B23A48"
sf6_dlpfc_col <- "#1B4F9C"

sf6_group_colors <- c(
  "Non-client mitochondrial proteins" = sf6_background_col,
  "Hsp60/10 clients" = sf6_hsp_col
)

sf6_group_fills <- c(
  "Non-client mitochondrial proteins" = sf6_background_fill,
  "Hsp60/10 clients" = sf6_hsp_fill
)

sf6_region_colors <- c(
  "DLPFC" = sf6_dlpfc_col,
  "STG" = sf6_stg_col
)

sf6_theme <- function(base_size = 8.5) {
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
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey94", color = "grey80", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", color = "grey10")
    )
}

sf6_panel_label <- function(plot, label) {
  plot +
    ggplot2::labs(tag = label) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 9.5),
      plot.tag.position = c(0, 1)
    )
}

# ============================================================
# Input helpers
# ============================================================

sf6_get_input <- function(candidate_names, required = TRUE, label = "input") {
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

sf6_clean_gene <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "\\|.*$", "")
  clean_gene(x)
}

sf6_validate_tbl <- function(tbl, required_cols, label) {
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

sf6_fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("P = NA")
  if (p < 0.001) return(paste0("P = ", formatC(p, format = "e", digits = 2)))
  paste0("P = ", signif(p, 3))
}

sf6_wilcox_p <- function(df, value_col) {
  df2 <- df |>
    dplyr::filter(.data$group %in% names(sf6_group_colors)) |>
    dplyr::select("group", value = dplyr::all_of(value_col)) |>
    tidyr::drop_na()

  if (dplyr::n_distinct(df2$group) < 2) return(NA_real_)

  tryCatch(
    stats::wilcox.test(value ~ group, data = df2, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

# ============================================================
# Data preparation
# ============================================================

sf6_prepare_ad_tbl <- function(regional_protein_screen, hsp_genes, background_genes) {
  regional_protein_screen <- tibble::as_tibble(regional_protein_screen)
  
  gene_col <- dplyr::case_when(
    "gene" %in% colnames(regional_protein_screen) ~ "gene",
    "gene_symbol" %in% colnames(regional_protein_screen) ~ "gene_symbol",
    "gene_raw" %in% colnames(regional_protein_screen) ~ "gene_raw",
    TRUE ~ NA_character_
  )
  
  if (is.na(gene_col)) {
    stop(
      "regional_protein_screen is missing a gene column. Tried: gene, gene_symbol, gene_raw\n",
      "Available columns: ", paste(colnames(regional_protein_screen), collapse = ", "),
      call. = FALSE
    )
  }
  
  ad_col <- dplyr::case_when(
    "adjusted_collapse_magnitude" %in% colnames(regional_protein_screen) ~ "adjusted_collapse_magnitude",
    "ad_minus_control" %in% colnames(regional_protein_screen) ~ "ad_minus_control",
    "beta_ad" %in% colnames(regional_protein_screen) ~ "beta_ad",
    TRUE ~ NA_character_
  )
  
  if (is.na(ad_col)) {
    stop(
      "regional_protein_screen is missing an AD effect column. Tried: adjusted_collapse_magnitude, ad_minus_control, beta_ad\n",
      "Available columns: ", paste(colnames(regional_protein_screen), collapse = ", "),
      call. = FALSE
    )
  }
  
  regional_protein_screen |>
    dplyr::mutate(
      gene_symbol = sf6_clean_gene(.data[[gene_col]]),
      region_short = as.character(.data$region_short),
      ad_effect_raw = suppressWarnings(as.numeric(.data[[ad_col]])),
      ad_associated_decline = if (ad_col == "adjusted_collapse_magnitude") {
        .data$ad_effect_raw
      } else if (ad_col %in% c("ad_minus_control", "beta_ad")) {
        -.data$ad_effect_raw
      } else {
        .data$ad_effect_raw
      },
      group = dplyr::case_when(
        .data$gene_symbol %in% hsp_genes ~ "Hsp60/10 clients",
        .data$gene_symbol %in% background_genes ~ "Non-client mitochondrial proteins",
        TRUE ~ NA_character_
      ),
      region_short = factor(.data$region_short, levels = c("DLPFC", "STG")),
      group = factor(.data$group, levels = c("Non-client mitochondrial proteins", "Hsp60/10 clients"))
    ) |>
    dplyr::filter(
      .data$region_short %in% c("DLPFC", "STG"),
      !is.na(.data$group),
      is.finite(.data$ad_associated_decline)
    )
}

sf6_prepare_braak_tbl <- function(braak_models, hsp_genes, background_genes) {
  braak_models |>
    sf6_validate_tbl(
      required_cols = c("gene_symbol", "region_short", "beta_braak"),
      label = "braak_models"
    ) |>
    dplyr::mutate(
      gene_symbol = sf6_clean_gene(.data$gene_symbol),
      region_short = as.character(.data$region_short),
      beta_braak = suppressWarnings(as.numeric(.data$beta_braak)),
      inverse_braak_beta = -.data$beta_braak,
      group = dplyr::case_when(
        .data$gene_symbol %in% hsp_genes ~ "Hsp60/10 clients",
        .data$gene_symbol %in% background_genes ~ "Non-client mitochondrial proteins",
        TRUE ~ NA_character_
      ),
      region_short = factor(.data$region_short, levels = c("DLPFC", "STG")),
      group = factor(.data$group, levels = c("Non-client mitochondrial proteins", "Hsp60/10 clients"))
    ) |>
    dplyr::filter(
      .data$region_short %in% c("DLPFC", "STG"),
      !is.na(.data$group),
      is.finite(.data$inverse_braak_beta)
    )
}

sf6_prepare_difference_tbl <- function(hsp_region_compare) {
  hsp_region_compare |>
    sf6_validate_tbl(
      required_cols = c(
        "gene_symbol",
        "stg_minus_dlpfc_adjusted_collapse",
        "stg_minus_dlpfc_adjusted_inverse_braak"
      ),
      label = "hsp_region_compare"
    ) |>
    dplyr::transmute(
      gene_symbol = sf6_clean_gene(.data$gene_symbol),
      `AD-associated decline` = suppressWarnings(as.numeric(.data$stg_minus_dlpfc_adjusted_collapse)),
      `High-Braak coupling` = suppressWarnings(as.numeric(.data$stg_minus_dlpfc_adjusted_inverse_braak))
    ) |>
    tidyr::pivot_longer(
      cols = c("AD-associated decline", "High-Braak coupling"),
      names_to = "metric",
      values_to = "stg_minus_dlpfc"
    ) |>
    dplyr::mutate(
      metric = factor(.data$metric, levels = c("AD-associated decline", "High-Braak coupling"))
    ) |>
    dplyr::filter(is.finite(.data$stg_minus_dlpfc))
}

sf6_prepare_top_braak_tbl <- function(hsp_region_compare, n_top = 10) {
  hsp_region_compare |>
    sf6_validate_tbl(
      required_cols = c(
        "gene_symbol",
        "adjusted_inverse_braak_magnitude_DLPFC",
        "adjusted_inverse_braak_magnitude_STG",
        "stg_minus_dlpfc_adjusted_inverse_braak"
      ),
      label = "hsp_region_compare"
    ) |>
    dplyr::transmute(
      gene_symbol = sf6_clean_gene(.data$gene_symbol),
      DLPFC = suppressWarnings(as.numeric(.data$adjusted_inverse_braak_magnitude_DLPFC)),
      STG = suppressWarnings(as.numeric(.data$adjusted_inverse_braak_magnitude_STG)),
      stg_minus_dlpfc = suppressWarnings(as.numeric(.data$stg_minus_dlpfc_adjusted_inverse_braak))
    ) |>
    dplyr::filter(is.finite(.data$DLPFC), is.finite(.data$STG), is.finite(.data$stg_minus_dlpfc)) |>
    dplyr::arrange(dplyr::desc(.data$stg_minus_dlpfc)) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::mutate(gene_symbol = factor(.data$gene_symbol, levels = rev(.data$gene_symbol)))
}

# ============================================================
# Plot helpers
# ============================================================

sf6_distribution_plot <- function(plot_tbl, value_col, title, subtitle, ylab) {
  p_tbl <- plot_tbl |>
    dplyr::group_by(.data$region_short) |>
    dplyr::summarise(
      p_value = sf6_wilcox_p(dplyr::pick(dplyr::everything()), value_col),
      y = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = vapply(.data$p_value, sf6_fmt_p, character(1)),
      y = .data$y + 0.05 * diff(range(plot_tbl[[value_col]], na.rm = TRUE))
    )

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$group, y = .data[[value_col]], fill = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_violin(width = 0.78, alpha = 0.72, color = NA, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, color = "grey25", linewidth = 0.38, alpha = 0.95) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data$group),
      width = 0.09,
      height = 0,
      size = 0.68,
      alpha = 0.24,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = p_tbl,
      ggplot2::aes(x = 1.5, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      size = 2.35,
      fontface = "bold",
      color = "grey20"
    ) +
    ggplot2::facet_wrap(~ region_short, nrow = 1) +
    ggplot2::scale_fill_manual(values = sf6_group_fills, guide = "none") +
    ggplot2::scale_color_manual(values = sf6_group_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = c(
      "Non-client mitochondrial proteins" = "Non-client\nmito",
      "Hsp60/10 clients" = "Hsp60/10\nclients"
    )) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = ylab
    ) +
    sf6_theme(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", angle = 0, hjust = 0.5),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

sf6_difference_plot <- function(diff_tbl, paired_summary) {
  x_rng <- range(diff_tbl$stg_minus_dlpfc, na.rm = TRUE)
  ann_x <- x_rng[1] + 0.06 * diff(x_rng)

  ann_tbl <- diff_tbl |>
    dplyr::group_by(.data$metric) |>
    dplyr::summarise(
      median_diff = stats::median(.data$stg_minus_dlpfc, na.rm = TRUE),
      p_value = tryCatch(stats::wilcox.test(.data$stg_minus_dlpfc, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0(
        "median = ", signif(.data$median_diff, 3),
        "\n", vapply(.data$p_value, sf6_fmt_p, character(1))
      ),
      x = ann_x,
      y = .data$metric
    )

  ggplot2::ggplot(diff_tbl, ggplot2::aes(x = .data$stg_minus_dlpfc, y = .data$metric)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.55) +
    ggplot2::geom_jitter(height = 0.12, width = 0, alpha = 0.30, size = 1.0, color = "grey35") +
    ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA, fill = "white", color = "grey20", linewidth = 0.45, alpha = 0.85) +
    ggplot2::geom_text(
      data = ann_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.45,
      lineheight = 0.92,
      color = "grey15"
    ) +
    ggplot2::labs(
      title = "Regional difference",
      subtitle = NULL,
      x = "STG - DLPFC difference",
      y = NULL
    ) +
    sf6_theme(base_size = 8.5)
}

sf6_top_braak_plot <- function(top_tbl) {
  segment_tbl <- top_tbl |>
    dplyr::select("gene_symbol", "DLPFC", "STG")

  long_tbl <- top_tbl |>
    tidyr::pivot_longer(
      cols = c("DLPFC", "STG"),
      names_to = "region",
      values_to = "value"
    ) |>
    dplyr::mutate(region = factor(.data$region, levels = c("DLPFC", "STG")))

  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_segment(
      data = segment_tbl,
      ggplot2::aes(x = .data$DLPFC, xend = .data$STG, y = .data$gene_symbol, yend = .data$gene_symbol),
      color = "grey55",
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      data = long_tbl,
      ggplot2::aes(x = .data$value, y = .data$gene_symbol, shape = .data$region, color = .data$region),
      size = 2.7
    ) +
    ggplot2::scale_color_manual(values = sf6_region_colors, name = "Region") +
    ggplot2::scale_shape_manual(values = c("DLPFC" = 16, "STG" = 17), name = "Region") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
    ggplot2::labs(
      title = "Top STG-coupled clients",
      subtitle = NULL,
      x = "Inverse Braak magnitude",
      y = NULL
    ) +
    sf6_theme(base_size = 8.5) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(face = "bold")
    )
}

# ============================================================
# Main builder
# ============================================================

make_supfig6_regional_proteomics_validation <- function(inputs) {
  message("Starting Supplementary Figure 6: regional proteomics validation.")

  regional_input <- sf6_get_input(
    c("regional_protein_screen", "regional_protein_results", "ampad_regional_results"),
    required = TRUE,
    label = "regional protein screen"
  )

  braak_models_input <- sf6_get_input(
    c("braak_models", "regional_braak_models", "adjusted_braak_models"),
    required = TRUE,
    label = "regional adjusted Braak models"
  )

  hsp_region_input <- sf6_get_input(
    c("hsp_region_compare", "regional_hsp_compare", "hsp_dlpfc_stg_compare"),
    required = TRUE,
    label = "Hsp60/10 DLPFC/STG comparison table"
  )

  background_input <- sf6_get_input(
    c("background_null_pool", "mitochondrial_background_null_pool", "background_specificity_pool"),
    required = TRUE,
    label = "non-Hsp60/10 mitochondrial background table"
  )

  hsp_region_compare <- tibble::as_tibble(hsp_region_input$object)
  background_tbl <- tibble::as_tibble(background_input$object)

  hsp_genes <- hsp_region_compare$gene_symbol |> sf6_clean_gene() |> unique()
  background_genes <- background_tbl$gene |> sf6_clean_gene() |> setdiff(hsp_genes)

  ad_tbl <- sf6_prepare_ad_tbl(tibble::as_tibble(regional_input$object), hsp_genes, background_genes)
  braak_tbl <- sf6_prepare_braak_tbl(tibble::as_tibble(braak_models_input$object), hsp_genes, background_genes)
  diff_tbl <- sf6_prepare_difference_tbl(hsp_region_compare)
  top_braak_tbl <- sf6_prepare_top_braak_tbl(hsp_region_compare, n_top = 10)

  if (nrow(ad_tbl) == 0) {
    stop(
      "Panel A table is empty after matching regional_protein_screen genes to Hsp60/10/background genes. ",
      "Most likely cause: regional_protein_screen$gene contains values like SYMBOL|UNIPROT and needs symbol stripping.",
      call. = FALSE
    )
  }
  if (nrow(braak_tbl) == 0) {
    stop(
      "Panel B table is empty after matching braak_models genes to Hsp60/10/background genes.",
      call. = FALSE
    )
  }
  if (nrow(diff_tbl) == 0) {
    stop("Panel C difference table is empty.", call. = FALSE)
  }
  if (nrow(top_braak_tbl) == 0) {
    stop("Panel D top Braak-coupled client table is empty.", call. = FALSE)
  }

  readr::write_csv(ad_tbl, file.path(audits_dir, "SuppFig6A_ad_associated_regional_values.csv"))
  readr::write_csv(braak_tbl, file.path(audits_dir, "SuppFig6B_adjusted_braak_regional_values.csv"))
  readr::write_csv(diff_tbl, file.path(audits_dir, "SuppFig6C_stg_minus_dlpfc_values.csv"))
  readr::write_csv(top_braak_tbl, file.path(audits_dir, "SuppFig6D_top_stg_braak_coupled_clients.csv"))

  ad_summary <- ad_tbl |>
    dplyr::group_by(.data$region_short, .data$group) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_decline = stats::median(.data$ad_associated_decline, na.rm = TRUE),
      .groups = "drop"
    )

  braak_summary <- braak_tbl |>
    dplyr::group_by(.data$region_short, .data$group) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_inverse_braak_beta = stats::median(.data$inverse_braak_beta, na.rm = TRUE),
      .groups = "drop"
    )

  diff_summary <- diff_tbl |>
    dplyr::group_by(.data$metric) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_stg_minus_dlpfc = stats::median(.data$stg_minus_dlpfc, na.rm = TRUE),
      p_value = tryCatch(stats::wilcox.test(.data$stg_minus_dlpfc, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    )

  readr::write_csv(ad_summary, file.path(audits_dir, "SuppFig6A_ad_associated_regional_summary.csv"))
  readr::write_csv(braak_summary, file.path(audits_dir, "SuppFig6B_adjusted_braak_regional_summary.csv"))
  readr::write_csv(diff_summary, file.path(audits_dir, "SuppFig6C_regional_difference_summary.csv"))

  panel_a <- sf6_panel_label(
    sf6_distribution_plot(
      ad_tbl,
      value_col = "ad_associated_decline",
      title = "AD-associated decline",
      subtitle = NULL,
      ylab = "AD decline\n-(AD - Control)"
    ),
    "A"
  )

  panel_b <- sf6_panel_label(
    sf6_distribution_plot(
      braak_tbl,
      value_col = "inverse_braak_beta",
      title = "Adjusted Braak decline",
      subtitle = NULL,
      ylab = "Adjusted inverse Braak"
    ),
    "B"
  )

  panel_c <- sf6_panel_label(sf6_difference_plot(diff_tbl, NULL), "C")
  panel_d <- sf6_panel_label(sf6_top_braak_plot(top_braak_tbl), "D")

  panel_paths <- c(
    save_panel_set(panel_a, "SuppFig6A_ad_associated_protein_decline", width = 6.2, height = 4.2, output_dir = panels_dir),
    save_panel_set(panel_b, "SuppFig6B_adjusted_braak_decline", width = 6.2, height = 4.2, output_dir = panels_dir),
    save_panel_set(panel_c, "SuppFig6C_stg_minus_dlpfc_difference", width = 6.2, height = 4.2, output_dir = panels_dir),
    save_panel_set(panel_d, "SuppFig6D_top_stg_coupled_clients", width = 6.2, height = 4.2, output_dir = panels_dir)
  )

  composite <- (panel_a | panel_b) /
    (panel_c | panel_d) +
    patchwork::plot_layout(heights = c(1, 1))

  save_plot_set(
    composite,
    "Supplementary_Figure_6_regional_proteomics_validation",
    width = 12.6,
    height = 8.8,
    output_dir = figures_dir
  )

  run_summary <- c(
    "Supplementary Figure 6 run summary",
    paste0("generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "figure: regional proteomics validation",
    paste0("regional protein input: ", regional_input$name),
    paste0("Braak model input: ", braak_models_input$name),
    paste0("Hsp60/10 regional comparison input: ", hsp_region_input$name),
    paste0("background input: ", background_input$name),
    paste0("Hsp60/10 genes in regional comparison table: ", length(hsp_genes)),
    paste0("background genes available: ", length(background_genes)),
    "AD-associated decline summary:",
    paste0(
      "  - ", ad_summary$region_short, " / ", ad_summary$group,
      ": n=", ad_summary$n,
      ", median=", signif(ad_summary$median_decline, 4),
      collapse = "\n"
    ),
    "Adjusted Braak decline summary:",
    paste0(
      "  - ", braak_summary$region_short, " / ", braak_summary$group,
      ": n=", braak_summary$n,
      ", median=", signif(braak_summary$median_inverse_braak_beta, 4),
      collapse = "\n"
    ),
    "Regional difference summary:",
    paste0(
      "  - ", diff_summary$metric,
      ": n=", diff_summary$n,
      ", median STG-DLPFC=", signif(diff_summary$median_stg_minus_dlpfc, 4),
      ", ", vapply(diff_summary$p_value, sf6_fmt_p, character(1)),
      collapse = "\n"
    )
  )

  writeLines(run_summary, file.path(audits_dir, "SuppFig6_run_summary.txt"))

  message("Supplementary Figure 6 complete.")
  print(ad_summary)
  print(braak_summary)
  print(diff_summary)

  invisible(list(
    plot = composite,
    ad_tbl = ad_tbl,
    braak_tbl = braak_tbl,
    diff_tbl = diff_tbl,
    top_braak_tbl = top_braak_tbl,
    panels = list(A = panel_a, B = panel_b, C = panel_c, D = panel_d),
    panel_paths = panel_paths
  ))
}

supfig6_outputs <- make_supfig6_regional_proteomics_validation(inputs)
