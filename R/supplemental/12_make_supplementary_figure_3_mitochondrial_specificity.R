# ============================================================
# Supplementary Figure 3: specificity against mitochondrial background
# ============================================================

# Purpose:
# Test whether Hsp60/10 clients show stronger collapse, inverse Braak coupling,
# cognition association, and AGORA support than abundance-matched
# non-Hsp60/10 mitochondrial background proteins.

if (!exists("inputs")) {
  stop("Run 01_supplemental_load_inputs.R before 12_make_supplementary_figure_3_mitochondrial_specificity.R.", call. = FALSE)
}
if (!exists("pick_col")) {
  stop("Run 02_supplemental_helper_functions.R before 12_make_supplementary_figure_3_mitochondrial_specificity.R.", call. = FALSE)
}

# ============================================================
# Main Figure style settings
# ============================================================

sf3_hsp_col <- "#B23A48"
sf3_background_col <- "#4C78A8"
sf3_null_fill <- "#D8DEE6"
sf3_null_line <- "grey35"
sf3_text <- "#111827"

sf3_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey25", linewidth = 0.45),
      axis.line = ggplot2::element_line(color = "grey25", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "grey25", linewidth = 0.35),
      axis.text = ggplot2::element_text(color = "grey15"),
      axis.title = ggplot2::element_text(color = "grey10", face = "bold"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

sf3_panel_label <- function(plot, label) {
  plot +
    ggplot2::labs(tag = label) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 9.5),
      plot.tag.position = c(0, 1)
    )
}

sf3_get_input <- function(candidate_names, required = TRUE, label = "input") {
  selected_name <- candidate_names[candidate_names %in% names(inputs)][1]

  if (!is.na(selected_name)) {
    message("Detected ", label, ": ", selected_name)
    return(list(name = selected_name, object = inputs[[selected_name]]))
  }

  # Fallback only helps if the object is already in the session.
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

sf3_validate_tbl <- function(tbl, label) {
  required_cols <- c(
    "gene",
    "protein_collapse_magnitude",
    "inverse_braak_magnitude",
    "cognition_composite_score",
    "agora_nominated_target",
    "abundance_bin"
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

  tbl |>
    tibble::as_tibble() |>
    dplyr::mutate(
      gene = clean_gene(.data$gene),
      abundance_bin = as.character(.data$abundance_bin),
      agora_nominated_target = dplyr::case_when(
        is.logical(.data$agora_nominated_target) ~ as.numeric(.data$agora_nominated_target),
        TRUE ~ suppressWarnings(as.numeric(.data$agora_nominated_target))
      )
    ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene))
}

############################################################
## Stand-alone AGORA target loading / annotation
############################################################

sf3_agora_path <- Sys.getenv(
  "AGORA_TARGET_FILE",
  unset = if (exists("agora_target_file")) agora_target_file else file.path(external_data_dir, "AGORA_nominated_targets.csv")
)

sf3_load_agora_targets <- function(agora_path = sf3_agora_path) {
  if (!file.exists(agora_path)) {
    stop(
      "AGORA target file not found:\n",
      agora_path,
      "\nUpdate sf3_agora_path in 12_make_supplementary_figure_3_mitochondrial_specificity.R.",
      call. = FALSE
    )
  }
  
  agora_raw <- readr::read_csv(agora_path, show_col_types = FALSE)
  
  gene_col <- dplyr::case_when(
    "Gene Symbol" %in% colnames(agora_raw) ~ "Gene Symbol",
    "gene_symbol" %in% colnames(agora_raw) ~ "gene_symbol",
    "gene" %in% colnames(agora_raw) ~ "gene",
    "Gene" %in% colnames(agora_raw) ~ "Gene",
    TRUE ~ NA_character_
  )
  
  if (is.na(gene_col)) {
    stop(
      "AGORA file lacks a recognizable gene column. Available columns: ",
      paste(colnames(agora_raw), collapse = ", "),
      call. = FALSE
    )
  }
  
  nomination_col <- dplyr::case_when(
    "Nominations" %in% colnames(agora_raw) ~ "Nominations",
    "nominations" %in% colnames(agora_raw) ~ "nominations",
    "agora_nominated_target" %in% colnames(agora_raw) ~ "agora_nominated_target",
    TRUE ~ NA_character_
  )
  
  agora_targets <- agora_raw |>
    dplyr::transmute(
      gene = clean_gene(.data[[gene_col]]),
      nomination_value = if (!is.na(nomination_col)) {
        suppressWarnings(as.numeric(.data[[nomination_col]]))
      } else {
        1
      }
    ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(
      agora_nominated_target = as.numeric(any(.data$nomination_value > 0, na.rm = TRUE)),
      .groups = "drop"
    )
  
  message(
    "Loaded AGORA nominated targets from: ", agora_path,
    "\nAGORA target genes loaded: ", sum(agora_targets$agora_nominated_target == 1, na.rm = TRUE),
    " / ", nrow(agora_targets)
  )
  
  agora_targets
}

sf3_reannotate_agora_from_csv <- function(tbl, agora_targets) {
  tbl |>
    tibble::as_tibble() |>
    dplyr::mutate(gene = clean_gene(.data$gene)) |>
    dplyr::select(-dplyr::any_of("agora_nominated_target")) |>
    dplyr::left_join(agora_targets, by = "gene") |>
    dplyr::mutate(
      agora_nominated_target = dplyr::coalesce(.data$agora_nominated_target, 0),
      agora_nominated_target = as.numeric(.data$agora_nominated_target)
    )
}

############################################################
## AGORA re-annotation helper
############################################################

sf3_reannotate_agora <- function(tbl, inputs) {
  tbl <- tibble::as_tibble(tbl) |>
    dplyr::mutate(gene = clean_gene(.data$gene))
  
  agora_candidates <- c(
    "agora_targets",
    "agora_target_tbl",
    "agora_raw",
    "agora_tbl",
    "agora_nominated_targets",
    "agora_target_table"
  )
  
  agora_obj_name <- agora_candidates[
    agora_candidates %in% names(inputs) |
      vapply(agora_candidates, exists, logical(1), envir = .GlobalEnv)
  ][1]
  
  if (is.na(agora_obj_name)) {
    warning(
      "No AGORA target table found in inputs/global environment. ",
      "Keeping existing agora_nominated_target column.",
      call. = FALSE
    )
    return(tbl)
  }
  
  agora_obj <- if (agora_obj_name %in% names(inputs)) {
    inputs[[agora_obj_name]]
  } else {
    get(agora_obj_name, envir = .GlobalEnv)
  }
  
  agora_tbl <- tibble::as_tibble(agora_obj)
  
  gene_col <- dplyr::case_when(
    "gene" %in% colnames(agora_tbl) ~ "gene",
    "gene_symbol" %in% colnames(agora_tbl) ~ "gene_symbol",
    "Gene Symbol" %in% colnames(agora_tbl) ~ "Gene Symbol",
    "Gene" %in% colnames(agora_tbl) ~ "Gene",
    TRUE ~ NA_character_
  )
  
  if (is.na(gene_col)) {
    stop(
      "AGORA object found but no recognizable gene column. Available columns: ",
      paste(colnames(agora_tbl), collapse = ", "),
      call. = FALSE
    )
  }
  
  nomination_col <- dplyr::case_when(
    "Nominations" %in% colnames(agora_tbl) ~ "Nominations",
    "nominations" %in% colnames(agora_tbl) ~ "nominations",
    "agora_nominated_target" %in% colnames(agora_tbl) ~ "agora_nominated_target",
    TRUE ~ NA_character_
  )
  
  agora_targets <- agora_tbl |>
    dplyr::transmute(
      gene = clean_gene(.data[[gene_col]]),
      agora_nominated_target_new = dplyr::case_when(
        !is.na(nomination_col) ~ suppressWarnings(as.numeric(.data[[nomination_col]])) > 0,
        TRUE ~ TRUE
      )
    ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(
      agora_nominated_target_new = as.numeric(any(.data$agora_nominated_target_new, na.rm = TRUE)),
      .groups = "drop"
    )
  
  out <- tbl |>
    dplyr::select(-dplyr::any_of("agora_nominated_target")) |>
    dplyr::left_join(agora_targets, by = "gene") |>
    dplyr::mutate(
      agora_nominated_target = dplyr::coalesce(.data$agora_nominated_target_new, 0),
      agora_nominated_target = as.numeric(.data$agora_nominated_target)
    ) |>
    dplyr::select(-.data$agora_nominated_target_new)
  
  message(
    "Re-annotated AGORA targets using ", agora_obj_name,
    ": ", sum(out$agora_nominated_target == 1, na.rm = TRUE),
    " / ", nrow(out), " genes flagged."
  )
  
  out
}

# ============================================================
# Matched null generation
# ============================================================

sf3_metric_tbl <- tibble::tibble(
  metric = c(
    "protein_collapse_magnitude",
    "inverse_braak_magnitude",
    "cognition_composite_score",
    "agora_nominated_target"
  ),
  metric_label = c(
    "Late-stage protein collapse",
    "Inverse Braak association",
    "Cognition association score",
    "AGORA target fraction"
  ),
  metric_type = c("continuous", "continuous", "continuous", "fraction")
)

sf3_format_empirical_p <- function(extreme_count, n_perm) {
  if (is.na(extreme_count) || is.na(n_perm) || n_perm <= 0) {
    return("empirical P unavailable")
  }
  if (extreme_count == 0) {
    return(paste0("empirical P < ", formatC(1 / n_perm, format = "e", digits = 1)))
  }
  p <- extreme_count / n_perm
  if (p < 0.001) {
    return(paste0("empirical P = ", formatC(p, format = "e", digits = 1)))
  }
  paste0("empirical P = ", signif(p, 3))
}

sf3_format_null_position <- function(extreme_count, n_perm) {
  if (is.na(extreme_count) || is.na(n_perm) || n_perm <= 0) {
    return("not estimable")
  }
  if (extreme_count == 0) {
    return(paste0("0 / ", scales::comma(n_perm), " null sets exceeded observed"))
  }
  paste0(extreme_count, " / ", scales::comma(n_perm), " null sets exceeded observed")
}


sf3_format_emp_p <- function(extreme_count, n_perm) {
  if (is.na(extreme_count) || is.na(n_perm) || n_perm <= 0) {
    return("NA")
  }
  if (extreme_count == 0) {
    return(paste0("<1/", scales::comma(n_perm)))
  }
  paste0(extreme_count, "/", scales::comma(n_perm))
}

sf3_format_emp_p_decimal <- function(extreme_count, n_perm) {
  if (is.na(extreme_count) || is.na(n_perm) || n_perm <= 0) {
    return(NA_character_)
  }
  if (extreme_count == 0) {
    return(paste0("<", formatC(1 / n_perm, format = "e", digits = 2)))
  }
  p <- extreme_count / n_perm
  if (p < 0.001) {
    formatC(p, format = "e", digits = 2)
  } else {
    as.character(signif(p, 3))
  }
}

sf3_clean_sci_label <- function(x, digits = 1) {
  out <- formatC(x, format = "e", digits = digits)
  out <- sub("\\.0e", "e", out)
  out <- sub("e\\+0*", "e", out)
  sub("e-0*", "e-", out)
}

sf3_format_compact_p <- function(extreme_count, n_perm) {
  if (is.na(extreme_count) || is.na(n_perm) || n_perm <= 0) {
    return("p = NA")
  }
  if (extreme_count == 0) {
    return(paste0("p < ", sf3_clean_sci_label(1 / n_perm, digits = 1)))
  }
  p <- extreme_count / n_perm
  if (p < 0.001) {
    return(paste0("p = ", sf3_clean_sci_label(p, digits = 1)))
  }
  paste0("p = ", signif(p, 3))
}

sf3_null_position_label <- function(extreme_count, n_perm) {
  if (is.na(extreme_count) || is.na(n_perm) || n_perm <= 0) {
    return("not estimable")
  }
  if (extreme_count == 0) {
    return(paste0("observed > all ", scales::comma(n_perm), " null sets"))
  }
  pct <- 100 * (1 - extreme_count / n_perm)
  paste0("observed at ", signif(pct, 3), "th null percentile")
}


sf3_make_matched_null <- function(hsp_tbl, background_tbl, n_perm = 5000, seed = 1300) {
  set.seed(seed)

  hsp_bin_counts <- hsp_tbl |>
    dplyr::filter(!is.na(.data$abundance_bin)) |>
    dplyr::count(.data$abundance_bin, name = "n_hsp")

  background_by_bin <- split(background_tbl, background_tbl$abundance_bin)

  missing_bins <- setdiff(hsp_bin_counts$abundance_bin, names(background_by_bin))
  if (length(missing_bins) > 0) {
    stop(
      "Background null pool lacks abundance bins needed to match Hsp60/10 clients: ",
      paste(missing_bins, collapse = ", "),
      call. = FALSE
    )
  }

  null_sets <- purrr::map_dfr(seq_len(n_perm), function(i) {
    sampled_tbl <- purrr::map_dfr(seq_len(nrow(hsp_bin_counts)), function(j) {
      bin_j <- hsp_bin_counts$abundance_bin[j]
      n_j <- hsp_bin_counts$n_hsp[j]
      pool_j <- background_by_bin[[bin_j]]

      # Use replacement only when the matched pool in that bin is smaller than required.
      sampled_idx <- sample(seq_len(nrow(pool_j)), size = n_j, replace = nrow(pool_j) < n_j)
      pool_j[sampled_idx, , drop = FALSE]
    })

    sf3_metric_tbl |>
      dplyr::mutate(
        perm = i,
        null_mean = purrr::map_dbl(.data$metric, ~ mean(sampled_tbl[[.x]], na.rm = TRUE)),
        null_n = nrow(sampled_tbl)
      )
  })

  null_sets
}

sf3_observed_summary <- function(hsp_tbl) {
  sf3_metric_tbl |>
    dplyr::mutate(
      observed_mean = purrr::map_dbl(.data$metric, ~ mean(hsp_tbl[[.x]], na.rm = TRUE)),
      observed_n = purrr::map_int(.data$metric, ~ sum(!is.na(hsp_tbl[[.x]])))
    )
}

sf3_null_summary <- function(null_tbl, observed_tbl) {
  null_tbl2 <- null_tbl |>
    dplyr::rename(null_perm_mean = .data$null_mean)
  
  null_tbl2 |>
    dplyr::group_by(.data$metric, .data$metric_label, .data$metric_type) |>
    dplyr::summarise(
      null_mean = mean(.data$null_perm_mean, na.rm = TRUE),
      null_sd = stats::sd(.data$null_perm_mean, na.rm = TRUE),
      null_q025 = stats::quantile(.data$null_perm_mean, 0.025, na.rm = TRUE),
      null_q975 = stats::quantile(.data$null_perm_mean, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(observed_tbl, by = c("metric", "metric_label", "metric_type")) |>
    dplyr::mutate(
      observed_to_null_ratio = .data$observed_mean / .data$null_mean,
      n_perm = purrr::map_int(.data$metric, function(m) {
        vals <- null_tbl2$null_perm_mean[null_tbl2$metric == m]
        sum(is.finite(vals))
      }),
      extreme_count_greater = purrr::map_int(.data$metric, function(m) {
        obs <- observed_tbl$observed_mean[observed_tbl$metric == m][1]
        vals <- null_tbl2$null_perm_mean[null_tbl2$metric == m]
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0 || !is.finite(obs)) return(NA_integer_)
        sum(vals >= obs)
      }),
      empirical_p_greater = (.data$extreme_count_greater + 1) / (.data$n_perm + 1),
      empirical_p_label = purrr::map2_chr(
        .data$extreme_count_greater,
        .data$n_perm,
        sf3_format_empirical_p
      ),
      null_position_label = purrr::map2_chr(
        .data$extreme_count_greater,
        .data$n_perm,
        sf3_format_null_position
      )
    )
}

# ============================================================
# Plot helpers
# ============================================================

sf3_make_schematic <- function() {
  nodes <- tibble::tibble(
    label = c(
      "Observed\nHsp60/10\nclients",
      "Match background\nproteins by\nabundance bin",
      "Generate matched\nmitochondrial\nnull sets",
      "Compare\nobserved to\nnull expectation"
    ),
    x = c(0.95, 3.55, 6.35, 9.55),
    y = 1,
    fill = c("#F7DEE2", "#E7EDF5", "#E7EDF5", "#F2F2F2")
  ) |>
    dplyr::mutate(
      width = dplyr::case_when(
        .data$label == "Match background\nproteins by\nabundance bin" ~ 0.95,
        .data$label == "Generate matched\nmitochondrial\nnull sets" ~ 1.18,
        .data$label == "Compare\nobserved to\nnull expectation" ~ 1.10,
        TRUE ~ 0.78
      ),
      height = 0.33
    )

  edges <- tibble::tibble(
    x = nodes$x[-nrow(nodes)] + nodes$width[-nrow(nodes)] + 0.18,
    xend = nodes$x[-1] - nodes$width[-1] - 0.20,
    y = 1,
    yend = 1
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      arrow = grid::arrow(length = grid::unit(0.075, "inches"), type = "closed"),
      linewidth = 0.42,
      color = "grey35"
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(
        xmin = .data$x - .data$width,
        xmax = .data$x + .data$width,
        ymin = .data$y - .data$height,
        ymax = .data$y + .data$height,
        fill = .data$fill
      ),
      color = "grey25",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(nodes, .data$label != "Match background\nproteins by\nabundance bin"),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.55,
      lineheight = 0.9,
      color = sf3_text
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(nodes, .data$label == "Match background\nproteins by\nabundance bin"),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.4,
      lineheight = 0.9,
      color = sf3_text
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.05, 10.85), ylim = c(0.45, 1.55), clip = "off") +
    ggplot2::labs(title = "Matched mitochondrial null design") +
    sf3_theme(base_size = 8.5) +
    ggplot2::theme_void(base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(8, 12, 4, 12)
    )
}

sf3_make_null_hist <- function(null_tbl, summary_tbl, metric_name, title, xlab) {
  plot_null <- null_tbl |>
    dplyr::filter(.data$metric == metric_name)

  plot_sum <- summary_tbl |>
    dplyr::filter(.data$metric == metric_name)

  obs <- plot_sum$observed_mean[1]
  null_mean <- plot_sum$null_mean[1]
  label_txt <- paste0(
    "observed = ", signif(obs, 3),
    "\nnull mean = ", signif(null_mean, 3),
    "\n", sf3_format_compact_p(plot_sum$extreme_count_greater[1], plot_sum$n_perm[1])
  )
  x_rng <- range(c(plot_null$null_mean, obs, null_mean), na.rm = TRUE)
  ann_shift <- 0.11 * diff(x_rng)
  ann_x <- if (obs > null_mean) obs - ann_shift else obs + ann_shift

  ggplot2::ggplot(plot_null, ggplot2::aes(x = .data$null_mean)) +
    ggplot2::geom_histogram(
      bins = 42,
      fill = sf3_null_fill,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(xintercept = null_mean, color = sf3_null_line, linewidth = 0.75) +
    ggplot2::geom_vline(xintercept = obs, color = sf3_hsp_col, linewidth = 0.95, linetype = "dashed") +
    ggplot2::annotate(
      "label",
      x = ann_x,
      y = Inf,
      label = label_txt,
      hjust = ifelse(obs > null_mean, 1, 0),
      vjust = 1.22,
      size = 3.1,
      label.size = 0.25,
      fill = "white",
      color = "grey10"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = NULL,
      x = xlab,
      y = NULL
    ) +
    sf3_theme(base_size = 8.5)
}

sf3_make_ratio_summary <- function(summary_tbl) {
  plot_tbl <- summary_tbl |>
    dplyr::mutate(
      metric_label = factor(
        .data$metric_label,
        levels = rev(c(
          "Late-stage protein collapse",
          "Inverse Braak association",
          "Cognition association score",
          "AGORA target fraction"
        ))
      ),
      label = paste0(
        signif(.data$observed_to_null_ratio, 2),
        "x, ",
        purrr::map2_chr(.data$extreme_count_greater, .data$n_perm, sf3_format_compact_p)
      )
    )

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$observed_to_null_ratio, y = .data$metric_label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.55) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 1, xend = .data$observed_to_null_ratio, yend = .data$metric_label),
      color = "grey60",
      linewidth = 0.75
    ) +
    ggplot2::geom_point(size = 3.2, color = sf3_hsp_col) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      nudge_x = 0.03,
      hjust = -0.05,
      size = 3.0,
      color = "grey15"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.28))) +
    ggplot2::labs(
      title = "Specificity summary",
      subtitle = NULL,
      x = "Observed / matched-null mean",
      y = NULL
    ) +
    sf3_theme(base_size = 8.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.25)
    )
}

# ============================================================
# Main builder
# ============================================================

make_supfig3_mitochondrial_specificity <- function(inputs, n_perm = 10000, seed = 1300) {
  message("Starting Supplementary Figure 3: mitochondrial specificity.")

  hsp_input <- sf3_get_input(
    c("hsp_null_tbl", "hsp60_10_null_tbl", "hsp_specificity_tbl"),
    required = TRUE,
    label = "observed Hsp60/10 null comparison table"
  )

  background_input <- sf3_get_input(
    c("background_null_pool", "mitochondrial_background_null_pool", "background_specificity_pool"),
    required = TRUE,
    label = "non-Hsp60/10 mitochondrial background null pool"
  )

  agora_targets <- sf3_load_agora_targets()
  
  hsp_tbl <- sf3_validate_tbl(hsp_input$object, "hsp_null_tbl") |>
    sf3_reannotate_agora_from_csv(agora_targets)
  
  background_tbl <- sf3_validate_tbl(background_input$object, "background_null_pool") |>
    sf3_reannotate_agora_from_csv(agora_targets) |>
    dplyr::filter(!.data$gene %in% hsp_tbl$gene)

  message("Observed Hsp60/10 table: ", nrow(hsp_tbl), " genes")
  message("Background null pool after excluding Hsp60/10 genes: ", nrow(background_tbl), " genes")
  message("Null permutations: ", n_perm)

  bin_audit <- dplyr::full_join(
    hsp_tbl |> dplyr::count(.data$abundance_bin, name = "n_hsp"),
    background_tbl |> dplyr::count(.data$abundance_bin, name = "n_background"),
    by = "abundance_bin"
  ) |>
    dplyr::arrange(.data$abundance_bin)

  readr::write_csv(bin_audit, file.path(audits_dir, "SuppFig3_abundance_bin_matching_audit.csv"))

  observed_tbl <- sf3_observed_summary(hsp_tbl)
  null_tbl <- sf3_make_matched_null(hsp_tbl, background_tbl, n_perm = n_perm, seed = seed)
  summary_tbl <- sf3_null_summary(null_tbl, observed_tbl)

  readr::write_csv(observed_tbl, file.path(audits_dir, "SuppFig3_observed_hsp60_10_summary.csv"))
  readr::write_csv(null_tbl, file.path(audits_dir, "SuppFig3_matched_null_permutation_values.csv"))
  readr::write_csv(summary_tbl, file.path(audits_dir, "SuppFig3_specificity_summary.csv"))

  panel_a <- sf3_panel_label(sf3_make_schematic(), "A")
  panel_b <- sf3_panel_label(
    sf3_make_null_hist(
      null_tbl,
      summary_tbl,
      metric_name = "protein_collapse_magnitude",
      title = "Late-stage collapse specificity",
      xlab = "Mean collapse magnitude"
    ),
    "B"
  )
  panel_c <- sf3_panel_label(
    sf3_make_null_hist(
      null_tbl,
      summary_tbl,
      metric_name = "inverse_braak_magnitude",
      title = "Inverse Braak specificity",
      xlab = "Mean inverse Braak magnitude"
    ),
    "C"
  )
  panel_d <- sf3_panel_label(sf3_make_ratio_summary(summary_tbl), "D")

  panel_paths <- c(
    save_panel_set(panel_a, "SuppFig3A_matched_null_design", width = 7.2, height = 2.4, output_dir = panels_dir),
    save_panel_set(panel_b, "SuppFig3B_collapse_null", width = 5.9, height = 4.0, output_dir = panels_dir),
    save_panel_set(panel_c, "SuppFig3C_inverse_braak_null", width = 5.9, height = 4.0, output_dir = panels_dir),
    save_panel_set(panel_d, "SuppFig3D_specificity_summary", width = 6.3, height = 4.0, output_dir = panels_dir)
  )

  composite <- panel_a /
    (panel_b | panel_c) /
    panel_d +
    patchwork::plot_layout(heights = c(0.55, 1, 0.88))

  save_plot_set(
    composite,
    "Supplementary_Figure_3_mitochondrial_specificity",
    width = 12.4,
    height = 10.2,
    output_dir = figures_dir
  )


  run_summary <- c(
    "Supplementary Figure 3 run summary",
    paste0("generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "figure: specificity against mitochondrial background proteins",
    paste0("observed Hsp60/10 genes: ", nrow(hsp_tbl)),
    paste0("background pool genes: ", nrow(background_tbl)),
    paste0("null permutations: ", n_perm),
    "summary:",
    paste0(
      "  - ", summary_tbl$metric_label,
      ": observed=", signif(summary_tbl$observed_mean, 4),
      ", null mean=", signif(summary_tbl$null_mean, 4),
      ", obs/null=", signif(summary_tbl$observed_to_null_ratio, 4),
      ", empirical_p_label=", summary_tbl$empirical_p_label,
      ", null_position=", summary_tbl$null_position_label,
      collapse = "\n"
    )
  )

  writeLines(run_summary, file.path(audits_dir, "SuppFig3_run_summary.txt"))

  message("Supplementary Figure 3 complete.")
  print(summary_tbl)

  invisible(list(
    hsp_tbl = hsp_tbl,
    background_tbl = background_tbl,
    observed_tbl = observed_tbl,
    null_tbl = null_tbl,
    summary_tbl = summary_tbl,
    panels = list(A = panel_a, B = panel_b, C = panel_c, D = panel_d),
    panel_paths = panel_paths
  ))
}

supfig3_outputs <- make_supfig3_mitochondrial_specificity(inputs)
