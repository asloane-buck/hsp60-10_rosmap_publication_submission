# S6: plot frozen formal paired DLPFC/STG interaction results.
# Genome-wide FDR values are read unchanged; no regional models are refitted.
if (!exists("project_dir") || !exists("save_plot_set")) {
  stop("Run supplemental configuration and helpers first.", call. = FALSE)
}
sf6_root <- file.path(project_dir, "outputs", "reviewer_revisions", "regional_formal_interaction")
sf6_paths <- c(
  AD = file.path(sf6_root, "validated_final", "FINAL_region_x_AD_interaction_Hsp60_clients.csv"),
  Braak = file.path(sf6_root, "validated_final", "FINAL_region_x_Braak_continuous_interaction_Hsp60_clients.csv"),
  summary = file.path(sf6_root, "validated_final", "FINAL_regional_interaction_summary.csv"),
  pathway = file.path(sf6_root, "hsp_conclusion_audit", "Hsp_pathway_score_primary_interactions.csv")
)
if (any(!file.exists(sf6_paths))) {
  stop("S6 frozen source files are missing:\n", paste(sf6_paths[!file.exists(sf6_paths)], collapse = "\n"), call. = FALSE)
}
sf6_sources <- lapply(sf6_paths, readr::read_csv, show_col_types = FALSE)
sf6_summary <- sf6_sources$summary
sf6_expected <- c(AD_models = 9005, AD_genomewide_FDR_lt_0.05 = 253,
  AD_Hsp_tested = 300, AD_Hsp_genomewide_FDR_lt_0.05 = 6,
  Braak_continuous_models = 9002, Braak_continuous_genomewide_FDR_lt_0.05 = 77,
  Braak_continuous_Hsp_tested = 300, Braak_continuous_Hsp_genomewide_FDR_lt_0.05 = 1,
  Braak_bin3_models = 9001, Braak_bin3_global_FDR_lt_0.05 = 53,
  Braak_bin3_Hsp_tested = 300)
if (!all(c("metric", "value") %in% names(sf6_summary)) || anyDuplicated(sf6_summary$metric)) {
  stop("S6 summary must have unique metric/value rows.", call. = FALSE)
}
sf6_actual <- setNames(sf6_summary$value, sf6_summary$metric)[names(sf6_expected)]
if (anyNA(sf6_actual) || any(sf6_actual != sf6_expected)) {
  stop("S6 summary does not match the verified final regional analysis.", call. = FALSE)
}
sf6_prepare <- function(tbl, n_sig) {
  cols <- c("accession", "gene_symbol", "beta_region_x_predictor", "p_region_x_predictor_fdr")
  if (!all(cols %in% names(tbl)) || nrow(tbl) != 300 || anyDuplicated(tbl$accession) ||
      any(!is.finite(tbl$beta_region_x_predictor)) || any(!is.finite(tbl$p_region_x_predictor_fdr)) ||
      any(tbl$p_region_x_predictor_fdr < 0 | tbl$p_region_x_predictor_fdr > 1) ||
      sum(tbl$p_region_x_predictor_fdr < .05) != n_sig) {
    stop("S6 client interaction table failed verification.", call. = FALSE)
  }
  tbl |>
    dplyr::mutate(
      q = .data$p_region_x_predictor_fdr,
      significance = factor(ifelse(.data$q < .05, "Genome-wide FDR < 0.05", "Other clients"),
        levels = c("Other clients", "Genome-wide FDR < 0.05")),
      minus_log10_q = -log10(pmax(.data$q, .Machine$double.xmin))
    )
}
sf6_ad <- sf6_prepare(sf6_sources$AD, 6)
sf6_braak <- sf6_prepare(sf6_sources$Braak, 1)
sf6_theme <- function() {
  ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "grey35", linewidth = .35),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 10),
      strip.background = ggplot2::element_rect(fill = "grey96", color = "grey40", linewidth = .3),
      strip.text = ggplot2::element_text(face = "bold", size = 8)
    )
}
sf6_volcano <- function(tbl, tag, xlab, n_total, n_total_sig) {
  # Label all FDR-significant clients; use repulsion if already installed.
  selected <- tbl[tbl$q < .05, ]
  p <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$beta_region_x_predictor, y = .data$minus_log10_q)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = .3) +
    ggplot2::geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "grey50", linewidth = .35) +
    ggplot2::geom_point(ggplot2::aes(color = .data$significance), size = 1.2, alpha = .75) +
    ggplot2::scale_color_manual(values = c("Other clients" = "#A9AFB8", "Genome-wide FDR < 0.05" = "#B23A48")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.03, .38))) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -.03, vjust = 1.1, size = 2.5,
      label = paste0("Clients: ", nrow(selected), "/300 at FDR < 0.05\nGenome-wide: ",
        n_total_sig, "/", scales::comma(n_total))) +
    ggplot2::labs(tag = tag, x = xlab, y = expression(-log[10]("genome-wide FDR"))) + sf6_theme()
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = selected, ggplot2::aes(label = .data$gene_symbol),
      size = 2.5, max.overlaps = Inf, min.segment.length = 0, seed = 1300, show.legend = FALSE)
  } else {
    p <- p + ggplot2::geom_text(data = selected, ggplot2::aes(label = .data$gene_symbol),
      size = 2.5, vjust = -1, check_overlap = TRUE, show.legend = FALSE)
  }
  p
}
sf6_pathway <- sf6_sources$pathway
sf6_columns <- c("analysis", "n_people", "beta_region_x_predictor", "ci_low_region_x_predictor",
  "ci_high_region_x_predictor", "p_region_x_predictor")
if (!all(sf6_columns %in% names(sf6_pathway))) stop("S6 pathway table lacks required columns.", call. = FALSE)
sf6_pathway <- sf6_pathway |>
  dplyr::filter(.data$analysis %in% c("Hsp_score_region_x_AD", "Hsp_score_region_x_Braak_continuous")) |>
  dplyr::mutate(model = factor(ifelse(.data$analysis == "Hsp_score_region_x_AD", "AD vs control", "Continuous Braak"),
    levels = c("AD vs control", "Continuous Braak")), label = sprintf("n = %d paired participants; P = %.3f", .data$n_people, .data$p_region_x_predictor)) |>
  dplyr::arrange(.data$model)
if (nrow(sf6_pathway) != 2 || anyDuplicated(sf6_pathway$analysis) ||
    any(sf6_pathway$n_people != c(138, 134)) ||
    any(abs(sf6_pathway$beta_region_x_predictor - c(.155131839749539, .021931406970494)) > 1e-10) ||
    any(abs(sf6_pathway$p_region_x_predictor - c(.175045202908122, .468183445138642)) > 1e-10)) {
  stop("S6 pathway results differ from the verified final analysis.", call. = FALSE)
}
sf6_panel_a <- sf6_volcano(sf6_ad, "A", "Region x AD interaction beta", 9005, 253)
sf6_panel_b <- sf6_volcano(sf6_braak, "B", "Region x Braak interaction beta", 9002, 77)
sf6_panel_c <- ggplot2::ggplot(sf6_pathway, ggplot2::aes(y = 1, x = .data$beta_region_x_predictor)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = .35) +
  ggplot2::geom_segment(ggplot2::aes(x = .data$ci_low_region_x_predictor,
    xend = .data$ci_high_region_x_predictor, yend = 1), linewidth = .6, color = "#B23A48") +
  ggplot2::geom_point(size = 2.3, color = "#B23A48") +
  ggplot2::geom_text(ggplot2::aes(label = .data$label), y = 1.38, size = 2.5) +
  ggplot2::facet_wrap(~model, nrow = 1, scales = "free_x") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2)) +
  ggplot2::coord_cartesian(ylim = c(.7, 1.6)) +
  ggplot2::labs(tag = "C", x = "Client-network interaction beta (95% CI)", y = NULL) +
  sf6_theme() + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank())
sf6_final <- (sf6_panel_a | sf6_panel_b) / sf6_panel_c +
  patchwork::plot_layout(heights = c(1.5, 1), guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
save_plot_set(sf6_final, "Supplementary_Figure_6_regional_proteomics_validation",
  width = 170/25.4, height = 140/25.4, output_dir = figures_dir)
readr::write_csv(dplyr::bind_rows(dplyr::mutate(sf6_ad, model = "AD"),
  dplyr::mutate(sf6_braak, model = "continuous_Braak")),
  file.path(audits_dir, "SuppFig6_formal_client_interactions_plotted.csv"))
readr::write_csv(sf6_pathway, file.path(audits_dir, "SuppFig6_formal_pathway_interactions_plotted.csv"))
readr::write_csv(sf6_summary, file.path(audits_dir, "SuppFig6_formal_regional_summary_used.csv"))
readr::write_csv(tibble::tibble(source = unname(sf6_paths), md5 = unname(tools::md5sum(sf6_paths))),
  file.path(audits_dir, "SuppFig6_frozen_source_manifest.csv"))
supfig6_outputs <- list(plot = sf6_final, AD = sf6_ad, Braak = sf6_braak,
  pathway = sf6_pathway, summary = sf6_summary, panels = list(A = sf6_panel_a, B = sf6_panel_b, C = sf6_panel_c))
message("Supplementary Figure 6 complete: formal paired interaction sources; no models refitted.")
