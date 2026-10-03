############################################################
## 11_make_main_figure_2_collapse_heterogeneity.R
## Literal migration of the final old plotting block.
## Inputs are supplied by adjusted upstream scripts.
############################################################

## MAIN FIGURE 2 — ALL HSP60/10 CLIENTS
## Heterogeneous late-stage decline and network centrality
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
})

############################################################
## 0. Output setup
############################################################

fig_id <- "main_fig2_all_clients"

## Keep Figure 2 outputs inside the covariate-adjusted modular output root.
fig2_all_dir <- file.path(cfg$output_root, fig_id)
fig2_all_plot_dir <- file.path(cfg$plot_dir, fig_id)
fig2_all_table_dir <- file.path(cfg$table_dir, fig_id)

dir.create(fig2_all_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig2_all_table_dir, recursive = TRUE, showWarnings = FALSE)

write_fig2_all_table <- function(x, name) {
  out <- file.path(fig2_all_table_dir, paste0(name, ".csv"))
  readr::write_csv(x, out)
  message("Wrote table: ", out)
  invisible(x)
}

save_fig2_all_plot <- function(plot, name, width = 15, height = 9, dpi = 600) {
  pdf_out <- file.path(fig2_all_plot_dir, paste0(name, ".pdf"))
  png_out <- file.path(fig2_all_plot_dir, paste0(name, ".png"))
  svg_out <- file.path(fig2_all_plot_dir, paste0(name, ".svg"))

  ggsave(
    pdf_out,
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    bg = "white",
    limitsize = FALSE,
    useDingbats = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_out,
      plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_png,
      bg = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_out,
      plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      svg_out,
      plot,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )
  }

  message("Saved PDF: ", pdf_out)
  message("Saved PNG: ", png_out)
  if (file.exists(svg_out)) message("Saved SVG: ", svg_out)

  invisible(plot)
}

############################################################
## 1. Helpers
############################################################

clean_gene <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\\s+", "") %>%
    toupper()
}

p_to_label <- function(p) {
  case_when(
    is.na(p) ~ "P = NA",
    p < 0.001 ~ paste0("P = ", formatC(p, format = "e", digits = 2)),
    TRUE ~ paste0("P = ", signif(p, 2))
  )
}

star_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

theme_fig2 <- function(base_size = 11) {
  if (exists("paper_theme_rebuilt")) {
    paper_theme_rebuilt(base_size)
  } else {
    theme_classic(base_size = base_size) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      )
  }
}

############################################################
## 2. Build all-client Figure 2 data
############################################################

require_objects("all_hsp60_10_client_tbl", context = "11_make_main_figure_2_collapse_heterogeneity.R")
require_columns(
  all_hsp60_10_client_tbl,
  c(
    "gene",
    "detected_in_protein",
    "detected_in_rna",
    "protein_late_decline_magnitude",
    "hub_mean_abs_cor"
  ),
  "Figure 2 source all-client table"
)

fig2_tbl <- all_hsp60_10_client_tbl %>%
  mutate(
    gene = clean_gene(gene),
    detected_in_protein = as.logical(detected_in_protein),
    detected_in_rna = as.logical(detected_in_rna),
    
    late_decline_value = protein_late_decline_magnitude,
    
    centrality_value = hub_mean_abs_cor
  ) %>%
  filter(
    detected_in_protein == TRUE,
    is.finite(late_decline_value)
  ) %>%
  arrange(desc(late_decline_value)) %>%
  mutate(
    late_decline_rank = row_number(),
    n_detected_clients = n(),

    ## Top quartile among detected protein clients.
    late_decline_group = if_else(
      late_decline_rank <= ceiling(0.25 * n_detected_clients),
      "Top-quartile late decline",
      "Other detected Hsp60/10 clients"
    ),

    late_decline_group = factor(
      late_decline_group,
      levels = c("Top-quartile late decline", "Other detected Hsp60/10 clients")
    ),

     label_gene = if_else(
      late_decline_rank <= 10 |
        gene %in% c("SHMT2", "IDH3A", "TUFM", "ECHS1", "DLD", "DLAT", "LRPPRC"),
      gene,
      NA_character_
    )
  )

fig2_label_tbl <- fig2_tbl %>%
  filter(!is.na(label_gene)) %>%
  arrange(late_decline_rank) %>%
  slice_head(n = 17) %>%
  mutate(
    repel_nudge_x = case_when(
      label_gene %in% c("MRPL47", "MRPS33", "DAP3", "MRPS16", "MRPS9", "MRPS23") ~ -6,
      label_gene %in% c("PTCD3", "MRPS35", "NDUFAF7", "MRPS22") ~ 18,
      label_gene %in% c("LRPPRC", "TUFM", "IDH3A", "SHMT2") ~ -10,
      label_gene %in% c("DLD", "DLAT", "ECHS1") ~ 12,
      TRUE ~ 8
    ),
    repel_nudge_y = case_when(
      label_gene == "MRPL47" ~ 0.022,
      label_gene == "MRPS33" ~ 0.016,
      label_gene == "DAP3" ~ 0.010,
      label_gene == "MRPS16" ~ 0.004,
      label_gene == "MRPS9" ~ -0.004,
      label_gene == "MRPS23" ~ -0.012,
      label_gene == "PTCD3" ~ 0.020,
      label_gene == "MRPS35" ~ 0.012,
      label_gene == "NDUFAF7" ~ 0.004,
      label_gene == "MRPS22" ~ -0.008,
      label_gene == "LRPPRC" ~ -0.016,
      label_gene == "TUFM" ~ -0.006,
      label_gene == "IDH3A" ~ -0.012,
      label_gene == "SHMT2" ~ -0.018,
      label_gene == "DLAT" ~ 0.007,
      label_gene == "DLD" ~ -0.004,
      label_gene == "ECHS1" ~ -0.010,
      TRUE ~ 0
    )
  )

fig2_ranked_output_tbl <- fig2_tbl %>%
  select(
    gene,
    late_decline_rank,
    late_decline_group,
    detected_in_protein,
    detected_in_rna,
    protein_late_decline_magnitude,
    late_decline_value,
    hub_mean_abs_cor,
    centrality_value,
    label_gene
  )

write_fig2_all_table(
  fig2_ranked_output_tbl,
  "Fig2_all_clients_ranked_client_table"
)

late_decline_group_counts <- fig2_tbl %>%
  count(late_decline_group, name = "n_clients")

write_fig2_all_table(late_decline_group_counts, "Fig2_all_clients_late_decline_group_counts")

print(late_decline_group_counts)

############################################################
## 3. Stats for Panels B and C
############################################################

centrality_stats_tbl <- fig2_tbl %>%
  filter(is.finite(centrality_value)) %>%
  summarise(
    comparison = "Top-quartile late decline vs other detected clients",
    metric = "Network centrality",
    n_top = sum(late_decline_group == "Top-quartile late decline"),
    n_other = sum(late_decline_group == "Other detected Hsp60/10 clients"),
    median_top = median(centrality_value[late_decline_group == "Top-quartile late decline"], na.rm = TRUE),
    median_other = median(centrality_value[late_decline_group == "Other detected Hsp60/10 clients"], na.rm = TRUE),
    p_value = suppressWarnings(
      wilcox.test(
        centrality_value[late_decline_group == "Top-quartile late decline"],
        centrality_value[late_decline_group == "Other detected Hsp60/10 clients"],
        exact = FALSE
      )$p.value
    ),
    p_label = p_to_label(p_value),
    star = star_label(p_value)
  )

fig2_stats_tbl <- centrality_stats_tbl

write_fig2_all_table(fig2_stats_tbl, "Fig2_all_clients_panel_stats")

print(fig2_stats_tbl)

############################################################
## 4. Plot colors
############################################################

late_decline_colors <- c(
  "Top-quartile late decline" = "#B22222",
  "Other detected Hsp60/10 clients" = "grey78"
)

late_decline_point_colors <- c(
  "Top-quartile late decline" = "#B22222",
  "Other detected Hsp60/10 clients" = "grey70"
)

############################################################
## 5. Panel A — all detected clients ranked by late-stage decline
############################################################

require_columns(
  fig2_tbl,
  c("gene", "late_decline_rank", "late_decline_value", "late_decline_group", "label_gene"),
  "Figure 2 Panel A ranked table"
)

pA <- ggplot(
  fig2_tbl,
  aes(x = late_decline_rank, y = late_decline_value)
) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      xend = late_decline_rank,
      y = 0,
      yend = late_decline_value,
      color = late_decline_group
    ),
    linewidth = 0.28,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = late_decline_group),
    shape = 21,
    size = 1.85,
    color = "grey35",
    stroke = 0.18
  ) +
  ggrepel::geom_text_repel(
    data = fig2_label_tbl,
    aes(label = label_gene),
    size = 2.55,
    fontface = "bold",
    color = "#7F1D1D",
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.22,
    box.padding = 0.68,
    point.padding = 0.42,
    force = 9.0,
    force_pull = 0.12,
    nudge_x = fig2_label_tbl$repel_nudge_x,
    nudge_y = fig2_label_tbl$repel_nudge_y,
    direction = "both",
    max.iter = 20000,
    max.time = 2,
    seed = 1
  ) +
  scale_fill_manual(values = late_decline_point_colors, name = "Late-decline group") +
  scale_color_manual(values = late_decline_point_colors, guide = "none") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.01, 0.035))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.48))
  ) +
  labs(
    title = "A. Late-stage decline is heterogeneous across the Hsp60/10 client network",
    subtitle = "All detected Hsp60/10 clients are ranked by late-stage protein decline; red marks the top quartile.",
    x = "Rank among detected Hsp60/10 clients",
    y = "Late-stage protein decline"
  ) +
  theme_fig2(11) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9.5),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

############################################################
############################################################
## 6. Panel B — centrality comparison
############################################################
############################################################

centrality_plot_tbl <- fig2_tbl %>%
  filter(is.finite(centrality_value)) %>%
  mutate(
      late_decline_group_short = case_when(
    late_decline_group == "Top-quartile late decline" ~ paste0("Top quartile\n(n=", sum(late_decline_group == "Top-quartile late decline"), ")"),
    late_decline_group == "Other detected Hsp60/10 clients" ~ paste0("Other clients\n(n=", sum(late_decline_group == "Other detected Hsp60/10 clients"), ")")
  ),
  late_decline_group_short = factor(
    late_decline_group_short,
    levels = c(
      paste0("Top quartile\n(n=", sum(late_decline_group == "Top-quartile late decline"), ")"),
      paste0("Other clients\n(n=", sum(late_decline_group == "Other detected Hsp60/10 clients"), ")")
    )
)
  )

require_columns(
  centrality_plot_tbl,
  c("late_decline_group_short", "centrality_value", "late_decline_group"),
  "Figure 2 Panel B centrality table"
)
require_values(centrality_plot_tbl, "late_decline_group", levels(fig2_tbl$late_decline_group), "Figure 2 Panel B")

centrality_ymax <- max(centrality_plot_tbl$centrality_value, na.rm = TRUE)
centrality_yrange <- diff(range(centrality_plot_tbl$centrality_value, na.rm = TRUE))
if (!is.finite(centrality_yrange) || centrality_yrange == 0) centrality_yrange <- 0.1

pB <- ggplot(
  centrality_plot_tbl,
  aes(x = late_decline_group_short, y = centrality_value, fill = late_decline_group)
) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    alpha = 0.75,
    color = "grey30",
    linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = late_decline_group),
    width = 0.13,
    size = 1.35,
    alpha = 0.45,
    show.legend = FALSE
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = centrality_ymax + centrality_yrange * 0.28,
    yend = centrality_ymax + centrality_yrange * 0.28,
    color = "black",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 1.5,
    y = centrality_ymax + centrality_yrange * 0.51,
    label = paste0(
      centrality_stats_tbl$star,
      "\n",
      centrality_stats_tbl$p_label
    ),
    size = 2.8,
    fontface = "bold",
    lineheight = 0.95
  ) +
  scale_fill_manual(values = late_decline_colors, guide = "none") +
  scale_color_manual(values = late_decline_point_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.70))) +
  labs(
    title = "B. Network centrality",
    subtitle = "Top-quartile late-declining clients occupy more central positions in the client network.",
    x = NULL,
    y = "Network centrality"
  ) +
  theme_fig2(11) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9.5),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

############################################################
## 7. Assemble Figure 2
############################################################




fig2_all_clients <- pA / pB +
  plot_layout(heights = c(1.35, 0.85)) +
  plot_annotation(
    title = "Hsp60/10 client vulnerability is selective across the full client network",
    subtitle = "Late-stage protein decline is heterogeneous; strongly late-declining clients show greater network centrality.",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey35")
    )
  )

save_fig2_all_plot(
  fig2_all_clients,
  "Main_Fig2_all_clients_late_decline_centrality",
  width = 15.5,
  height = 10.2
)

fig2_all_clients

fig2_missingness_audit <- fig2_tbl %>%
  summarise(
    n_total_detected_protein = n(),
    n_with_centrality = sum(is.finite(centrality_value)),
    n_missing_centrality = sum(!is.finite(centrality_value))
  )

write_fig2_all_table(fig2_missingness_audit, "Fig2_all_clients_metric_missingness_audit")
print(fig2_missingness_audit)

fig2_metric_missing_genes <- fig2_tbl %>%
  filter(!is.finite(centrality_value)) %>%
  select(
    any_of(c(
      "gene",
      "late_decline_rank",
      "late_decline_group",
      "detected_in_protein",
      "detected_in_rna",
      "protein_late_decline_magnitude",
      "centrality_value",
      "hub_mean_abs_cor"
    ))
  ) %>%
  arrange(late_decline_rank)

write_fig2_all_table(fig2_metric_missing_genes, "Fig2_all_clients_metric_missing_genes")
print(fig2_metric_missing_genes, n = 50)

############################################################
