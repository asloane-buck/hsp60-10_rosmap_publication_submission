############################################################
## 12_make_main_figure_3_pathology_coupling.R
## MAIN FIGURE 3 — ALL HSP60/10 CLIENTS
##
## Reframed pathology-coupling figure:
##   A. Paired inverse Braak/tau vs inverse CERAD/amyloid coupling
##   B. Late-stage decline vs inverse Braak/tau coupling
##   C. Distribution of Braak-minus-CERAD pathology-coupling bias
##
## Key interpretation:
##   Protein-level pathology metrics are covariate-adjusted upstream
##   and modeled separately. Braak/tau coupling is direct and is not
##   CERAD-residualized.
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

fig_id <- "main_fig3_all_clients"

fig3_all_dir <- file.path(cfg$output_root, fig_id)
fig3_all_plot_dir <- file.path(cfg$plot_dir, fig_id)
fig3_all_table_dir <- file.path(cfg$table_dir, fig_id)

dir.create(fig3_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig3_all_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig3_all_table_dir, recursive = TRUE, showWarnings = FALSE)

write_fig3_all_table <- function(x, name) {
  out <- file.path(fig3_all_table_dir, paste0(name, ".csv"))
  readr::write_csv(x, out)
  message("Wrote table: ", out)
  invisible(x)
}

save_fig3_all_plot <- function(plot, name, width = 15.5, height = 9.2, dpi = 600) {
  pdf_out <- file.path(fig3_all_plot_dir, paste0(name, ".pdf"))
  png_out <- file.path(fig3_all_plot_dir, paste0(name, ".png"))
  svg_out <- file.path(fig3_all_plot_dir, paste0(name, ".svg"))
  
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

safe_wilcox_paired <- function(x, y) {
  suppressWarnings(
    wilcox.test(
      x,
      y,
      paired = TRUE,
      exact = FALSE
    )$p.value
  )
}

safe_wilcox_unpaired <- function(x, y) {
  suppressWarnings(
    wilcox.test(
      x,
      y,
      paired = FALSE,
      exact = FALSE
    )$p.value
  )
}

theme_fig3 <- function(base_size = 11) {
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

theme_fig3_clean <- function(base_size = 11) {
  theme_fig3(base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "grey35", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1, color = "grey20"),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.key = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

############################################################
## 2. Build Figure 3 table
############################################################

require_objects(
  "all_hsp60_10_client_tbl",
  context = "12_make_main_figure_3_pathology_coupling.R"
)

require_columns(
  all_hsp60_10_client_tbl,
  c(
    "gene",
    "detected_in_protein",
    "protein_late_decline_magnitude",
    "braak_rho",
    "cerad_rho"
  ),
  "Figure 3 source all-client table"
)

fig3_tbl <- all_hsp60_10_client_tbl %>%
  mutate(
    gene = clean_gene(gene),
    detected_in_protein = as.logical(detected_in_protein),
    
    ## Late-stage protein decline:
    ## more positive = stronger MCI-to-AD protein decline.
    late_decline_value = suppressWarnings(as.numeric(protein_late_decline_magnitude)),
    
    ## Direct inverse Braak/tau association:
    ## more positive = stronger negative protein-Braak association.
    ## Prefer existing inverse_braak_magnitude if available.
    inverse_braak_value = case_when(
      "inverse_braak_magnitude" %in% names(.) ~ suppressWarnings(as.numeric(inverse_braak_magnitude)),
      is.finite(suppressWarnings(as.numeric(braak_rho))) &
        suppressWarnings(as.numeric(braak_rho)) < 0 ~ abs(suppressWarnings(as.numeric(braak_rho))),
      TRUE ~ 0
    ),
    
    ## Direct inverse CERAD/amyloid association:
    ## more positive = stronger negative protein-CERAD association.
    ## Prefer existing inverse_cerad_magnitude if available.
    inverse_cerad_value = case_when(
      "inverse_cerad_magnitude" %in% names(.) ~ suppressWarnings(as.numeric(inverse_cerad_magnitude)),
      is.finite(suppressWarnings(as.numeric(cerad_rho))) &
        suppressWarnings(as.numeric(cerad_rho)) < 0 ~ abs(suppressWarnings(as.numeric(cerad_rho))),
      TRUE ~ 0
    ),
    
    ## Plotting axes.
    braak_axis = inverse_braak_value,
    cerad_axis = inverse_cerad_value,
    braak_minus_cerad_inverse = braak_axis - cerad_axis
  ) %>%
  filter(
    detected_in_protein == TRUE,
    is.finite(late_decline_value),
    is.finite(braak_axis),
    is.finite(cerad_axis),
    is.finite(braak_minus_cerad_inverse)
  ) %>%
  arrange(desc(late_decline_value)) %>%
  mutate(
    late_decline_rank = row_number(),
    n_detected_clients = n(),
    late_decline_group = if_else(
      late_decline_rank <= ceiling(0.25 * n_detected_clients),
      "Top-quartile late decline",
      "Other detected Hsp60/10 clients"
    ),
    late_decline_group = factor(
      late_decline_group,
      levels = c(
        "Top-quartile late decline",
        "Other detected Hsp60/10 clients"
      )
    )
  )

write_fig3_all_table(fig3_tbl, "Fig3_all_clients_pathology_table")

cat("\nFigure 3 table dimensions:\n")
print(dim(fig3_tbl))

cat("\nFigure 3 late-decline group counts:\n")
print(table(fig3_tbl$late_decline_group))

############################################################
## 3. Statistics
############################################################

pathology_pair_tbl <- fig3_tbl %>%
  filter(
    is.finite(braak_axis),
    is.finite(cerad_axis),
    is.finite(braak_minus_cerad_inverse)
  )

pathology_axis_stats <- pathology_pair_tbl %>%
  summarise(
    n_genes = n(),
    median_inverse_braak = median(braak_axis, na.rm = TRUE),
    median_inverse_cerad = median(cerad_axis, na.rm = TRUE),
    mean_inverse_braak = mean(braak_axis, na.rm = TRUE),
    mean_inverse_cerad = mean(cerad_axis, na.rm = TRUE),
    median_braak_minus_cerad = median(braak_minus_cerad_inverse, na.rm = TRUE),
    mean_braak_minus_cerad = mean(braak_minus_cerad_inverse, na.rm = TRUE),
    n_braak_gt_cerad = sum(braak_minus_cerad_inverse > 0, na.rm = TRUE),
    n_equal = sum(braak_minus_cerad_inverse == 0, na.rm = TRUE),
    n_cerad_gt_braak = sum(braak_minus_cerad_inverse < 0, na.rm = TRUE),
    p_value = safe_wilcox_paired(braak_axis, cerad_axis),
    p_label = p_to_label(p_value),
    star = star_label(p_value)
  )

late_decline_braak_tbl <- fig3_tbl %>%
  filter(
    is.finite(late_decline_value),
    is.finite(braak_axis)
  )

late_decline_braak_cor <- suppressWarnings(
  cor.test(
    late_decline_braak_tbl$late_decline_value,
    late_decline_braak_tbl$braak_axis,
    method = "spearman",
    exact = FALSE
  )
)

late_decline_braak_stats <- tibble(
  n_genes = nrow(late_decline_braak_tbl),
  spearman_rho = unname(late_decline_braak_cor$estimate),
  p_value = late_decline_braak_cor$p.value,
  p_label = p_to_label(p_value),
  star = star_label(p_value)
)

braak_group_stats <- fig3_tbl %>%
  filter(is.finite(braak_axis)) %>%
  summarise(
    n_top = sum(late_decline_group == "Top-quartile late decline"),
    n_other = sum(late_decline_group == "Other detected Hsp60/10 clients"),
    median_top = median(braak_axis[late_decline_group == "Top-quartile late decline"], na.rm = TRUE),
    median_other = median(braak_axis[late_decline_group == "Other detected Hsp60/10 clients"], na.rm = TRUE),
    p_value = safe_wilcox_unpaired(
      braak_axis[late_decline_group == "Top-quartile late decline"],
      braak_axis[late_decline_group == "Other detected Hsp60/10 clients"]
    ),
    p_label = p_to_label(p_value),
    star = star_label(p_value)
  )

delta_summary_tbl <- pathology_pair_tbl %>%
  summarise(
    n_genes = n(),
    median_braak_minus_cerad = median(braak_minus_cerad_inverse, na.rm = TRUE),
    mean_braak_minus_cerad = mean(braak_minus_cerad_inverse, na.rm = TRUE),
    n_braak_gt_cerad = sum(braak_minus_cerad_inverse > 0, na.rm = TRUE),
    n_equal = sum(braak_minus_cerad_inverse == 0, na.rm = TRUE),
    n_cerad_gt_braak = sum(braak_minus_cerad_inverse < 0, na.rm = TRUE),
    p_value = pathology_axis_stats$p_value,
    p_label = pathology_axis_stats$p_label,
    star = pathology_axis_stats$star
  )

fig3_stats_tbl <- bind_rows(
  pathology_axis_stats %>%
    mutate(test = "paired_inverse_braak_vs_inverse_cerad"),
  late_decline_braak_stats %>%
    mutate(test = "late_decline_vs_inverse_braak_spearman"),
  braak_group_stats %>%
    mutate(test = "inverse_braak_top_quartile_vs_other"),
  delta_summary_tbl %>%
    mutate(test = "braak_minus_cerad_delta_summary")
)

write_fig3_all_table(pathology_axis_stats, "Fig3A_inverse_braak_vs_inverse_cerad_stats")
write_fig3_all_table(late_decline_braak_stats, "Fig3B_late_decline_vs_inverse_braak_correlation_stats")
write_fig3_all_table(braak_group_stats, "Fig3_inverse_braak_group_stats_audit")
write_fig3_all_table(delta_summary_tbl, "Fig3C_braak_minus_cerad_delta_summary")
write_fig3_all_table(fig3_stats_tbl, "Fig3_all_clients_stats_combined")

print(fig3_stats_tbl)

############################################################
## 4. Labels and colors
############################################################

late_decline_colors <- c(
  "Top-quartile late decline" = "#9E4A4A",
  "Other detected Hsp60/10 clients" = "grey75"
)

late_decline_edge_colors <- c(
  "Top-quartile late decline" = "#9E4A4A",
  "Other detected Hsp60/10 clients" = "grey55"
)

## Panel B labels: requested high-priority genes only.
fig3_label_genes_B <- c(
  "DAP3", "MRPS35", "MRPS22", "MRPS16", "MRPS9", "MRPL48",
  "MRPS23", "NDUFAF7", "PTCD3", "LRPPRC", "MRPS33"
)

fig3_label_tbl_B <- fig3_tbl %>%
  filter(gene %in% fig3_label_genes_B) %>%
  mutate(
    label_gene_B = gene,
    label_x = case_when(
      gene == "MRPS33" ~ 0.096,
      gene == "LRPPRC" ~ 0.086,
      gene == "DAP3" ~ 0.104,
      gene == "MRPS35" ~ 0.098,
      gene == "PTCD3" ~ 0.101,
      gene == "NDUFAF7" ~ 0.096,
      gene == "MRPS22" ~ 0.086,
      gene == "MRPS16" ~ 0.058,
      gene == "MRPS9" ~ 0.038,
      gene == "MRPL48" ~ 0.058,
      gene == "MRPS23" ~ 0.058,
      TRUE ~ late_decline_value
    ),
    label_y = case_when(
      gene == "MRPS33" ~ 0.415,
      gene == "LRPPRC" ~ 0.382,
      gene == "DAP3" ~ 0.338,
      gene == "MRPS35" ~ 0.298,
      gene == "PTCD3" ~ 0.258,
      gene == "NDUFAF7" ~ 0.096,
      gene == "MRPS22" ~ 0.158,
      gene == "MRPS16" ~ 0.336,
      gene == "MRPS9" ~ 0.282,
      gene == "MRPL48" ~ 0.202,
      gene == "MRPS23" ~ 0.132,
      TRUE ~ braak_axis
    ),
    label_hjust = if_else(label_x < late_decline_value, 1, 0)
  ) %>%
  arrange(match(gene, fig3_label_genes_B))

############################################################
## 5. Panel A — paired summary with all gene-level points
############################################################

panelA_long <- fig3_tbl %>%
  select(gene, late_decline_group, braak_axis, cerad_axis) %>%
  pivot_longer(
    cols = c(cerad_axis, braak_axis),
    names_to = "pathology_axis",
    values_to = "inverse_association"
  ) %>%
  mutate(
    pathology_axis = recode(
      pathology_axis,
      cerad_axis = "CERAD/amyloid",
      braak_axis = "Braak/tau"
    ),
    pathology_axis = factor(
      pathology_axis,
      levels = c("CERAD/amyloid", "Braak/tau")
    )
  )

panelA_label <- paste0(
  "Paired Wilcoxon\n",
  pathology_axis_stats$p_label,
  "\nn = ",
  pathology_axis_stats$n_genes
)

pA <- ggplot(
  panelA_long,
  aes(x = pathology_axis, y = inverse_association)
) +
  geom_violin(
    width = 0.72,
    fill = "grey92",
    color = "grey70",
    linewidth = 0.30,
    alpha = 0.45,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    color = "grey30",
    linewidth = 0.40
  ) +
  geom_point(
    aes(fill = late_decline_group),
    position = position_jitter(width = 0.05, height = 0, seed = 1),
    shape = 21,
    size = 1.7,
    stroke = 0.18,
    color = "grey35",
    alpha = 0.65
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = panelA_label,
    hjust = -0.02,
    vjust = 1.10,
    size = 2.7,
    linewidth = 0.25,
    fill = "white",
    color = "grey20"
  ) +
  scale_fill_manual(
    values = late_decline_colors,
    name = "Late-decline group",
    labels = c("Top-quartile late decline", "Other clients")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.18))
  ) +
  labs(
    title = "A. Inverse pathology coupling is stronger for Braak/tau than CERAD/amyloid",
    subtitle = "Covariate-adjusted Braak/tau and CERAD/amyloid associations were modeled separately for each client.",
    x = NULL,
    y = "Inverse pathology association"
  ) +
  theme_fig3_clean(11) +
  theme(
    legend.position = "top"
  )

############################################################
## 6. Panel B — collapse vs inverse Braak/tau coupling
############################################################

cor_label_B <- paste0(
  "Spearman rho = ",
  signif(late_decline_braak_stats$spearman_rho, 2),
  "\n",
  late_decline_braak_stats$p_label,
  "\nn = ",
  late_decline_braak_stats$n_genes
)

pB <- ggplot(
  late_decline_braak_tbl,
  aes(x = late_decline_value, y = braak_axis)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "grey35",
    fill = "grey82",
    linewidth = 0.75
  ) +
  geom_point(
    aes(fill = late_decline_group),
    shape = 21,
    color = "grey35",
    alpha = 0.78,
    size = 2.15,
    stroke = 0.22
  ) +
  geom_segment(
    data = fig3_label_tbl_B,
    aes(x = late_decline_value, y = braak_axis, xend = label_x, yend = label_y),
    inherit.aes = FALSE,
    color = "#7F1D1D",
    linewidth = 0.18,
    alpha = 0.95
  ) +
  geom_text(
    data = fig3_label_tbl_B,
    aes(x = label_x, y = label_y, label = label_gene_B, hjust = label_hjust),
    inherit.aes = FALSE,
    size = 2.55,
    fontface = "bold",
    color = "#7F1D1D",
    show.legend = FALSE
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = cor_label_B,
    hjust = -0.02,
    vjust = 1.10,
    size = 2.7,
    linewidth = 0.25,
    fill = "white",
    color = "grey20"
  ) +
  scale_fill_manual(values = late_decline_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.34))) +
  labs(
    title = "B. Late-stage decline is continuously associated with inverse Braak/tau coupling",
    subtitle = "Greater late-stage protein decline aligns with stronger direct inverse Braak association.",
    x = "Late-stage protein decline",
    y = "Inverse Braak/tau association"
  ) +
  theme_fig3_clean(11) +
  theme(
    legend.position = "none"
  )

############################################################
## 7. Panel C — distribution of Braak-minus-CERAD pathology bias
############################################################

delta_label <- paste0(
  "Median Braak - CERAD = ",
  signif(delta_summary_tbl$median_braak_minus_cerad, 2),
  "\n",
  delta_summary_tbl$p_label,
  "\nn = ",
  delta_summary_tbl$n_genes
)

pC <- ggplot(
  pathology_pair_tbl,
  aes(x = braak_minus_cerad_inverse)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.45
  ) +
  geom_histogram(
    bins = 38,
    fill = "grey78",
    color = "white",
    linewidth = 0.20
  ) +
  geom_rug(
    aes(color = late_decline_group),
    alpha = 0.45,
    sides = "b",
    linewidth = 0.35,
    show.legend = TRUE
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    hjust = 1.02,
    vjust = 1.10,
    label = delta_label,
    size = 3.0,
    linewidth = 0.22,
    fill = "white",
    color = "grey20"
  ) +
  scale_color_manual(
    values = late_decline_edge_colors,
    name = "Late-decline group",
    labels = c("Top-quartile late decline", "Other clients")
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  labs(
    title = "C. Braak/tau coupling exceeds CERAD/amyloid coupling across the client network",
    subtitle = "Positive values indicate stronger inverse Braak/tau coupling than inverse CERAD/amyloid coupling for the same client.",
    x = "Inverse Braak/tau association - inverse CERAD/amyloid association",
    y = "Number of Hsp60/10 clients"
  ) +
  theme_fig3_clean(11) +
  theme(
    legend.position = "top"
  )

############################################################
## 8. Assemble and save Figure 3
############################################################

fig3_all_clients <- (pA | pB) / pC +
  plot_layout(
    heights = c(1.00, 0.92),
    widths = c(1.0, 1.0),
    guides = "collect"
  ) +
  plot_annotation(
    title = "Hsp60/10 client vulnerability is preferentially coupled to Braak/tau pathology",
    subtitle = "Covariate-adjusted Braak/tau and CERAD/amyloid associations were modeled separately across detected Hsp60/10 clients.",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey35"),
      plot.margin = margin(8, 12, 8, 12)
    )
  ) &
  theme(
    legend.position = "top"
  )

save_fig3_all_plot(
  fig3_all_clients,
  "Main_Fig3_all_clients_adjusted_pathology_coupling_REFRAMED",
  width = 15.5,
  height = 9.2
)

fig3_all_clients

message("Loaded 12_make_main_figure_3_pathology_coupling.R")
