############################################################
## 12_make_main_figure_3_pathology_coupling.R
## MAIN FIGURE 3 — ALL HSP60/10 CLIENTS
##
## Corrected pathology-coupling figure:
##   A. Paired signed Braak/tau vs CERAD/amyloid pathology associations
##   B. Late-stage decline vs signed Braak/tau-aligned association
##   C. Late-stage decline vs signed CERAD/amyloid-aligned association
##
## Positive pathology-aligned beta means lower protein abundance
## with worse pathology for BOTH endpoints:
##   Braak: -braak_beta
##   CERAD: +cerad_beta
##
## Braak and CERAD are modeled separately using the same protein
## covariate structure. No zero-truncation is used in this figure.
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
    "braak_pathology_aligned_beta",
    "cerad_pathology_aligned_beta"
  ),
  "Figure 3 source all-client table"
)

fig3_tbl <- all_hsp60_10_client_tbl %>%
  mutate(
    gene = clean_gene(gene),

    detected_in_protein =
      as.logical(detected_in_protein),

    ## Late-stage protein decline:
    ## more positive = stronger MCI-to-AD protein decline.
    late_decline_value =
      suppressWarnings(
        as.numeric(
          protein_late_decline_magnitude
        )
      ),

    ## FULL SIGNED pathology-aligned coefficients.
    ##
    ## Positive values have the same biological direction
    ## for both endpoints:
    ## lower protein abundance with worse pathology.
    braak_aligned_beta =
      suppressWarnings(
        as.numeric(
          braak_pathology_aligned_beta
        )
      ),

    cerad_aligned_beta =
      suppressWarnings(
        as.numeric(
          cerad_pathology_aligned_beta
        )
      ),

    braak_minus_cerad_beta =
      braak_aligned_beta -
        cerad_aligned_beta
  ) %>%
  filter(
    detected_in_protein == TRUE,
    is.finite(late_decline_value),
    is.finite(braak_aligned_beta),
    is.finite(cerad_aligned_beta)
  ) %>%
  arrange(
    desc(late_decline_value)
  ) %>%
  mutate(
    late_decline_rank =
      row_number(),

    n_detected_clients =
      n(),

    late_decline_group =
      if_else(
        late_decline_rank <=
          ceiling(
            0.25 *
              n_detected_clients
          ),
        "Top-quartile late decline",
        "Other detected Hsp60/10 clients"
      ),

    late_decline_group =
      factor(
        late_decline_group,
        levels = c(
          "Top-quartile late decline",
          "Other detected Hsp60/10 clients"
        )
      )
  )

write_fig3_all_table(
  fig3_tbl,
  "Fig3_all_clients_pathology_table"
)

cat("\nFigure 3 table dimensions:\n")
print(dim(fig3_tbl))

cat("\nFigure 3 late-decline group counts:\n")
print(
  table(
    fig3_tbl$late_decline_group
  )
)

############################################################
## 3. Statistics
############################################################

pathology_pair_tbl <- fig3_tbl %>%
  filter(
    is.finite(braak_aligned_beta),
    is.finite(cerad_aligned_beta),
    is.finite(braak_minus_cerad_beta)
  )

pathology_axis_stats <- pathology_pair_tbl %>%
  summarise(
    n_genes = n(),

    median_braak_aligned_beta =
      median(
        braak_aligned_beta,
        na.rm = TRUE
      ),

    median_cerad_aligned_beta =
      median(
        cerad_aligned_beta,
        na.rm = TRUE
      ),

    mean_braak_aligned_beta =
      mean(
        braak_aligned_beta,
        na.rm = TRUE
      ),

    mean_cerad_aligned_beta =
      mean(
        cerad_aligned_beta,
        na.rm = TRUE
      ),

    median_braak_minus_cerad =
      median(
        braak_minus_cerad_beta,
        na.rm = TRUE
      ),

    mean_braak_minus_cerad =
      mean(
        braak_minus_cerad_beta,
        na.rm = TRUE
      ),

    n_braak_positive =
      sum(
        braak_aligned_beta > 0,
        na.rm = TRUE
      ),

    n_braak_negative =
      sum(
        braak_aligned_beta < 0,
        na.rm = TRUE
      ),

    n_cerad_positive =
      sum(
        cerad_aligned_beta > 0,
        na.rm = TRUE
      ),

    n_cerad_negative =
      sum(
        cerad_aligned_beta < 0,
        na.rm = TRUE
      ),

    n_braak_gt_cerad =
      sum(
        braak_minus_cerad_beta > 0,
        na.rm = TRUE
      ),

    n_equal =
      sum(
        braak_minus_cerad_beta == 0,
        na.rm = TRUE
      ),

    n_cerad_gt_braak =
      sum(
        braak_minus_cerad_beta < 0,
        na.rm = TRUE
      ),

    p_value =
      safe_wilcox_paired(
        braak_aligned_beta,
        cerad_aligned_beta
      ),

    p_label =
      p_to_label(p_value),

    star =
      star_label(p_value)
  )

############################################################
## Late decline vs Braak
############################################################

late_decline_braak_tbl <- fig3_tbl %>%
  filter(
    is.finite(late_decline_value),
    is.finite(braak_aligned_beta)
  )

late_decline_braak_cor <- suppressWarnings(
  cor.test(
    late_decline_braak_tbl$late_decline_value,
    late_decline_braak_tbl$braak_aligned_beta,
    method = "spearman",
    exact = FALSE
  )
)

late_decline_braak_stats <- tibble(
  n_genes =
    nrow(late_decline_braak_tbl),

  spearman_rho =
    unname(
      late_decline_braak_cor$estimate
    ),

  p_value =
    late_decline_braak_cor$p.value,

  p_label =
    p_to_label(p_value),

  star =
    star_label(p_value)
)

############################################################
## Late decline vs CERAD
############################################################

late_decline_cerad_tbl <- fig3_tbl %>%
  filter(
    is.finite(late_decline_value),
    is.finite(cerad_aligned_beta)
  )

late_decline_cerad_cor <- suppressWarnings(
  cor.test(
    late_decline_cerad_tbl$late_decline_value,
    late_decline_cerad_tbl$cerad_aligned_beta,
    method = "spearman",
    exact = FALSE
  )
)

late_decline_cerad_stats <- tibble(
  n_genes =
    nrow(late_decline_cerad_tbl),

  spearman_rho =
    unname(
      late_decline_cerad_cor$estimate
    ),

  p_value =
    late_decline_cerad_cor$p.value,

  p_label =
    p_to_label(p_value),

  star =
    star_label(p_value)
)

############################################################
## Top-quartile late-decline descriptive comparisons
############################################################

braak_group_stats <- fig3_tbl %>%
  summarise(
    n_top =
      sum(
        late_decline_group ==
          "Top-quartile late decline"
      ),

    n_other =
      sum(
        late_decline_group ==
          "Other detected Hsp60/10 clients"
      ),

    median_top =
      median(
        braak_aligned_beta[
          late_decline_group ==
            "Top-quartile late decline"
        ],
        na.rm = TRUE
      ),

    median_other =
      median(
        braak_aligned_beta[
          late_decline_group ==
            "Other detected Hsp60/10 clients"
        ],
        na.rm = TRUE
      ),

    p_value =
      safe_wilcox_unpaired(
        braak_aligned_beta[
          late_decline_group ==
            "Top-quartile late decline"
        ],
        braak_aligned_beta[
          late_decline_group ==
            "Other detected Hsp60/10 clients"
        ]
      ),

    p_label =
      p_to_label(p_value),

    star =
      star_label(p_value)
  )

cerad_group_stats <- fig3_tbl %>%
  summarise(
    n_top =
      sum(
        late_decline_group ==
          "Top-quartile late decline"
      ),

    n_other =
      sum(
        late_decline_group ==
          "Other detected Hsp60/10 clients"
      ),

    median_top =
      median(
        cerad_aligned_beta[
          late_decline_group ==
            "Top-quartile late decline"
        ],
        na.rm = TRUE
      ),

    median_other =
      median(
        cerad_aligned_beta[
          late_decline_group ==
            "Other detected Hsp60/10 clients"
        ],
        na.rm = TRUE
      ),

    p_value =
      safe_wilcox_unpaired(
        cerad_aligned_beta[
          late_decline_group ==
            "Top-quartile late decline"
        ],
        cerad_aligned_beta[
          late_decline_group ==
            "Other detected Hsp60/10 clients"
        ]
      ),

    p_label =
      p_to_label(p_value),

    star =
      star_label(p_value)
  )

fig3_stats_tbl <- bind_rows(
  pathology_axis_stats %>%
    mutate(
      test =
        "paired_signed_braak_vs_cerad"
    ),

  late_decline_braak_stats %>%
    mutate(
      test =
        "late_decline_vs_signed_braak_spearman"
    ),

  late_decline_cerad_stats %>%
    mutate(
      test =
        "late_decline_vs_signed_cerad_spearman"
    ),

  braak_group_stats %>%
    mutate(
      test =
        "signed_braak_top_late_decline_vs_other"
    ),

  cerad_group_stats %>%
    mutate(
      test =
        "signed_cerad_top_late_decline_vs_other"
    )
)

write_fig3_all_table(
  pathology_axis_stats,
  "Fig3A_signed_braak_vs_cerad_stats"
)

write_fig3_all_table(
  late_decline_braak_stats,
  "Fig3B_late_decline_vs_signed_braak_correlation_stats"
)

write_fig3_all_table(
  late_decline_cerad_stats,
  "Fig3C_late_decline_vs_signed_cerad_correlation_stats"
)

write_fig3_all_table(
  braak_group_stats,
  "Fig3_signed_braak_group_stats_audit"
)

write_fig3_all_table(
  cerad_group_stats,
  "Fig3_signed_cerad_group_stats_audit"
)

write_fig3_all_table(
  fig3_stats_tbl,
  "Fig3_all_clients_stats_combined"
)

cat("\n===== FIGURE 3 PATHOLOGY STATISTICS =====\n")
print(
  pathology_axis_stats,
  width = Inf
)

cat("\n===== LATE DECLINE VS BRAAK =====\n")
print(
  late_decline_braak_stats,
  width = Inf
)

cat("\n===== LATE DECLINE VS CERAD =====\n")
print(
  late_decline_cerad_stats,
  width = Inf
)

cat("\n===== BRAAK TOP-LATE-DECLINE GROUP =====\n")
print(
  braak_group_stats,
  width = Inf
)

cat("\n===== CERAD TOP-LATE-DECLINE GROUP =====\n")
print(
  cerad_group_stats,
  width = Inf
)

############################################################
## 4. Labels, colors, and shared pathology scale
############################################################

late_decline_colors <- c(
  "Top-quartile late decline" =
    "#9E4A4A",

  "Other detected Hsp60/10 clients" =
    "grey75"
)

fig3_label_genes <- c(
  "DAP3",
  "LRPPRC",
  "PDHA1",
  "TUFM",
  "MRPL48",
  "NDUFAF7",
  "PTCD3",
  "MRPS35"
)

fig3_label_tbl <- fig3_tbl %>%
  filter(
    gene %in% fig3_label_genes
  )

pathology_range <- range(
  c(
    fig3_tbl$braak_aligned_beta,
    fig3_tbl$cerad_aligned_beta
  ),
  na.rm = TRUE
)

pathology_pad <-
  diff(pathology_range) * 0.08

if (
  !is.finite(pathology_pad) ||
  pathology_pad <= 0
) {
  pathology_pad <- 0.05
}

pathology_y_limits <- c(
  pathology_range[[1]] -
    pathology_pad,
  pathology_range[[2]] +
    pathology_pad
)

############################################################
## 5. Panel A — paired signed pathology associations
############################################################

panelA_long <- fig3_tbl %>%
  select(
    gene,
    late_decline_group,
    braak_aligned_beta,
    cerad_aligned_beta
  ) %>%
  pivot_longer(
    cols = c(
      braak_aligned_beta,
      cerad_aligned_beta
    ),
    names_to =
      "pathology_axis",
    values_to =
      "pathology_aligned_beta"
  ) %>%
  mutate(
    pathology_axis =
      recode(
        pathology_axis,
        braak_aligned_beta =
          "Braak/tau",
        cerad_aligned_beta =
          "CERAD/amyloid"
      ),

    pathology_axis =
      factor(
        pathology_axis,
        levels = c(
          "Braak/tau",
          "CERAD/amyloid"
        )
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
  aes(
    x = pathology_axis,
    y = pathology_aligned_beta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.4
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
    aes(
      fill = late_decline_group
    ),
    position =
      position_jitter(
        width = 0.05,
        height = 0,
        seed = 1
      ),
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
    size = 2.35,
    linewidth = 0.25,
    fill = "white",
    color = "grey20"
  ) +
  scale_fill_manual(
    values =
      late_decline_colors,
    name =
      "Late-decline group",
    labels = c(
      "Top-quartile late decline",
      "Other clients"
    )
  ) +
  coord_cartesian(
    ylim =
      pathology_y_limits,
    clip = "off"
  ) +
  labs(
    title =
      "A. Pathology-aligned protein associations span both AD neuropathologic measures",

    subtitle =
      "Positive beta indicates lower protein abundance with worse pathology.",

    x = NULL,

    y =
      "Pathology-aligned beta"
  ) +
  theme_fig3_clean(11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 9.2)
  )

############################################################
## 6. Panel B — late decline vs signed Braak/tau
############################################################

cor_label_B <- paste0(
  "Spearman \u03c1 = ",
  signif(
    late_decline_braak_stats$spearman_rho,
    2
  ),
  "\n",
  late_decline_braak_stats$p_label,
  "\nn = ",
  late_decline_braak_stats$n_genes
)

pB <- ggplot(
  late_decline_braak_tbl,
  aes(
    x = late_decline_value,
    y = braak_aligned_beta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey35",
    fill = "grey82",
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      fill =
        late_decline_group
    ),
    shape = 21,
    color = "grey35",
    alpha = 0.78,
    size = 2.15,
    stroke = 0.22
  ) +
  ggrepel::geom_text_repel(
    data =
      fig3_label_tbl,
    aes(
      x =
        late_decline_value,
      y =
        braak_aligned_beta,
      label =
        gene
    ),
    inherit.aes = FALSE,
    seed = 12,
    size = 2.15,
    fontface = "bold",
    color = "#7F1D1D",
    min.segment.length = 0,
    segment.size = 0.18,
    box.padding = 0.22,
    point.padding = 0.12,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = cor_label_B,
    hjust = -0.02,
    vjust = 1.10,
    size = 2.35,
    linewidth = 0.25,
    fill = "white",
    color = "grey20"
  ) +
  scale_fill_manual(
    values =
      late_decline_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    expand =
      expansion(
        mult = c(
          0.08,
          0.12
        )
      )
  ) +
  coord_cartesian(
    ylim =
      pathology_y_limits,
    clip = "off"
  ) +
  labs(
    title =
      "B. Late-stage decline is associated with Braak/tau-aligned protein loss",

    subtitle =
      "Higher pathology-aligned beta indicates lower protein abundance with worse Braak stage.",

    x =
      "Late-stage protein decline",

    y =
      "Braak/tau-aligned beta"
  ) +
  theme_fig3_clean(11) +
  theme(
    legend.position = "none"
  )

############################################################
## 7. Panel C — late decline vs signed CERAD/amyloid
############################################################

cor_label_C <- paste0(
  "Spearman \u03c1 = ",
  signif(
    late_decline_cerad_stats$spearman_rho,
    2
  ),
  "\n",
  late_decline_cerad_stats$p_label,
  "\nn = ",
  late_decline_cerad_stats$n_genes
)

pC <- ggplot(
  late_decline_cerad_tbl,
  aes(
    x = late_decline_value,
    y = cerad_aligned_beta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey35",
    fill = "grey82",
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      fill =
        late_decline_group
    ),
    shape = 21,
    color = "grey35",
    alpha = 0.78,
    size = 2.15,
    stroke = 0.22
  ) +
  ggrepel::geom_text_repel(
    data =
      fig3_label_tbl,
    aes(
      x =
        late_decline_value,
      y =
        cerad_aligned_beta,
      label =
        gene
    ),
    inherit.aes = FALSE,
    seed = 13,
    size = 2.15,
    fontface = "bold",
    color = "#7F1D1D",
    min.segment.length = 0,
    segment.size = 0.18,
    box.padding = 0.22,
    point.padding = 0.12,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = cor_label_C,
    hjust = -0.02,
    vjust = 1.10,
    size = 2.35,
    linewidth = 0.25,
    fill = "white",
    color = "grey20"
  ) +
  scale_fill_manual(
    values =
      late_decline_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    expand =
      expansion(
        mult = c(
          0.08,
          0.12
        )
      )
  ) +
  coord_cartesian(
    ylim =
      pathology_y_limits,
    clip = "off"
  ) +
  labs(
    title =
      "C. Late-stage decline is associated with CERAD/amyloid-aligned protein loss",

    subtitle =
      "Because lower CERAD scores indicate worse pathology, positive beta is the pathology-aligned decline direction.",

    x =
      "Late-stage protein decline",

    y =
      "CERAD/amyloid-aligned beta"
  ) +
  theme_fig3_clean(11) +
  theme(
    legend.position = "none"
  )

############################################################
## 8. Assemble and save Figure 3
############################################################

fig3_all_clients <-
  (pA | pB) /
    pC +
    plot_layout(
      heights = c(
        1.00,
        0.95
      ),
      widths = c(
        1.0,
        1.0
      ),
      guides = "collect"
    ) +
    plot_annotation(
      title =
        "Late-stage Hsp60/10 client decline is associated with both tau and amyloid pathology",

      subtitle =
        paste0(
          "Braak/tau and CERAD/amyloid were modeled separately. ",
          "Positive pathology-aligned beta indicates lower protein abundance with worse pathology."
        ),

      theme =
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),

          plot.subtitle = element_text(hjust = 0.5, size = 10.5, color = "grey35"),

          plot.margin =
            margin(
              8,
              12,
              8,
              12
            )
        )
    ) &
    theme(
      legend.position = "top"
    )

save_fig3_all_plot(
  fig3_all_clients,
  "Main_Fig3_all_clients_signed_braak_cerad_pathology_coupling",
  width = 7.1,
  height = 5.6
)

fig3_all_clients

message(
  "Loaded 12_make_main_figure_3_pathology_coupling.R"
)
