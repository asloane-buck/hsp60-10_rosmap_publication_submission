############################################################
## 10_make_main_figure_1_pathway_remodeling.R
## Literal migration of the final old plotting block.
## Inputs are supplied by adjusted upstream scripts.
############################################################

## MAIN FIGURE 1 — ALL HSP60/10 CLIENTS
## Clean version with:
## - all 320+ Hsp60/10 clients
## - broad/non-overlapping comparator pathways
## - RNA vs protein significance in Panel A
## - category-colored facet strips
## - cleaned Panel B labels
## - bootstrap CIs/significance for Panel B
## - no legend "a" artifact
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
})

if (!requireNamespace("ggtext", quietly = TRUE)) {
  stop("Please install ggtext first: install.packages('ggtext')", call. = FALSE)
}

require_objects(
  c("cfg", "prot_mat", "rna_mat", "prot_scores", "rna_scores"),
  context = "10_make_main_figure_1_pathway_remodeling.R"
)

############################################################
## 0. Output setup
############################################################

fig_id <- "main_fig1_all_clients"

## Keep Figure 1 outputs inside the covariate-adjusted modular output root.
fig1_all_dir <- file.path(cfg$output_root, fig_id)
fig1_all_plot_dir <- file.path(cfg$plot_dir, fig_id)
fig1_all_table_dir <- file.path(cfg$table_dir, fig_id)

dir.create(fig1_all_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig1_all_table_dir, recursive = TRUE, showWarnings = FALSE)

write_fig1_all_table <- function(x, name) {
  out <- file.path(fig1_all_table_dir, paste0(name, ".csv"))
  readr::write_csv(x, out)
  message("Wrote table: ", out)
  invisible(x)
}

save_fig1_all_plot <- function(plot, name, width = 18, height = 12, dpi = 600) {
  pdf_out <- file.path(fig1_all_plot_dir, paste0(name, ".pdf"))
  png_out <- file.path(fig1_all_plot_dir, paste0(name, ".png"))
  svg_out <- file.path(fig1_all_plot_dir, paste0(name, ".svg"))

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
## 1. Load clean pathway sets
############################################################

pathway_set_file <- first_existing(c(
  file.path(
    cfg$pathway_dir,
    "pathway_gene_sets_all_clients_NO_Hsp60_10_overlap.rds"
  )
), "clean Figure 1 pathway set RDS")

pathway_gene_sets_all <- readRDS(pathway_set_file)

############################################################
## 2. Figure 1 pathway selection
############################################################

figure1_pathways_all <- c(
  "Hsp60_10_all_clients",
  "Broad_MitoCarta_non_Hsp60_10",
  "Mito_translation_non_Hsp60_10",
  "OXPHOS_ETC_non_Hsp60_10",
  "TCA_pyruvate_metabolism_non_Hsp60_10",
  "FAO_metabolism_non_Hsp60_10",
  "Mito_protein_quality_control_non_Hsp60_10",
  "Proteasome_core_non_Hsp60_10",
  "Lysosome_core_non_Hsp60_10",
  "UPRmt_core_non_Hsp60_10"
)

figure1_pathways_all <- intersect(
  figure1_pathways_all,
  names(pathway_gene_sets_all)
)

if (length(figure1_pathways_all) < 3) {
  stop(
    "Too few Figure 1 pathways found. Found: ",
    paste(figure1_pathways_all, collapse = ", "),
    call. = FALSE
  )
}

fig1_pathway_sets <- pathway_gene_sets_all[figure1_pathways_all]

############################################################
## 3. Labels and colors
############################################################

fig1_panelA_plain_labels <- c(
  Hsp60_10_all_clients = "Hsp60/10 clients",
  Broad_MitoCarta_non_Hsp60_10 = "Broad mitochondrial\nnon-Hsp60/10",
  Mito_translation_non_Hsp60_10 = "Mito translation\nnon-Hsp60/10",
  OXPHOS_ETC_non_Hsp60_10 = "OXPHOS/ETC\nnon-Hsp60/10",
  TCA_pyruvate_metabolism_non_Hsp60_10 = "TCA/pyruvate\nnon-Hsp60/10",
  FAO_metabolism_non_Hsp60_10 = "FAO/metabolism\nnon-Hsp60/10",
  Mito_protein_quality_control_non_Hsp60_10 = "Mito protein QC\nnon-Hsp60/10",
  Proteasome_core_non_Hsp60_10 = "Proteasome\ncore",
  Lysosome_core_non_Hsp60_10 = "Lysosome\ncore",
  UPRmt_core_non_Hsp60_10 = "UPRmt\ncore"
)

fig1_panelB_labels <- c(
  Hsp60_10_all_clients = "Hsp60/10 clients",
  Broad_MitoCarta_non_Hsp60_10 = "Broad mito",
  Mito_translation_non_Hsp60_10 = "Mito translation",
  OXPHOS_ETC_non_Hsp60_10 = "OXPHOS/ETC",
  TCA_pyruvate_metabolism_non_Hsp60_10 = "TCA/pyruvate",
  FAO_metabolism_non_Hsp60_10 = "FAO/metabolism",
  Mito_protein_quality_control_non_Hsp60_10 = "Mito protein QC",
  Proteasome_core_non_Hsp60_10 = "Proteasome",
  Lysosome_core_non_Hsp60_10 = "Lysosome",
  UPRmt_core_non_Hsp60_10 = "UPRmt"
)

fig1_pathway_colors <- c(
  Hsp60_10_all_clients = "#B23A48",
  Broad_MitoCarta_non_Hsp60_10 = "#4C78A8",
  Mito_translation_non_Hsp60_10 = "#7A3E9D",
  OXPHOS_ETC_non_Hsp60_10 = "#2E86DE",
  TCA_pyruvate_metabolism_non_Hsp60_10 = "#2CA02C",
  FAO_metabolism_non_Hsp60_10 = "#E67E22",
  Mito_protein_quality_control_non_Hsp60_10 = "#A65E2E",
  Proteasome_core_non_Hsp60_10 = "#7D3C98",
  Lysosome_core_non_Hsp60_10 = "#1E8449",
  UPRmt_core_non_Hsp60_10 = "#F39C12"
)

modality_colors <- c(
  "Protein" = "#8E1B1B",
  "RNA" = "#1B4F9C",
  "Protein vs RNA" = "black"
)

fig1_strip_labels <- c(
  Hsp60_10_all_clients =
    "<span style='color:#B23A48'><b>Hsp60/10<br>clients</b></span>",

  Broad_MitoCarta_non_Hsp60_10 =
    "<span style='color:#4C78A8'><b>Broad<br>mitochondrial<br>non-Hsp60/10</b></span>",

  Mito_translation_non_Hsp60_10 =
    "<span style='color:#7A3E9D'>Mito<br>translation<br>non-Hsp60/10</span>",

  OXPHOS_ETC_non_Hsp60_10 =
    "<span style='color:#2E86DE'>OXPHOS/ETC<br>non-Hsp60/10</span>",

  TCA_pyruvate_metabolism_non_Hsp60_10 =
    "<span style='color:#2CA02C'>TCA/<br>pyruvate<br>non-Hsp60/10</span>",

  FAO_metabolism_non_Hsp60_10 =
    "<span style='color:#E67E22'>FAO/<br>metabolism<br>non-Hsp60/10</span>",

  Mito_protein_quality_control_non_Hsp60_10 =
    "<span style='color:#A65E2E'>Mito protein<br>QC<br>non-Hsp60/10</span>",

  Proteasome_core_non_Hsp60_10 =
    "<span style='color:#7D3C98'>Proteasome<br>core</span>",

  Lysosome_core_non_Hsp60_10 =
    "<span style='color:#1E8449'>Lysosome<br>core</span>",

  UPRmt_core_non_Hsp60_10 =
    "<span style='color:#F39C12'>UPRmt<br>core</span>"
)

missing_fig1_label_paths <- setdiff(figure1_pathways_all, names(fig1_panelA_plain_labels))
if (length(missing_fig1_label_paths) > 0) {
  stop(
    "Figure 1 pathway label map is missing expected pathway label(s): ",
    paste(missing_fig1_label_paths, collapse = ", "),
    call. = FALSE
  )
}

############################################################
## 4. Coverage audit
############################################################

fig1_coverage_protein <- summarize_pathway_coverage(
  prot_mat,
  fig1_pathway_sets
) %>%
  mutate(Modality = "Protein")

fig1_coverage_rna <- summarize_pathway_coverage(
  rna_mat,
  fig1_pathway_sets
) %>%
  mutate(Modality = "RNA")

fig1_coverage <- bind_rows(fig1_coverage_protein, fig1_coverage_rna) %>%
  mutate(
    PathwayLabel = fig1_panelA_plain_labels[Pathway]
  )

write_fig1_all_table(fig1_coverage, "Fig1_all_clients_pathway_coverage")
print(fig1_coverage)

############################################################
## 5. Compute pathway scores
############################################################

prot_scores_all_clients <- compute_pathway_scores(
  prot_mat,
  fig1_pathway_sets,
  method = "mean_z",
  sample_col = "SampleID"
)

rna_scores_all_clients <- compute_pathway_scores(
  rna_mat,
  fig1_pathway_sets,
  method = "mean_z",
  sample_col = "sample_id"
)

require_columns(prot_scores, c("SampleID", "EmoryStrictDx.2019"), "Figure 1 protein stage metadata")
prot_stage_meta <- prot_scores %>%
  select(SampleID, EmoryStrictDx.2019) %>%
  distinct()

require_columns(rna_scores, c("sample_id", "diagnosis_stage"), "Figure 1 RNA stage metadata")
rna_stage_meta <- rna_scores %>%
  select(sample_id, diagnosis_stage) %>%
  distinct()

prot_scores_all_clients <- prot_scores_all_clients %>%
  checked_left_join(prot_stage_meta, by = "SampleID", label = "Figure 1 protein pathway scores to stage metadata")

rna_scores_all_clients <- rna_scores_all_clients %>%
  checked_left_join(rna_stage_meta, by = "sample_id", label = "Figure 1 RNA pathway scores to stage metadata")

############################################################
## 6. Long format and summaries
############################################################

overlay_long_all_clients <- bind_rows(
  prot_scores_all_clients %>%
    pivot_longer(
      cols = all_of(figure1_pathways_all),
      names_to = "Pathway",
      values_to = "Score"
    ) %>%
    transmute(
      sample_id = SampleID,
      Modality = "Protein",
      Stage = recode(
        as.character(EmoryStrictDx.2019),
        Control = "NCI",
        AsymAD = "MCI",
        AD = "AD"
      ),
      Pathway,
      PathwayLabel = fig1_panelA_plain_labels[Pathway],
      PathwayShort = fig1_panelB_labels[Pathway],
      PathwayStrip = fig1_strip_labels[Pathway],
      Score
    ),

  rna_scores_all_clients %>%
    pivot_longer(
      cols = all_of(figure1_pathways_all),
      names_to = "Pathway",
      values_to = "Score"
    ) %>%
    transmute(
      sample_id = sample_id,
      Modality = "RNA",
      Stage = recode(
        as.character(diagnosis_stage),
        Control = "NCI",
        Early_AD = "MCI",
        AD = "AD"
      ),
      Pathway,
      PathwayLabel = fig1_panelA_plain_labels[Pathway],
      PathwayShort = fig1_panelB_labels[Pathway],
      PathwayStrip = fig1_strip_labels[Pathway],
      Score
    )
) %>%
  filter(!is.na(Stage)) %>%
  mutate(
    Stage = factor(Stage, levels = c("NCI", "MCI", "AD")),
    Modality = factor(Modality, levels = c("Protein", "RNA")),
    Pathway = factor(Pathway, levels = figure1_pathways_all),
    PathwayLabel = factor(
      PathwayLabel,
      levels = fig1_panelA_plain_labels[figure1_pathways_all]
    ),
    PathwayShort = factor(
      PathwayShort,
      levels = fig1_panelB_labels[figure1_pathways_all]
    ),
    PathwayStrip = factor(
      PathwayStrip,
      levels = fig1_strip_labels[figure1_pathways_all]
    )
  )

overlay_summary_all_clients <- overlay_long_all_clients %>%
  group_by(Pathway, PathwayLabel, PathwayShort, PathwayStrip, Modality, Stage) %>%
  summarise(
    n = sum(is.finite(Score)),
    mean_score = mean(Score, na.rm = TRUE),
    sem = sd(Score, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

write_fig1_all_table(overlay_long_all_clients, "Fig1A_all_clients_overlay_long")
write_fig1_all_table(overlay_summary_all_clients, "Fig1A_all_clients_overlay_summary")

############################################################
## 7. Panel A significance annotations
## Colored brackets = within-modality stage changes
## Black brackets = Protein vs RNA within same stage
############################################################
make_fig1_all_sig_annotations <- function(overlay_long, overlay_summary) {
  stage_levels <- c("NCI", "MCI", "AD")

  facet_range_tbl <- overlay_summary %>%
    group_by(Pathway, PathwayLabel, PathwayStrip) %>%
    summarise(
      ymax = max(mean_score + sem, na.rm = TRUE),
      ymin = min(mean_score - sem, na.rm = TRUE),
      yrange = ymax - ymin,
      .groups = "drop"
    ) %>%
    mutate(
      yrange = if_else(!is.finite(yrange) | yrange <= 0, 0.1, yrange)
    )

  within_grid <- tidyr::expand_grid(
    pathway_i = unique(overlay_long$Pathway),
    modality_i = c("Protein", "RNA"),
    comparison_i = c("NCI_vs_MCI", "MCI_vs_AD")
  ) %>%
    mutate(
      group1_i = if_else(comparison_i == "NCI_vs_MCI", "NCI", "MCI"),
      group2_i = if_else(comparison_i == "NCI_vs_MCI", "MCI", "AD")
    )

  within_tbl <- purrr::pmap_dfr(
    within_grid,
    function(pathway_i, modality_i, comparison_i, group1_i, group2_i) {
      sub <- overlay_long %>%
        filter(
          as.character(Pathway) == as.character(pathway_i),
          as.character(Modality) == as.character(modality_i)
        )

      x <- sub %>%
        filter(as.character(Stage) == group1_i) %>%
        pull(Score)

      y <- sub %>%
        filter(as.character(Stage) == group2_i) %>%
        pull(Score)

      p <- if (sum(is.finite(x)) >= 2 && sum(is.finite(y)) >= 2) {
        suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
      } else {
        NA_real_
      }

      tibble(
        Pathway = pathway_i,
        Modality = modality_i,
        comparison = comparison_i,
        group1 = group1_i,
        group2 = group2_i,
        n1 = sum(is.finite(x)),
        n2 = sum(is.finite(y)),
        mean1 = mean(x, na.rm = TRUE),
        mean2 = mean(y, na.rm = TRUE),
        p_value = p
      )
    }
  ) %>%
    mutate(
      star = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(star != "") %>%
    checked_left_join(facet_range_tbl, by = "Pathway", label = "Figure 1 within-stage annotations to facet ranges") %>%
    mutate(
      x = as.numeric(factor(group1, levels = stage_levels)),
      xend = as.numeric(factor(group2, levels = stage_levels)),
    base_offset = case_when(
      Modality == "Protein" & comparison == "NCI_vs_MCI" ~ 0.20,
      Modality == "Protein" & comparison == "MCI_vs_AD" ~ 0.28,
      Modality == "RNA" & comparison == "NCI_vs_MCI" ~ 0.36,
      Modality == "RNA" & comparison == "MCI_vs_AD" ~ 0.44,
      TRUE ~ 0.25
    ),
      y = ymax + yrange * base_offset,
      tick_y = y - yrange * 0.03,
      text_y = y + yrange * 0.015,
      annotation_class = as.character(Modality)
    ) %>%
    select(
      Pathway, PathwayLabel, PathwayStrip, annotation_class,
      comparison, group1, group2, n1, n2, mean1, mean2,
      p_value, star, x, xend, y, tick_y, text_y
    )

  between_grid <- tidyr::expand_grid(
    pathway_i = unique(overlay_long$Pathway),
    stage_i = stage_levels
  )

  between_tbl <- purrr::pmap_dfr(
    between_grid,
    function(pathway_i, stage_i) {
      sub <- overlay_long %>%
        filter(
          as.character(Pathway) == as.character(pathway_i),
          as.character(Stage) == stage_i
        )

      x <- sub %>%
        filter(as.character(Modality) == "Protein") %>%
        pull(Score)

      y <- sub %>%
        filter(as.character(Modality) == "RNA") %>%
        pull(Score)

      p <- if (sum(is.finite(x)) >= 2 && sum(is.finite(y)) >= 2) {
        suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
      } else {
        NA_real_
      }

      tibble(
        Pathway = pathway_i,
        Stage = stage_i,
        n_protein = sum(is.finite(x)),
        n_rna = sum(is.finite(y)),
        mean_protein = mean(x, na.rm = TRUE),
        mean_rna = mean(y, na.rm = TRUE),
        p_value = p
      )
    }
  ) %>%
    mutate(
      star = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(star != "") %>%
    checked_left_join(facet_range_tbl, by = "Pathway", label = "Figure 1 between-modality annotations to facet ranges") %>%
    mutate(
      x = as.numeric(factor(Stage, levels = stage_levels)) - 0.18,
      xend = as.numeric(factor(Stage, levels = stage_levels)) + 0.18,
    
      ## Put black Protein-vs-RNA bars just above the higher point at that stage,
      ## instead of floating above all within-stage red/blue brackets.
      y_stage_max = pmax(mean_protein, mean_rna, na.rm = TRUE),
      y = y_stage_max + yrange * 0.23,
      tick_y = y - yrange * 0.025,
      text_y = y + yrange * 0.012,
    
      annotation_class = "Protein vs RNA",
      comparison = paste0("Protein_vs_RNA_at_", Stage)
    ) %>%
    select(
      Pathway, PathwayLabel, PathwayStrip, annotation_class,
      comparison, Stage, n_protein, n_rna, mean_protein, mean_rna,
      p_value, star, x, xend, y, tick_y, text_y
    )

  bind_rows(within_tbl, between_tbl)
}

fig1_sig_ann_all_clients <- make_fig1_all_sig_annotations(
  overlay_long_all_clients,
  overlay_summary_all_clients
)

write_fig1_all_table(
  fig1_sig_ann_all_clients,
  "Fig1A_all_clients_significance_annotations_FIXED"
)

fig1_sig_ann_all_clients %>%
  select(Pathway, annotation_class, comparison, p_value, star) %>%
  arrange(annotation_class, comparison, Pathway) %>%
  print(n = 100)

############################################################
## 8. Panel B effect-difference table
############################################################

effect_diff_all_clients <- overlay_summary_all_clients %>%
  select(Pathway, PathwayLabel, PathwayShort, Modality, Stage, mean_score) %>%
  pivot_wider(
    names_from = c(Modality, Stage),
    values_from = mean_score
  ) %>%
  mutate(
    protein_early_effect = `Protein_MCI` - `Protein_NCI`,
    protein_late_effect = `Protein_AD` - `Protein_MCI`,
    rna_early_effect = `RNA_MCI` - `RNA_NCI`,
    rna_late_effect = `RNA_AD` - `RNA_MCI`,

    early_protein_minus_rna_magnitude = abs(protein_early_effect) - abs(rna_early_effect),
    late_protein_minus_rna_magnitude = abs(protein_late_effect) - abs(rna_late_effect)
  ) %>%
  pivot_longer(
    cols = c(early_protein_minus_rna_magnitude, late_protein_minus_rna_magnitude),
    names_to = "Shift",
    values_to = "protein_minus_rna_magnitude"
  ) %>%
  mutate(
    Shift = recode(
      Shift,
      early_protein_minus_rna_magnitude = "NCI to MCI",
      late_protein_minus_rna_magnitude = "MCI to AD"
    ),
    Shift = factor(
      Shift,
      levels = c("NCI to MCI", "MCI to AD")
    ),
    PathwayShort = factor(
      PathwayShort,
      levels = rev(fig1_panelB_labels[figure1_pathways_all])
    )
  )

write_fig1_all_table(
  effect_diff_all_clients,
  "Fig1B_all_clients_protein_minus_rna_effect_difference"
)

############################################################
## 9. Clean Panel B table
## No bootstrap CIs here. Panel B is an effect-size summary.
## Protein-vs-RNA significance is shown in Panel A.
############################################################

calc_panelB_bootstrap <- function(overlay_long, B = 1000, seed = 1) {
  set.seed(seed)

  stage_pairs <- tribble(
    ~Shift, ~s1, ~s2,
    "NCI to MCI", "NCI", "MCI",
    "MCI to AD", "MCI", "AD"
  )

  boot_grid <- crossing(
    pathway_i = unique(overlay_long$Pathway),
    Shift = stage_pairs$Shift
  ) %>%
    checked_left_join(stage_pairs, by = "Shift", label = "Figure 1 bootstrap grid to stage pairs")

  purrr::pmap_dfr(
    boot_grid,
    function(pathway_i, Shift, s1, s2) {
      sub <- overlay_long %>%
        filter(as.character(Pathway) == as.character(pathway_i))

      p1 <- sub %>% filter(as.character(Modality) == "Protein", as.character(Stage) == s1) %>% pull(Score)
      p2 <- sub %>% filter(as.character(Modality) == "Protein", as.character(Stage) == s2) %>% pull(Score)
      r1 <- sub %>% filter(as.character(Modality) == "RNA", as.character(Stage) == s1) %>% pull(Score)
      r2 <- sub %>% filter(as.character(Modality) == "RNA", as.character(Stage) == s2) %>% pull(Score)

      if (
        sum(is.finite(p1)) < 2 ||
        sum(is.finite(p2)) < 2 ||
        sum(is.finite(r1)) < 2 ||
        sum(is.finite(r2)) < 2
      ) {
        return(tibble(
          Pathway = pathway_i,
          Shift = Shift,
          delta = NA_real_,
          ci_lo = NA_real_,
          ci_hi = NA_real_,
          p_value = NA_real_,
          star = ""
        ))
      }

      p1 <- p1[is.finite(p1)]
      p2 <- p2[is.finite(p2)]
      r1 <- r1[is.finite(r1)]
      r2 <- r2[is.finite(r2)]

      boot_delta <- replicate(B, {
        bp1 <- mean(sample(p1, length(p1), replace = TRUE))
        bp2 <- mean(sample(p2, length(p2), replace = TRUE))
        br1 <- mean(sample(r1, length(r1), replace = TRUE))
        br2 <- mean(sample(r2, length(r2), replace = TRUE))

        abs(bp2 - bp1) - abs(br2 - br1)
      })

      delta_hat <- mean(boot_delta, na.rm = TRUE)
      p_val <- 2 * min(
        mean(boot_delta >= 0, na.rm = TRUE),
        mean(boot_delta <= 0, na.rm = TRUE)
      )

      tibble(
        Pathway = pathway_i,
        Shift = Shift,
        delta = delta_hat,
        ci_lo = as.numeric(quantile(boot_delta, 0.025, na.rm = TRUE)),
        ci_hi = as.numeric(quantile(boot_delta, 0.975, na.rm = TRUE)),
        p_value = p_val,
        star = case_when(
          is.na(p_val) ~ "",
          p_val < 0.001 ~ "***",
          p_val < 0.01 ~ "**",
          p_val < 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    }
  )
}

panelB_sig_tbl <- calc_panelB_bootstrap(
  overlay_long_all_clients,
  B = 1000,
  seed = 1
)

write_fig1_all_table(
  panelB_sig_tbl,
  "Fig1B_all_clients_bootstrap_significance_FIXED"
)

panelB_sig_tbl %>%
  arrange(Shift, Pathway) %>%
  print(n = 100)

############################################################
## 11. Plot Panel B — cleaner effect-size summary
############################################################

############################################################
## 11. Plot Panel B — remodeling magnitude summary
############################################################

pB_tbl <- effect_diff_all_clients %>%
  mutate(
    PathwayShort = factor(
      PathwayShort,
      levels = rev(fig1_panelB_labels[figure1_pathways_all])
    ),
    value_label = sprintf(
      "%.2f",
      if_else(abs(protein_minus_rna_magnitude) < 0.005, 0, protein_minus_rna_magnitude)
    ),
    label_x = case_when(
      protein_minus_rna_magnitude > 0 ~ protein_minus_rna_magnitude + 0.004,
      protein_minus_rna_magnitude < 0 ~ protein_minus_rna_magnitude - 0.004,
      TRUE ~ 0.004
    ),
    label_hjust = case_when(
      protein_minus_rna_magnitude > 0 ~ 0,
      protein_minus_rna_magnitude < 0 ~ 1,
      TRUE ~ 0
    )
  )

x_min_b <- min(pB_tbl$protein_minus_rna_magnitude, na.rm = TRUE)
x_max_b <- max(pB_tbl$protein_minus_rna_magnitude, na.rm = TRUE)

write_fig1_all_table(
  pB_tbl,
  "Fig1B_all_clients_clean_effect_summary"
)

############################################################
## Define paper theme if missing
############################################################

if (!exists("paper_theme_rebuilt")) {
  paper_theme_rebuilt <- function(base_size = 11) {
    theme_classic(base_size = base_size) +
      theme(
        plot.title = element_text(
          face = "bold",
          color = "black",
          size = base_size + 3,
          hjust = 0
        ),
        plot.subtitle = element_text(
          color = "grey35",
          size = base_size - 1,
          hjust = 0
        ),
        axis.title = element_text(
          face = "bold",
          color = "black"
        ),
        axis.text = element_text(
          color = "black"
        ),
        strip.background = element_rect(
          fill = "grey95",
          color = "grey70",
          linewidth = 0.35
        ),
        strip.text = element_text(
          face = "bold",
          color = "black"
        ),
        legend.title = element_text(
          face = "bold",
          color = "black"
        ),
        legend.text = element_text(
          color = "black"
        ),
        panel.grid.major.y = element_line(
          color = "grey90",
          linewidth = 0.25
        ),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
      )
  }
}

############################################################
## 10. Plot Panel A
############################################################

require_columns(
  overlay_summary_all_clients,
  c("Pathway", "PathwayStrip", "Stage", "Modality", "mean_score", "sem"),
  "Figure 1 Panel A overlay summary"
)
require_columns(
  fig1_sig_ann_all_clients,
  c("Pathway", "annotation_class", "x", "xend", "y", "tick_y", "text_y", "star"),
  "Figure 1 Panel A significance annotations"
)
require_values(overlay_summary_all_clients, "Pathway", figure1_pathways_all, "Figure 1 Panel A")

pA <- ggplot(
  overlay_summary_all_clients,
  aes(x = Stage, y = mean_score, color = Modality, group = Modality)
) +
  geom_hline(yintercept = 0, color = "grey82", linewidth = 0.4) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.35) +
  geom_errorbar(
    aes(ymin = mean_score - sem, ymax = mean_score + sem),
    width = 0.06,
    linewidth = 0.42
  ) +
  geom_segment(
    data = fig1_sig_ann_all_clients,
    aes(x = x, xend = xend, y = y, yend = y, color = annotation_class),
    inherit.aes = FALSE,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_segment(
    data = fig1_sig_ann_all_clients,
    aes(x = x, xend = x, y = tick_y, yend = y, color = annotation_class),
    inherit.aes = FALSE,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_segment(
    data = fig1_sig_ann_all_clients,
    aes(x = xend, xend = xend, y = tick_y, yend = y, color = annotation_class),
    inherit.aes = FALSE,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = fig1_sig_ann_all_clients,
    aes(x = (x + xend) / 2, y = text_y, label = star, color = annotation_class),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 2.7,
    vjust = 0,
    show.legend = FALSE
  ) +
  facet_wrap(~ PathwayStrip, nrow = 2, scales = "free_y") +
  scale_color_manual(
    values = modality_colors,
    breaks = c("Protein", "RNA"),
    name = "Modality"
  ) +
  scale_x_discrete(expand = expansion(add = 0.30)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.52))) +
  labs(
    title = "A. All-client pathway remodeling across AD stage",
    subtitle = "Colored brackets show within-modality stage changes; black brackets show Protein vs RNA differences within each stage.",
    x = NULL,
    y = "Mean standardized pathway score"
  ) +
  paper_theme_rebuilt(10.5) +
  theme(
    legend.position = "top",
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(hjust = 0, size = 14),
    plot.subtitle = element_text(hjust = 0, size = 10),
    strip.text = ggtext::element_markdown(size = 8.7, lineheight = 0.90),
    axis.text.x = element_text(size = 7.0, angle = 35, hjust = 1, vjust = 1, margin = margin(t = 1)),
    axis.text.y = element_text(size = 8.6),
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(1.0, "lines"),
    plot.margin = margin(10, 14, 18, 10)
  ) +
  coord_cartesian(clip = "off")

############################################################
## 11. Plot Panel B — clean effect-size summary
############################################################

require_columns(
  pB_tbl,
  c("Pathway", "PathwayShort", "Shift", "protein_minus_rna_magnitude", "value_label", "label_x", "label_hjust"),
  "Figure 1 Panel B effect table"
)
require_values(pB_tbl, "Pathway", figure1_pathways_all, "Figure 1 Panel B")

pB <- ggplot(
  pB_tbl,
  aes(x = protein_minus_rna_magnitude, y = PathwayShort, fill = Pathway)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.45
  ) +
  geom_col(
    width = 0.62,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      x = label_x,
      label = value_label,
      hjust = label_hjust
    ),
    size = 2.55,
    color = "grey20",
    na.rm = TRUE
  ) +
  facet_wrap(~ Shift, ncol = 1) +
  scale_fill_manual(
    values = fig1_pathway_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(x_min_b - 0.025, x_max_b + 0.03),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B. Protein-biased remodeling magnitude",
    subtitle = "Positive values indicate a larger protein-level stage change than RNA-level stage change.",
    x = "Protein remodeling magnitude − RNA remodeling magnitude",
    y = NULL
  ) +
  paper_theme_rebuilt(10.3) +
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0, size = 10),
    strip.text = element_text(face = "bold", size = 9.2),
    axis.text.y = element_text(size = 7.9, lineheight = 0.88),
    axis.text.x = element_text(size = 8.0),
    axis.title.x = element_text(size = 9.4, face = "bold"),
    panel.spacing.y = unit(0.90, "lines"),
    plot.margin = margin(8, 18, 10, 10)
  ) +
  coord_cartesian(clip = "off")

############################################################
## 12. Assemble and save Figure 1
############################################################



fig1_all_clients <- pA / pB +
  plot_layout(heights = c(1.36, 0.92)) +
  plot_annotation(
    title = "All-client Hsp60/10 analysis resolves protein-biased mitochondrial remodeling in AD",
    subtitle = "Hsp60/10 clients are evaluated against broad MitoCarta and non-overlapping mitochondrial/comparator gene sets.",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey35")
    )
  )

save_fig1_all_plot(
  fig1_all_clients,
  "Main_Fig1_all_clients_pathway_remodeling",
  width = 18,
  height = 12.6
)

fig1_all_clients

############################################################
