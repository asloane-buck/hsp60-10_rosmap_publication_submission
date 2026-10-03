############################################################
## 15_make_main_figure_6_candidate_classification.R
##
## Figure 6 candidate-classification synthesis using MitoCarta.
##
## Publication version notes:
## - Panel B simplified:
##     * removes centrality size mapping to reduce clutter
##     * uses a single readable late-decline–Braak landscape
##     * fewer labels
## - Panel C:
##     * uses a muted 50–100% fill range (with squish)
##     * keeps transparent module ranking
## - Panel D:
##     * shows UP TO 6 genes per candidate layer
##     * no misleading claim that each layer has exactly the same count
##     * uses muted 50–100% fill range (with squish)
## - Saves an audit table of how many genes were available/shown per
##   candidate layer.
##
## Recommended run:
##   source("00_config.R")
##   source("01_utils.R")
##   source("15_make_main_figure_6_candidate_classification.R")
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(readxl)
})

############################################################
## 0. Input/output setup
############################################################

if (!exists("cfg")) {
  stop("cfg does not exist. Run source('00_config.R') first.", call. = FALSE)
}

if (!exists("priority_tbl")) {
  priority_file <- file.path(cfg$table_dir, "priority_tbl_COVARIATE_ADJUSTED.csv")
  if (!file.exists(priority_file)) {
    stop(
      "priority_tbl is not in memory and could not find:\n",
      priority_file,
      "\nRun 90_run_main_pipeline.R first, or provide priority_tbl.",
      call. = FALSE
    )
  }
  priority_tbl <- readr::read_csv(priority_file, show_col_types = FALSE)
}

mitocarta_path <- cfg$mitocarta_xls

if (!file.exists(mitocarta_path)) {
  stop(
    "MitoCarta file not found:\n",
    mitocarta_path,
    "\nCheck the path and filename.",
    call. = FALSE
  )
}

out_dir <- file.path(cfg$plot_dir, "main_fig6_candidate_classification")
tab_dir <- file.path(cfg$table_dir, "main_fig6_candidate_classification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

write_fig6_table <- function(x, name) {
  out <- file.path(tab_dir, paste0(name, ".csv"))
  readr::write_csv(x, out)
  message("Wrote table: ", out)
  invisible(x)
}

save_fig6_plot <- function(plot, name, width = 17.8, height = 11.2, dpi = 600) {
  pdf_out <- file.path(out_dir, paste0(name, ".pdf"))
  png_out <- file.path(out_dir, paste0(name, ".png"))
  svg_out <- file.path(out_dir, paste0(name, ".svg"))

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
  toupper(trimws(as.character(x)))
}

clean_col <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
}

as_fraction <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite <- is.finite(x)
  if (any(abs(x[finite]) > 1, na.rm = TRUE)) x / 100 else x
}

rank_to_percentile <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  if (sum(finite) == 0) return(out)
  if (dplyr::n_distinct(x[finite]) <= 1) {
    out[finite] <- 0.5
  } else {
    out[finite] <- dplyr::percent_rank(x[finite])
  }
  out
}

require_any_col <- function(tbl, candidates, label) {
  hit <- intersect(candidates, colnames(tbl))[1]
  if (is.na(hit)) {
    stop(
      "Could not find ", label, " column. Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(tbl), collapse = ", "),
      call. = FALSE
    )
  }
  hit
}

pick_optional_col <- function(tbl, candidates) {
  hit <- intersect(candidates, colnames(tbl))[1]
  if (is.na(hit)) return(NA_character_)
  hit
}

pareto_frontier_flag <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  complete <- stats::complete.cases(mat)
  out <- rep(FALSE, nrow(mat))
  idx <- which(complete)
  if (length(idx) == 0) return(out)

  mm <- mat[idx, , drop = FALSE]

  for (ii in seq_along(idx)) {
    current <- mm[ii, ]
    dominated <- apply(mm, 1, function(candidate) {
      all(candidate >= current, na.rm = TRUE) &&
        any(candidate > current, na.rm = TRUE)
    })
    dominated[ii] <- FALSE
    out[idx[ii]] <- !any(dominated)
  }

  out
}

muted_full_fill_scale <- function(name = "Percentile /\nfraction") {
  scale_fill_gradientn(
    colors = c("grey96", "#EADDD5", "#C28F7B", "#7A4536"),
    limits = c(0.00, 1.00),
    oob = scales::squish,
    labels = percent_format(accuracy = 1),
    name = name
  )
}

muted_high_fill_scale <- function(name = "Percentile") {
  scale_fill_gradientn(
    colors = c("grey96", "#E8DED6", "#C9A28E", "#8C5A44"),
    limits = c(0.50, 1.00),
    oob = scales::squish,
    labels = percent_format(accuracy = 1),
    name = name
  )
}

############################################################
## 2. Read and audit MitoCarta file
############################################################

sheets <- readxl::excel_sheets(mitocarta_path)

sheet_audit <- purrr::map_dfr(sheets, function(sh) {
  tmp <- tryCatch(
    readxl::read_excel(
      mitocarta_path,
      sheet = sh,
      n_max = 5,
      col_types = "text",
      na = c("", "NA")
    ),
    error = function(e) NULL
  )

  if (is.null(tmp)) {
    return(tibble(sheet = sh, n_preview_rows = NA_integer_, n_cols = NA_integer_, columns = NA_character_))
  }

  tibble(
    sheet = sh,
    n_preview_rows = nrow(tmp),
    n_cols = ncol(tmp),
    columns = paste(colnames(tmp), collapse = "; ")
  )
})

write_fig6_table(sheet_audit, "Fig6_MitoCarta_sheet_audit")

read_mitocarta_sheet <- function(sh) {
  raw <- readxl::read_excel(
    mitocarta_path,
    sheet = sh,
    col_types = "text",
    na = c("", "NA")
  )
  names(raw) <- clean_col(names(raw))
  raw
}

candidate_sheets <- purrr::map(sheets, function(sh) {
  tmp <- tryCatch(read_mitocarta_sheet(sh), error = function(e) NULL)
  if (is.null(tmp)) return(NULL)

  cols <- colnames(tmp)

  gene_col <- intersect(
    c("symbol", "gene_symbol", "genesymbol", "hgnc_symbol", "gene", "gene_name"),
    cols
  )[1]

  pathway_col <- intersect(
    c(
      "mitopathways",
      "mito_pathways",
      "mitocarta3_0_mitopathways",
      "mitocarta3_mitopathways",
      "mitocarta_pathways",
      "mitocarta_pathway",
      "pathways",
      "pathway",
      "category",
      "categories"
    ),
    cols
  )[1]

  list(
    sheet = sh,
    table = tmp,
    gene_col = gene_col,
    pathway_col = pathway_col,
    n_rows = nrow(tmp),
    n_cols = ncol(tmp)
  )
})

candidate_sheets <- candidate_sheets[!purrr::map_lgl(candidate_sheets, is.null)]

sheet_choice_tbl <- purrr::map_dfr(candidate_sheets, function(x) {
  tibble(
    sheet = x$sheet,
    n_rows = x$n_rows,
    n_cols = x$n_cols,
    gene_col = ifelse(is.na(x$gene_col), NA_character_, x$gene_col),
    pathway_col = ifelse(is.na(x$pathway_col), NA_character_, x$pathway_col),
    score = as.integer(!is.na(x$gene_col)) * 10 +
      as.integer(!is.na(x$pathway_col)) * 5 +
      log10(max(x$n_rows, 1))
  )
}) %>%
  arrange(desc(score), desc(n_rows))

write_fig6_table(sheet_choice_tbl, "Fig6_MitoCarta_sheet_choice_audit")

chosen_sheet <- sheet_choice_tbl$sheet[1]
chosen <- candidate_sheets[[which(purrr::map_chr(candidate_sheets, "sheet") == chosen_sheet)[1]]]

if (is.na(chosen$gene_col)) {
  stop(
    "Could not identify a gene symbol column in MitoCarta file.\n",
    "See Fig6_MitoCarta_sheet_choice_audit.csv for sheet/column audit.",
    call. = FALSE
  )
}

mito_raw <- chosen$table
gene_col <- chosen$gene_col
pathway_col <- chosen$pathway_col
desc_col <- intersect(c("description", "protein_description", "name", "gene_description"), colnames(mito_raw))[1]
sub_loc_col <- intersect(
  c(
    "sub_mito_localization",
    "submito_localization",
    "mitocarta3_0_submito_localization",
    "submitochondrial_location",
    "localization"
  ),
  colnames(mito_raw)
)[1]

mitocarta_tbl <- mito_raw %>%
  transmute(
    gene = clean_gene(.data[[gene_col]]),
    mitocarta_sheet = chosen_sheet,
    mitocarta_pathways_raw = if (!is.na(pathway_col)) as.character(.data[[pathway_col]]) else NA_character_,
    mitocarta_description = if (!is.na(desc_col)) as.character(.data[[desc_col]]) else NA_character_,
    mitocarta_submito = if (!is.na(sub_loc_col)) as.character(.data[[sub_loc_col]]) else NA_character_
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)

write_fig6_table(mitocarta_tbl, "Fig6_MitoCarta_join_table_cleaned")

############################################################
## 3. Functional class from MitoCarta pathways
############################################################

classify_mitocarta <- function(gene, pathways, description = NA_character_) {
  g <- clean_gene(gene)
  x <- paste(pathways, description, sep = " | ")
  x <- stringr::str_to_lower(x)
  x[is.na(x)] <- ""

  dplyr::case_when(
    stringr::str_detect(x, "ribosom|translation|mitoribosom|mitochondrial central dogma|mrna|trna|rrna|rna metabolism|transcription") |
      stringr::str_detect(g, "^MRPL|^MRPS") |
      g %in% c("TUFM", "TSFM", "GFM1", "GFM2", "DAP3", "LRPPRC", "SLIRP", "PTCD3") ~
      "Mito translation / ribosome",

    stringr::str_detect(x, "tca|tricarboxylic|pyruvate|carbohydrate|ketoglutarate|dehydrogenase|nad|redox|fe-s|heme|iron-sulfur|organic acid") |
      g %in% c(
        "PDHA1", "PDHB", "PDHX", "DLAT", "DLD", "DLST", "OGDH",
        "IDH3A", "IDH3B", "IDH3G", "IDH2", "MDH2", "CS", "FH", "ACO2",
        "SUCLA2", "SUCLG1", "SUCLG2", "GOT2", "SHMT2", "GLUD1", "GLUD2",
        "OAT", "ALDH2", "ME2"
      ) ~
      "TCA / pyruvate / redox metabolism",

    stringr::str_detect(x, "oxidative phosphorylation|oxphos|respiratory chain|electron transport|complex i|complex ii|complex iii|complex iv|complex v|atp synthase") |
      stringr::str_detect(g, "^NDUF|^UQCR|^COX|^SDH|^ATP5|^MT-") ~
      "OXPHOS / ETC",

    stringr::str_detect(x, "fatty acid|lipid|beta.oxidation|β.oxidation|acyl|carnitine") |
      g %in% c("ECHS1", "ECH1", "HADHA", "HADHB", "ACAA2", "ACADVL", "ETFA", "ETFB", "ETFDH") ~
      "FAO / lipid metabolism",

    stringr::str_detect(x, "protein import|protein sorting|protein homeostasis|proteostasis|protein folding|chaperone|protease|quality control|stress|detoxification|antioxidant|ros") |
      g %in% c("HSPA9", "TRAP1", "CLPX", "CLPP", "LONP1", "AFG3L2", "SPG7", "PHB", "PHB2", "PRDX3", "SOD2", "TXN2") ~
      "Proteostasis / stress",

    stringr::str_detect(x, "transport|carrier|membrane|import|export|transloc|dynamics|fission|fusion|signaling|apoptosis|calcium") ~
      "Transport / membrane / signaling",

    stringr::str_detect(x, "metabolism|amino acid|nucleotide|one-carbon|folate|cofactor|vitamin|porphyrin|ubiquinone|coq") ~
      "Other mitochondrial metabolism",

    x != "" ~ "Other MitoCarta-annotated clients",
    TRUE ~ "No MitoCarta annotation"
  )
}

functional_levels <- c(
  "Mito translation / ribosome",
  "TCA / pyruvate / redox metabolism",
  "OXPHOS / ETC",
  "FAO / lipid metabolism",
  "Proteostasis / stress",
  "Transport / membrane / signaling",
  "Other mitochondrial metabolism",
  "Other MitoCarta-annotated clients",
  "No MitoCarta annotation"
)

functional_pal <- c(
  "Mito translation / ribosome" = "#7A2E55",          # muted plum
  "TCA / pyruvate / redox metabolism" = "#C8962E",   # muted gold; distinct from OXPHOS
  "OXPHOS / ETC" = "#3B6EA8",                        # clear muted blue
  "FAO / lipid metabolism" = "#7A8F2A",              # olive
  "Proteostasis / stress" = "#4F8A5B",               # green
  "Transport / membrane / signaling" = "#7B6EA8",    # lavender
  "Other mitochondrial metabolism" = "#8E6C8A",      # mauve
  "Other MitoCarta-annotated clients" = "grey58",
  "No MitoCarta annotation" = "grey82"
)

functional_legend_labels <- c(
  "Mito translation / ribosome" = "Translation / ribosome",
  "TCA / pyruvate / redox metabolism" = "TCA / pyruvate / redox",
  "OXPHOS / ETC" = "OXPHOS / ETC",
  "FAO / lipid metabolism" = "FAO / lipid",
  "Proteostasis / stress" = "Proteostasis / stress",
  "Transport / membrane / signaling" = "Transport / membrane / signaling",
  "Other mitochondrial metabolism" = "Other mito metabolism",
  "Other MitoCarta-annotated clients" = "Other MitoCarta",
  "No MitoCarta annotation" = "No annotation"
)

############################################################
## 4. Join priority table to MitoCarta
############################################################

late_decline_col <- require_any_col(
  priority_tbl,
  c("late_decline_percentile", "late_decline_pct"),
  "late-stage decline percentile"
)
joint_pathology_col <- require_any_col(
  priority_tbl,
  c(
    "joint_pathology_percentile",
    "joint_pathology_pct"
  ),
  "joint AD-pathology percentile"
)

strict_joint_pathology_col <- require_any_col(
  priority_tbl,
  c(
    "strict_joint_pathology_percentile",
    "strict_joint_pathology_pct"
  ),
  "strict joint AD-pathology percentile"
)

## Legacy Braak-only axis retained solely for migration
## validation against the frozen Figure 6 frontier.
legacy_braak_col <- require_any_col(
  priority_tbl,
  c(
    "inverse_braak_percentile",
    "braak_pct"
  ),
  "legacy inverse Braak percentile"
)
centrality_col <- require_any_col(
  priority_tbl,
  c("centrality_percentile", "centrality_pct"),
  "centrality percentile"
)

agora_col <- pick_optional_col(
  priority_tbl,
  c("agora_target", "agora_flag", "agora_nominated_target", "AGORA_target", "is_agora_target")
)

if (is.na(agora_col)) {
  priority_tbl$agora_target_tmp <- 0
  agora_col <- "agora_target_tmp"
}

fig6_tbl <- priority_tbl %>%
  mutate(gene = clean_gene(gene)) %>%
  left_join(mitocarta_tbl, by = "gene") %>%
  mutate(
    late_decline_pct =
      as_fraction(
        .data[[late_decline_col]]
      ),

    ## Primary pathology dimension:
    ## one equal-weight joint Braak/CERAD percentile axis.
    joint_pathology_pct =
      as_fraction(
        .data[[joint_pathology_col]]
      ),

    ## Strict convergence sensitivity:
    ## percentile of the weaker Braak/CERAD pathology score.
    strict_joint_pathology_pct =
      as_fraction(
        .data[[strict_joint_pathology_col]]
      ),

    ## Legacy Braak-only dimension retained for exact migration
    ## comparison; it is not used in the revised primary figure.
    legacy_braak_pct =
      as_fraction(
        .data[[legacy_braak_col]]
      ),

    centrality_pct =
      as_fraction(
        .data[[centrality_col]]
      ),
    agora_binary = as.integer(as.numeric(.data[[agora_col]]) > 0),
    functional_class = classify_mitocarta(gene, mitocarta_pathways_raw, mitocarta_description),
    functional_class = factor(functional_class, levels = functional_levels)
  )

############################################################
## 4b. Add client-level cognition support
##
## Preferred source is the integrated Figure 4 output. If priority_tbl
## already contains a cognition percentile/score, that is used instead.
############################################################

cognition_pct_col <- pick_optional_col(
  fig6_tbl,
  c("cognition_percentile", "cognition_pct", "cognitive_percentile", "cognitive_pct")
)

cognition_score_col <- pick_optional_col(
  fig6_tbl,
  c("cognition_priority_score", "cognition_support_score", "cognition_preservation_score")
)

if (!is.na(cognition_pct_col)) {
  fig6_tbl <- fig6_tbl %>%
    mutate(
      cognition_pct = as_fraction(.data[[cognition_pct_col]]),
      cognition_source = cognition_pct_col
    )
} else if (!is.na(cognition_score_col)) {
  fig6_tbl <- fig6_tbl %>%
    mutate(
      cognition_raw_score = suppressWarnings(as.numeric(.data[[cognition_score_col]])),
      cognition_pct = rank_to_percentile(cognition_raw_score),
      cognition_source = cognition_score_col
    )
} else {
  fig4_cognition_candidates <- c(
    file.path(cfg$table_dir, "main_fig4_cognition", "Fig4_client_level_cognition_summary_by_gene.csv"),
    file.path(cfg$plot_dir, "client_level_cognition_prioritization", "client_level_cognition_summary_by_gene.csv")
  )
  fig4_cognition_file <- fig4_cognition_candidates[file.exists(fig4_cognition_candidates)][1]

  if (is.na(fig4_cognition_file)) {
    stop(
      "Figure 6 needs client-level cognition output, but no cognition columns were found in priority_tbl ",
      "and no Figure 4 cognition summary file was found.
Tried:
",
      paste(fig4_cognition_candidates, collapse = "
"),
      call. = FALSE
    )
  }

  cognition_tbl <- readr::read_csv(fig4_cognition_file, show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    mutate(gene = clean_gene(gene))

  cognition_score_col_file <- pick_optional_col(
    cognition_tbl,
    c("cognition_priority_score", "cognition_support_score", "cognition_preservation_score")
  )

  if (is.na(cognition_score_col_file)) {
    stop(
      "Found cognition summary file but could not identify a cognition score column: ",
      fig4_cognition_file,
      "
Available columns: ", paste(colnames(cognition_tbl), collapse = ", "),
      call. = FALSE
    )
  }

  n_fdr_col_file <- pick_optional_col(cognition_tbl, c("n_fdr_better_cognition"))

  cognition_tbl <- cognition_tbl %>%
    transmute(
      gene,
      cognition_raw_score = suppressWarnings(as.numeric(.data[[cognition_score_col_file]])),
      n_fdr_better_cognition = if (!is.na(n_fdr_col_file)) {
        suppressWarnings(as.numeric(.data[[n_fdr_col_file]]))
      } else {
        NA_real_
      }
    ) %>%
    mutate(cognition_pct = rank_to_percentile(cognition_raw_score))

  fig6_tbl <- fig6_tbl %>%
    left_join(cognition_tbl, by = "gene") %>%
    mutate(cognition_source = basename(fig4_cognition_file))
}

fig6_tbl <- fig6_tbl %>%
  mutate(
    cognition_pct = as.numeric(cognition_pct),
    cognition_pct = if_else(is.na(cognition_pct), 0, cognition_pct)
  )

## Primary Figure 6 Pareto frontier.
fig6_tbl$pareto_frontier <- pareto_frontier_flag(
  fig6_tbl %>%
    select(
      late_decline_pct,
      joint_pathology_pct,
      cognition_pct,
      centrality_pct
    )
)

## Strict joint-pathology sensitivity frontier.
fig6_tbl$strict_joint_pareto_frontier <- pareto_frontier_flag(
  fig6_tbl %>%
    select(
      late_decline_pct,
      strict_joint_pathology_pct,
      cognition_pct,
      centrality_pct
    )
)

## Frozen Braak-only frontier retained for migration validation.
fig6_tbl$legacy_braak_pareto_frontier <- pareto_frontier_flag(
  fig6_tbl %>%
    select(
      late_decline_pct,
      legacy_braak_pct,
      cognition_pct,
      centrality_pct
    )
)

fig6_tbl <- fig6_tbl %>%
  mutate(
    mean_axis_percentile_for_display_only =
      rowMeans(
        cbind(
          late_decline_pct,
          joint_pathology_pct,
          cognition_pct,
          centrality_pct
        ),
        na.rm = TRUE
      ),
    candidate_layer = case_when(
      pareto_frontier & agora_binary == 1 ~ "Pareto frontier + AGORA",
      pareto_frontier ~ "Pareto frontier",
      agora_binary == 1 ~ "AGORA-supported",
      TRUE ~ "Other classified client"
    ),
    candidate_layer = factor(
      candidate_layer,
      levels = c("Pareto frontier + AGORA", "Pareto frontier", "AGORA-supported", "Other classified client")
    )
  )

write_fig6_table(fig6_tbl, "Fig6_full_candidate_classification")

join_audit <- fig6_tbl %>%
  summarise(
    n_priority_genes = n(),
    n_joined_to_mitocarta = sum(!is.na(mitocarta_pathways_raw) | !is.na(mitocarta_description)),
    fraction_joined = n_joined_to_mitocarta / n_priority_genes,
    n_no_mitocarta = sum(functional_class == "No MitoCarta annotation", na.rm = TRUE),
    mitocarta_sheet_used = chosen_sheet,
    gene_column_used = gene_col,
    pathway_column_used = ifelse(is.na(pathway_col), "NONE", pathway_col)
  )

write_fig6_table(join_audit, "Fig6_MitoCarta_join_audit")

############################################################
## 5. Summaries and display tables
############################################################

category_summary <- fig6_tbl %>%
  group_by(functional_class) %>%
  summarise(
    n_genes = n(),
    median_late_decline =
      median(
        late_decline_pct,
        na.rm = TRUE
      ),

    median_joint_pathology =
      median(
        joint_pathology_pct,
        na.rm = TRUE
      ),

    median_cognition =
      median(
        cognition_pct,
        na.rm = TRUE
      ),
    median_centrality = median(centrality_pct, na.rm = TRUE),
    pareto_fraction = mean(pareto_frontier, na.rm = TRUE),
    agora_fraction = mean(agora_binary == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_genes > 0) %>%
  arrange(
    desc(pareto_fraction),
    desc(median_cognition),
    desc(median_joint_pathology),
    desc(median_late_decline),
    desc(median_centrality)
  ) %>%
  mutate(module_rank = row_number())

write_fig6_table(category_summary, "Fig6_MitoCarta_functional_class_summary_ranked")

heatmap_tbl <- category_summary %>%
  select(
    module_rank, functional_class, n_genes,
    median_late_decline, median_joint_pathology, median_cognition, median_centrality,
    pareto_fraction, agora_fraction
  ) %>%
  pivot_longer(
    cols = c(
      median_late_decline,
      median_joint_pathology,
      median_cognition,
      median_centrality,
      pareto_fraction,
      agora_fraction
    ),
    names_to = "evidence_layer_raw",
    values_to = "value"
  ) %>%
  mutate(
    evidence_layer = dplyr::recode(
      evidence_layer_raw,
      median_late_decline = "Late decline",
      median_joint_pathology = "Joint pathology",
      median_cognition = "Cognition",
      median_centrality = "Centrality",
      pareto_fraction = "Frontier\nfraction",
      agora_fraction = "AGORA\nfraction"
    ),
    evidence_layer = factor(
      evidence_layer,
      levels = c("Late decline", "Joint pathology", "Cognition", "Centrality", "Frontier\nfraction", "AGORA\nfraction")
    ),
    row_label = paste0(module_rank, ". ", functional_class, " (n=", n_genes, ")"),
    row_label = factor(row_label, levels = rev(unique(row_label)))
  )

candidate_rank_table <- fig6_tbl %>%
  arrange(candidate_layer, desc(mean_axis_percentile_for_display_only), desc(agora_binary), desc(cognition_pct), desc(centrality_pct)) %>%
  group_by(candidate_layer) %>%
  mutate(display_rank_within_layer = row_number()) %>%
  ungroup() %>%
  select(
    display_rank_within_layer, gene, candidate_layer, functional_class,
    late_decline_pct,
    joint_pathology_pct,
    strict_joint_pathology_pct,
    legacy_braak_pct,
    cognition_pct,
    centrality_pct,
    pareto_frontier,
    strict_joint_pareto_frontier,
    legacy_braak_pareto_frontier,
    agora_binary,
    mean_axis_percentile_for_display_only,
    mitocarta_pathways_raw, mitocarta_description, mitocarta_submito, everything()
  )

write_fig6_table(candidate_rank_table, "Fig6_full_ranked_candidate_table_for_supplement")

display_n_per_layer <- 6

candidate_layer_counts <- candidate_rank_table %>%
  filter(candidate_layer != "Other classified client") %>%
  count(candidate_layer, name = "available_n") %>%
  mutate(shown_n = pmin(available_n, display_n_per_layer))

write_fig6_table(candidate_layer_counts, "Fig6_candidate_layer_display_counts")


candidate_layer_class_counts <- candidate_rank_table %>%
  filter(candidate_layer != "Other classified client") %>%
  count(candidate_layer, functional_class, name = "n") %>%
  group_by(candidate_layer) %>%
  mutate(
    layer_total = sum(n),
    fraction_within_layer = n / layer_total
  ) %>%
  ungroup() %>%
  mutate(
    candidate_layer = factor(
      candidate_layer,
      levels = c("Pareto frontier + AGORA", "Pareto frontier", "AGORA-supported")
    ),
    functional_class = factor(functional_class, levels = functional_levels)
  )

write_fig6_table(candidate_layer_class_counts, "Fig6_candidate_layer_by_MitoCarta_class_counts")


candidate_display_tbl <- candidate_rank_table %>%
  filter(candidate_layer != "Other classified client") %>%
  group_by(candidate_layer) %>%
  slice_head(n = display_n_per_layer) %>%
  ungroup() %>%
  arrange(candidate_layer, display_rank_within_layer) %>%
  mutate(
    row_label = gene,
    row_label = factor(row_label, levels = rev(unique(row_label)))
  )

write_fig6_table(candidate_display_tbl, "Fig6_panelD_display_candidates_rule_based")

candidate_heatmap_tbl <- candidate_display_tbl %>%
  select(
    row_label,
    candidate_layer,
    functional_class,
    late_decline_pct,
    joint_pathology_pct,
    cognition_pct,
    centrality_pct
  ) %>%
  pivot_longer(
    cols = c(
      late_decline_pct,
      joint_pathology_pct,
      cognition_pct,
      centrality_pct
    ),
    names_to = "axis",
    values_to = "percentile"
  ) %>%
  mutate(
    axis = recode(
      axis,
      late_decline_pct = "Late decline",
      joint_pathology_pct = "Joint pathology",
      cognition_pct = "Cognition",
      centrality_pct = "Centrality"
    ),
    axis = factor(axis, levels = c("Late decline", "Joint pathology", "Cognition", "Centrality"))
  )

candidate_class_dot_tbl <- candidate_display_tbl %>%
  distinct(row_label, candidate_layer, functional_class) %>%
  mutate(axis = factor("Class", levels = c("Class", "Late decline", "Joint pathology", "Cognition", "Centrality")))

candidate_layer_class_counts <- fig6_tbl %>%
  filter(candidate_layer != "Other classified client") %>%
  count(functional_class, candidate_layer, name = "n") %>%
  group_by(functional_class) %>%
  mutate(total_candidate_layer_n = sum(n)) %>%
  ungroup() %>%
  mutate(functional_class = factor(functional_class, levels = rev(functional_levels)))

write_fig6_table(candidate_layer_class_counts, "Fig6_candidate_layer_counts_by_MitoCarta_class")

############################################################
## 6. Plot helpers
############################################################

theme_fig6 <- function(base_size = 10.0) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size - 0.7, color = "grey35"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(6, 7, 6, 6)
    )
}

## Fewer labels to calm Panel B.
label_tbl <- bind_rows(
  candidate_rank_table %>%
    filter(candidate_layer == "Pareto frontier + AGORA") %>%
    slice_head(n = 6),
  candidate_rank_table %>%
    filter(candidate_layer == "Pareto frontier") %>%
    slice_head(n = 4)
) %>%
  distinct(gene, .keep_all = TRUE)

############################################################
## 7. Panel A schematic
############################################################

schem_tbl <- tibble(
  xmin = c(0.06, 1.42, 3.98, 5.62),
  xmax = c(1.14, 3.74, 5.38, 6.86),
  ymin = c(0.44, 0.44, 0.44, 0.44),
  ymax = c(1.66, 1.66, 1.66, 1.66),
  header = c("Input", "Four prioritization axes", "Support/context layers", "Output"),
  body = c(
    "Detected Hsp60/10\nclients",
    "Late-stage protein decline\nJoint Braak/CERAD pathology\nCognition preservation\nNetwork centrality",
    "Pareto/frontier status\nAGORA support\nMitoCarta class",
    "Rule-based\ncandidate classes\nfor mechanistic\nfollow-up"
  ),
  fill = c("#F4F4F4", "#EEF1E6", "#F2ECE6", "#EAF1F6"),
  header_fill = c("#E3E3E3", "#DCE2CF", "#E6DCD3", "#DCE8F0"),
  body_size = c(2.25, 2.02, 2.04, 1.95)
)

arrow_tbl <- tibble(
  x = schem_tbl$xmax[-nrow(schem_tbl)] + 0.08,
  xend = schem_tbl$xmin[-1] - 0.08,
  y = 1.05,
  yend = 1.05
)

pA <- ggplot() +
  xlim(-0.02, 6.92) +
  ylim(0.36, 1.74) +
  geom_rect(
    data = schem_tbl,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = "black",
    linewidth = 0.24
  ) +
  geom_rect(
    data = schem_tbl,
    aes(xmin = xmin, xmax = xmax, ymin = ymax - 0.34, ymax = ymax, fill = header_fill),
    color = "black",
    linewidth = 0.24
  ) +
  geom_text(
    data = schem_tbl,
    aes(x = (xmin + xmax) / 2, y = ymax - 0.17, label = header),
    size = 2.12,
    fontface = "bold",
    lineheight = 0.86
  ) +
  geom_text(
    data = schem_tbl,
    aes(x = (xmin + xmax) / 2, y = ymin + 0.42, label = body, size = body_size),
    lineheight = 0.88
  ) +
  scale_fill_identity() +
  scale_size_identity() +
  geom_segment(
    data = arrow_tbl,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.15, "cm")),
    linewidth = 0.42
  ) +
  labs(
    title = "A. Candidate classification framework",
    subtitle = "Rule-based prioritization keeps quantitative axes separate and adds support/context layers for follow-up."
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 9.8) +
  theme(
    plot.title = element_text(face = "bold", size = 12.8),
    plot.subtitle = element_text(size = 8.7, color = "grey35"),
    plot.margin = margin(2, 10, 2, 10)
  )

############################################################
## 8. Panel B simplified collapse-Braak landscape
############################################################

pB <- ggplot(
  fig6_tbl,
  aes(
    x = late_decline_pct,
    y = joint_pathology_pct
  )
) +
  geom_point(
    data = fig6_tbl %>% filter(candidate_layer == "Other classified client"),
    color = "grey83",
    shape = 16,
    size = 1.10,
    alpha = 0.36
  ) +
  geom_point(
    data = fig6_tbl %>% filter(candidate_layer != "Other classified client"),
    aes(color = functional_class),
    shape = 16,
    size = 2.25,
    alpha = 0.82
  ) +
  geom_point(
    data = fig6_tbl %>% filter(agora_binary == 1),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 0.62,
    size = 3.25,
    alpha = 0.95
  ) +
  geom_point(
    data = fig6_tbl %>% filter(pareto_frontier),
    shape = 21,
    fill = NA,
    color = "grey15",
    stroke = 0.90,
    size = 3.85,
    alpha = 0.95
  ) +
  scale_color_manual(values = functional_pal, drop = FALSE, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "B. Late-decline / joint-pathology target landscape",
    subtitle = "Candidate-layer genes are colored by MitoCarta class; unlabeled background clients are grey. \nPareto rings use late-stage decline, joint pathology, cognition, and centrality axes.",
    x = "Late-stage protein decline\npercentile",
    y = "Joint AD pathology\npercentile"
  ) +
  theme_fig6(10.2) +
  theme(
    plot.title = element_text(face = "bold", size = 13.8),
    plot.subtitle = element_text(size = 9.1, color = "grey35"),
    axis.text = element_text(size = 6.7),
    axis.title.x = element_text(face = "bold", size = 8.5, lineheight = 0.88),
    axis.title.y = element_text(face = "bold", size = 7.1, lineheight = 0.78, margin = margin(r = 2)),
    plot.margin = margin(4, 14, 5, 9)
  )

############################################################
## 9. Panel C ranked functional-class summary heatmap
############################################################

pC <- ggplot(heatmap_tbl, aes(x = evidence_layer, y = row_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(
    aes(label = percent(value, accuracy = 1)),
    size = 2.15,
    fontface = "bold"
  ) +
  muted_full_fill_scale(name = "Percentile /\nfraction") +
  labs(
    title = "C. Ranked MitoCarta functional-class evidence summary",
    subtitle = "Rows are ordered by frontier fraction, then cognition, joint pathology, late-stage decline, and centrality.",
    x = NULL, y = NULL
  ) +
  theme_fig6(9.8) +
  guides(
    fill = guide_colorbar(
      barheight = unit(1.0, "cm"),
      barwidth = unit(0.18, "cm"),
      title.position = "top"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1, face = "bold", size = 6.1, lineheight = 0.88),
    axis.text.y = element_text(face = "bold", size = 6.6),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 6.9, lineheight = 0.88),
    legend.text = element_text(size = 6.2),
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 8.9, color = "grey35"),
    plot.margin = margin(4, 8, 5, 8)
  )

############################################################
## 10. Panel D candidate evidence table/heatmap
############################################################

pD <- ggplot() +
  geom_tile(
    data = candidate_class_dot_tbl,
    aes(x = axis, y = row_label),
    fill = "grey96",
    color = "white",
    linewidth = 0.36
  ) +
  geom_point(
    data = candidate_class_dot_tbl,
    aes(x = axis, y = row_label, color = functional_class),
    size = 2.35,
    alpha = 0.95
  ) +
  geom_tile(
    data = candidate_heatmap_tbl,
    aes(x = axis, y = row_label, fill = percentile),
    color = "white",
    linewidth = 0.36
  ) +
  geom_text(
    data = candidate_heatmap_tbl,
    aes(x = axis, y = row_label, label = percent(percentile, accuracy = 1)),
    size = 1.86
  ) +
  facet_grid(
    candidate_layer ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  muted_high_fill_scale(name = "Percentile") +
  scale_color_manual(values = functional_pal, drop = FALSE, guide = "none") +
  labs(
    title = "D. Rule-based candidate-layer evidence table",
    subtitle = paste0(
      "Rows show up to ", display_n_per_layer,
      " genes per candidate layer when available. Class dots use the same MitoCarta colors as Panel C/B. \nOrdering uses the unweighted mean of late-stage decline, joint pathology, cognition, and centrality percentiles for display only; \nfull ranked table is saved."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_fig6(9.8) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey92", color = "grey80"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 6.3, margin = margin(2, 2, 2, 2)),
    axis.text.x = element_text(angle = 34, hjust = 1, face = "bold", size = 6.2),
    axis.text.y = element_text(face = "bold", size = 5.65),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 8.7, color = "grey35"),
    panel.spacing.y = unit(0.62, "lines"),
    plot.margin = margin(4, 10, 8, 8)
  )

############################################################
## 10b. Panel E candidate-layer composition
############################################################

pE <- ggplot(
  candidate_layer_class_counts,
  aes(x = candidate_layer, y = functional_class)
) +
  geom_point(
    aes(size = n, color = functional_class),
    alpha = 0.88
  ) +
  geom_text(
    aes(label = n),
    size = 2.8,
    fontface = "bold",
    color = "white"
  ) +
  scale_color_manual(values = functional_pal, drop = FALSE, guide = "none") +
  scale_size_continuous(range = c(3.0, 8.5), guide = "none") +
  labs(
    title = "E. Candidate-layer composition by MitoCarta class",
    subtitle = "Counts show how many genes in each functional class are Pareto+AGORA, Pareto-only, or AGORA-supported.",
    x = NULL,
    y = NULL
  ) +
  theme_fig6(9.5) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", size = 7.8),
    axis.text.y = element_text(face = "bold", size = 7.0),
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 8.2, color = "grey35"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25)
  )


############################################################
## 10b. Panel E candidate-layer composition by MitoCarta class
############################################################

candidate_layer_class_counts <- candidate_layer_class_counts %>%
  mutate(
    candidate_layer = factor(
      candidate_layer,
      levels = c("Pareto frontier + AGORA", "Pareto frontier", "AGORA-supported")
    )
  )

pE <- ggplot(
  candidate_layer_class_counts,
  aes(y = candidate_layer, x = n, fill = functional_class)
) +
  geom_col(width = 0.64, color = "white", linewidth = 0.25) +
  geom_text(
    data = candidate_layer_class_counts %>% filter(n >= 3),
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 2.7,
    fontface = "bold",
    color = "white"
  ) +
  scale_fill_manual(
    values = functional_pal,
    drop = FALSE,
    name = "MitoCarta class",
    labels = function(x) stringr::str_wrap(unname(functional_legend_labels[x]), width = 22)
  ) +
  scale_y_discrete(
    labels = c(
      "Pareto frontier + AGORA" = "Pareto + AGORA",
      "Pareto frontier" = "Pareto frontier",
      "AGORA-supported" = "AGORA-supported"
    )
  ) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(linewidth = 0.2))) +
  labs(
    title = "E. Candidate-layer composition",
    subtitle = "Counts show which MitoCarta classes contribute to each candidate layer.",
    x = "Number of candidates",
    y = NULL
  ) +
  theme_fig6(9.8) +
  theme(
    axis.text.x = element_text(size = 7.9),
    axis.text.y = element_text(face = "bold", size = 7.1, margin = margin(r = 2)),
    axis.title.x = element_text(face = "bold", size = 8.6),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 6.4),
    legend.text = element_text(size = 5.3, lineheight = 0.82),
    legend.key.size = unit(0.17, "cm"),
    legend.spacing.x = unit(3, "pt"),
    legend.spacing.y = unit(0, "pt"),
    plot.title = element_text(face = "bold", size = 13.0),
    plot.subtitle = element_text(size = 8.5, color = "grey35"),
    panel.grid = element_blank(),
    plot.margin = margin(4, 6, 3, 24)
  )

############################################################
## 11. Assemble
############################################################

fig6_candidate_classification <- patchwork::free(pA) / patchwork::free(pB) / pC / pD / patchwork::free(pE) +
  plot_layout(
    heights = c(0.72, 1.14, 0.80, 1.52, 1.25),
    guides = "keep"
  ) +
  plot_annotation(
    title = "Integrated Hsp60/10 client evidence defines candidate classes for mechanistic follow-up",
    subtitle = "Late-stage decline, joint Braak/CERAD pathology, cognition support, and centrality are continuous axes; AGORA support and MitoCarta-derived functional classes annotate rule-based candidate layers.",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17.2),
      plot.subtitle = element_text(hjust = 0.5, size = 11.0, color = "grey35"),
      plot.margin = margin(8, 12, 8, 12)
    )
  )

save_fig6_plot(
  fig6_candidate_classification,
  "Figure6_candidate_classification",
  width = 17.8,
  height = 11.2
)

fig6_candidate_classification

message("\n============================================================")
message("Figure 6 candidate-classification with cognition and MitoCarta annotations complete.")
message("MitoCarta sheet used: ", chosen_sheet)
message("Priority genes joined to MitoCarta: ", join_audit$n_joined_to_mitocarta, " / ", join_audit$n_priority_genes)
message("Displayed candidate counts by layer:")
print(candidate_layer_counts)
message("Outputs written to:")
message(out_dir)
message(tab_dir)
message("============================================================")
