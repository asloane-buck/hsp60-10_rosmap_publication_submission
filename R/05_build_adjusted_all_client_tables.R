############################################################
## 05_build_adjusted_all_client_tables.R
## Build source-of-truth all_hsp60_10_client_tbl using
## covariate-adjusted RNA/protein inputs for stage/late-decline/centrality
## and processed non-residualized protein input for pathology beta models.
##
## Inputs expected from scripts 00-04:
##   - cfg
##   - prot_mat, rna_mat                    # adjusted matrices from 03
##   - prot_meta_adj, rna_meta_adj
##   - protein_covars
##   - all_hsp60_10_clients
##   - helper functions from 01_utils.R
##
## Outputs:
##   - all_hsp60_10_client_tbl
##   - all_client_stage_long
##   - priority_tbl
##   - saved CSV/RDS outputs
############################################################

required_objects <- c(
  "cfg", "prot_mat", "rna_mat", "prot_mat_raw", "prot_meta_adj", "rna_meta_adj",
  "protein_covars", "all_hsp60_10_clients"
)
require_objects(required_objects, context = "05_build_adjusted_all_client_tables.R")

required_functions <- c(
  "clean_gene_symbols", "canonical_gene_symbol", "display_gene_symbol",
  "fit_adjusted_braak_beta", "fit_adjusted_cerad_beta", "percentile01",
  "write_tbl", "save_obj"
)
missing_functions <- required_functions[!vapply(required_functions, exists, logical(1))]
if (length(missing_functions) > 0) {
  stop(
    "05_build_adjusted_all_client_tables.R is missing required helper functions: ",
    paste(missing_functions, collapse = ", "),
    "\nRun 01_utils.R first.",
    call. = FALSE
  )
}

message("\n============================================================")
message("05: Building all-Hsp60/10 client tables")
message("Braak model CERAD adjustment enabled: ", cfg$adjust_braak_for_cerad)
message("============================================================")

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

make_stage_meta <- function(meta, sample_col, stage_col) {
  require_columns(
    meta,
    c(sample_col, stage_col),
    "canonical clinical-stage metadata"
  )

  out <- meta |>
    dplyr::transmute(
      sample_id_for_model = as.character(.data[[sample_col]]),
      stage = factor(
        as.character(.data[[stage_col]]),
        levels = c("NCI", "MCI", "AD")
      )
    ) |>
    dplyr::filter(
      !is.na(.data$sample_id_for_model),
      nzchar(.data$sample_id_for_model),
      !is.na(.data$stage)
    )

  if (anyDuplicated(out$sample_id_for_model) > 0) {
    stop("Canonical stage metadata contain duplicate sample IDs.", call. = FALSE)
  }

  out
}

build_stage_summary <- function(mat, meta, genes, prefix) {
  genes <- intersect(clean_gene_symbols(genes), rownames(mat))

  purrr::map_dfr(genes, function(g) {
    samples <- intersect(meta$sample_id_for_model, colnames(mat))
    meta2 <- meta[match(samples, meta$sample_id_for_model), , drop = FALSE]
    y <- as.numeric(mat[g, samples])
    tmp <- tibble::tibble(stage = meta2$stage, y = y)

    nci <- tmp$y[tmp$stage == "NCI"]
    mci <- tmp$y[tmp$stage == "MCI"]
    ad  <- tmp$y[tmp$stage == "AD"]

    tibble::tibble(
      gene = g,
      "{prefix}_nci_mean"     := safe_mean(nci),
      "{prefix}_mci_mean"     := safe_mean(mci),
      "{prefix}_ad_mean"      := safe_mean(ad),
      "{prefix}_nci_sem"      := safe_sem(nci),
      "{prefix}_mci_sem"      := safe_sem(mci),
      "{prefix}_ad_sem"       := safe_sem(ad),
      "{prefix}_n_nci"        := sum(is.finite(nci)),
      "{prefix}_n_mci"        := sum(is.finite(mci)),
      "{prefix}_n_ad"         := sum(is.finite(ad)),
      "{prefix}_early_effect" := safe_mean(mci) - safe_mean(nci),
      "{prefix}_late_effect"  := safe_mean(ad) - safe_mean(mci),
      "{prefix}_total_effect" := safe_mean(ad) - safe_mean(nci)
    )
  })
}


assign_informative_function <- function(gene) {
  dplyr::case_when(
    gene %in% c("TUFM", "TSFM", "GFM1", "GFM2", "DAP3", "LRPPRC", "MRPS23") ~ "Mitochondrial translation",
    stringr::str_detect(gene, "^MRPL|^MRPS") ~ "Mitochondrial translation",
    stringr::str_detect(gene, "^NDUF|^UQCR|^COX|^SDH") ~ "OXPHOS / ETC assembly",
    stringr::str_detect(gene, "^ATP5|MT-ATP") ~ "ATP synthase / energy coupling",
    gene %in% c("PDHA1", "PDHB", "PDHX", "DLAT", "DLD", "DLST", "OGDH", "IDH2", "IDH3A", "MDH2", "CS", "FH", "ACO2") ~ "Pyruvate / TCA metabolism",
    gene %in% c("ECHS1", "ECH1", "HADHA", "HADHB", "ACAA2", "ACADVL", "ETFA", "ETFB", "ETFDH") ~ "FAO / amino-acid catabolism",
    gene %in% c("SHMT2", "MTHFD1L", "MTHFD2", "GLDC", "GOT2", "OAT", "HSD17B10") ~ "One-carbon / amino-acid metabolism",
    gene %in% c("TRAP1", "HSPA9", "GRPEL1", "DNAJA3", "LONP1", "CLPP", "CLPX", "PMPCB", "PMPCA") ~ "Mitochondrial proteostasis",
    gene %in% c("PRDX3", "PRDX4", "SOD2", "GLRX5", "GSR", "GSTK1") ~ "Redox / Fe-S / heme biology",
    stringr::str_detect(gene, "^TIMM|^TOMM|^SLC25") ~ "Mitochondrial membrane / import",
    TRUE ~ "Other / unassigned client"
  )
}

load_agora_targets <- function() {
  if (!"agora_target_file" %in% names(cfg)) {
    warning("cfg$agora_target_file is not defined. AGORA support will be set to 0.", call. = FALSE)
    return(tibble::tibble(gene = character(), agora_nominations = numeric(), agora_nominated_target = logical()))
  }

  if (!file.exists(cfg$agora_target_file)) {
    warning("AGORA target file not found: ", cfg$agora_target_file, ". AGORA support will be set to 0.", call. = FALSE)
    return(tibble::tibble(gene = character(), agora_nominations = numeric(), agora_nominated_target = logical()))
  }

  agora_raw <- readr::read_csv(cfg$agora_target_file, show_col_types = FALSE)
  message("AGORA raw dimensions: ", paste(dim(agora_raw), collapse = " x "))
  message("AGORA raw columns: ", paste(colnames(agora_raw), collapse = ", "))

  gene_col_candidates <- c(
    "Gene Symbol", "gene symbol", "GeneSymbol", "GENE_SYMBOL",
    "gene_symbol", "Gene", "gene", "SYMBOL", "Symbol", "symbol",
    "hgnc_symbol", "HGNC Symbol", "target", "Target"
  )

  agora_gene_col <- intersect(gene_col_candidates, colnames(agora_raw))[1]
  if (is.na(agora_gene_col)) {
    stop(
      "Could not identify AGORA gene-symbol column. Available columns: ",
      paste(colnames(agora_raw), collapse = ", "),
      call. = FALSE
    )
  }

  nomination_col <- intersect(c("Nominations", "nominations", "Nomination", "nomination_count"), colnames(agora_raw))[1]

  message("Using AGORA gene column: ", agora_gene_col)
  if (!is.na(nomination_col)) message("Using AGORA nominations column: ", nomination_col)

  out <- agora_raw |>
    dplyr::transmute(
      gene = canonical_gene_symbol(.data[[agora_gene_col]]),
      agora_nominations = if (!is.na(nomination_col)) safe_num(.data[[nomination_col]]) else 1,
      agora_nominated_target = TRUE
    ) |>
    dplyr::filter(!is.na(gene), gene != "") |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      agora_nominations = max(agora_nominations, na.rm = TRUE),
      agora_nominated_target = TRUE,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      agora_nominations = dplyr::if_else(is.finite(agora_nominations), agora_nominations, 1)
    )

  message("AGORA nominated target genes loaded: ", nrow(out))
  out
}

############################################################
## 1. Stage metadata and stage summaries from adjusted matrices
############################################################

prot_stage_meta <- make_stage_meta(
  prot_meta_adj,
  sample_col = "SampleID",
  stage_col = "clinical_stage"
) |>
  dplyr::filter(.data$sample_id_for_model %in% colnames(prot_mat))

rna_stage_meta <- make_stage_meta(
  rna_meta_adj,
  sample_col = "sample_id",
  stage_col = "clinical_stage"
) |>
  dplyr::filter(.data$sample_id_for_model %in% colnames(rna_mat))

stage_analysis_sample_counts <- dplyr::bind_rows(
  prot_stage_meta |>
    dplyr::count(.data$stage, name = "n") |>
    dplyr::mutate(modality = "Protein"),
  rna_stage_meta |>
    dplyr::count(.data$stage, name = "n") |>
    dplyr::mutate(modality = "RNA")
) |>
  dplyr::select(modality, stage, n)

write_tbl(
  stage_analysis_sample_counts,
  "clinical_stage_sample_counts_before_abundance_missingness"
)

all_detected_or_supplied <- clean_gene_symbols(all_hsp60_10_clients)
prot_genes <- intersect(all_detected_or_supplied, rownames(prot_mat))
rna_genes  <- intersect(all_detected_or_supplied, rownames(rna_mat))

## Pathology beta models use the processed TMT matrix before additional
## covariate residualization, with covariates included directly in each
## gene-level regression model. This avoids residualizing and then
## adjusting for the same protein covariates a second time.
prot_genes_pathology <- Reduce(intersect, list(all_detected_or_supplied, rownames(prot_mat), rownames(prot_mat_raw)))

message("Hsp60/10 clients supplied: ", length(all_detected_or_supplied))
message("Detected in adjusted protein matrix: ", length(prot_genes))
message("Detected in adjusted RNA matrix: ", length(rna_genes))
message("Detected in processed protein matrix for pathology models: ", length(prot_genes_pathology))

protein_stage_tbl <- build_stage_summary(prot_mat, prot_stage_meta, all_detected_or_supplied, "protein")
rna_stage_tbl     <- build_stage_summary(rna_mat, rna_stage_meta, all_detected_or_supplied, "rna")

############################################################
## 2. Adjusted pathology effects
##
## Use prot_mat_raw here: pathology effects are estimated from
## processed protein abundance with protein covariates included
## directly in the gene-level model. Stage/collapse and centrality
## continue to use covariate-residualized prot_mat.
############################################################

adjusted_braak_tbl <- fit_adjusted_braak_beta(
  mat = prot_mat_raw,
  meta_df = prot_meta_adj,
  sample_col = "SampleID",
  genes = prot_genes_pathology,
  covars = protein_covars,
  adjust_for_cerad = cfg$adjust_braak_for_cerad,
  min_n = cfg$min_n_gene_model
)

adjusted_cerad_tbl <- fit_adjusted_cerad_beta(
  mat = prot_mat_raw,
  meta_df = prot_meta_adj,
  sample_col = "SampleID",
  genes = prot_genes_pathology,
  covars = protein_covars,
  min_n = cfg$min_n_gene_model
)

############################################################
## 3. Network centrality from adjusted protein matrix
############################################################

if (length(prot_genes) < 3) {
  stop("Too few Hsp60/10 client proteins detected to compute centrality.", call. = FALSE)
}

prot_cor_mat_all <- suppressWarnings(
  stats::cor(
    t(prot_mat[prot_genes, , drop = FALSE]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
)
diag(prot_cor_mat_all) <- NA_real_

centrality_tbl <- tibble::tibble(
  gene = colnames(prot_cor_mat_all),
  hub_mean_abs_cor = apply(abs(prot_cor_mat_all), 2, function(x) mean(x, na.rm = TRUE)),
  hub_degree = apply(abs(prot_cor_mat_all), 2, function(x) sum(x >= 0.30, na.rm = TRUE))
)

############################################################
## 4. AGORA support
############################################################

agora_tbl <- load_agora_targets()
message("Overlap between AGORA and supplied Hsp60/10 clients: ", length(intersect(all_detected_or_supplied, agora_tbl$gene)))
if (length(intersect(all_detected_or_supplied, agora_tbl$gene)) > 0) {
  message("Overlapping AGORA/Hsp60 genes: ", paste(sort(intersect(all_detected_or_supplied, agora_tbl$gene)), collapse = ", "))
}

############################################################
## 5. Source-of-truth all-client table
############################################################

all_hsp60_10_client_tbl <- tibble::tibble(gene = all_detected_or_supplied) |>
  dplyr::mutate(
    detected_in_protein = gene %in% prot_genes,
    detected_in_rna = gene %in% rna_genes,
    display_gene = display_gene_symbol(gene)
  ) |>
  checked_left_join(protein_stage_tbl, by = "gene", label = "all clients to protein stage summary") |>
  checked_left_join(rna_stage_tbl, by = "gene", label = "all clients to RNA stage summary") |>
  checked_left_join(adjusted_braak_tbl, by = "gene", label = "all clients to adjusted Braak table") |>
  checked_left_join(adjusted_cerad_tbl, by = "gene", label = "all clients to adjusted CERAD table") |>
  checked_left_join(centrality_tbl, by = "gene", label = "all clients to centrality table") |>
  checked_left_join(agora_tbl, by = "gene", label = "all clients to AGORA table") |>
  dplyr::mutate(
    agora_nominated_target = dplyr::coalesce(agora_nominated_target, FALSE),
    agora_target = as.integer(agora_nominated_target),
    agora_flag = agora_target,
    agora_nominations = dplyr::coalesce(agora_nominations, 0),

    functional_class = assign_informative_function(gene),

    protein_late_decline_magnitude = dplyr::case_when(
      !is.finite(protein_late_effect) ~ NA_real_,
      protein_late_effect < 0 ~ abs(protein_late_effect),
      TRUE ~ 0
    ),
    late_effect_absolute_difference = abs(
      protein_late_effect - rna_late_effect
    ),

    ########################################################
    ## Canonical pathology variables
    ##
    ## Positive pathology-aligned beta means lower protein
    ## abundance with WORSE pathology for both endpoints.
    ##
    ## Braak:
    ##   higher score = worse pathology
    ##   therefore pathology-aligned beta = -braak_beta
    ##
    ## CERAD:
    ##   lower score = worse pathology
    ##   therefore pathology-aligned beta = +cerad_beta
    ########################################################

    braak_pathology_aligned_beta = dplyr::if_else(
      is.finite(braak_beta),
      -braak_beta,
      NA_real_
    ),

    cerad_pathology_aligned_beta = dplyr::if_else(
      is.finite(cerad_beta),
      cerad_beta,
      NA_real_
    ),

    ## One-sided vulnerability magnitudes.
    braak_pathology_magnitude = dplyr::if_else(
      is.finite(braak_pathology_aligned_beta),
      pmax(braak_pathology_aligned_beta, 0),
      NA_real_
    ),

    cerad_pathology_magnitude = dplyr::if_else(
      is.finite(cerad_pathology_aligned_beta),
      pmax(cerad_pathology_aligned_beta, 0),
      NA_real_
    ),

    ## Raw standardized-coefficient joint pathology metrics.
    ## These require both pathology models to be available.
    joint_pathology_magnitude = dplyr::if_else(
      is.finite(braak_pathology_magnitude) &
        is.finite(cerad_pathology_magnitude),
      (
        braak_pathology_magnitude +
          cerad_pathology_magnitude
      ) / 2,
      NA_real_
    ),

    strict_joint_pathology_magnitude = dplyr::if_else(
      is.finite(braak_pathology_magnitude) &
        is.finite(cerad_pathology_magnitude),
      pmin(
        braak_pathology_magnitude,
        cerad_pathology_magnitude
      ),
      NA_real_
    ),

    late_decline_percentile = percentile01(
      protein_late_decline_magnitude
    ),

    braak_pathology_percentile = percentile01(
      braak_pathology_magnitude
    ),

    cerad_pathology_percentile = percentile01(
      cerad_pathology_magnitude
    ),

    ## Equal-weight single AD-neuropathology score.
    ##
    ## This is the MEAN OF TWO PERCENTILES and therefore is
    ## on a 0-100 scale, but is not itself a percentile rank.
    joint_pathology_score = dplyr::if_else(
      is.finite(braak_pathology_percentile) &
        is.finite(cerad_pathology_percentile),
      (
        braak_pathology_percentile +
          cerad_pathology_percentile
      ) / 2,
      NA_real_
    ),

    ## Strict continuous convergence score: a gene is limited
    ## by its weaker pathology dimension.
    strict_joint_pathology_score = dplyr::if_else(
      is.finite(braak_pathology_percentile) &
        is.finite(cerad_pathology_percentile),
      pmin(
        braak_pathology_percentile,
        cerad_pathology_percentile
      ),
      NA_real_
    ),

    ## Proper percentile ranks of the combined scores.
    joint_pathology_percentile = percentile01(
      joint_pathology_score
    ),

    strict_joint_pathology_percentile = percentile01(
      strict_joint_pathology_score
    ),

    ## Descriptive binary concordance criterion.
    dual_pathology_topq = dplyr::case_when(
      is.finite(braak_pathology_percentile) &
        is.finite(cerad_pathology_percentile) ~
        as.integer(
          braak_pathology_percentile >= 75 &
            cerad_pathology_percentile >= 75
        ),
      TRUE ~ NA_integer_
    ),

    ########################################################
    ## Backward-compatible aliases
    ##
    ## Retain temporarily while downstream scripts are
    ## migrated. New code should use the canonical names.
    ########################################################

    inverse_braak_magnitude =
      braak_pathology_magnitude,

    inverse_cerad_magnitude =
      cerad_pathology_magnitude,

    inverse_braak_percentile =
      braak_pathology_percentile,

    inverse_cerad_percentile =
      cerad_pathology_percentile,

    centrality_percentile = percentile01(
      hub_mean_abs_cor
    ),

    late_effect_difference_percentile = percentile01(
      late_effect_absolute_difference
    ),

    ## Legacy Braak-only score retained solely for migration
    ## validation. Do not use as the revised primary pathology
    ## framework.
    legacy_braak_pathology_vulnerability_score = sqrt(
      late_decline_percentile *
        inverse_braak_percentile
    ),

    pathology_vulnerability_score =
      legacy_braak_pathology_vulnerability_score,

    ## Generic four-axis score used by some intermediate scripts.
    ## The RNA-protein late-effect difference is descriptive only.
    priority_score_internal_only = rowMeans(
      cbind(
        late_decline_percentile,
        inverse_braak_percentile,
        centrality_percentile,
        late_effect_difference_percentile
      ),
      na.rm = TRUE
    ),
    vulnerability_rank = rank(-priority_score_internal_only, ties.method = "min", na.last = "keep"),

    braak_metric = "covariate-adjusted inverse Braak beta",

    ## Backward-compatible columns required by exact migrated Figure 3 code.
    braak_rho = -inverse_braak_magnitude,
    rho_braak = braak_rho,
    cerad_rho = -inverse_cerad_magnitude,
    rho_cerad = cerad_rho,
    p_braak = braak_p,
    p_cerad = cerad_p,
    n_braak = braak_n,
    n_cerad = cerad_n,
    padj_braak = braak_padj,
    padj_cerad = cerad_padj
  ) |>
  dplyr::arrange(vulnerability_rank, gene)

message("AGORA target counts in all_hsp60_10_client_tbl:")
print(table(all_hsp60_10_client_tbl$agora_target, useNA = "ifany"))

############################################################
## 6. Long-format stage table for Figure 2-style trajectories
############################################################

all_client_stage_long <- dplyr::bind_rows(
  all_hsp60_10_client_tbl |>
    dplyr::filter(.data$detected_in_protein) |>
    dplyr::select(
      gene,
      protein_nci_mean,
      protein_mci_mean,
      protein_ad_mean,
      protein_nci_sem,
      protein_mci_sem,
      protein_ad_sem
    ) |>
    tidyr::pivot_longer(
      cols = c(protein_nci_mean, protein_mci_mean, protein_ad_mean),
      names_to = "stage_key",
      values_to = "mean_expr"
    ) |>
    dplyr::mutate(
      Modality = "Protein",
      Stage = dplyr::recode(
        .data$stage_key,
        protein_nci_mean = "NCI",
        protein_mci_mean = "MCI",
        protein_ad_mean = "AD"
      ),
      sem = dplyr::case_when(
        .data$stage_key == "protein_nci_mean" ~ .data$protein_nci_sem,
        .data$stage_key == "protein_mci_mean" ~ .data$protein_mci_sem,
        .data$stage_key == "protein_ad_mean" ~ .data$protein_ad_sem,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::select(gene, Modality, Stage, mean_expr, sem),

  all_hsp60_10_client_tbl |>
    dplyr::filter(.data$detected_in_rna) |>
    dplyr::select(
      gene,
      rna_nci_mean,
      rna_mci_mean,
      rna_ad_mean,
      rna_nci_sem,
      rna_mci_sem,
      rna_ad_sem
    ) |>
    tidyr::pivot_longer(
      cols = c(rna_nci_mean, rna_mci_mean, rna_ad_mean),
      names_to = "stage_key",
      values_to = "mean_expr"
    ) |>
    dplyr::mutate(
      Modality = "RNA",
      Stage = dplyr::recode(
        .data$stage_key,
        rna_nci_mean = "NCI",
        rna_mci_mean = "MCI",
        rna_ad_mean = "AD"
      ),
      sem = dplyr::case_when(
        .data$stage_key == "rna_nci_mean" ~ .data$rna_nci_sem,
        .data$stage_key == "rna_mci_mean" ~ .data$rna_mci_sem,
        .data$stage_key == "rna_ad_mean" ~ .data$rna_ad_sem,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::select(gene, Modality, Stage, mean_expr, sem)
) |>
  dplyr::mutate(
    Stage = factor(.data$Stage, levels = c("NCI", "MCI", "AD")),
    Modality = factor(.data$Modality, levels = c("Protein", "RNA"))
  )

############################################################
## 7. Priority table required by Figure 6
############################################################

priority_input_tbl <- all_hsp60_10_client_tbl |>
  dplyr::filter(detected_in_protein)

message("AGORA-like columns in priority_input_tbl: ", paste(grep("agora|target|nom", colnames(priority_input_tbl), value = TRUE, ignore.case = TRUE), collapse = ", "))

priority_tbl <- priority_input_tbl |>
  dplyr::mutate(
    late_decline_topq = as.integer(
      late_decline_percentile >= 75
    ),

    ## Legacy Braak-only flag retained during migration.
    braak_topq = as.integer(
      inverse_braak_percentile >= 75
    ),

    ## Canonical pathology support flags.
    cerad_topq = as.integer(
      cerad_pathology_percentile >= 75
    ),

    joint_pathology_topq = as.integer(
      joint_pathology_percentile >= 75
    ),

    strict_joint_pathology_topq = as.integer(
      strict_joint_pathology_percentile >= 75
    ),

    centrality_topq = as.integer(
      centrality_percentile >= 75
    ),

    ## LEGACY classification field retained temporarily for
    ## exact migration validation. It double-uses pathology
    ## in the old priority-score architecture and will not be
    ## the revised Figure 6 primary framework.
    clinical_topq = as.integer(
      late_decline_percentile >= 75 |
        inverse_braak_percentile >= 75
    ),

    n_total_axes = late_decline_topq + braak_topq + centrality_topq + clinical_topq + agora_target,

    proteomic_vulnerability_axis = rowMeans(
      cbind(late_decline_percentile, inverse_braak_percentile, centrality_percentile),
      na.rm = TRUE
    ),

    clinical_external_support_axis = rowMeans(
      cbind(clinical_topq * 100, agora_target * 100),
      na.rm = TRUE
    ),

    ## Main Figure 6 priority score.
    priority_score = rowMeans(
      cbind(
        late_decline_percentile,
        inverse_braak_percentile,
        centrality_percentile,
        clinical_topq * 100,
        agora_target * 100
      ),
      na.rm = TRUE
    ),

    nomination_class = dplyr::case_when(
      agora_target == 1 & priority_score >= 75 ~ "AGORA-supported high-priority",
      agora_target == 0 & priority_score >= 75 ~ "Newly nominated high-priority",
      priority_score >= 60 ~ "Emerging candidate",
      agora_target == 1 ~ "AGORA-supported contextual",
      TRUE ~ "Lower-priority / background"
    )
  ) |>
  dplyr::arrange(dplyr::desc(priority_score))

message("Priority table built successfully: ", nrow(priority_tbl), " rows x ", ncol(priority_tbl), " columns")
message("AGORA target counts in priority_tbl:")
print(table(priority_tbl$agora_target, useNA = "ifany"))
message("Nomination class counts:")
print(table(priority_tbl$nomination_class, useNA = "ifany"))

############################################################
## 8. Save source-of-truth outputs
############################################################

write_tbl(all_hsp60_10_client_tbl, "all_hsp60_10_client_tbl_COVARIATE_ADJUSTED")
write_tbl(all_client_stage_long, "all_hsp60_10_client_stage_long_COVARIATE_ADJUSTED")
write_tbl(priority_tbl, "priority_tbl_COVARIATE_ADJUSTED")

save_obj(all_hsp60_10_client_tbl, "all_hsp60_10_client_tbl")
save_obj(priority_tbl, "priority_tbl")
save_obj(all_client_stage_long, "all_client_stage_long")

message("Loaded 05_build_adjusted_all_client_tables.R")
message("============================================================\n")
