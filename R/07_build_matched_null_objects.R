############################################################
## 07_build_matched_null_objects.R
## Build abundance-matched mitochondrial null objects from the
## corrected full proteomics cohort.
##
## Key invariants:
## - Late-stage decline uses adjusted protein abundance and the canonical
##   cogdx-derived NCI/MCI/AD stage.
## - Abundance matching uses processed non-residualized protein abundance.
## - Braak effects use processed abundance with covariates in the model.
## - No AsymAD-to-MCI mapping and no "collapse" terminology.
############################################################

require_objects(
  c(
    "cfg",
    "all_hsp60_10_client_tbl",
    "prot_mat",
    "prot_mat_raw",
    "prot_meta_adj",
    "protein_covars",
    "mitocarta_genes"
  ),
  context = "07_build_matched_null_objects.R"
)

if (!file.exists(cfg$agora_target_file)) {
  stop(
    "Missing AGORA nominated target file required for matched null: ",
    cfg$agora_target_file,
    call. = FALSE
  )
}

############################################################
## 1. Gene universes
############################################################

agora_targets <- readr::read_csv(
  cfg$agora_target_file,
  show_col_types = FALSE
) |>
  janitor::clean_names()

gene_col_agora <- intersect(
  c(
    "gene",
    "gene_symbol",
    "symbol",
    "hgnc_symbol",
    "target",
    "target_gene",
    "gene_name"
  ),
  colnames(agora_targets)
)[1]

if (is.na(gene_col_agora)) {
  stop(
    "Could not find gene column in AGORA target file. Available columns: ",
    available_cols_msg(agora_targets),
    call. = FALSE
  )
}

agora_target_symbols <- agora_targets |>
  dplyr::transmute(
    gene = canonical_gene_symbol(.data[[gene_col_agora]])
  ) |>
  dplyr::filter(!is.na(.data$gene), .data$gene != "") |>
  dplyr::distinct(.data$gene) |>
  dplyr::pull(.data$gene)

hsp_symbols <- clean_gene_symbols(all_hsp60_10_client_tbl$gene)

hsp_symbols_detected <- intersect(
  hsp_symbols,
  rownames(prot_mat)
)

background_symbols <- setdiff(
  intersect(
    clean_gene_symbols(mitocarta_genes),
    rownames(prot_mat)
  ),
  hsp_symbols
)

all_null_genes <- union(
  hsp_symbols_detected,
  background_symbols
)

############################################################
## 2. Canonical clinical-stage metadata
############################################################

require_columns(
  prot_meta_adj,
  c("SampleID", "clinical_stage"),
  "matched-null protein metadata"
)

prot_stage_meta_null <- prot_meta_adj |>
  dplyr::transmute(
    SampleID = as.character(.data$SampleID),
    stage = factor(
      as.character(.data$clinical_stage),
      levels = c("NCI", "MCI", "AD")
    )
  ) |>
  dplyr::filter(
    .data$SampleID %in% colnames(prot_mat),
    !is.na(.data$stage)
  )

if (anyDuplicated(prot_stage_meta_null$SampleID) > 0) {
  stop("Duplicate SampleID in matched-null stage metadata.", call. = FALSE)
}

############################################################
## 3. Stage-associated late decline and matching abundance
############################################################

calc_background_metric_one <- function(g) {
  adjusted_y <- as.numeric(
    prot_mat[g, prot_stage_meta_null$SampleID]
  )

  raw_y <- as.numeric(
    prot_mat_raw[g, prot_stage_meta_null$SampleID]
  )

  tmp <- tibble::tibble(
    stage = prot_stage_meta_null$stage,
    adjusted_y = adjusted_y,
    raw_y = raw_y
  )

  mci <- tmp$adjusted_y[tmp$stage == "MCI"]
  ad <- tmp$adjusted_y[tmp$stage == "AD"]

  protein_late_effect <-
    mean(ad, na.rm = TRUE) -
    mean(mci, na.rm = TRUE)

  protein_late_decline_magnitude <- dplyr::case_when(
    !is.finite(protein_late_effect) ~ NA_real_,
    protein_late_effect < 0 ~ abs(protein_late_effect),
    TRUE ~ 0
  )

  tibble::tibble(
    gene = g,
    protein_late_effect = protein_late_effect,
    protein_late_decline_magnitude =
      protein_late_decline_magnitude,
    matching_abundance = mean(tmp$raw_y, na.rm = TRUE)
  )
}

background_metrics <- purrr::map_dfr(
  all_null_genes,
  calc_background_metric_one
)

############################################################
## 4. Adjusted Braak and CERAD pathology associations
############################################################

## Braak:
## higher score = worse pathology.
## Negative abundance beta therefore indicates protein loss
## with worse pathology.

background_braak <- fit_adjusted_braak_beta(
  mat = prot_mat_raw,
  meta_df = prot_meta_adj,
  sample_col = "SampleID",
  genes = all_null_genes,
  covars = protein_covars,
  adjust_for_cerad = cfg$adjust_braak_for_cerad,
  min_n = cfg$min_n_gene_model
) |>
  dplyr::transmute(
    gene,
    braak_pathology_magnitude =
      adjusted_inverse_braak_beta,

    ## Legacy alias retained for exact migration validation.
    inverse_braak_magnitude =
      adjusted_inverse_braak_beta
  )

## CERAD:
## lower score = worse pathology in this ROSMAP coding.
## Positive abundance beta therefore indicates protein loss
## with worse pathology.

background_cerad <- fit_adjusted_cerad_beta(
  mat = prot_mat_raw,
  meta_df = prot_meta_adj,
  sample_col = "SampleID",
  genes = all_null_genes,
  covars = protein_covars,
  min_n = cfg$min_n_gene_model
) |>
  dplyr::transmute(
    gene,
    cerad_pathology_magnitude =
      adjusted_inverse_cerad_beta
  )

############################################################
## 5. Combined null input
############################################################

null_input_tbl <- tibble::tibble(
  gene = all_null_genes
) |>
  dplyr::mutate(
    is_hsp60_10_client = .data$gene %in% hsp_symbols_detected
  ) |>
  checked_left_join(
    background_metrics,
    by = "gene",
    label = "null genes to corrected late-decline metrics"
  ) |>
  checked_left_join(
    background_braak,
    by = "gene",
    label = "null genes to adjusted Braak metrics"
  ) |>
  checked_left_join(
    background_cerad,
    by = "gene",
    label = "null genes to adjusted CERAD metrics"
  ) |>
  dplyr::mutate(
    agora_nominated_target =
      .data$gene %in% agora_target_symbols,

    ########################################################
    ## Rank-based legacy quantities
    ########################################################

    late_decline_rank =
      percentile01(
        .data$protein_late_decline_magnitude
      ) / 100,

    braak_pathology_rank =
      percentile01(
        .data$braak_pathology_magnitude
      ) / 100,

    ## Legacy alias.
    inverse_braak_rank =
      .data$braak_pathology_rank,

    cerad_pathology_rank =
      percentile01(
        .data$cerad_pathology_magnitude
      ) / 100,

    ########################################################
    ## Primary raw standardized joint-pathology metrics
    ########################################################

    joint_pathology_magnitude =
      (
        .data$braak_pathology_magnitude +
          .data$cerad_pathology_magnitude
      ) / 2,

    strict_joint_pathology_magnitude =
      pmin(
        .data$braak_pathology_magnitude,
        .data$cerad_pathology_magnitude
      ),

    ########################################################
    ## Rank-based joint pathology sensitivity quantities
    ########################################################

    joint_pathology_score =
      (
        .data$braak_pathology_rank +
          .data$cerad_pathology_rank
      ) / 2,

    strict_joint_pathology_score =
      pmin(
        .data$braak_pathology_rank,
        .data$cerad_pathology_rank
      ),

    ########################################################
    ## LEGACY Figure 5 composite.
    ##
    ## Retained only so the frozen null can be reproduced
    ## exactly during migration. It will no longer be the
    ## revised primary Figure 5 pathology metric.
    ########################################################

    pathology_vulnerability_score =
      rowMeans(
        cbind(
          .data$late_decline_rank,
          .data$inverse_braak_rank
        ),
        na.rm = TRUE
      )
  )

message(
  "Null input dimensions before complete-metric filter: ",
  paste(dim(null_input_tbl), collapse = " x ")
)

null_input_tbl <- null_input_tbl |>
  dplyr::filter(
    !is.na(.data$protein_late_decline_magnitude),
    !is.na(.data$inverse_braak_magnitude),
    !is.na(.data$cerad_pathology_magnitude),
    !is.na(.data$joint_pathology_magnitude),
    !is.na(.data$strict_joint_pathology_magnitude),
    !is.na(.data$pathology_vulnerability_score),
    !is.na(.data$matching_abundance)
  )

message(
  "Null input dimensions after complete-metric filter: ",
  paste(dim(null_input_tbl), collapse = " x ")
)

hsp_null_tbl <- null_input_tbl |>
  dplyr::filter(.data$is_hsp60_10_client)

background_null_pool <- null_input_tbl |>
  dplyr::filter(!.data$is_hsp60_10_client)

cat("\nNull input summary:\n")
cat("Hsp60/10 clients:", nrow(hsp_null_tbl), "\n")
cat(
  "Background non-Hsp mitochondrial proteins:",
  nrow(background_null_pool),
  "\n"
)
cat(
  "AGORA targets in Hsp60/10:",
  sum(hsp_null_tbl$agora_nominated_target),
  "\n"
)
cat(
  "AGORA targets in background:",
  sum(background_null_pool$agora_nominated_target),
  "\n"
)

if (nrow(hsp_null_tbl) < 50) {
  stop("Too few Hsp60/10 clients for null.", call. = FALSE)
}

if (nrow(background_null_pool) < 50) {
  stop(
    "Too few non-Hsp mitochondrial background proteins for null.",
    call. = FALSE
  )
}

############################################################
## 6. Verify Hsp60/10 late effects against source-of-truth table
############################################################

hsp_metric_agreement <- hsp_null_tbl |>
  dplyr::select(
    gene,
    null_late_effect = protein_late_effect,
    null_late_decline = protein_late_decline_magnitude
  ) |>
  dplyr::inner_join(
    all_hsp60_10_client_tbl |>
      dplyr::filter(.data$detected_in_protein) |>
      dplyr::select(
        gene,
        source_late_effect = protein_late_effect,
        source_late_decline = protein_late_decline_magnitude
      ),
    by = "gene"
  ) |>
  dplyr::mutate(
    abs_effect_difference = abs(
      .data$null_late_effect - .data$source_late_effect
    ),
    abs_decline_difference = abs(
      .data$null_late_decline - .data$source_late_decline
    )
  )

if (
  nrow(hsp_metric_agreement) != length(hsp_symbols_detected) ||
  max(hsp_metric_agreement$abs_effect_difference, na.rm = TRUE) > 1e-10 ||
  max(hsp_metric_agreement$abs_decline_difference, na.rm = TRUE) > 1e-10
) {
  stop(
    "Matched-null Hsp60/10 late-decline metrics do not reproduce the ",
    "source-of-truth client table.",
    call. = FALSE
  )
}

write_tbl(
  hsp_metric_agreement,
  "figure5_hsp_late_decline_source_agreement"
)

############################################################
## 7. Abundance-matched mitochondrial null
############################################################

set.seed(1300)

n_iter <- cfg$null_n_iter

null_input_tbl <- null_input_tbl |>
  dplyr::mutate(
    abundance_bin = dplyr::ntile(.data$matching_abundance, 10)
  )

hsp_null_tbl <- null_input_tbl |>
  dplyr::filter(.data$is_hsp60_10_client)

background_null_pool <- null_input_tbl |>
  dplyr::filter(!.data$is_hsp60_10_client)

n_hsp <- nrow(hsp_null_tbl)

cat("\nRunning abundance-matched mitochondrial null...\n")
cat("Iterations:", n_iter, "\n")
cat("Observed Hsp60/10 n:", n_hsp, "\n")
cat(
  "Background pool n:",
  nrow(background_null_pool),
  "\n"
)

observed_stats <- hsp_null_tbl |>
  dplyr::summarise(
    observed_mean_late_decline = mean(
      .data$protein_late_decline_magnitude,
      na.rm = TRUE
    ),
    observed_mean_inverse_braak = mean(
      .data$inverse_braak_magnitude,
      na.rm = TRUE
    ),

    observed_mean_cerad_pathology = mean(
      .data$cerad_pathology_magnitude,
      na.rm = TRUE
    ),

    observed_mean_joint_pathology = mean(
      .data$joint_pathology_magnitude,
      na.rm = TRUE
    ),

    observed_mean_strict_joint_pathology = mean(
      .data$strict_joint_pathology_magnitude,
      na.rm = TRUE
    ),

    ## Legacy Braak + late-decline composite.
    observed_mean_pathology_score = mean(
      .data$pathology_vulnerability_score,
      na.rm = TRUE
    ),
    observed_agora_fraction = mean(
      .data$agora_nominated_target,
      na.rm = TRUE
    ),
    observed_agora_n = sum(
      .data$agora_nominated_target,
      na.rm = TRUE
    ),
    n_hsp = dplyr::n()
  )

sample_abundance_matched_null <- function(hsp_tbl, bg_tbl) {
  hsp_tbl |>
    dplyr::count(.data$abundance_bin, name = "n_needed") |>
    dplyr::group_split(.data$abundance_bin) |>
    purrr::map_dfr(function(bin_df) {
      this_bin <- bin_df$abundance_bin[[1]]
      n_needed <- bin_df$n_needed[[1]]

      bg_bin <- bg_tbl |>
        dplyr::filter(.data$abundance_bin == this_bin)

      if (nrow(bg_bin) == 0) {
        stop(
          "No background proteins available in abundance bin ",
          this_bin,
          ".",
          call. = FALSE
        )
      }

      bg_bin |>
        dplyr::slice_sample(
          n = n_needed,
          replace = nrow(bg_bin) < n_needed
        )
    })
}

null_results <- purrr::map_dfr(
  seq_len(n_iter),
  function(i) {
    sampled_bg <- sample_abundance_matched_null(
      hsp_tbl = hsp_null_tbl,
      bg_tbl = background_null_pool
    )

    tibble::tibble(
      iter = i,
      null_mean_late_decline = mean(
        sampled_bg$protein_late_decline_magnitude,
        na.rm = TRUE
      ),
      null_mean_inverse_braak = mean(
        sampled_bg$inverse_braak_magnitude,
        na.rm = TRUE
      ),

      null_mean_cerad_pathology = mean(
        sampled_bg$cerad_pathology_magnitude,
        na.rm = TRUE
      ),

      null_mean_joint_pathology = mean(
        sampled_bg$joint_pathology_magnitude,
        na.rm = TRUE
      ),

      null_mean_strict_joint_pathology = mean(
        sampled_bg$strict_joint_pathology_magnitude,
        na.rm = TRUE
      ),

      ## Legacy migration-control distribution.
      null_mean_pathology_score = mean(
        sampled_bg$pathology_vulnerability_score,
        na.rm = TRUE
      ),
      null_agora_fraction = mean(
        sampled_bg$agora_nominated_target,
        na.rm = TRUE
      ),
      null_agora_n = sum(
        sampled_bg$agora_nominated_target,
        na.rm = TRUE
      ),
      null_n = nrow(sampled_bg)
    )
  }
)

obs <- observed_stats

null_summary <- tibble::tibble(
  metric = c(
    "protein_late_decline_magnitude",

    ## Legacy migration-control metrics.
    "inverse_braak_magnitude",
    "pathology_vulnerability_score",

    ## Corrected pathology metrics.
    "cerad_pathology_magnitude",
    "joint_pathology_magnitude",
    "strict_joint_pathology_magnitude",

    "agora_fraction"
  ),

  observed = c(
    obs$observed_mean_late_decline,

    obs$observed_mean_inverse_braak,
    obs$observed_mean_pathology_score,

    obs$observed_mean_cerad_pathology,
    obs$observed_mean_joint_pathology,
    obs$observed_mean_strict_joint_pathology,

    obs$observed_agora_fraction
  ),

  null_mean = c(
    mean(
      null_results$null_mean_late_decline,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_inverse_braak,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_pathology_score,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_cerad_pathology,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_joint_pathology,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_strict_joint_pathology,
      na.rm = TRUE
    ),

    mean(
      null_results$null_agora_fraction,
      na.rm = TRUE
    )
  ),

  null_sd = c(
    stats::sd(
      null_results$null_mean_late_decline,
      na.rm = TRUE
    ),

    stats::sd(
      null_results$null_mean_inverse_braak,
      na.rm = TRUE
    ),

    stats::sd(
      null_results$null_mean_pathology_score,
      na.rm = TRUE
    ),

    stats::sd(
      null_results$null_mean_cerad_pathology,
      na.rm = TRUE
    ),

    stats::sd(
      null_results$null_mean_joint_pathology,
      na.rm = TRUE
    ),

    stats::sd(
      null_results$null_mean_strict_joint_pathology,
      na.rm = TRUE
    ),

    stats::sd(
      null_results$null_agora_fraction,
      na.rm = TRUE
    )
  ),

  empirical_p_greater = c(
    mean(
      null_results$null_mean_late_decline >=
        obs$observed_mean_late_decline,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_inverse_braak >=
        obs$observed_mean_inverse_braak,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_pathology_score >=
        obs$observed_mean_pathology_score,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_cerad_pathology >=
        obs$observed_mean_cerad_pathology,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_joint_pathology >=
        obs$observed_mean_joint_pathology,
      na.rm = TRUE
    ),

    mean(
      null_results$null_mean_strict_joint_pathology >=
        obs$observed_mean_strict_joint_pathology,
      na.rm = TRUE
    ),

    mean(
      null_results$null_agora_fraction >=
        obs$observed_agora_fraction,
      na.rm = TRUE
    )
  )
) |>
  dplyr::mutate(
    z_score =
      (
        .data$observed -
          .data$null_mean
      ) /
        .data$null_sd,

    empirical_p_greater =
      pmax(
        .data$empirical_p_greater,
        1 / n_iter
      )
  )

############################################################
## 8. Save source objects
############################################################

output_dir <- file.path(
  cfg$plot_dir,
  "main_fig5_matched_null_specificity_exact"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tbl(
  null_summary,
  "figure5_null_summary_COVARIATE_ADJUSTED"
)
write_tbl(
  null_results,
  "figure5_null_results_COVARIATE_ADJUSTED"
)
write_tbl(
  null_input_tbl,
  "figure5_null_input_COVARIATE_ADJUSTED"
)

save_obj(null_results, "figure5_null_results")
save_obj(null_summary, "figure5_null_summary")
save_obj(observed_stats, "figure5_observed_stats")

message("Loaded 07_build_matched_null_objects.R")
