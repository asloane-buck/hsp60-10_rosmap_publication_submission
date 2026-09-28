############################################################
## 50_run_APOE_sensitivity.R
## Reviewer 3.6: APOE-e4 sensitivity for key Hsp60/10
## protein stage/pathology analyses.
##
## Design:
##   1) full original cohort/model
##   2) APOE-complete subset, original covariates
##   3) same APOE-complete subset + APOE-e4 carrier
##
## This separates cohort-restriction effects from APOE
## adjustment effects.
############################################################

options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R",
  "R/04_build_pathway_sets_all_clients.R"
)) {
  if (!file.exists(f)) stop("Missing required production script: ", f)
  source(f)
}

OUTDIR <- file.path("outputs", "reviewer_revisions", "APOE_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

required_objects <- c(
  "prot_mat_raw", "prot_mat", "prot_meta_adj",
  "protein_covars", "all_hsp60_10_clients"
)
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects) > 0) {
  stop("Missing required objects: ", paste(missing_objects, collapse = ", "))
}

expected_covars <- c("age_num", "sex_factor", "pmi_num", "batch_factor")
if (!identical(protein_covars, expected_covars)) {
  stop("Protein covariates changed: ", paste(protein_covars, collapse = ";"))
}
if (isTRUE(cfg$adjust_braak_for_cerad)) {
  stop("Expected cfg$adjust_braak_for_cerad = FALSE.")
}
if (!identical(cfg$min_n_gene_model, 30)) {
  stop("Expected cfg$min_n_gene_model = 30.")
}

############################################################
## APOE lookup from script 49
############################################################

apoe_file <- file.path(OUTDIR, "APOE_participant_lookup.csv")
if (!file.exists(apoe_file)) {
  stop("Missing APOE lookup from script 49: ", apoe_file)
}

apoe_lookup <- read.csv(apoe_file, stringsAsFactors = FALSE)
required_apoe_cols <- c(
  "individual_id", "apoe_genotype",
  "apoe4_carrier", "apoe4_allele_count"
)
missing_apoe_cols <- setdiff(required_apoe_cols, names(apoe_lookup))
if (length(missing_apoe_cols) > 0) {
  stop("APOE lookup missing: ", paste(missing_apoe_cols, collapse = ", "))
}
if (anyDuplicated(apoe_lookup$individual_id)) {
  stop("APOE lookup is not unique by individual_id.")
}

prot_meta_full <- prot_meta_adj |>
  dplyr::select(
    -dplyr::any_of(c(
      "apoe_genotype", "apoe4_carrier",
      "apoe4_allele_count", "apoe4_carrier_num"
    ))
  ) |>
  dplyr::left_join(
    apoe_lookup |>
      dplyr::select(
        individual_id, apoe_genotype,
        apoe4_carrier, apoe4_allele_count
      ),
    by = "individual_id"
  ) |>
  dplyr::mutate(
    apoe4_carrier_num = dplyr::case_when(
      .data$apoe4_carrier %in% TRUE ~ 1,
      .data$apoe4_carrier %in% FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    apoe4_carrier_factor = factor(
      .data$apoe4_carrier_num,
      levels = c(0, 1),
      labels = c("noncarrier", "carrier")
    )
  )

if (nrow(prot_meta_full) != 400L) {
  stop("Expected full protein metadata n=400; found ", nrow(prot_meta_full))
}
if (anyDuplicated(prot_meta_full$SampleID)) {
  stop("Protein SampleID is not unique.")
}

prot_meta_full <- prot_meta_full[
  match(colnames(prot_mat_raw), prot_meta_full$SampleID),
  ,
  drop = FALSE
]
if (
  any(is.na(prot_meta_full$SampleID)) ||
  !all(colnames(prot_mat_raw) == prot_meta_full$SampleID)
) {
  stop("Could not align protein metadata to prot_mat_raw.")
}

apoe_complete <- is.finite(prot_meta_full$apoe4_carrier_num)
prot_meta_apoe <- prot_meta_full[apoe_complete, , drop = FALSE]
prot_raw_apoe <- prot_mat_raw[, prot_meta_apoe$SampleID, drop = FALSE]

if (nrow(prot_meta_apoe) != 338L) {
  stop("Expected APOE-complete protein cohort n=338; found ", nrow(prot_meta_apoe))
}
carrier_counts <- table(prot_meta_apoe$apoe4_carrier_num)
if (!identical(as.integer(carrier_counts[c("0", "1")]), c(271L, 67L))) {
  stop("Unexpected APOE-e4 carrier counts.")
}

hsp_genes <- intersect(
  clean_gene_symbols(all_hsp60_10_clients),
  rownames(prot_mat_raw)
)
if (length(hsp_genes) != 306L) {
  stop("Expected 306 detected Hsp60/10 proteins; found ", length(hsp_genes))
}

############################################################
## Same APOE-complete subset, residualized two ways
############################################################

prot_adj_subset_control <- residualize_matrix(
  prot_raw_apoe,
  prot_meta_apoe,
  sample_col = "SampleID",
  covars = protein_covars
)

prot_adj_subset_apoe <- residualize_matrix(
  prot_raw_apoe,
  prot_meta_apoe,
  sample_col = "SampleID",
  covars = c(protein_covars, "apoe4_carrier_factor")
)

if (!identical(dim(prot_adj_subset_control), dim(prot_adj_subset_apoe))) {
  stop("Subset residualized matrices differ in dimension.")
}
if (!identical(colnames(prot_adj_subset_control), colnames(prot_adj_subset_apoe))) {
  stop("Subset residualized matrices differ in samples.")
}

make_stage_meta_local <- function(meta_df) {
  meta_df |>
    dplyr::transmute(
      sample_id = as.character(.data$SampleID),
      stage = factor(
        as.character(.data$clinical_stage),
        levels = c("NCI", "MCI", "AD")
      )
    ) |>
    dplyr::filter(!is.na(.data$stage))
}

stage_meta_full <- make_stage_meta_local(prot_meta_full)
stage_meta_apoe <- make_stage_meta_local(prot_meta_apoe)

## Phenotype-eligible APOE-complete pure-stage counts are 149/86/83.
## One APOE-complete NCI participant is nuisance-incomplete, so the
## effective residualized stage cohort is 148/86/83 = 317.
eligible_stage_counts <- table(
  factor(stage_meta_apoe$stage, levels = c("NCI", "MCI", "AD"))
)
if (!identical(as.integer(eligible_stage_counts), c(149L, 86L, 83L))) {
  stop(
    "Unexpected APOE-complete phenotype-eligible stage counts: ",
    paste(as.integer(eligible_stage_counts), collapse = "/")
  )
}

nuisance_complete_apoe <- stats::complete.cases(
  prot_meta_apoe[, c(protein_covars, "clinical_stage"), drop = FALSE]
)
effective_stage_counts <- table(
  factor(
    prot_meta_apoe$clinical_stage[nuisance_complete_apoe],
    levels = c("NCI", "MCI", "AD")
  )
)
if (!identical(as.integer(effective_stage_counts), c(148L, 86L, 83L))) {
  stop(
    "Unexpected APOE-complete nuisance-complete stage counts: ",
    paste(as.integer(effective_stage_counts), collapse = "/")
  )
}

############################################################
## Late-stage AD-MCI effects
############################################################

fit_late_stage <- function(mat, stage_meta, genes, model_name) {
  purrr::map_dfr(genes, function(g) {
    samples <- intersect(stage_meta$sample_id, colnames(mat))
    meta2 <- stage_meta[
      match(samples, stage_meta$sample_id),
      ,
      drop = FALSE
    ]
    df <- tibble::tibble(
      y = as.numeric(mat[g, samples]),
      stage = factor(as.character(meta2$stage), levels = c("MCI", "AD"))
    ) |>
      dplyr::filter(!is.na(.data$stage), is.finite(.data$y))

    n_mci <- sum(df$stage == "MCI")
    n_ad <- sum(df$stage == "AD")
    if (n_mci < 2L || n_ad < 2L) {
      return(tibble::tibble(
        gene = g, model = model_name, n_mci = n_mci, n_ad = n_ad,
        effect = NA_real_, se = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, p_value = NA_real_
      ))
    }

    mean_difference <-
      mean(df$y[df$stage == "AD"], na.rm = TRUE) -
      mean(df$y[df$stage == "MCI"], na.rm = TRUE)

    hit <- broom::tidy(stats::lm(y ~ stage, data = df), conf.int = TRUE) |>
      dplyr::filter(.data$term == "stageAD")

    if (nrow(hit) != 1L) stop("Could not extract AD-MCI effect for ", g)
    if (!isTRUE(all.equal(
      as.numeric(hit$estimate),
      mean_difference,
      tolerance = 1e-10
    ))) {
      stop("AD-MCI coefficient != production mean difference for ", g)
    }

    tibble::tibble(
      gene = g,
      model = model_name,
      n_mci = n_mci,
      n_ad = n_ad,
      effect = as.numeric(hit$estimate),
      se = as.numeric(hit$std.error),
      ci_low = as.numeric(hit$conf.low),
      ci_high = as.numeric(hit$conf.high),
      p_value = as.numeric(hit$p.value)
    )
  }) |>
    dplyr::group_by(.data$model) |>
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = "BH"),
      late_decline_magnitude = dplyr::case_when(
        !is.finite(.data$effect) ~ NA_real_,
        .data$effect < 0 ~ abs(.data$effect),
        TRUE ~ 0
      )
    ) |>
    dplyr::ungroup()
}

late_results <- dplyr::bind_rows(
  fit_late_stage(
    prot_mat, stage_meta_full, hsp_genes, "full_original"
  ),
  fit_late_stage(
    prot_adj_subset_control,
    stage_meta_apoe,
    hsp_genes,
    "APOE_complete_original_covariates"
  ),
  fit_late_stage(
    prot_adj_subset_apoe,
    stage_meta_apoe,
    hsp_genes,
    "APOE_complete_plus_e4"
  )
)

############################################################
## Direct Braak/CERAD models with SE/CI
############################################################

fit_pathology_augmented <- function(
  mat, meta_df, genes, predictor, covars,
  model_name, endpoint_name, min_n = 30
) {
  mat <- as.matrix(mat)
  meta_df <- meta_df[
    match(colnames(mat), meta_df$SampleID),
    ,
    drop = FALSE
  ]
  if (
    any(is.na(meta_df$SampleID)) ||
    !all(colnames(mat) == meta_df$SampleID)
  ) {
    stop("Pathology-model metadata alignment failed.")
  }

  covars <- intersect(covars, names(meta_df))

  purrr::map_dfr(genes, function(g) {
    y <- safe_z(mat[g, ])
    df <- tibble::tibble(
      abundance = y,
      predictor_value = meta_df[[predictor]]
    ) |>
      dplyr::bind_cols(
        meta_df |>
          dplyr::select(dplyr::all_of(covars))
      ) |>
      dplyr::mutate(
        dplyr::across(where(is.character), as.factor)
      )

    keep <- stats::complete.cases(
      df[, c("abundance", "predictor_value", covars), drop = FALSE]
    )
    n_model <- sum(keep)

    if (
      n_model < min_n ||
      length(unique(df$predictor_value[keep])) < 3L
    ) {
      return(tibble::tibble(
        gene = g, endpoint = endpoint_name, model = model_name, n = n_model,
        effect = NA_real_, se = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, p_value = NA_real_
      ))
    }

    form <- stats::as.formula(
      paste(
        "abundance ~",
        paste(c("predictor_value", covars), collapse = " + ")
      )
    )
    fit <- tryCatch(
      stats::lm(form, data = df[keep, , drop = FALSE]),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble::tibble(
        gene = g, endpoint = endpoint_name, model = model_name, n = n_model,
        effect = NA_real_, se = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, p_value = NA_real_
      ))
    }

    hit <- broom::tidy(fit, conf.int = TRUE) |>
      dplyr::filter(.data$term == "predictor_value")

    if (nrow(hit) != 1L) {
      return(tibble::tibble(
        gene = g, endpoint = endpoint_name, model = model_name, n = n_model,
        effect = NA_real_, se = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, p_value = NA_real_
      ))
    }

    tibble::tibble(
      gene = g,
      endpoint = endpoint_name,
      model = model_name,
      n = n_model,
      effect = as.numeric(hit$estimate),
      se = as.numeric(hit$std.error),
      ci_low = as.numeric(hit$conf.low),
      ci_high = as.numeric(hit$conf.high),
      p_value = as.numeric(hit$p.value)
    )
  }) |>
    dplyr::group_by(.data$endpoint, .data$model) |>
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = "BH"),
      inverse_magnitude = dplyr::case_when(
        !is.finite(.data$effect) ~ NA_real_,
        .data$effect < 0 ~ abs(.data$effect),
        TRUE ~ 0
      )
    ) |>
    dplyr::ungroup()
}

braak_results <- dplyr::bind_rows(
  fit_pathology_augmented(
    prot_mat_raw, prot_meta_full, hsp_genes, "braak_num_std",
    protein_covars, "full_original", "Braak", cfg$min_n_gene_model
  ),
  fit_pathology_augmented(
    prot_raw_apoe, prot_meta_apoe, hsp_genes, "braak_num_std",
    protein_covars, "APOE_complete_original_covariates",
    "Braak", cfg$min_n_gene_model
  ),
  fit_pathology_augmented(
    prot_raw_apoe, prot_meta_apoe, hsp_genes, "braak_num_std",
    c(protein_covars, "apoe4_carrier_factor"),
    "APOE_complete_plus_e4", "Braak", cfg$min_n_gene_model
  )
)

cerad_results <- dplyr::bind_rows(
  fit_pathology_augmented(
    prot_mat_raw, prot_meta_full, hsp_genes, "cerad_num_std",
    protein_covars, "full_original", "CERAD", cfg$min_n_gene_model
  ),
  fit_pathology_augmented(
    prot_raw_apoe, prot_meta_apoe, hsp_genes, "cerad_num_std",
    protein_covars, "APOE_complete_original_covariates",
    "CERAD", cfg$min_n_gene_model
  ),
  fit_pathology_augmented(
    prot_raw_apoe, prot_meta_apoe, hsp_genes, "cerad_num_std",
    c(protein_covars, "apoe4_carrier_factor"),
    "APOE_complete_plus_e4", "CERAD", cfg$min_n_gene_model
  )
)

############################################################
## Validate effect/P/n against production pathology helpers
############################################################

validate_pathology <- function(
  endpoint, mat, meta_df, covars, augmented_tbl, model_name
) {
  if (endpoint == "Braak") {
    prod <- fit_adjusted_braak_beta(
      mat = mat,
      meta_df = meta_df,
      sample_col = "SampleID",
      genes = hsp_genes,
      covars = covars,
      adjust_for_cerad = FALSE,
      min_n = cfg$min_n_gene_model
    ) |>
      dplyr::transmute(
        gene,
        prod_effect = braak_beta,
        prod_p = braak_p,
        prod_n = braak_n
      )
  } else {
    prod <- fit_adjusted_cerad_beta(
      mat = mat,
      meta_df = meta_df,
      sample_col = "SampleID",
      genes = hsp_genes,
      covars = covars,
      min_n = cfg$min_n_gene_model
    ) |>
      dplyr::transmute(
        gene,
        prod_effect = cerad_beta,
        prod_p = cerad_p,
        prod_n = cerad_n
      )
  }

  chk <- augmented_tbl |>
    dplyr::filter(.data$model == model_name) |>
    dplyr::select(gene, effect, p_value, n) |>
    dplyr::left_join(prod, by = "gene")

  finite_effect <- is.finite(chk$effect) & is.finite(chk$prod_effect)
  finite_p <- is.finite(chk$p_value) & is.finite(chk$prod_p)

  max_effect_diff <- if (any(finite_effect)) {
    max(abs(chk$effect[finite_effect] - chk$prod_effect[finite_effect]))
  } else {
    NA_real_
  }
  max_p_diff <- if (any(finite_p)) {
    max(abs(chk$p_value[finite_p] - chk$prod_p[finite_p]))
  } else {
    NA_real_
  }
  n_mismatch <- sum(chk$n != chk$prod_n, na.rm = TRUE)

  if (is.finite(max_effect_diff) && max_effect_diff > 1e-10) {
    stop(endpoint, " effect differs from production helper in ", model_name)
  }
  if (is.finite(max_p_diff) && max_p_diff > 1e-10) {
    stop(endpoint, " P differs from production helper in ", model_name)
  }
  if (n_mismatch > 0) {
    stop(endpoint, " n differs from production helper in ", model_name)
  }

  tibble::tibble(
    endpoint = endpoint,
    model = model_name,
    max_effect_difference = max_effect_diff,
    max_p_difference = max_p_diff,
    n_mismatches = n_mismatch
  )
}

production_validation <- dplyr::bind_rows(
  validate_pathology(
    "Braak", prot_mat_raw, prot_meta_full, protein_covars,
    braak_results, "full_original"
  ),
  validate_pathology(
    "Braak", prot_raw_apoe, prot_meta_apoe, protein_covars,
    braak_results, "APOE_complete_original_covariates"
  ),
  validate_pathology(
    "Braak", prot_raw_apoe, prot_meta_apoe,
    c(protein_covars, "apoe4_carrier_factor"),
    braak_results, "APOE_complete_plus_e4"
  ),
  validate_pathology(
    "CERAD", prot_mat_raw, prot_meta_full, protein_covars,
    cerad_results, "full_original"
  ),
  validate_pathology(
    "CERAD", prot_raw_apoe, prot_meta_apoe, protein_covars,
    cerad_results, "APOE_complete_original_covariates"
  ),
  validate_pathology(
    "CERAD", prot_raw_apoe, prot_meta_apoe,
    c(protein_covars, "apoe4_carrier_factor"),
    cerad_results, "APOE_complete_plus_e4"
  )
)

############################################################
## Derived late-decline/Braak vulnerability
############################################################

build_vulnerability <- function(model_name) {
  late_results |>
    dplyr::filter(.data$model == model_name) |>
    dplyr::select(gene, late_effect = effect, late_decline_magnitude) |>
    dplyr::inner_join(
      braak_results |>
        dplyr::filter(.data$model == model_name) |>
        dplyr::select(
          gene,
          braak_beta = effect,
          inverse_braak_magnitude = inverse_magnitude
        ),
      by = "gene"
    ) |>
    dplyr::mutate(
      late_decline_percentile = percentile01(.data$late_decline_magnitude),
      inverse_braak_percentile = percentile01(.data$inverse_braak_magnitude),
      pathology_vulnerability_score = sqrt(
        .data$late_decline_percentile *
          .data$inverse_braak_percentile
      ),
      model = model_name
    )
}

vulnerability_results <- dplyr::bind_rows(
  build_vulnerability("full_original"),
  build_vulnerability("APOE_complete_original_covariates"),
  build_vulnerability("APOE_complete_plus_e4")
)

############################################################
## Same-subset effect comparisons
############################################################

compare_effect_models <- function(tbl, endpoint_name) {
  a <- tbl |>
    dplyr::filter(
      .data$model == "APOE_complete_original_covariates"
    ) |>
    dplyr::select(
      gene,
      effect_without_APOE = effect,
      p_without_APOE = p_value,
      fdr_without_APOE = fdr
    )

  b <- tbl |>
    dplyr::filter(.data$model == "APOE_complete_plus_e4") |>
    dplyr::select(
      gene,
      effect_with_APOE = effect,
      p_with_APOE = p_value,
      fdr_with_APOE = fdr
    )

  paired <- a |>
    dplyr::inner_join(b, by = "gene") |>
    dplyr::filter(
      is.finite(.data$effect_without_APOE),
      is.finite(.data$effect_with_APOE)
    ) |>
    dplyr::mutate(
      effect_change = .data$effect_with_APOE - .data$effect_without_APOE,
      absolute_effect_change = abs(.data$effect_change),
      same_direction = sign(.data$effect_without_APOE) ==
        sign(.data$effect_with_APOE)
    )

  summary <- tibble::tibble(
    endpoint = endpoint_name,
    n_shared = nrow(paired),
    pearson_effect = stats::cor(
      paired$effect_without_APOE,
      paired$effect_with_APOE,
      method = "pearson"
    ),
    spearman_effect = stats::cor(
      paired$effect_without_APOE,
      paired$effect_with_APOE,
      method = "spearman"
    ),
    direction_concordance = mean(paired$same_direction),
    mean_effect_without_APOE = mean(paired$effect_without_APOE),
    mean_effect_with_APOE = mean(paired$effect_with_APOE),
    median_absolute_effect_change = stats::median(
      paired$absolute_effect_change
    ),
    max_absolute_effect_change = max(paired$absolute_effect_change),
    n_nominal_without_APOE = sum(
      paired$p_without_APOE < 0.05,
      na.rm = TRUE
    ),
    n_nominal_with_APOE = sum(
      paired$p_with_APOE < 0.05,
      na.rm = TRUE
    ),
    n_FDR_without_APOE = sum(
      paired$fdr_without_APOE < 0.05,
      na.rm = TRUE
    ),
    n_FDR_with_APOE = sum(
      paired$fdr_with_APOE < 0.05,
      na.rm = TRUE
    )
  )

  list(
    paired = paired |>
      dplyr::mutate(endpoint = endpoint_name, .before = 1),
    summary = summary
  )
}

late_compare <- compare_effect_models(late_results, "Late_AD_minus_MCI")
braak_compare <- compare_effect_models(braak_results, "Braak")
cerad_compare <- compare_effect_models(cerad_results, "CERAD")

paired_effect_comparisons <- dplyr::bind_rows(
  late_compare$paired,
  braak_compare$paired,
  cerad_compare$paired
)
effect_summary <- dplyr::bind_rows(
  late_compare$summary,
  braak_compare$summary,
  cerad_compare$summary
)

############################################################
## Vulnerability concordance and top-quartile retention
############################################################

vulnerability_paired <- vulnerability_results |>
  dplyr::filter(
    .data$model == "APOE_complete_original_covariates"
  ) |>
  dplyr::select(
    gene,
    late_no_APOE = late_decline_percentile,
    braak_no_APOE = inverse_braak_percentile,
    vulnerability_no_APOE = pathology_vulnerability_score
  ) |>
  dplyr::inner_join(
    vulnerability_results |>
      dplyr::filter(.data$model == "APOE_complete_plus_e4") |>
      dplyr::select(
        gene,
        late_with_APOE = late_decline_percentile,
        braak_with_APOE = inverse_braak_percentile,
        vulnerability_with_APOE = pathology_vulnerability_score
      ),
    by = "gene"
  ) |>
  dplyr::mutate(
    late_topq_no_APOE = .data$late_no_APOE >= 75,
    late_topq_with_APOE = .data$late_with_APOE >= 75,
    braak_topq_no_APOE = .data$braak_no_APOE >= 75,
    braak_topq_with_APOE = .data$braak_with_APOE >= 75
  )

vulnerability_summary <- tibble::tibble(
  metric = c(
    "late_decline_percentile",
    "inverse_braak_percentile",
    "pathology_vulnerability_score"
  ),
  spearman = c(
    stats::cor(
      vulnerability_paired$late_no_APOE,
      vulnerability_paired$late_with_APOE,
      method = "spearman",
      use = "complete.obs"
    ),
    stats::cor(
      vulnerability_paired$braak_no_APOE,
      vulnerability_paired$braak_with_APOE,
      method = "spearman",
      use = "complete.obs"
    ),
    stats::cor(
      vulnerability_paired$vulnerability_no_APOE,
      vulnerability_paired$vulnerability_with_APOE,
      method = "spearman",
      use = "complete.obs"
    )
  ),
  pearson = c(
    stats::cor(
      vulnerability_paired$late_no_APOE,
      vulnerability_paired$late_with_APOE,
      method = "pearson",
      use = "complete.obs"
    ),
    stats::cor(
      vulnerability_paired$braak_no_APOE,
      vulnerability_paired$braak_with_APOE,
      method = "pearson",
      use = "complete.obs"
    ),
    stats::cor(
      vulnerability_paired$vulnerability_no_APOE,
      vulnerability_paired$vulnerability_with_APOE,
      method = "pearson",
      use = "complete.obs"
    )
  )
)

topq_overlap <- tibble::tibble(
  axis = c("late_decline", "inverse_braak"),
  n_topq_without_APOE = c(
    sum(vulnerability_paired$late_topq_no_APOE, na.rm = TRUE),
    sum(vulnerability_paired$braak_topq_no_APOE, na.rm = TRUE)
  ),
  n_topq_with_APOE = c(
    sum(vulnerability_paired$late_topq_with_APOE, na.rm = TRUE),
    sum(vulnerability_paired$braak_topq_with_APOE, na.rm = TRUE)
  ),
  n_overlap = c(
    sum(
      vulnerability_paired$late_topq_no_APOE &
        vulnerability_paired$late_topq_with_APOE,
      na.rm = TRUE
    ),
    sum(
      vulnerability_paired$braak_topq_no_APOE &
        vulnerability_paired$braak_topq_with_APOE,
      na.rm = TRUE
    )
  )
) |>
  dplyr::mutate(
    pct_original_topq_retained =
      100 * .data$n_overlap / .data$n_topq_without_APOE
  )

############################################################
## Full original vs APOE-complete subset (no APOE adjustment)
############################################################

compare_full_to_subset <- function(tbl, endpoint_name) {
  x <- tbl |>
    dplyr::filter(.data$model == "full_original") |>
    dplyr::select(gene, effect_full = effect) |>
    dplyr::inner_join(
      tbl |>
        dplyr::filter(
          .data$model == "APOE_complete_original_covariates"
        ) |>
        dplyr::select(gene, effect_subset = effect),
      by = "gene"
    ) |>
    dplyr::filter(
      is.finite(.data$effect_full),
      is.finite(.data$effect_subset)
    )

  tibble::tibble(
    endpoint = endpoint_name,
    n_shared = nrow(x),
    pearson = stats::cor(
      x$effect_full, x$effect_subset, method = "pearson"
    ),
    spearman = stats::cor(
      x$effect_full, x$effect_subset, method = "spearman"
    ),
    direction_concordance = mean(
      sign(x$effect_full) == sign(x$effect_subset)
    )
  )
}

full_vs_subset_summary <- dplyr::bind_rows(
  compare_full_to_subset(late_results, "Late_AD_minus_MCI"),
  compare_full_to_subset(braak_results, "Braak"),
  compare_full_to_subset(cerad_results, "CERAD")
)

############################################################
## Write outputs and provenance
############################################################

cohort_audit <- tibble::tibble(
  metric = c(
    "full_protein_samples",
    "APOE_complete_protein_samples",
    "APOE_e4_noncarriers",
    "APOE_e4_carriers",
    "APOE_complete_stage_phenotype_eligible",
    "APOE_complete_stage_nuisance_complete",
    "Hsp60_10_client_proteins",
    "primary_protein_covariates",
    "APOE_sensitivity_covariate",
    "Braak_adjusted_for_CERAD"
  ),
  value = c(
    "400",
    "338",
    "271",
    "67",
    "NCI=149;MCI=86;AD=83",
    "NCI=148;MCI=86;AD=83",
    as.character(length(hsp_genes)),
    paste(protein_covars, collapse = ";"),
    "apoe4_carrier_factor",
    as.character(cfg$adjust_braak_for_cerad)
  )
)

outputs <- list(
  "APOE_sensitivity_cohort_audit.csv" = cohort_audit,
  "APOE_sensitivity_late_stage_clients.csv" = late_results,
  "APOE_sensitivity_Braak_clients.csv" = braak_results,
  "APOE_sensitivity_CERAD_clients.csv" = cerad_results,
  "APOE_sensitivity_production_model_validation.csv" =
    production_validation,
  "APOE_sensitivity_same_subset_paired_effects.csv" =
    paired_effect_comparisons,
  "APOE_sensitivity_effect_summary.csv" = effect_summary,
  "APOE_sensitivity_vulnerability_results.csv" =
    vulnerability_results,
  "APOE_sensitivity_vulnerability_paired.csv" =
    vulnerability_paired,
  "APOE_sensitivity_vulnerability_summary.csv" =
    vulnerability_summary,
  "APOE_sensitivity_top_quartile_overlap.csv" = topq_overlap,
  "APOE_sensitivity_full_vs_subset_summary.csv" =
    full_vs_subset_summary
)
for (nm in names(outputs)) {
  utils::write.csv(
    outputs[[nm]],
    file.path(OUTDIR, nm),
    row.names = FALSE
  )
}

git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

provenance <- tibble::tibble(
  item = c(
    "analysis",
    "APOE_encoding",
    "APOE_carrier_definition",
    "full_protein_n",
    "APOE_complete_n",
    "APOE_complete_stage_nuisance_complete",
    "Hsp_client_n",
    "stage_primary_covariates",
    "stage_APOE_covariates",
    "Braak_primary_model",
    "Braak_APOE_model",
    "CERAD_primary_model",
    "CERAD_APOE_model",
    "Braak_adjusted_for_CERAD",
    "multiple_testing",
    "git_HEAD"
  ),
  value = c(
    "Reviewer 3.6 APOE-e4 sensitivity",
    "22;23;24;33;34;44",
    "24,34,44",
    "400",
    "338",
    "NCI=148;MCI=86;AD=83",
    as.character(length(hsp_genes)),
    paste(protein_covars, collapse = ";"),
    paste(c(protein_covars, "apoe4_carrier_factor"), collapse = ";"),
    paste0(
      "safe_z(protein) ~ braak_num_std + ",
      paste(protein_covars, collapse = " + ")
    ),
    paste0(
      "safe_z(protein) ~ braak_num_std + ",
      paste(c(protein_covars, "apoe4_carrier_factor"), collapse = " + ")
    ),
    paste0(
      "safe_z(protein) ~ cerad_num_std + ",
      paste(protein_covars, collapse = " + ")
    ),
    paste0(
      "safe_z(protein) ~ cerad_num_std + ",
      paste(c(protein_covars, "apoe4_carrier_factor"), collapse = " + ")
    ),
    "FALSE",
    "Benjamini-Hochberg within endpoint/model across 306 Hsp60/10 clients",
    paste(git_head, collapse = ";")
  )
)
utils::write.csv(
  provenance,
  file.path(OUTDIR, "APOE_sensitivity_provenance.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "APOE_sensitivity_sessionInfo.txt")
)

############################################################
## Console report
############################################################

cat("\n============================================================\n")
cat("R3.6 APOE-e4 SENSITIVITY\n")
cat("============================================================\n\n")
cat("Protein cohort:\n")
cat("  Full: 400\n")
cat("  APOE complete: 338\n")
cat("  e4 non-carrier: 271\n")
cat("  e4 carrier: 67\n")
cat("  APOE-complete pure-stage phenotype eligible: NCI=149, MCI=86, AD=83\n")
cat("  APOE-complete pure-stage nuisance complete: NCI=148, MCI=86, AD=83\n\n")

cat("Production-model validation:\n")
print(production_validation, n = Inf, width = Inf)

cat("\nSame-subset effect comparison:\n")
print(effect_summary, n = Inf, width = Inf)

cat("\nVulnerability-score concordance:\n")
print(vulnerability_summary, n = Inf, width = Inf)

cat("\nTop-quartile retention:\n")
print(topq_overlap, n = Inf, width = Inf)

cat("\nFull original vs APOE-complete subset (WITHOUT APOE adjustment):\n")
print(full_vs_subset_summary, n = Inf, width = Inf)

cat("\nOutputs written to:\n", OUTDIR, "\n", sep = "")
cat("\n============================================================\n")
cat("APOE SENSITIVITY ANALYSIS COMPLETE — VALIDATION STILL REQUIRED\n")
cat("============================================================\n")
