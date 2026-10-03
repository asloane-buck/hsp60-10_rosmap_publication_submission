############################################################
## 47_run_alternative_mechanism_marker_stage_analysis.R
##
## Reviewer 3.1
## Prespecified alternative-mechanism TMT marker panel:
## NCI / MCI / AD differential abundance using the same
## protein-stage limma model as the validated conventional
## differential analysis.
############################################################

options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  if (!file.exists(f)) {
    stop("Missing required production script: ", f)
  }
  source(f)
}

require_objects(
  c(
    "prot_mat_raw",
    "prot_meta_adj",
    "protein_covars"
  ),
  context = "47_run_alternative_mechanism_marker_stage_analysis.R"
)

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required.")
}

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "alternative_mechanism_marker_panel"
)

manifest_file <- file.path(
  OUTDIR,
  "prespecified_marker_detectability_manifest.csv"
)

if (!file.exists(manifest_file)) {
  stop(
    "Missing prespecified marker manifest: ",
    manifest_file
  )
}

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE
)

required_manifest_cols <- c(
  "category",
  "gene",
  "detected_raw_TMT"
)

if (!all(required_manifest_cols %in% colnames(manifest))) {
  stop("Marker manifest is missing required columns.")
}

## =========================================================
## 1. Canonical protein matrix
## =========================================================

protein_mat <- as.matrix(prot_mat_raw)
mode(protein_mat) <- "numeric"

if (ncol(protein_mat) != 400L) {
  stop(
    "Expected 400 columns in production protein matrix; found ",
    ncol(protein_mat)
  )
}

protein_symbols <- canonical_gene_symbol(
  rownames(protein_mat)
)

if (anyDuplicated(protein_symbols) > 0) {
  stop(
    "Duplicate canonical protein symbols detected."
  )
}

rownames(protein_mat) <- protein_symbols

prot_meta_all <- prot_meta_adj[
  match(
    colnames(protein_mat),
    prot_meta_adj$SampleID
  ),
  ,
  drop = FALSE
]

if (
  !all(
    colnames(protein_mat) ==
      prot_meta_all$SampleID
  )
) {
  stop(
    "Protein matrix/metadata alignment failed."
  )
}


## =========================================================
## 2. Exact production covariates
## =========================================================

protein_covars_local <- unique(
  as.character(protein_covars)
)

expected_covars <- c(
  "age_num",
  "sex_factor",
  "pmi_num",
  "batch_factor"
)

if (
  !identical(
    protein_covars_local,
    expected_covars
  )
) {
  stop(
    "Protein covariate specification changed. Observed: ",
    paste(
      protein_covars_local,
      collapse = ", "
    )
  )
}

stage_levels <- c(
  "NCI",
  "MCI",
  "AD"
)

prot_meta_model <- prot_meta_all |>
  dplyr::filter(
    .data$clinical_stage %in%
      stage_levels
  ) |>
  dplyr::mutate(
    SampleID = as.character(
      .data$SampleID
    ),
    clinical_stage = factor(
      as.character(
        .data$clinical_stage
      ),
      levels = stage_levels
    )
  )

protein_complete_covars <- complete.cases(
  prot_meta_model[
    ,
    c(
      protein_covars_local,
      "clinical_stage"
    ),
    drop = FALSE
  ]
)

prot_meta_model <- prot_meta_model[
  protein_complete_covars,
  ,
  drop = FALSE
]

expected_counts <- c(
  NCI = 167L,
  MCI = 96L,
  AD = 109L
)

observed_counts <- table(
  prot_meta_model$clinical_stage
)

if (
  !identical(
    as.integer(
      observed_counts[
        stage_levels
      ]
    ),
    as.integer(
      expected_counts[
        stage_levels
      ]
    )
  )
) {
  stop(
    "Protein nuisance-complete stage counts changed."
  )
}

for (nm in protein_covars_local) {
  if (
    is.character(
      prot_meta_model[[nm]]
    )
  ) {
    prot_meta_model[[nm]] <- factor(
      prot_meta_model[[nm]]
    )
  }

  if (
    is.factor(
      prot_meta_model[[nm]]
    )
  ) {
    prot_meta_model[[nm]] <- droplevels(
      prot_meta_model[[nm]]
    )
  }
}


## =========================================================
## 3. Exact production design / contrasts
## =========================================================

protein_formula <- stats::as.formula(
  paste(
    "~ 0 +",
    paste(
      c(
        "clinical_stage",
        protein_covars_local
      ),
      collapse = " + "
    )
  )
)

protein_design <- stats::model.matrix(
  protein_formula,
  data = prot_meta_model
)

if (
  qr(protein_design)$rank !=
    ncol(protein_design)
) {
  stop(
    "Protein limma design matrix is not full rank."
  )
}

needed_stage_coef <- paste0(
  "clinical_stage",
  stage_levels
)

if (
  !all(
    needed_stage_coef %in%
      colnames(protein_design)
  )
) {
  stop(
    "Protein design missing stage coefficients."
  )
}

contrast_matrix <- limma::makeContrasts(
  MCI_vs_NCI =
    clinical_stageMCI -
    clinical_stageNCI,

  AD_vs_NCI =
    clinical_stageAD -
    clinical_stageNCI,

  AD_vs_MCI =
    clinical_stageAD -
    clinical_stageMCI,

  levels = protein_design
)


## =========================================================
## 4. Detected prespecified markers
## =========================================================

detected_manifest <- manifest |>
  dplyr::filter(
    .data$detected_raw_TMT
  )

if (
  any(
    !detected_manifest$gene %in%
      rownames(protein_mat)
  )
) {
  stop(
    "Manifest marks genes detected that are absent ",
    "from prot_mat_raw."
  )
}

marker_mat <- protein_mat[
  detected_manifest$gene,
  prot_meta_model$SampleID,
  drop = FALSE
]

if (
  !all(
    colnames(marker_mat) ==
      prot_meta_model$SampleID
  )
) {
  stop(
    "Marker matrix/sample alignment failed."
  )
}

cat(
  "Prespecified markers:",
  nrow(manifest),
  "\n"
)

cat(
  "Detected markers:",
  nrow(detected_manifest),
  "\n"
)

cat(
  "Protein model samples:",
  ncol(marker_mat),
  "\n"
)


## =========================================================
## 5. Marker-specific observation counts
## =========================================================

count_marker_stage <- function(gene) {
  y <- marker_mat[
    gene,
    ,
    drop = TRUE
  ]

  finite <- is.finite(y)

  data.frame(
    gene = gene,
    n_total = sum(finite),

    n_NCI = sum(
      finite &
        prot_meta_model$clinical_stage ==
          "NCI"
    ),

    n_MCI = sum(
      finite &
        prot_meta_model$clinical_stage ==
          "MCI"
    ),

    n_AD = sum(
      finite &
        prot_meta_model$clinical_stage ==
          "AD"
    ),

    stringsAsFactors = FALSE
  )
}

observation_counts <- dplyr::bind_rows(
  lapply(
    rownames(marker_mat),
    count_marker_stage
  )
)

observation_counts <- detected_manifest |>
  dplyr::select(
    category,
    gene
  ) |>
  dplyr::left_join(
    observation_counts,
    by = "gene"
  )


## =========================================================
## 6. Fit same limma model used by production analysis
## =========================================================

warning_messages <- character(0)

fit <- withCallingHandlers(
  limma::lmFit(
    marker_mat,
    protein_design
  ),
  warning = function(w) {
    warning_messages <<- c(
      warning_messages,
      conditionMessage(w)
    )
    invokeRestart(
      "muffleWarning"
    )
  }
)

coefficient_na <- is.na(
  fit$coefficients
)

estimability <- data.frame(
  gene = rownames(
    fit$coefficients
  ),

  any_na_coefficient =
    rowSums(
      coefficient_na
    ) > 0,

  any_na_stage_coefficient =
    rowSums(
      coefficient_na[
        ,
        needed_stage_coef,
        drop = FALSE
      ]
    ) > 0,

  stringsAsFactors = FALSE
)


## =========================================================
## 6b. Explicit partial-NA coefficient audit
## =========================================================
##
## Some markers are absent from entire TMT batches.
##
## With treatment-coded batch_factor, if the absent batch is
## the reference level, the resulting aliased coefficient need
## not carry the name of that missing reference batch. Therefore
## the valid structural check is gene-level:
##
##   * all NA coefficients must be batch_factor coefficients;
##   * no stage/age/sex/PMI coefficient may be NA;
##   * for each affected marker, the number of NA batch
##     coefficients must equal the number of TMT batches with
##     zero finite observations for that marker.

partial_na_details <- dplyr::bind_rows(
  lapply(
    rownames(fit$coefficients),
    function(g) {

      bad_coef <- colnames(fit$coefficients)[
        is.na(fit$coefficients[g, ])
      ]

      if (length(bad_coef) == 0) {
        return(NULL)
      }

      data.frame(
        gene = g,
        coefficient = bad_coef,
        is_batch_coefficient = grepl(
          "^batch_factor",
          bad_coef
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

if (nrow(partial_na_details) == 0) {
  partial_na_details <- data.frame(
    gene = character(),
    coefficient = character(),
    is_batch_coefficient = logical(),
    stringsAsFactors = FALSE
  )
}

if (
  any(
    estimability$any_na_stage_coefficient
  )
) {
  stop(
    "At least one marker has a non-estimable ",
    "clinical-stage coefficient."
  )
}

if (
  nrow(partial_na_details) > 0 &&
  any(
    !partial_na_details$is_batch_coefficient
  )
) {
  bad <- partial_na_details[
    !partial_na_details$is_batch_coefficient,
    ,
    drop = FALSE
  ]

  stop(
    "Unexpected non-batch NA coefficient(s): ",
    paste(
      paste0(
        bad$gene,
        ":",
        bad$coefficient
      ),
      collapse = "; "
    )
  )
}

batch_levels <- levels(
  factor(
    prot_meta_model$batch_factor
  )
)

partial_na_gene_audit <- dplyr::bind_rows(
  lapply(
    rownames(marker_mat),
    function(g) {

      y <- marker_mat[
        g,
        ,
        drop = TRUE
      ]

      zero_batches <- batch_levels[
        vapply(
          batch_levels,
          function(b) {
            idx <- as.character(
              prot_meta_model$batch_factor
            ) == b

            !any(
              is.finite(
                y[idx]
              )
            )
          },
          logical(1)
        )
      ]

      na_coefficients <- partial_na_details$coefficient[
        partial_na_details$gene == g
      ]

      data.frame(
        gene = g,

        n_zero_data_batches =
          length(zero_batches),

        zero_data_batches =
          paste(
            zero_batches,
            collapse = ";"
          ),

        n_na_batch_coefficients =
          length(na_coefficients),

        na_batch_coefficients =
          paste(
            na_coefficients,
            collapse = ";"
          ),

        counts_match =
          length(zero_batches) ==
            length(na_coefficients),

        stringsAsFactors = FALSE
      )
    }
  )
)

if (
  any(
    !partial_na_gene_audit$counts_match
  )
) {
  bad <- partial_na_gene_audit[
    !partial_na_gene_audit$counts_match,
    ,
    drop = FALSE
  ]

  stop(
    "Partial-NA coefficient count does not match ",
    "zero-data TMT-batch count for: ",
    paste(
      bad$gene,
      collapse = ", "
    )
  )
}

unique_model_warnings <- unique(
  warning_messages
)

unexpected_model_warnings <- unique_model_warnings[
  !grepl(
    "^Partial NA coefficients for [0-9]+ probe\\(s\\)$",
    unique_model_warnings
  )
]

if (
  length(
    unexpected_model_warnings
  ) > 0
) {
  stop(
    "Unexpected limma warning(s): ",
    paste(
      unexpected_model_warnings,
      collapse = " | "
    )
  )
}

fit2 <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- limma::eBayes(
  fit2
)


## =========================================================
## 7. Extract all three contrasts
## =========================================================

contrast_names <- colnames(
  contrast_matrix
)

results <- dplyr::bind_rows(
  lapply(
    contrast_names,
    function(contrast_i) {
      tt <- limma::topTable(
        fit2,
        coef = contrast_i,
        number = Inf,
        sort.by = "none",
        adjust.method = "BH",
        confint = 0.95
      ) |>
        tibble::rownames_to_column(
          "gene"
        )

      ## Recalculate BH explicitly across the detected
      ## prespecified marker family for this contrast.
      tt$marker_panel_fdr <- p.adjust(
        tt$P.Value,
        method = "BH"
      )

      tt |>
        dplyr::transmute(
          gene = canonical_gene_symbol(
            gene
          ),

          contrast = contrast_i,

          effect =
            as.numeric(
              .data$logFC
            ),

          effect_label =
            "limma log2 TMT abundance difference",

          std_error =
            abs(
              as.numeric(
                .data$logFC
              ) /
                as.numeric(
                  .data$t
                )
            ),

          conf_low =
            as.numeric(
              .data$CI.L
            ),

          conf_high =
            as.numeric(
              .data$CI.R
            ),

          statistic =
            as.numeric(
              .data$t
            ),

          p_value =
            as.numeric(
              .data$P.Value
            ),

          marker_panel_fdr =
            as.numeric(
              marker_panel_fdr
            )
        )
    }
  )
)

results <- detected_manifest |>
  dplyr::select(
    category,
    gene
  ) |>
  dplyr::left_join(
    results,
    by = "gene"
  ) |>
  dplyr::left_join(
    observation_counts,
    by = c(
      "category",
      "gene"
    )
  ) |>
  dplyr::left_join(
    estimability,
    by = "gene"
  ) |>
  dplyr::mutate(
    significant_marker_panel_fdr05 =
      is.finite(
        marker_panel_fdr
      ) &
        marker_panel_fdr <
          0.05
  )


## =========================================================
## 8. Per-category summary
## =========================================================

category_contrast_summary <- results |>
  dplyr::group_by(
    category,
    contrast
  ) |>
  dplyr::summarise(
    n_markers =
      dplyr::n(),

    n_estimable =
      sum(
        is.finite(
          effect
        )
      ),

    n_decreased =
      sum(
        is.finite(
          effect
        ) &
          effect < 0
      ),

    n_increased =
      sum(
        is.finite(
          effect
        ) &
          effect > 0
      ),

    median_effect =
      median(
        effect,
        na.rm = TRUE
      ),

    mean_effect =
      mean(
        effect,
        na.rm = TRUE
      ),

    n_nominal_p05 =
      sum(
        is.finite(
          p_value
        ) &
          p_value < 0.05
      ),

    n_fdr05 =
      sum(
        .data$significant_marker_panel_fdr05,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


## =========================================================
## 9. Unavailable requested markers
## =========================================================

unavailable <- manifest |>
  dplyr::filter(
    !.data$detected_raw_TMT
  ) |>
  dplyr::select(
    category,
    gene
  )


## =========================================================
## 10. Hard estimability audit
## =========================================================

estimability_summary <- estimability |>
  dplyr::summarise(
    n_detected_markers =
      dplyr::n(),

    n_any_na_coefficient =
      sum(
        .data$any_na_coefficient
      ),

    n_na_stage_coefficient =
      sum(
        .data$any_na_stage_coefficient
      )
  )

non_estimable_results <- results |>
  dplyr::filter(
    !is.finite(
      effect
    ) |
      !is.finite(
        p_value
      )
  )


## =========================================================
## 11. Write outputs
## =========================================================

write.csv(
  results,
  file.path(
    OUTDIR,
    "marker_stage_limma_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  observation_counts,
  file.path(
    OUTDIR,
    "marker_stage_observation_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  category_contrast_summary,
  file.path(
    OUTDIR,
    "marker_category_stage_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  unavailable,
  file.path(
    OUTDIR,
    "undetected_prespecified_markers.csv"
  ),
  row.names = FALSE
)

write.csv(
  estimability,
  file.path(
    OUTDIR,
    "marker_model_estimability.csv"
  ),
  row.names = FALSE
)

write.csv(
  estimability_summary,
  file.path(
    OUTDIR,
    "marker_model_estimability_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  non_estimable_results,
  file.path(
    OUTDIR,
    "non_estimable_marker_contrasts.csv"
  ),
  row.names = FALSE
)

warning_audit <- data.frame(
  n_warning_messages =
    length(
      unique(
        warning_messages
      )
    ),

  warning_messages =
    paste(
      unique(
        warning_messages
      ),
      collapse = " | "
    ),

  stringsAsFactors = FALSE
)

write.csv(
  warning_audit,
  file.path(
    OUTDIR,
    "marker_model_warning_audit.csv"
  ),
  row.names = FALSE
)


## =========================================================

write.csv(
  partial_na_details,
  file.path(
    OUTDIR,
    "marker_partial_na_coefficient_details.csv"
  ),
  row.names = FALSE
)


write.csv(
  partial_na_gene_audit,
  file.path(
    OUTDIR,
    "marker_partial_na_gene_audit.csv"
  ),
  row.names = FALSE
)

## 12. Provenance
## =========================================================

git_head <- tryCatch(
  system2(
    "git",
    c(
      "rev-parse",
      "HEAD"
    ),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) NA_character_
)

provenance <- data.frame(
  item = c(
    "protein_expression_object",
    "starting_protein_samples",
    "model_samples",
    "stage_counts",
    "protein_covariates",
    "design",
    "contrasts",
    "prespecified_markers",
    "detected_markers",
    "multiple_testing_family",
    "multiple_testing_method",
    "git_HEAD"
  ),

  value = c(
    "prot_mat_raw",
    ncol(protein_mat),
    nrow(prot_meta_model),

    paste(
      names(expected_counts),
      expected_counts,
      sep = "=",
      collapse = ";"
    ),

    paste(
      protein_covars_local,
      collapse = ";"
    ),

    deparse(
      protein_formula
    ),

    paste(
      contrast_names,
      collapse = ";"
    ),

    nrow(manifest),
    nrow(detected_manifest),

    "all detected prespecified R3.1 markers within each stage contrast",
    "Benjamini-Hochberg",

    paste(
      git_head,
      collapse = ";"
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(
    OUTDIR,
    "marker_stage_analysis_provenance.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "marker_stage_sessionInfo.txt"
  )
)


## =========================================================
## 13. Console report
## =========================================================

cat(
  "\n============================================================\n"
)

cat(
  "R3.1 ALTERNATIVE-MECHANISM MARKER STAGE ANALYSIS\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Protein cohort:",
  paste(
    names(expected_counts),
    expected_counts,
    sep = "=",
    collapse = ", "
  ),
  "\n"
)

cat(
  "Prespecified markers:",
  nrow(manifest),
  "\n"
)

cat(
  "Detected markers:",
  nrow(detected_manifest),
  "\n"
)

cat(
  "Undetected markers:",
  nrow(unavailable),
  "\n\n"
)

cat(
  "Model warning messages:",
  nrow(
    unique(
      data.frame(
        warning_messages
      )
    )
  ),
  "\n"
)

cat(
  "Markers with non-estimable stage coefficients:",
  estimability_summary$n_na_stage_coefficient,
  "\n"
)

cat(
  "Non-estimable marker/contrast results:",
  nrow(non_estimable_results),
  "\n\n"
)

cat(
  "Partial-NA coefficients explained by TMT-batch missingness:",
  (
    nrow(partial_na_details) > 0 &&
    all(partial_na_details$is_batch_coefficient) &&
    all(partial_na_gene_audit$counts_match) &&
    !any(estimability$any_na_stage_coefficient)
  ),
  "\n\n"
)

cat(
  "Category/contrast summary:\n"
)

print(
  category_contrast_summary,
  n = Inf
)

cat(
  "\nFDR-significant marker results:\n"
)

sig_results <- results |>
  dplyr::filter(
    significant_marker_panel_fdr05
  ) |>
  dplyr::select(
    category,
    gene,
    contrast,
    effect,
    conf_low,
    conf_high,
    p_value,
    marker_panel_fdr,
    n_total,
    n_NCI,
    n_MCI,
    n_AD
  )

print(
  tibble::as_tibble(sig_results),
  n = Inf,
  width = Inf
)

cat(
  "\nUndetected prespecified markers:\n"
)

print(
  unavailable,
  row.names = FALSE
)

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)
