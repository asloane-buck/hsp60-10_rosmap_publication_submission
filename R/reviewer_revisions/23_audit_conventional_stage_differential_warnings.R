#!/usr/bin/env Rscript

############################################################
## 23_audit_conventional_stage_differential_warnings.R
##
## READ-ONLY audit of the first conventional stage differential
## analysis run. This script does not overwrite analytical results.
##
## It specifically addresses:
##   1. DESeq2 numeric-covariate scale message.
##   2. limma "Partial NA coefficients" warning.
##   3. Ensembl -> gene-symbol one-to-many mapping message.
##   4. Whether saved stage-contrast results remain fully estimable.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
  library(readr)
})

source("R/00_config.R")
source("R/01_utils.R")

out_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential"
)

audit_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "conventional_stage_differential_warning_audit"
)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

write_audit <- function(x, filename) {
  path <- file.path(audit_dir, filename)
  readr::write_csv(x, path, na = "")
  message("Wrote: ", path)
  invisible(x)
}

load_required_object <- function(name) {
  path <- file.path(cfg$object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing production object: ", path, call. = FALSE)
  }
  readRDS(path)
}

required_result_files <- c(
  "02_RNA_DESeq2_all_contrasts.csv",
  "03_protein_limma_all_contrasts.csv",
  "04_reviewer_contrast_summary.csv",
  "11_sample_stage_audit.csv",
  "12_model_specification_audit.csv",
  "15_validation_checks.csv"
)

missing_results <- required_result_files[
  !file.exists(file.path(out_dir, required_result_files))
]
if (length(missing_results) > 0) {
  stop(
    "Missing first-run result file(s): ",
    paste(missing_results, collapse = ", "),
    call. = FALSE
  )
}

rna_results <- readr::read_csv(
  file.path(out_dir, "02_RNA_DESeq2_all_contrasts.csv"),
  show_col_types = FALSE
)
protein_results_saved <- readr::read_csv(
  file.path(out_dir, "03_protein_limma_all_contrasts.csv"),
  show_col_types = FALSE
)
sample_audit_saved <- readr::read_csv(
  file.path(out_dir, "11_sample_stage_audit.csv"),
  show_col_types = FALSE
)
internal_validation <- readr::read_csv(
  file.path(out_dir, "15_validation_checks.csv"),
  show_col_types = FALSE
)

rna_meta_adj <- load_required_object("rna_meta_adj")
rna_covars <- as.character(load_required_object("rna_covars"))
prot_mat_raw <- load_required_object("prot_mat_raw")
prot_meta_adj <- load_required_object("prot_meta_adj")
protein_covars <- as.character(load_required_object("protein_covars"))
all_hsp60_10_client_tbl <- load_required_object("all_hsp60_10_client_tbl")

stage_levels <- c("NCI", "MCI", "AD")
contrast_levels <- c("MCI_vs_NCI", "AD_vs_NCI", "AD_vs_MCI")

# ============================================================
# 1. Confirm first-run validation and canonical sample counts
# ============================================================

if (!all(internal_validation$passed)) {
  stop(
    "The first-run internal validation table contains failed checks. ",
    "Do not interpret the differential results.",
    call. = FALSE
  )
}

expected_sample_counts <- tibble::tribble(
  ~modality, ~clinical_stage, ~expected_n,
  "RNA", "NCI", 200L,
  "RNA", "MCI", 158L,
  "RNA", "AD", 220L,
  "Protein", "NCI", 167L,
  "Protein", "MCI", 96L,
  "Protein", "AD", 109L
)

sample_check <- expected_sample_counts |>
  dplyr::left_join(
    sample_audit_saved |>
      dplyr::select("modality", "clinical_stage", observed_n = "n"),
    by = c("modality", "clinical_stage")
  ) |>
  dplyr::mutate(
    passed = .data$observed_n == .data$expected_n
  )

write_audit(sample_check, "01_sample_count_check.csv")

if (any(!sample_check$passed)) {
  stop("Canonical sample-count audit failed.", call. = FALSE)
}

# ============================================================
# 2. Audit DESeq2 covariate scale message
# ============================================================

missing_rna_covars <- setdiff(rna_covars, colnames(rna_meta_adj))
if (length(missing_rna_covars) > 0) {
  stop(
    "RNA covariates missing from production metadata: ",
    paste(missing_rna_covars, collapse = ", "),
    call. = FALSE
  )
}

rna_covariate_scale_audit <- dplyr::bind_rows(
  lapply(rna_covars, function(nm) {
    x <- rna_meta_adj[[nm]]

    if (is.numeric(x) || is.integer(x)) {
      tibble::tibble(
        covariate = nm,
        storage_class = class(x)[1],
        is_numeric = TRUE,
        n_non_missing = sum(is.finite(x)),
        n_unique = dplyr::n_distinct(x[is.finite(x)]),
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE),
        triggers_deseq_scale_message =
          abs(mean(x, na.rm = TRUE)) > 5 ||
          stats::sd(x, na.rm = TRUE) > 5
      )
    } else {
      tibble::tibble(
        covariate = nm,
        storage_class = class(x)[1],
        is_numeric = FALSE,
        n_non_missing = sum(!is.na(x)),
        n_unique = dplyr::n_distinct(x[!is.na(x)]),
        mean = NA_real_,
        sd = NA_real_,
        min = NA_real_,
        max = NA_real_,
        triggers_deseq_scale_message = FALSE
      )
    }
  })
)

write_audit(
  rna_covariate_scale_audit,
  "02_RNA_covariate_scale_audit.csv"
)

# ============================================================
# 3. Refit exact protein model and audit partial NA coefficients
# ============================================================

require_columns(
  prot_meta_adj,
  c("SampleID", "clinical_stage"),
  "protein production metadata"
)

missing_prot_covars <- setdiff(
  protein_covars,
  colnames(prot_meta_adj)
)
if (length(missing_prot_covars) > 0) {
  stop(
    "Protein covariates missing from production metadata: ",
    paste(missing_prot_covars, collapse = ", "),
    call. = FALSE
  )
}

protein_mat <- as.matrix(prot_mat_raw)
mode(protein_mat) <- "numeric"
rownames(protein_mat) <- canonical_gene_symbol(rownames(protein_mat))

if (anyDuplicated(rownames(protein_mat)) > 0) {
  stop(
    "Duplicate canonical protein gene symbols detected in prot_mat_raw.",
    call. = FALSE
  )
}

prot_meta_model <- prot_meta_adj |>
  dplyr::filter(.data$clinical_stage %in% stage_levels) |>
  dplyr::mutate(
    SampleID = as.character(.data$SampleID),
    clinical_stage = factor(
      as.character(.data$clinical_stage),
      levels = stage_levels
    )
  )

protein_complete <- complete.cases(
  prot_meta_model[
    ,
    c(protein_covars, "clinical_stage"),
    drop = FALSE
  ]
)

prot_meta_model <- prot_meta_model[
  protein_complete,
  ,
  drop = FALSE
]

if (!all(prot_meta_model$SampleID %in% colnames(protein_mat))) {
  stop(
    "Protein model metadata contains samples absent from prot_mat_raw.",
    call. = FALSE
  )
}

protein_mat_model <- protein_mat[
  ,
  prot_meta_model$SampleID,
  drop = FALSE
]

if (!all(colnames(protein_mat_model) == prot_meta_model$SampleID)) {
  stop("Protein matrix/metadata alignment failed.", call. = FALSE)
}

for (nm in protein_covars) {
  if (is.character(prot_meta_model[[nm]])) {
    prot_meta_model[[nm]] <- factor(prot_meta_model[[nm]])
  }
  if (is.factor(prot_meta_model[[nm]])) {
    prot_meta_model[[nm]] <- droplevels(prot_meta_model[[nm]])
  }
}

protein_formula <- stats::as.formula(
  paste(
    "~ 0 +",
    paste(
      c("clinical_stage", protein_covars),
      collapse = " + "
    )
  )
)

protein_design <- stats::model.matrix(
  protein_formula,
  data = prot_meta_model
)

if (qr(protein_design)$rank != ncol(protein_design)) {
  stop("Protein design matrix is not full rank.", call. = FALSE)
}

needed_stage_coef <- paste0(
  "clinical_stage",
  stage_levels
)

if (!all(needed_stage_coef %in% colnames(protein_design))) {
  stop(
    "Protein design missing expected stage coefficients.",
    call. = FALSE
  )
}

protein_fit <- suppressWarnings(
  limma::lmFit(
    protein_mat_model,
    protein_design
  )
)

coef_na_matrix <- is.na(protein_fit$coefficients)

protein_coefficient_audit <- tibble::tibble(
  coefficient = colnames(protein_fit$coefficients),
  coefficient_role = dplyr::if_else(
    colnames(protein_fit$coefficients) %in% needed_stage_coef,
    "stage",
    "nuisance"
  ),
  n_proteins_na = colSums(coef_na_matrix),
  n_proteins_finite = colSums(
    is.finite(protein_fit$coefficients)
  )
)

write_audit(
  protein_coefficient_audit,
  "03_protein_coefficient_estimability_by_term.csv"
)

protein_probe_estimability <- tibble::tibble(
  gene_symbol = rownames(protein_fit$coefficients),
  n_observed_samples = rowSums(
    is.finite(protein_mat_model)
  ),
  any_na_coefficient = rowSums(coef_na_matrix) > 0,
  any_na_stage_coefficient =
    rowSums(
      coef_na_matrix[
        ,
        needed_stage_coef,
        drop = FALSE
      ]
    ) > 0,
  any_na_nuisance_coefficient =
    rowSums(
      coef_na_matrix[
        ,
        setdiff(
          colnames(protein_fit$coefficients),
          needed_stage_coef
        ),
        drop = FALSE
      ]
    ) > 0
)

write_audit(
  protein_probe_estimability,
  "04_protein_probe_estimability.csv"
)

protein_estimability_summary <- tibble::tibble(
  metric = c(
    "proteins modeled",
    "proteins with any NA coefficient",
    "proteins with any NA stage coefficient",
    "proteins with NA nuisance coefficient but all stage coefficients finite",
    "minimum observed samples per protein",
    "median observed samples per protein",
    "maximum observed samples per protein"
  ),
  value = c(
    nrow(protein_probe_estimability),
    sum(protein_probe_estimability$any_na_coefficient),
    sum(protein_probe_estimability$any_na_stage_coefficient),
    sum(
      protein_probe_estimability$any_na_nuisance_coefficient &
        !protein_probe_estimability$any_na_stage_coefficient
    ),
    min(protein_probe_estimability$n_observed_samples),
    stats::median(protein_probe_estimability$n_observed_samples),
    max(protein_probe_estimability$n_observed_samples)
  )
)

write_audit(
  protein_estimability_summary,
  "05_protein_estimability_summary.csv"
)

if (any(
  protein_probe_estimability$any_na_stage_coefficient
)) {
  stop(
    "At least one protein has a non-estimable stage coefficient. ",
    "Do not interpret protein stage contrasts yet.",
    call. = FALSE
  )
}

protein_contrast_matrix <- limma::makeContrasts(
  MCI_vs_NCI =
    clinical_stageMCI - clinical_stageNCI,
  AD_vs_NCI =
    clinical_stageAD - clinical_stageNCI,
  AD_vs_MCI =
    clinical_stageAD - clinical_stageMCI,
  levels = protein_design
)

protein_fit2 <- limma::contrasts.fit(
  protein_fit,
  protein_contrast_matrix
)
protein_fit2 <- limma::eBayes(protein_fit2)

protein_refit_results <- lapply(
  contrast_levels,
  function(contrast_i) {
    tt <- limma::topTable(
      protein_fit2,
      coef = contrast_i,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH",
      confint = 0.95
    ) |>
      tibble::rownames_to_column("gene_symbol")

    tibble::tibble(
      contrast = contrast_i,
      gene_symbol =
        canonical_gene_symbol(tt$gene_symbol),
      effect_refit = as.numeric(tt$logFC),
      t_refit = as.numeric(tt$t),
      p_value_refit = as.numeric(tt$P.Value),
      fdr_refit = as.numeric(tt$adj.P.Val)
    )
  }
) |>
  dplyr::bind_rows()

protein_contrast_estimability <- protein_refit_results |>
  dplyr::group_by(.data$contrast) |>
  dplyr::summarise(
    n_features = dplyr::n(),
    n_effect_finite = sum(is.finite(.data$effect_refit)),
    n_t_finite = sum(is.finite(.data$t_refit)),
    n_p_finite = sum(is.finite(.data$p_value_refit)),
    n_fdr_finite = sum(is.finite(.data$fdr_refit)),
    .groups = "drop"
  )

write_audit(
  protein_contrast_estimability,
  "06_protein_stage_contrast_estimability.csv"
)

if (any(
  protein_contrast_estimability$n_effect_finite !=
    protein_contrast_estimability$n_features
) ||
  any(
    protein_contrast_estimability$n_p_finite !=
      protein_contrast_estimability$n_features
  ) ||
  any(
    protein_contrast_estimability$n_fdr_finite !=
      protein_contrast_estimability$n_features
  )) {
  stop(
    "One or more protein stage contrasts are not fully estimable.",
    call. = FALSE
  )
}

protein_saved_compare <- protein_results_saved |>
  dplyr::select(
    "contrast",
    "gene_symbol",
    effect_saved = "effect",
    p_value_saved = "p_value",
    fdr_saved = "fdr"
  ) |>
  dplyr::inner_join(
    protein_refit_results,
    by = c("contrast", "gene_symbol")
  ) |>
  dplyr::mutate(
    abs_effect_difference =
      abs(.data$effect_saved - .data$effect_refit),
    abs_p_difference =
      abs(.data$p_value_saved - .data$p_value_refit),
    abs_fdr_difference =
      abs(.data$fdr_saved - .data$fdr_refit)
  )

protein_refit_comparison_summary <- tibble::tibble(
  metric = c(
    "saved/refit joined rows",
    "maximum absolute effect difference",
    "maximum absolute p-value difference",
    "maximum absolute FDR difference"
  ),
  value = c(
    nrow(protein_saved_compare),
    max(
      protein_saved_compare$abs_effect_difference,
      na.rm = TRUE
    ),
    max(
      protein_saved_compare$abs_p_difference,
      na.rm = TRUE
    ),
    max(
      protein_saved_compare$abs_fdr_difference,
      na.rm = TRUE
    )
  )
)

write_audit(
  protein_refit_comparison_summary,
  "07_protein_saved_vs_refit_summary.csv"
)

if (
  nrow(protein_saved_compare) !=
    nrow(protein_results_saved) ||
  max(
    protein_saved_compare$abs_effect_difference,
    na.rm = TRUE
  ) > 1e-10 ||
  max(
    protein_saved_compare$abs_p_difference,
    na.rm = TRUE
  ) > 1e-10 ||
  max(
    protein_saved_compare$abs_fdr_difference,
    na.rm = TRUE
  ) > 1e-10
) {
  stop(
    "Protein warning-audit refit does not exactly reproduce saved results.",
    call. = FALSE
  )
}

# ============================================================
# 4. Audit Ensembl -> SYMBOL ambiguity explicitly
# ============================================================

client_gene_col <- intersect(
  c("gene", "gene_symbol", "symbol"),
  colnames(all_hsp60_10_client_tbl)
)[1]

if (is.na(client_gene_col)) {
  stop(
    "Could not identify Hsp60/10 client gene column.",
    call. = FALSE
  )
}

hsp_clients <- all_hsp60_10_client_tbl[[client_gene_col]] |>
  as.character() |>
  canonical_gene_symbol() |>
  unique()

rna_feature_map_base <- rna_results |>
  dplyr::distinct(
    .data$feature_id,
    .data$ensembl_id
  ) |>
  dplyr::filter(
    !is.na(.data$ensembl_id),
    nzchar(.data$ensembl_id)
  )

annotation_rows <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(rna_feature_map_base$ensembl_id),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
)

annotation_summary <- annotation_rows |>
  dplyr::mutate(
    SYMBOL = as.character(.data$SYMBOL)
  ) |>
  dplyr::group_by(.data$ENSEMBL) |>
  dplyr::summarise(
    n_symbol_mappings = dplyr::n_distinct(
      .data$SYMBOL[
        !is.na(.data$SYMBOL) &
          nzchar(.data$SYMBOL)
      ]
    ),
    symbol_candidates = paste(
      sort(
        unique(
          .data$SYMBOL[
            !is.na(.data$SYMBOL) &
              nzchar(.data$SYMBOL)
          ]
        )
      ),
      collapse = ";"
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    uniquely_mapped_symbol = dplyr::if_else(
      .data$n_symbol_mappings == 1L,
      canonical_gene_symbol(.data$symbol_candidates),
      NA_character_
    )
  )

mapping_audit <- rna_feature_map_base |>
  dplyr::left_join(
    annotation_summary,
    by = c("ensembl_id" = "ENSEMBL")
  ) |>
  dplyr::mutate(
    n_symbol_mappings = dplyr::coalesce(
      .data$n_symbol_mappings,
      0L
    ),
    mapping_status = dplyr::case_when(
      .data$n_symbol_mappings == 0L ~ "unmapped",
      .data$n_symbol_mappings == 1L ~ "unique",
      .data$n_symbol_mappings > 1L ~ "ambiguous",
      TRUE ~ "unknown"
    ),
    any_candidate_is_hsp60_10_client =
      vapply(
        strsplit(
          dplyr::coalesce(
            .data$symbol_candidates,
            ""
          ),
          ";",
          fixed = TRUE
        ),
        function(z) {
          z <- canonical_gene_symbol(z)
          any(z %in% hsp_clients)
        },
        logical(1)
      )
  )

write_audit(
  mapping_audit,
  "08_RNA_ensembl_symbol_mapping_audit.csv"
)

mapping_summary <- mapping_audit |>
  dplyr::count(
    .data$mapping_status,
    name = "n_features"
  )

write_audit(
  mapping_summary,
  "09_RNA_mapping_summary.csv"
)

ambiguous_hsp <- mapping_audit |>
  dplyr::filter(
    .data$mapping_status == "ambiguous",
    .data$any_candidate_is_hsp60_10_client
  )

write_audit(
  ambiguous_hsp,
  "10_ambiguous_mapping_Hsp60_10_audit.csv"
)

unique_symbol_feature_counts <- mapping_audit |>
  dplyr::filter(
    .data$mapping_status == "unique",
    !is.na(.data$uniquely_mapped_symbol)
  ) |>
  dplyr::count(
    .data$uniquely_mapped_symbol,
    name = "n_ensembl_features"
  )

duplicated_hsp_symbols <- unique_symbol_feature_counts |>
  dplyr::filter(
    .data$uniquely_mapped_symbol %in% hsp_clients,
    .data$n_ensembl_features > 1
  )

write_audit(
  duplicated_hsp_symbols,
  "11_Hsp60_10_symbols_with_multiple_ensembl_features.csv"
)

core_mapping <- mapping_audit |>
  dplyr::filter(
    grepl(
      "(^|;)(HSPD1|HSPE1)(;|$)",
      dplyr::coalesce(
        .data$symbol_candidates,
        ""
      )
    )
  )

write_audit(
  core_mapping,
  "12_HSPD1_HSPE1_mapping_audit.csv"
)

# ============================================================
# 5. Compact decision table
# ============================================================

stage_na_count <- sum(
  protein_probe_estimability$any_na_stage_coefficient
)

protein_all_contrasts_finite <- all(
  protein_contrast_estimability$n_effect_finite ==
    protein_contrast_estimability$n_features &
  protein_contrast_estimability$n_p_finite ==
    protein_contrast_estimability$n_features &
  protein_contrast_estimability$n_fdr_finite ==
    protein_contrast_estimability$n_features
)

decision <- tibble::tribble(
  ~check, ~observed, ~required, ~passed,
  "First-run internal validation",
  as.character(all(internal_validation$passed)),
  "TRUE",
  all(internal_validation$passed),

  "Canonical sample counts",
  as.character(all(sample_check$passed)),
  "TRUE",
  all(sample_check$passed),

  "Protein stage coefficients with NA",
  as.character(stage_na_count),
  "0",
  stage_na_count == 0,

  "Protein stage contrasts fully finite",
  as.character(protein_all_contrasts_finite),
  "TRUE",
  protein_all_contrasts_finite,

  "Protein warning-audit refit reproduces saved results",
  as.character(
    nrow(protein_saved_compare) ==
      nrow(protein_results_saved) &&
      max(
        protein_saved_compare$abs_effect_difference,
        na.rm = TRUE
      ) <= 1e-10 &&
      max(
        protein_saved_compare$abs_p_difference,
        na.rm = TRUE
      ) <= 1e-10 &&
      max(
        protein_saved_compare$abs_fdr_difference,
        na.rm = TRUE
      ) <= 1e-10
  ),
  "TRUE",
  nrow(protein_saved_compare) ==
    nrow(protein_results_saved) &&
    max(
      protein_saved_compare$abs_effect_difference,
      na.rm = TRUE
    ) <= 1e-10 &&
    max(
      protein_saved_compare$abs_p_difference,
      na.rm = TRUE
    ) <= 1e-10 &&
    max(
      protein_saved_compare$abs_fdr_difference,
      na.rm = TRUE
    ) <= 1e-10,

  "Ambiguous Ensembl mappings involving an Hsp60/10 client",
  as.character(nrow(ambiguous_hsp)),
  "0",
  nrow(ambiguous_hsp) == 0,

  "Hsp60/10 symbols represented by >1 uniquely mapped Ensembl feature",
  as.character(nrow(duplicated_hsp_symbols)),
  "0",
  nrow(duplicated_hsp_symbols) == 0
)

write_audit(
  decision,
  "00_warning_audit_decision_table.csv"
)

cat("\n============================================================\n")
cat("CONVENTIONAL DIFFERENTIAL ANALYSIS WARNING AUDIT\n")
cat("============================================================\n\n")

cat("RNA production covariates and scale:\n")
print(rna_covariate_scale_audit, n = Inf)

cat("\nProtein coefficient estimability:\n")
print(protein_estimability_summary, n = Inf)

cat("\nProtein stage-contrast estimability:\n")
print(protein_contrast_estimability, n = Inf)

cat("\nRNA Ensembl -> SYMBOL mapping summary:\n")
print(mapping_summary, n = Inf)

cat("\nDecision table:\n")
print(decision, n = Inf)

if (any(!decision$passed)) {
  stop(
    "WARNING AUDIT FOUND A BLOCKER. ",
    "Do not treat the conventional differential results as final.",
    call. = FALSE
  )
}

cat(
  "\nWARNING AUDIT PASSED: the limma partial-NA warning does not ",
  "affect stage-coefficient estimability, saved protein contrasts ",
  "are exactly reproduced, and Hsp60/10 symbol mapping is unambiguous.\n"
)
