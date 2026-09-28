#!/usr/bin/env Rscript

############################################################
## 24_run_conventional_stage_differential_FINAL.R
##
## Reviewer revision: conventional differential gene and
## protein abundance analyses across canonical NCI / MCI / AD.
##
## RNA:
##   DESeq2 on raw bulk RNA-seq counts.
##
## Protein PRIMARY:
##   limma on processed, non-residualized log2 TMT abundance,
##   restricted to proteins observed in ALL 400 TMT participants.
##
## Protein SENSITIVITY 1:
##   proteins absent from no more than 2 whole TMT batches.
##
## Protein SENSITIVITY 2:
##   proteins with <=20% participant-level missingness.
##
## Contrasts:
##   MCI vs NCI
##   AD  vs NCI
##   AD  vs MCI
##
## No RNA-protein effect subtraction is performed.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
  library(readr)
})

source("R/00_config.R")
source("R/01_utils.R")

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[1] == length(args)) stop("Missing value after ", flag, call. = FALSE)
  args[hit[1] + 1]
}

counts_path <- get_arg("--counts", Sys.getenv("ROSMAP_RAW_COUNTS_FILE", unset = ""))
if (!nzchar(counts_path)) {
  stop("Supply raw RNA-seq counts with --counts /path/to/count_matrix.rds", call. = FALSE)
}
counts_path <- path.expand(counts_path)
if (!file.exists(counts_path)) stop("Raw RNA-seq count file not found: ", counts_path, call. = FALSE)

out_dir <- get_arg(
  "--outdir",
  file.path(cfg$project_dir, "outputs", "reviewer_revisions", "conventional_stage_differential_final")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_out <- function(x, filename) {
  path <- file.path(out_dir, filename)
  readr::write_csv(x, path, na = "")
  message("Wrote: ", path)
  invisible(x)
}

load_required_object <- function(name) {
  path <- file.path(cfg$object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) stop("Missing production object: ", path, call. = FALSE)
  readRDS(path)
}

stage_levels <- c("NCI", "MCI", "AD")

contrast_specs <- tibble::tribble(
  ~contrast, ~numerator, ~denominator,
  "MCI_vs_NCI", "MCI", "NCI",
  "AD_vs_NCI",  "AD",  "NCI",
  "AD_vs_MCI",  "AD",  "MCI"
)

rna_meta_adj <- load_required_object("rna_meta_adj")
rna_covars <- unique(as.character(load_required_object("rna_covars")))
prot_mat_raw <- load_required_object("prot_mat_raw")
prot_meta_adj <- load_required_object("prot_meta_adj")
protein_covars <- unique(as.character(load_required_object("protein_covars")))
all_hsp60_10_client_tbl <- load_required_object("all_hsp60_10_client_tbl")

require_columns(rna_meta_adj, c("sample_id", "clinical_stage"), "RNA production metadata")
require_columns(prot_meta_adj, c("SampleID", "clinical_stage", "batch_factor"), "protein production metadata")

missing_rna_covars <- setdiff(rna_covars, colnames(rna_meta_adj))
missing_prot_covars <- setdiff(protein_covars, colnames(prot_meta_adj))
if (length(missing_rna_covars) > 0) stop("Missing RNA production covariates: ", paste(missing_rna_covars, collapse = ", "), call. = FALSE)
if (length(missing_prot_covars) > 0) stop("Missing protein production covariates: ", paste(missing_prot_covars, collapse = ", "), call. = FALSE)

client_gene_col <- intersect(c("gene", "gene_symbol", "symbol"), colnames(all_hsp60_10_client_tbl))[1]
if (is.na(client_gene_col)) stop("Could not identify Hsp60/10 client gene column.", call. = FALSE)

hsp_clients <- all_hsp60_10_client_tbl[[client_gene_col]] |>
  as.character() |>
  canonical_gene_symbol() |>
  unique()
hsp_clients <- hsp_clients[!is.na(hsp_clients) & nzchar(hsp_clients)]
if (length(hsp_clients) != 321L) stop("Expected 321 canonical Hsp60/10 clients; observed ", length(hsp_clients), ".", call. = FALSE)

message("\n============================================================")
message("RNA DIFFERENTIAL EXPRESSION: DESeq2")
message("============================================================")

rna_counts <- readRDS(counts_path)
if (inherits(rna_counts, "data.frame")) rna_counts <- as.matrix(rna_counts)
if (is.null(rownames(rna_counts)) || is.null(colnames(rna_counts))) stop("RNA count matrix must have row/column names.", call. = FALSE)

normalized_sample_ids <- basename(colnames(rna_counts))
normalized_sample_ids <- sub("\\.bam$", "", normalized_sample_ids)
if (anyDuplicated(normalized_sample_ids) > 0) stop("Raw-count sample-ID normalization created duplicates.", call. = FALSE)
colnames(rna_counts) <- normalized_sample_ids

rna_meta_model <- rna_meta_adj |>
  dplyr::filter(.data$clinical_stage %in% stage_levels) |>
  dplyr::mutate(
    sample_id = as.character(.data$sample_id),
    clinical_stage = factor(as.character(.data$clinical_stage), levels = stage_levels)
  )

if (anyDuplicated(rna_meta_model$sample_id) > 0) stop("Duplicate RNA sample_id.", call. = FALSE)

rna_overlap <- sum(rna_meta_model$sample_id %in% colnames(rna_counts))
message("RNA production sample IDs matched after harmonization: ", rna_overlap, "/", nrow(rna_meta_model))
if (rna_overlap != nrow(rna_meta_model)) stop("Not all production RNA samples occur in raw counts.", call. = FALSE)

rna_counts <- rna_counts[, rna_meta_model$sample_id, drop = FALSE]
if (!all(colnames(rna_counts) == rna_meta_model$sample_id)) stop("RNA count/metadata alignment failed.", call. = FALSE)

expected_rna_counts <- c(NCI = 200L, MCI = 158L, AD = 219L)
observed_rna_counts <- table(rna_meta_model$clinical_stage)
if (!identical(as.integer(observed_rna_counts[stage_levels]), as.integer(expected_rna_counts[stage_levels]))) {
  stop("RNA canonical stage counts changed.", call. = FALSE)
}

if (anyNA(rna_counts) || any(rna_counts < 0)) stop("RNA counts contain NA or negative values.", call. = FALSE)
max_noninteger <- max(abs(rna_counts - round(rna_counts)))
if (!is.finite(max_noninteger) || max_noninteger > 1e-6) stop("RNA count matrix is not integer-like.", call. = FALSE)

rna_model_covars <- character(0)
transform_rows <- list()

for (nm in rna_covars) {
  x <- rna_meta_model[[nm]]
  if (is.numeric(x) || is.integer(x)) {
    new_nm <- paste0(nm, "_z")
    if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) stop("RNA numeric covariate cannot be standardized: ", nm, call. = FALSE)
    rna_meta_model[[new_nm]] <- as.numeric(scale(x))
    rna_model_covars <- c(rna_model_covars, new_nm)
    transform_rows[[length(transform_rows) + 1]] <- tibble::tibble(
      original_covariate = nm,
      model_covariate = new_nm,
      transformation = "z-score",
      original_mean = mean(x, na.rm = TRUE),
      original_sd = stats::sd(x, na.rm = TRUE),
      transformed_mean = mean(rna_meta_model[[new_nm]], na.rm = TRUE),
      transformed_sd = stats::sd(rna_meta_model[[new_nm]], na.rm = TRUE)
    )
  } else {
    if (is.character(x)) rna_meta_model[[nm]] <- factor(x)
    if (is.factor(rna_meta_model[[nm]])) rna_meta_model[[nm]] <- droplevels(rna_meta_model[[nm]])
    rna_model_covars <- c(rna_model_covars, nm)
    transform_rows[[length(transform_rows) + 1]] <- tibble::tibble(
      original_covariate = nm,
      model_covariate = nm,
      transformation = "none; categorical",
      original_mean = NA_real_,
      original_sd = NA_real_,
      transformed_mean = NA_real_,
      transformed_sd = NA_real_
    )
  }
}

rna_covariate_transform_audit <- dplyr::bind_rows(transform_rows)
write_out(rna_covariate_transform_audit, "01_RNA_covariate_transform_audit.csv")

if (!all(complete.cases(rna_meta_model[, c(rna_model_covars, "clinical_stage"), drop = FALSE]))) {
  stop("RNA model has incomplete covariates.", call. = FALSE)
}

rna_formula <- stats::as.formula(paste("~", paste(c(rna_model_covars, "clinical_stage"), collapse = " + ")))
rna_mm <- stats::model.matrix(rna_formula, data = rna_meta_model)
if (qr(rna_mm)$rank != ncol(rna_mm)) stop("RNA DESeq2 design matrix is not full rank.", call. = FALSE)

rna_counts_integer <- round(rna_counts)
storage.mode(rna_counts_integer) <- "integer"

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = rna_counts_integer,
  colData = rna_meta_model,
  design = rna_formula
)

rna_genes_before_filter <- nrow(dds)
rna_keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
dds <- dds[rna_keep, ]
rna_genes_after_filter <- nrow(dds)

message("RNA genes before filter: ", rna_genes_before_filter, "; after count>=10 in >=3 samples: ", rna_genes_after_filter)
dds <- DESeq2::DESeq(dds, quiet = TRUE)

ensembl_id_version <- rownames(dds)
ensembl_id <- sub("\\..*$", "", ensembl_id_version)

annotation_long <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(ensembl_id),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
)

annotation_summary <- annotation_long |>
  dplyr::mutate(SYMBOL = as.character(.data$SYMBOL)) |>
  dplyr::group_by(.data$ENSEMBL) |>
  dplyr::summarise(
    n_symbol_mappings = dplyr::n_distinct(.data$SYMBOL[!is.na(.data$SYMBOL) & nzchar(.data$SYMBOL)]),
    symbol_candidates = paste(sort(unique(.data$SYMBOL[!is.na(.data$SYMBOL) & nzchar(.data$SYMBOL)])), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    gene_symbol = dplyr::if_else(.data$n_symbol_mappings == 1L, canonical_gene_symbol(.data$symbol_candidates), NA_character_),
    mapping_status = dplyr::case_when(
      .data$n_symbol_mappings == 0L ~ "unmapped",
      .data$n_symbol_mappings == 1L ~ "unique",
      .data$n_symbol_mappings > 1L ~ "ambiguous",
      TRUE ~ "unknown"
    )
  )

rna_gene_map <- tibble::tibble(
  ensembl_id_version = ensembl_id_version,
  ensembl_id = ensembl_id
) |>
  dplyr::left_join(annotation_summary, by = c("ensembl_id" = "ENSEMBL"))

write_out(rna_gene_map, "02_RNA_ensembl_symbol_mapping.csv")
write_out(rna_gene_map |> dplyr::count(.data$mapping_status, name = "n_features"), "03_RNA_mapping_summary.csv")

candidate_lists <- strsplit(dplyr::coalesce(rna_gene_map$symbol_candidates, ""), ";", fixed = TRUE)
ambiguous_hsp <- vapply(candidate_lists, function(z) {
  z <- canonical_gene_symbol(z)
  any(z %in% hsp_clients)
}, logical(1)) & rna_gene_map$mapping_status == "ambiguous"

if (any(ambiguous_hsp)) {
  write_out(rna_gene_map[ambiguous_hsp, ], "ERROR_ambiguous_Hsp60_10_RNA_mappings.csv")
  stop("At least one ambiguous Ensembl->SYMBOL mapping touches an Hsp60/10 client.", call. = FALSE)
}

duplicated_hsp_symbols <- rna_gene_map |>
  dplyr::filter(.data$mapping_status == "unique", .data$gene_symbol %in% hsp_clients) |>
  dplyr::count(.data$gene_symbol, name = "n_ensembl_features") |>
  dplyr::filter(.data$n_ensembl_features > 1)

if (nrow(duplicated_hsp_symbols) > 0) {
  write_out(duplicated_hsp_symbols, "ERROR_Hsp60_10_symbols_with_multiple_ensembl_features.csv")
  stop("At least one Hsp60/10 symbol maps uniquely from multiple modeled Ensembl features.", call. = FALSE)
}

extract_rna_contrast <- function(contrast_name, numerator, denominator) {
  res <- DESeq2::results(
    dds,
    contrast = c("clinical_stage", numerator, denominator),
    alpha = 0.05,
    independentFiltering = TRUE
  )

  as.data.frame(res) |>
    tibble::rownames_to_column("ensembl_id_version") |>
    dplyr::left_join(rna_gene_map, by = "ensembl_id_version") |>
    dplyr::transmute(
      modality = "RNA",
      universe = "DESeq2_count_filtered",
      contrast = contrast_name,
      feature_id = .data$ensembl_id_version,
      ensembl_id = .data$ensembl_id,
      gene_symbol = .data$gene_symbol,
      mapping_status = .data$mapping_status,
      base_mean = as.numeric(.data$baseMean),
      effect = as.numeric(.data$log2FoldChange),
      effect_label = "DESeq2 log2 fold change",
      std_error = as.numeric(.data$lfcSE),
      conf_low = .data$effect - 1.96 * .data$std_error,
      conf_high = .data$effect + 1.96 * .data$std_error,
      statistic = as.numeric(.data$stat),
      p_value = as.numeric(.data$pvalue),
      fdr = as.numeric(.data$padj),
      significant_fdr05 = is.finite(.data$fdr) & .data$fdr < 0.05,
      is_hsp60_10_client = .data$mapping_status == "unique" & .data$gene_symbol %in% hsp_clients
    )
}

rna_results <- dplyr::bind_rows(lapply(seq_len(nrow(contrast_specs)), function(i) {
  extract_rna_contrast(
    contrast_specs$contrast[i],
    contrast_specs$numerator[i],
    contrast_specs$denominator[i]
  )
}))

write_out(rna_results, "04_RNA_DESeq2_all_contrasts.csv")

message("\n============================================================")
message("PROTEIN DIFFERENTIAL ABUNDANCE: LIMMA")
message("============================================================")

protein_mat <- as.matrix(prot_mat_raw)
mode(protein_mat) <- "numeric"
if (ncol(protein_mat) != 400L) stop("Expected 400 columns in production protein matrix.", call. = FALSE)

protein_symbols <- canonical_gene_symbol(rownames(protein_mat))
if (anyDuplicated(protein_symbols) > 0) stop("Duplicate canonical protein symbols detected.", call. = FALSE)
rownames(protein_mat) <- protein_symbols

prot_meta_all <- prot_meta_adj[match(colnames(protein_mat), prot_meta_adj$SampleID), , drop = FALSE]
if (!all(colnames(protein_mat) == prot_meta_all$SampleID)) stop("Protein matrix/metadata alignment failed.", call. = FALSE)
if (anyNA(prot_meta_all$batch_factor)) stop("Protein batch_factor incomplete across 400 samples.", call. = FALSE)

batch_factor_all <- droplevels(factor(prot_meta_all$batch_factor))
batch_indices <- split(seq_along(batch_factor_all), batch_factor_all)

present_by_batch <- sapply(batch_indices, function(idx) {
  rowSums(is.finite(protein_mat[, idx, drop = FALSE])) > 0
})

n_missing_batches <- rowSums(!present_by_batch)
missing_fraction <- rowMeans(!is.finite(protein_mat))
complete_all_400 <- rowSums(is.finite(protein_mat)) == ncol(protein_mat)

protein_universe_tbl <- tibble::tibble(
  gene_symbol = rownames(protein_mat),
  complete_all_400 = complete_all_400,
  n_missing_batches = n_missing_batches,
  missing_fraction = missing_fraction,
  primary_complete400 = complete_all_400,
  sensitivity_missing_batches_le2 = n_missing_batches <= 2,
  sensitivity_missingness_le20 = missing_fraction <= 0.20
)

write_out(protein_universe_tbl, "05_protein_feature_universe_manifest.csv")

universe_counts <- tibble::tibble(
  universe = c(
    "primary_complete400",
    "sensitivity_missing_batches_le2",
    "sensitivity_missingness_le20"
  ),
  n_features = c(
    sum(protein_universe_tbl$primary_complete400),
    sum(protein_universe_tbl$sensitivity_missing_batches_le2),
    sum(protein_universe_tbl$sensitivity_missingness_le20)
  )
)

expected_universe_counts <- c(
  primary_complete400 = 5088L,
  sensitivity_missing_batches_le2 = 5899L,
  sensitivity_missingness_le20 = 7159L
)

for (nm in names(expected_universe_counts)) {
  observed <- universe_counts$n_features[universe_counts$universe == nm]
  if (length(observed) != 1L || observed != expected_universe_counts[[nm]]) {
    stop("Protein universe count changed for ", nm, ": observed ", observed, "; expected ", expected_universe_counts[[nm]], ".", call. = FALSE)
  }
}

write_out(universe_counts, "06_protein_feature_universe_counts.csv")

prot_meta_model <- prot_meta_all |>
  dplyr::filter(.data$clinical_stage %in% stage_levels) |>
  dplyr::mutate(
    SampleID = as.character(.data$SampleID),
    clinical_stage = factor(as.character(.data$clinical_stage), levels = stage_levels)
  )

protein_complete_covars <- complete.cases(
  prot_meta_model[, c(protein_covars, "clinical_stage"), drop = FALSE]
)
prot_meta_model <- prot_meta_model[protein_complete_covars, , drop = FALSE]

expected_protein_counts <- c(NCI = 167L, MCI = 96L, AD = 109L)
observed_protein_counts <- table(prot_meta_model$clinical_stage)
if (!identical(as.integer(observed_protein_counts[stage_levels]), as.integer(expected_protein_counts[stage_levels]))) {
  stop("Protein canonical nuisance-complete stage counts changed.", call. = FALSE)
}

for (nm in protein_covars) {
  if (is.character(prot_meta_model[[nm]])) prot_meta_model[[nm]] <- factor(prot_meta_model[[nm]])
  if (is.factor(prot_meta_model[[nm]])) prot_meta_model[[nm]] <- droplevels(prot_meta_model[[nm]])
}

protein_formula <- stats::as.formula(
  paste("~ 0 +", paste(c("clinical_stage", protein_covars), collapse = " + "))
)

protein_design <- stats::model.matrix(protein_formula, data = prot_meta_model)
if (qr(protein_design)$rank != ncol(protein_design)) stop("Protein limma design matrix is not full rank.", call. = FALSE)

needed_stage_coef <- paste0("clinical_stage", stage_levels)
if (!all(needed_stage_coef %in% colnames(protein_design))) stop("Protein design missing stage coefficients.", call. = FALSE)

protein_contrast_matrix <- limma::makeContrasts(
  MCI_vs_NCI = clinical_stageMCI - clinical_stageNCI,
  AD_vs_NCI = clinical_stageAD - clinical_stageNCI,
  AD_vs_MCI = clinical_stageAD - clinical_stageMCI,
  levels = protein_design
)

run_protein_universe <- function(universe_name, keep_vector, primary = FALSE) {
  mat <- protein_mat[keep_vector, prot_meta_model$SampleID, drop = FALSE]
  if (!all(colnames(mat) == prot_meta_model$SampleID)) stop("Protein alignment failed for universe: ", universe_name, call. = FALSE)

  warning_messages <- character(0)

  fit <- withCallingHandlers(
    limma::lmFit(mat, protein_design),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  coefficient_na <- is.na(fit$coefficients)

  n_na_stage <- sum(
    rowSums(coefficient_na[, needed_stage_coef, drop = FALSE]) > 0
  )

  fit2 <- limma::contrasts.fit(fit, protein_contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  all_contrasts <- dplyr::bind_rows(lapply(contrast_specs$contrast, function(contrast_i) {
    tt <- limma::topTable(
      fit2,
      coef = contrast_i,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH",
      confint = 0.95
    ) |>
      tibble::rownames_to_column("gene_symbol")

    tt |>
      dplyr::transmute(
        modality = "Protein",
        universe = universe_name,
        contrast = contrast_i,
        feature_id = canonical_gene_symbol(.data$gene_symbol),
        ensembl_id = NA_character_,
        gene_symbol = canonical_gene_symbol(.data$gene_symbol),
        mapping_status = "gene_symbol_primary",
        base_mean = rowMeans(mat[.data$gene_symbol, , drop = FALSE], na.rm = TRUE),
        effect = as.numeric(.data$logFC),
        effect_label = "limma log2 TMT abundance difference",
        std_error = abs(.data$effect / as.numeric(.data$t)),
        conf_low = as.numeric(.data$CI.L),
        conf_high = as.numeric(.data$CI.R),
        statistic = as.numeric(.data$t),
        p_value = as.numeric(.data$P.Value),
        fdr = as.numeric(.data$adj.P.Val),
        significant_fdr05 = is.finite(.data$fdr) & .data$fdr < 0.05,
        is_hsp60_10_client = .data$gene_symbol %in% hsp_clients
      )
  }))

  contrast_estimability <- all_contrasts |>
    dplyr::group_by(.data$universe, .data$contrast) |>
    dplyr::summarise(
      n_features = dplyr::n(),
      n_effect_finite = sum(is.finite(.data$effect)),
      n_p_finite = sum(is.finite(.data$p_value)),
      n_fdr_finite = sum(is.finite(.data$fdr)),
      .groups = "drop"
    )

  if (
    primary &&
    (
      n_na_stage != 0 ||
      any(contrast_estimability$n_effect_finite != contrast_estimability$n_features) ||
      any(contrast_estimability$n_p_finite != contrast_estimability$n_features) ||
      any(contrast_estimability$n_fdr_finite != contrast_estimability$n_features)
    )
  ) {
    stop("Primary protein universe has non-estimable stage coefficients or contrasts.", call. = FALSE)
  }

  warning_audit <- tibble::tibble(
    universe = universe_name,
    primary = primary,
    n_features = nrow(mat),
    n_features_with_any_na_coefficient = sum(rowSums(coefficient_na) > 0),
    n_features_with_na_stage_coefficient = n_na_stage,
    n_warning_messages = length(unique(warning_messages)),
    warning_messages = paste(unique(warning_messages), collapse = " | ")
  )

  list(
    results = all_contrasts,
    warning_audit = warning_audit,
    estimability = contrast_estimability
  )
}

protein_primary <- run_protein_universe(
  "primary_complete400",
  protein_universe_tbl$primary_complete400,
  primary = TRUE
)

protein_batch2 <- run_protein_universe(
  "sensitivity_missing_batches_le2",
  protein_universe_tbl$sensitivity_missing_batches_le2,
  primary = FALSE
)

protein_missing20 <- run_protein_universe(
  "sensitivity_missingness_le20",
  protein_universe_tbl$sensitivity_missingness_le20,
  primary = FALSE
)

write_out(protein_primary$results, "07_protein_PRIMARY_complete400_limma.csv")
write_out(protein_batch2$results, "08_protein_SENSITIVITY_missing_batches_le2_limma.csv")
write_out(protein_missing20$results, "09_protein_SENSITIVITY_missingness_le20_limma.csv")

protein_warning_audit <- dplyr::bind_rows(
  protein_primary$warning_audit,
  protein_batch2$warning_audit,
  protein_missing20$warning_audit
)

protein_estimability_audit <- dplyr::bind_rows(
  protein_primary$estimability,
  protein_batch2$estimability,
  protein_missing20$estimability
)

write_out(protein_warning_audit, "10_protein_limma_warning_audit.csv")
write_out(protein_estimability_audit, "11_protein_contrast_estimability_audit.csv")

summarize_results <- function(tbl) {
  tbl |>
    dplyr::group_by(.data$modality, .data$universe, .data$contrast) |>
    dplyr::summarise(
      n_features = dplyr::n(),
      n_fdr_testable = sum(is.finite(.data$fdr)),
      n_fdr_significant = sum(.data$significant_fdr05, na.rm = TRUE),
      n_significant_higher = sum(.data$significant_fdr05 & .data$effect > 0, na.rm = TRUE),
      n_significant_lower = sum(.data$significant_fdr05 & .data$effect < 0, na.rm = TRUE),
      n_hsp_clients_tested = dplyr::n_distinct(.data$gene_symbol[.data$is_hsp60_10_client & !is.na(.data$gene_symbol)]),
      n_hsp_clients_fdr_significant = dplyr::n_distinct(.data$gene_symbol[.data$is_hsp60_10_client & .data$significant_fdr05 & !is.na(.data$gene_symbol)]),
      n_hsp_clients_higher = dplyr::n_distinct(.data$gene_symbol[.data$is_hsp60_10_client & .data$significant_fdr05 & .data$effect > 0 & !is.na(.data$gene_symbol)]),
      n_hsp_clients_lower = dplyr::n_distinct(.data$gene_symbol[.data$is_hsp60_10_client & .data$significant_fdr05 & .data$effect < 0 & !is.na(.data$gene_symbol)]),
      .groups = "drop"
    )
}

primary_summary <- dplyr::bind_rows(
  summarize_results(rna_results),
  summarize_results(protein_primary$results)
)

protein_sensitivity_summary <- summarize_results(
  dplyr::bind_rows(
    protein_primary$results,
    protein_batch2$results,
    protein_missing20$results
  )
)

write_out(primary_summary, "12_PRIMARY_reviewer_contrast_summary.csv")
write_out(protein_sensitivity_summary, "13_protein_sensitivity_contrast_summary.csv")

primary_combined <- dplyr::bind_rows(
  rna_results,
  protein_primary$results
)

hsp_primary <- primary_combined |>
  dplyr::filter(.data$is_hsp60_10_client) |>
  dplyr::arrange(
    factor(.data$contrast, levels = contrast_specs$contrast),
    .data$modality,
    .data$fdr
  )

write_out(hsp_primary, "14_PRIMARY_Hsp60_10_client_differential_results.csv")

core_chaperonin <- primary_combined |>
  dplyr::filter(.data$gene_symbol %in% c("HSPD1", "HSPE1")) |>
  dplyr::arrange(
    factor(.data$contrast, levels = contrast_specs$contrast),
    .data$modality,
    .data$gene_symbol
  )

write_out(core_chaperonin, "15_PRIMARY_HSPD1_HSPE1_results.csv")

best_unique_symbol_row <- function(tbl, modality_prefix) {
  tbl |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol)) |>
    dplyr::group_by(.data$contrast, .data$gene_symbol) |>
    dplyr::arrange(
      dplyr::if_else(is.finite(.data$fdr), .data$fdr, Inf),
      dplyr::desc(.data$base_mean)
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(
      "contrast", "gene_symbol", "effect", "std_error",
      "conf_low", "conf_high", "p_value", "fdr", "significant_fdr05"
    ) |>
    dplyr::rename_with(
      ~ paste0(modality_prefix, "_", .x),
      -c("contrast", "gene_symbol")
    )
}

rna_by_symbol <- best_unique_symbol_row(rna_results, "rna")
protein_by_symbol <- best_unique_symbol_row(protein_primary$results, "protein")

cross_modal_primary <- dplyr::inner_join(
  rna_by_symbol,
  protein_by_symbol,
  by = c("contrast", "gene_symbol")
) |>
  dplyr::mutate(
    is_hsp60_10_client = .data$gene_symbol %in% hsp_clients,
    significance_class = dplyr::case_when(
      .data$rna_significant_fdr05 & .data$protein_significant_fdr05 ~ "FDR-significant in both",
      .data$rna_significant_fdr05 ~ "RNA only",
      .data$protein_significant_fdr05 ~ "Protein only",
      TRUE ~ "Neither"
    ),
    direction_concordant = sign(.data$rna_effect) == sign(.data$protein_effect)
  )

write_out(cross_modal_primary, "16_PRIMARY_cross_modal_shared_features_side_by_side_NO_SUBTRACTION.csv")

cross_modal_hsp_primary <- cross_modal_primary |>
  dplyr::filter(.data$is_hsp60_10_client)

write_out(cross_modal_hsp_primary, "17_PRIMARY_cross_modal_Hsp60_10_side_by_side_NO_SUBTRACTION.csv")

cross_modal_hsp_summary <- cross_modal_hsp_primary |>
  dplyr::group_by(.data$contrast) |>
  dplyr::summarise(
    n_shared_hsp_clients = dplyr::n(),
    n_rna_only_fdr = sum(.data$significance_class == "RNA only"),
    n_protein_only_fdr = sum(.data$significance_class == "Protein only"),
    n_both_fdr = sum(.data$significance_class == "FDR-significant in both"),
    n_neither_fdr = sum(.data$significance_class == "Neither"),
    n_direction_concordant = sum(.data$direction_concordant, na.rm = TRUE),
    pct_direction_concordant = 100 * mean(.data$direction_concordant, na.rm = TRUE),
    .groups = "drop"
  )

write_out(cross_modal_hsp_summary, "18_PRIMARY_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv")

# ============================================================
# 13b. Hsp60/10 primary-universe bookkeeping manifest
# ============================================================

rna_hsp_primary <- unique(
  rna_results$gene_symbol[
    rna_results$is_hsp60_10_client &
      !is.na(rna_results$gene_symbol)
  ]
)

protein_hsp_primary <- unique(
  protein_primary$results$gene_symbol[
    protein_primary$results$is_hsp60_10_client &
      !is.na(protein_primary$results$gene_symbol)
  ]
)

hsp_primary_universe_manifest <- tibble::tibble(
  gene_symbol = sort(hsp_clients),
  in_rna_primary =
    gene_symbol %in% rna_hsp_primary,
  in_protein_primary =
    gene_symbol %in% protein_hsp_primary
) |>
  dplyr::mutate(
    category = dplyr::case_when(
      .data$in_rna_primary &
        .data$in_protein_primary ~ "Both",
      .data$in_rna_primary ~ "RNA only",
      .data$in_protein_primary ~ "Protein only",
      TRUE ~ "Neither"
    )
  )

hsp_primary_universe_counts <- hsp_primary_universe_manifest |>
  dplyr::count(
    .data$category,
    name = "n"
  ) |>
  dplyr::mutate(
    category = factor(
      .data$category,
      levels = c(
        "Both",
        "RNA only",
        "Protein only",
        "Neither"
      )
    )
  ) |>
  dplyr::arrange(.data$category) |>
  dplyr::mutate(
    category = as.character(.data$category)
  )

expected_hsp_primary_counts <- c(
  "Both" = 257L,
  "RNA only" = 40L,
  "Protein only" = 17L,
  "Neither" = 7L
)

observed_hsp_primary_counts <- stats::setNames(
  hsp_primary_universe_counts$n,
  hsp_primary_universe_counts$category
)

for (nm in names(expected_hsp_primary_counts)) {
  observed <- observed_hsp_primary_counts[[nm]]

  if (
    is.null(observed) ||
    length(observed) != 1L ||
    observed != expected_hsp_primary_counts[[nm]]
  ) {
    stop(
      "Hsp60/10 primary-universe count changed for ",
      nm,
      ": observed ",
      ifelse(is.null(observed), "missing", observed),
      "; expected ",
      expected_hsp_primary_counts[[nm]],
      ".",
      call. = FALSE
    )
  }
}

if (nrow(hsp_primary_universe_manifest) != 321L) {
  stop(
    "Expected 321 rows in Hsp60/10 primary-universe manifest; observed ",
    nrow(hsp_primary_universe_manifest),
    ".",
    call. = FALSE
  )
}

if (length(rna_hsp_primary) != 297L) {
  stop(
    "Expected 297 Hsp60/10 clients in RNA primary universe; observed ",
    length(rna_hsp_primary),
    ".",
    call. = FALSE
  )
}

if (length(protein_hsp_primary) != 274L) {
  stop(
    "Expected 274 Hsp60/10 clients in protein primary universe; observed ",
    length(protein_hsp_primary),
    ".",
    call. = FALSE
  )
}

write_out(
  hsp_primary_universe_manifest,
  "25_Hsp60_10_primary_universe_manifest.csv"
)

write_out(
  hsp_primary_universe_counts,
  "26_Hsp60_10_primary_universe_counts.csv"
)

get_primary_reference <- function(sensitivity_results) {
  dplyr::inner_join(
    protein_primary$results |>
      dplyr::select("contrast", "gene_symbol", primary_effect = "effect", primary_fdr = "fdr", primary_sig = "significant_fdr05"),
    sensitivity_results |>
      dplyr::select("contrast", "gene_symbol", sensitivity_effect = "effect", sensitivity_fdr = "fdr", sensitivity_sig = "significant_fdr05"),
    by = c("contrast", "gene_symbol")
  )
}

stability_one <- function(sensitivity_results, sensitivity_name) {
  joined <- get_primary_reference(sensitivity_results)
  dplyr::bind_rows(lapply(contrast_specs$contrast, function(contrast_i) {
    x <- joined |>
      dplyr::filter(.data$contrast == contrast_i)

    tibble::tibble(
      sensitivity_universe = sensitivity_name,
      contrast = contrast_i,
      n_shared_primary_features = nrow(x),
      effect_spearman = suppressWarnings(stats::cor(x$primary_effect, x$sensitivity_effect, method = "spearman", use = "complete.obs")),
      direction_agreement = mean(sign(x$primary_effect) == sign(x$sensitivity_effect), na.rm = TRUE),
      primary_fdr_significant = sum(x$primary_sig, na.rm = TRUE),
      sensitivity_fdr_significant_among_shared = sum(x$sensitivity_sig, na.rm = TRUE),
      primary_significant_retained = sum(x$primary_sig & x$sensitivity_sig, na.rm = TRUE)
    )
  }))
}

protein_sensitivity_stability <- dplyr::bind_rows(
  stability_one(protein_batch2$results, "sensitivity_missing_batches_le2"),
  stability_one(protein_missing20$results, "sensitivity_missingness_le20")
)

write_out(protein_sensitivity_stability, "19_protein_sensitivity_stability_vs_primary.csv")

sample_audit <- dplyr::bind_rows(
  tibble::tibble(
    modality = "RNA",
    clinical_stage = names(expected_rna_counts),
    n = as.integer(expected_rna_counts),
    eligibility = "canonical stage + production RNA covariates complete"
  ),
  tibble::tibble(
    modality = "Protein",
    clinical_stage = names(expected_protein_counts),
    n = as.integer(expected_protein_counts),
    eligibility = "canonical stage + production protein covariates complete"
  )
)

write_out(sample_audit, "20_sample_stage_audit.csv")

policy_audit <- tibble::tibble(
  modality = c("RNA", "Protein", "Protein", "Protein"),
  role = c("primary", "primary", "sensitivity", "sensitivity"),
  universe = c(
    "DESeq2_count_filtered",
    "primary_complete400",
    "sensitivity_missing_batches_le2",
    "sensitivity_missingness_le20"
  ),
  definition = c(
    "raw counts; count >=10 in >=3 samples",
    "finite TMT abundance in all 400 participants",
    "absent from no more than 2 whole TMT batches",
    "participant-level missingness <=20%"
  ),
  n_features = c(
    rna_genes_after_filter,
    sum(protein_universe_tbl$primary_complete400),
    sum(protein_universe_tbl$sensitivity_missing_batches_le2),
    sum(protein_universe_tbl$sensitivity_missingness_le20)
  ),
  method = c("DESeq2", "limma", "limma", "limma")
)

write_out(policy_audit, "21_analysis_policy_audit.csv")

model_audit <- tibble::tibble(
  modality = c("RNA", "Protein"),
  formula = c(
    paste(deparse(rna_formula), collapse = ""),
    paste(deparse(protein_formula), collapse = "")
  ),
  biological_variable = "clinical_stage",
  stage_levels = "NCI;MCI;AD",
  contrasts = "MCI_vs_NCI;AD_vs_NCI;AD_vs_MCI",
  multiple_testing = "Benjamini-Hochberg FDR",
  significance = "FDR < 0.05",
  cross_modal_policy = "report side-by-side; no RNA-protein effect subtraction"
)

write_out(model_audit, "22_model_specification_audit.csv")

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g"))
}

report <- c(
  "CONVENTIONAL STAGE DIFFERENTIAL ANALYSIS",
  "========================================",
  "",
  "Purpose:",
  "  Conventional RNA differential-expression and protein differential-abundance",
  "  analyses across canonical NCI, MCI, and AD stages.",
  "",
  "Cross-modal policy:",
  "  RNA and protein effects are modeled and reported separately.",
  "  No RNA-minus-protein or protein-minus-RNA effect is calculated.",
  "",
  paste0("RNA primary universe: ", rna_genes_after_filter, " count-filtered features."),
  paste0("Protein primary universe: ", sum(protein_universe_tbl$primary_complete400), " proteins observed in all 400 participants."),
  paste0("Protein sensitivity universe 1: ", sum(protein_universe_tbl$sensitivity_missing_batches_le2), " proteins absent from <=2 whole TMT batches."),
  paste0("Protein sensitivity universe 2: ", sum(protein_universe_tbl$sensitivity_missingness_le20), " proteins with <=20% participant-level missingness."),
  "",
  "PRIMARY CONTRAST SUMMARY:"
)

for (i in seq_len(nrow(primary_summary))) {
  x <- primary_summary[i, ]
  report <- c(
    report,
    paste0(
      "  ", x$modality, " | ", x$contrast,
      " | features=", x$n_features,
      " | FDR-testable=", x$n_fdr_testable,
      " | FDR<0.05=", x$n_fdr_significant,
      " | Hsp60/10 FDR<0.05=",
      x$n_hsp_clients_fdr_significant, "/", x$n_hsp_clients_tested
    )
  )
}

report <- c(
  report,
  "",
  "HSPD1 / HSPE1 (study-specific anchor genes; not a reviewer-requested subset):"
)

core_display <- core_chaperonin |>
  dplyr::select("modality", "contrast", "gene_symbol", "effect", "conf_low", "conf_high", "p_value", "fdr")

for (i in seq_len(nrow(core_display))) {
  x <- core_display[i, ]
  report <- c(
    report,
    paste0(
      "  ", x$modality, " | ", x$contrast, " | ", x$gene_symbol,
      " | effect=", fmt(x$effect),
      " | 95% CI ", fmt(x$conf_low), " to ", fmt(x$conf_high),
      " | P=", fmt(x$p_value),
      " | FDR=", fmt(x$fdr)
    )
  )
}

report <- c(
  report,
  "",
  "Interpretation guardrails:",
  "  - Protein primary inference uses only proteins complete across all 400 TMT participants.",
  "  - Broader protein universes are sensitivity analyses.",
  "  - Cross-modal direction concordance is descriptive only.",
  "  - No cross-modal subtraction is used as an inferential statistic.",
  ""
)

report_path <- file.path(out_dir, "23_REVIEWER_REPORT_SUMMARY.txt")
writeLines(report, report_path)
message("Wrote: ", report_path)

validation <- tibble::tribble(
  ~check, ~observed, ~expected, ~passed,
  "RNA NCI", as.character(observed_rna_counts["NCI"]), "200", observed_rna_counts["NCI"] == 200,
  "RNA MCI", as.character(observed_rna_counts["MCI"]), "158", observed_rna_counts["MCI"] == 158,
  "RNA AD", as.character(observed_rna_counts["AD"]), "219", observed_rna_counts["AD"] == 219,
  "Protein NCI", as.character(observed_protein_counts["NCI"]), "167", observed_protein_counts["NCI"] == 167,
  "Protein MCI", as.character(observed_protein_counts["MCI"]), "96", observed_protein_counts["MCI"] == 96,
  "Protein AD", as.character(observed_protein_counts["AD"]), "109", observed_protein_counts["AD"] == 109,
  "TMT batches", as.character(length(batch_indices)), "50", length(batch_indices) == 50,
  "Primary protein universe", as.character(sum(protein_universe_tbl$primary_complete400)), "5088", sum(protein_universe_tbl$primary_complete400) == 5088,
  "Batch<=2 sensitivity universe", as.character(sum(protein_universe_tbl$sensitivity_missing_batches_le2)), "5899", sum(protein_universe_tbl$sensitivity_missing_batches_le2) == 5899,
  "<=20% missing sensitivity universe", as.character(sum(protein_universe_tbl$sensitivity_missingness_le20)), "7159", sum(protein_universe_tbl$sensitivity_missingness_le20) == 7159,
  "Primary protein NA stage coefficients", as.character(protein_primary$warning_audit$n_features_with_na_stage_coefficient), "0", protein_primary$warning_audit$n_features_with_na_stage_coefficient == 0,
  "Primary protein contrasts fully testable", as.character(all(protein_primary$estimability$n_fdr_finite == protein_primary$estimability$n_features)), "TRUE", all(protein_primary$estimability$n_fdr_finite == protein_primary$estimability$n_features),
  "RNA contrasts", as.character(dplyr::n_distinct(rna_results$contrast)), "3", dplyr::n_distinct(rna_results$contrast) == 3,
  "Primary protein contrasts", as.character(dplyr::n_distinct(protein_primary$results$contrast)), "3", dplyr::n_distinct(protein_primary$results$contrast) == 3,
  "No subtraction columns in cross-modal table", as.character(!any(grepl("minus|difference", colnames(cross_modal_primary), ignore.case = TRUE))), "TRUE", !any(grepl("minus|difference", colnames(cross_modal_primary), ignore.case = TRUE)),
  "Canonical Hsp60/10 manifest total", as.character(nrow(hsp_primary_universe_manifest)), "321", nrow(hsp_primary_universe_manifest) == 321,
  "RNA primary Hsp60/10 total", as.character(length(rna_hsp_primary)), "297", length(rna_hsp_primary) == 297,
  "Protein primary Hsp60/10 total", as.character(length(protein_hsp_primary)), "274", length(protein_hsp_primary) == 274,
  "Shared primary Hsp60/10 total", as.character(observed_hsp_primary_counts[["Both"]]), "257", observed_hsp_primary_counts[["Both"]] == 257,
  "RNA-only primary Hsp60/10 total", as.character(observed_hsp_primary_counts[["RNA only"]]), "40", observed_hsp_primary_counts[["RNA only"]] == 40,
  "Protein-only primary Hsp60/10 total", as.character(observed_hsp_primary_counts[["Protein only"]]), "17", observed_hsp_primary_counts[["Protein only"]] == 17,
  "Neither-primary Hsp60/10 total", as.character(observed_hsp_primary_counts[["Neither"]]), "7", observed_hsp_primary_counts[["Neither"]] == 7
)

write_out(validation, "24_validation_checks.csv")

if (any(!validation$passed)) stop("FINAL conventional differential-analysis validation failed.", call. = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

git_head <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE), error = function(e) NA_character_)
git_status <- tryCatch(system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE), error = function(e) NA_character_)

writeLines(
  c(
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Git HEAD: ", paste(git_head, collapse = "")),
    "Git status --short:",
    git_status,
    "",
    paste0("Controlled RNA counts input: ", counts_path),
    paste0("RNA formula: ", paste(deparse(rna_formula), collapse = "")),
    paste0("Protein formula: ", paste(deparse(protein_formula), collapse = "")),
    "Primary protein universe: complete in all 400 participants",
    "Sensitivity protein universe 1: <=2 whole TMT batches absent",
    "Sensitivity protein universe 2: <=20% participant-level missingness",
    "No imputation",
    "Multiple testing: Benjamini-Hochberg FDR",
    "Significance: FDR < 0.05",
    "Cross-modal policy: side-by-side reporting; no effect subtraction"
  ),
  file.path(out_dir, "run_provenance.txt")
)

message("\n============================================================")
message("FINAL CONVENTIONAL STAGE DIFFERENTIAL ANALYSIS COMPLETE")
message("============================================================")
message("\nPrimary reviewer summary:")
print(tibble::as_tibble(primary_summary), n = Inf)
message("\nProtein sensitivity summary:")
print(tibble::as_tibble(protein_sensitivity_summary), n = Inf)
message("\nProtein warning audit:")
print(tibble::as_tibble(protein_warning_audit), n = Inf)
message("\nHSPD1 / HSPE1:")
print(tibble::as_tibble(core_display), n = Inf)
message("\nCross-modal Hsp60/10 primary summary:")
print(tibble::as_tibble(cross_modal_hsp_summary), n = Inf)
message("\nReviewer report:")
message(report_path)
