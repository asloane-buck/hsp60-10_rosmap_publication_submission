#!/usr/bin/env Rscript

############################################################
## 21_run_conventional_stage_differential_analysis.R
##
## Reviewer-requested conventional differential analyses
## across canonical NCI / MCI / AD stages.
##
## RNA:
##   DESeq2 on raw bulk RNA-seq counts.
##
## Protein:
##   limma on processed log2 TMT abundance.
##
## Primary contrasts:
##   MCI vs NCI
##   AD  vs NCI
##   AD  vs MCI
##
## IMPORTANT:
## - Uses the current production clinical_stage definitions.
## - Uses the same nuisance-covariate sets selected by production.
## - RNA and protein are modeled separately.
## - No RNA-minus-protein or protein-minus-RNA subtraction is performed.
##
## Controlled input:
##   Raw ROSMAP RNA-seq counts are not redistributed.
##   Supply them with:
##
##   Rscript R/reviewer_revisions/21_run_conventional_stage_differential_analysis.R \
##     --counts ~/Downloads/count_matrix.rds
##
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

source("R/00_config.R")
source("R/01_utils.R")

# ------------------------------------------------------------
# 0. CLI and paths
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[1] == length(args)) {
    stop("Missing value after ", flag, call. = FALSE)
  }
  args[hit[1] + 1]
}

counts_path <- get_arg(
  "--counts",
  Sys.getenv("ROSMAP_RAW_COUNTS_FILE", unset = "")
)

if (!nzchar(counts_path)) {
  stop(
    "Raw RNA-seq counts are required. Supply --counts /path/to/count_matrix.rds ",
    "or set ROSMAP_RAW_COUNTS_FILE.",
    call. = FALSE
  )
}

counts_path <- path.expand(counts_path)

if (!file.exists(counts_path)) {
  stop("Raw RNA-seq counts file not found: ", counts_path, call. = FALSE)
}

out_dir <- get_arg(
  "--outdir",
  file.path(
    cfg$project_dir,
    "outputs",
    "reviewer_revisions",
    "conventional_stage_differential"
  )
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_out <- function(x, filename) {
  path <- file.path(out_dir, filename)
  readr::write_csv(x, path, na = "")
  message("Wrote: ", path)
  invisible(x)
}

# ------------------------------------------------------------
# 1. Load current production objects
# ------------------------------------------------------------

load_required_object <- function(name) {
  path <- file.path(cfg$object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing production object: ", path, call. = FALSE)
  }
  readRDS(path)
}

rna_meta_adj <- load_required_object("rna_meta_adj")
rna_covars <- load_required_object("rna_covars")
prot_mat_raw <- load_required_object("prot_mat_raw")
prot_meta_adj <- load_required_object("prot_meta_adj")
protein_covars <- load_required_object("protein_covars")
all_hsp60_10_client_tbl <- load_required_object("all_hsp60_10_client_tbl")

require_columns(
  rna_meta_adj,
  c("sample_id", "clinical_stage"),
  "RNA production metadata"
)

require_columns(
  prot_meta_adj,
  c("SampleID", "clinical_stage"),
  "protein production metadata"
)

rna_covars <- unique(as.character(rna_covars))
protein_covars <- unique(as.character(protein_covars))

if (length(rna_covars) == 0) {
  stop("Production RNA covariate set is empty.", call. = FALSE)
}
if (length(protein_covars) == 0) {
  stop("Production protein covariate set is empty.", call. = FALSE)
}

missing_rna_covars <- setdiff(rna_covars, colnames(rna_meta_adj))
missing_prot_covars <- setdiff(protein_covars, colnames(prot_meta_adj))

if (length(missing_rna_covars) > 0) {
  stop(
    "RNA production covariate(s) missing from rna_meta_adj: ",
    paste(missing_rna_covars, collapse = ", "),
    call. = FALSE
  )
}
if (length(missing_prot_covars) > 0) {
  stop(
    "Protein production covariate(s) missing from prot_meta_adj: ",
    paste(missing_prot_covars, collapse = ", "),
    call. = FALSE
  )
}

stage_levels <- c("NCI", "MCI", "AD")

# ------------------------------------------------------------
# 2. Canonical Hsp60/10 client set
# ------------------------------------------------------------

client_gene_col <- intersect(
  c("gene", "gene_symbol", "symbol"),
  colnames(all_hsp60_10_client_tbl)
)[1]

if (is.na(client_gene_col)) {
  stop(
    "Could not identify Hsp60/10 client gene column. Available columns: ",
    paste(colnames(all_hsp60_10_client_tbl), collapse = ", "),
    call. = FALSE
  )
}

hsp_clients <- all_hsp60_10_client_tbl[[client_gene_col]] |>
  as.character() |>
  canonical_gene_symbol() |>
  unique()

hsp_clients <- hsp_clients[!is.na(hsp_clients) & nzchar(hsp_clients)]

if (length(hsp_clients) != 321) {
  stop(
    "Expected 321 canonical Hsp60/10 clients; observed ",
    length(hsp_clients),
    ".",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 3. RNA: raw counts + canonical production metadata
# ------------------------------------------------------------

message("\n============================================================")
message("RNA DIFFERENTIAL EXPRESSION: DESeq2")
message("============================================================")

rna_counts <- readRDS(counts_path)

if (inherits(rna_counts, "data.frame")) {
  rna_counts <- as.matrix(rna_counts)
}

if (is.null(rownames(rna_counts)) || is.null(colnames(rna_counts))) {
  stop(
    "RNA count matrix must have gene row names and sample column names.",
    call. = FALSE
  )
}

## Harmonize raw-count sample identifiers exactly as in the original
## ROSMAP DESeq2 preprocessing: raw count columns may be full BAM paths
## and/or end in ".bam", whereas production metadata stores sample_id.
raw_count_sample_ids_original <- colnames(rna_counts)
raw_count_sample_ids_normalized <- basename(raw_count_sample_ids_original)
raw_count_sample_ids_normalized <- sub("\\.bam$", "", raw_count_sample_ids_normalized)

if (anyDuplicated(raw_count_sample_ids_normalized) > 0) {
  dup_ids <- unique(
    raw_count_sample_ids_normalized[
      duplicated(raw_count_sample_ids_normalized)
    ]
  )
  stop(
    "Normalizing raw-count sample IDs created duplicates: ",
    paste(head(dup_ids, 20), collapse = ", "),
    call. = FALSE
  )
}

colnames(rna_counts) <- raw_count_sample_ids_normalized

message(
  "Raw-count sample-ID harmonization: basename() + strip terminal .bam"
)

rna_meta_stage <- rna_meta_adj |>
  dplyr::filter(.data$clinical_stage %in% stage_levels) |>
  dplyr::mutate(
    sample_id = as.character(.data$sample_id),
    clinical_stage = factor(
      as.character(.data$clinical_stage),
      levels = stage_levels
    )
  )

if (anyDuplicated(rna_meta_stage$sample_id) > 0) {
  stop("Duplicate RNA sample_id in production metadata.", call. = FALSE)
}

expected_rna_stage_counts <- c(NCI = 200L, MCI = 158L, AD = 220L)
observed_rna_stage_counts <- table(rna_meta_stage$clinical_stage)

if (!identical(
  as.integer(observed_rna_stage_counts[stage_levels]),
  as.integer(expected_rna_stage_counts[stage_levels])
)) {
  stop(
    "RNA canonical stage counts changed. Observed: ",
    paste(
      stage_levels,
      as.integer(observed_rna_stage_counts[stage_levels]),
      sep = "=",
      collapse = ", "
    ),
    "; expected NCI=200, MCI=158, AD=220.",
    call. = FALSE
  )
}

col_overlap <- sum(rna_meta_stage$sample_id %in% colnames(rna_counts))
row_overlap <- sum(rna_meta_stage$sample_id %in% rownames(rna_counts))

message(
  "RNA production sample IDs matched after harmonization: ",
  col_overlap, "/", nrow(rna_meta_stage),
  " count-matrix columns; ",
  row_overlap, "/", nrow(rna_meta_stage),
  " count-matrix rows."
)

orientation <- NA_character_

if (col_overlap == nrow(rna_meta_stage) && row_overlap < col_overlap) {
  orientation <- "genes_rows_samples_columns"
} else if (row_overlap == nrow(rna_meta_stage) && col_overlap < row_overlap) {
  rna_counts <- t(rna_counts)
  orientation <- "samples_rows_genes_columns_transposed"
} else {
  stop(
    "Could not unambiguously align raw count matrix to the 578 production RNA samples. ",
    "Column overlap=", col_overlap,
    "; row overlap=", row_overlap,
    ".",
    call. = FALSE
  )
}

if (anyDuplicated(colnames(rna_counts)) > 0) {
  stop("Duplicate sample columns in raw RNA count matrix.", call. = FALSE)
}
if (anyDuplicated(rownames(rna_counts)) > 0) {
  stop("Duplicate gene rows in raw RNA count matrix.", call. = FALSE)
}

if (!all(rna_meta_stage$sample_id %in% colnames(rna_counts))) {
  stop(
    "One or more production RNA sample IDs are absent from the raw count matrix.",
    call. = FALSE
  )
}

rna_counts <- rna_counts[
  ,
  rna_meta_stage$sample_id,
  drop = FALSE
]

if (!all(colnames(rna_counts) == rna_meta_stage$sample_id)) {
  stop("Raw RNA counts and production metadata are not exactly aligned.", call. = FALSE)
}

if (anyNA(rna_counts)) {
  stop("Raw RNA count matrix contains NA values.", call. = FALSE)
}
if (any(rna_counts < 0)) {
  stop("Raw RNA count matrix contains negative values.", call. = FALSE)
}

max_noninteger <- max(abs(rna_counts - round(rna_counts)))
if (!is.finite(max_noninteger) || max_noninteger > 1e-6) {
  stop(
    "RNA count matrix is not integer-like; max distance to nearest integer = ",
    signif(max_noninteger, 6),
    ". DESeq2 requires raw counts.",
    call. = FALSE
  )
}

rna_complete <- complete.cases(
  rna_meta_stage[, c(rna_covars, "clinical_stage"), drop = FALSE]
)

if (!all(rna_complete)) {
  stop(
    "Current production RNA metadata has ",
    sum(!rna_complete),
    " incomplete rows for the production DESeq2 covariates. ",
    "Expected all 578 canonical-stage RNA samples to be complete.",
    call. = FALSE
  )
}

rna_meta_model <- rna_meta_stage

for (nm in rna_covars) {
  if (is.character(rna_meta_model[[nm]])) {
    rna_meta_model[[nm]] <- factor(rna_meta_model[[nm]])
  }
  if (is.factor(rna_meta_model[[nm]])) {
    rna_meta_model[[nm]] <- droplevels(rna_meta_model[[nm]])
  }
}

rna_formula <- stats::as.formula(
  paste(
    "~",
    paste(c(rna_covars, "clinical_stage"), collapse = " + ")
  )
)

rna_mm <- stats::model.matrix(rna_formula, data = rna_meta_model)

if (qr(rna_mm)$rank != ncol(rna_mm)) {
  stop(
    "RNA DESeq2 design matrix is not full rank. Formula: ",
    paste(deparse(rna_formula), collapse = ""),
    call. = FALSE
  )
}

rna_counts_model <- round(rna_counts)
storage.mode(rna_counts_model) <- "integer"

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = rna_counts_model,
  colData = rna_meta_model,
  design = rna_formula
)

rna_genes_before_filter <- nrow(dds)

rna_keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
dds <- dds[rna_keep, ]

rna_genes_after_filter <- nrow(dds)

if (rna_genes_after_filter < 10000) {
  stop(
    "Unexpectedly few RNA genes passed the prespecified count filter: ",
    rna_genes_after_filter,
    ".",
    call. = FALSE
  )
}

message(
  "RNA genes before filter: ", rna_genes_before_filter,
  "; after count>=10 in >=3 samples: ", rna_genes_after_filter
)

dds <- DESeq2::DESeq(dds, quiet = TRUE)

ensembl_id_version <- rownames(dds)
ensembl_id <- sub("\\..*$", "", ensembl_id_version)

gene_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) |>
  as.character()

gene_symbol <- canonical_gene_symbol(gene_symbol)

rna_gene_map <- tibble::tibble(
  ensembl_id_version = ensembl_id_version,
  ensembl_id = ensembl_id,
  gene_symbol = gene_symbol
)

write_out(
  rna_gene_map,
  "01_rna_ensembl_to_symbol_mapping.csv"
)

# ------------------------------------------------------------
# 4. Protein: limma on processed non-residualized TMT matrix
# ------------------------------------------------------------

message("\n============================================================")
message("PROTEIN DIFFERENTIAL ABUNDANCE: LIMMA")
message("============================================================")

protein_mat <- as.matrix(prot_mat_raw)
mode(protein_mat) <- "numeric"

if (is.null(rownames(protein_mat)) || is.null(colnames(protein_mat))) {
  stop("Protein matrix lacks row/column names.", call. = FALSE)
}

protein_genes <- canonical_gene_symbol(rownames(protein_mat))

if (anyDuplicated(protein_genes) > 0) {
  dup <- unique(protein_genes[duplicated(protein_genes)])
  stop(
    "Canonical protein gene symbols are duplicated: ",
    paste(head(dup, 20), collapse = ", "),
    call. = FALSE
  )
}

rownames(protein_mat) <- protein_genes

prot_meta_stage <- prot_meta_adj |>
  dplyr::filter(.data$clinical_stage %in% stage_levels) |>
  dplyr::mutate(
    SampleID = as.character(.data$SampleID),
    clinical_stage = factor(
      as.character(.data$clinical_stage),
      levels = stage_levels
    )
  )

if (anyDuplicated(prot_meta_stage$SampleID) > 0) {
  stop("Duplicate protein SampleID in production metadata.", call. = FALSE)
}

protein_complete <- complete.cases(
  prot_meta_stage[, c(protein_covars, "clinical_stage"), drop = FALSE]
)

prot_meta_model <- prot_meta_stage[protein_complete, , drop = FALSE]

observed_protein_stage_counts <- table(prot_meta_model$clinical_stage)
expected_protein_stage_counts <- c(NCI = 167L, MCI = 96L, AD = 109L)

if (!identical(
  as.integer(observed_protein_stage_counts[stage_levels]),
  as.integer(expected_protein_stage_counts[stage_levels])
)) {
  stop(
    "Protein nuisance-complete canonical stage counts changed. Observed: ",
    paste(
      stage_levels,
      as.integer(observed_protein_stage_counts[stage_levels]),
      sep = "=",
      collapse = ", "
    ),
    "; expected NCI=167, MCI=96, AD=109.",
    call. = FALSE
  )
}

if (!all(prot_meta_model$SampleID %in% colnames(protein_mat))) {
  stop("One or more modeled protein samples are missing from prot_mat_raw.", call. = FALSE)
}

protein_mat_model <- protein_mat[
  ,
  prot_meta_model$SampleID,
  drop = FALSE
]

if (!all(colnames(protein_mat_model) == prot_meta_model$SampleID)) {
  stop("Protein matrix and metadata are not exactly aligned.", call. = FALSE)
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
    paste(c("clinical_stage", protein_covars), collapse = " + ")
  )
)

protein_design <- stats::model.matrix(
  protein_formula,
  data = prot_meta_model
)

if (qr(protein_design)$rank != ncol(protein_design)) {
  stop(
    "Protein limma design matrix is not full rank. Formula: ",
    paste(deparse(protein_formula), collapse = ""),
    call. = FALSE
  )
}

needed_stage_coef <- paste0("clinical_stage", stage_levels)

if (!all(needed_stage_coef %in% colnames(protein_design))) {
  stop(
    "Protein design is missing stage coefficients: ",
    paste(
      setdiff(needed_stage_coef, colnames(protein_design)),
      collapse = ", "
    ),
    call. = FALSE
  )
}

protein_fit <- limma::lmFit(
  protein_mat_model,
  protein_design
)

protein_contrast_matrix <- limma::makeContrasts(
  MCI_vs_NCI = clinical_stageMCI - clinical_stageNCI,
  AD_vs_NCI = clinical_stageAD - clinical_stageNCI,
  AD_vs_MCI = clinical_stageAD - clinical_stageMCI,
  levels = protein_design
)

protein_fit2 <- limma::contrasts.fit(
  protein_fit,
  protein_contrast_matrix
)
protein_fit2 <- limma::eBayes(protein_fit2)

# ------------------------------------------------------------
# 5. Contrast extraction helpers
# ------------------------------------------------------------

contrast_specs <- tibble::tribble(
  ~contrast, ~numerator, ~denominator,
  "MCI_vs_NCI", "MCI", "NCI",
  "AD_vs_NCI",  "AD",  "NCI",
  "AD_vs_MCI",  "AD",  "MCI"
)

extract_rna_contrast <- function(contrast, numerator, denominator) {
  res <- DESeq2::results(
    dds,
    contrast = c(
      "clinical_stage",
      numerator,
      denominator
    ),
    alpha = 0.05,
    independentFiltering = TRUE
  )

  tbl <- as.data.frame(res) |>
    tibble::rownames_to_column("ensembl_id_version") |>
    dplyr::left_join(
      rna_gene_map,
      by = "ensembl_id_version"
    ) |>
    dplyr::transmute(
      modality = "RNA",
      contrast = contrast,
      feature_id = .data$ensembl_id_version,
      ensembl_id = .data$ensembl_id,
      gene_symbol = .data$gene_symbol,
      base_mean = as.numeric(.data$baseMean),
      effect = as.numeric(.data$log2FoldChange),
      effect_label = "DESeq2 log2 fold change",
      std_error = as.numeric(.data$lfcSE),
      conf_low = .data$effect - 1.96 * .data$std_error,
      conf_high = .data$effect + 1.96 * .data$std_error,
      statistic = as.numeric(.data$stat),
      p_value = as.numeric(.data$pvalue),
      fdr = as.numeric(.data$padj),
      significant_fdr05 =
        is.finite(.data$fdr) & .data$fdr < 0.05,
      direction = dplyr::case_when(
        !.data$significant_fdr05 ~ "Not FDR-significant",
        .data$effect > 0 ~ "Higher in numerator stage",
        .data$effect < 0 ~ "Lower in numerator stage",
        TRUE ~ "FDR-significant, zero estimate"
      ),
      is_hsp60_10_client = .data$gene_symbol %in% hsp_clients
    )

  tbl
}

extract_protein_contrast <- function(contrast) {
  tt <- limma::topTable(
    protein_fit2,
    coef = contrast,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH",
    confint = 0.95
  ) |>
    tibble::rownames_to_column("gene_symbol")

  tt |>
    dplyr::transmute(
      modality = "Protein",
      contrast = contrast,
      feature_id = canonical_gene_symbol(.data$gene_symbol),
      ensembl_id = NA_character_,
      gene_symbol = canonical_gene_symbol(.data$gene_symbol),
      base_mean = rowMeans(
        protein_mat_model[.data$gene_symbol, , drop = FALSE],
        na.rm = TRUE
      ),
      effect = as.numeric(.data$logFC),
      effect_label = "limma log2 TMT abundance difference",
      std_error = abs(.data$effect / as.numeric(.data$t)),
      conf_low = as.numeric(.data$CI.L),
      conf_high = as.numeric(.data$CI.R),
      statistic = as.numeric(.data$t),
      p_value = as.numeric(.data$P.Value),
      fdr = as.numeric(.data$adj.P.Val),
      significant_fdr05 =
        is.finite(.data$fdr) & .data$fdr < 0.05,
      direction = dplyr::case_when(
        !.data$significant_fdr05 ~ "Not FDR-significant",
        .data$effect > 0 ~ "Higher in numerator stage",
        .data$effect < 0 ~ "Lower in numerator stage",
        TRUE ~ "FDR-significant, zero estimate"
      ),
      is_hsp60_10_client = .data$gene_symbol %in% hsp_clients
    )
}

rna_results <- purrr::pmap_dfr(
  contrast_specs,
  extract_rna_contrast
)

protein_results <- purrr::map_dfr(
  contrast_specs$contrast,
  extract_protein_contrast
)

if (anyDuplicated(
  rna_results[, c("contrast", "feature_id")]
) > 0) {
  stop("Duplicate RNA contrast-feature rows.", call. = FALSE)
}

if (anyDuplicated(
  protein_results[, c("contrast", "feature_id")]
) > 0) {
  stop("Duplicate protein contrast-feature rows.", call. = FALSE)
}

write_out(
  rna_results,
  "02_RNA_DESeq2_all_contrasts.csv"
)

write_out(
  protein_results,
  "03_protein_limma_all_contrasts.csv"
)

# ------------------------------------------------------------
# 6. Reviewer-facing summaries
# ------------------------------------------------------------

summarise_contrast <- function(tbl) {
  tbl |>
    dplyr::group_by(.data$modality, .data$contrast) |>
    dplyr::summarise(
      n_features = dplyr::n(),
      n_fdr_testable = sum(is.finite(.data$fdr)),
      n_fdr_significant = sum(.data$significant_fdr05, na.rm = TRUE),
      n_significant_higher = sum(
        .data$significant_fdr05 & .data$effect > 0,
        na.rm = TRUE
      ),
      n_significant_lower = sum(
        .data$significant_fdr05 & .data$effect < 0,
        na.rm = TRUE
      ),
      n_hsp_clients_tested = dplyr::n_distinct(
        .data$gene_symbol[
          .data$is_hsp60_10_client &
            !is.na(.data$gene_symbol)
        ]
      ),
      n_hsp_clients_fdr_significant = dplyr::n_distinct(
        .data$gene_symbol[
          .data$is_hsp60_10_client &
            .data$significant_fdr05 &
            !is.na(.data$gene_symbol)
        ]
      ),
      n_hsp_clients_higher = dplyr::n_distinct(
        .data$gene_symbol[
          .data$is_hsp60_10_client &
            .data$significant_fdr05 &
            .data$effect > 0 &
            !is.na(.data$gene_symbol)
        ]
      ),
      n_hsp_clients_lower = dplyr::n_distinct(
        .data$gene_symbol[
          .data$is_hsp60_10_client &
            .data$significant_fdr05 &
            .data$effect < 0 &
            !is.na(.data$gene_symbol)
        ]
      ),
      .groups = "drop"
    )
}

contrast_summary <- dplyr::bind_rows(
  summarise_contrast(rna_results),
  summarise_contrast(protein_results)
) |>
  dplyr::arrange(
    factor(.data$contrast, levels = contrast_specs$contrast),
    .data$modality
  )

write_out(
  contrast_summary,
  "04_reviewer_contrast_summary.csv"
)

core_chaperonin_results <- dplyr::bind_rows(
  rna_results,
  protein_results
) |>
  dplyr::filter(.data$gene_symbol %in% c("HSPD1", "HSPE1")) |>
  dplyr::arrange(
    factor(.data$contrast, levels = contrast_specs$contrast),
    .data$modality,
    .data$gene_symbol
  )

write_out(
  core_chaperonin_results,
  "05_HSPD1_HSPE1_differential_results.csv"
)

hsp_client_results <- dplyr::bind_rows(
  rna_results,
  protein_results
) |>
  dplyr::filter(.data$is_hsp60_10_client) |>
  dplyr::arrange(
    factor(.data$contrast, levels = contrast_specs$contrast),
    .data$modality,
    .data$fdr,
    dplyr::desc(abs(.data$effect))
  )

write_out(
  hsp_client_results,
  "06_Hsp60_10_client_differential_results.csv"
)

top_hsp_hits <- hsp_client_results |>
  dplyr::filter(.data$significant_fdr05) |>
  dplyr::group_by(.data$modality, .data$contrast) |>
  dplyr::arrange(.data$fdr, dplyr::desc(abs(.data$effect))) |>
  dplyr::slice_head(n = 20) |>
  dplyr::ungroup()

write_out(
  top_hsp_hits,
  "07_top_FDR_significant_Hsp60_10_clients.csv"
)

# ------------------------------------------------------------
# 7. Side-by-side cross-modal results WITHOUT subtraction
# ------------------------------------------------------------

best_symbol_row <- function(tbl, modality_label) {
  tbl |>
    dplyr::filter(
      !is.na(.data$gene_symbol),
      nzchar(.data$gene_symbol)
    ) |>
    dplyr::group_by(.data$contrast, .data$gene_symbol) |>
    dplyr::arrange(
      dplyr::if_else(is.finite(.data$fdr), .data$fdr, Inf),
      dplyr::desc(.data$base_mean)
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(
      .data$contrast,
      .data$gene_symbol,
      effect,
      std_error,
      conf_low,
      conf_high,
      p_value,
      fdr,
      significant_fdr05
    ) |>
    dplyr::rename_with(
      ~ paste0(modality_label, "_", .x),
      -c("contrast", "gene_symbol")
    )
}

rna_by_symbol <- best_symbol_row(rna_results, "rna")
protein_by_symbol <- best_symbol_row(protein_results, "protein")

cross_modal <- dplyr::inner_join(
  rna_by_symbol,
  protein_by_symbol,
  by = c("contrast", "gene_symbol")
) |>
  dplyr::mutate(
    is_hsp60_10_client = .data$gene_symbol %in% hsp_clients,
    significance_class = dplyr::case_when(
      .data$rna_significant_fdr05 &
        .data$protein_significant_fdr05 ~ "FDR-significant in both",
      .data$rna_significant_fdr05 ~ "RNA only",
      .data$protein_significant_fdr05 ~ "Protein only",
      TRUE ~ "Neither"
    ),
    direction_concordant =
      sign(.data$rna_effect) == sign(.data$protein_effect)
  ) |>
  dplyr::arrange(
    factor(.data$contrast, levels = contrast_specs$contrast),
    dplyr::desc(.data$is_hsp60_10_client),
    .data$gene_symbol
  )

write_out(
  cross_modal,
  "08_cross_modal_shared_features_side_by_side_NO_SUBTRACTION.csv"
)

cross_modal_hsp <- cross_modal |>
  dplyr::filter(.data$is_hsp60_10_client)

write_out(
  cross_modal_hsp,
  "09_cross_modal_Hsp60_10_clients_side_by_side_NO_SUBTRACTION.csv"
)

cross_modal_hsp_summary <- cross_modal_hsp |>
  dplyr::group_by(.data$contrast) |>
  dplyr::summarise(
    n_shared_hsp_clients = dplyr::n(),
    n_rna_only_fdr = sum(.data$significance_class == "RNA only"),
    n_protein_only_fdr = sum(.data$significance_class == "Protein only"),
    n_both_fdr = sum(.data$significance_class == "FDR-significant in both"),
    n_neither_fdr = sum(.data$significance_class == "Neither"),
    n_direction_concordant = sum(
      .data$direction_concordant,
      na.rm = TRUE
    ),
    pct_direction_concordant = 100 * mean(
      .data$direction_concordant,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_out(
  cross_modal_hsp_summary,
  "10_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv"
)

# ------------------------------------------------------------
# 8. Stage/sample/model audit
# ------------------------------------------------------------

sample_audit <- dplyr::bind_rows(
  rna_meta_model |>
    dplyr::count(.data$clinical_stage, name = "n") |>
    dplyr::mutate(
      modality = "RNA",
      eligibility = "canonical stage + production RNA covariates complete"
    ),
  prot_meta_model |>
    dplyr::count(.data$clinical_stage, name = "n") |>
    dplyr::mutate(
      modality = "Protein",
      eligibility = "canonical stage + production protein covariates complete"
    )
) |>
  dplyr::select(
    .data$modality,
    .data$clinical_stage,
    .data$n,
    .data$eligibility
  )

write_out(
  sample_audit,
  "11_sample_stage_audit.csv"
)

model_audit <- tibble::tibble(
  modality = c("RNA", "Protein"),
  method = c("DESeq2", "limma"),
  abundance_input = c(
    "raw integer RNA-seq counts",
    "processed non-residualized log2 TMT abundance"
  ),
  formula = c(
    paste(deparse(rna_formula), collapse = ""),
    paste(deparse(protein_formula), collapse = "")
  ),
  production_covariates = c(
    paste(rna_covars, collapse = ";"),
    paste(protein_covars, collapse = ";")
  ),
  n_samples = c(
    nrow(rna_meta_model),
    nrow(prot_meta_model)
  ),
  n_features_before_or_available = c(
    rna_genes_before_filter,
    nrow(protein_mat_model)
  ),
  n_features_modeled = c(
    rna_genes_after_filter,
    nrow(protein_mat_model)
  )
)

write_out(
  model_audit,
  "12_model_specification_audit.csv"
)

input_audit <- tibble::tibble(
  item = c(
    "raw_counts_path",
    "raw_counts_orientation",
    "raw_counts_rows_original",
    "raw_counts_columns_original",
    "rna_samples_aligned",
    "hsp60_10_reference_clients"
  ),
  value = c(
    counts_path,
    orientation,
    as.character(nrow(rna_counts)),
    as.character(ncol(rna_counts)),
    as.character(nrow(rna_meta_model)),
    as.character(length(hsp_clients))
  )
)

write_out(
  input_audit,
  "13_input_audit.csv"
)

# ------------------------------------------------------------
# 9. Reviewer-facing plain-text report
# ------------------------------------------------------------

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, digits = digits, format = "g")
  )
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      formatC(x, format = "e", digits = 2),
      formatC(x, format = "f", digits = 4)
    )
  )
}

report_lines <- c(
  "CONVENTIONAL STAGE DIFFERENTIAL ANALYSIS — REVIEWER REPORT",
  "=========================================================",
  "",
  "Reviewer-facing purpose:",
  "  Perform conventional differential gene and protein expression analyses",
  "  across canonical NCI, MCI, and AD stages as an alternative to direct",
  "  RNA-protein effect subtraction.",
  "",
  paste0(
    "RNA method: DESeq2 on raw counts; design ",
    paste(deparse(rna_formula), collapse = "")
  ),
  paste0(
    "Protein method: limma on processed log2 TMT abundance; design ",
    paste(deparse(protein_formula), collapse = "")
  ),
  "",
  "No RNA-minus-protein or protein-minus-RNA effect was calculated.",
  "",
  "Sample counts:",
  paste0(
    "  RNA: NCI=", expected_rna_stage_counts["NCI"],
    ", MCI=", expected_rna_stage_counts["MCI"],
    ", AD=", expected_rna_stage_counts["AD"]
  ),
  paste0(
    "  Protein nuisance-complete: NCI=",
    expected_protein_stage_counts["NCI"],
    ", MCI=", expected_protein_stage_counts["MCI"],
    ", AD=", expected_protein_stage_counts["AD"]
  ),
  ""
)

for (contrast_i in contrast_specs$contrast) {
  report_lines <- c(
    report_lines,
    paste0("CONTRAST: ", contrast_i),
    paste(rep("-", nchar(contrast_i) + 10), collapse = "")
  )

  for (modality_i in c("RNA", "Protein")) {
    s <- contrast_summary |>
      dplyr::filter(
        .data$contrast == contrast_i,
        .data$modality == modality_i
      )

    report_lines <- c(
      report_lines,
      paste0(
        "  ", modality_i,
        ": ", s$n_fdr_significant,
        " FDR-significant features (",
        s$n_significant_higher, " higher; ",
        s$n_significant_lower, " lower in numerator stage)."
      ),
      paste0(
        "    Hsp60/10 clients: ",
        s$n_hsp_clients_fdr_significant,
        " FDR-significant of ",
        s$n_hsp_clients_tested,
        " tested (",
        s$n_hsp_clients_higher, " higher; ",
        s$n_hsp_clients_lower, " lower)."
      )
    )
  }

  core_i <- core_chaperonin_results |>
    dplyr::filter(.data$contrast == contrast_i)

  for (gene_i in c("HSPD1", "HSPE1")) {
    for (modality_i in c("RNA", "Protein")) {
      g <- core_i |>
        dplyr::filter(
          .data$gene_symbol == gene_i,
          .data$modality == modality_i
        )

      if (nrow(g) == 1) {
        report_lines <- c(
          report_lines,
          paste0(
            "    ", gene_i, " [", modality_i, "] effect=",
            fmt_num(g$effect),
            ", 95% CI ",
            fmt_num(g$conf_low), " to ", fmt_num(g$conf_high),
            ", P=", fmt_p(g$p_value),
            ", FDR=", fmt_p(g$fdr)
          )
        )
      }
    }
  }

  hsp_top_i <- top_hsp_hits |>
    dplyr::filter(.data$contrast == contrast_i) |>
    dplyr::group_by(.data$modality) |>
    dplyr::slice_head(n = 5) |>
    dplyr::ungroup()

  if (nrow(hsp_top_i) > 0) {
    report_lines <- c(
      report_lines,
      "  Top FDR-significant Hsp60/10 clients:"
    )

    for (j in seq_len(nrow(hsp_top_i))) {
      report_lines <- c(
        report_lines,
        paste0(
          "    ", hsp_top_i$modality[j], " — ",
          hsp_top_i$gene_symbol[j],
          ": effect=", fmt_num(hsp_top_i$effect[j]),
          ", FDR=", fmt_p(hsp_top_i$fdr[j])
        )
      )
    }
  } else {
    report_lines <- c(
      report_lines,
      "  No Hsp60/10 clients reached FDR < 0.05 in this contrast."
    )
  }

  cm <- cross_modal_hsp_summary |>
    dplyr::filter(.data$contrast == contrast_i)

  report_lines <- c(
    report_lines,
    paste0(
      "  Shared Hsp60/10 clients analyzed side-by-side: ",
      cm$n_shared_hsp_clients,
      "; RNA-only FDR=", cm$n_rna_only_fdr,
      "; protein-only FDR=", cm$n_protein_only_fdr,
      "; both FDR=", cm$n_both_fdr,
      "; neither=", cm$n_neither_fdr,
      "."
    ),
    paste0(
      "  Direction concordance across shared client effects: ",
      fmt_num(cm$pct_direction_concordant, digits = 4),
      "%."
    ),
    ""
  )
}

report_lines <- c(
  report_lines,
  "Interpretation guardrails:",
  "  - RNA and protein effects are reported separately.",
  "  - Cross-modal tables are side-by-side only; effects are not subtracted.",
  "  - FDR < 0.05 defines differential-expression/abundance significance.",
  "  - Direction concordance is descriptive and is not evidence of",
  "    post-transcriptional regulation.",
  ""
)

report_path <- file.path(out_dir, "14_REVIEWER_REPORT_SUMMARY.txt")
writeLines(report_lines, report_path)
message("Wrote: ", report_path)

# ------------------------------------------------------------
# 10. Validation table and provenance
# ------------------------------------------------------------

validation <- tibble::tribble(
  ~check, ~observed, ~expected, ~passed,
  "RNA canonical stage NCI", as.character(observed_rna_stage_counts["NCI"]), "200",
  observed_rna_stage_counts["NCI"] == 200,
  "RNA canonical stage MCI", as.character(observed_rna_stage_counts["MCI"]), "158",
  observed_rna_stage_counts["MCI"] == 158,
  "RNA canonical stage AD", as.character(observed_rna_stage_counts["AD"]), "220",
  observed_rna_stage_counts["AD"] == 220,
  "Protein nuisance-complete NCI", as.character(observed_protein_stage_counts["NCI"]), "167",
  observed_protein_stage_counts["NCI"] == 167,
  "Protein nuisance-complete MCI", as.character(observed_protein_stage_counts["MCI"]), "96",
  observed_protein_stage_counts["MCI"] == 96,
  "Protein nuisance-complete AD", as.character(observed_protein_stage_counts["AD"]), "109",
  observed_protein_stage_counts["AD"] == 109,
  "RNA contrasts", as.character(dplyr::n_distinct(rna_results$contrast)), "3",
  dplyr::n_distinct(rna_results$contrast) == 3,
  "Protein contrasts", as.character(dplyr::n_distinct(protein_results$contrast)), "3",
  dplyr::n_distinct(protein_results$contrast) == 3,
  "Hsp60/10 reference clients", as.character(length(hsp_clients)), "321",
  length(hsp_clients) == 321,
  "HSPD1/HSPE1 rows", as.character(nrow(core_chaperonin_results)), "12",
  nrow(core_chaperonin_results) == 12,
  "No cross-modal subtraction columns",
  as.character(
    !any(
      grepl(
        "minus|difference",
        colnames(cross_modal),
        ignore.case = TRUE
      )
    )
  ),
  "TRUE",
  !any(
    grepl(
      "minus|difference",
      colnames(cross_modal),
      ignore.case = TRUE
    )
  )
)

write_out(
  validation,
  "15_validation_checks.csv"
)

if (any(!validation$passed)) {
  stop(
    "Conventional differential-analysis validation failed. ",
    "See 15_validation_checks.csv.",
    call. = FALSE
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

git_head <- tryCatch(
  system2(
    "git",
    c("rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) NA_character_
)

git_status <- tryCatch(
  system2(
    "git",
    c("status", "--short"),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) NA_character_
)

writeLines(
  c(
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Git HEAD: ", paste(git_head, collapse = "")),
    "Git status --short:",
    git_status,
    "",
    paste0("Raw counts controlled input: ", counts_path),
    paste0("RNA formula: ", paste(deparse(rna_formula), collapse = "")),
    paste0("Protein formula: ", paste(deparse(protein_formula), collapse = "")),
    "RNA filter: count >= 10 in >= 3 samples",
    "Multiple testing: Benjamini-Hochberg FDR",
    "Significance threshold: FDR < 0.05",
    "Cross-modal policy: report side-by-side; do not subtract effects"
  ),
  file.path(out_dir, "run_provenance.txt")
)

message("\n============================================================")
message("CONVENTIONAL STAGE DIFFERENTIAL ANALYSIS COMPLETE")
message("============================================================")
print(sample_audit, n = Inf)
message("\nContrast summary:")
print(contrast_summary, n = Inf)
message("\nHSPD1 / HSPE1:")
print(
  core_chaperonin_results |>
    dplyr::select(
      .data$modality,
      .data$contrast,
      .data$gene_symbol,
      .data$effect,
      .data$std_error,
      .data$conf_low,
      .data$conf_high,
      .data$p_value,
      .data$fdr
    ),
  n = Inf
)
message("\nCross-modal Hsp60/10 side-by-side summary:")
print(cross_modal_hsp_summary, n = Inf)
message("\nReviewer report:")
message(report_path)
