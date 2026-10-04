############################################################
## Reviewer revision audit 05
## Compare primary-pure vs broad ROSMAP cogdx stage definitions
## on the correctly ingested full TMT cohort
##
## Rationale
## ---------
## ROSMAP cogdx codes distinguish:
##   1 = NCI
##   2 = MCI, no other contributing cause
##   3 = MCI + another contributing cause
##   4 = AD dementia, no other contributing cause
##   5 = AD dementia + another contributing cause
##   6 = other dementia
##
## The submitted RNA analysis effectively used the narrower 1/2/4
## subset, whereas the source ROSMAP proteomics publication describes
## clinical NCI/MCI/AD groups based on cogdx more broadly.
##
## This audit compares:
##
##   PRIMARY_PURE:
##     1 -> NCI
##     2 -> MCI
##     4 -> AD
##     3/5/6 excluded
##
##   BROAD_CONSENSUS:
##     1 -> NCI
##     2/3 -> MCI
##     4/5 -> AD
##     6 excluded
##
## It also audits agreement between final consensus cogdx and the
## last-valid clinical diagnosis variable (dcfdx_lv) when available.
##
## This script is audit-only. It does not modify production objects.
##
## Run from repository root:
##   Rscript R/reviewer_revisions/05_audit_cogdx_stage_definition_sensitivity.R
############################################################

options(stringsAsFactors = FALSE)

required_scripts <- file.path("R", c("00_config.R", "01_utils.R"))
missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "Run this script from the repository root. Missing:\n",
    paste(missing_scripts, collapse = "\n"),
    call. = FALSE
  )
}

for (script_i in required_scripts) source(script_i)

require_objects("cfg", context = "reviewer audit 05")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "cogdx_stage_definition_sensitivity"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_review_csv <- function(x, filename) {
  out <- file.path(out_dir, filename)
  readr::write_csv(x, out)
  message("Wrote: ", out)
  invisible(x)
}

canonical_batch_channel <- function(x) {
  x <- stringr::str_trim(as.character(x))
  hit <- stringr::str_match(
    x,
    stringr::regex("^b0*([0-9]+)\\.(.+)$", ignore_case = TRUE)
  )
  out <- x
  ok <- !is.na(hit[, 1])
  out[ok] <- paste0("b", as.integer(hit[ok, 2]), ".", hit[ok, 3])
  out
}

stage_primary_pure <- function(x) {
  x <- safe_num(x)
  dplyr::case_when(
    x == 1 ~ "NCI",
    x == 2 ~ "MCI",
    x == 4 ~ "AD",
    TRUE ~ NA_character_
  )
}

stage_broad_consensus <- function(x) {
  x <- safe_num(x)
  dplyr::case_when(
    x == 1 ~ "NCI",
    x %in% c(2, 3) ~ "MCI",
    x %in% c(4, 5) ~ "AD",
    TRUE ~ NA_character_
  )
}

is_present_chr <- function(x) {
  !is.na(x) & nzchar(stringr::str_trim(as.character(x)))
}

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

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

# ------------------------------------------------------------
# 1. Resolve source files
# ------------------------------------------------------------

analysis_meta_file <- first_existing(
  c(
    file.path(cfg$metadata_dir, "Analysis_Meta_Merged.csv"),
    file.path(cfg$metadata_dir, "Analysis_Meta_Merged_copy.csv"),
    file.path(cfg$derived_dir, "Analysis_Meta_Merged_copy.csv")
  ),
  "Analysis_Meta_Merged metadata"
)

rna_file <- first_existing(
  c(
    file.path(cfg$derived_dir, "ROSMAP_vst_gene_symbol_matrix_STAGE.csv"),
    file.path(
      cfg$project_dir,
      "ROSMAP_pathway_scores_robustness",
      "ROSMAP_vst_gene_symbol_matrix_STAGE.csv"
    )
  ),
  "ROSMAP VST gene-symbol RNA matrix"
)

protein_file <- first_existing(
  c(
    file.path(
      cfg$derived_dir,
      "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_copy.csv"
    ),
    file.path(
      cfg$proteomics_input_dir,
      "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"
    ),
    file.path(
      cfg$proteomics_input_dir,
      "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_copy.csv"
    )
  ),
  "full ROSMAP TMT protein matrix"
)

protein_meta_file <- first_existing(
  c(
    file.path(cfg$derived_dir, "matched_metadata_copy.csv"),
    file.path(cfg$proteomics_input_dir, "matched_metadata.csv"),
    file.path(cfg$proteomics_input_dir, "matched_metadata_copy.csv")
  ),
  "ROSMAP TMT matched metadata"
)

clinical_file <- cfg$rosmap_clinical_file
if (!file.exists(clinical_file)) {
  stop("ROSMAP clinical file not found: ", clinical_file, call. = FALSE)
}

write_review_csv(
  tibble::tibble(
    input = c(
      "Analysis_Meta_Merged",
      "RNA VST matrix",
      "Full TMT matrix",
      "TMT matched metadata",
      "ROSMAP clinical"
    ),
    path = c(
      analysis_meta_file,
      rna_file,
      protein_file,
      protein_meta_file,
      clinical_file
    )
  ),
  "01_input_source_paths.csv"
)

# ------------------------------------------------------------
# 2. Clinical metadata and stage definitions
# ------------------------------------------------------------

clinical_raw <- readr::read_csv(
  clinical_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (!"individualID" %in% colnames(clinical_raw)) {
  stop(
    "ROSMAP_clinical.csv must contain individualID for the proteomics join. ",
    "Available columns: ",
    paste(colnames(clinical_raw), collapse = ", "),
    call. = FALSE
  )
}

cogdx_col <- require_alias_col(
  clinical_raw,
  c("cogdx", "COGDX"),
  "cogdx",
  context = "ROSMAP clinical"
)

dcfdx_lv_col <- pick_existing(
  clinical_raw,
  c("dcfdx_lv", "dcfdx_l", "dcfdx")
)

age_col <- require_alias_col(
  clinical_raw,
  alias_sets$age_death,
  "age at death",
  context = "ROSMAP clinical"
)

sex_col <- require_alias_col(
  clinical_raw,
  alias_sets$sex,
  "sex",
  context = "ROSMAP clinical"
)

pmi_col <- require_alias_col(
  clinical_raw,
  alias_sets$pmi,
  "PMI",
  context = "ROSMAP clinical"
)

braak_col <- require_alias_col(
  clinical_raw,
  alias_sets$braak,
  "Braak",
  context = "ROSMAP clinical"
)

cerad_col <- require_alias_col(
  clinical_raw,
  alias_sets$cerad,
  "CERAD",
  context = "ROSMAP clinical"
)

clinical <- tibble::tibble(
  individual_id = as.character(clinical_raw$individualID),
  cogdx = safe_num(clinical_raw[[cogdx_col]]),
  dcfdx_lv = if (!is.na(dcfdx_lv_col)) {
    safe_num(clinical_raw[[dcfdx_lv_col]])
  } else {
    NA_real_
  },
  age_death = safe_num(clinical_raw[[age_col]]),
  sex = as.character(clinical_raw[[sex_col]]),
  pmi = safe_num(clinical_raw[[pmi_col]]),
  braak = safe_num(clinical_raw[[braak_col]]),
  cerad = safe_num(clinical_raw[[cerad_col]])
) |>
  dplyr::filter(is_present_chr(.data$individual_id)) |>
  dplyr::mutate(
    stage_primary_pure = stage_primary_pure(.data$cogdx),
    stage_broad_consensus = stage_broad_consensus(.data$cogdx),
    dcfdx_stage_broad = stage_broad_consensus(.data$dcfdx_lv)
  )

if (anyDuplicated(clinical$individual_id) > 0) {
  stop("ROSMAP clinical individualID is not unique.", call. = FALSE)
}

stage_definition_key <- tibble::tibble(
  cogdx = 1:6,
  codebook_interpretation = c(
    "NCI",
    "MCI; no other contributing cause",
    "MCI; another contributing cause",
    "AD dementia; no other contributing cause",
    "AD dementia; another contributing cause",
    "Other dementia"
  ),
  primary_pure_stage = stage_primary_pure(1:6),
  broad_consensus_stage = stage_broad_consensus(1:6)
)

write_review_csv(stage_definition_key, "02_stage_definition_key.csv")

# ------------------------------------------------------------
# 3. RNA sample counts under both definitions
# ------------------------------------------------------------

analysis_meta <- readr::read_csv(
  analysis_meta_file,
  show_col_types = FALSE
) |>
  janitor::clean_names() |>
  standardize_common_aliases(
    context = "Analysis_Meta_Merged metadata",
    require_ids = TRUE,
    require_diagnosis = TRUE
  ) |>
  dplyr::mutate(
    sample_id = as.character(.data$sample_id),
    individual_id = as.character(.data$individual_id),
    cogdx = safe_num(.data$diagnosis),
    stage_primary_pure = stage_primary_pure(.data$cogdx),
    stage_broad_consensus = stage_broad_consensus(.data$cogdx)
  )

rna_mat_full <- read_gene_matrix_csv(rna_file, gene_col = "gene_symbol")
rna_sample_ids <- colnames(rna_mat_full)

rna_stage_manifest <- analysis_meta |>
  dplyr::filter(.data$sample_id %in% rna_sample_ids) |>
  dplyr::distinct(.data$sample_id, .keep_all = TRUE) |>
  dplyr::select(
    "sample_id",
    "individual_id",
    "cogdx",
    "stage_primary_pure",
    "stage_broad_consensus"
  )

write_review_csv(rna_stage_manifest, "03_rna_stage_manifest.csv")

rna_stage_counts <- dplyr::bind_rows(
  rna_stage_manifest |>
    dplyr::filter(!is.na(.data$stage_primary_pure)) |>
    dplyr::count(stage = .data$stage_primary_pure, name = "n") |>
    dplyr::mutate(definition = "PRIMARY_PURE"),
  rna_stage_manifest |>
    dplyr::filter(!is.na(.data$stage_broad_consensus)) |>
    dplyr::count(stage = .data$stage_broad_consensus, name = "n") |>
    dplyr::mutate(definition = "BROAD_CONSENSUS")
) |>
  dplyr::select("definition", "stage", "n")

write_review_csv(rna_stage_counts, "04_rna_stage_counts_by_definition.csv")

# ------------------------------------------------------------
# 4. Correct full 400-sample protein ingestion
# ------------------------------------------------------------

protein_mat <- read_gene_matrix_csv(protein_file, gene_col = 1)

protein_meta <- readr::read_csv(
  protein_meta_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  janitor::clean_names()

names(protein_meta) <- stringr::str_replace_all(
  names(protein_meta),
  "\\.",
  "_"
)

protein_meta <- protein_meta |>
  add_alias_column(
    "batch_channel",
    c("batch_channel", "batch.channel", "batch", "channel", "tmt_channel"),
    context = "protein matched metadata"
  ) |>
  add_alias_column(
    "sample_id",
    alias_sets$sample_id,
    context = "protein matched metadata"
  ) |>
  add_alias_column(
    "individual_id",
    alias_sets$individual_id,
    context = "protein matched metadata"
  ) |>
  add_alias_column(
    "legacy_dx",
    c(
      "emory_strict_dx_2019",
      "emorystrictdx_2019",
      "EmoryStrictDx.2019",
      "diagnosis"
    ),
    context = "protein matched metadata"
  ) |>
  dplyr::mutate(
    batch_channel_raw = as.character(.data$batch_channel),
    batch_channel_canonical = canonical_batch_channel(.data$batch_channel),
    sample_id = as.character(.data$sample_id),
    individual_id = as.character(.data$individual_id),
    legacy_dx = as.character(.data$legacy_dx)
  )

if (anyDuplicated(protein_meta$batch_channel_canonical) > 0) {
  stop("Duplicate canonical batch/channel IDs in protein metadata.", call. = FALSE)
}

matrix_channel_raw <- colnames(protein_mat)
matrix_channel_canonical <- canonical_batch_channel(matrix_channel_raw)

if (anyDuplicated(matrix_channel_canonical) > 0) {
  stop("Duplicate canonical batch/channel IDs in protein matrix.", call. = FALSE)
}

idx <- match(matrix_channel_canonical, protein_meta$batch_channel_canonical)

if (any(is.na(idx))) {
  stop(
    "At least one protein matrix column does not map to matched metadata ",
    "after canonical batch/channel normalization.",
    call. = FALSE
  )
}

protein_meta_aligned <- protein_meta[idx, , drop = FALSE]

if (nrow(protein_meta_aligned) != ncol(protein_mat)) {
  stop("Protein metadata does not align one-to-one to matrix columns.", call. = FALSE)
}

# Preserve original matrix order, but use biological sample IDs downstream.
colnames(protein_mat) <- protein_meta_aligned$sample_id

# Choose a true batch term. Prefer an explicit batch column that is not the
# batch_channel identifier; otherwise derive the batch number from canonical
# batch/channel.
batch_col_candidates <- c(
  "batch",
  "batch_id",
  "study_batch",
  "tmt_batch",
  "plex",
  "plex_id"
)
batch_col <- pick_existing(protein_meta_aligned, batch_col_candidates)

batch_factor_raw <- if (!is.na(batch_col)) {
  as.character(protein_meta_aligned[[batch_col]])
} else {
  stringr::str_remove(
    protein_meta_aligned$batch_channel_canonical,
    "\\..*$"
  )
}

protein_manifest <- tibble::tibble(
  SampleID = protein_meta_aligned$sample_id,
  IndividualID = protein_meta_aligned$individual_id,
  batch_channel_raw = protein_meta_aligned$batch_channel_raw,
  batch_channel_canonical = protein_meta_aligned$batch_channel_canonical,
  batch_for_adjustment = batch_factor_raw,
  legacy_dx = protein_meta_aligned$legacy_dx
) |>
  dplyr::left_join(
    clinical |>
      dplyr::rename(IndividualID = .data$individual_id),
    by = "IndividualID"
  ) |>
  dplyr::mutate(
    stage_primary_pure = factor(
      .data$stage_primary_pure,
      levels = c("NCI", "MCI", "AD")
    ),
    stage_broad_consensus = factor(
      .data$stage_broad_consensus,
      levels = c("NCI", "MCI", "AD")
    ),
    sex_factor = factor(.data$sex),
    batch_factor = factor(.data$batch_for_adjustment)
  )

write_review_csv(protein_manifest, "05_protein_stage_manifest.csv")

protein_stage_counts <- dplyr::bind_rows(
  protein_manifest |>
    dplyr::filter(!is.na(.data$stage_primary_pure)) |>
    dplyr::count(stage = .data$stage_primary_pure, name = "n") |>
    dplyr::mutate(definition = "PRIMARY_PURE"),
  protein_manifest |>
    dplyr::filter(!is.na(.data$stage_broad_consensus)) |>
    dplyr::count(stage = .data$stage_broad_consensus, name = "n") |>
    dplyr::mutate(definition = "BROAD_CONSENSUS")
) |>
  dplyr::select("definition", "stage", "n")

write_review_csv(
  protein_stage_counts,
  "06_protein_stage_counts_by_definition.csv"
)

# ------------------------------------------------------------
# 5. Final consensus vs last-valid diagnosis audit
# ------------------------------------------------------------

diagnosis_agreement <- protein_manifest |>
  dplyr::transmute(
    SampleID,
    IndividualID,
    cogdx,
    dcfdx_lv,
    final_stage = as.character(.data$stage_broad_consensus),
    last_valid_stage = .data$dcfdx_stage_broad,
    both_available = !is.na(.data$final_stage) & !is.na(.data$last_valid_stage),
    stage_agrees = dplyr::if_else(
      .data$both_available,
      .data$final_stage == .data$last_valid_stage,
      NA
    )
  )

write_review_csv(
  diagnosis_agreement,
  "07_final_cogdx_vs_last_valid_diagnosis.csv"
)

diagnosis_agreement_summary <- tibble::tibble(
  metric = c(
    "Protein samples with final broad cogdx stage",
    "Protein samples with last-valid broad stage",
    "Protein samples comparable on both",
    "Broad stage agreements",
    "Broad stage disagreements"
  ),
  value = c(
    sum(!is.na(diagnosis_agreement$final_stage)),
    sum(!is.na(diagnosis_agreement$last_valid_stage)),
    sum(diagnosis_agreement$both_available),
    sum(diagnosis_agreement$stage_agrees %in% TRUE, na.rm = TRUE),
    sum(diagnosis_agreement$stage_agrees %in% FALSE, na.rm = TRUE)
  )
)

write_review_csv(
  diagnosis_agreement_summary,
  "08_final_vs_last_valid_diagnosis_summary.csv"
)

# ------------------------------------------------------------
# 6. Nuisance adjustment on the complete full protein cohort
# ------------------------------------------------------------

protein_adjust_meta <- protein_manifest |>
  dplyr::transmute(
    SampleID,
    age_num = .data$age_death,
    sex_factor = .data$sex_factor,
    pmi_num = .data$pmi,
    batch_factor = .data$batch_factor
  )

protein_covars <- select_usable_covars(
  protein_adjust_meta,
  c("age_num", "sex_factor", "pmi_num", "batch_factor")
)

if (!all(c("age_num", "sex_factor", "pmi_num") %in% protein_covars)) {
  stop(
    "Age, sex, and PMI were not all selected as usable protein covariates. ",
    "Selected: ",
    paste(protein_covars, collapse = ", "),
    call. = FALSE
  )
}

protein_mat_adj <- residualize_matrix(
  mat = protein_mat,
  meta_df = protein_adjust_meta,
  sample_col = "SampleID",
  covars = protein_covars
)

n_adjusted_samples <- sum(
  colSums(is.finite(protein_mat_adj)) > 0
)

adjustment_summary <- tibble::tibble(
  metric = c(
    "Raw protein matrix samples",
    "Samples with at least one finite adjusted protein",
    "Protein covariates used"
  ),
  value = c(
    as.character(ncol(protein_mat)),
    as.character(n_adjusted_samples),
    paste(protein_covars, collapse = ";")
  )
)

write_review_csv(
  adjustment_summary,
  "09_protein_adjustment_summary.csv"
)

# ------------------------------------------------------------
# 7. Hsp60/10 client inventory
# ------------------------------------------------------------

if (!file.exists(cfg$hsp60_client_file)) {
  stop(
    "Hsp60/10 client inventory not found: ",
    cfg$hsp60_client_file,
    call. = FALSE
  )
}

client_raw <- readxl::read_excel(
  cfg$hsp60_client_file,
  sheet = cfg$hsp60_client_sheet,
  skip = cfg$hsp60_client_skip
)

gene_candidates <- c(
  "Gene name", "gene", "Gene", "SYMBOL", "symbol", "gene_symbol"
)
client_gene_col <- gene_candidates[
  gene_candidates %in% colnames(client_raw)
][1]

if (is.na(client_gene_col)) {
  stop(
    "No recognizable Hsp60/10 client gene column. Available: ",
    paste(colnames(client_raw), collapse = ", "),
    call. = FALSE
  )
}

hsp_clients <- client_raw |>
  dplyr::transmute(
    gene = canonical_gene_symbol(.data[[client_gene_col]])
  ) |>
  dplyr::filter(
    !is.na(.data$gene),
    nzchar(.data$gene),
    !.data$gene %in% c("N/A", "NA", "HSPD1", "HSPE1", "HSPE1-MOB4")
  ) |>
  dplyr::distinct(.data$gene) |>
  dplyr::pull(.data$gene)

detected_clients <- intersect(hsp_clients, rownames(protein_mat_adj))

if (length(detected_clients) < 10) {
  stop("Too few Hsp60/10 clients detected in corrected protein matrix.", call. = FALSE)
}

# ------------------------------------------------------------
# 8. Pathway-score trajectories under both cogdx definitions
# ------------------------------------------------------------

client_score <- score_pathway_mean_z(
  protein_mat_adj,
  detected_clients,
  min_genes = 2
)

score_tbl <- tibble::tibble(
  SampleID = colnames(protein_mat_adj),
  Hsp60_10_client_score = as.numeric(client_score)
) |>
  dplyr::left_join(
    protein_manifest |>
      dplyr::select(
        "SampleID",
        "stage_primary_pure",
        "stage_broad_consensus"
      ),
    by = "SampleID"
  )

score_long <- dplyr::bind_rows(
  score_tbl |>
    dplyr::transmute(
      SampleID,
      definition = "PRIMARY_PURE",
      stage = as.character(.data$stage_primary_pure),
      score = .data$Hsp60_10_client_score
    ),
  score_tbl |>
    dplyr::transmute(
      SampleID,
      definition = "BROAD_CONSENSUS",
      stage = as.character(.data$stage_broad_consensus),
      score = .data$Hsp60_10_client_score
    )
) |>
  dplyr::filter(
    .data$stage %in% c("NCI", "MCI", "AD"),
    is.finite(.data$score)
  ) |>
  dplyr::mutate(
    stage = factor(.data$stage, levels = c("NCI", "MCI", "AD"))
  )

pathway_summary <- score_long |>
  dplyr::group_by(.data$definition, .data$stage) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_score = mean(.data$score),
    sem = safe_sem(.data$score),
    .groups = "drop"
  )

write_review_csv(
  pathway_summary,
  "10_hsp60_pathway_stage_summary.csv"
)

pathway_effects <- pathway_summary |>
  dplyr::select("definition", "stage", "mean_score") |>
  tidyr::pivot_wider(
    names_from = "stage",
    values_from = "mean_score"
  ) |>
  dplyr::mutate(
    early_MCI_minus_NCI = .data$MCI - .data$NCI,
    late_AD_minus_MCI = .data$AD - .data$MCI,
    total_AD_minus_NCI = .data$AD - .data$NCI
  )

write_review_csv(
  pathway_effects,
  "11_hsp60_pathway_transition_effects.csv"
)

# ------------------------------------------------------------
# 9. Client gene-level effect sensitivity
# ------------------------------------------------------------

make_stage_meta <- function(stage_col) {
  protein_manifest |>
    dplyr::transmute(
      SampleID,
      stage = .data[[stage_col]]
    )
}

gene_effects_pure <- fit_stage_effects(
  expr_mat = protein_mat_adj[detected_clients, , drop = FALSE],
  meta_df = make_stage_meta("stage_primary_pure"),
  sample_col = "SampleID",
  stage_col = "stage",
  stage_levels = c("NCI", "MCI", "AD"),
  modality_name = "pure",
  min_n = 10
)

gene_effects_broad <- fit_stage_effects(
  expr_mat = protein_mat_adj[detected_clients, , drop = FALSE],
  meta_df = make_stage_meta("stage_broad_consensus"),
  sample_col = "SampleID",
  stage_col = "stage",
  stage_levels = c("NCI", "MCI", "AD"),
  modality_name = "broad",
  min_n = 10
)

client_effect_comparison <- dplyr::inner_join(
  gene_effects_pure,
  gene_effects_broad,
  by = "gene"
) |>
  dplyr::mutate(
    late_decline_pure = dplyr::case_when(
      !is.finite(.data$late_shift__pure) ~ NA_real_,
      .data$late_shift__pure < 0 ~ abs(.data$late_shift__pure),
      TRUE ~ 0
    ),
    late_decline_broad = dplyr::case_when(
      !is.finite(.data$late_shift__broad) ~ NA_real_,
      .data$late_shift__broad < 0 ~ abs(.data$late_shift__broad),
      TRUE ~ 0
    )
  )

write_review_csv(
  client_effect_comparison,
  "12_hsp60_client_gene_effects_pure_vs_broad.csv"
)

effect_concordance <- tibble::tibble(
  metric = c(
    "Detected Hsp60/10 clients",
    "Spearman rho early shift: pure vs broad",
    "Spearman rho late shift: pure vs broad",
    "Spearman rho total shift: pure vs broad",
    "Spearman rho late decline magnitude: pure vs broad",
    "Fraction negative late shift: pure",
    "Fraction negative late shift: broad",
    "Mean late decline magnitude: pure",
    "Mean late decline magnitude: broad"
  ),
  value = c(
    length(detected_clients),
    safe_cor(
      client_effect_comparison$early_shift__pure,
      client_effect_comparison$early_shift__broad
    ),
    safe_cor(
      client_effect_comparison$late_shift__pure,
      client_effect_comparison$late_shift__broad
    ),
    safe_cor(
      client_effect_comparison$total_shift__pure,
      client_effect_comparison$total_shift__broad
    ),
    safe_cor(
      client_effect_comparison$late_decline_pure,
      client_effect_comparison$late_decline_broad
    ),
    mean(
      client_effect_comparison$late_shift__pure < 0,
      na.rm = TRUE
    ),
    mean(
      client_effect_comparison$late_shift__broad < 0,
      na.rm = TRUE
    ),
    mean(
      client_effect_comparison$late_decline_pure,
      na.rm = TRUE
    ),
    mean(
      client_effect_comparison$late_decline_broad,
      na.rm = TRUE
    )
  )
)

write_review_csv(
  effect_concordance,
  "00_stage_definition_effect_sensitivity_summary.csv"
)

# ------------------------------------------------------------
# 10. Provenance
# ------------------------------------------------------------

git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

git_branch <- tryCatch(
  system2("git", c("branch", "--show-current"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

provenance <- c(
  "Reviewer revision audit 05 provenance",
  paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0(
    "Repository root: ",
    normalizePath(cfg$project_dir, mustWork = FALSE)
  ),
  paste0("Git branch: ", paste(git_branch, collapse = " ")),
  paste0("Git HEAD: ", paste(git_head, collapse = " ")),
  "",
  "Stage definitions compared:",
  "  PRIMARY_PURE: 1=NCI, 2=MCI, 4=AD; 3/5/6 excluded",
  "  BROAD_CONSENSUS: 1=NCI, 2/3=MCI, 4/5=AD; 6 excluded",
  "",
  paste0("Detected Hsp60/10 clients: ", length(detected_clients)),
  paste0("Protein covariates used: ", paste(protein_covars, collapse = ";")),
  paste0("Adjusted protein samples with finite data: ", n_adjusted_samples),
  paste0(
    "Last-valid diagnosis column used: ",
    ifelse(is.na(dcfdx_lv_col), "NONE", dcfdx_lv_col)
  ),
  "",
  "This script is audit-only and does not modify production objects."
)

writeLines(
  provenance,
  file.path(out_dir, "13_run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "14_sessionInfo.txt")
)

# ------------------------------------------------------------
# 11. Terminal summary
# ------------------------------------------------------------

message("\n============================================================")
message("COGDX STAGE-DEFINITION SENSITIVITY AUDIT")
message("============================================================\n")

message("RNA STAGE COUNTS")
print(rna_stage_counts, n = Inf)

message("\nPROTEIN STAGE COUNTS")
print(protein_stage_counts, n = Inf)

message("\nFINAL COGDX VS LAST-VALID DIAGNOSIS")
print(diagnosis_agreement_summary, n = Inf)

message("\nHsp60/10 PATHWAY TRANSITION EFFECTS")
print(pathway_effects, n = Inf)

message("\nGENE-LEVEL EFFECT CONCORDANCE")
print(effect_concordance, n = Inf)

message("\nPrimary files to inspect:")
message(file.path(out_dir, "00_stage_definition_effect_sensitivity_summary.csv"))
message(file.path(out_dir, "04_rna_stage_counts_by_definition.csv"))
message(file.path(out_dir, "06_protein_stage_counts_by_definition.csv"))
message(file.path(out_dir, "08_final_vs_last_valid_diagnosis_summary.csv"))
message(file.path(out_dir, "11_hsp60_pathway_transition_effects.csv"))
message(file.path(out_dir, "12_hsp60_client_gene_effects_pure_vs_broad.csv"))
message("\nAudit complete.")
