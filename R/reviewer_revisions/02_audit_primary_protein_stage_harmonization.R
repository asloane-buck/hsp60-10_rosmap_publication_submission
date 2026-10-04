############################################################
## Reviewer revision audit 02
## Primary proteomics stage harmonization impact audit
##
## Purpose:
## 1. Reassign ROSMAP proteomics samples using the same participant-level
##    cogdx-derived clinical stage used for RNA (NCI/MCI/AD).
## 2. Compare those assignments with the legacy proteomics
##    Control/AsymAD/AD grouping used in the submitted analysis.
## 3. Quantify the impact on pathway stage-transition effects and
##    gene-level late-stage protein decline.
## 4. Write explicit manifests, effect tables, and run provenance.
##
## Reproducibility:
## - Run from repository root.
## - Uses only repository-relative paths through R/00_config.R.
## - Uses the canonical loaders and preprocessing functions from the
##   manuscript pipeline.
## - Does not modify raw data or existing manuscript figures.
## - No stochastic procedures are used.
##
## Run:
##   Rscript R/reviewer_revisions/02_audit_primary_protein_stage_harmonization.R
############################################################

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Repository/root checks and canonical pipeline objects
# ------------------------------------------------------------

required_core_scripts <- file.path(
  "R",
  c(
    "00_config.R",
    "01_utils.R",
    "02_load_data.R",
    "03_build_adjusted_core_objects.R",
    "04_build_pathway_sets_all_clients.R"
  )
)

missing_core_scripts <- required_core_scripts[!file.exists(required_core_scripts)]
if (length(missing_core_scripts) > 0) {
  stop(
    "Run this script from the repository root. Missing required script(s):\n",
    paste(missing_core_scripts, collapse = "\n"),
    call. = FALSE
  )
}

for (script_i in required_core_scripts) {
  source(script_i)
}

require_objects(
  c(
    "cfg",
    "analysis_meta",
    "prot_meta_adj",
    "prot_mat_adj",
    "pathway_gene_sets",
    "all_hsp60_10_clients"
  ),
  context = "reviewer audit 02"
)

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "primary_protein_stage_harmonization"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_review_csv <- function(x, filename) {
  out <- file.path(out_dir, filename)
  readr::write_csv(x, out)
  message("Wrote: ", out)
  invisible(x)
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- sort(unique(x))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

harmonize_common_stage <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("Control", "NCI") ~ "NCI",
    x %in% c("Early_AD", "MCI", "Intermediate") ~ "MCI",
    x == "AD" ~ "AD",
    TRUE ~ NA_character_
  )
}

harmonize_legacy_protein_stage <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "Control" ~ "NCI",
    x == "AsymAD" ~ "MCI",
    x == "AD" ~ "AD",
    TRUE ~ NA_character_
  )
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

safe_pearson <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = "pearson"))
}

# ------------------------------------------------------------
# 1. Build one canonical cogdx-derived stage per participant
# ------------------------------------------------------------

require_columns(
  analysis_meta,
  c("individual_id", "diagnosis", "diagnosis_stage"),
  "analysis_meta"
)

common_stage_raw <- analysis_meta |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    diagnosis_raw = as.character(.data$diagnosis),
    diagnosis_stage_raw = as.character(.data$diagnosis_stage),
    common_stage = harmonize_common_stage(.data$diagnosis_stage)
  ) |>
  dplyr::filter(
    !is.na(.data$individual_id),
    nzchar(.data$individual_id)
  )

stage_conflicts <- common_stage_raw |>
  dplyr::filter(!is.na(.data$common_stage)) |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    n_common_stages = dplyr::n_distinct(.data$common_stage),
    common_stage_values = collapse_unique(.data$common_stage),
    .groups = "drop"
  ) |>
  dplyr::filter(.data$n_common_stages > 1)

write_review_csv(stage_conflicts, "01_common_stage_conflicts.csv")

if (nrow(stage_conflicts) > 0) {
  stop(
    "Participant-level cogdx stage conflicts detected. ",
    "Inspect 01_common_stage_conflicts.csv before continuing.",
    call. = FALSE
  )
}

common_stage_by_id <- common_stage_raw |>
  dplyr::group_by(.data$individual_id) |>
  dplyr::summarise(
    diagnosis_raw = collapse_unique(.data$diagnosis_raw),
    diagnosis_stage_raw = collapse_unique(.data$diagnosis_stage_raw),
    common_stage = {
      vals <- unique(.data$common_stage[!is.na(.data$common_stage)])
      if (length(vals) == 0) NA_character_ else vals[[1]]
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    common_stage_status = dplyr::case_when(
      !is.na(.data$common_stage) ~ "Primary cogdx stage (1/2/4)",
      !is.na(.data$diagnosis_raw) ~ "Excluded non-primary cogdx",
      TRUE ~ "Missing cogdx"
    )
  )

# ------------------------------------------------------------
# 2. Build an explicit matrix-backed proteomics sample manifest
# ------------------------------------------------------------

require_columns(
  prot_meta_adj,
  c("SampleID", "IndividualID", "EmoryStrictDx.2019"),
  "prot_meta_adj"
)

protein_manifest <- prot_meta_adj |>
  dplyr::transmute(
    protein_sample_id = as.character(.data$SampleID),
    individual_id = as.character(.data$IndividualID),
    protein_legacy_dx = as.character(.data$EmoryStrictDx.2019),
    legacy_stage = harmonize_legacy_protein_stage(.data$EmoryStrictDx.2019)
  ) |>
  dplyr::filter(.data$protein_sample_id %in% colnames(prot_mat_adj)) |>
  dplyr::left_join(common_stage_by_id, by = "individual_id") |>
  dplyr::arrange(.data$protein_sample_id)

if (nrow(protein_manifest) != ncol(prot_mat_adj)) {
  stop(
    "Proteomics manifest does not contain exactly one row per matrix column: ",
    nrow(protein_manifest), " manifest rows versus ",
    ncol(prot_mat_adj), " matrix columns.",
    call. = FALSE
  )
}

if (anyDuplicated(protein_manifest$protein_sample_id) > 0) {
  stop("Duplicate protein SampleID values detected in matrix-backed manifest.", call. = FALSE)
}

if (anyDuplicated(protein_manifest$individual_id) > 0) {
  stop(
    "More than one matrix-backed proteomics sample exists for at least one participant. ",
    "Resolve repeated participants before stage-level analysis.",
    call. = FALSE
  )
}

if (!setequal(protein_manifest$protein_sample_id, colnames(prot_mat_adj))) {
  stop("Proteomics manifest SampleID set does not match the adjusted protein matrix.", call. = FALSE)
}

protein_manifest <- protein_manifest |>
  dplyr::mutate(
    common_stage = factor(.data$common_stage, levels = c("NCI", "MCI", "AD")),
    legacy_stage = factor(.data$legacy_stage, levels = c("NCI", "MCI", "AD")),
    legacy_matches_common = dplyr::if_else(
      !is.na(.data$common_stage) & !is.na(.data$legacy_stage),
      as.character(.data$common_stage) == as.character(.data$legacy_stage),
      NA
    )
  )

write_review_csv(protein_manifest, "02_protein_sample_stage_manifest.csv")

# ------------------------------------------------------------
# 3. Counts and stage crosswalk
# ------------------------------------------------------------

legacy_counts <- protein_manifest |>
  dplyr::count(stage = .data$legacy_stage, name = "n") |>
  dplyr::mutate(definition = "Legacy protein label: Control/AsymAD/AD")

common_counts <- protein_manifest |>
  dplyr::filter(!is.na(.data$common_stage)) |>
  dplyr::count(stage = .data$common_stage, name = "n") |>
  dplyr::mutate(definition = "Common clinical stage: cogdx 1/2/4")

excluded_counts <- protein_manifest |>
  dplyr::filter(is.na(.data$common_stage)) |>
  dplyr::count(
    stage = .data$common_stage_status,
    name = "n"
  ) |>
  dplyr::mutate(definition = "Not included in common cogdx 1/2/4 stage analysis")

stage_counts <- dplyr::bind_rows(
  legacy_counts |>
    dplyr::mutate(stage = as.character(.data$stage)),
  common_counts |>
    dplyr::mutate(stage = as.character(.data$stage)),
  excluded_counts |>
    dplyr::mutate(stage = as.character(.data$stage))
)

write_review_csv(stage_counts, "03_primary_protein_stage_counts.csv")

stage_crosswalk <- protein_manifest |>
  dplyr::mutate(
    common_stage_display = dplyr::if_else(
      is.na(.data$common_stage),
      .data$common_stage_status,
      as.character(.data$common_stage)
    )
  ) |>
  dplyr::count(
    .data$common_stage_display,
    .data$protein_legacy_dx,
    .data$legacy_stage,
    name = "n"
  ) |>
  dplyr::arrange(.data$common_stage_display, .data$protein_legacy_dx)

write_review_csv(stage_crosswalk, "04_common_vs_legacy_protein_stage_crosswalk.csv")

primary_stage_manifest <- protein_manifest |>
  dplyr::filter(!is.na(.data$common_stage)) |>
  dplyr::select(
    .data$protein_sample_id,
    .data$individual_id,
    .data$diagnosis_raw,
    .data$common_stage,
    .data$protein_legacy_dx,
    .data$legacy_stage,
    .data$legacy_matches_common
  )

write_review_csv(primary_stage_manifest, "05_common_stage_primary_protein_manifest.csv")

# ------------------------------------------------------------
# 4. Pathway abundance-score impact
# ------------------------------------------------------------

required_pathways <- c(
  "Hsp60_10_all_clients",
  "Broad_MitoCarta_non_Hsp60_10",
  "TCA_pyruvate_metabolism_non_Hsp60_10"
)

missing_pathways <- setdiff(required_pathways, names(pathway_gene_sets))
if (length(missing_pathways) > 0) {
  stop(
    "Required pathway set(s) missing: ",
    paste(missing_pathways, collapse = ", "),
    call. = FALSE
  )
}

too_small <- required_pathways[
  vapply(pathway_gene_sets[required_pathways], length, integer(1)) < 2
]
if (length(too_small) > 0) {
  stop(
    "Required pathway set(s) contain fewer than two genes: ",
    paste(too_small, collapse = ", "),
    ". Check controlled/public annotation inputs.",
    call. = FALSE
  )
}

protein_pathway_scores <- compute_pathway_scores(
  prot_mat_adj,
  pathway_gene_sets,
  sample_col = "protein_sample_id",
  method = "mean_z"
) |>
  dplyr::left_join(
    protein_manifest |>
      dplyr::select(
        .data$protein_sample_id,
        .data$legacy_stage,
        .data$common_stage
      ),
    by = "protein_sample_id"
  )

pathway_long <- protein_pathway_scores |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(pathway_gene_sets)),
    names_to = "pathway",
    values_to = "score"
  )

pathway_by_definition <- dplyr::bind_rows(
  pathway_long |>
    dplyr::transmute(
      protein_sample_id = .data$protein_sample_id,
      pathway = .data$pathway,
      score = .data$score,
      stage_definition = "legacy_Control_AsymAD_AD",
      stage = as.character(.data$legacy_stage)
    ),
  pathway_long |>
    dplyr::transmute(
      protein_sample_id = .data$protein_sample_id,
      pathway = .data$pathway,
      score = .data$score,
      stage_definition = "common_cogdx_1_2_4",
      stage = as.character(.data$common_stage)
    )
) |>
  dplyr::filter(
    .data$stage %in% c("NCI", "MCI", "AD"),
    is.finite(.data$score)
  ) |>
  dplyr::mutate(
    stage = factor(.data$stage, levels = c("NCI", "MCI", "AD"))
  )

pathway_group_summary <- pathway_by_definition |>
  dplyr::group_by(
    .data$stage_definition,
    .data$pathway,
    .data$stage
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_score = mean(.data$score),
    sd = stats::sd(.data$score),
    sem = sem(.data$score),
    .groups = "drop"
  )

write_review_csv(pathway_group_summary, "06_pathway_stage_group_summary.csv")

pathway_transition_effects <- pathway_group_summary |>
  dplyr::select(
    .data$stage_definition,
    .data$pathway,
    .data$stage,
    .data$mean_score
  ) |>
  tidyr::pivot_wider(
    names_from = .data$stage,
    values_from = .data$mean_score
  ) |>
  dplyr::mutate(
    early_shift_MCI_minus_NCI = .data$MCI - .data$NCI,
    late_shift_AD_minus_MCI = .data$AD - .data$MCI,
    total_shift_AD_minus_NCI = .data$AD - .data$NCI
  )

write_review_csv(
  pathway_transition_effects,
  "07_pathway_transition_effects_by_stage_definition.csv"
)

pathway_effect_comparison <- pathway_transition_effects |>
  dplyr::select(
    .data$stage_definition,
    .data$pathway,
    .data$early_shift_MCI_minus_NCI,
    .data$late_shift_AD_minus_MCI,
    .data$total_shift_AD_minus_NCI
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(c(
      "early_shift_MCI_minus_NCI",
      "late_shift_AD_minus_MCI",
      "total_shift_AD_minus_NCI"
    )),
    names_to = "transition",
    values_to = "effect"
  ) |>
  tidyr::pivot_wider(
    names_from = .data$stage_definition,
    values_from = .data$effect
  ) |>
  dplyr::mutate(
    common_minus_legacy = .data$common_cogdx_1_2_4 -
      .data$legacy_Control_AsymAD_AD
  )

write_review_csv(
  pathway_effect_comparison,
  "08_pathway_effect_common_vs_legacy.csv"
)

# ------------------------------------------------------------
# 5. Gene-level stage effects and late-stage decline impact
# ------------------------------------------------------------

legacy_meta <- protein_manifest |>
  dplyr::transmute(
    SampleID = .data$protein_sample_id,
    stage = .data$legacy_stage
  )

common_meta <- protein_manifest |>
  dplyr::transmute(
    SampleID = .data$protein_sample_id,
    stage = .data$common_stage
  )

legacy_gene_effects <- fit_stage_effects(
  expr_mat = prot_mat_adj,
  meta_df = legacy_meta,
  sample_col = "SampleID",
  stage_col = "stage",
  stage_levels = c("NCI", "MCI", "AD"),
  modality_name = "legacy",
  min_n = cfg$min_n_stage_model
)

common_gene_effects <- fit_stage_effects(
  expr_mat = prot_mat_adj,
  meta_df = common_meta,
  sample_col = "SampleID",
  stage_col = "stage",
  stage_levels = c("NCI", "MCI", "AD"),
  modality_name = "common",
  min_n = cfg$min_n_stage_model
)

gene_effect_comparison <- dplyr::full_join(
  legacy_gene_effects,
  common_gene_effects,
  by = "gene"
) |>
  dplyr::mutate(
    is_hsp60_10_client = .data$gene %in% all_hsp60_10_clients,
    late_decline_magnitude__legacy = dplyr::case_when(
      !is.finite(.data$late_shift__legacy) ~ NA_real_,
      .data$late_shift__legacy < 0 ~ abs(.data$late_shift__legacy),
      TRUE ~ 0
    ),
    late_decline_magnitude__common = dplyr::case_when(
      !is.finite(.data$late_shift__common) ~ NA_real_,
      .data$late_shift__common < 0 ~ abs(.data$late_shift__common),
      TRUE ~ 0
    ),
    late_shift_difference_common_minus_legacy =
      .data$late_shift__common - .data$late_shift__legacy,
    late_decline_difference_common_minus_legacy =
      .data$late_decline_magnitude__common - .data$late_decline_magnitude__legacy
  )

write_review_csv(
  gene_effect_comparison,
  "09_gene_stage_effects_common_vs_legacy.csv"
)

client_effect_comparison <- gene_effect_comparison |>
  dplyr::filter(.data$is_hsp60_10_client)

write_review_csv(
  client_effect_comparison,
  "10_hsp60_client_stage_effects_common_vs_legacy.csv"
)

# ------------------------------------------------------------
# 6. Concordance summaries
# ------------------------------------------------------------

stage_comparable <- protein_manifest |>
  dplyr::filter(
    !is.na(.data$legacy_stage),
    !is.na(.data$common_stage)
  )

n_stage_comparable <- nrow(stage_comparable)
n_stage_concordant <- sum(stage_comparable$legacy_matches_common)
n_stage_discordant <- n_stage_comparable - n_stage_concordant

client_q75_legacy <- stats::quantile(
  client_effect_comparison$late_decline_magnitude__legacy,
  probs = 0.75,
  na.rm = TRUE,
  names = FALSE
)

client_q75_common <- stats::quantile(
  client_effect_comparison$late_decline_magnitude__common,
  probs = 0.75,
  na.rm = TRUE,
  names = FALSE
)

top_legacy <- client_effect_comparison |>
  dplyr::filter(.data$late_decline_magnitude__legacy >= client_q75_legacy) |>
  dplyr::pull(.data$gene) |>
  unique()

top_common <- client_effect_comparison |>
  dplyr::filter(.data$late_decline_magnitude__common >= client_q75_common) |>
  dplyr::pull(.data$gene) |>
  unique()

top_intersection <- intersect(top_legacy, top_common)
top_union <- union(top_legacy, top_common)

concordance_summary <- tibble::tibble(
  metric = c(
    "Protein matrix samples",
    "Protein participants",
    "Protein samples with primary common cogdx stage",
    "Protein samples excluded from common cogdx 1/2/4 stage analysis",
    "Samples comparable under both stage definitions",
    "Stage-concordant samples",
    "Stage-discordant samples",
    "Stage discordance fraction",
    "All-protein Spearman rho: late shift legacy vs common",
    "All-protein Pearson r: late shift legacy vs common",
    "Hsp60/10-client Spearman rho: late shift legacy vs common",
    "Hsp60/10-client Pearson r: late shift legacy vs common",
    "Hsp60/10-client Spearman rho: decline magnitude legacy vs common",
    "Legacy client top-quartile threshold",
    "Common-stage client top-quartile threshold",
    "Legacy client top-quartile set size",
    "Common-stage client top-quartile set size",
    "Top-quartile client intersection size",
    "Top-quartile client union size",
    "Top-quartile client Jaccard"
  ),
  value = c(
    ncol(prot_mat_adj),
    dplyr::n_distinct(protein_manifest$individual_id),
    sum(!is.na(protein_manifest$common_stage)),
    sum(is.na(protein_manifest$common_stage)),
    n_stage_comparable,
    n_stage_concordant,
    n_stage_discordant,
    if (n_stage_comparable > 0) n_stage_discordant / n_stage_comparable else NA_real_,
    safe_spearman(
      gene_effect_comparison$late_shift__legacy,
      gene_effect_comparison$late_shift__common
    ),
    safe_pearson(
      gene_effect_comparison$late_shift__legacy,
      gene_effect_comparison$late_shift__common
    ),
    safe_spearman(
      client_effect_comparison$late_shift__legacy,
      client_effect_comparison$late_shift__common
    ),
    safe_pearson(
      client_effect_comparison$late_shift__legacy,
      client_effect_comparison$late_shift__common
    ),
    safe_spearman(
      client_effect_comparison$late_decline_magnitude__legacy,
      client_effect_comparison$late_decline_magnitude__common
    ),
    client_q75_legacy,
    client_q75_common,
    length(top_legacy),
    length(top_common),
    length(top_intersection),
    length(top_union),
    if (length(top_union) > 0) length(top_intersection) / length(top_union) else NA_real_
  )
)

write_review_csv(concordance_summary, "00_primary_stage_harmonization_summary.csv")

# Membership table for the two client top-quartile sets.
top_quartile_overlap <- tibble::tibble(gene = sort(top_union)) |>
  dplyr::mutate(
    top_quartile_legacy = .data$gene %in% top_legacy,
    top_quartile_common = .data$gene %in% top_common,
    in_both = .data$gene %in% top_intersection
  )

write_review_csv(
  top_quartile_overlap,
  "11_hsp60_client_top_quartile_overlap.csv"
)

# ------------------------------------------------------------
# 7. Run provenance
# ------------------------------------------------------------

git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

git_branch <- tryCatch(
  system2("git", c("branch", "--show-current"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

git_status <- tryCatch(
  system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
  error = function(e) "git status unavailable"
)

provenance_lines <- c(
  "Reviewer revision audit 02 provenance",
  paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Repository root: ", normalizePath(cfg$project_dir, mustWork = FALSE)),
  paste0("Git branch: ", paste(git_branch, collapse = " ")),
  paste0("Git HEAD: ", paste(git_head, collapse = " ")),
  paste0("Working tree dirty: ", length(git_status) > 0),
  "",
  "Canonical scripts sourced:",
  paste0("  ", required_core_scripts),
  "",
  "Key analysis settings:",
  paste0("  adjust_main_figures = ", cfg$adjust_main_figures),
  paste0("  min_n_stage_model = ", cfg$min_n_stage_model),
  "",
  "Git status --porcelain:",
  if (length(git_status) == 0) "  <clean>" else paste0("  ", git_status)
)

writeLines(
  provenance_lines,
  con = file.path(out_dir, "12_run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "13_sessionInfo.txt")
)

# ------------------------------------------------------------
# 8. Terminal summary
# ------------------------------------------------------------

message("\n============================================================")
message("PRIMARY PROTEOMICS STAGE HARMONIZATION IMPACT AUDIT")
message("============================================================\n")

message("STAGE COUNTS")
print(stage_counts, n = Inf)

message("\nCOMMON cogdx STAGE VS LEGACY PROTEIN LABEL")
print(stage_crosswalk, n = Inf)

message("\nCONCORDANCE / IMPACT SUMMARY")
print(concordance_summary, n = Inf)

message("\nPrimary files to inspect next:")
message(file.path(out_dir, "00_primary_stage_harmonization_summary.csv"))
message(file.path(out_dir, "03_primary_protein_stage_counts.csv"))
message(file.path(out_dir, "04_common_vs_legacy_protein_stage_crosswalk.csv"))
message(file.path(out_dir, "08_pathway_effect_common_vs_legacy.csv"))
message(file.path(out_dir, "10_hsp60_client_stage_effects_common_vs_legacy.csv"))
message("\nAudit complete.")
