# Run from the repository root:
# Rscript R/supplemental/91_export_supplemental_figures_manuscript_ready.R
# Replot frozen estimates, matching the revised supplement; do not refit analyses.
options(stringsAsFactors = FALSE)
if (!file.exists(file.path("R", "supplemental", "90_run_supplemental_pipeline.R"))) {
  stop("Run this command from the production repository root.", call. = FALSE)
}
SUPP_USE_FROZEN_MODELS <- TRUE
set.seed(1300)
root <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(HSP60_ROSMAP_PROJECT_DIR = root)
base <- file.path(root, "outputs", "supplemental_figures", "audits")
revision <- file.path(root, "outputs", "reviewer_revisions")
# Preflight every numerical source used by the new and corrected figures,
# before any figure script can overwrite an old export.
required_sources <- c(
  file.path(root, "outputs", "main_figures", "tables", "main_fig1_all_clients", "Fig1A_all_clients_overlay_summary.csv"),
  file.path(root, "outputs", "main_figures", "tables", "figure5_null_input_COVARIATE_ADJUSTED.csv"),
  file.path(root, "outputs", "main_figures", "tables", "main_fig5_matched_null_specificity", "Fig5_gene_universe_cognition_summary_by_gene.csv"),
  file.path(base, c("SuppFig3_matched_null_permutation_values.csv",
    "SuppFig4_plotted_gene_level_values.csv", "SuppFig5_cross_cohort_gene_values.csv")),
  file.path(revision, "regional_formal_interaction", "validated_final", c(
    "FINAL_region_x_AD_interaction_Hsp60_clients.csv",
    "FINAL_region_x_Braak_continuous_interaction_Hsp60_clients.csv", "FINAL_regional_interaction_summary.csv")),
  file.path(revision, "regional_formal_interaction", "hsp_conclusion_audit", "Hsp_pathway_score_primary_interactions.csv"),
  file.path(revision, "conventional_stage_differential_final", c(
    "17_PRIMARY_cross_modal_Hsp60_10_side_by_side_NO_SUBTRACTION.csv",
    "18_PRIMARY_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv")),
  file.path(revision, "Hsp_pathway_PC1_sensitivity", c(
    "participant_mean_z_and_PC1_scores.csv", "stage_score_summary.csv",
    "stage_Wilcoxon_mean_z_vs_PC1.csv", "mean_z_vs_PC1_concordance.csv")),
  file.path(revision, "variancePartition_source_of_variation", "MASTER_variance_fraction_summary.csv"),
  file.path(revision, "alternative_mechanism_marker_panel", c(
    "marker_stage_limma_results.csv", "marker_category_stage_summary.csv", "undetected_prespecified_markers.csv"))
)
missing_sources <- required_sources[!file.exists(required_sources)]
if (length(missing_sources)) {
  stop("Frozen supplemental sources are missing; nothing has been regenerated:\n",
    paste0(" - ", missing_sources, collapse = "\n"),
    "\nDo not rerun the statistical analyses just to bypass this check. Locate the frozen files first.", call. = FALSE)
}
source_md5_before <- tools::md5sum(required_sources)
# Models and frozen matched null draws remain read-only in this export mode.
protected <- file.path(base, c("SuppFig3_matched_null_permutation_values.csv",
  "SuppFig4_plotted_gene_level_values.csv", "SuppFig5_cross_cohort_gene_values.csv"))
source(file.path(root, "R", "supplemental", "90_run_supplemental_pipeline.R"))
source_md5_after <- tools::md5sum(required_sources)
if (any(source_md5_before[protected] != source_md5_after[protected])) {
  stop("A protected frozen source changed unexpectedly during manuscript export.", call. = FALSE)
}
readr::write_csv(tibble::tibble(source = required_sources,
  md5_before = unname(source_md5_before), md5_after = unname(source_md5_after)),
  file.path(audits_dir, "supplemental_manuscript_export_source_manifest.csv"))
message("DONE: S1-S10 regenerated at 170 mm wide, <=225 mm high; PDF + 600 dpi PNG.")
