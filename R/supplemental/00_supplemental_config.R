# ============================================================
# Supplemental figure pipeline configuration
#
# Publication-ready configuration. Controlled-access ROSMAP/MSBB data,
# individual-level metadata, and generated intermediate tables are not
# redistributed with this repository. Place local inputs under data/raw/,
# data/external/, or data/supplemental_inputs/, or override paths with
# environment variables.
# ============================================================

locate_repo_root <- function() {
  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )

  if (is.na(this_file) && requireNamespace("rstudioapi", quietly = TRUE)) {
    this_file <- tryCatch(
      normalizePath(rstudioapi::getActiveDocumentContext()$path, mustWork = TRUE),
      error = function(e) NA_character_
    )
  }

  if (!is.na(this_file)) {
    ## This file lives in R/supplemental/.
    return(normalizePath(file.path(dirname(this_file), "..", ".."), mustWork = FALSE))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

project_dir <- Sys.getenv("HSP60_ROSMAP_PROJECT_DIR", unset = locate_repo_root())
data_dir <- file.path(project_dir, "data")
raw_data_dir <- file.path(data_dir, "raw")
metadata_dir <- file.path(raw_data_dir, "metadata")
proteomics_dir <- file.path(raw_data_dir, "proteomics")
derived_dir <- file.path(raw_data_dir, "derived")
external_data_dir <- file.path(data_dir, "external")
supplemental_input_dir <- file.path(data_dir, "supplemental_inputs")
main_outputs_dir <- file.path(project_dir, "outputs", "main_figures")

supplemental_dir <- file.path(project_dir, "outputs", "supplemental_figures")
outputs_dir <- supplemental_dir
figures_dir <- file.path(outputs_dir, "figures")
panels_dir <- file.path(outputs_dir, "panels")
tables_dir <- file.path(outputs_dir, "tables")
audits_dir <- file.path(outputs_dir, "audits")
manuscript_ready_pdf_dir <- file.path(outputs_dir, "manuscript_ready_supplemental_pdfs")
manuscript_ready_png_dir <- file.path(outputs_dir, "manuscript_ready_supplemental_png_qc")

# Journal supplemental figure export constraints.
SUPP_FULL_WIDTH_MM <- 170
SUPP_MAX_HEIGHT_MM <- 225
MM_TO_IN <- 1 / 25.4
SUPP_FULL_WIDTH_IN <- SUPP_FULL_WIDTH_MM * MM_TO_IN

SUPP_FINAL_FIGURE_HEIGHT_MM <- c(
  Supplementary_Figure_1_cohort_detection = 139,
  Supplementary_Figure_2_matched_individual_sensitivity = 85,
  Supplementary_Figure_3_mitochondrial_specificity = 140,
  Supplementary_Figure_4_pathology_model_robustness = 119,
  Supplementary_Figure_5_msbb_cross_cohort_validation = 119,
  Supplementary_Figure_6_regional_proteomics_validation = 120
)

if (any(SUPP_FINAL_FIGURE_HEIGHT_MM > SUPP_MAX_HEIGHT_MM)) {
  stop("Configured supplemental figure height exceeds SUPP_MAX_HEIGHT_MM.", call. = FALSE)
}

output_dirs <- c(
  supplemental_dir,
  outputs_dir,
  figures_dir,
  panels_dir,
  tables_dir,
  audits_dir,
  manuscript_ready_pdf_dir,
  manuscript_ready_png_dir,
  metadata_dir,
  proteomics_dir,
  derived_dir,
  external_data_dir,
  supplemental_input_dir
)

for (dir_i in output_dirs) {
  if (!dir.exists(dir_i)) {
    dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
  }
}

# ============================================================
# Package setup
# ============================================================

required_packages <- c(
  "tidyverse",
  "ggplot2",
  "patchwork",
  "readr",
  "stringr",
  "scales"
)

optional_packages <- c("cowplot", "janitor")

missing_required_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_required_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(scales)
})

for (pkg_i in optional_packages) {
  if (requireNamespace(pkg_i, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg_i, character.only = TRUE))
  } else {
    message("Optional package not installed: ", pkg_i)
  }
}

# ============================================================
# Editable input candidates
# ============================================================

# Set overwrite_inputs <- TRUE if a later source should replace an earlier one.
overwrite_inputs <- FALSE

# Optional workspaces are loaded into a temporary environment, and only object
# names listed in candidate_object_names are retained in inputs. Do not commit
# these workspaces.
candidate_rdata_paths <- c(
  main_analysis_workspace = Sys.getenv(
    "SUPP_MAIN_ANALYSIS_RDATA",
    unset = file.path(supplemental_input_dir, "main_analysis_workspace.RData")
  ),
  robustness_workspace = Sys.getenv(
    "SUPP_ROBUSTNESS_RDATA",
    unset = file.path(supplemental_input_dir, "robustness_workspace.RData")
  )
)

candidate_rds_paths <- c(
  rna_vst = Sys.getenv(
    "SUPP_RNA_VST_RDS",
    unset = file.path(supplemental_input_dir, "ROSMAP_vsd.rds")
  ),
  rna_mat = Sys.getenv(
    "SUPP_RNA_MAT_RDS",
    unset = file.path(main_outputs_dir, "objects", "rna_mat_adj.rds")
  ),
  rna_meta = Sys.getenv(
    "SUPP_RNA_META_RDS",
    unset = file.path(main_outputs_dir, "objects", "rna_meta_adj.rds")
  ),
  prot_mat = Sys.getenv(
    "SUPP_PROT_MAT_RDS",
    unset = file.path(main_outputs_dir, "objects", "prot_mat_adj.rds")
  ),
  prot_mat_raw = Sys.getenv(
    "SUPP_PROT_MAT_RAW_RDS",
    unset = file.path(main_outputs_dir, "objects", "prot_mat_raw.rds")
  ),
  prot_meta = Sys.getenv(
    "SUPP_PROT_META_RDS",
    unset = file.path(main_outputs_dir, "objects", "prot_meta_adj.rds")
  )
)

agora_target_file <- Sys.getenv(
  "AGORA_TARGET_FILE",
  unset = file.path(external_data_dir, "AGORA_nominated_targets.csv")
)

msbb_validation_paths <- c(
  Sys.getenv("MSBB_VALIDATION_FILE", unset = ""),
  file.path(supplemental_input_dir, "MSBB_covariate_adjusted_braak_and_collapse.csv"),
  file.path(supplemental_input_dir, "MSBB_standardized_braak_collapse_percentiles.csv"),
  file.path(supplemental_input_dir, "MSBB_covariate_adjusted_braak_collapse_percentiles.csv"),
  file.path(supplemental_input_dir, "COMBINED_ROSMAP_MSBB_validation_results.csv")
)

candidate_csv_paths <- list(
  hsp60_hsp10_interactor_inventory = c(
    file.path(main_outputs_dir, "pathway_gene_sets_all_clients", "Hsp60_10_all_clients.csv"),
    file.path(supplemental_input_dir, "Hsp60_10_all_clients.csv"),
    file.path(supplemental_input_dir, "all_hsp60_10_clients_from_supplement_cleaned.csv")
  ),
  mito_background_tbl = c(
    file.path(main_outputs_dir, "pathway_gene_sets_all_clients", "Broad_MitoCarta_non_Hsp60_10.csv"),
    file.path(supplemental_input_dir, "Broad_MitoCarta_non_Hsp60_10.csv"),
    file.path(supplemental_input_dir, "Broad_MitoCarta_non_Hsp60_10_NO_Hsp60_10_overlap.csv")
  ),
  tca_pyruvate_tbl = c(
    file.path(main_outputs_dir, "pathway_gene_sets_all_clients", "TCA_pyruvate_metabolism_non_Hsp60_10.csv"),
    file.path(supplemental_input_dir, "TCA_pyruvate_metabolism_non_Hsp60_10.csv")
  ),
  prot_mat = character(),
  protein_matrix = character(),
  prot_meta = character(),
  protein_meta = character(),
  protein_sample_meta = c(
    file.path(proteomics_dir, "rosmap_50batch_specimen_metadata_for_batch_correction.csv"),
    file.path(proteomics_dir, "matched_metadata.csv"),
    file.path(supplemental_input_dir, "matched_metadata_copy.csv")
  ),
  protein_clinical_meta = c(
    file.path(metadata_dir, "ROSMAP_clinical.csv")
  ),
  protein_biospecimen_meta = c(
    file.path(metadata_dir, "ROSMAP_biospecimen_metadata.csv")
  ),
  regional_protein_screen = c(
    file.path(supplemental_input_dir, "AMP_covariate_adjusted_regional_screen_labeled_Hsp60_mito.csv")
  ),
  braak_models = c(
    file.path(supplemental_input_dir, "AMP_covariate_adjusted_Braak_models_by_region.csv")
  ),
  hsp_region_compare = c(
    file.path(supplemental_input_dir, "AMP_covariate_adjusted_Hsp60_clients_DLPFC_vs_STG.csv")
  ),
  msbb_braak_effects = c(
    file.path(supplemental_input_dir, "MSBB_braak_effects.csv")
  ),
  msbb_collapse = c(
    file.path(supplemental_input_dir, "MSBB_collapse.csv")
  ),
  rna_mat = character(),
  rna_meta = character(),
  per_gene_results = c(
    file.path(main_outputs_dir, "tables", "main_fig4_cognition", "client_level_cognition_gene_summary.csv"),
    file.path(supplemental_input_dir, "per_gene_results.csv")
  ),
  priority_tbl = c(
    file.path(main_outputs_dir, "tables", "priority_tbl_COVARIATE_ADJUSTED.csv"),
    file.path(supplemental_input_dir, "priority_tbl_COVARIATE_ADJUSTED.csv")
  )
)

candidate_object_names <- c(
  "hsp60_hsp10_interactor_inventory",
  "hsp60_client_tbl",
  "hsp_client_tbl",
  "all_hsp60_10_client_tbl",
  "hsp_clients_all",
  "hsp60_clients",
  "prot_mat",
  "protein_matrix",
  "protein_mat",
  "tmt_mat",
  "prot_expr_mat",
  "prot_mat_in",
  "prot_mat_raw",
  "prot_meta",
  "prot_meta_path",
  "meta_protein_sub",
  "meta_protein",
  "protein_meta",
  "protein_meta_raw",
  "protein_metadata",
  "prot_metadata",
  "protein_covariates",
  "rna_mat",
  "vst_mat",
  "rna_vst",
  "rna_matrix",
  "vst_symbol_mat",
  "rna_meta",
  "meta_rna_sub",
  "rna_metadata",
  "mito_background_tbl",
  "mitochondrial_background_tbl",
  "per_gene_results",
  "priority_tbl",
  "background_null_pool",
  "hsp_null_tbl",
  "null_results_tbl",
  "matched_null_results",
  "mito_null_results",
  "regional_protein_screen",
  "braak_models",
  "hsp_region_compare",
  "msbb_braak_effects",
  "msbb_collapse"
)

message("Supplemental figure directory: ", supplemental_dir)
message("Output directory: ", outputs_dir)
