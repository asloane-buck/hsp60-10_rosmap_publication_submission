############################################################
## 90_run_supplemental_pipeline.R
##
## Deterministic supplemental-figure runner.
##
## Prerequisite:
##   Run the corrected main pipeline through Main Figure 5 first.
##
## The supplemental pipeline consumes the canonical production
## Figure 5 null universe and gene-level cognition summary directly.
## It does NOT reconstruct Hsp60/10/background metrics from the
## AMP-AD regional screen, and it does NOT rediscover cognition files
## by scanning the repository.
############################################################

options(stringsAsFactors = FALSE)

message("Starting corrected supplemental figure pipeline (S1-S10).")
supplemental_regenerated_figures <- character()

script_dir <- normalizePath(
  file.path(getwd(), "R", "supplemental"),
  mustWork = TRUE
)

required_scripts <- c(
  "00_supplemental_config.R",
  "01_supplemental_load_inputs.R",
  "02_supplemental_helper_functions.R",
  "10_make_supplementary_figure_1_cohort_detection.R",
  "11_make_supplementary_figure_2_matched_individual_sensitivity.R",
  "12_make_supplementary_figure_3_mitochondrial_specificity.R",
  "13_make_supplementary_figure_4_pathology_model_robustness.R",
  "14_make_supplementary_figure_5_msbb_cross_cohort_validation.R",
  "15_make_supplementary_figure_6_regional_proteomics_validation.R",
  "16_make_supplementary_figure_7_conventional_differential.R",
  "17_make_supplementary_figure_8_PC1_sensitivity.R",
  "18_make_supplementary_figure_9_variance_partition.R",
  "19_make_supplementary_figure_10_marker_sensitivity.R"
)

missing_scripts <- required_scripts[
  !file.exists(file.path(script_dir, required_scripts))
]

if (length(missing_scripts) > 0) {
  stop(
    "Missing supplemental script(s): ",
    paste(missing_scripts, collapse = ", "),
    call. = FALSE
  )
}

source(file.path(script_dir, "00_supplemental_config.R"))
source(file.path(script_dir, "01_supplemental_load_inputs.R"))
source(file.path(script_dir, "02_supplemental_helper_functions.R"))

production_null_path <- file.path(
  main_outputs_dir,
  "tables",
  "figure5_null_input_COVARIATE_ADJUSTED.csv"
)

production_cognition_path <- file.path(
  main_outputs_dir,
  "tables",
  "main_fig5_matched_null_specificity",
  "Fig5_gene_universe_cognition_summary_by_gene.csv"
)

required_production_files <- c(
  production_null_path,
  production_cognition_path
)

missing_production_files <- required_production_files[
  !file.exists(required_production_files)
]

if (length(missing_production_files) > 0) {
  stop(
    "Corrected main-pipeline prerequisite output(s) are missing:\n",
    paste0(" - ", missing_production_files, collapse = "\n"),
    "\nRun the main pipeline through Main Figure 5 before the supplemental runner.",
    call. = FALSE
  )
}

production_null_tbl <- readr::read_csv(
  production_null_path,
  show_col_types = FALSE
)

production_cognition_tbl <- readr::read_csv(
  production_cognition_path,
  show_col_types = FALSE
)

required_null_cols <- c(
  "gene",
  "is_hsp60_10_client",
  "protein_late_decline_magnitude",
  "braak_pathology_magnitude",
  "cerad_pathology_magnitude",
  "joint_pathology_magnitude",
  "strict_joint_pathology_magnitude",
  "matching_abundance",
  "abundance_bin",
  "agora_nominated_target"
)

missing_null_cols <- setdiff(
  required_null_cols,
  colnames(production_null_tbl)
)

if (length(missing_null_cols) > 0) {
  stop(
    "Production Figure 5 null input is missing corrected column(s): ",
    paste(missing_null_cols, collapse = ", "),
    call. = FALSE
  )
}

required_cognition_cols <- c(
  "gene",
  "cognition_priority_score"
)

missing_cognition_cols <- setdiff(
  required_cognition_cols,
  colnames(production_cognition_tbl)
)

if (length(missing_cognition_cols) > 0) {
  stop(
    "Production Figure 5 cognition summary is missing column(s): ",
    paste(missing_cognition_cols, collapse = ", "),
    call. = FALSE
  )
}

production_null_tbl <- production_null_tbl |>
  dplyr::mutate(
    gene_original_production = as.character(.data$gene),
    gene = clean_gene(.data$gene)
  )

production_cognition_tbl <- production_cognition_tbl |>
  dplyr::transmute(
    gene_original_cognition = as.character(.data$gene),
    gene = clean_gene(.data$gene),
    cognition_priority_score =
      as.numeric(.data$cognition_priority_score)
  )

if (anyDuplicated(production_null_tbl$gene) > 0) {
  stop(
    "Production null input has duplicate genes after supplemental gene normalization.",
    call. = FALSE
  )
}

if (anyDuplicated(production_cognition_tbl$gene) > 0) {
  stop(
    "Production cognition summary has duplicate genes after supplemental gene normalization.",
    call. = FALSE
  )
}

supplemental_null_tbl <- production_null_tbl |>
  dplyr::left_join(
    production_cognition_tbl |>
      dplyr::select("gene", "cognition_priority_score"),
    by = "gene"
  )

if (nrow(supplemental_null_tbl) != 915) {
  stop(
    "Corrected supplemental null universe must contain 915 proteins; observed ",
    nrow(supplemental_null_tbl),
    ".",
    call. = FALSE
  )
}

if (sum(supplemental_null_tbl$is_hsp60_10_client, na.rm = TRUE) != 306) {
  stop(
    "Corrected supplemental null universe must contain 306 Hsp60/10 clients.",
    call. = FALSE
  )
}

if (sum(!supplemental_null_tbl$is_hsp60_10_client, na.rm = TRUE) != 609) {
  stop(
    "Corrected supplemental null universe must contain 609 non-client mitochondrial proteins.",
    call. = FALSE
  )
}

n_hsp_cognition <- sum(
  supplemental_null_tbl$is_hsp60_10_client &
    is.finite(supplemental_null_tbl$cognition_priority_score),
  na.rm = TRUE
)

n_background_cognition <- sum(
  !supplemental_null_tbl$is_hsp60_10_client &
    is.finite(supplemental_null_tbl$cognition_priority_score),
  na.rm = TRUE
)

if (n_hsp_cognition < 5 || n_background_cognition < 5) {
  stop(
    "Too few corrected Figure 5 genes have cognition_priority_score: ",
    "Hsp60/10=", n_hsp_cognition,
    ", background=", n_background_cognition,
    ". Re-run corrected Main Figure 5.",
    call. = FALSE
  )
}

hsp_null_tbl <- supplemental_null_tbl |>
  dplyr::filter(.data$is_hsp60_10_client)

background_null_pool <- supplemental_null_tbl |>
  dplyr::filter(!.data$is_hsp60_10_client)

inputs$hsp_null_tbl <- hsp_null_tbl
inputs$background_null_pool <- background_null_pool

if ("protein_meta" %in% names(inputs)) protein_meta <- inputs$protein_meta
if ("prot_meta" %in% names(inputs)) prot_meta <- inputs$prot_meta
if ("prot_meta_aligned" %in% names(inputs)) prot_meta_aligned <- inputs$prot_meta_aligned
if ("prot_mat" %in% names(inputs)) prot_mat <- inputs$prot_mat
if ("prot_mat_raw" %in% names(inputs)) prot_mat_raw <- inputs$prot_mat_raw
if ("rna_meta" %in% names(inputs)) rna_meta <- inputs$rna_meta
if ("rna_mat" %in% names(inputs)) rna_mat <- inputs$rna_mat

readr::write_csv(
  supplemental_null_tbl,
  file.path(
    audits_dir,
    "CANONICAL_PRODUCTION_NULL_INPUT_USED_BY_SUPPLEMENTAL_PIPELINE.csv"
  )
)

message("Canonical supplemental null inputs loaded from corrected production outputs:")
message("  Hsp60/10 clients: ", nrow(hsp_null_tbl))
message("  non-client mitochondrial proteins: ", nrow(background_null_pool))
message(
  "  cognition scores available: Hsp60/10=",
  n_hsp_cognition,
  "; background=",
  n_background_cognition
)

source(file.path(
  script_dir,
  "10_make_supplementary_figure_1_cohort_detection.R"
))

source(file.path(
  script_dir,
  "11_make_supplementary_figure_2_matched_individual_sensitivity.R"
))

# Freeze the revised cohort/detection and strictly matched-panel populations.
if (isTRUE(get0("SUPP_USE_FROZEN_MODELS", ifnotfound = FALSE))) {
  detected <- readr::read_csv(file.path(audits_dir, "SuppFig1_interactor_detection_audit.csv"), show_col_types = FALSE)
  if (nrow(detected) != 321 || sum(detected$detected_rna, na.rm = TRUE) != 297 ||
      sum(detected$detected_protein, na.rm = TRUE) != 306 ||
      sum(detected$detected_rna & detected$detected_protein, na.rm = TRUE) != 285) {
    stop("S1 detected coverage differs from the verified 321-reference/297-RNA/306-protein inventory.", call. = FALSE)
  }
  stages <- readr::read_csv(file.path(audits_dir, "SuppFig1_plotted_primary_stage_counts.csv"), show_col_types = FALSE)
  expected <- c(RNA_NCI = 200, RNA_MCI = 158, RNA_AD = 219,
    Protein_NCI = 167, Protein_MCI = 96, Protein_AD = 109)
  actual <- setNames(stages$n_samples, paste(stages$modality, stages$stage, sep = "_"))[names(expected)]
  if (anyNA(actual) || any(actual != expected)) stop("S1 primary stage counts differ from the frozen analysis.", call. = FALSE)
  matched <- supfig2_outputs$matched_summary
  expected_n <- c(NCI = 94, MCI = 56, AD = 48)
  if (length(supfig2_outputs$overlap_ids) != 198 || nrow(matched) != 18 ||
      anyNA(matched$n) || any(matched$n != expected_n[as.character(matched$Stage)])) {
    stop("S2 must plot the same 198 matched participants (94 NCI/56 MCI/48 AD) in every pathway and modality.", call. = FALSE)
  }
}

source(file.path(
  script_dir,
  "12_make_supplementary_figure_3_mitochondrial_specificity.R"
))

source(file.path(
  script_dir,
  "13_make_supplementary_figure_4_pathology_model_robustness.R"
))

source(file.path(
  script_dir,
  "14_make_supplementary_figure_5_msbb_cross_cohort_validation.R"
))

source(file.path(
  script_dir,
  "15_make_supplementary_figure_6_regional_proteomics_validation.R"
))

source(file.path(script_dir, "16_make_supplementary_figure_7_conventional_differential.R"))
source(file.path(script_dir, "17_make_supplementary_figure_8_PC1_sensitivity.R"))
source(file.path(script_dir, "18_make_supplementary_figure_9_variance_partition.R"))
source(file.path(script_dir, "19_make_supplementary_figure_10_marker_sensitivity.R"))

required_output_objects <- c(
  "supfig1_outputs",
  "supfig2_outputs",
  "supfig3_outputs",
  "supfig4_outputs",
  "supfig5_outputs",
  "supfig6_outputs", "supfig7_outputs", "supfig8_outputs",
  "supfig9_outputs", "supfig10_outputs"
)

missing_output_objects <- required_output_objects[
  !vapply(
    required_output_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_output_objects) > 0) {
  stop(
    "Supplemental figure object(s) missing after run: ",
    paste(missing_output_objects, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(supfig3_outputs$hsp_tbl) != 306) {
  stop("Supp Fig 3 did not retain 306 Hsp60/10 clients.", call. = FALSE)
}

if (nrow(supfig3_outputs$background_tbl) != 609) {
  stop("Supp Fig 3 did not retain 609 background proteins.", call. = FALSE)
}

if (nrow(supfig5_outputs$rosmap_tbl) != 915) {
  stop("Supp Fig 5 ROSMAP discovery table is not 915 genes.", call. = FALSE)
}

copy_manuscript_ready_supplemental_figures <- function() {
  expected <- names(SUPP_FINAL_FIGURE_HEIGHT_MM)
  if (!setequal(expected, supplemental_regenerated_figures)) {
    stop("All ten final figures must be freshly exported during this run; missing: ",
      paste(setdiff(expected, supplemental_regenerated_figures), collapse = ", "), call. = FALSE)
  }
  ensure_dir(manuscript_ready_pdf_dir)
  ensure_dir(manuscript_ready_png_dir)
  pdf_names <- paste0(expected, ".pdf")
  png_names <- paste0(expected, ".png")
  source_pdf <- file.path(figures_dir, pdf_names)
  source_png <- file.path(figures_dir, png_names)
  if (any(!file.exists(c(source_pdf, source_png))) ||
      any(file.info(c(source_pdf, source_png))$size <= 0)) {
    stop("Missing or empty newly exported supplemental files.", call. = FALSE)
  }
  # PNG IHDR: verify actual pixel dimensions independently of export arguments.
  read_png_size <- function(path) {
    con <- file(path, "rb")
    on.exit(close(con))
    header <- readBin(con, "raw", n = 24)
    if (length(header) != 24 || !identical(header[1:8], as.raw(c(137,80,78,71,13,10,26,10)))) {
      stop("Invalid PNG: ", path, call. = FALSE)
    }
    c(sum(as.numeric(header[17:20]) * 256^(3:0)),
      sum(as.numeric(header[21:24]) * 256^(3:0)))
  }
  dims <- vapply(source_png, read_png_size, numeric(2))
  expected_width <- 600 * SUPP_FULL_WIDTH_MM / 25.4
  expected_height <- 600 * unname(SUPP_FINAL_FIGURE_HEIGHT_MM) / 25.4
  if (any(abs(dims[1,] - expected_width) > 2) || any(abs(dims[2,] - expected_height) > 2)) {
    stop("Supplemental PNG pixel dimensions do not match 170 mm / 600 dpi exports.", call. = FALSE)
  }
  pdf_path <- file.path(manuscript_ready_pdf_dir, pdf_names)
  png_path <- file.path(manuscript_ready_png_dir, png_names)
  # Overwrite only these ten expected filenames; preserve unrelated files.
  if (!all(file.copy(source_pdf, pdf_path, overwrite = TRUE)) ||
      !all(file.copy(source_png, png_path, overwrite = TRUE))) {
    stop("Failed to copy manuscript-ready supplemental exports.", call. = FALSE)
  }
  manifest <- tibble::tibble(
    supplemental_figure = paste0("S", seq_along(expected)),
    pdf_path = pdf_path, png_path = png_path,
    width_mm = SUPP_FULL_WIDTH_MM, height_mm = unname(SUPP_FINAL_FIGURE_HEIGHT_MM),
    png_dpi = 600, png_width_px = dims[1,], png_height_px = dims[2,],
    pdf_md5 = unname(tools::md5sum(pdf_path)), png_md5 = unname(tools::md5sum(png_path))
  )
  readr::write_csv(manifest, file.path(manuscript_ready_pdf_dir, "manifest.csv"))
  readr::write_csv(manifest, file.path(manuscript_ready_png_dir, "manifest.csv"))
  invisible(manifest)
}
manuscript_ready_pdf_manifest <- copy_manuscript_ready_supplemental_figures()

if (exists(
  "audit_manuscript_ready_supplemental_pdfs",
  mode = "function"
)) {
  manuscript_ready_pdf_audit <-
    audit_manuscript_ready_supplemental_pdfs()
}

git_branch <- tryCatch(
  system2(
    "git",
    c("branch", "--show-current"),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) NA_character_
)

git_sha <- tryCatch(
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

provenance_lines <- c(
  "Corrected supplemental figure pipeline provenance",
  paste0("generated_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("project_dir: ", project_dir),
  paste0("git_branch: ", paste(git_branch, collapse = " ")),
  paste0("git_sha: ", paste(git_sha, collapse = " ")),
  paste0(
    "working_tree_dirty: ",
    ifelse(length(git_status) > 0 && any(nzchar(git_status)), "TRUE", "FALSE")
  ),
  paste0("production_null_input: ", production_null_path),
  paste0("production_cognition_input: ", production_cognition_path),
  "supplemental_null_universe: 915",
  "hsp60_10_clients: 306",
  "non_client_mitochondrial_background: 609"
)

writeLines(
  provenance_lines,
  file.path(
    audits_dir,
    "supplemental_pipeline_provenance.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    audits_dir,
    "supplemental_pipeline_sessionInfo.txt"
  )
)

message("Supplemental figure pipeline complete.")
message("Figures: ", figures_dir)
message("Panels: ", panels_dir)
message("Audits: ", audits_dir)
message("Tables: ", tables_dir)
message("Manuscript-ready PDFs: ", manuscript_ready_pdf_dir)
message("Manuscript-ready 600 dpi PNGs: ", manuscript_ready_png_dir)
