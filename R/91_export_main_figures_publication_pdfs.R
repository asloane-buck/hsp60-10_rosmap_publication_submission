############################################################
## 91_export_main_figures_publication_pdfs.R
## Fresh-session publication PDF export pipeline.
##
## Run from the repository root with:
##   source("R/91_export_main_figures_publication_pdfs.R")
############################################################

message("\n============================================================")
message("Running ROSMAP publication PDF export pipeline")
message("============================================================")

locate_script_dir <- function() {
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

  if (!is.na(this_file)) dirname(this_file) else normalizePath("R", mustWork = FALSE)
}

script_dir <- locate_script_dir()
old_wd <- getwd()

message("Previous working directory: ", old_wd)
message("Script directory: ", script_dir)
message("Current working directory: ", getwd())

setup_scripts <- c(
  "00_config.R",
  "80_publication_export_helpers.R"
)

core_scripts <- c(
  "01_utils.R",
  "02_load_data.R",
  "03_build_adjusted_core_objects.R",
  "04_build_pathway_sets_all_clients.R",
  "05_build_adjusted_all_client_tables.R",
  "06_build_cognition_objects.R",
  "07_build_matched_null_objects.R"
)

figure_specs <- data.frame(
  figure = 1:6,
  script = c(
    "10_make_main_figure_1_pathway_remodeling.R",
    "11_make_main_figure_2_collapse_heterogeneity.R",
    "12_make_main_figure_3_pathology_coupling.R",
    "13_make_main_figure_4_cognition.R",
    "14_make_main_figure_5_matched_null_specificity.R",
    "15_make_main_figure_6_candidate_classification.R"
  ),
  object_name = c(
    "fig1_all_clients",
    "fig2_all_clients",
    "fig3_all_clients",
    "fig4_integrated",
    "final_fig",
    "fig6_candidate_classification"
  ),
  file = c(
    "Figure_01_manuscript_ready_full_width.pdf",
    "Figure_02_manuscript_ready_full_width.pdf",
    "Figure_03_manuscript_ready_full_width.pdf",
    "Figure_04_manuscript_ready_full_width.pdf",
    "Figure_05_manuscript_ready_full_width.pdf",
    "Figure_06_manuscript_ready_full_width.pdf"
  ),
  width_mm = rep(170, 6),
  height_mm = c(220, 160, 165, 220, 210, 225),
  stringsAsFactors = FALSE
)

required_scripts <- unique(c(setup_scripts, core_scripts, figure_specs$script))
script_paths <- file.path(script_dir, required_scripts)
names(script_paths) <- required_scripts

missing_scripts <- required_scripts[!file.exists(script_paths)]

if (length(missing_scripts) > 0) {
  stop(
    "Missing required script(s) in script directory:\n",
    paste(missing_scripts, collapse = "\n"),
    "\n\nScript directory searched:\n",
    script_dir,
    call. = FALSE
  )
}

message("\nFinal plot objects identified from original functional scripts:")
for (i in seq_len(nrow(figure_specs))) {
  message(
    "  Figure ", figure_specs$figure[[i]], ": ",
    figure_specs$object_name[[i]], " from ", figure_specs$script[[i]]
  )
}

source(script_paths[["00_config.R"]], local = .GlobalEnv)
source(script_paths[["80_publication_export_helpers.R"]], local = .GlobalEnv)

original_cfg_dirs <- mr_capture_output_dirs(cfg)
manuscript_output_dir <- mr_manuscript_ready_output_dir(cfg)
staging_root <- mr_new_staging_root()

mr_sync_global_cfg <- function(new_cfg) {
  assign("cfg", new_cfg, envir = .GlobalEnv)
  invisible(new_cfg)
}

message("\nOriginal manuscript-ready PDF output directory:")
message(manuscript_output_dir)
message("\nTemporary side-effect staging directory for original scripts:")
message(staging_root)

cfg <- mr_use_staged_output_dirs(cfg, staging_root)
mr_sync_global_cfg(cfg)

for (script in core_scripts) {
  mr_source_script(script_paths[[script]], env = .GlobalEnv)
}

ggsave_state <- mr_install_ggsave_skipper(env = .GlobalEnv)
audit_rows <- list()

for (i in seq_len(nrow(figure_specs))) {
  spec <- figure_specs[i, ]

  cfg <- mr_use_staged_output_dirs(cfg, staging_root)
  mr_sync_global_cfg(cfg)
  mr_source_script(script_paths[[spec$script]], env = .GlobalEnv)

  if (!exists(spec$object_name, envir = .GlobalEnv, inherits = FALSE)) {
    stop(
      "Expected final plot object not found after sourcing ",
      spec$script,
      ": ",
      spec$object_name,
      call. = FALSE
    )
  }

  final_plot <- get(spec$object_name, envir = .GlobalEnv, inherits = FALSE)
  manuscript_plot <- mr_make_manuscript_ready_plot(final_plot, tag_levels = "A", figure = spec$figure)

  cfg <- mr_restore_output_dirs(cfg, original_cfg_dirs)
  mr_sync_global_cfg(cfg)
  manuscript_output_dir <- mr_manuscript_ready_output_dir(cfg)

  message("\nExporting manuscript-ready Figure ", spec$figure, ":")
  message("  Object: ", spec$object_name)
  message("  File:   ", file.path(manuscript_output_dir, spec$file))
  message("  Size:   ", spec$width_mm, " mm x ", spec$height_mm, " mm")

  audit_rows[[length(audit_rows) + 1]] <- mr_export_pdf(
    plot = manuscript_plot,
    filename = spec$file,
    output_dir = manuscript_output_dir,
    figure = spec$figure,
    object_name = spec$object_name,
    width_mm = spec$width_mm,
    height_mm = spec$height_mm,
    embed_fonts = TRUE
  )
}

cfg <- mr_restore_output_dirs(cfg, original_cfg_dirs)
mr_sync_global_cfg(cfg)
mr_restore_ggsave(ggsave_state, env = .GlobalEnv)

audit_df <- do.call(rbind, audit_rows)
mr_write_audit_csv(audit_df, manuscript_output_dir)

message("\n============================================================")
message("Manuscript-ready PDF export complete.")
message("============================================================")

print(
  audit_df[
    ,
    c(
      "figure",
      "object_name",
      "path",
      "width_mm",
      "height_mm",
      "file_size_mb",
      "under_10mb",
      "font_embedding_status",
      "pdffonts_message"
    )
  ],
  row.names = FALSE
)

if (any(!audit_df$under_10mb)) {
  warning("One or more manuscript-ready PDFs are >= 10 MB. See manuscript_ready_pdf_audit.csv.", call. = FALSE)
}

message("\nNo copied figure scripts were required for manuscript-ready export.")
message("Original ggsave figure exports were suppressed during this manuscript-ready run.")
