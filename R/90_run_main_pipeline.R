############################################################
## 90_run_main_pipeline.R
## Fresh-session main analysis and figure pipeline.
##
## Run from the repository root with:
##   source("R/90_run_main_pipeline.R")
############################################################

message("\n============================================================")
message("Running ROSMAP main figure pipeline")
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

required_scripts <- c(
  "00_config.R",
  "01_utils.R",
  "02_load_data.R",
  "03_build_adjusted_core_objects.R",
  "04_build_pathway_sets_all_clients.R",
  "05_build_adjusted_all_client_tables.R",
  "06_build_cognition_objects.R",
  "07_build_matched_null_objects.R",
  "10_make_main_figure_1_pathway_remodeling.R",
  "11_make_main_figure_2_collapse_heterogeneity.R",
  "12_make_main_figure_3_pathology_coupling.R",
  "13_make_main_figure_4_cognition.R",
  "14_make_main_figure_5_matched_null_specificity.R",
  "15_make_main_figure_6_candidate_classification.R"
)

script_paths <- file.path(script_dir, required_scripts)
names(script_paths) <- required_scripts

missing_scripts <- required_scripts[!file.exists(script_paths)]

if (length(missing_scripts) > 0) {
  stop(
    "Missing required script(s) in script directory:\n",
    paste(missing_scripts, collapse = "\n"),
    "\n\nScript directory searched:\n",
    script_dir,
    "\n\nAvailable .R files here:\n",
    paste(list.files(script_dir, pattern = "\\.R$"), collapse = "\n"),
    call. = FALSE
  )
}

for (script in required_scripts) {
  message("\n------------------------------------------------------------")
  message("Sourcing: ", script)
  message("------------------------------------------------------------")

  tryCatch(
    source(script_paths[[script]], local = .GlobalEnv),
    error = function(e) {
      stop(
        "\nPipeline failed while sourcing: ", script,
        "\n\nOriginal error:\n",
        conditionMessage(e),
        "\n\nCurrent working directory:\n",
        getwd(),
        "\n\nAvailable .R files:\n",
        paste(list.files(script_dir, pattern = "\\.R$"), collapse = "\n"),
        call. = FALSE
      )
    }
  )
}

message("\n============================================================")
message("All six main figure scripts completed.")
message("============================================================")
