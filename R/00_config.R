############################################################
## 00_config.R
## Project paths, run flags, output folders
##
## Publication-ready configuration.
## Controlled-access ROSMAP data are not redistributed with this repository.
## Place required input files under data/raw/ and data/external/ as documented
## in data/README.md, or override paths using environment variables.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(readxl)
  library(broom)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

set.seed(123)
options(stringsAsFactors = FALSE)

locate_project_dir <- function() {
  find_repo_root <- function(start_path) {
    if (is.na(start_path) || !nzchar(start_path)) {
      return(NA_character_)
    }

    start_path <- normalizePath(
      start_path,
      mustWork = FALSE
    )

    if (!dir.exists(start_path)) {
      start_path <- dirname(start_path)
    }

    current <- start_path

    repeat {
      if (
        file.exists(file.path(current, "R", "00_config.R")) &&
        file.exists(file.path(current, "R", "01_utils.R"))
      ) {
        return(normalizePath(current, mustWork = FALSE))
      }

      parent <- dirname(current)

      if (identical(parent, current)) {
        break
      }

      current <- parent
    }

    NA_character_
  }

  candidates <- character()

  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )

  if (!is.na(this_file)) {
    candidates <- c(candidates, this_file)
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    active_file <- tryCatch(
      normalizePath(
        rstudioapi::getActiveDocumentContext()$path,
        mustWork = TRUE
      ),
      error = function(e) NA_character_
    )

    if (!is.na(active_file)) {
      candidates <- c(candidates, active_file)
    }
  }

  candidates <- c(candidates, getwd())

  for (candidate in unique(candidates)) {
    root <- find_repo_root(candidate)

    if (!is.na(root)) {
      return(root)
    }
  }

  stop(
    paste0(
      "Could not locate repository root. Expected to find ",
      "R/00_config.R and R/01_utils.R in a parent directory. ",
      "Set HSP60_ROSMAP_PROJECT_DIR explicitly if needed."
    )
  )
}

cfg <- list()

## Repository root. Override with HSP60_ROSMAP_PROJECT_DIR if needed.
cfg$project_dir <- Sys.getenv("HSP60_ROSMAP_PROJECT_DIR", unset = locate_project_dir())

## Input folders. Controlled-access files should be placed locally only and
## remain ignored by Git.
cfg$data_dir <- file.path(cfg$project_dir, "data")
cfg$raw_data_dir <- file.path(cfg$data_dir, "raw")
cfg$metadata_dir <- file.path(cfg$raw_data_dir, "metadata")
cfg$derived_dir <- file.path(cfg$raw_data_dir, "derived")
cfg$proteomics_input_dir <- file.path(cfg$raw_data_dir, "proteomics")
cfg$external_data_dir <- file.path(cfg$data_dir, "external")

## Output folders.
cfg$output_root <- file.path(cfg$project_dir, "outputs", "main_figures")
cfg$table_dir <- file.path(cfg$output_root, "tables")
cfg$plot_dir <- file.path(cfg$output_root, "plots")
cfg$object_dir <- file.path(cfg$output_root, "objects")
cfg$pathway_dir <- file.path(cfg$output_root, "pathway_gene_sets_all_clients")

for (path in c(
  cfg$table_dir,
  cfg$plot_dir,
  cfg$object_dir,
  cfg$pathway_dir,
  cfg$metadata_dir,
  cfg$derived_dir,
  cfg$proteomics_input_dir,
  cfg$external_data_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

## External/public annotation files. These are not automatically redistributed;
## see data/README.md for expected filenames and source notes.
cfg$hsp60_client_file <- Sys.getenv(
  "HSP60_CLIENT_XLSX",
  unset = file.path(cfg$external_data_dir, "1-s2.0-S1355814523012117-MOESM4_ESM.xlsx")
)
cfg$hsp60_client_sheet <- "Supplementary Table S1"
cfg$hsp60_client_skip <- 2

cfg$mitocarta_xls <- Sys.getenv(
  "MITOCARTA_XLS",
  unset = file.path(cfg$external_data_dir, "Human.MitoCarta3.0.xls")
)

## Controlled-access ROSMAP and derived inputs. Do not commit these files.
cfg$rosmap_clinical_file <- Sys.getenv(
  "ROSMAP_CLINICAL_FILE",
  unset = file.path(cfg$metadata_dir, "ROSMAP_clinical.csv")
)

## External target nomination file used in candidate-prioritization analyses.
cfg$agora_target_file <- Sys.getenv(
  "AGORA_TARGET_FILE",
  unset = file.path(cfg$external_data_dir, "AGORA_nominated_targets.csv")
)

## Main statistical choices.
cfg$adjust_main_figures <- TRUE
cfg$adjust_braak_for_cerad <- FALSE
cfg$min_n_gene_model <- 30
cfg$min_n_stage_model <- 10
cfg$null_n_iter <- 10000

## Figure export.
cfg$png_dpi <- 600

write_tbl <- function(x, name, subdir = cfg$table_dir) {
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(subdir, paste0(name, ".csv"))
  readr::write_csv(x, out)
  message("Wrote table: ", out)
  invisible(x)
}

save_obj <- function(x, name) {
  out <- file.path(cfg$object_dir, paste0(name, ".rds"))
  saveRDS(x, out)
  message("Wrote object: ", out)
  invisible(x)
}

read_obj <- function(name) {
  readRDS(file.path(cfg$object_dir, paste0(name, ".rds")))
}

save_plot <- function(plot, name, width = 12, height = 8, subdir = cfg$plot_dir, dpi = cfg$png_dpi) {
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

  pdf_out <- file.path(subdir, paste0(name, ".pdf"))
  png_out <- file.path(subdir, paste0(name, ".png"))

  ggsave(pdf_out, plot, width = width, height = height, units = "in",
         bg = "white", device = grDevices::pdf, limitsize = FALSE, useDingbats = FALSE)

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_out, plot, width = width, height = height, units = "in",
           dpi = dpi, bg = "white", device = ragg::agg_png, limitsize = FALSE)
  } else {
    ggsave(png_out, plot, width = width, height = height, units = "in",
           dpi = dpi, bg = "white", limitsize = FALSE)
  }

  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_out <- file.path(subdir, paste0(name, ".svg"))
    ggsave(svg_out, plot, width = width, height = height, units = "in",
           bg = "white", device = svglite::svglite, limitsize = FALSE)
  }

  message("Saved plot: ", name)
  invisible(plot)
}

message("Loaded 00_config.R")
