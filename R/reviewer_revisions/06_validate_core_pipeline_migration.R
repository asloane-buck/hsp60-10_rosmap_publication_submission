############################################################
## 06_validate_core_pipeline_migration.R
##
## Checkpoint after the upstream ROSMAP protein ingestion/stage migration.
## This script intentionally sources ONLY R/00 through R/03.
## It must pass before any downstream manuscript analysis is rerun.
############################################################

options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")
source("R/02_load_data.R")
source("R/03_build_adjusted_core_objects.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "core_pipeline_migration"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

failures <- character()

expect_equal <- function(label, observed, expected) {
  ok <- identical(as.integer(observed), as.integer(expected))
  if (!ok) {
    failures <<- c(
      failures,
      paste0(label, ": observed ", observed, ", expected ", expected)
    )
  }
  tibble::tibble(
    check = label,
    observed = as.character(observed),
    expected = as.character(expected),
    passed = ok
  )
}

expect_true <- function(label, condition, observed = condition) {
  ok <- isTRUE(condition)
  if (!ok) {
    failures <<- c(
      failures,
      paste0(label, ": failed (observed ", observed, ")")
    )
  }
  tibble::tibble(
    check = label,
    observed = as.character(observed),
    expected = "TRUE",
    passed = ok
  )
}

protein_primary_counts <- prot_meta_path |>
  dplyr::filter(!is.na(.data$clinical_stage)) |>
  dplyr::count(.data$clinical_stage, name = "n")

rna_primary_counts <- rna_meta |>
  dplyr::count(.data$clinical_stage, name = "n")

protein_has_adjusted_data <- colSums(is.finite(prot_mat_adj)) > 0
protein_nuisance_complete <- stats::complete.cases(
  prot_meta_adj[, protein_covars, drop = FALSE]
)
protein_primary_analyzable <-
  !is.na(prot_meta_adj$clinical_stage) &
  protein_nuisance_complete &
  protein_has_adjusted_data

protein_analyzable_counts <- tibble::tibble(
  stage = prot_meta_adj$clinical_stage,
  analyzable = protein_primary_analyzable
) |>
  dplyr::filter(!is.na(.data$stage)) |>
  dplyr::group_by(.data$stage, .drop = FALSE) |>
  dplyr::summarise(
    n = sum(.data$analyzable),
    .groups = "drop"
  )

get_n <- function(tbl, stage_name) {
  stage_col <- colnames(tbl)[1]
  keep <- as.character(tbl[[stage_col]]) == stage_name
  hit <- tbl$n[keep]
  if (length(hit) == 0) 0L else as.integer(hit[[1]])
}

checks <- dplyr::bind_rows(
  expect_equal("Full TMT matrix samples", ncol(prot_mat_raw), 400),
  expect_equal("Protein metadata rows", nrow(prot_meta_path), 400),
  expect_equal(
    "Unique protein participants",
    dplyr::n_distinct(prot_meta_path$IndividualID),
    400
  ),
  expect_equal(
    "Protein participants with cogdx",
    sum(!is.na(prot_meta_path$cogdx_num)),
    400
  ),
  expect_equal(
    "Protein primary-stage NCI eligible",
    get_n(protein_primary_counts, "NCI"),
    168
  ),
  expect_equal(
    "Protein primary-stage MCI eligible",
    get_n(protein_primary_counts, "MCI"),
    97
  ),
  expect_equal(
    "Protein primary-stage AD eligible",
    get_n(protein_primary_counts, "AD"),
    109
  ),
  expect_equal(
    "RNA primary-stage NCI",
    get_n(rna_primary_counts, "NCI"),
    200
  ),
  expect_equal(
    "RNA primary-stage MCI",
    get_n(rna_primary_counts, "MCI"),
    158
  ),
  expect_equal(
    "RNA primary-stage AD",
    get_n(rna_primary_counts, "AD"),
    220
  ),
  expect_equal(
    "Protein samples with adjusted data",
    sum(protein_has_adjusted_data),
    398
  ),
  expect_equal(
    "Protein primary-stage NCI analyzable",
    get_n(protein_analyzable_counts, "NCI"),
    167
  ),
  expect_equal(
    "Protein primary-stage MCI analyzable",
    get_n(protein_analyzable_counts, "MCI"),
    96
  ),
  expect_equal(
    "Protein primary-stage AD analyzable",
    get_n(protein_analyzable_counts, "AD"),
    109
  ),
  expect_equal(
    "TMT batch levels",
    dplyr::n_distinct(prot_meta_path$tmt_batch),
    50
  ),
  expect_true(
    "Required protein covariates present",
    all(
      c(
        "age_num",
        "sex_factor",
        "pmi_num",
        "batch_factor"
      ) %in% protein_covars
    ),
    paste(protein_covars, collapse = ";")
  ),
  expect_true(
    "Protein matrix and metadata remain aligned",
    identical(colnames(prot_mat_raw), prot_meta_path$SampleID)
  ),
  expect_true(
    "Adjusted protein matrix and metadata remain aligned",
    identical(colnames(prot_mat_adj), prot_meta_adj$SampleID)
  ),
  expect_true(
    "Primary clinical stage is derived from cogdx only",
    all(
      as.character(prot_meta_path$clinical_stage) ==
        clinical_stage_primary(prot_meta_path$cogdx_num) |
        (
          is.na(prot_meta_path$clinical_stage) &
          is.na(clinical_stage_primary(prot_meta_path$cogdx_num))
        )
    )
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "core_pipeline_migration_checks.csv")
)

readr::write_csv(
  protein_primary_counts,
  file.path(out_dir, "protein_primary_stage_counts.csv")
)

readr::write_csv(
  protein_analyzable_counts,
  file.path(out_dir, "protein_primary_stage_analyzable_counts.csv")
)

git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)
git_branch <- tryCatch(
  system2("git", c("branch", "--show-current"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

writeLines(
  c(
    "Core pipeline migration checkpoint",
    paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Git branch: ", paste(git_branch, collapse = " ")),
    paste0("Git HEAD: ", paste(git_head, collapse = " ")),
    paste0("Protein covariates: ", paste(protein_covars, collapse = ";")),
    paste0("All checks passed: ", length(failures) == 0)
  ),
  file.path(out_dir, "run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("\n============================================================")
message("CORE PIPELINE MIGRATION CHECKPOINT")
message("============================================================\n")
print(checks, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))
  stop(
    "Core migration checkpoint FAILED. Do not run R/04 or later.",
    call. = FALSE
  )
}

message("\nALL CORE MIGRATION CHECKS PASSED.")
message("It is now safe to begin migrating downstream scripts 04+.")
