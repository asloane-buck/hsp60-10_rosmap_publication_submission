############################################################
## 15_validate_sensitivity_migration.R
##
## Fast smoke validation for migrated R/92.
## Uses 100 permutations and no plots; this is NOT the final archival run.
############################################################

options(stringsAsFactors = FALSE)

script_path <- "R/92_pre_submission_sensitivity_checks_PIPELINE_OBJECTS_FAST.R"

if (!file.exists(script_path)) {
  stop("Missing migrated R/92 sensitivity script.", call. = FALSE)
}

script_text <- paste(
  readLines(script_path, warn = FALSE),
  collapse = "\n"
)

forbidden_patterns <- c(
  "EmoryStrictDx\\.2019",
  "AsymAD",
  "protein_collapse_magnitude",
  "collapse_percentile",
  "SENSITIVITY_TOP_N_COLLAPSE",
  "top_late_collapse_clients",
  "collapse_central_clients"
)

forbidden_hits <- forbidden_patterns[
  vapply(
    forbidden_patterns,
    grepl,
    logical(1),
    x = script_text
  )
]

if (length(forbidden_hits) > 0) {
  stop(
    "Migrated R/92 still contains stale identifier(s): ",
    paste(forbidden_hits, collapse = ", "),
    call. = FALSE
  )
}

required_patterns <- c(
  "clinical_stage",
  'levels = c("NCI", "MCI", "AD")',
  "protein_late_decline_magnitude",
  "late_decline_percentile",
  "SENSITIVITY_TOP_N_LATE_DECLINE",
  "sensitivity_provenance.txt",
  "sensitivity_sessionInfo.txt"
)

missing_required <- required_patterns[
  !vapply(
    required_patterns,
    grepl,
    logical(1),
    x = script_text,
    fixed = TRUE
  )
]

if (length(missing_required) > 0) {
  stop(
    "Migrated R/92 is missing corrected identifier(s): ",
    paste(missing_required, collapse = ", "),
    call. = FALSE
  )
}

source("R/00_config.R")

smoke_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "sensitivity_migration_smoke"
)

dir.create(
  smoke_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

Sys.setenv(
  SENSITIVITY_OUTPUT_DIR = smoke_dir,
  SENSITIVITY_N_PERM = "100",
  SENSITIVITY_N_PERM_CHECK1 = "100",
  SENSITIVITY_N_PERM_CHECK2 = "100",
  SENSITIVITY_WITHIN_STAGE_NULLS = "FALSE",
  SENSITIVITY_SAVE_PLOTS = "FALSE",
  SENSITIVITY_WRITE_FULL_NULL_DISTRIBUTIONS = "FALSE",
  SENSITIVITY_SOURCE_UPSTREAM_IF_MISSING = "TRUE"
)

source(script_path, local = .GlobalEnv)

failures <- character()

check_one <- function(label, observed, expected, ok) {
  if (!isTRUE(ok)) {
    failures <<- c(
      failures,
      paste0(
        label,
        ": observed=",
        paste(observed, collapse = ";"),
        " expected=",
        paste(expected, collapse = ";")
      )
    )
  }

  tibble::tibble(
    check = label,
    observed = paste(observed, collapse = ";"),
    expected = paste(expected, collapse = ";"),
    passed = isTRUE(ok)
  )
}

stage_counts <- make_protein_stage_meta() |>
  dplyr::count(.data$stage, name = "n") |>
  dplyr::mutate(stage = as.character(.data$stage))

get_stage_n <- function(stage_name) {
  x <- stage_counts$n[stage_counts$stage == stage_name]
  if (length(x) == 0) 0L else as.integer(x[[1]])
}

check1_metric_names <- results$check1$residual_centrality_tests$predictor
check2_metric_names <- unique(results$check2$overall_matched$metric)

checks <- dplyr::bind_rows(
  check_one(
    "Sensitivity client proteins",
    nrow(client_metrics),
    306,
    nrow(client_metrics) == 306
  ),
  check_one(
    "Primary stage NCI",
    get_stage_n("NCI"),
    168,
    get_stage_n("NCI") == 168
  ),
  check_one(
    "Primary stage MCI",
    get_stage_n("MCI"),
    97,
    get_stage_n("MCI") == 97
  ),
  check_one(
    "Primary stage AD",
    get_stage_n("AD"),
    109,
    get_stage_n("AD") == 109
  ),
  check_one(
    "Check 1 completed",
    !is.null(results$check1),
    TRUE,
    !is.null(results$check1)
  ),
  check_one(
    "Check 2 completed",
    !is.null(results$check2),
    TRUE,
    !is.null(results$check2)
  ),
  check_one(
    "Check 1 late-decline magnitude predictor",
    "protein_late_decline_magnitude" %in% check1_metric_names,
    TRUE,
    "protein_late_decline_magnitude" %in% check1_metric_names
  ),
  check_one(
    "Check 1 late-decline percentile predictor",
    "late_decline_percentile" %in% check1_metric_names,
    TRUE,
    "late_decline_percentile" %in% check1_metric_names
  ),
  check_one(
    "Check 2 late-decline metric",
    "protein_late_decline_magnitude" %in% check2_metric_names,
    TRUE,
    "protein_late_decline_magnitude" %in% check2_metric_names
  ),
  check_one(
    "Check 2 inverse-Braak metric",
    "inverse_braak_magnitude" %in% check2_metric_names,
    TRUE,
    "inverse_braak_magnitude" %in% check2_metric_names
  ),
  check_one(
    "Check 2 pathology score metric",
    "pathology_vulnerability_score" %in% check2_metric_names,
    TRUE,
    "pathology_vulnerability_score" %in% check2_metric_names
  ),
  check_one(
    "Smoke Check 1 permutations",
    unique(results$check1$residual_summary$n_perm),
    100,
    all(results$check1$residual_summary$n_perm == 100)
  ),
  check_one(
    "Smoke Check 2 permutations",
    unique(results$check2$overall_matched$n_perm),
    100,
    all(results$check2$overall_matched$n_perm == 100)
  ),
  check_one(
    "Stage audit written",
    file.exists(
      file.path(
        smoke_dir,
        "check1_primary_clinical_stage_sample_counts.csv"
      )
    ),
    TRUE,
    file.exists(
      file.path(
        smoke_dir,
        "check1_primary_clinical_stage_sample_counts.csv"
      )
    )
  ),
  check_one(
    "Sensitivity provenance written",
    file.exists(
      file.path(
        smoke_dir,
        "sensitivity_provenance.txt"
      )
    ),
    TRUE,
    file.exists(
      file.path(
        smoke_dir,
        "sensitivity_provenance.txt"
      )
    )
  ),
  check_one(
    "Sensitivity sessionInfo written",
    file.exists(
      file.path(
        smoke_dir,
        "sensitivity_sessionInfo.txt"
      )
    ),
    TRUE,
    file.exists(
      file.path(
        smoke_dir,
        "sensitivity_sessionInfo.txt"
      )
    )
  )
)

readr::write_csv(
  checks,
  file.path(
    smoke_dir,
    "sensitivity_migration_smoke_checks.csv"
  )
)

message("\n============================================================")
message("R/92 SENSITIVITY MIGRATION SMOKE CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nCHECK 1 STAGE-RESIDUALIZED SUMMARY")
print(results$check1$residual_summary, n = Inf)

message("\nCHECK 1 WITHIN-STAGE SUMMARY")
print(results$check1$within_stage_summary, n = Inf)

message("\nCHECK 2 OVERALL FUNCTION/ABUNDANCE-MATCHED SUMMARY")
print(results$check2$overall_matched, n = Inf)

message("\nCHECK 2 FUNCTION-ADJUSTED REGRESSION SUMMARY")
print(results$check2$regression_tests, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))
  stop(
    "R/92 sensitivity migration smoke validation FAILED.",
    call. = FALSE
  )
}

message("\nALL R/92 SENSITIVITY MIGRATION SMOKE CHECKS PASSED.")
message(
  "This used 100 permutations for structural validation only. ",
  "Do not use these smoke-run p-values in the manuscript."
)
