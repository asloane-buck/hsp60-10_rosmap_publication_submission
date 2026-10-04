############################################################
## 08_validate_cognition_null_migration.R
##
## Sources production scripts R/00 through R/07 and checks that
## cognition/null objects use the corrected full TMT cohort and
## canonical clinical-stage definitions.
############################################################

options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")
source("R/02_load_data.R")
source("R/03_build_adjusted_core_objects.R")
source("R/04_build_pathway_sets_all_clients.R")
source("R/05_build_adjusted_all_client_tables.R")
source("R/06_build_cognition_objects.R")
source("R/07_build_matched_null_objects.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "cognition_null_migration"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

score_compare <- cognition_model_df |>
  dplyr::select(
    SampleID,
    cognition_score = Hsp60_client_score
  ) |>
  dplyr::inner_join(
    prot_scores |>
      dplyr::select(
        SampleID,
        pathway_score = Hsp60_10_all_clients
      ),
    by = "SampleID"
  ) |>
  dplyr::mutate(
    abs_difference = abs(
      .data$cognition_score - .data$pathway_score
    )
  )

hsp_null_compare <- hsp_null_tbl |>
  dplyr::select(
    gene,
    null_late_effect = protein_late_effect,
    null_late_decline = protein_late_decline_magnitude
  ) |>
  dplyr::inner_join(
    all_hsp60_10_client_tbl |>
      dplyr::filter(.data$detected_in_protein) |>
      dplyr::select(
        gene,
        source_late_effect = protein_late_effect,
        source_late_decline = protein_late_decline_magnitude
      ),
    by = "gene"
  ) |>
  dplyr::mutate(
    effect_difference = abs(
      .data$null_late_effect - .data$source_late_effect
    ),
    decline_difference = abs(
      .data$null_late_decline - .data$source_late_decline
    )
  )

script06 <- paste(
  readLines("R/06_build_cognition_objects.R", warn = FALSE),
  collapse = "\n"
)

script07 <- paste(
  readLines("R/07_build_matched_null_objects.R", warn = FALSE),
  collapse = "\n"
)

checks <- dplyr::bind_rows(
  check_one(
    "Cognition model rows",
    nrow(cognition_model_df),
    400,
    nrow(cognition_model_df) == 400
  ),
  check_one(
    "Unique cognition SampleID",
    dplyr::n_distinct(cognition_model_df$SampleID),
    400,
    dplyr::n_distinct(cognition_model_df$SampleID) == 400
  ),
  check_one(
    "Unique cognition individualID",
    dplyr::n_distinct(cognition_model_df$individualID),
    400,
    dplyr::n_distinct(cognition_model_df$individualID) == 400
  ),
  check_one(
    "Finite Hsp60 cognition scores",
    sum(is.finite(cognition_model_df$Hsp60_client_score)),
    398,
    sum(is.finite(cognition_model_df$Hsp60_client_score)) == 398
  ),
  check_one(
    "Cognition cogdx available",
    sum(is.finite(cognition_model_df$cogdx)),
    400,
    sum(is.finite(cognition_model_df$cogdx)) == 400
  ),
  check_one(
    "Cognition score exactly equals prot_scores",
    max(score_compare$abs_difference, na.rm = TRUE),
    "<1e-12",
    max(score_compare$abs_difference, na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Matched-null Hsp clients",
    nrow(hsp_null_tbl),
    306,
    nrow(hsp_null_tbl) == 306
  ),
  check_one(
    "Matched-null non-Hsp mitochondrial background",
    nrow(background_null_pool),
    609,
    nrow(background_null_pool) == 609
  ),
  check_one(
    "Matched-null Hsp late effects reproduce source table",
    max(hsp_null_compare$effect_difference, na.rm = TRUE),
    "<1e-10",
    max(hsp_null_compare$effect_difference, na.rm = TRUE) < 1e-10
  ),
  check_one(
    "Matched-null Hsp late-decline magnitudes reproduce source table",
    max(hsp_null_compare$decline_difference, na.rm = TRUE),
    "<1e-10",
    max(hsp_null_compare$decline_difference, na.rm = TRUE) < 1e-10
  ),
  check_one(
    "Observed matched-null mean late decline",
    round(observed_stats$observed_mean_late_decline, 4),
    "0.0230 +/- 0.005",
    abs(observed_stats$observed_mean_late_decline - 0.0230) <= 0.005
  ),
  check_one(
    "Null iterations",
    nrow(null_results),
    cfg$null_n_iter,
    nrow(null_results) == cfg$null_n_iter
  ),
  check_one(
    "Null summary uses late-decline terminology",
    "protein_late_decline_magnitude" %in% null_summary$metric,
    TRUE,
    "protein_late_decline_magnitude" %in% null_summary$metric
  ),
  check_one(
    "R/06 has no invalid legacy-stage mapping",
    stringr::str_detect(
      script06,
      "AsymAD.{0,80}(MCI|Early AD)|(MCI|Early AD).{0,80}AsymAD"
    ),
    FALSE,
    !stringr::str_detect(
      script06,
      "AsymAD.{0,80}(MCI|Early AD)|(MCI|Early AD).{0,80}AsymAD"
    )
  ),
  check_one(
    "R/07 has no invalid legacy-stage mapping",
    stringr::str_detect(
      script07,
      "AsymAD.{0,80}(MCI|Early AD)|(MCI|Early AD).{0,80}AsymAD"
    ),
    FALSE,
    !stringr::str_detect(
      script07,
      "AsymAD.{0,80}(MCI|Early AD)|(MCI|Early AD).{0,80}AsymAD"
    )
  ),
  check_one(
    "R/07 has no legacy collapse identifiers",
    stringr::str_detect(
      script07,
      "protein_collapse|collapse_rank|null_mean_collapse|observed_mean_collapse"
    ),
    FALSE,
    !stringr::str_detect(
      script07,
      "protein_collapse|collapse_rank|null_mean_collapse|observed_mean_collapse"
    )
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "cognition_null_migration_checks.csv")
)

readr::write_csv(
  cognition_availability,
  file.path(out_dir, "cognition_variable_availability.csv")
)

readr::write_csv(
  null_summary,
  file.path(out_dir, "corrected_null_summary.csv")
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
    "Cognition/null migration checkpoint",
    paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Git branch: ", paste(git_branch, collapse = " ")),
    paste0("Git HEAD: ", paste(git_head, collapse = " ")),
    paste0("Null iterations: ", cfg$null_n_iter),
    paste0("All checks passed: ", length(failures) == 0)
  ),
  file.path(out_dir, "run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("\n============================================================")
message("COGNITION / MATCHED-NULL MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nCOGNITION VARIABLE AVAILABILITY")
print(cognition_availability, n = Inf)

message("\nCORRECTED MATCHED-NULL SUMMARY")
print(null_summary, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))

  stop(
    "Cognition/null migration FAILED. Do not run figure scripts.",
    call. = FALSE
  )
}

message("\nALL COGNITION / MATCHED-NULL MIGRATION CHECKS PASSED.")
message("Production object builders R/00-R/07 are now internally consistent.")
