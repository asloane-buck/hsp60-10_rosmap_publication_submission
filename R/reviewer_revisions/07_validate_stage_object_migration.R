############################################################
## 07_validate_stage_object_migration.R
##
## Sources ONLY production scripts R/00 through R/05.
## Validates canonical NCI/MCI/AD stage propagation before
## cognition/null/figure scripts are migrated.
############################################################

options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")
source("R/02_load_data.R")
source("R/03_build_adjusted_core_objects.R")
source("R/04_build_pathway_sets_all_clients.R")
source("R/05_build_adjusted_all_client_tables.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "stage_object_migration"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

failures <- character()

record_check <- function(label, observed, expected, ok) {
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

count_stage <- function(df, stage_col, score_col = NULL) {
  x <- df

  if (!is.null(score_col)) {
    x <- x |>
      dplyr::filter(is.finite(.data[[score_col]]))
  }

  x |>
    dplyr::filter(!is.na(.data[[stage_col]])) |>
    dplyr::count(stage = .data[[stage_col]], name = "n")
}

get_stage_n <- function(tbl, stage_name) {
  hit <- tbl$n[as.character(tbl$stage) == stage_name]
  if (length(hit) == 0) 0L else as.integer(hit[[1]])
}

protein_score_counts <- count_stage(
  prot_scores,
  "clinical_stage",
  "Hsp60_10_all_clients"
)

rna_score_counts <- count_stage(
  rna_scores,
  "clinical_stage",
  "Hsp60_10_all_clients"
)

required_new_columns <- c(
  "protein_nci_mean",
  "protein_mci_mean",
  "protein_ad_mean",
  "protein_early_effect",
  "protein_late_effect",
  "protein_total_effect",
  "protein_late_decline_magnitude",
  "late_decline_percentile",
  "rna_nci_mean",
  "rna_mci_mean",
  "rna_ad_mean"
)

forbidden_old_columns <- c(
  "protein_control_mean",
  "protein_early_mean",
  "protein_collapse_magnitude",
  "collapse_percentile",
  "rna_control_mean",
  "rna_early_mean"
)

protein_client_rows <- all_hsp60_10_client_tbl |>
  dplyr::filter(.data$detected_in_protein)

trajectory <- prot_scores |>
  dplyr::filter(
    !is.na(.data$clinical_stage),
    is.finite(.data$Hsp60_10_all_clients)
  ) |>
  dplyr::group_by(.data$clinical_stage, .drop = FALSE) |>
  dplyr::summarise(
    mean_score = mean(.data$Hsp60_10_all_clients),
    n = dplyr::n(),
    .groups = "drop"
  )

late_fraction_negative <- mean(
  protein_client_rows$protein_late_effect < 0,
  na.rm = TRUE
)

mean_late_decline <- mean(
  protein_client_rows$protein_late_decline_magnitude,
  na.rm = TRUE
)

script04 <- paste(
  readLines("R/04_build_pathway_sets_all_clients.R", warn = FALSE),
  collapse = "\n"
)

script05 <- paste(
  readLines("R/05_build_adjusted_all_client_tables.R", warn = FALSE),
  collapse = "\n"
)

invalid_mapping_regex <- paste0(
  "AsymAD.{0,80}(MCI|Early AD)|",
  "(MCI|Early AD).{0,80}AsymAD"
)

checks <- dplyr::bind_rows(
  record_check(
    "Protein pathway-score rows",
    nrow(prot_scores),
    400,
    nrow(prot_scores) == 400
  ),
  record_check(
    "RNA pathway-score rows",
    nrow(rna_scores),
    578,
    nrow(rna_scores) == 578
  ),
  record_check(
    "Protein Hsp score NCI finite n",
    get_stage_n(protein_score_counts, "NCI"),
    167,
    get_stage_n(protein_score_counts, "NCI") == 167
  ),
  record_check(
    "Protein Hsp score MCI finite n",
    get_stage_n(protein_score_counts, "MCI"),
    96,
    get_stage_n(protein_score_counts, "MCI") == 96
  ),
  record_check(
    "Protein Hsp score AD finite n",
    get_stage_n(protein_score_counts, "AD"),
    109,
    get_stage_n(protein_score_counts, "AD") == 109
  ),
  record_check(
    "RNA Hsp score NCI finite n",
    get_stage_n(rna_score_counts, "NCI"),
    200,
    get_stage_n(rna_score_counts, "NCI") == 200
  ),
  record_check(
    "RNA Hsp score MCI finite n",
    get_stage_n(rna_score_counts, "MCI"),
    158,
    get_stage_n(rna_score_counts, "MCI") == 158
  ),
  record_check(
    "RNA Hsp score AD finite n",
    get_stage_n(rna_score_counts, "AD"),
    220,
    get_stage_n(rna_score_counts, "AD") == 220
  ),
  record_check(
    "Reference Hsp60/10 client rows",
    nrow(all_hsp60_10_client_tbl),
    321,
    nrow(all_hsp60_10_client_tbl) == 321
  ),
  record_check(
    "Detected Hsp60/10 proteins",
    sum(all_hsp60_10_client_tbl$detected_in_protein),
    306,
    sum(all_hsp60_10_client_tbl$detected_in_protein) == 306
  ),
  record_check(
    "Detected Hsp60/10 RNA genes",
    sum(all_hsp60_10_client_tbl$detected_in_rna),
    294,
    sum(all_hsp60_10_client_tbl$detected_in_rna) == 294
  ),
  record_check(
    "Required canonical stage columns present",
    all(required_new_columns %in% colnames(all_hsp60_10_client_tbl)),
    TRUE,
    all(required_new_columns %in% colnames(all_hsp60_10_client_tbl))
  ),
  record_check(
    "Old stage/collapse columns absent",
    any(forbidden_old_columns %in% colnames(all_hsp60_10_client_tbl)),
    FALSE,
    !any(forbidden_old_columns %in% colnames(all_hsp60_10_client_tbl))
  ),
  record_check(
    "Long table stages",
    sort(unique(as.character(all_client_stage_long$Stage))),
    c("AD", "MCI", "NCI"),
    identical(
      sort(unique(as.character(all_client_stage_long$Stage))),
      c("AD", "MCI", "NCI")
    )
  ),
  record_check(
    "R/04 invalid AsymAD mapping absent",
    stringr::str_detect(script04, invalid_mapping_regex),
    FALSE,
    !stringr::str_detect(script04, invalid_mapping_regex)
  ),
  record_check(
    "R/05 invalid AsymAD mapping absent",
    stringr::str_detect(script05, invalid_mapping_regex),
    FALSE,
    !stringr::str_detect(script05, invalid_mapping_regex)
  ),
  record_check(
    "Late-shift negative fraction approximately audit-05 value",
    round(late_fraction_negative, 4),
    "0.794 +/- 0.01",
    abs(late_fraction_negative - 0.794) <= 0.01
  ),
  record_check(
    "Mean late-decline magnitude approximately audit-05 value",
    round(mean_late_decline, 4),
    "0.0230 +/- 0.005",
    abs(mean_late_decline - 0.0230) <= 0.005
  )
)

expected_trajectory <- tibble::tibble(
  clinical_stage = factor(
    c("NCI", "MCI", "AD"),
    levels = c("NCI", "MCI", "AD")
  ),
  expected = c(0.0934, 0.0370, -0.132)
)

trajectory_check <- trajectory |>
  dplyr::left_join(
    expected_trajectory,
    by = "clinical_stage"
  ) |>
  dplyr::mutate(
    abs_difference = abs(.data$mean_score - .data$expected),
    passed = .data$abs_difference <= 0.01
  )

if (any(!trajectory_check$passed)) {
  failures <- c(
    failures,
    "Protein Hsp60/10 pathway trajectory differs materially from audit 05."
  )
}

readr::write_csv(
  checks,
  file.path(out_dir, "stage_object_migration_checks.csv")
)

readr::write_csv(
  trajectory,
  file.path(out_dir, "protein_hsp60_stage_trajectory.csv")
)

readr::write_csv(
  trajectory_check,
  file.path(out_dir, "protein_hsp60_stage_trajectory_check.csv")
)

readr::write_csv(
  stage_analysis_sample_counts,
  file.path(out_dir, "stage_analysis_sample_counts.csv")
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
    "Stage-object migration checkpoint",
    paste0("Run time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Git branch: ", paste(git_branch, collapse = " ")),
    paste0("Git HEAD: ", paste(git_head, collapse = " ")),
    paste0("All checks passed: ", length(failures) == 0)
  ),
  file.path(out_dir, "run_provenance.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("\n============================================================")
message("STAGE-OBJECT MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nPROTEIN Hsp60/10 PATHWAY TRAJECTORY")
print(trajectory_check, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))

  stop(
    "Stage-object migration FAILED. Do not run R/06 or later.",
    call. = FALSE
  )
}

message("\nALL STAGE-OBJECT MIGRATION CHECKS PASSED.")
message("R/04 and R/05 now use canonical cogdx-derived clinical stages.")
