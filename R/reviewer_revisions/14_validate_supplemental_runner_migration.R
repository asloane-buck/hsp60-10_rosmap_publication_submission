############################################################
## 14_validate_supplemental_runner_migration.R
############################################################

options(stringsAsFactors = FALSE)

runner_path <- "R/supplemental/90_run_supplemental_pipeline.R"

if (!file.exists(runner_path)) {
  stop("Missing supplemental runner: ", runner_path, call. = FALSE)
}

runner_text <- paste(
  readLines(runner_path, warn = FALSE),
  collapse = "\n"
)

forbidden_runner_patterns <- c(
  "protein_collapse_magnitude",
  "cognition_composite_score",
  "BASE_hsp_null_tbl_built_from_regional_DLPFC",
  "sf99_scan_project_for_cognition_tables",
  "regional_base <-",
  "EmoryStrictDx\\.2019",
  "AsymAD"
)

forbidden_hits <- forbidden_runner_patterns[
  vapply(
    forbidden_runner_patterns,
    grepl,
    logical(1),
    x = runner_text
  )
]

if (length(forbidden_hits) > 0) {
  stop(
    "Corrected supplemental runner still contains stale pattern(s): ",
    paste(forbidden_hits, collapse = ", "),
    call. = FALSE
  )
}

required_runner_patterns <- c(
  "figure5_null_input_COVARIATE_ADJUSTED.csv",
  "Fig5_gene_universe_cognition_summary_by_gene.csv",
  "protein_late_decline_magnitude",
  "cognition_priority_score",
  "915",
  "306",
  "609",
  "supplemental_pipeline_sessionInfo.txt",
  "supplemental_pipeline_provenance.txt"
)

missing_required_patterns <- required_runner_patterns[
  !vapply(
    required_runner_patterns,
    grepl,
    logical(1),
    x = runner_text,
    fixed = TRUE
  )
]

if (length(missing_required_patterns) > 0) {
  stop(
    "Corrected supplemental runner is missing required invariant(s): ",
    paste(missing_required_patterns, collapse = ", "),
    call. = FALSE
  )
}

source(runner_path)

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

source("R/00_config.R")
source("R/01_utils.R")

external_key <- function(x) {
  toupper(canonical_gene_symbol(x))
}

prod <- readr::read_csv(
  file.path(
    cfg$table_dir,
    "figure5_null_input_COVARIATE_ADJUSTED.csv"
  ),
  show_col_types = FALSE
) |>
  dplyr::transmute(
    gene_key = external_key(.data$gene),
    group = dplyr::if_else(
      .data$is_hsp60_10_client,
      "Hsp60/10 clients",
      "Non-client mitochondrial proteins"
    ),
    inverse_braak = as.numeric(.data$inverse_braak_magnitude),
    late_decline = as.numeric(.data$protein_late_decline_magnitude)
  )

supp <- supfig5_outputs$rosmap_tbl |>
  dplyr::transmute(
    gene_key = external_key(.data$gene),
    group = as.character(.data$group),
    inverse_braak = as.numeric(.data$rosmap_inverse_braak),
    late_decline = as.numeric(.data$rosmap_late_decline)
  )

if (anyDuplicated(prod$gene_key) > 0 ||
    anyDuplicated(supp$gene_key) > 0) {
  stop(
    "Duplicate case-insensitive canonical keys in production/Supp Fig 5.",
    call. = FALSE
  )
}

agreement <- prod |>
  dplyr::inner_join(
    supp,
    by = "gene_key",
    suffix = c("_production", "_supplement")
  ) |>
  dplyr::mutate(
    group_agrees =
      .data$group_production == .data$group_supplement,
    braak_abs_difference = abs(
      .data$inverse_braak_production -
        .data$inverse_braak_supplement
    ),
    late_abs_difference = abs(
      .data$late_decline_production -
        .data$late_decline_supplement
    )
  )

checks <- dplyr::bind_rows(
  check_one(
    "Runner canonical null universe",
    nrow(supplemental_null_tbl),
    915,
    nrow(supplemental_null_tbl) == 915
  ),
  check_one(
    "Runner Hsp60/10 clients",
    nrow(hsp_null_tbl),
    306,
    nrow(hsp_null_tbl) == 306
  ),
  check_one(
    "Runner mitochondrial background",
    nrow(background_null_pool),
    609,
    nrow(background_null_pool) == 609
  ),
  check_one(
    "Runner Hsp cognition coverage matches joined production table",
    sum(
      supplemental_null_tbl$is_hsp60_10_client &
        is.finite(supplemental_null_tbl$cognition_priority_score),
      na.rm = TRUE
    ),
    n_hsp_cognition,
    sum(
      supplemental_null_tbl$is_hsp60_10_client &
        is.finite(supplemental_null_tbl$cognition_priority_score),
      na.rm = TRUE
    ) == n_hsp_cognition
  ),
  check_one(
    "Runner background cognition coverage matches joined production table",
    sum(
      !supplemental_null_tbl$is_hsp60_10_client &
        is.finite(supplemental_null_tbl$cognition_priority_score),
      na.rm = TRUE
    ),
    n_background_cognition,
    sum(
      !supplemental_null_tbl$is_hsp60_10_client &
        is.finite(supplemental_null_tbl$cognition_priority_score),
      na.rm = TRUE
    ) == n_background_cognition
  ),
  check_one(
    "Supp Fig 1 protein samples",
    ncol(supfig1_outputs$protein_info$matrix),
    400,
    ncol(supfig1_outputs$protein_info$matrix) == 400
  ),
  check_one(
    "Supp Fig 2 matched participants",
    nrow(supfig2_outputs$matched_stage_crosswalk),
    198,
    nrow(supfig2_outputs$matched_stage_crosswalk) == 198
  ),
  check_one(
    "Supp Fig 2 cross-modal stage agreement",
    all(supfig2_outputs$matched_stage_crosswalk$stage_agrees),
    TRUE,
    all(supfig2_outputs$matched_stage_crosswalk$stage_agrees)
  ),
  check_one(
    "Supp Fig 3 Hsp60/10 clients",
    nrow(supfig3_outputs$hsp_tbl),
    306,
    nrow(supfig3_outputs$hsp_tbl) == 306
  ),
  check_one(
    "Supp Fig 3 mitochondrial background",
    nrow(supfig3_outputs$background_tbl),
    609,
    nrow(supfig3_outputs$background_tbl) == 609
  ),
  check_one(
    "Supp Fig 4 Hsp60/10 modeled",
    sum(supfig4_outputs$plot_tbl$group == "Hsp60/10 clients"),
    306,
    sum(supfig4_outputs$plot_tbl$group == "Hsp60/10 clients") == 306
  ),
  check_one(
    "Supp Fig 4 background modeled",
    sum(
      supfig4_outputs$plot_tbl$group ==
        "Non-client mitochondrial proteins"
    ),
    609,
    sum(
      supfig4_outputs$plot_tbl$group ==
        "Non-client mitochondrial proteins"
    ) == 609
  ),
  check_one(
    "Supp Fig 5 full production agreement",
    nrow(agreement),
    915,
    nrow(agreement) == 915
  ),
  check_one(
    "Supp Fig 5 group agreement",
    all(agreement$group_agrees),
    TRUE,
    all(agreement$group_agrees)
  ),
  check_one(
    "Supp Fig 5 Braak agreement",
    max(agreement$braak_abs_difference, na.rm = TRUE),
    "<1e-12",
    max(agreement$braak_abs_difference, na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Supp Fig 5 late-decline agreement",
    max(agreement$late_abs_difference, na.rm = TRUE),
    "<1e-12",
    max(agreement$late_abs_difference, na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Supp Fig 6 difference metrics",
    sort(unique(as.character(supfig6_outputs$diff_tbl$metric))),
    sort(c("AD-associated decline", "High-Braak coupling")),
    identical(
      sort(unique(as.character(supfig6_outputs$diff_tbl$metric))),
      sort(c("AD-associated decline", "High-Braak coupling"))
    )
  ),
  check_one(
    "Six manuscript-ready supplemental PDFs",
    nrow(manuscript_ready_pdf_manifest),
    6,
    nrow(manuscript_ready_pdf_manifest) == 6 &&
      all(file.exists(manuscript_ready_pdf_manifest$pdf_path))
  ),
  check_one(
    "Runner provenance file written",
    file.exists(
      file.path(
        audits_dir,
        "supplemental_pipeline_provenance.txt"
      )
    ),
    TRUE,
    file.exists(
      file.path(
        audits_dir,
        "supplemental_pipeline_provenance.txt"
      )
    )
  ),
  check_one(
    "Runner sessionInfo file written",
    file.exists(
      file.path(
        audits_dir,
        "supplemental_pipeline_sessionInfo.txt"
      )
    ),
    TRUE,
    file.exists(
      file.path(
        audits_dir,
        "supplemental_pipeline_sessionInfo.txt"
      )
    )
  )
)

validation_dir <- file.path(
  project_dir,
  "outputs",
  "reviewer_revisions",
  "supplemental_runner_migration"
)

dir.create(
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  checks,
  file.path(
    validation_dir,
    "supplemental_runner_migration_checks.csv"
  )
)

readr::write_csv(
  agreement,
  file.path(
    validation_dir,
    "suppfig5_runner_vs_production_agreement.csv"
  )
)

message("\n============================================================")
message("SUPPLEMENTAL RUNNER MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))
  stop(
    "Supplemental runner migration FAILED.",
    call. = FALSE
  )
}

message("\nALL SUPPLEMENTAL RUNNER MIGRATION CHECKS PASSED.")
