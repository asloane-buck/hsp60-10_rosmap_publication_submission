############################################################
## 13_validate_supplemental_figures_5_6_migration.R
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

production_hsp_null_tbl <- hsp_null_tbl
production_background_null_pool <- background_null_pool

source("R/supplemental/00_supplemental_config.R")
source("R/supplemental/01_supplemental_load_inputs.R")
source("R/supplemental/02_supplemental_helper_functions.R")

## Explicitly inject corrected ROSMAP discovery objects.
inputs$hsp_null_tbl <- production_hsp_null_tbl
inputs$background_null_pool <- production_background_null_pool

source("R/supplemental/14_make_supplementary_figure_5_msbb_cross_cohort_validation.R")
source("R/supplemental/15_make_supplementary_figure_6_regional_proteomics_validation.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "supplemental_figures_5_6_migration"
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

############################################################
## Supp Fig 5 agreement with corrected ROSMAP source objects
############################################################

production_rosmap <- dplyr::bind_rows(
  production_hsp_null_tbl |>
    dplyr::transmute(
      gene,
      production_group = "Hsp60/10 clients",
      production_inverse_braak = inverse_braak_magnitude,
      production_late_decline = protein_late_decline_magnitude
    ),
  production_background_null_pool |>
    dplyr::transmute(
      gene,
      production_group = "Non-client mitochondrial proteins",
      production_inverse_braak = inverse_braak_magnitude,
      production_late_decline = protein_late_decline_magnitude
    )
) |>
  dplyr::distinct(.data$gene, .keep_all = TRUE)

sf5_agreement <- supfig5_outputs$rosmap_tbl |>
  dplyr::mutate(group_chr = as.character(.data$group)) |>
  dplyr::inner_join(production_rosmap, by = "gene") |>
  dplyr::mutate(
    group_agrees = .data$group_chr == .data$production_group,
    braak_abs_difference = abs(
      .data$rosmap_inverse_braak - .data$production_inverse_braak
    ),
    decline_abs_difference = abs(
      .data$rosmap_late_decline - .data$production_late_decline
    )
  )

script14 <- paste(
  readLines(
    "R/supplemental/14_make_supplementary_figure_5_msbb_cross_cohort_validation.R",
    warn = FALSE
  ),
  collapse = "\n"
)

script15 <- paste(
  readLines(
    "R/supplemental/15_make_supplementary_figure_6_regional_proteomics_validation.R",
    warn = FALSE
  ),
  collapse = "\n"
)

sf6_metric_levels <- levels(supfig6_outputs$diff_tbl$metric)

checks <- dplyr::bind_rows(
  check_one(
    "Supp Fig 5 ROSMAP discovery genes",
    nrow(supfig5_outputs$rosmap_tbl),
    915,
    nrow(supfig5_outputs$rosmap_tbl) == 915
  ),
  check_one(
    "Supp Fig 5 ROSMAP Hsp60/10 clients",
    sum(supfig5_outputs$rosmap_tbl$group == "Hsp60/10 clients"),
    306,
    sum(supfig5_outputs$rosmap_tbl$group == "Hsp60/10 clients") == 306
  ),
  check_one(
    "Supp Fig 5 ROSMAP mitochondrial background",
    sum(
      supfig5_outputs$rosmap_tbl$group ==
        "Non-client mitochondrial proteins"
    ),
    609,
    sum(
      supfig5_outputs$rosmap_tbl$group ==
        "Non-client mitochondrial proteins"
    ) == 609
  ),
  check_one(
    "Supp Fig 5 ROSMAP genes equal production universe",
    nrow(sf5_agreement),
    915,
    nrow(sf5_agreement) == 915 &&
      setequal(
        supfig5_outputs$rosmap_tbl$gene,
        production_rosmap$gene
      )
  ),
  check_one(
    "Supp Fig 5 ROSMAP group assignments reproduce production",
    all(sf5_agreement$group_agrees),
    TRUE,
    all(sf5_agreement$group_agrees)
  ),
  check_one(
    "Supp Fig 5 inverse Braak reproduces production",
    max(sf5_agreement$braak_abs_difference, na.rm = TRUE),
    "<1e-12",
    max(sf5_agreement$braak_abs_difference, na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Supp Fig 5 late decline reproduces production",
    max(sf5_agreement$decline_abs_difference, na.rm = TRUE),
    "<1e-12",
    max(sf5_agreement$decline_abs_difference, na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Supp Fig 5 MSBB matched genes nonempty",
    nrow(supfig5_outputs$msbb_tbl),
    ">0",
    nrow(supfig5_outputs$msbb_tbl) > 0
  ),
  check_one(
    "Supp Fig 5 cross-cohort Hsp overlap adequate",
    sum(supfig5_outputs$cross_tbl$group == "Hsp60/10 clients"),
    ">=10",
    sum(supfig5_outputs$cross_tbl$group == "Hsp60/10 clients") >= 10
  ),
  check_one(
    "Supp Fig 5 old ROSMAP collapse identifier absent",
    grepl("protein_collapse_magnitude|rosmap_collapse", script14),
    FALSE,
    !grepl("protein_collapse_magnitude|rosmap_collapse", script14)
  ),
  check_one(
    "Supp Fig 6 regional AD table nonempty",
    nrow(supfig6_outputs$ad_tbl),
    ">0",
    nrow(supfig6_outputs$ad_tbl) > 0
  ),
  check_one(
    "Supp Fig 6 regional Braak table nonempty",
    nrow(supfig6_outputs$braak_tbl),
    ">0",
    nrow(supfig6_outputs$braak_tbl) > 0
  ),
  check_one(
    "Supp Fig 6 regional difference metrics",
    sf6_metric_levels,
    c("AD-associated decline", "High-Braak coupling"),
    identical(
      sf6_metric_levels,
      c("AD-associated decline", "High-Braak coupling")
    )
  ),
  check_one(
    "Supp Fig 6 displayed collapse terminology absent",
    grepl(
      'AD-associated collapse|ylab = "AD decline\\\\n-\\\\(AD - Control\\\\)"',
      script15
    ),
    FALSE,
    !grepl(
      'AD-associated collapse|ylab = "AD decline\\\\n-\\\\(AD - Control\\\\)"',
      script15
    )
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "supplemental_figures_5_6_migration_checks.csv")
)

readr::write_csv(
  sf5_agreement,
  file.path(out_dir, "suppfig5_rosmap_vs_production_agreement.csv")
)

readr::write_csv(
  supfig5_outputs$rosmap_summary,
  file.path(out_dir, "suppfig5_rosmap_summary.csv")
)

readr::write_csv(
  supfig5_outputs$msbb_summary,
  file.path(out_dir, "suppfig5_msbb_summary.csv")
)

readr::write_csv(
  supfig5_outputs$concordance_summary,
  file.path(out_dir, "suppfig5_concordance_summary.csv")
)

readr::write_csv(
  supfig6_outputs$ad_tbl,
  file.path(out_dir, "suppfig6_regional_ad_values.csv")
)

readr::write_csv(
  supfig6_outputs$braak_tbl,
  file.path(out_dir, "suppfig6_regional_braak_values.csv")
)

readr::write_csv(
  supfig6_outputs$diff_tbl,
  file.path(out_dir, "suppfig6_regional_difference_values.csv")
)

message("\n============================================================")
message("SUPPLEMENTAL FIGURES 5-6 MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nSUPP FIG 5 ROSMAP SUMMARY")
print(supfig5_outputs$rosmap_summary, n = Inf)

message("\nSUPP FIG 5 MSBB SUMMARY")
print(supfig5_outputs$msbb_summary, n = Inf)

message("\nSUPP FIG 5 CROSS-COHORT CONCORDANCE")
print(supfig5_outputs$concordance_summary, n = Inf)

message("\nSUPP FIG 5 VALIDATION TESTS")
print(supfig5_outputs$validation_test_summary, n = Inf)

message("\nSUPP FIG 6 REGIONAL DIFFERENCE SUMMARY")
print(
  supfig6_outputs$diff_tbl |>
    dplyr::group_by(.data$metric) |>
    dplyr::summarise(
      n = dplyr::n(),
      median_stg_minus_dlpfc = stats::median(
        .data$stg_minus_dlpfc,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),
  n = Inf
)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))

  stop(
    "Supplemental Figures 5-6 migration FAILED. ",
    "Do not migrate the full supplemental runner yet.",
    call. = FALSE
  )
}

message("\nALL SUPPLEMENTAL FIGURE 5-6 MIGRATION CHECKS PASSED.")
