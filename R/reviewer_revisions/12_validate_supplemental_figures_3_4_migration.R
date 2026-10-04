############################################################
## 12_validate_supplemental_figures_3_4_migration.R
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

## Main Figure 5 computes cognition_priority_score for the complete
## Hsp60/10 + non-client mitochondrial gene universe.
source("R/14_make_main_figure_5_matched_null_specificity.R")

main_fig5_null_summary <- null_summary
main_fig5_hsp_null_tbl <- hsp_null_tbl
main_fig5_background_null_pool <- background_null_pool

source("R/supplemental/00_supplemental_config.R")
source("R/supplemental/01_supplemental_load_inputs.R")
source("R/supplemental/02_supplemental_helper_functions.R")

## Explicitly provide corrected production objects. This intentionally
## bypasses the old supplemental runner's regional-screen reconstruction.
inputs$hsp_null_tbl <- main_fig5_hsp_null_tbl
inputs$background_null_pool <- main_fig5_background_null_pool
inputs$prot_mat_raw <- prot_mat_raw
inputs$prot_mat <- prot_mat
inputs$prot_meta <- prot_meta_adj
inputs$protein_meta <- prot_meta_adj
inputs$prot_meta_aligned <- prot_meta_adj

source("R/supplemental/12_make_supplementary_figure_3_mitochondrial_specificity.R")
source("R/supplemental/13_make_supplementary_figure_4_pathology_model_robustness.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "supplemental_figures_3_4_migration"
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

get_metric_value <- function(tbl, metric_name, value_col) {
  hit <- tbl |>
    dplyr::filter(.data$metric == metric_name) |>
    dplyr::pull(dplyr::all_of(value_col))

  if (length(hit) != 1) NA_real_ else as.numeric(hit[[1]])
}

sf3_summary <- supfig3_outputs$summary_tbl

main_late_obs <- get_metric_value(
  main_fig5_null_summary,
  "protein_late_decline_magnitude",
  "observed"
)
main_braak_obs <- get_metric_value(
  main_fig5_null_summary,
  "inverse_braak_magnitude",
  "observed"
)
main_cog_obs <- get_metric_value(
  main_fig5_null_summary,
  "cognition_priority_score",
  "observed"
)
main_agora_obs <- get_metric_value(
  main_fig5_null_summary,
  "agora_fraction",
  "observed"
)

sf3_late_obs <- get_metric_value(
  sf3_summary,
  "protein_late_decline_magnitude",
  "observed_mean"
)
sf3_braak_obs <- get_metric_value(
  sf3_summary,
  "inverse_braak_magnitude",
  "observed_mean"
)
sf3_cog_obs <- get_metric_value(
  sf3_summary,
  "cognition_priority_score",
  "observed_mean"
)
sf3_agora_obs <- get_metric_value(
  sf3_summary,
  "agora_nominated_target",
  "observed_mean"
)

production_braak <- dplyr::bind_rows(
  main_fig5_hsp_null_tbl |>
    dplyr::select(gene, production_inverse_braak = inverse_braak_magnitude),
  main_fig5_background_null_pool |>
    dplyr::select(gene, production_inverse_braak = inverse_braak_magnitude)
) |>
  dplyr::distinct(.data$gene, .keep_all = TRUE)

sf4_compare <- supfig4_outputs$adjusted_tbl |>
  dplyr::select(
    gene,
    supplemental_inverse_braak = inverse_braak_beta_adjusted
  ) |>
  dplyr::inner_join(production_braak, by = "gene") |>
  dplyr::mutate(
    abs_difference = abs(
      .data$supplemental_inverse_braak -
        .data$production_inverse_braak
    )
  )

sf4_metadata_audit <- readr::read_csv(
  file.path(audits_dir, "SuppFig4_metadata_column_audit.csv"),
  show_col_types = FALSE
)

batch_detected <- sf4_metadata_audit |>
  dplyr::filter(.data$conceptual_variable == "tmt_batch") |>
  dplyr::pull(.data$detected_column)

script12 <- paste(
  readLines(
    "R/supplemental/12_make_supplementary_figure_3_mitochondrial_specificity.R",
    warn = FALSE
  ),
  collapse = "\n"
)

script13 <- paste(
  readLines(
    "R/supplemental/13_make_supplementary_figure_4_pathology_model_robustness.R",
    warn = FALSE
  ),
  collapse = "\n"
)

checks <- dplyr::bind_rows(
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
    "Supp Fig 3 null permutations per metric",
    paste(sort(unique(sf3_summary$n_perm)), collapse = ";"),
    10000,
    all(sf3_summary$n_perm == 10000)
  ),
  check_one(
    "Supp Fig 3 metric set",
    sort(as.character(sf3_summary$metric)),
    sort(c(
      "protein_late_decline_magnitude",
      "inverse_braak_magnitude",
      "cognition_priority_score",
      "agora_nominated_target"
    )),
    setequal(
      sf3_summary$metric,
      c(
        "protein_late_decline_magnitude",
        "inverse_braak_magnitude",
        "cognition_priority_score",
        "agora_nominated_target"
      )
    )
  ),
  check_one(
    "Supp Fig 3 observed late decline matches Main Fig 5",
    abs(sf3_late_obs - main_late_obs),
    "<1e-12",
    is.finite(sf3_late_obs) &&
      is.finite(main_late_obs) &&
      abs(sf3_late_obs - main_late_obs) < 1e-12
  ),
  check_one(
    "Supp Fig 3 observed inverse Braak matches Main Fig 5",
    abs(sf3_braak_obs - main_braak_obs),
    "<1e-12",
    is.finite(sf3_braak_obs) &&
      is.finite(main_braak_obs) &&
      abs(sf3_braak_obs - main_braak_obs) < 1e-12
  ),
  check_one(
    "Supp Fig 3 observed cognition score matches Main Fig 5",
    abs(sf3_cog_obs - main_cog_obs),
    "<1e-12",
    is.finite(sf3_cog_obs) &&
      is.finite(main_cog_obs) &&
      abs(sf3_cog_obs - main_cog_obs) < 1e-12
  ),
  check_one(
    "Supp Fig 3 observed AGORA fraction matches Main Fig 5",
    abs(sf3_agora_obs - main_agora_obs),
    "<1e-12",
    is.finite(sf3_agora_obs) &&
      is.finite(main_agora_obs) &&
      abs(sf3_agora_obs - main_agora_obs) < 1e-12
  ),
  check_one(
    "Supp Fig 3 old collapse identifier absent",
    grepl("protein_collapse_magnitude", script12),
    FALSE,
    !grepl("protein_collapse_magnitude", script12)
  ),
  check_one(
    "Supp Fig 4 raw protein matrix prioritized",
    grepl(
      'c\\("prot_mat_raw", "prot_mat"',
      script13
    ),
    TRUE,
    grepl(
      'c\\("prot_mat_raw", "prot_mat"',
      script13
    )
  ),
  check_one(
    "Supp Fig 4 Hsp60/10 clients modeled",
    sum(supfig4_outputs$plot_tbl$group == "Hsp60/10 clients"),
    306,
    sum(supfig4_outputs$plot_tbl$group == "Hsp60/10 clients") == 306
  ),
  check_one(
    "Supp Fig 4 mitochondrial background modeled",
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
    "Supp Fig 4 TMT batch covariate detected",
    paste(batch_detected, collapse = ";"),
    "batch_factor",
    length(batch_detected) == 1 &&
      batch_detected[[1]] == "batch_factor"
  ),
  check_one(
    "Supp Fig 4 adjusted Braak reproduces production",
    max(sf4_compare$abs_difference, na.rm = TRUE),
    "<1e-10",
    nrow(sf4_compare) == 915 &&
      max(sf4_compare$abs_difference, na.rm = TRUE) < 1e-10
  ),
  check_one(
    "Supp Fig 4 old collapse/legacy diagnosis identifiers absent",
    grepl(
      "protein_collapse_magnitude|EmoryStrictDx\\.2019|diagnosis_model",
      script13
    ),
    FALSE,
    !grepl(
      "protein_collapse_magnitude|EmoryStrictDx\\.2019|diagnosis_model",
      script13
    )
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "supplemental_figures_3_4_migration_checks.csv")
)

readr::write_csv(
  sf3_summary,
  file.path(out_dir, "suppfig3_corrected_specificity_summary.csv")
)

readr::write_csv(
  supfig4_outputs$summary_tbl,
  file.path(out_dir, "suppfig4_corrected_pathology_summary.csv")
)

readr::write_csv(
  sf4_compare,
  file.path(out_dir, "suppfig4_vs_production_braak_agreement.csv")
)

message("\n============================================================")
message("SUPPLEMENTAL FIGURES 3-4 MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nSUPP FIG 3 CORRECTED SPECIFICITY SUMMARY")
print(
  sf3_summary |>
    dplyr::select(
      metric,
      observed_mean,
      null_mean,
      observed_to_null_ratio,
      empirical_p_greater
    ),
  n = Inf
)

message("\nSUPP FIG 4 PATHOLOGY ROBUSTNESS SUMMARY")
print(supfig4_outputs$summary_tbl, n = Inf)

message("\nSUPP FIG 4 VS PRODUCTION BRAAK AGREEMENT")
print(
  sf4_compare |>
    dplyr::summarise(
      n = dplyr::n(),
      max_abs_difference = max(.data$abs_difference, na.rm = TRUE),
      mean_abs_difference = mean(.data$abs_difference, na.rm = TRUE)
    )
)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))

  stop(
    "Supplemental Figures 3-4 migration FAILED. ",
    "Do not proceed to Supplemental Figures 5-6.",
    call. = FALSE
  )
}

message("\nALL SUPPLEMENTAL FIGURE 3-4 MIGRATION CHECKS PASSED.")
