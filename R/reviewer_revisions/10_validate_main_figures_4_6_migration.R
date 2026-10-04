############################################################
## 10_validate_main_figures_4_6_migration.R
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

source("R/13_make_main_figure_4_cognition.R")
source("R/14_make_main_figure_5_matched_null_specificity.R")
source("R/15_make_main_figure_6_candidate_classification.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "main_figures_4_6_migration"
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

broad_counts <- plot_df |>
  dplyr::filter(!is.na(.data$cog_group)) |>
  dplyr::count(.data$cog_group, name = "n")

get_group_n <- function(group_name) {
  hit <- broad_counts$n[
    as.character(broad_counts$cog_group) == group_name
  ]
  if (length(hit) == 0) 0L else as.integer(hit[[1]])
}

expected_broad <- clinical_stage_broad(cognition_model_df$cogdx)

broad_group_matches_upstream <- all(
  as.character(plot_df$cog_group) == as.character(expected_broad) |
    (
      is.na(plot_df$cog_group) &
      is.na(expected_broad)
    )
)

late_row <- null_summary |>
  dplyr::filter(.data$metric == "protein_late_decline_magnitude")

fig5_metrics <- as.character(figure_tbl$metric)

priority_late_fraction <- as.numeric(priority_tbl$late_decline_percentile)
if (any(abs(priority_late_fraction[is.finite(priority_late_fraction)]) > 1)) {
  priority_late_fraction <- priority_late_fraction / 100
}

fig6_source <- priority_late_fraction[
  match(fig6_tbl$gene, priority_tbl$gene)
]

script14 <- paste(
  readLines("R/14_make_main_figure_5_matched_null_specificity.R", warn = FALSE),
  collapse = "\n"
)
script15 <- paste(
  readLines("R/15_make_main_figure_6_candidate_classification.R", warn = FALSE),
  collapse = "\n"
)

checks <- dplyr::bind_rows(
  check_one(
    "Figure 4 broad cognition NCI n",
    get_group_n("NCI"),
    168,
    get_group_n("NCI") == 168
  ),
  check_one(
    "Figure 4 broad cognition MCI n",
    get_group_n("MCI"),
    101,
    get_group_n("MCI") == 101
  ),
  check_one(
    "Figure 4 broad cognition AD n",
    get_group_n("AD"),
    123,
    get_group_n("AD") == 123
  ),
  check_one(
    "Figure 4 broad groups match canonical helper",
    broad_group_matches_upstream,
    TRUE,
    broad_group_matches_upstream
  ),
  check_one(
    "Figure 4 network cognition models",
    nrow(network_cognition_tbl),
    3,
    nrow(network_cognition_tbl) == 3
  ),
  check_one(
    "Figure 4 client genes modeled",
    dplyr::n_distinct(client_cognition_all$gene),
    306,
    dplyr::n_distinct(client_cognition_all$gene) == 306
  ),
  check_one(
    "Figure 5 Hsp null clients",
    nrow(hsp_null_tbl),
    306,
    nrow(hsp_null_tbl) == 306
  ),
  check_one(
    "Figure 5 mitochondrial background",
    nrow(background_null_pool),
    609,
    nrow(background_null_pool) == 609
  ),
  check_one(
    "Figure 5 null iterations",
    nrow(null_results),
    10000,
    nrow(null_results) == 10000
  ),
  check_one(
    "Figure 5 includes corrected late-decline metric",
    "protein_late_decline_magnitude" %in% fig5_metrics,
    TRUE,
    "protein_late_decline_magnitude" %in% fig5_metrics
  ),
  check_one(
    "Figure 5 observed late decline preserved",
    round(late_row$observed, 4),
    "0.0230 +/- 0.005",
    nrow(late_row) == 1 &&
      abs(late_row$observed - 0.0230) <= 0.005
  ),
  check_one(
    "Figure 5 old collapse source identifier absent",
    grepl("protein_collapse_magnitude", script14),
    FALSE,
    !grepl("protein_collapse_magnitude", script14)
  ),
  check_one(
    "Figure 6 candidate rows",
    nrow(fig6_tbl),
    306,
    nrow(fig6_tbl) == 306
  ),
  check_one(
    "Figure 6 late-decline percentile reproduces priority table",
    max(
      abs(fig6_tbl$late_decline_pct - fig6_source),
      na.rm = TRUE
    ),
    "<1e-12",
    max(
      abs(fig6_tbl$late_decline_pct - fig6_source),
      na.rm = TRUE
    ) < 1e-12
  ),
  check_one(
    "Figure 6 old collapse analytical identifiers absent",
    grepl(
      "collapse_percentile|collapse_pct|collapse_col|median_collapse",
      script15
    ),
    FALSE,
    !grepl(
      "collapse_percentile|collapse_pct|collapse_col|median_collapse",
      script15
    )
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "main_figures_4_6_migration_checks.csv")
)

readr::write_csv(
  broad_counts,
  file.path(out_dir, "figure4_broad_cognition_group_counts.csv")
)

readr::write_csv(
  network_cognition_tbl,
  file.path(out_dir, "figure4_network_cognition_models.csv")
)

readr::write_csv(
  figure_tbl,
  file.path(out_dir, "figure5_corrected_summary.csv")
)

readr::write_csv(
  candidate_layer_counts,
  file.path(out_dir, "figure6_candidate_layer_counts.csv")
)

message("\n============================================================")
message("MAIN FIGURES 4-6 MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nFIGURE 4 NETWORK COGNITION MODELS")
print(network_cognition_tbl, n = Inf)

message("\nFIGURE 5 CORRECTED MATCHED-NULL SUMMARY")
print(
  figure_tbl |>
    dplyr::select(
      metric,
      observed,
      null_mean,
      fold_enrichment,
      empirical_p_greater
    ),
  n = Inf
)

message("\nFIGURE 6 CANDIDATE-LAYER COUNTS")
print(candidate_layer_counts, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))

  stop(
    "Main Figure 4-6 migration FAILED. Do not proceed to supplemental figures.",
    call. = FALSE
  )
}

message("\nALL MAIN FIGURE 4-6 MIGRATION CHECKS PASSED.")
