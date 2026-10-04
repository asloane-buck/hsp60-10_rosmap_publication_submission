############################################################
## 09_validate_main_figures_1_3_migration.R
############################################################

options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")
source("R/02_load_data.R")
source("R/03_build_adjusted_core_objects.R")
source("R/04_build_pathway_sets_all_clients.R")
source("R/05_build_adjusted_all_client_tables.R")

source("R/10_make_main_figure_1_pathway_remodeling.R")
source("R/11_make_main_figure_2_collapse_heterogeneity.R")
source("R/12_make_main_figure_3_pathology_coupling.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "main_figures_1_3_migration"
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

fig1_counts <- overlay_long_all_clients |>
  dplyr::filter(.data$Pathway == "Hsp60_10_all_clients") |>
  dplyr::group_by(.data$Modality, .data$Stage) |>
  dplyr::summarise(n = sum(is.finite(.data$Score)), .groups = "drop")

get_fig1_n <- function(modality, stage) {
  hit <- fig1_counts$n[
    as.character(fig1_counts$Modality) == modality &
    as.character(fig1_counts$Stage) == stage
  ]
  if (length(hit) == 0) 0L else as.integer(hit[[1]])
}

script10 <- paste(
  readLines("R/10_make_main_figure_1_pathway_remodeling.R", warn = FALSE),
  collapse = "\n"
)
script11 <- paste(
  readLines("R/11_make_main_figure_2_collapse_heterogeneity.R", warn = FALSE),
  collapse = "\n"
)
script12 <- paste(
  readLines("R/12_make_main_figure_3_pathology_coupling.R", warn = FALSE),
  collapse = "\n"
)

fig2_source <- all_hsp60_10_client_tbl$protein_late_decline_magnitude[
  match(fig2_tbl$gene, all_hsp60_10_client_tbl$gene)
]
fig3_source <- all_hsp60_10_client_tbl$protein_late_decline_magnitude[
  match(fig3_tbl$gene, all_hsp60_10_client_tbl$gene)
]

checks <- dplyr::bind_rows(
  check_one(
    "Figure 1 stage levels",
    levels(overlay_long_all_clients$Stage),
    c("NCI", "MCI", "AD"),
    identical(levels(overlay_long_all_clients$Stage), c("NCI", "MCI", "AD"))
  ),
  check_one("Figure 1 protein Hsp NCI n", get_fig1_n("Protein", "NCI"), 167, get_fig1_n("Protein", "NCI") == 167),
  check_one("Figure 1 protein Hsp MCI n", get_fig1_n("Protein", "MCI"), 96, get_fig1_n("Protein", "MCI") == 96),
  check_one("Figure 1 protein Hsp AD n", get_fig1_n("Protein", "AD"), 109, get_fig1_n("Protein", "AD") == 109),
  check_one("Figure 1 RNA Hsp NCI n", get_fig1_n("RNA", "NCI"), 200, get_fig1_n("RNA", "NCI") == 200),
  check_one("Figure 1 RNA Hsp MCI n", get_fig1_n("RNA", "MCI"), 158, get_fig1_n("RNA", "MCI") == 158),
  check_one("Figure 1 RNA Hsp AD n", get_fig1_n("RNA", "AD"), 219, get_fig1_n("RNA", "AD") == 219),
  check_one(
    "Figure 1 invalid AsymAD-to-MCI mapping absent",
    grepl('AsymAD\\s*=\\s*"MCI"', script10),
    FALSE,
    !grepl('AsymAD\\s*=\\s*"MCI"', script10)
  ),
  check_one("Figure 2 source clients", nrow(fig2_tbl), 306, nrow(fig2_tbl) == 306),
  check_one(
    "Figure 2 top-quartile late-decline n",
    sum(fig2_tbl$late_decline_group == "Top-quartile late decline"),
    ceiling(0.25 * 306),
    sum(fig2_tbl$late_decline_group == "Top-quartile late decline") == ceiling(0.25 * 306)
  ),
  check_one(
    "Figure 2 late-decline metric reproduces source",
    max(abs(fig2_tbl$late_decline_value - fig2_source), na.rm = TRUE),
    "<1e-12",
    max(abs(fig2_tbl$late_decline_value - fig2_source), na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Figure 2 legacy source identifiers absent",
    grepl(
      "protein_collapse_magnitude|abs_delta_late|protein_minus_rna_collapse",
      script11
    ),
    FALSE,
    !grepl(
      "protein_collapse_magnitude|abs_delta_late|protein_minus_rna_collapse",
      script11
    )
  ),
  check_one("Figure 3 source clients", nrow(fig3_tbl), 306, nrow(fig3_tbl) == 306),
  check_one(
    "Figure 3 late-decline metric reproduces source",
    max(abs(fig3_tbl$late_decline_value - fig3_source), na.rm = TRUE),
    "<1e-12",
    max(abs(fig3_tbl$late_decline_value - fig3_source), na.rm = TRUE) < 1e-12
  ),
  check_one(
    "Figure 3 legacy source identifier absent",
    grepl("protein_collapse_magnitude", script12),
    FALSE,
    !grepl("protein_collapse_magnitude", script12)
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "main_figures_1_3_migration_checks.csv")
)
readr::write_csv(
  fig1_counts,
  file.path(out_dir, "figure1_hsp_stage_counts.csv")
)
readr::write_csv(
  fig2_stats_tbl,
  file.path(out_dir, "figure2_corrected_stats.csv")
)
readr::write_csv(
  fig3_stats_tbl,
  file.path(out_dir, "figure3_corrected_stats.csv")
)

message("\n============================================================")
message("MAIN FIGURES 1-3 MIGRATION CHECKPOINT")
message("============================================================\n")
print(checks, n = Inf)

message("\nFIGURE 2 CORRECTED STATS")
print(fig2_stats_tbl, n = Inf)

message("\nFIGURE 3 CORRECTED STATS")
print(fig3_stats_tbl, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))
  stop(
    "Main Figure 1-3 migration FAILED. Do not proceed to Figures 4-6.",
    call. = FALSE
  )
}

message("\nALL MAIN FIGURE 1-3 MIGRATION CHECKS PASSED.")
