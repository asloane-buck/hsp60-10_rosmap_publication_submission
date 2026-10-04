############################################################
## 11_validate_supplemental_figures_1_2_migration.R
############################################################

options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")
source("R/02_load_data.R")
source("R/03_build_adjusted_core_objects.R")
source("R/04_build_pathway_sets_all_clients.R")
source("R/05_build_adjusted_all_client_tables.R")

source("R/supplemental/00_supplemental_config.R")
source("R/supplemental/01_supplemental_load_inputs.R")
source("R/supplemental/02_supplemental_helper_functions.R")
source("R/supplemental/10_make_supplementary_figure_1_cohort_detection.R")
source("R/supplemental/11_make_supplementary_figure_2_matched_individual_sensitivity.R")

out_dir <- file.path(
  cfg$output_root,
  "reviewer_revisions",
  "supplemental_figures_1_2_migration"
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

canonical_rna <- rna_meta_adj |>
  dplyr::filter(
    .data$sample_id %in% colnames(rna_mat),
    !is.na(.data$clinical_stage)
  ) |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    stage = as.character(.data$clinical_stage)
  ) |>
  dplyr::distinct()

canonical_protein <- prot_meta_adj |>
  dplyr::filter(
    .data$SampleID %in% colnames(prot_mat),
    !is.na(.data$clinical_stage)
  ) |>
  dplyr::transmute(
    individual_id = as.character(.data$IndividualID),
    stage = as.character(.data$clinical_stage)
  ) |>
  dplyr::distinct()

expected_crosswalk <- canonical_protein |>
  dplyr::rename(protein_stage = stage) |>
  dplyr::inner_join(
    canonical_rna |>
      dplyr::rename(rna_stage = stage),
    by = "individual_id"
  )

expected_stage_agreement <- all(
  expected_crosswalk$protein_stage == expected_crosswalk$rna_stage
)

sf1_stage_audit <- readr::read_csv(
  file.path(audits_dir, "SuppFig1_stage_sample_size_audit.csv"),
  show_col_types = FALSE
)

client_genes <- clean_gene(supfig1_outputs$interactor_tbl$gene)
rna_genes <- rownames(supfig1_outputs$rna_info$matrix)
protein_genes <- rownames(supfig1_outputs$protein_info$matrix)

n_both <- sum(client_genes %in% rna_genes & client_genes %in% protein_genes)
n_rna_only <- sum(client_genes %in% rna_genes & !client_genes %in% protein_genes)
n_protein_only <- sum(!client_genes %in% rna_genes & client_genes %in% protein_genes)
n_neither <- sum(!client_genes %in% rna_genes & !client_genes %in% protein_genes)

get_sf1_stage_n <- function(modality, stage) {
  hit <- sf1_stage_audit$n_samples[
    sf1_stage_audit$modality == modality &
      as.character(sf1_stage_audit$stage) == stage
  ]
  if (length(hit) == 0) 0L else as.integer(hit[[1]])
}

sf2_crosswalk <- supfig2_outputs$matched_stage_crosswalk
sf2_stages <- sort(unique(as.character(supfig2_outputs$matched_long$Stage)))

script01 <- paste(
  readLines("R/supplemental/01_supplemental_load_inputs.R", warn = FALSE),
  collapse = "\n"
)
script10 <- paste(
  readLines("R/supplemental/10_make_supplementary_figure_1_cohort_detection.R", warn = FALSE),
  collapse = "\n"
)
script11 <- paste(
  readLines("R/supplemental/11_make_supplementary_figure_2_matched_individual_sensitivity.R", warn = FALSE),
  collapse = "\n"
)

checks <- dplyr::bind_rows(
  check_one(
    "Supplemental protein matrix samples",
    ncol(supfig1_outputs$protein_info$matrix),
    400,
    ncol(supfig1_outputs$protein_info$matrix) == 400
  ),
  check_one(
    "Supplemental RNA matrix samples",
    ncol(supfig1_outputs$rna_info$matrix),
    577,
    ncol(supfig1_outputs$rna_info$matrix) == 577
  ),
  check_one(
    "Supp Fig 1 client inventory",
    length(client_genes),
    321,
    length(client_genes) == 321
  ),
  check_one(
    "Supp Fig 1 detected both",
    n_both,
    285,
    n_both == 285
  ),
  check_one(
    "Supp Fig 1 RNA only",
    n_rna_only,
    12,
    n_rna_only == 12
  ),
  check_one(
    "Supp Fig 1 protein only",
    n_protein_only,
    21,
    n_protein_only == 21
  ),
  check_one(
    "Supp Fig 1 neither",
    n_neither,
    3,
    n_neither == 3
  ),
  check_one(
    "Supp Fig 1 RNA NCI",
    get_sf1_stage_n("RNA", "NCI"),
    200,
    get_sf1_stage_n("RNA", "NCI") == 200
  ),
  check_one(
    "Supp Fig 1 RNA MCI",
    get_sf1_stage_n("RNA", "MCI"),
    158,
    get_sf1_stage_n("RNA", "MCI") == 158
  ),
  check_one(
    "Supp Fig 1 RNA AD",
    get_sf1_stage_n("RNA", "AD"),
    219,
    get_sf1_stage_n("RNA", "AD") == 219
  ),
  check_one(
    "Supp Fig 1 protein NCI phenotype-eligible",
    get_sf1_stage_n("Protein", "NCI"),
    168,
    get_sf1_stage_n("Protein", "NCI") == 168
  ),
  check_one(
    "Supp Fig 1 protein MCI phenotype-eligible",
    get_sf1_stage_n("Protein", "MCI"),
    97,
    get_sf1_stage_n("Protein", "MCI") == 97
  ),
  check_one(
    "Supp Fig 1 protein AD phenotype-eligible",
    get_sf1_stage_n("Protein", "AD"),
    109,
    get_sf1_stage_n("Protein", "AD") == 109
  ),
  check_one(
    "Production cross-modal stage agreement",
    expected_stage_agreement,
    TRUE,
    expected_stage_agreement
  ),
  check_one(
    "Supp Fig 2 matched participants equal production intersection",
    nrow(sf2_crosswalk),
    nrow(expected_crosswalk),
    nrow(sf2_crosswalk) == nrow(expected_crosswalk) &&
      setequal(sf2_crosswalk$individual_id, expected_crosswalk$individual_id)
  ),
  check_one(
    "Supp Fig 2 cross-modal stage agreement",
    all(sf2_crosswalk$stage_agrees),
    TRUE,
    all(sf2_crosswalk$stage_agrees)
  ),
  check_one(
    "Supp Fig 2 stage levels",
    sf2_stages,
    c("AD", "MCI", "NCI"),
    identical(sf2_stages, c("AD", "MCI", "NCI"))
  ),
  check_one(
    "Supplemental loader prioritizes clinical_stage",
    grepl(
      'protein_stage_candidates <- c\\(\\n  "clinical_stage"',
      script01
    ),
    TRUE,
    grepl(
      'protein_stage_candidates <- c\\(\\n  "clinical_stage"',
      script01
    )
  ),
  check_one(
    "Supp Fig 1 no legacy Emory stage source",
    grepl('"EmoryStrictDx.2019"', script10),
    FALSE,
    !grepl('"EmoryStrictDx.2019"', script10)
  ),
  check_one(
    "Supp Fig 2 no legacy Emory stage source",
    grepl('"EmoryStrictDx.2019"', script11),
    FALSE,
    !grepl('"EmoryStrictDx.2019"', script11)
  )
)

readr::write_csv(
  checks,
  file.path(out_dir, "supplemental_figures_1_2_migration_checks.csv")
)

readr::write_csv(
  sf1_stage_audit,
  file.path(out_dir, "suppfig1_stage_sample_size_audit.csv")
)

readr::write_csv(
  sf2_crosswalk,
  file.path(out_dir, "suppfig2_matched_stage_crosswalk.csv")
)

readr::write_csv(
  supfig2_outputs$matched_summary,
  file.path(out_dir, "suppfig2_matched_summary.csv")
)

message("\n============================================================")
message("SUPPLEMENTAL FIGURES 1-2 MIGRATION CHECKPOINT")
message("============================================================\n")

print(checks, n = Inf)

message("\nSUPP FIG 1 STAGE COUNTS")
print(sf1_stage_audit, n = Inf)

message("\nSUPP FIG 2 MATCHED STAGE COUNTS")
print(
  sf2_crosswalk |>
    dplyr::count(.data$protein_stage, name = "n"),
  n = Inf
)

message("\nSUPP FIG 2 MATCHED PATHWAY SUMMARY")
print(supfig2_outputs$matched_summary, n = Inf)

if (length(failures) > 0) {
  message("\nFAILED CHECKS:")
  message(paste0(" - ", failures, collapse = "\n"))

  stop(
    "Supplemental Figures 1-2 migration FAILED. ",
    "Do not proceed to Supplemental Figures 3-6.",
    call. = FALSE
  )
}

message("\nALL SUPPLEMENTAL FIGURE 1-2 MIGRATION CHECKS PASSED.")
