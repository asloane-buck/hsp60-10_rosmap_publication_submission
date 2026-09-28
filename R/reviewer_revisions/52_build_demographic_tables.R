############################################################
## 52_build_demographic_tables.R
##
## Canonical demographic tables for reviewer revision.
##
## Cohorts:
##   1. Corrected RNA cohort: n = 577
##   2. Full protein cohort: n = 400
##   3. Pure-stage protein cohort: n = 374
##   4. Nuisance-complete pure-stage protein cohort: n = 372
##   5. Validated paired DLPFC/STG regional cohort: n = 215
##
## Outputs are analysis/reviewer-ready. Manuscript formatting
## can be performed later without recomputing cohort summaries.
############################################################

options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  if (!file.exists(f)) {
    stop("Missing required production script: ", f)
  }
  source(f)
}

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "demographic_tables"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Required files / objects
############################################################

required_objects <- c(
  "rna_meta_adj",
  "prot_meta_adj",
  "prot_mat_adj",
  "protein_covars"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0L) {
  stop(
    "Missing required production objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

apoe_file <- file.path(
  "outputs",
  "reviewer_revisions",
  "APOE_sensitivity",
  "APOE_participant_lookup.csv"
)

regional_file <- file.path(
  "outputs",
  "reviewer_revisions",
  "regional_formal_interaction",
  "paired_participant_analysis_manifest.csv"
)

for (f in c(apoe_file, regional_file)) {
  if (!file.exists(f)) {
    stop("Missing required input file: ", f)
  }
}

############################################################
## 2. Helper functions
############################################################

normalize_sex <- function(x) {
  z <- trimws(as.character(x))
  dplyr::case_when(
    z %in% c("1", "Male", "male", "M", "m") ~ "Male",
    z %in% c("0", "Female", "female", "F", "f") ~ "Female",
    is.na(x) | z == "" | z %in% c("NA", "NaN") ~ NA_character_,
    TRUE ~ z
  )
}

num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

summarize_continuous <- function(
  df,
  cohort,
  group,
  variable,
  label
) {
  x <- num(df[[variable]])
  n_total <- length(x)
  n_nonmissing <- sum(is.finite(x))
  n_missing <- n_total - n_nonmissing

  if (n_nonmissing == 0L) {
    return(
      tibble::tibble(
        cohort = cohort,
        group = group,
        variable = variable,
        label = label,
        n_total = n_total,
        n_nonmissing = 0L,
        n_missing = n_missing,
        mean = NA_real_,
        sd = NA_real_,
        median = NA_real_,
        q1 = NA_real_,
        q3 = NA_real_,
        min = NA_real_,
        max = NA_real_
      )
    )
  }

  xf <- x[is.finite(x)]

  tibble::tibble(
    cohort = cohort,
    group = group,
    variable = variable,
    label = label,
    n_total = n_total,
    n_nonmissing = n_nonmissing,
    n_missing = n_missing,
    mean = mean(xf),
    sd = if (length(xf) > 1L) stats::sd(xf) else NA_real_,
    median = stats::median(xf),
    q1 = as.numeric(stats::quantile(xf, 0.25, names = FALSE)),
    q3 = as.numeric(stats::quantile(xf, 0.75, names = FALSE)),
    min = min(xf),
    max = max(xf)
  )
}

summarize_sex <- function(df, cohort, group) {
  x <- normalize_sex(df$sex_display)
  n_total <- length(x)
  n_available <- sum(!is.na(x))

  tibble::tibble(
    cohort = cohort,
    group = group,
    n_total = n_total,
    n_sex_available = n_available,
    n_sex_missing = n_total - n_available,
    n_female = sum(x == "Female", na.rm = TRUE),
    n_male = sum(x == "Male", na.rm = TRUE),
    pct_female_among_available = if (n_available > 0L) {
      100 * sum(x == "Female", na.rm = TRUE) / n_available
    } else {
      NA_real_
    },
    pct_male_among_available = if (n_available > 0L) {
      100 * sum(x == "Male", na.rm = TRUE) / n_available
    } else {
      NA_real_
    }
  )
}

summarize_apoe <- function(df, cohort, group) {
  carrier <- df$apoe4_carrier
  n_total <- nrow(df)
  n_available <- sum(!is.na(carrier))

  tibble::tibble(
    cohort = cohort,
    group = group,
    n_total = n_total,
    n_apoe_available = n_available,
    n_apoe_missing = n_total - n_available,
    n_e4_noncarrier = sum(carrier %in% FALSE, na.rm = TRUE),
    n_e4_carrier = sum(carrier %in% TRUE, na.rm = TRUE),
    pct_e4_carrier_among_available = if (n_available > 0L) {
      100 * sum(carrier %in% TRUE, na.rm = TRUE) / n_available
    } else {
      NA_real_
    }
  )
}

mean_sd_string <- function(x) {
  if (nrow(x) != 1L || !is.finite(x$mean)) return("NA")
  sprintf("%.1f (%.1f)", x$mean, x$sd)
}

median_iqr_string <- function(x) {
  if (nrow(x) != 1L || !is.finite(x$median)) return("NA")
  sprintf("%.1f [%.1f–%.1f]", x$median, x$q1, x$q3)
}

############################################################
## 3. APOE participant lookup
############################################################

apoe_lookup <- utils::read.csv(
  apoe_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_apoe_cols <- c(
  "individual_id",
  "apoe_genotype",
  "apoe4_carrier",
  "apoe4_allele_count"
)

missing_apoe_cols <- setdiff(
  required_apoe_cols,
  names(apoe_lookup)
)

if (length(missing_apoe_cols) > 0L) {
  stop(
    "APOE lookup missing required columns: ",
    paste(missing_apoe_cols, collapse = ", ")
  )
}

if (anyDuplicated(apoe_lookup$individual_id)) {
  stop("APOE lookup is not unique by individual_id.")
}

############################################################
## 4. Corrected RNA cohort
############################################################

required_rna_cols <- c(
  "sample_id",
  "individual_id",
  "clinical_stage",
  "age_num",
  "sex_label",
  "pmi_num",
  "rin_num",
  "braak_num",
  "cerad_num"
)

missing_rna_cols <- setdiff(
  required_rna_cols,
  names(rna_meta_adj)
)

if (length(missing_rna_cols) > 0L) {
  stop(
    "RNA metadata missing required columns: ",
    paste(missing_rna_cols, collapse = ", ")
  )
}

rna_demo <- rna_meta_adj |>
  dplyr::transmute(
    sample_id = as.character(.data$sample_id),
    individual_id = as.character(.data$individual_id),
    stage = as.character(.data$clinical_stage),
    age_death = num(.data$age_num),
    sex_display = normalize_sex(.data$sex_label),
    pmi = num(.data$pmi_num),
    rin = num(.data$rin_num),
    braak = num(.data$braak_num),
    cerad = num(.data$cerad_num)
  ) |>
  dplyr::left_join(
    apoe_lookup |>
      dplyr::select(
        individual_id,
        apoe_genotype,
        apoe4_carrier,
        apoe4_allele_count
      ),
    by = "individual_id"
  )

if (nrow(rna_demo) != 577L) {
  stop("Expected corrected RNA cohort n=577; found ", nrow(rna_demo))
}

if (anyDuplicated(rna_demo$individual_id)) {
  stop("Corrected RNA cohort is not unique by individual_id.")
}

rna_stage_counts <- table(
  factor(rna_demo$stage, levels = c("NCI", "MCI", "AD"))
)

if (!identical(as.integer(rna_stage_counts), c(200L, 158L, 219L))) {
  stop(
    "Corrected RNA stage counts changed: ",
    paste(as.integer(rna_stage_counts), collapse = "/")
  )
}

############################################################
## 5. Protein cohorts
############################################################

required_protein_cols <- c(
  "SampleID",
  "individual_id",
  "clinical_stage",
  "age_num",
  "sex_label",
  "pmi_num",
  "braak_num",
  "cerad_num"
)

missing_protein_cols <- setdiff(
  required_protein_cols,
  names(prot_meta_adj)
)

if (length(missing_protein_cols) > 0L) {
  stop(
    "Protein metadata missing required columns: ",
    paste(missing_protein_cols, collapse = ", ")
  )
}

protein_demo_all <- prot_meta_adj |>
  dplyr::transmute(
    sample_id = as.character(.data$SampleID),
    individual_id = as.character(.data$individual_id),
    stage = as.character(.data$clinical_stage),
    age_death = num(.data$age_num),
    sex_display = normalize_sex(.data$sex_label),
    pmi = num(.data$pmi_num),
    braak = num(.data$braak_num),
    cerad = num(.data$cerad_num)
  ) |>
  dplyr::left_join(
    apoe_lookup |>
      dplyr::select(
        individual_id,
        apoe_genotype,
        apoe4_carrier,
        apoe4_allele_count
      ),
    by = "individual_id"
  )

if (nrow(protein_demo_all) != 400L) {
  stop("Expected full protein cohort n=400; found ", nrow(protein_demo_all))
}

if (anyDuplicated(protein_demo_all$individual_id)) {
  stop("Full protein cohort is not unique by individual_id.")
}

protein_primary <- protein_demo_all |>
  dplyr::filter(.data$stage %in% c("NCI", "MCI", "AD"))

if (nrow(protein_primary) != 374L) {
  stop(
    "Expected pure-stage protein cohort n=374; found ",
    nrow(protein_primary)
  )
}

protein_primary_counts <- table(
  factor(protein_primary$stage, levels = c("NCI", "MCI", "AD"))
)

if (!identical(
  as.integer(protein_primary_counts),
  c(168L, 97L, 109L)
)) {
  stop(
    "Pure-stage protein counts changed: ",
    paste(as.integer(protein_primary_counts), collapse = "/")
  )
}

## Reproduce the canonical analyzability definition from R/03.
meta_index <- match(
  protein_demo_all$sample_id,
  prot_meta_adj$SampleID
)

if (any(is.na(meta_index))) {
  stop("Could not align protein demographic table to prot_meta_adj.")
}

nuisance_complete <- stats::complete.cases(
  prot_meta_adj[
    meta_index,
    protein_covars,
    drop = FALSE
  ]
)

adjusted_data_available_by_sample <- setNames(
  colSums(is.finite(prot_mat_adj)) > 0,
  colnames(prot_mat_adj)
)

has_adjusted_data <- unname(
  adjusted_data_available_by_sample[
    protein_demo_all$sample_id
  ]
)

if (any(is.na(has_adjusted_data))) {
  stop("Could not align adjusted protein matrix availability.")
}

protein_demo_all$nuisance_complete <- nuisance_complete
protein_demo_all$adjusted_data_available <- has_adjusted_data

protein_stage_complete <- protein_demo_all |>
  dplyr::filter(
    .data$stage %in% c("NCI", "MCI", "AD"),
    .data$nuisance_complete,
    .data$adjusted_data_available
  )

if (nrow(protein_stage_complete) != 372L) {
  stop(
    "Expected nuisance-complete pure-stage protein cohort n=372; found ",
    nrow(protein_stage_complete)
  )
}

protein_complete_counts <- table(
  factor(
    protein_stage_complete$stage,
    levels = c("NCI", "MCI", "AD")
  )
)

if (!identical(
  as.integer(protein_complete_counts),
  c(167L, 96L, 109L)
)) {
  stop(
    "Nuisance-complete protein stage counts changed: ",
    paste(as.integer(protein_complete_counts), collapse = "/")
  )
}

############################################################
## 6. Validated paired regional cohort
############################################################

regional_raw <- utils::read.csv(
  regional_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_regional_cols <- c(
  "individual_id",
  "ad_status",
  "age_death_num",
  "sex",
  "pmi_num",
  "Braak"
)

missing_regional_cols <- setdiff(
  required_regional_cols,
  names(regional_raw)
)

if (length(missing_regional_cols) > 0L) {
  stop(
    "Regional manifest missing required columns: ",
    paste(missing_regional_cols, collapse = ", ")
  )
}

regional_demo <- regional_raw |>
  dplyr::transmute(
    individual_id = as.character(.data$individual_id),
    ad_status = as.character(.data$ad_status),
    age_death = num(.data$age_death_num),
    sex_display = normalize_sex(.data$sex),
    pmi = num(.data$pmi_num),
    braak = num(.data$Braak)
  ) |>
  dplyr::left_join(
    apoe_lookup |>
      dplyr::select(
        individual_id,
        apoe_genotype,
        apoe4_carrier,
        apoe4_allele_count
      ),
    by = "individual_id"
  )

if (nrow(regional_demo) != 215L) {
  stop(
    "Expected paired regional cohort n=215; found ",
    nrow(regional_demo)
  )
}

if (anyDuplicated(regional_demo$individual_id)) {
  stop("Paired regional cohort is not unique by individual_id.")
}

############################################################
## 7. Cohort/group registry
############################################################

cohort_groups <- list(
  list(
    cohort = "RNA_corrected_577",
    group = "Overall",
    data = rna_demo
  ),
  list(
    cohort = "RNA_corrected_577",
    group = "NCI",
    data = dplyr::filter(rna_demo, .data$stage == "NCI")
  ),
  list(
    cohort = "RNA_corrected_577",
    group = "MCI",
    data = dplyr::filter(rna_demo, .data$stage == "MCI")
  ),
  list(
    cohort = "RNA_corrected_577",
    group = "AD",
    data = dplyr::filter(rna_demo, .data$stage == "AD")
  ),
  list(
    cohort = "Protein_full_400",
    group = "Overall",
    data = protein_demo_all
  ),
  list(
    cohort = "Protein_primary_374",
    group = "Overall",
    data = protein_primary
  ),
  list(
    cohort = "Protein_primary_374",
    group = "NCI",
    data = dplyr::filter(protein_primary, .data$stage == "NCI")
  ),
  list(
    cohort = "Protein_primary_374",
    group = "MCI",
    data = dplyr::filter(protein_primary, .data$stage == "MCI")
  ),
  list(
    cohort = "Protein_primary_374",
    group = "AD",
    data = dplyr::filter(protein_primary, .data$stage == "AD")
  ),
  list(
    cohort = "Protein_stage_nuisance_complete_372",
    group = "Overall",
    data = protein_stage_complete
  ),
  list(
    cohort = "Protein_stage_nuisance_complete_372",
    group = "NCI",
    data = dplyr::filter(protein_stage_complete, .data$stage == "NCI")
  ),
  list(
    cohort = "Protein_stage_nuisance_complete_372",
    group = "MCI",
    data = dplyr::filter(protein_stage_complete, .data$stage == "MCI")
  ),
  list(
    cohort = "Protein_stage_nuisance_complete_372",
    group = "AD",
    data = dplyr::filter(protein_stage_complete, .data$stage == "AD")
  ),
  list(
    cohort = "Regional_paired_DLPFC_STG_215",
    group = "Overall",
    data = regional_demo
  )
)

regional_status_levels <- sort(unique(
  regional_demo$ad_status[
    !is.na(regional_demo$ad_status) &
      nzchar(regional_demo$ad_status)
  ]
))

for (lev in regional_status_levels) {
  cohort_groups[[length(cohort_groups) + 1L]] <- list(
    cohort = "Regional_paired_DLPFC_STG_215",
    group = paste0("AD_status=", lev),
    data = dplyr::filter(regional_demo, .data$ad_status == lev)
  )
}

cohort_counts <- purrr::map_dfr(
  cohort_groups,
  function(x) {
    tibble::tibble(
      cohort = x$cohort,
      group = x$group,
      n = nrow(x$data)
    )
  }
)

############################################################
## 8. Continuous summaries
############################################################

continuous_summaries <- list()

for (cg in cohort_groups) {
  df <- cg$data

  variable_map <- c(
    age_death = "Age at death, years",
    pmi = "Postmortem interval",
    braak = "Braak stage"
  )

  ## CERAD is available in the canonical RNA/protein cohorts.
  if ("cerad" %in% names(df)) {
    variable_map <- c(
      variable_map,
      cerad = "CERAD score"
    )
  }

  ## RIN is RNA-specific.
  if ("rin" %in% names(df)) {
    variable_map <- c(
      variable_map,
      rin = "RNA integrity number"
    )
  }

  for (v in names(variable_map)) {
    continuous_summaries[[length(continuous_summaries) + 1L]] <-
      summarize_continuous(
        df = df,
        cohort = cg$cohort,
        group = cg$group,
        variable = v,
        label = unname(variable_map[[v]])
      )
  }
}

continuous_summary <- dplyr::bind_rows(continuous_summaries)

############################################################
## 9. Sex and APOE summaries
############################################################

sex_summary <- purrr::map_dfr(
  cohort_groups,
  function(x) {
    summarize_sex(
      x$data,
      x$cohort,
      x$group
    )
  }
)

apoe_summary <- purrr::map_dfr(
  cohort_groups,
  function(x) {
    summarize_apoe(
      x$data,
      x$cohort,
      x$group
    )
  }
)

## The validated regional manifest uses a different participant-ID
## namespace from the canonical APOE lookup (0/215 IDs overlap).
## Therefore APOE is not linkable for this cohort; this is not
## biological APOE missingness.
regional_apoe_rows <- (
  apoe_summary$cohort == "Regional_paired_DLPFC_STG_215"
)

apoe_summary$n_apoe_available[regional_apoe_rows] <- NA_integer_
apoe_summary$n_apoe_missing[regional_apoe_rows] <- NA_integer_
apoe_summary$n_e4_noncarrier[regional_apoe_rows] <- NA_integer_
apoe_summary$n_e4_carrier[regional_apoe_rows] <- NA_integer_
apoe_summary$pct_e4_carrier_among_available[regional_apoe_rows] <- NA_real_

############################################################
## 10. Stage/group counts
############################################################

stage_counts <- dplyr::bind_rows(
  rna_demo |>
    dplyr::count(.data$stage, name = "n") |>
    dplyr::mutate(
      cohort = "RNA_corrected_577",
      grouping_variable = "clinical_stage",
      group = .data$stage,
      .before = 1
    ) |>
    dplyr::select(-stage),

  protein_primary |>
    dplyr::count(.data$stage, name = "n") |>
    dplyr::mutate(
      cohort = "Protein_primary_374",
      grouping_variable = "clinical_stage",
      group = .data$stage,
      .before = 1
    ) |>
    dplyr::select(-stage),

  protein_stage_complete |>
    dplyr::count(.data$stage, name = "n") |>
    dplyr::mutate(
      cohort = "Protein_stage_nuisance_complete_372",
      grouping_variable = "clinical_stage",
      group = .data$stage,
      .before = 1
    ) |>
    dplyr::select(-stage),

  regional_demo |>
    dplyr::count(.data$ad_status, name = "n") |>
    dplyr::mutate(
      cohort = "Regional_paired_DLPFC_STG_215",
      grouping_variable = "regional_AD_status",
      group = as.character(.data$ad_status),
      .before = 1
    ) |>
    dplyr::select(-ad_status)
)

############################################################
## 11. Regional brain-region coverage table
############################################################

regional_region_counts <- tibble::tibble(
  cohort = "Regional_paired_DLPFC_STG_215",
  region = c("DLPFC", "STG"),
  n_participants = c(215L, 215L),
  paired_design = TRUE,
  note = c(
    "Each participant contributes DLPFC in the validated paired regional analysis.",
    "Each participant contributes STG in the validated paired regional analysis."
  )
)

############################################################
## 12. Variable availability / interpretation table
############################################################

variable_availability <- tibble::tribble(
  ~cohort, ~variable, ~available, ~note,
  "RNA_corrected_577", "Age at death", TRUE,
  "Canonical corrected RNA metadata: age_num.",
  "RNA_corrected_577", "Sex", TRUE,
  "Canonical corrected RNA metadata: sex_label.",
  "RNA_corrected_577", "PMI", TRUE,
  "Canonical corrected RNA metadata: pmi_num.",
  "RNA_corrected_577", "RIN", TRUE,
  "RNA-specific technical/quality variable: rin_num.",
  "RNA_corrected_577", "Braak", TRUE,
  "Canonical corrected RNA metadata: braak_num.",
  "RNA_corrected_577", "CERAD", TRUE,
  "Canonical corrected RNA metadata: cerad_num.",
  "RNA_corrected_577", "APOE e4", TRUE,
  "Joined by individual_id from frozen APOE participant lookup.",
  "RNA_corrected_577", "Education", FALSE,
  "Not present in the canonical main RNA metadata used for this table.",

  "Protein_full_400", "Age at death", TRUE,
  "Canonical protein metadata: age_num.",
  "Protein_full_400", "Sex", TRUE,
  "Canonical protein metadata: sex_label.",
  "Protein_full_400", "PMI", TRUE,
  "Canonical protein metadata: pmi_num.",
  "Protein_full_400", "RIN", FALSE,
  "RIN is not a protein nuisance covariate and is not reported for protein demographics.",
  "Protein_full_400", "Braak", TRUE,
  "Canonical protein metadata: braak_num.",
  "Protein_full_400", "CERAD", TRUE,
  "Canonical protein metadata: cerad_num.",
  "Protein_full_400", "APOE e4", TRUE,
  "Joined by individual_id from frozen APOE participant lookup.",
  "Protein_full_400", "Education", FALSE,
  "Not present in the canonical main protein metadata used for this table.",

  "Regional_paired_DLPFC_STG_215", "Age at death", TRUE,
  "Validated paired regional manifest: age_death_num.",
  "Regional_paired_DLPFC_STG_215", "Sex", TRUE,
  "Validated paired regional manifest: sex.",
  "Regional_paired_DLPFC_STG_215", "PMI", TRUE,
  "Validated paired regional manifest: pmi_num.",
  "Regional_paired_DLPFC_STG_215", "Braak", TRUE,
  "Validated paired regional manifest: Braak.",
  "Regional_paired_DLPFC_STG_215", "CERAD", FALSE,
  "Not labeled as CERAD in the validated paired regional manifest; not inferred from bScore.",
  "Regional_paired_DLPFC_STG_215", "APOE e4", FALSE,
  "Unavailable from the canonical APOE lookup because regional and canonical participant-ID namespaces have 0/215 overlap.",
  "Regional_paired_DLPFC_STG_215", "Education", FALSE,
  "Not present in the validated paired regional manifest."
)

############################################################
## 13. Analysis-ready manuscript display table
##
## This is a formatting convenience only. Long-format source
## tables above remain the analytical source of truth.
############################################################

display_groups <- cohort_counts |>
  dplyr::filter(
    .data$group %in% c("Overall", "NCI", "MCI", "AD")
  ) |>
  dplyr::mutate(
    column_id = paste(.data$cohort, .data$group, sep = "__")
  )

display_rows <- c(
  "N",
  "Age at death, mean (SD)",
  "Female, n (%)",
  "PMI, median [IQR]",
  "RIN, mean (SD)",
  "Braak, median [IQR]",
  "CERAD, median [IQR]",
  "APOE e4 carrier, n/N (%)",
  "APOE missing, n"
)

display_long <- purrr::map_dfr(
  seq_len(nrow(display_groups)),
  function(i) {
    cohort_i <- display_groups$cohort[i]
    group_i <- display_groups$group[i]
    column_i <- display_groups$column_id[i]

    get_cont <- function(v) {
      continuous_summary |>
        dplyr::filter(
          .data$cohort == cohort_i,
          .data$group == group_i,
          .data$variable == v
        )
    }

    sex_i <- sex_summary |>
      dplyr::filter(
        .data$cohort == cohort_i,
        .data$group == group_i
      )

    apoe_i <- apoe_summary |>
      dplyr::filter(
        .data$cohort == cohort_i,
        .data$group == group_i
      )

    n_i <- display_groups$n[i]

    age_i <- get_cont("age_death")
    pmi_i <- get_cont("pmi")
    rin_i <- get_cont("rin")
    braak_i <- get_cont("braak")
    cerad_i <- get_cont("cerad")

    female_string <- if (
      nrow(sex_i) == 1L &&
      sex_i$n_sex_available > 0L
    ) {
      sprintf(
        "%d (%.1f%%)",
        sex_i$n_female,
        sex_i$pct_female_among_available
      )
    } else {
      "NA"
    }

    apoe_carrier_string <- if (
      nrow(apoe_i) == 1L &&
      !is.na(apoe_i$n_apoe_available) &&
      apoe_i$n_apoe_available > 0L
    ) {
      sprintf(
        "%d/%d (%.1f%%)",
        apoe_i$n_e4_carrier,
        apoe_i$n_apoe_available,
        apoe_i$pct_e4_carrier_among_available
      )
    } else {
      "NA"
    }

    tibble::tibble(
      column_id = column_i,
      row = display_rows,
      value = c(
        as.character(n_i),
        mean_sd_string(age_i),
        female_string,
        median_iqr_string(pmi_i),
        mean_sd_string(rin_i),
        median_iqr_string(braak_i),
        median_iqr_string(cerad_i),
        apoe_carrier_string,
        if (nrow(apoe_i) == 1L) {
          as.character(apoe_i$n_apoe_missing)
        } else {
          "NA"
        }
      )
    )
  }
)

display_table <- display_long |>
  tidyr::pivot_wider(
    names_from = .data$column_id,
    values_from = .data$value
  )

############################################################
## 14. Participant manifests for reproducibility
############################################################

rna_manifest <- rna_demo |>
  dplyr::mutate(cohort = "RNA_corrected_577", .before = 1)

protein_full_manifest <- protein_demo_all |>
  dplyr::mutate(cohort = "Protein_full_400", .before = 1)

protein_primary_manifest <- protein_primary |>
  dplyr::mutate(cohort = "Protein_primary_374", .before = 1)

protein_complete_manifest <- protein_stage_complete |>
  dplyr::mutate(
    cohort = "Protein_stage_nuisance_complete_372",
    .before = 1
  )

regional_manifest <- regional_demo |>
  dplyr::mutate(
    cohort = "Regional_paired_DLPFC_STG_215",
    .before = 1
  )

############################################################
## 15. Write outputs
############################################################

outputs <- list(
  "demographics_cohort_counts.csv" = cohort_counts,
  "demographics_stage_group_counts.csv" = stage_counts,
  "demographics_continuous_summary.csv" = continuous_summary,
  "demographics_sex_summary.csv" = sex_summary,
  "demographics_APOE_summary.csv" = apoe_summary,
  "demographics_regional_region_counts.csv" = regional_region_counts,
  "demographics_variable_availability.csv" = variable_availability,
  "demographics_manuscript_ready_table.csv" = display_table,
  "demographics_RNA_manifest.csv" = rna_manifest,
  "demographics_protein_full_manifest.csv" = protein_full_manifest,
  "demographics_protein_primary_manifest.csv" = protein_primary_manifest,
  "demographics_protein_stage_complete_manifest.csv" =
    protein_complete_manifest,
  "demographics_regional_paired_manifest.csv" = regional_manifest
)

for (nm in names(outputs)) {
  utils::write.csv(
    outputs[[nm]],
    file.path(OUTDIR, nm),
    row.names = FALSE
  )
}

############################################################
## 16. Provenance
############################################################

git_head <- tryCatch(
  system2(
    "git",
    c("rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) NA_character_
)

provenance <- tibble::tibble(
  item = c(
    "analysis",
    "RNA_cohort",
    "RNA_stage_counts",
    "protein_full_cohort",
    "protein_primary_cohort",
    "protein_primary_stage_counts",
    "protein_stage_complete_cohort",
    "protein_stage_complete_counts",
    "regional_paired_cohort",
    "regional_regions",
    "APOE_source",
    "education_policy",
    "regional_CERAD_policy",
    "regional_APOE_policy",
    "protein_nuisance_covariates",
    "git_HEAD"
  ),
  value = c(
    "Canonical reviewer demographic tables",
    "577",
    "NCI=200;MCI=158;AD=219",
    "400",
    "374",
    "NCI=168;MCI=97;AD=109",
    "372",
    "NCI=167;MCI=96;AD=109",
    "215",
    "DLPFC=215;STG=215;paired=TRUE",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_participant_lookup.csv",
    "Education not inferred because absent from canonical main metadata/regional manifest.",
    "bScore is not relabeled as CERAD without explicit source support.",
    "Regional APOE unavailable: regional participant-ID namespace has 0/215 overlap with canonical APOE lookup; not treated as APOE missingness.",
    paste(protein_covars, collapse = ";"),
    paste(git_head, collapse = ";")
  )
)

utils::write.csv(
  provenance,
  file.path(OUTDIR, "demographics_provenance.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "demographics_sessionInfo.txt")
)

############################################################
## 17. Console report
############################################################

cat("\n============================================================\n")
cat("CANONICAL DEMOGRAPHIC TABLES\n")
cat("============================================================\n\n")

cat("Core cohort counts:\n")
print(
  cohort_counts |>
    dplyr::filter(
      .data$group %in% c("Overall", "NCI", "MCI", "AD")
    ),
  n = Inf,
  width = Inf
)

cat("\nAPOE summary for core overall cohorts:\n")
print(
  apoe_summary |>
    dplyr::filter(.data$group == "Overall"),
  n = Inf,
  width = Inf
)

cat("\nRegional brain-region coverage:\n")
print(regional_region_counts, n = Inf, width = Inf)

cat("\nOutputs written to:\n", OUTDIR, "\n", sep = "")
cat("\n============================================================\n")
cat("DEMOGRAPHIC TABLES BUILT — VALIDATION STILL REQUIRED\n")
cat("============================================================\n")
