############################################################
## 58_run_remaining_analysis_gates.R
##
## Consolidated audit of the remaining reviewer-analysis gates:
##   A. Main ROSMAP Braak/CERAD categorical-vs-continuous sensitivity
##   B. Figure 1B numerical/scale audit
##   C. Matched-cohort + multiple-testing/FDR consistency audit
##   D. RNA-protein subtraction inventory/policy audit
##
## This script is analytical only. It does not edit the manuscript,
## figures, or reviewer-response text.
############################################################

options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R",
  "R/04_build_pathway_sets_all_clients.R"
)) {
  if (!file.exists(f)) stop("Missing required production script: ", f)
  source(f)
}

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "remaining_analysis_gates"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

tol <- 1e-10

max_abs <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(0)
  max(abs(x))
}

pick_first <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

## =========================================================
## A. BRAAK / CERAD: CATEGORICAL VS CONTINUOUS
## =========================================================

path_file <- file.path(
  "outputs",
  "main_figures",
  "tables",
  "main_fig3_all_clients",
  "Fig3_all_clients_pathology_table.csv"
)

if (!file.exists(path_file)) stop("Missing canonical pathology table: ", path_file)

path_tbl <- read.csv(
  path_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(path_tbl) != 306L) {
  stop("Expected 306 Hsp60/10 clients in canonical pathology table.")
}

genes <- as.character(path_tbl$gene)

if (!all(genes %in% rownames(prot_mat_raw))) {
  stop("Not all canonical pathology genes are present in prot_mat_raw.")
}

if (!identical(protein_covars, c(
  "age_num", "sex_factor", "pmi_num", "batch_factor"
))) {
  stop("Canonical protein covariates changed.")
}

idx <- match(colnames(prot_mat_raw), prot_meta_adj$SampleID)
if (anyNA(idx)) stop("Could not align prot_mat_raw to prot_meta_adj.")
meta <- prot_meta_adj[idx, , drop = FALSE]

fit_factor_sensitivity <- function(
  predictor_raw,
  predictor_std,
  endpoint_label,
  primary_beta_col,
  primary_p_col,
  primary_fdr_col
) {
  raw <- suppressWarnings(as.numeric(as.character(meta[[predictor_raw]])))
  std <- suppressWarnings(as.numeric(as.character(meta[[predictor_std]])))

  levs <- sort(unique(raw[is.finite(raw)]))
  if (length(levs) < 3L) {
    stop(endpoint_label, " has fewer than three observed categories.")
  }

  counts <- as.data.frame(table(
    factor(raw, levels = levs),
    useNA = "no"
  ))
  names(counts) <- c("numeric_category", "n")
  counts$endpoint <- endpoint_label
  counts$numeric_category <- as.character(counts$numeric_category)

  results <- lapply(genes, function(g) {
    y <- safe_z(as.numeric(prot_mat_raw[g, ]))

    d <- data.frame(
      abundance = y,
      pathology_raw = raw,
      pathology_std = std,
      meta[, protein_covars, drop = FALSE],
      check.names = FALSE
    )
    d <- d[complete.cases(d), , drop = FALSE]

    if (nrow(d) < 30L || length(unique(d$pathology_raw)) < 3L) {
      return(data.frame(
        gene = g, n = nrow(d),
        continuous_beta = NA_real_,
        continuous_p = NA_real_,
        categorical_omnibus_p = NA_real_,
        high_vs_low_numeric_effect = NA_real_,
        high_vs_low_numeric_se = NA_real_,
        high_vs_low_numeric_ci_low = NA_real_,
        high_vs_low_numeric_ci_high = NA_real_,
        nonlinearity_p = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    d$pathology_factor <- factor(
      d$pathology_raw,
      levels = levs
    )

    continuous_formula <- as.formula(
      paste(
        "abundance ~ pathology_std +",
        paste(protein_covars, collapse = " + ")
      )
    )

    categorical_formula <- as.formula(
      paste(
        "abundance ~ pathology_factor +",
        paste(protein_covars, collapse = " + ")
      )
    )

    reduced_formula <- as.formula(
      paste(
        "abundance ~",
        paste(protein_covars, collapse = " + ")
      )
    )

    fit_cont <- lm(continuous_formula, data = d)
    fit_cat <- lm(categorical_formula, data = d)
    fit_red <- lm(reduced_formula, data = d)

    tt_cont <- summary(fit_cont)$coefficients
    beta <- unname(tt_cont["pathology_std", "Estimate"])
    p_cont <- unname(tt_cont["pathology_std", "Pr(>|t|)"])

    omnibus <- anova(fit_red, fit_cat)
    p_omnibus <- omnibus$`Pr(>F)`[2]

    ## Exact factor model is a superset of a linear trend model.
    nonlinear <- anova(fit_cont, fit_cat)
    p_nonlinear <- nonlinear$`Pr(>F)`[2]

    ## Lowest numeric category is the factor reference. Therefore the
    ## coefficient for the highest category is highest - lowest.
    highest <- tail(levs, 1)
    term_candidates <- paste0(
      "pathology_factor",
      c(
        as.character(highest),
        make.names(as.character(highest))
      )
    )
    term <- term_candidates[term_candidates %in% names(coef(fit_cat))][1]

    if (is.na(term)) {
      ## Robust fallback: identify the non-reference factor coefficient
      ## whose suffix corresponds to the highest numeric level.
      cn <- names(coef(fit_cat))
      cand <- cn[grepl("^pathology_factor", cn)]
      term <- cand[
        grepl(
          paste0(gsub("\\.", "\\\\.", as.character(highest)), "$"),
          cand
        )
      ][1]
    }

    if (is.na(term) || !nzchar(term)) {
      stop(
        "Could not identify highest-vs-lowest factor contrast for ",
        endpoint_label, " / ", g
      )
    }

    tt_cat <- summary(fit_cat)$coefficients
    eff <- unname(tt_cat[term, "Estimate"])
    se <- unname(tt_cat[term, "Std. Error"])
    df <- df.residual(fit_cat)
    crit <- qt(0.975, df)

    data.frame(
      gene = g,
      n = nrow(d),
      continuous_beta = beta,
      continuous_p = p_cont,
      categorical_omnibus_p = p_omnibus,
      high_vs_low_numeric_effect = eff,
      high_vs_low_numeric_se = se,
      high_vs_low_numeric_ci_low = eff - crit * se,
      high_vs_low_numeric_ci_high = eff + crit * se,
      nonlinearity_p = p_nonlinear,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results)
  results$endpoint <- endpoint_label
  results$continuous_fdr_bh <- p.adjust(
    results$continuous_p,
    method = "BH"
  )
  results$categorical_omnibus_fdr_bh <- p.adjust(
    results$categorical_omnibus_p,
    method = "BH"
  )
  results$nonlinearity_fdr_bh <- p.adjust(
    results$nonlinearity_p,
    method = "BH"
  )

  primary <- path_tbl[
    match(results$gene, path_tbl$gene),
    c(
      "gene",
      primary_beta_col,
      primary_p_col,
      primary_fdr_col
    ),
    drop = FALSE
  ]

  names(primary)[2:4] <- c(
    "production_beta",
    "production_p",
    "production_fdr"
  )

  results <- merge(
    results,
    primary,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  results$continuous_effect_difference_vs_production <- (
    results$continuous_beta - results$production_beta
  )
  results$continuous_p_difference_vs_production <- (
    results$continuous_p - results$production_p
  )
  results$continuous_fdr_difference_vs_production <- (
    results$continuous_fdr_bh - results$production_fdr
  )

  finite <- is.finite(results$continuous_beta) &
    is.finite(results$high_vs_low_numeric_effect)

  summary <- data.frame(
    endpoint = endpoint_label,
    n_genes = nrow(results),
    numeric_levels = paste(levs, collapse = ";"),
    n_continuous_fdr05 = sum(
      results$continuous_fdr_bh < 0.05,
      na.rm = TRUE
    ),
    n_categorical_omnibus_fdr05 = sum(
      results$categorical_omnibus_fdr_bh < 0.05,
      na.rm = TRUE
    ),
    n_nonlinearity_fdr05 = sum(
      results$nonlinearity_fdr_bh < 0.05,
      na.rm = TRUE
    ),
    pearson_continuous_vs_high_low = cor(
      results$continuous_beta[finite],
      results$high_vs_low_numeric_effect[finite],
      method = "pearson"
    ),
    spearman_continuous_vs_high_low = cor(
      results$continuous_beta[finite],
      results$high_vs_low_numeric_effect[finite],
      method = "spearman"
    ),
    sign_concordance_continuous_vs_high_low = mean(
      sign(results$continuous_beta[finite]) ==
        sign(results$high_vs_low_numeric_effect[finite])
    ),
    max_abs_beta_difference_vs_production = max_abs(
      results$continuous_effect_difference_vs_production
    ),
    max_abs_p_difference_vs_production = max_abs(
      results$continuous_p_difference_vs_production
    ),
    max_abs_fdr_difference_vs_production = max_abs(
      results$continuous_fdr_difference_vs_production
    ),
    stringsAsFactors = FALSE
  )

  list(
    counts = counts,
    results = results,
    summary = summary
  )
}

braak_sens <- fit_factor_sensitivity(
  predictor_raw = "braak_num",
  predictor_std = "braak_num_std",
  endpoint_label = "Braak",
  primary_beta_col = "braak_beta",
  primary_p_col = "braak_p",
  primary_fdr_col = "braak_padj"
)

cerad_sens <- fit_factor_sensitivity(
  predictor_raw = "cerad_num",
  predictor_std = "cerad_num_std",
  endpoint_label = "CERAD",
  primary_beta_col = "cerad_beta",
  primary_p_col = "cerad_p",
  primary_fdr_col = "cerad_padj"
)

category_counts <- rbind(
  braak_sens$counts,
  cerad_sens$counts
)

categorical_results <- rbind(
  braak_sens$results,
  cerad_sens$results
)

categorical_summary <- rbind(
  braak_sens$summary,
  cerad_sens$summary
)

## Braak 3-bin descriptive sensitivity if the canonical levels are 1:6.
braak_levels <- sort(unique(
  suppressWarnings(as.numeric(as.character(meta$braak_num)))
))
braak_levels <- braak_levels[is.finite(braak_levels)]

braak_bin3_summary <- data.frame(
  available = FALSE,
  level_I_II_n = NA_integer_,
  level_III_IV_n = NA_integer_,
  level_V_VI_n = NA_integer_,
  note = "Not calculated because Braak numeric categories were not exactly 1:6.",
  stringsAsFactors = FALSE
)

if (identical(braak_levels, 1:6)) {
  b <- suppressWarnings(as.numeric(as.character(meta$braak_num)))
  bin <- cut(
    b,
    breaks = c(0.5, 2.5, 4.5, 6.5),
    labels = c("I-II", "III-IV", "V-VI"),
    include.lowest = TRUE
  )
  tb <- table(bin, useNA = "no")

  braak_bin3_summary <- data.frame(
    available = TRUE,
    level_I_II_n = unname(tb["I-II"]),
    level_III_IV_n = unname(tb["III-IV"]),
    level_V_VI_n = unname(tb["V-VI"]),
    note = "Three-bin sensitivity uses canonical numeric Braak 1-2 / 3-4 / 5-6.",
    stringsAsFactors = FALSE
  )
}

## =========================================================
## B. FIGURE 1B NUMERICAL / SCALE AUDIT
## =========================================================

fig1_summary_file <- file.path(
  "outputs", "main_figures", "tables",
  "main_fig1_all_clients",
  "Fig1B_all_clients_clean_effect_summary.csv"
)

fig1_boot_file <- file.path(
  "outputs", "main_figures", "tables",
  "main_fig1_all_clients",
  "Fig1B_all_clients_bootstrap_significance_FIXED.csv"
)

fig1_minus_file <- file.path(
  "outputs", "main_figures", "tables",
  "main_fig1_all_clients",
  "Fig1B_all_clients_protein_minus_rna_effect_difference.csv"
)

for (f in c(fig1_summary_file, fig1_boot_file, fig1_minus_file)) {
  if (!file.exists(f)) stop("Missing Fig1B source file: ", f)
}

fig1 <- read.csv(
  fig1_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
fig1_boot <- read.csv(
  fig1_boot_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
fig1_minus <- read.csv(
  fig1_minus_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_fig1 <- c(
  "Pathway", "Shift",
  "Protein_NCI", "Protein_MCI", "Protein_AD",
  "RNA_NCI", "RNA_MCI", "RNA_AD",
  "protein_early_effect", "protein_late_effect",
  "rna_early_effect", "rna_late_effect",
  "protein_minus_rna_magnitude"
)
if (!all(required_fig1 %in% names(fig1))) {
  stop("Fig1B clean effect summary schema changed.")
}

fig1$protein_early_recalc <- fig1$Protein_MCI - fig1$Protein_NCI
fig1$protein_late_recalc <- fig1$Protein_AD - fig1$Protein_MCI
fig1$rna_early_recalc <- fig1$RNA_MCI - fig1$RNA_NCI
fig1$rna_late_recalc <- fig1$RNA_AD - fig1$RNA_MCI

fig1$expected_magnitude_difference <- ifelse(
  fig1$Shift == "NCI to MCI",
  abs(fig1$protein_early_effect) - abs(fig1$rna_early_effect),
  ifelse(
    fig1$Shift == "MCI to AD",
    abs(fig1$protein_late_effect) - abs(fig1$rna_late_effect),
    NA_real_
  )
)

fig1$arithmetic_difference <- (
  fig1$protein_minus_rna_magnitude -
    fig1$expected_magnitude_difference
)

boot_key <- paste(fig1_boot$Pathway, fig1_boot$Shift, sep = "||")
fig_key <- paste(fig1$Pathway, fig1$Shift, sep = "||")
boot_idx <- match(fig_key, boot_key)

if (anyNA(boot_idx)) {
  stop("Could not align Fig1B clean summary to bootstrap table.")
}

fig1$bootstrap_delta <- fig1_boot$delta[boot_idx]
fig1$bootstrap_ci_lo <- fig1_boot$ci_lo[boot_idx]
fig1$bootstrap_ci_hi <- fig1_boot$ci_hi[boot_idx]
fig1$bootstrap_p <- fig1_boot$p_value[boot_idx]
fig1$bootstrap_minus_point_difference <- (
  fig1$bootstrap_delta -
    fig1$protein_minus_rna_magnitude
)

fig1_scale_audit <- data.frame(
  check = c(
    "protein_early_effect_arithmetic",
    "protein_late_effect_arithmetic",
    "RNA_early_effect_arithmetic",
    "RNA_late_effect_arithmetic",
    "magnitude_difference_arithmetic",
    "bootstrap_rows_align",
    "bootstrap_CI_ordered",
    "score_scale_policy"
  ),
  max_abs_difference = c(
    max_abs(fig1$protein_early_effect - fig1$protein_early_recalc),
    max_abs(fig1$protein_late_effect - fig1$protein_late_recalc),
    max_abs(fig1$rna_early_effect - fig1$rna_early_recalc),
    max_abs(fig1$rna_late_effect - fig1$rna_late_recalc),
    max_abs(fig1$arithmetic_difference),
    0,
    0,
    NA_real_
  ),
  passed = c(
    max_abs(fig1$protein_early_effect - fig1$protein_early_recalc) < tol,
    max_abs(fig1$protein_late_effect - fig1$protein_late_recalc) < tol,
    max_abs(fig1$rna_early_effect - fig1$rna_early_recalc) < tol,
    max_abs(fig1$rna_late_effect - fig1$rna_late_recalc) < tol,
    max_abs(fig1$arithmetic_difference) < tol,
    all(!is.na(boot_idx)),
    all(fig1$bootstrap_ci_lo <= fig1$bootstrap_ci_hi),
    TRUE
  ),
  note = c(
    "Protein NCI->MCI equals MCI mean minus NCI mean.",
    "Protein MCI->AD equals AD mean minus MCI mean.",
    "RNA NCI->MCI equals MCI mean minus NCI mean.",
    "RNA MCI->AD equals AD mean minus MCI mean.",
    "Stored protein_minus_rna_magnitude equals |protein transition| - |RNA transition|.",
    "Bootstrap table aligns one-to-one by Pathway + Shift.",
    "Bootstrap confidence intervals are numerically ordered.",
    paste(
      "RNA and protein pathway scores are treated as modality-specific",
      "dimensionless z-score-derived pathway summaries.",
      "The cross-modal magnitude difference is therefore descriptive,",
      "not a molecular-scale fold-change comparison."
    )
  ),
  stringsAsFactors = FALSE
)

fig1_scale_summary <- data.frame(
  n_rows = nrow(fig1),
  max_abs_point_arithmetic_error = max_abs(fig1$arithmetic_difference),
  max_abs_bootstrap_delta_minus_point = max_abs(
    fig1$bootstrap_minus_point_difference
  ),
  mean_abs_bootstrap_delta_minus_point = mean(
    abs(fig1$bootstrap_minus_point_difference),
    na.rm = TRUE
  ),
  interpretation = paste(
    "Fig1B arithmetic is internally consistent.",
    "Protein-vs-RNA magnitude differences are descriptive comparisons",
    "of standardized pathway-score changes and should not be interpreted",
    "as subtraction of directly commensurate molecular abundance scales."
  ),
  stringsAsFactors = FALSE
)

## =========================================================
## C1. CANONICAL MATCHED RNA/PROTEIN COHORT
## =========================================================
##
## IMPORTANT:
## The reviewer-facing matched-individual sensitivity is
## Supplementary Figure 2. Its canonical cohort is the exact
## RNA/protein participant intersection after matrix-backed
## primary-stage filtering, with identical canonical stage
## assignments across modalities.
##
## Do NOT use the older stage_harmonization "strict_pair"
## manifest here; that artifact represents a different,
## more restrictive audit definition.

matched_file <- file.path(
  "outputs",
  "supplemental_figures",
  "audits",
  "SuppFig2_matched_stage_crosswalk.csv"
)

if (!file.exists(matched_file)) {
  stop(
    "Missing canonical Supplementary Figure 2 matched crosswalk: ",
    matched_file,
    ". Run the validated supplemental pipeline first."
  )
}

matched <- read.csv(
  matched_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_matched_cols <- c(
  "individual_id",
  "protein_sample_id",
  "protein_stage",
  "rna_sample_id",
  "rna_stage",
  "stage_agrees"
)

missing_matched_cols <- setdiff(
  required_matched_cols,
  names(matched)
)

if (length(missing_matched_cols)) {
  stop(
    "Canonical SuppFig2 matched crosswalk is missing columns: ",
    paste(missing_matched_cols, collapse = ", ")
  )
}

matched$individual_id <- as.character(
  matched$individual_id
)
matched$protein_sample_id <- as.character(
  matched$protein_sample_id
)
matched$rna_sample_id <- as.character(
  matched$rna_sample_id
)

protein_stage <- toupper(
  trimws(as.character(matched$protein_stage))
)
rna_stage <- toupper(
  trimws(as.character(matched$rna_stage))
)

if (anyDuplicated(matched$individual_id)) {
  stop(
    "Canonical SuppFig2 matched crosswalk contains duplicate participants."
  )
}

if (anyDuplicated(matched$protein_sample_id)) {
  stop(
    "Canonical SuppFig2 matched crosswalk contains duplicate protein samples."
  )
}

if (anyDuplicated(matched$rna_sample_id)) {
  stop(
    "Canonical SuppFig2 matched crosswalk contains duplicate RNA samples."
  )
}

if (!all(matched$stage_agrees %in% TRUE)) {
  stop(
    "Canonical SuppFig2 matched crosswalk contains RNA/protein stage disagreement."
  )
}

if (!all(protein_stage == rna_stage)) {
  stop(
    "Canonical SuppFig2 protein_stage and rna_stage differ."
  )
}

allowed_stage <- c("NCI", "MCI", "AD")

if (!all(protein_stage %in% allowed_stage)) {
  stop(
    "Canonical SuppFig2 contains an unexpected clinical stage: ",
    paste(
      setdiff(unique(protein_stage), allowed_stage),
      collapse = ", "
    )
  )
}

matched_counts <- as.data.frame(
  table(
    factor(
      protein_stage,
      levels = allowed_stage
    )
  )
)

names(matched_counts) <- c(
  "clinical_stage",
  "n"
)

matched_summary <- data.frame(
  manifest_file = matched_file,
  cohort_definition = paste(
    "Exact participant intersection with matrix-backed RNA and protein",
    "and identical canonical NCI/MCI/AD stage assignments"
  ),
  participant_key = "individual_id",
  n_rows = nrow(matched),
  n_unique_participants = length(
    unique(matched$individual_id)
  ),
  n_unique_RNA_samples = length(
    unique(matched$rna_sample_id)
  ),
  n_unique_protein_samples = length(
    unique(matched$protein_sample_id)
  ),
  duplicate_participant_rows = sum(
    duplicated(matched$individual_id)
  ),
  stage_disagreements = sum(
    protein_stage != rna_stage
  ),
  n_NCI = matched_counts$n[
    matched_counts$clinical_stage == "NCI"
  ],
  n_MCI = matched_counts$n[
    matched_counts$clinical_stage == "MCI"
  ],
  n_AD = matched_counts$n[
    matched_counts$clinical_stage == "AD"
  ],
  expected_n = 198L,
  expected_NCI = 94L,
  expected_MCI = 56L,
  expected_AD = 48L,
  counts_match_expected = (
    nrow(matched) == 198L &&
      length(unique(matched$individual_id)) == 198L &&
      length(unique(matched$rna_sample_id)) == 198L &&
      length(unique(matched$protein_sample_id)) == 198L &&
      matched_counts$n[
        matched_counts$clinical_stage == "NCI"
      ] == 94L &&
      matched_counts$n[
        matched_counts$clinical_stage == "MCI"
      ] == 56L &&
      matched_counts$n[
        matched_counts$clinical_stage == "AD"
      ] == 48L
  ),
  stringsAsFactors = FALSE
)

stage_consistency <- data.frame(
  metric = c(
    "participants",
    "unique_RNA_samples",
    "unique_protein_samples",
    "RNA_protein_stage_disagreements",
    "NCI",
    "MCI",
    "AD"
  ),
  value = c(
    nrow(matched),
    length(unique(matched$rna_sample_id)),
    length(unique(matched$protein_sample_id)),
    sum(protein_stage != rna_stage),
    matched_summary$n_NCI,
    matched_summary$n_MCI,
    matched_summary$n_AD
  ),
  expected = c(
    198L,
    198L,
    198L,
    0L,
    94L,
    56L,
    48L
  ),
  passed = c(
    nrow(matched) == 198L,
    length(unique(matched$rna_sample_id)) == 198L,
    length(unique(matched$protein_sample_id)) == 198L,
    sum(protein_stage != rna_stage) == 0L,
    matched_summary$n_NCI == 94L,
    matched_summary$n_MCI == 56L,
    matched_summary$n_AD == 48L
  ),
  stringsAsFactors = FALSE
)

## C2. FDR / MULTIPLE-TESTING NUMERICAL CONSISTENCY
## =========================================================

fdr_rows <- list()

add_fdr_audit <- function(
  family,
  file,
  p_col,
  fdr_col,
  group_col = NULL
) {
  if (!file.exists(file)) {
    stop("Missing FDR audit file: ", file)
  }

  d <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c(p_col, fdr_col) %in% names(d))) {
    stop(
      "Missing p/FDR columns for ", family, ": ",
      p_col, " / ", fdr_col, " in ", file
    )
  }

  if (is.null(group_col)) {
    recomputed <- p.adjust(d[[p_col]], method = "BH")
    maxdiff <- max_abs(recomputed - d[[fdr_col]])
    groups <- 1L
  } else {
    if (!group_col %in% names(d)) {
      stop("Missing grouping column ", group_col, " in ", file)
    }

    recomputed <- rep(NA_real_, nrow(d))
    lev <- unique(as.character(d[[group_col]]))

    for (g in lev) {
      ii <- which(as.character(d[[group_col]]) == g)
      recomputed[ii] <- p.adjust(d[[p_col]][ii], method = "BH")
    }

    maxdiff <- max_abs(recomputed - d[[fdr_col]])
    groups <- length(lev)
  }

  fdr_rows[[length(fdr_rows) + 1L]] <<- data.frame(
    family = family,
    file = file,
    p_column = p_col,
    fdr_column = fdr_col,
    group_column = ifelse(is.null(group_col), "", group_col),
    n_rows = nrow(d),
    n_groups = groups,
    max_abs_fdr_difference = maxdiff,
    passed = maxdiff < tol,
    stringsAsFactors = FALSE
  )
}

## Pathology: endpoint-specific BH across 306 clients.
path_resolution <- file.path(
  "outputs", "reviewer_revisions",
  "reporting_requirements_resolution",
  "02_pathology_Braak_CERAD_effect_SE_CI_n_P_FDR.csv"
)

add_fdr_audit(
  "Hsp client pathology models",
  path_resolution,
  "p_value",
  "fdr_bh",
  "endpoint"
)

## Network cognition: BH across three prespecified outcomes.
network_resolution <- file.path(
  "outputs", "reviewer_revisions",
  "reporting_requirements_resolution",
  "06_network_cognition_reporting_complete.csv"
)
add_fdr_audit(
  "Network cognition outcomes",
  network_resolution,
  "p.value",
  "fdr_bh_three_network_outcomes"
)

## Matched/null metrics.
null_resolution <- file.path(
  "outputs", "reviewer_revisions",
  "reporting_requirements_resolution",
  "08_cognition_null_reporting_complete.csv"
)
add_fdr_audit(
  "Corrected cognition null metrics",
  null_resolution,
  "empirical_p_greater",
  "fdr_bh_across_null_metrics"
)

## Regional Hsp two primary tests.
regional_hsp <- file.path(
  "outputs", "reviewer_revisions",
  "regional_formal_interaction", "hsp_conclusion_audit",
  "Hsp_pathway_score_primary_interactions.csv"
)
add_fdr_audit(
  "Regional Hsp pathway primary interactions",
  regional_hsp,
  "p_region_x_predictor",
  "BH_FDR_two_primary_tests"
)

## Regional genome-wide AD and continuous Braak.
regional_ad <- file.path(
  "outputs", "reviewer_revisions",
  "regional_formal_interaction", "validated_final",
  "FINAL_region_x_AD_interaction_all_proteins.csv"
)
regional_braak <- file.path(
  "outputs", "reviewer_revisions",
  "regional_formal_interaction", "validated_final",
  "FINAL_region_x_Braak_continuous_interaction_all_proteins.csv"
)
regional_bin3 <- file.path(
  "outputs", "reviewer_revisions",
  "regional_formal_interaction", "validated_final",
  "FINAL_region_x_Braak_bin3_interaction_all_proteins.csv"
)

add_fdr_audit(
  "Regional genome-wide AD interaction",
  regional_ad,
  "p_region_x_predictor",
  "p_region_x_predictor_fdr"
)
add_fdr_audit(
  "Regional genome-wide continuous Braak interaction",
  regional_braak,
  "p_region_x_predictor",
  "p_region_x_predictor_fdr"
)
add_fdr_audit(
  "Regional genome-wide categorical Braak interaction",
  regional_bin3,
  "global_region_x_braak_bin3_p",
  "global_region_x_braak_bin3_fdr"
)

## PC1 sensitivity: two contrasts within each score type.
pc1_file <- file.path(
  "outputs", "reviewer_revisions",
  "Hsp_pathway_PC1_sensitivity",
  "stage_Wilcoxon_mean_z_vs_PC1.csv"
)
add_fdr_audit(
  "RNA pathway PC1 sensitivity",
  pc1_file,
  "p_value",
  "p_fdr_bh_within_score",
  "score_type"
)

## Alternative-mechanism panel: 39 markers within each contrast.
marker_file <- file.path(
  "outputs", "reviewer_revisions",
  "alternative_mechanism_marker_panel",
  "marker_stage_limma_results.csv"
)
if (file.exists(marker_file)) {
  md <- read.csv(marker_file, nrows = 2, check.names = FALSE)
  p_marker <- pick_first(
    md,
    c("p_value", "p.value", "p", "P.Value")
  )
  f_marker <- pick_first(
    md,
    c("marker_panel_fdr", "fdr", "adj.P.Val")
  )
  contrast_marker <- pick_first(
    md,
    c("contrast", "comparison")
  )

  if (
    !is.na(p_marker) &&
    !is.na(f_marker) &&
    !is.na(contrast_marker)
  ) {
    add_fdr_audit(
      "Alternative-mechanism marker panel",
      marker_file,
      p_marker,
      f_marker,
      contrast_marker
    )
  }
}

fdr_audit <- do.call(rbind, fdr_rows)

## =========================================================
## D. RNA-PROTEIN SUBTRACTION INVENTORY / POLICY
## =========================================================

r_files <- list.files(
  "R",
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

subtraction_patterns <- c(
  "protein_minus_rna",
  "rna_minus_protein",
  "minus_rna",
  "minus_protein",
  "subtraction"
)

scan_rows <- list()

for (f in r_files) {
  txt <- readLines(f, warn = FALSE)
  hit <- which(vapply(
    txt,
    function(line) {
      any(vapply(
        subtraction_patterns,
        function(pat) grepl(
          pat,
          tolower(line),
          fixed = TRUE
        ),
        logical(1)
      ))
    },
    logical(1)
  ))

  if (length(hit)) {
    scan_rows[[length(scan_rows) + 1L]] <- data.frame(
      file = f,
      line = hit,
      text = trimws(txt[hit]),
      stringsAsFactors = FALSE
    )
  }
}

if (length(scan_rows)) {
  subtraction_code_inventory <- do.call(rbind, scan_rows)
} else {
  subtraction_code_inventory <- data.frame(
    file = character(),
    line = integer(),
    text = character(),
    stringsAsFactors = FALSE
  )
}

output_candidates <- list.files(
  "outputs",
  pattern = "protein_minus_rna|rna_minus_protein|NO_SUBTRACTION",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

subtraction_output_inventory <- data.frame(
  file = output_candidates,
  classification = ifelse(
    grepl("NO_SUBTRACTION", output_candidates, fixed = TRUE),
    "preferred_side_by_side_no_subtraction",
    ifelse(
      grepl(
        "main_fig1_all_clients/Fig1B",
        output_candidates,
        fixed = TRUE
      ) &
        !grepl(
          "RNA_batch_canonicalization",
          output_candidates,
          fixed = TRUE
        ),
      "canonical_Fig1B_descriptive_subtraction_to_deemphasize",
      "legacy_or_snapshot_subtraction_related"
    )
  ),
  stringsAsFactors = FALSE
)

preferred_files <- c(
  "outputs/reviewer_revisions/conventional_stage_differential_final/16_PRIMARY_cross_modal_shared_features_side_by_side_NO_SUBTRACTION.csv",
  "outputs/reviewer_revisions/conventional_stage_differential_final/17_PRIMARY_cross_modal_Hsp60_10_side_by_side_NO_SUBTRACTION.csv",
  "outputs/reviewer_revisions/conventional_stage_differential_final/18_PRIMARY_cross_modal_Hsp60_10_summary_NO_SUBTRACTION.csv"
)

if (!all(file.exists(preferred_files))) {
  stop("One or more canonical NO_SUBTRACTION conventional-DE outputs are missing.")
}

subtraction_policy <- data.frame(
  analysis_component = c(
    "Conventional RNA-vs-protein differential analysis",
    "Figure 1B protein-minus-RNA magnitude difference",
    "Cross-modal biological interpretation"
  ),
  analytical_policy = c(
    "Use side-by-side modality-specific effects and FDR; no RNA-protein subtraction.",
    paste(
      "Arithmetic retained only as a descriptive standardized-score comparison;",
      "do not treat as a molecular-scale effect."
    ),
    paste(
      "Do not infer translation/protein-specific regulation from direct subtraction",
      "of RNA and protein effect magnitudes."
    )
  ),
  manuscript_action_later = c(
    "Use canonical NO_SUBTRACTION outputs.",
    "De-emphasize/remove inferential framing when manuscript/figure is edited.",
    "Describe discordance qualitatively or with side-by-side modality-specific estimates."
  ),
  stringsAsFactors = FALSE
)

## =========================================================
## FINAL GATE STATUS
## =========================================================

gate_status <- data.frame(
  gate = c(
    "Braak_CERAD_categorical_vs_continuous",
    "Fig1B_numerical_and_scale",
    "Matched_cohort_and_FDR",
    "Subtraction_analysis"
  ),
  analytical_status = c(
    ifelse(
      all(
        categorical_summary$max_abs_beta_difference_vs_production < tol,
        categorical_summary$max_abs_p_difference_vs_production < tol,
        categorical_summary$max_abs_fdr_difference_vs_production < tol
      ),
      "ANALYTICALLY_REVALIDATED",
      "FAILED"
    ),
    ifelse(
      all(fig1_scale_audit$passed),
      "ANALYTICALLY_REVALIDATED",
      "FAILED"
    ),
    ifelse(
      isTRUE(matched_summary$counts_match_expected) &&
        all(fdr_audit$passed),
      "ANALYTICALLY_REVALIDATED",
      "FAILED"
    ),
    "ANALYTICALLY_RESOLVED_MANUSCRIPT_EDIT_PENDING"
  ),
  note = c(
    paste(
      "Exact categorical-factor sensitivity completed for Braak and CERAD;",
      "continuous models reproduce production exactly."
    ),
    paste(
      "Transition arithmetic, cross-modal magnitude arithmetic, bootstrap alignment,",
      "and standardized-score interpretation audited."
    ),
    paste(
      "Strict matched cohort and major BH/FDR families audited numerically."
    ),
    paste(
      "Canonical conventional DE uses NO_SUBTRACTION outputs.",
      "Fig1B subtraction-derived descriptive metric remains to be de-emphasized/removed",
      "during manuscript/figure editing."
    )
  ),
  stringsAsFactors = FALSE
)

if (any(gate_status$analytical_status == "FAILED")) {
  print(gate_status)
  stop("One or more consolidated remaining analysis gates failed.")
}

## =========================================================
## WRITE OUTPUTS
## =========================================================

write.csv(
  category_counts,
  file.path(OUTDIR, "01_pathology_numeric_category_counts.csv"),
  row.names = FALSE
)
write.csv(
  categorical_results,
  file.path(OUTDIR, "02_pathology_categorical_vs_continuous_gene_results.csv"),
  row.names = FALSE
)
write.csv(
  categorical_summary,
  file.path(OUTDIR, "03_pathology_categorical_vs_continuous_summary.csv"),
  row.names = FALSE
)
write.csv(
  braak_bin3_summary,
  file.path(OUTDIR, "04_Braak_three_bin_summary.csv"),
  row.names = FALSE
)
write.csv(
  fig1,
  file.path(OUTDIR, "05_Fig1B_numerical_audit_rows.csv"),
  row.names = FALSE
)
write.csv(
  fig1_scale_audit,
  file.path(OUTDIR, "06_Fig1B_scale_audit.csv"),
  row.names = FALSE
)
write.csv(
  fig1_scale_summary,
  file.path(OUTDIR, "07_Fig1B_scale_summary.csv"),
  row.names = FALSE
)
write.csv(
  matched_summary,
  file.path(OUTDIR, "08_matched_cohort_summary.csv"),
  row.names = FALSE
)
write.csv(
  matched_counts,
  file.path(OUTDIR, "09_matched_cohort_stage_counts.csv"),
  row.names = FALSE
)
write.csv(
  stage_consistency,
  file.path(OUTDIR, "10_matched_cohort_stage_column_audit.csv"),
  row.names = FALSE
)
write.csv(
  fdr_audit,
  file.path(OUTDIR, "11_FDR_numerical_consistency_audit.csv"),
  row.names = FALSE
)
write.csv(
  subtraction_code_inventory,
  file.path(OUTDIR, "12_subtraction_code_inventory.csv"),
  row.names = FALSE
)
write.csv(
  subtraction_output_inventory,
  file.path(OUTDIR, "13_subtraction_output_inventory.csv"),
  row.names = FALSE
)
write.csv(
  subtraction_policy,
  file.path(OUTDIR, "14_subtraction_policy.csv"),
  row.names = FALSE
)
write.csv(
  gate_status,
  file.path(OUTDIR, "15_remaining_gate_status.csv"),
  row.names = FALSE
)

git_head <- tryCatch(
  system2(
    "git",
    c("rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) NA_character_
)

provenance <- data.frame(
  item = c(
    "analysis",
    "Braak_adjusted_for_CERAD",
    "protein_covariates",
    "Hsp_pathology_gene_n",
    "strict_matched_expected_n",
    "strict_matched_expected_stage_counts",
    "FDR_families_reaudited",
    "git_HEAD"
  ),
  value = c(
    "Consolidated remaining reviewer-analysis gates",
    as.character(cfg$adjust_braak_for_cerad),
    paste(protein_covars, collapse = ";"),
    "306",
    "198",
    "NCI=94;MCI=56;AD=48",
    as.character(nrow(fdr_audit)),
    paste(git_head, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(OUTDIR, "16_provenance.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "17_sessionInfo.txt")
)

## =========================================================
## CONSOLE
## =========================================================

cat("\n============================================================\n")
cat("CONSOLIDATED REMAINING ANALYSIS GATES\n")
cat("============================================================\n\n")

cat("PATHOLOGY CATEGORICAL-VS-CONTINUOUS SUMMARY\n")
print(categorical_summary, row.names = FALSE)

cat("\nBRAAK 3-BIN SUMMARY\n")
print(braak_bin3_summary, row.names = FALSE)

cat("\nFIG1B SCALE SUMMARY\n")
print(fig1_scale_summary, row.names = FALSE)

cat("\nSTRICT MATCHED COHORT\n")
print(matched_summary, row.names = FALSE)

cat("\nFDR NUMERICAL CONSISTENCY\n")
print(
  fdr_audit[
    ,
    c(
      "family", "n_rows", "n_groups",
      "max_abs_fdr_difference", "passed"
    )
  ],
  row.names = FALSE
)

cat("\nSUBTRACTION POLICY\n")
print(subtraction_policy, row.names = FALSE)

cat("\nFINAL GATE STATUS\n")
print(gate_status, row.names = FALSE)

cat("\nOutputs written to:\n", OUTDIR, "\n", sep = "")

cat("\n============================================================\n")
cat("CONSOLIDATED ANALYSIS GATES BUILT — VALIDATION REQUIRED\n")
cat("============================================================\n")
