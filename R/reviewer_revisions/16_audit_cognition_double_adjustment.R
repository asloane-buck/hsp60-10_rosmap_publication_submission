#!/usr/bin/env Rscript

############################################################
## 16_audit_cognition_double_adjustment.R
##
## Purpose
## -------
## Audit the remaining cognition-method blocker before changing any
## production Figure 4/5 scripts.
##
## OLD specification:
##   - abundance from prot_mat_adj (already residualized for
##     age + sex + PMI + TMT batch)
##   - cognition regression additionally includes
##     age + sex + education + PMI + Braak + CERAD
##
## CORRECTED specification:
##   - abundance from processed/non-residualized prot_mat_raw
##   - cognition regression directly includes
##     age + sex + education + PMI + Braak + CERAD + TMT batch
##
## This script:
##   1. hard-confirms the old double-adjustment structure;
##   2. compares old vs corrected network-level cognition models;
##   3. compares all 306 client-level cognition models;
##   4. compares Figure 5 cognition scores across the 915-gene null universe;
##   5. recomputes the 10,000-iteration cognition-only abundance-matched null
##      under both specifications using the same seed.
##
## It does NOT modify production scripts or production outputs.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

source("R/00_config.R")
source("R/01_utils.R")

out_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "cognition_double_adjustment_audit"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- function(name) {
  path <- file.path(cfg$object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing required production object: ", path, call. = FALSE)
  }
  readRDS(path)
}

prot_mat_adj <- obj("prot_mat_adj")
prot_mat_raw <- obj("prot_mat_raw")
prot_meta_adj <- obj("prot_meta_adj")
cognition_model_df <- obj("cognition_model_df")
all_hsp60_10_client_tbl <- obj("all_hsp60_10_client_tbl")

null_path <- file.path(
  cfg$table_dir,
  "figure5_null_input_COVARIATE_ADJUSTED.csv"
)
if (!file.exists(null_path)) {
  stop("Missing Figure 5 null input: ", null_path, call. = FALSE)
}
null_input_tbl <- readr::read_csv(null_path, show_col_types = FALSE)

############################################################
## 1. Structural invariants
############################################################

required_meta_cols <- c("SampleID", "batch_factor")
missing_meta <- setdiff(required_meta_cols, colnames(prot_meta_adj))
if (length(missing_meta) > 0) {
  stop(
    "prot_meta_adj missing required field(s): ",
    paste(missing_meta, collapse = ", "),
    call. = FALSE
  )
}

required_cognition_cols <- c(
  "SampleID", "individualID",
  "cogdx", "dcfdx_lv", "mmse_last_valid",
  "age_death", "educ", "sex", "pmi", "Braak", "CERAD",
  "Hsp60_client_score_z"
)
missing_cog <- setdiff(required_cognition_cols, colnames(cognition_model_df))
if (length(missing_cog) > 0) {
  stop(
    "cognition_model_df missing required field(s): ",
    paste(missing_cog, collapse = ", "),
    call. = FALSE
  )
}

if (!identical(colnames(prot_mat_adj), colnames(prot_mat_raw))) {
  stop("Adjusted/raw protein matrices do not have identical sample columns.", call. = FALSE)
}

canonicalize_matrix_rows <- function(mat, label) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  new_names <- canonical_gene_symbol(rownames(mat))

  if (anyDuplicated(new_names) > 0) {
    dup <- unique(new_names[duplicated(new_names)])
    stop(
      label,
      " has duplicate canonical gene names: ",
      paste(head(dup, 20), collapse = ", "),
      call. = FALSE
    )
  }

  rownames(mat) <- new_names
  mat
}

prot_mat_adj <- canonicalize_matrix_rows(prot_mat_adj, "prot_mat_adj")
prot_mat_raw <- canonicalize_matrix_rows(prot_mat_raw, "prot_mat_raw")

shared_samples <- intersect(colnames(prot_mat_adj), colnames(prot_mat_raw))
if (length(shared_samples) != 400) {
  stop(
    "Expected 400 shared protein samples; observed ",
    length(shared_samples),
    ".",
    call. = FALSE
  )
}

batch_tbl <- prot_meta_adj |>
  transmute(
    SampleID = as.character(.data$SampleID),
    batch_factor = as.factor(.data$batch_factor)
  ) |>
  distinct(.data$SampleID, .keep_all = TRUE)

cog <- cognition_model_df |>
  mutate(
    SampleID = as.character(.data$SampleID),
    individualID = as.character(.data$individualID)
  ) |>
  select(-any_of("batch_factor")) |>
  checked_left_join(
    batch_tbl,
    by = "SampleID",
    label = "cognition metadata to TMT batch"
  ) |>
  mutate(
    cogdx_severity = case_when(
      .data$cogdx == 1 ~ 0,
      .data$cogdx %in% c(2, 3) ~ 1,
      .data$cogdx %in% c(4, 5) ~ 2,
      TRUE ~ NA_real_
    ),
    dcfdx_lv_severity = case_when(
      .data$dcfdx_lv == 1 ~ 0,
      .data$dcfdx_lv %in% c(2, 3) ~ 1,
      .data$dcfdx_lv %in% c(4, 5) ~ 2,
      TRUE ~ NA_real_
    ),
    cogdx_better = case_when(
      .data$cogdx_severity == 0 ~ 2,
      .data$cogdx_severity == 1 ~ 1,
      .data$cogdx_severity == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    dcfdx_lv_better = case_when(
      .data$dcfdx_lv_severity == 0 ~ 2,
      .data$dcfdx_lv_severity == 1 ~ 1,
      .data$dcfdx_lv_severity == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    cogdx_better_z = as.numeric(scale(.data$cogdx_better)),
    dcfdx_lv_better_z = as.numeric(scale(.data$dcfdx_lv_better)),
    mmse_last_valid_z = as.numeric(scale(.data$mmse_last_valid)),
    cog_group = case_when(
      .data$cogdx == 1 ~ "NCI",
      .data$cogdx %in% c(2, 3) ~ "MCI",
      .data$cogdx %in% c(4, 5) ~ "AD",
      TRUE ~ NA_character_
    ),
    cog_group = factor(.data$cog_group, levels = c("NCI", "MCI", "AD")),
    sex = as.factor(.data$sex),
    batch_factor = as.factor(.data$batch_factor)
  )

if (sum(!is.na(cog$batch_factor)) != 400) {
  stop(
    "Expected batch_factor for all 400 protein samples; observed ",
    sum(!is.na(cog$batch_factor)),
    ".",
    call. = FALSE
  )
}

############################################################
## 2. Helpers
############################################################

row_z <- function(mat) {
  out <- t(scale(t(mat)))
  out[rowSums(is.finite(out)) > 0, , drop = FALSE]
}

z_or_na <- function(x) {
  if (sum(is.finite(x)) < 3 || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

run_one_lm <- function(data, outcome, predictor, corrected = FALSE) {
  covars <- c("age_death", "sex", "educ", "pmi", "Braak", "CERAD")
  if (corrected) covars <- c(covars, "batch_factor")

  keep <- c(outcome, predictor, covars)
  dat <- data |>
    select(all_of(keep)) |>
    drop_na()

  if (nrow(dat) < 40 ||
      sd(dat[[outcome]], na.rm = TRUE) == 0 ||
      sd(dat[[predictor]], na.rm = TRUE) == 0) {
    return(tibble(
      n = nrow(dat),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      model_status = "skipped_low_n_or_no_variance"
    ))
  }

  form <- reformulate(
    c(predictor, covars),
    response = outcome
  )

  fit <- tryCatch(lm(form, data = dat), error = function(e) e)

  if (inherits(fit, "error")) {
    return(tibble(
      n = nrow(dat),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      model_status = paste0("model_error: ", fit$message)
    ))
  }

  ans <- broom::tidy(fit, conf.int = TRUE) |>
    filter(.data$term == predictor)

  if (nrow(ans) != 1) {
    return(tibble(
      n = nobs(fit),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      model_status = "missing_predictor_term"
    ))
  }

  ans |>
    transmute(
      n = nobs(fit),
      estimate = .data$estimate,
      std.error = .data$std.error,
      statistic = .data$statistic,
      p.value = .data$p.value,
      conf.low = .data$conf.low,
      conf.high = .data$conf.high,
      model_status = "ok"
    )
}

run_gene_models <- function(mat, genes, spec_label, corrected) {
  genes <- intersect(genes, rownames(mat))
  mat_z <- row_z(mat[genes, , drop = FALSE])

  long <- as.data.frame(mat_z) |>
    rownames_to_column("gene") |>
    pivot_longer(
      cols = -.data$gene,
      names_to = "SampleID",
      values_to = "abundance_z"
    ) |>
    mutate(
      SampleID = as.character(.data$SampleID),
      abundance_z = as.numeric(.data$abundance_z)
    ) |>
    checked_left_join(
      cog |>
        select(
          .data$SampleID,
          .data$cogdx_better_z,
          .data$dcfdx_lv_better_z,
          .data$mmse_last_valid_z,
          .data$age_death,
          .data$sex,
          .data$educ,
          .data$pmi,
          .data$Braak,
          .data$CERAD,
          .data$batch_factor
        ),
      by = "SampleID",
      label = paste0(spec_label, " abundance to cognition metadata")
    )

  outcomes <- tribble(
    ~outcome, ~outcome_label,
    "cogdx_better_z", "Final cognitive diagnosis",
    "dcfdx_lv_better_z", "Last-valid cognitive diagnosis",
    "mmse_last_valid_z", "Last-valid MMSE"
  )

  res <- purrr::map_dfr(seq_len(nrow(outcomes)), function(i) {
    out <- outcomes$outcome[[i]]
    lab <- outcomes$outcome_label[[i]]

    message(
      spec_label, " gene cognition model: ", lab,
      " (", length(genes), " genes)"
    )

    long |>
      group_by(.data$gene) |>
      group_modify(~ {
        dat <- .x
        fit <- run_one_lm(
          data = dat,
          outcome = out,
          predictor = "abundance_z",
          corrected = corrected
        )
        fit
      }) |>
      ungroup() |>
      mutate(
        outcome = lab,
        specification = spec_label
      )
  }) |>
    group_by(.data$outcome) |>
    mutate(
      q.value = p.adjust(.data$p.value, method = "BH"),
      fdr_positive = .data$model_status == "ok" &
        .data$q.value < 0.05 &
        .data$estimate > 0,
      nominal_positive = .data$model_status == "ok" &
        .data$p.value < 0.05 &
        .data$estimate > 0
    ) |>
    ungroup()

  res
}

summarize_gene_priority <- function(model_tbl) {
  model_tbl |>
    filter(.data$model_status == "ok") |>
    group_by(.data$gene) |>
    summarise(
      n_outcomes_tested = n(),
      n_nominal_better_cognition =
        sum(.data$nominal_positive, na.rm = TRUE),
      n_fdr_better_cognition =
        sum(.data$fdr_positive, na.rm = TRUE),
      best_p = min(.data$p.value, na.rm = TRUE),
      best_q = min(.data$q.value, na.rm = TRUE),
      mean_beta = mean(.data$estimate, na.rm = TRUE),
      median_beta = median(.data$estimate, na.rm = TRUE),
      mmse_beta =
        .data$estimate[.data$outcome == "Last-valid MMSE"][1],
      mmse_p =
        .data$p.value[.data$outcome == "Last-valid MMSE"][1],
      mmse_q =
        .data$q.value[.data$outcome == "Last-valid MMSE"][1],
      .groups = "drop"
    ) |>
    mutate(
      cognition_priority_score =
        .data$n_fdr_better_cognition * 3 +
        .data$n_nominal_better_cognition +
        if_else(
          !is.na(.data$mmse_beta) & .data$mmse_beta > 0,
          0.5,
          0
        )
    ) |>
    arrange(
      desc(.data$cognition_priority_score),
      .data$best_q,
      desc(.data$mean_beta)
    )
}

############################################################
## 3. Explicitly confirm the structural concern
############################################################

structural_audit <- tibble(
  item = c(
    "prot_mat_adj nuisance covariates",
    "old cognition regression covariates",
    "overlap causing repeated adjustment",
    "corrected abundance matrix",
    "corrected direct cognition covariates"
  ),
  value = c(
    "age + sex + PMI + TMT batch",
    "age + sex + education + PMI + Braak + CERAD",
    "age + sex + PMI",
    "prot_mat_raw (processed, non-residualized)",
    "age + sex + education + PMI + Braak + CERAD + TMT batch"
  )
)

readr::write_csv(
  structural_audit,
  file.path(out_dir, "01_structural_double_adjustment_audit.csv")
)

############################################################
## 4. Network-level score and models
############################################################

hsp_genes <- all_hsp60_10_client_tbl |>
  filter(.data$detected_in_protein) |>
  transmute(gene = canonical_gene_symbol(.data$gene)) |>
  distinct() |>
  pull(.data$gene)

hsp_genes_old <- intersect(hsp_genes, rownames(prot_mat_adj))
hsp_genes_new <- intersect(hsp_genes, rownames(prot_mat_raw))

if (length(hsp_genes_old) != 306 || length(hsp_genes_new) != 306) {
  stop(
    "Expected 306 detected Hsp60/10 proteins in both matrices; observed old=",
    length(hsp_genes_old),
    ", corrected=",
    length(hsp_genes_new),
    ".",
    call. = FALSE
  )
}

old_network_score <- colMeans(
  row_z(prot_mat_adj[hsp_genes_old, , drop = FALSE]),
  na.rm = TRUE
)

new_network_score <- colMeans(
  row_z(prot_mat_raw[hsp_genes_new, , drop = FALSE]),
  na.rm = TRUE
)

network_df <- cog |>
  select(
    .data$SampleID,
    .data$cogdx_better_z,
    .data$dcfdx_lv_better_z,
    .data$mmse_last_valid_z,
    .data$cog_group,
    .data$age_death,
    .data$sex,
    .data$educ,
    .data$pmi,
    .data$Braak,
    .data$CERAD,
    .data$batch_factor
  ) |>
  mutate(
    old_score_z = z_or_na(old_network_score[
      match(.data$SampleID, names(old_network_score))
    ]),
    corrected_score_z = z_or_na(new_network_score[
      match(.data$SampleID, names(new_network_score))
    ])
  )

network_outcomes <- tribble(
  ~outcome, ~label,
  "cogdx_better_z", "Final cognitive diagnosis",
  "dcfdx_lv_better_z", "Last-valid cognitive diagnosis",
  "mmse_last_valid_z", "Last-valid MMSE"
)

network_comparison <- purrr::map_dfr(
  seq_len(nrow(network_outcomes)),
  function(i) {
    out <- network_outcomes$outcome[[i]]
    lab <- network_outcomes$label[[i]]

    old <- run_one_lm(
      network_df,
      out,
      "old_score_z",
      corrected = FALSE
    ) |>
      mutate(
        outcome = lab,
        specification = "old_adjusted_matrix_plus_covariates"
      )

    new <- run_one_lm(
      network_df,
      out,
      "corrected_score_z",
      corrected = TRUE
    ) |>
      mutate(
        outcome = lab,
        specification = "corrected_raw_matrix_plus_covariates_and_batch"
      )

    bind_rows(old, new)
  }
)

network_wide <- network_comparison |>
  select(
    .data$outcome,
    .data$specification,
    .data$n,
    .data$estimate,
    .data$std.error,
    .data$p.value,
    .data$conf.low,
    .data$conf.high
  ) |>
  pivot_wider(
    names_from = .data$specification,
    values_from = c(
      .data$n,
      .data$estimate,
      .data$std.error,
      .data$p.value,
      .data$conf.low,
      .data$conf.high
    )
  )

run_group_trend <- function(score_col, corrected) {
  covars <- c("age_death", "sex", "educ", "pmi", "Braak", "CERAD")
  if (corrected) covars <- c(covars, "batch_factor")

  dat <- network_df |>
    mutate(
      cog_group_num = case_when(
        .data$cog_group == "NCI" ~ 0,
        .data$cog_group == "MCI" ~ 1,
        .data$cog_group == "AD" ~ 2,
        TRUE ~ NA_real_
      )
    ) |>
    select(
      all_of(c(score_col, "cog_group_num", covars))
    ) |>
    drop_na()

  form <- reformulate(
    c("cog_group_num", covars),
    response = score_col
  )

  fit <- lm(form, data = dat)

  broom::tidy(fit, conf.int = TRUE) |>
    filter(.data$term == "cog_group_num") |>
    transmute(
      n = nobs(fit),
      estimate = .data$estimate,
      std.error = .data$std.error,
      statistic = .data$statistic,
      p.value = .data$p.value,
      conf.low = .data$conf.low,
      conf.high = .data$conf.high
    )
}

group_trend_comparison <- bind_rows(
  run_group_trend("old_score_z", FALSE) |>
    mutate(specification = "old_adjusted_matrix_plus_covariates"),
  run_group_trend("corrected_score_z", TRUE) |>
    mutate(specification = "corrected_raw_matrix_plus_covariates_and_batch")
)

readr::write_csv(
  network_comparison,
  file.path(out_dir, "02_network_level_old_vs_corrected_models.csv")
)

readr::write_csv(
  network_wide,
  file.path(out_dir, "03_network_level_old_vs_corrected_wide.csv")
)

readr::write_csv(
  group_trend_comparison,
  file.path(out_dir, "04_diagnosis_group_trend_old_vs_corrected.csv")
)

############################################################
## 5. Client-level Figure 4 cognition audit
############################################################

client_old <- run_gene_models(
  prot_mat_adj,
  hsp_genes_old,
  spec_label = "old_adjusted_matrix_plus_covariates",
  corrected = FALSE
)

client_new <- run_gene_models(
  prot_mat_raw,
  hsp_genes_new,
  spec_label = "corrected_raw_matrix_plus_covariates_and_batch",
  corrected = TRUE
)

client_pair <- client_old |>
  select(
    gene,
    outcome,
    old_n = n,
    old_estimate = estimate,
    old_p = p.value,
    old_q = q.value,
    old_fdr_positive = fdr_positive
  ) |>
  inner_join(
    client_new |>
      select(
        gene,
        outcome,
        corrected_n = n,
        corrected_estimate = estimate,
        corrected_p = p.value,
        corrected_q = q.value,
        corrected_fdr_positive = fdr_positive
      ),
    by = c("gene", "outcome")
  )

client_outcome_summary <- client_pair |>
  group_by(.data$outcome) |>
  summarise(
    n_genes = n(),
    spearman_beta_old_vs_corrected =
      suppressWarnings(
        cor(
          .data$old_estimate,
          .data$corrected_estimate,
          method = "spearman",
          use = "complete.obs"
        )
      ),
    sign_agreement =
      mean(
        sign(.data$old_estimate) ==
          sign(.data$corrected_estimate),
        na.rm = TRUE
      ),
    old_fdr_positive =
      sum(.data$old_fdr_positive, na.rm = TRUE),
    corrected_fdr_positive =
      sum(.data$corrected_fdr_positive, na.rm = TRUE),
    old_and_corrected_fdr_positive =
      sum(
        .data$old_fdr_positive &
          .data$corrected_fdr_positive,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

client_old_priority <- summarize_gene_priority(client_old) |>
  mutate(old_rank = row_number())

client_new_priority <- summarize_gene_priority(client_new) |>
  mutate(corrected_rank = row_number())

client_priority_comparison <- client_old_priority |>
  select(
    gene,
    old_rank,
    old_cognition_priority_score = cognition_priority_score,
    old_n_fdr = n_fdr_better_cognition,
    old_n_nominal = n_nominal_better_cognition,
    old_best_q = best_q
  ) |>
  full_join(
    client_new_priority |>
      select(
        gene,
        corrected_rank,
        corrected_cognition_priority_score =
          cognition_priority_score,
        corrected_n_fdr = n_fdr_better_cognition,
        corrected_n_nominal = n_nominal_better_cognition,
        corrected_best_q = best_q
      ),
    by = "gene"
  )

top25_old <- head(client_old_priority$gene, 25)
top25_new <- head(client_new_priority$gene, 25)

client_priority_summary <- tibble(
  metric = c(
    "Spearman priority score",
    "Top-25 overlap",
    "Top-25 Jaccard",
    "Genes FDR-positive for all 3 outcomes: old",
    "Genes FDR-positive for all 3 outcomes: corrected"
  ),
  value = c(
    suppressWarnings(
      cor(
        client_priority_comparison$old_cognition_priority_score,
        client_priority_comparison$corrected_cognition_priority_score,
        method = "spearman",
        use = "complete.obs"
      )
    ),
    length(intersect(top25_old, top25_new)),
    length(intersect(top25_old, top25_new)) /
      length(union(top25_old, top25_new)),
    sum(client_old_priority$n_fdr_better_cognition == 3),
    sum(client_new_priority$n_fdr_better_cognition == 3)
  )
)

readr::write_csv(
  bind_rows(client_old, client_new),
  file.path(out_dir, "05_client_level_all_models_old_vs_corrected.csv")
)
readr::write_csv(
  client_pair,
  file.path(out_dir, "06_client_level_gene_outcome_comparison.csv")
)
readr::write_csv(
  client_outcome_summary,
  file.path(out_dir, "07_client_level_outcome_summary.csv")
)
readr::write_csv(
  client_priority_comparison,
  file.path(out_dir, "08_client_priority_old_vs_corrected.csv")
)
readr::write_csv(
  client_priority_summary,
  file.path(out_dir, "09_client_priority_summary.csv")
)

############################################################
## 6. Figure 5 915-gene cognition audit
############################################################

fig5_genes <- canonical_gene_symbol(null_input_tbl$gene)

if (length(fig5_genes) != 915 ||
    anyDuplicated(fig5_genes) > 0) {
  stop(
    "Expected 915 unique Figure 5 genes after canonicalization.",
    call. = FALSE
  )
}

if (!all(fig5_genes %in% rownames(prot_mat_adj)) ||
    !all(fig5_genes %in% rownames(prot_mat_raw))) {
  stop(
    "Not all 915 Figure 5 genes are present in both protein matrices.",
    call. = FALSE
  )
}

fig5_old <- run_gene_models(
  prot_mat_adj,
  fig5_genes,
  spec_label = "old_adjusted_matrix_plus_covariates",
  corrected = FALSE
)

fig5_new <- run_gene_models(
  prot_mat_raw,
  fig5_genes,
  spec_label = "corrected_raw_matrix_plus_covariates_and_batch",
  corrected = TRUE
)

fig5_old_priority <- summarize_gene_priority(fig5_old)
fig5_new_priority <- summarize_gene_priority(fig5_new)

fig5_priority_pair <- fig5_old_priority |>
  select(
    gene,
    old_score = cognition_priority_score,
    old_n_fdr = n_fdr_better_cognition,
    old_n_nominal = n_nominal_better_cognition
  ) |>
  inner_join(
    fig5_new_priority |>
      select(
        gene,
        corrected_score = cognition_priority_score,
        corrected_n_fdr = n_fdr_better_cognition,
        corrected_n_nominal = n_nominal_better_cognition
      ),
    by = "gene"
  )

null_for_cognition <- null_input_tbl |>
  mutate(gene = canonical_gene_symbol(.data$gene)) |>
  select(
    .data$gene,
    .data$is_hsp60_10_client,
    .data$abundance_bin
  ) |>
  left_join(
    fig5_priority_pair,
    by = "gene"
  )

if (any(!is.finite(null_for_cognition$old_score)) ||
    any(!is.finite(null_for_cognition$corrected_score))) {
  stop(
    "One or more Figure 5 genes lack old/corrected cognition priority score.",
    call. = FALSE
  )
}

sample_matched_rows <- function() {
  hsp <- null_for_cognition |>
    filter(.data$is_hsp60_10_client)
  bg <- null_for_cognition |>
    filter(!.data$is_hsp60_10_client)

  needed <- hsp |>
    count(.data$abundance_bin, name = "n_needed")

  purrr::map_dfr(
    seq_len(nrow(needed)),
    function(i) {
      bin_i <- needed$abundance_bin[[i]]
      n_i <- needed$n_needed[[i]]

      pool <- bg |>
        filter(.data$abundance_bin == bin_i)

      if (nrow(pool) == 0) {
        stop(
          "No background proteins available for abundance bin ",
          bin_i,
          ".",
          call. = FALSE
        )
      }

      pool |>
        slice_sample(
          n = n_i,
          replace = nrow(pool) < n_i
        )
    }
  )
}

set.seed(1300)
n_perm <- 10000L

message("Running cognition-only old/corrected matched nulls: ", n_perm)

old_null <- numeric(n_perm)
new_null <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  if (i == 1 || i %% 1000 == 0 || i == n_perm) {
    message("  cognition null iteration ", i, "/", n_perm)
  }

  ## Use the exact same abundance-matched sampled genes for the old and
  ## corrected cognition specifications within each iteration.
  sampled <- sample_matched_rows()
  old_null[[i]] <- mean(sampled$old_score, na.rm = TRUE)
  new_null[[i]] <- mean(sampled$corrected_score, na.rm = TRUE)
}

obs_old <- mean(
  null_for_cognition$old_score[
    null_for_cognition$is_hsp60_10_client
  ]
)

obs_new <- mean(
  null_for_cognition$corrected_score[
    null_for_cognition$is_hsp60_10_client
  ]
)

emp_p <- function(null, observed) {
  (sum(null >= observed) + 1) / (length(null) + 1)
}

fig5_null_comparison <- tibble(
  specification = c(
    "old_adjusted_matrix_plus_covariates",
    "corrected_raw_matrix_plus_covariates_and_batch"
  ),
  observed_hsp_mean = c(obs_old, obs_new),
  null_mean = c(mean(old_null), mean(new_null)),
  null_sd = c(sd(old_null), sd(new_null)),
  observed_to_null_ratio =
    c(obs_old / mean(old_null), obs_new / mean(new_null)),
  empirical_p_greater =
    c(emp_p(old_null, obs_old), emp_p(new_null, obs_new)),
  n_perm = n_perm
)

readr::write_csv(
  bind_rows(fig5_old, fig5_new),
  file.path(out_dir, "10_fig5_915gene_all_models_old_vs_corrected.csv")
)
readr::write_csv(
  fig5_priority_pair,
  file.path(out_dir, "11_fig5_915gene_priority_old_vs_corrected.csv")
)
readr::write_csv(
  fig5_null_comparison,
  file.path(out_dir, "12_fig5_cognition_null_old_vs_corrected.csv")
)

############################################################
## 7. Concise decision checkpoint
############################################################

network_decision <- network_comparison |>
  select(
    .data$outcome,
    .data$specification,
    .data$n,
    .data$estimate,
    .data$p.value
  )

decision_summary <- c(
  "COGNITION DOUBLE-ADJUSTMENT AUDIT",
  "",
  "Structural finding:",
  "  OLD: prot_mat_adj already removes age + sex + PMI + batch, then cognition models adjust age + sex + education + PMI + Braak + CERAD.",
  "  Therefore age, sex, and PMI are adjusted twice in the old cognition workflow.",
  "  CORRECTED: prot_mat_raw + age + sex + education + PMI + Braak + CERAD + batch in each cognition model.",
  "",
  "Network-level corrected models:",
  paste(
    capture.output(print(
      network_decision |>
        filter(
          .data$specification ==
            "corrected_raw_matrix_plus_covariates_and_batch"
        )
    )),
    collapse = "\n"
  ),
  "",
  "Client-level old-vs-corrected summary:",
  paste(capture.output(print(client_outcome_summary)), collapse = "\n"),
  "",
  "Client cognition-priority summary:",
  paste(capture.output(print(client_priority_summary)), collapse = "\n"),
  "",
  "Figure 5 cognition matched-null comparison:",
  paste(capture.output(print(fig5_null_comparison)), collapse = "\n")
)

writeLines(
  decision_summary,
  file.path(out_dir, "COGNITION_AUDIT_DECISION_SUMMARY.txt")
)

message("\n============================================================")
message("COGNITION DOUBLE-ADJUSTMENT AUDIT COMPLETE")
message("============================================================\n")

message("STRUCTURAL FINDING")
print(structural_audit, n = Inf)

message("\nNETWORK-LEVEL OLD VS CORRECTED")
print(network_comparison, n = Inf)

message("\nDIAGNOSIS-GROUP TREND OLD VS CORRECTED")
print(group_trend_comparison, n = Inf)

message("\nCLIENT-LEVEL OUTCOME SUMMARY")
print(client_outcome_summary, n = Inf)

message("\nCLIENT PRIORITY SUMMARY")
print(client_priority_summary, n = Inf)

message("\nFIGURE 5 COGNITION NULL OLD VS CORRECTED")
print(fig5_null_comparison, n = Inf)

message("\nAudit outputs written to: ", out_dir)
message(
  "No production scripts or production outputs were modified."
)
