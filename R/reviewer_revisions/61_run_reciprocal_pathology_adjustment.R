#!/usr/bin/env Rscript

############################################################
## 61_run_reciprocal_pathology_adjustment.R
##
## Reciprocal pathology-model robustness analysis.
##
## Primary models remain unchanged:
##   abundance ~ Braak + production covariates
##   abundance ~ CERAD + production covariates
##
## Robustness models add the other pathology measure:
##   abundance ~ Braak + CERAD + production covariates
##   abundance ~ CERAD + Braak + production covariates
##
## Uses processed, non-residualized TMT abundance (prot_mat_raw)
## and the same production protein covariates as the manuscript.
##
## Runs on the same 915-protein Figure 5 universe:
##   306 Hsp60/10 clients
##   609 non-client mitochondrial proteins
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(broom)
})

source("R/00_config.R")
source("R/01_utils.R")

load_required_object <- function(name) {
  path <- file.path(cfg$object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing production object: ", path, call. = FALSE)
  }
  readRDS(path)
}

prot_mat_raw <- load_required_object("prot_mat_raw")
prot_meta_adj <- load_required_object("prot_meta_adj")
protein_covars <- unique(as.character(load_required_object("protein_covars")))

require_columns(
  prot_meta_adj,
  c("SampleID", "braak_num_std", "cerad_num_std"),
  "protein production metadata"
)

missing_covars <- setdiff(protein_covars, colnames(prot_meta_adj))
if (length(missing_covars) > 0) {
  stop(
    "Production protein covariate(s) missing from prot_meta_adj: ",
    paste(missing_covars, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Recover the exact manuscript mitochondrial universe.
# ------------------------------------------------------------

null_input_path <- file.path(
  cfg$project_dir,
  "outputs",
  "main_figures",
  "tables",
  "figure5_null_input_COVARIATE_ADJUSTED.csv"
)

if (!file.exists(null_input_path)) {
  stop("Missing canonical Figure 5 null input: ", null_input_path, call. = FALSE)
}

universe <- readr::read_csv(
  null_input_path,
  show_col_types = FALSE
) |>
  transmute(
    gene = as.character(.data$gene),
    is_hsp60_10_client = as.logical(.data$is_hsp60_10_client)
  ) |>
  distinct(.data$gene, .keep_all = TRUE)

if (nrow(universe) != 915L) {
  stop("Expected 915 proteins in manuscript universe; observed ", nrow(universe), ".", call. = FALSE)
}
if (sum(universe$is_hsp60_10_client) != 306L) {
  stop("Expected 306 clients; observed ", sum(universe$is_hsp60_10_client), ".", call. = FALSE)
}
if (sum(!universe$is_hsp60_10_client) != 609L) {
  stop("Expected 609 non-client mitochondrial proteins.", call. = FALSE)
}

genes <- intersect(universe$gene, rownames(prot_mat_raw))
if (length(genes) != 915L) {
  stop(
    "All 915 manuscript-universe proteins must be present in prot_mat_raw; observed ",
    length(genes), ".",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Generic gene-level pathology model.
# ------------------------------------------------------------

fit_pathology_model <- function(gene, predictor, extra_pathology = NULL, min_n = 30L) {
  y <- safe_z(prot_mat_raw[gene, ])

  vars <- unique(c(predictor, extra_pathology, protein_covars))
  vars <- vars[!is.na(vars) & nzchar(vars)]

  dat <- tibble(
    abundance = y,
    SampleID = colnames(prot_mat_raw)
  ) |>
    left_join(
      prot_meta_adj |>
        select(any_of(c("SampleID", vars))),
      by = "SampleID"
    ) |>
    mutate(across(where(is.character), as.factor))

  keep_vars <- c("abundance", vars)
  keep <- complete.cases(dat[, keep_vars, drop = FALSE])

  if (sum(keep) < min_n || length(unique(dat[[predictor]][keep])) < 3L) {
    return(tibble(
      gene = gene,
      predictor = predictor,
      adjusted_for = ifelse(is.null(extra_pathology), "none", extra_pathology),
      beta = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p = NA_real_,
      n = sum(keep)
    ))
  }

  form <- as.formula(
    paste("abundance ~", paste(vars, collapse = " + "))
  )

  fit <- tryCatch(
    lm(form, data = dat[keep, , drop = FALSE]),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      gene = gene,
      predictor = predictor,
      adjusted_for = ifelse(is.null(extra_pathology), "none", extra_pathology),
      beta = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p = NA_real_,
      n = sum(keep)
    ))
  }

  hit <- broom::tidy(fit, conf.int = TRUE) |>
    filter(.data$term == predictor)

  if (nrow(hit) != 1L) {
    stop("Could not uniquely extract ", predictor, " for ", gene, call. = FALSE)
  }

  tibble(
    gene = gene,
    predictor = predictor,
    adjusted_for = ifelse(is.null(extra_pathology), "none", extra_pathology),
    beta = hit$estimate,
    se = hit$std.error,
    ci_low = hit$conf.low,
    ci_high = hit$conf.high,
    p = hit$p.value,
    n = sum(keep)
  )
}

message("Fitting primary and reciprocal pathology models for 915 proteins...")

results <- bind_rows(
  map_dfr(genes, fit_pathology_model, predictor = "braak_num_std"),
  map_dfr(genes, fit_pathology_model,
          predictor = "braak_num_std", extra_pathology = "cerad_num_std"),
  map_dfr(genes, fit_pathology_model, predictor = "cerad_num_std"),
  map_dfr(genes, fit_pathology_model,
          predictor = "cerad_num_std", extra_pathology = "braak_num_std")
) |>
  left_join(universe, by = "gene") |>
  group_by(.data$predictor, .data$adjusted_for, .data$is_hsp60_10_client) |>
  mutate(q_bh = p.adjust(.data$p, method = "BH")) |>
  ungroup() |>
  mutate(
    protein_set = if_else(
      .data$is_hsp60_10_client,
      "Hsp60/10 clients",
      "Non-client mitochondrial proteins"
    ),
    pathology = case_when(
      .data$predictor == "braak_num_std" ~ "Braak",
      .data$predictor == "cerad_num_std" ~ "CERAD",
      TRUE ~ NA_character_
    ),
    specification = case_when(
      .data$predictor == "braak_num_std" & .data$adjusted_for == "none" ~ "Braak primary",
      .data$predictor == "braak_num_std" & .data$adjusted_for == "cerad_num_std" ~ "Braak adjusted for CERAD",
      .data$predictor == "cerad_num_std" & .data$adjusted_for == "none" ~ "CERAD primary",
      .data$predictor == "cerad_num_std" & .data$adjusted_for == "braak_num_std" ~ "CERAD adjusted for Braak",
      TRUE ~ NA_character_
    ),
    # Pathology-aligned signed coefficient used by the manuscript:
    # higher Braak = worse pathology => aligned = -beta
    # lower CERAD = worse pathology => aligned = +beta
    aligned_beta = if_else(
      .data$pathology == "Braak",
      -.data$beta,
      .data$beta
    ),
    directional_magnitude = pmax(.data$aligned_beta, 0)
  )

# ------------------------------------------------------------
# Summary tables.
# ------------------------------------------------------------

set_summary <- results |>
  group_by(.data$specification, .data$pathology, .data$protein_set) |>
  summarise(
    n_models = sum(is.finite(.data$beta)),
    median_beta = median(.data$beta, na.rm = TRUE),
    median_aligned_beta = median(.data$aligned_beta, na.rm = TRUE),
    median_directional_magnitude = median(.data$directional_magnitude, na.rm = TRUE),
    n_fdr = sum(.data$q_bh < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

client_wide <- results |>
  filter(.data$is_hsp60_10_client) |>
  select(.data$gene, .data$specification, .data$beta, .data$aligned_beta, .data$q_bh) |>
  tidyr::pivot_wider(
    names_from = .data$specification,
    values_from = c(.data$beta, .data$aligned_beta, .data$q_bh),
    names_sep = "__"
  )

compare_one <- function(pathology) {
  if (pathology == "Braak") {
    p_spec <- "Braak primary"
    a_spec <- "Braak adjusted for CERAD"
  } else {
    p_spec <- "CERAD primary"
    a_spec <- "CERAD adjusted for Braak"
  }

  x <- results |>
    filter(
      .data$is_hsp60_10_client,
      .data$specification %in% c(p_spec, a_spec)
    ) |>
    select(.data$gene, .data$specification, .data$beta, .data$aligned_beta, .data$q_bh) |>
    tidyr::pivot_wider(
      names_from = .data$specification,
      values_from = c(.data$beta, .data$aligned_beta, .data$q_bh),
      names_sep = "__"
    )

  beta_primary <- x[[paste0("beta__", p_spec)]]
  beta_adjusted <- x[[paste0("beta__", a_spec)]]
  q_primary <- x[[paste0("q_bh__", p_spec)]]
  q_adjusted <- x[[paste0("q_bh__", a_spec)]]

  tibble(
    pathology = pathology,
    spearman_primary_vs_adjusted = suppressWarnings(
      cor(beta_primary, beta_adjusted, method = "spearman", use = "complete.obs")
    ),
    direction_agreement = mean(
      sign(beta_primary) == sign(beta_adjusted),
      na.rm = TRUE
    ),
    n_fdr_primary = sum(q_primary < 0.05, na.rm = TRUE),
    n_fdr_adjusted = sum(q_adjusted < 0.05, na.rm = TRUE)
  )
}

robustness_summary <- bind_rows(
  compare_one("Braak"),
  compare_one("CERAD")
)

# Client-vs-background comparisons for the directional robustness metric.
set_tests <- results |>
  group_by(.data$specification, .data$pathology) |>
  group_modify(~{
    dat <- .x
    a <- dat$directional_magnitude[dat$protein_set == "Hsp60/10 clients"]
    b <- dat$directional_magnitude[dat$protein_set == "Non-client mitochondrial proteins"]
    wt <- suppressWarnings(wilcox.test(a, b, exact = FALSE))
    tibble(
      client_median = median(a, na.rm = TRUE),
      background_median = median(b, na.rm = TRUE),
      wilcoxon_p = wt$p.value
    )
  }) |>
  ungroup()

out_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_revisions",
  "reciprocal_pathology_adjustment"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  results,
  file.path(out_dir, "01_all_915_primary_and_reciprocal_pathology_models.csv"),
  na = ""
)
readr::write_csv(
  set_summary,
  file.path(out_dir, "02_set_level_pathology_summary.csv"),
  na = ""
)
readr::write_csv(
  robustness_summary,
  file.path(out_dir, "03_client_primary_vs_adjusted_robustness_summary.csv"),
  na = ""
)
readr::write_csv(
  set_tests,
  file.path(out_dir, "04_client_vs_nonclient_directional_magnitude_tests.csv"),
  na = ""
)

cat("\n============================================================\n")
cat("RECIPROCAL PATHOLOGY ROBUSTNESS COMPLETE\n")
cat("============================================================\n\n")
print(robustness_summary, n = Inf)
cat("\nSet-level summaries:\n")
print(set_summary, n = Inf)
cat("\nClient-vs-background tests:\n")
print(set_tests, n = Inf)
cat("\nOutputs written to:\n", out_dir, "\n")
