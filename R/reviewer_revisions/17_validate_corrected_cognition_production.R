#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(tidyverse))

source("R/00_config.R")
source("R/01_utils.R")

read_script <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
r06_txt <- read_script("R/06_build_cognition_objects.R")
r13_txt <- read_script("R/13_make_main_figure_4_cognition.R")
r14_txt <- read_script("R/14_make_main_figure_5_matched_null_specificity.R")

static_checks <- tibble::tribble(
  ~check, ~passed,
  "R/06 uses prot_mat_raw", grepl("protein_mat <- prot_mat_raw", r06_txt, fixed = TRUE),
  "R/06 joins TMT batch", grepl("client score to TMT batch", r06_txt, fixed = TRUE),
  "R/13 client models use prot_mat_raw", grepl("protein_mat <- prot_mat_raw", r13_txt, fixed = TRUE),
  "R/13 network models include batch",
  grepl("Hsp60_client_score_z + age_death + sex + educ + pmi + Braak + CERAD + batch_factor", r13_txt, fixed = TRUE),
  "R/13 client models include batch",
  grepl("client_abundance_z + age_death + sex + educ + pmi + Braak + CERAD + batch_factor", r13_txt, fixed = TRUE),
  "R/14 uses prot_mat_raw", grepl('get("prot_mat_raw"', r14_txt, fixed = TRUE),
  "R/14 cognition models include batch",
  grepl("abundance_z + age_death + sex + educ + pmi + Braak + CERAD + batch_factor", r14_txt, fixed = TRUE)
)

if (any(!static_checks$passed)) {
  print(static_checks, n = Inf)
  stop("Static corrected-cognition implementation checks failed.", call. = FALSE)
}

load_obj <- function(name, assign_as = name) {
  path <- file.path(cfg$object_dir, paste0(name, ".rds"))
  if (!file.exists(path)) stop("Missing saved production object: ", path, call. = FALSE)
  assign(assign_as, readRDS(path), envir = .GlobalEnv)
}

for (nm in c(
  "prot_mat_raw", "prot_mat_adj", "prot_meta_adj", "prot_meta_path",
  "analysis_meta", "prot_scores", "all_hsp60_10_client_tbl", "priority_tbl"
)) load_obj(nm)

prot_mat <- prot_mat_adj

source("R/06_build_cognition_objects.R", local = .GlobalEnv)
source("R/13_make_main_figure_4_cognition.R", local = .GlobalEnv)

audit_dir <- file.path(cfg$project_dir, "outputs", "reviewer_revisions", "cognition_double_adjustment_audit")
audit_network <- readr::read_csv(file.path(audit_dir, "02_network_level_old_vs_corrected_models.csv"), show_col_types = FALSE) |>
  filter(.data$specification == "corrected_raw_matrix_plus_covariates_and_batch") |>
  arrange(.data$outcome)
audit_priority <- readr::read_csv(file.path(audit_dir, "08_client_priority_old_vs_corrected.csv"), show_col_types = FALSE)
audit_fig5 <- readr::read_csv(file.path(audit_dir, "12_fig5_cognition_null_old_vs_corrected.csv"), show_col_types = FALSE) |>
  filter(.data$specification == "corrected_raw_matrix_plus_covariates_and_batch")

prod_network <- network_cognition_tbl |> arrange(.data$outcome)

if (!identical(as.character(prod_network$outcome), as.character(audit_network$outcome))) {
  stop("Production/audit network cognition outcomes do not align.", call. = FALSE)
}

max_network_beta_diff <- max(abs(prod_network$estimate - audit_network$estimate), na.rm = TRUE)
max_network_p_diff <- max(abs(prod_network$p.value - audit_network$p.value), na.rm = TRUE)

prod_priority <- client_cognition_summary |>
  transmute(
    gene = as.character(.data$gene),
    production_score = as.numeric(.data$cognition_priority_score),
    production_n_fdr = as.integer(.data$n_fdr_better_cognition)
  )

priority_compare <- audit_priority |> inner_join(prod_priority, by = "gene")
if (nrow(priority_compare) != 306) stop("Expected 306 client genes in priority comparison.", call. = FALSE)

max_priority_score_diff <- max(
  abs(priority_compare$corrected_cognition_priority_score - priority_compare$production_score),
  na.rm = TRUE
)
n_fdr_mismatch <- sum(priority_compare$corrected_n_fdr != priority_compare$production_n_fdr, na.rm = TRUE)

load_obj("figure5_null_results", "null_results")
load_obj("figure5_null_summary", "null_summary")
load_obj("figure5_observed_stats", "observed_stats")

null_input_path <- file.path(cfg$table_dir, "figure5_null_input_COVARIATE_ADJUSTED.csv")
null_input_tbl <- readr::read_csv(null_input_path, show_col_types = FALSE)
hsp_null_tbl <- null_input_tbl |> filter(.data$is_hsp60_10_client)
background_null_pool <- null_input_tbl |> filter(!.data$is_hsp60_10_client)

output_dir <- file.path(cfg$plot_dir, "main_fig5_matched_null_specificity_exact")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("R/14_make_main_figure_5_matched_null_specificity.R", local = .GlobalEnv)

prod_cog_row <- null_summary |> filter(.data$metric == "cognition_priority_score")
if (nrow(prod_cog_row) != 1) stop("Expected exactly one production cognition null-summary row.", call. = FALSE)

observed_cognition_diff <- abs(prod_cog_row$observed - audit_fig5$observed_hsp_mean)
null_mean_diff <- abs(prod_cog_row$null_mean - audit_fig5$null_mean)

checks <- bind_rows(
  static_checks |>
    transmute(check = .data$check, observed = as.character(.data$passed), expected = "TRUE", passed = .data$passed),
  tibble(
    check = c(
      "cognition_model_df participants",
      "cognition_model_df batch complete",
      "Figure 4 network beta reproduces corrected audit",
      "Figure 4 network p reproduces corrected audit",
      "Figure 4 client priority score reproduces corrected audit",
      "Figure 4 client FDR counts reproduce corrected audit",
      "Figure 5 observed cognition mean reproduces corrected audit",
      "Figure 5 cognition null mean consistent with corrected audit",
      "Figure 5 cognition empirical p respects +1 floor"
    ),
    observed = c(
      nrow(cognition_model_df),
      sum(!is.na(cognition_model_df$batch_factor)),
      signif(max_network_beta_diff, 5),
      signif(max_network_p_diff, 5),
      signif(max_priority_score_diff, 5),
      n_fdr_mismatch,
      signif(observed_cognition_diff, 5),
      signif(null_mean_diff, 5),
      signif(prod_cog_row$empirical_p_greater, 7)
    ),
    expected = c("400","400","<1e-10","<1e-10","<1e-10","0","<1e-10","<0.03",">=9.999e-5"),
    passed = c(
      nrow(cognition_model_df) == 400,
      sum(!is.na(cognition_model_df$batch_factor)) == 400,
      max_network_beta_diff < 1e-10,
      max_network_p_diff < 1e-10,
      max_priority_score_diff < 1e-10,
      n_fdr_mismatch == 0,
      observed_cognition_diff < 1e-10,
      null_mean_diff < 0.03,
      prod_cog_row$empirical_p_greater >= 1/10001
    )
  )
)

validation_dir <- file.path(cfg$project_dir, "outputs", "reviewer_revisions", "corrected_cognition_production_validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(checks, file.path(validation_dir, "corrected_cognition_production_checks.csv"))
readr::write_csv(network_cognition_tbl, file.path(validation_dir, "production_network_cognition_models.csv"))
readr::write_csv(client_cognition_summary, file.path(validation_dir, "production_client_cognition_summary.csv"))
readr::write_csv(null_summary, file.path(validation_dir, "production_figure5_null_summary_with_corrected_cognition.csv"))

message("\n============================================================")
message("CORRECTED COGNITION PRODUCTION VALIDATION")
message("============================================================\n")
print(checks, n = Inf)
message("\nCORRECTED NETWORK COGNITION MODELS")
print(network_cognition_tbl, n = Inf)
message("\nCORRECTED FIGURE 5 COGNITION NULL")
print(prod_cog_row, n = Inf)

if (any(!checks$passed)) stop("Corrected cognition production validation FAILED.", call. = FALSE)

message("\nALL CORRECTED COGNITION PRODUCTION CHECKS PASSED.")
