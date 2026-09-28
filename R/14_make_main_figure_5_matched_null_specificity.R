############################################################
## 14_main_fig5_matched_null_specificity.R
## Literal migration of the final old plotting block.
## Inputs are supplied by adjusted upstream scripts.
############################################################

## Hsp60/10 clients define a selectively AD- and cognition-relevant
## mitochondrial subnetwork
##
## Requires:
##   null_results, null_summary, observed_stats,
##   hsp_null_tbl, background_null_pool, output_dir
############################################################

library(dplyr)
library(tidyr)
library(readr)
library(broom)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)

############################################################
## 0. Checks
############################################################

req_objs <- c(
  "null_results", "null_summary", "observed_stats",
  "hsp_null_tbl", "background_null_pool", "output_dir"
)

missing_objs <- req_objs[!vapply(req_objs, exists, logical(1))]
if (length(missing_objs) > 0) {
  stop(
    "Figure 5 missing required objects: ", paste(missing_objs, collapse = ", "),
    ". Available objects: ", available_objects_msg(),
    call. = FALSE
  )
}



############################################################
## 0b. Add cognition metric to matched-null objects if absent
##
## The original matched-null objects were built before the
## client-level cognition screen existed. The client-level cognition
## summary from Figure 4 only covers Hsp60/10 clients, so it cannot
## be used directly for non-client mitochondrial null sets.
##
## This block therefore computes the same covariate-adjusted
## cognition score for the full Figure 5 gene universe:
##   Hsp60/10 clients + non-client mitochondrial background proteins.
##
## It then joins cognition_priority_score onto hsp_null_tbl and
## background_null_pool and builds a matched null cognition distribution.
############################################################

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

require_gene_column <- function(tbl, label) {
  candidates <- c("gene", "Gene", "GENE", "symbol", "Symbol", "Gene Symbol", "gene_symbol", "hgnc_symbol")
  hit <- candidates[candidates %in% colnames(tbl)][1]
  if (length(hit) == 0 || is.na(hit)) {
    stop(
      "Could not find gene column in ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(tbl), collapse = ", "),
      call. = FALSE
    )
  }
  hit
}

## Reconstruct cognition outcomes from cognition_model_df.
prepare_cognition_df_for_fig5 <- function(cognition_model_df) {
  required <- c(
    "SampleID", "individualID", "cogdx", "dcfdx_lv", "mmse_last_valid",
    "age_death", "sex", "educ", "pmi", "Braak", "CERAD", "batch_factor"
  )

  missing <- setdiff(required, colnames(cognition_model_df))
  if (length(missing) > 0) {
    stop(
      "cognition_model_df is missing required columns for Figure 5 cognition scoring: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cognition_model_df %>%
    mutate(
      SampleID = as.character(SampleID),
      individualID = as.character(individualID),

      cogdx_severity = case_when(
        cogdx == 1 ~ 0,
        cogdx %in% c(2, 3) ~ 1,
        cogdx %in% c(4, 5) ~ 2,
        TRUE ~ NA_real_
      ),

      dcfdx_lv_severity = case_when(
        dcfdx_lv == 1 ~ 0,
        dcfdx_lv %in% c(2, 3) ~ 1,
        dcfdx_lv %in% c(4, 5) ~ 2,
        TRUE ~ NA_real_
      ),

      ## Reoriented so higher = better cognition.
      cogdx_better = case_when(
        cogdx_severity == 0 ~ 2,
        cogdx_severity == 1 ~ 1,
        cogdx_severity == 2 ~ 0,
        TRUE ~ NA_real_
      ),

      dcfdx_lv_better = case_when(
        dcfdx_lv_severity == 0 ~ 2,
        dcfdx_lv_severity == 1 ~ 1,
        dcfdx_lv_severity == 2 ~ 0,
        TRUE ~ NA_real_
      ),

      cogdx_better_z = as.numeric(scale(cogdx_better)),
      dcfdx_lv_better_z = as.numeric(scale(dcfdx_lv_better)),
      mmse_last_valid_z = as.numeric(scale(mmse_last_valid)),
      sex = as.factor(sex)
    ) %>%
    select(
      SampleID,
      individualID,
      cogdx_better_z,
      dcfdx_lv_better_z,
      mmse_last_valid_z,
      age_death,
      sex,
      educ,
      pmi,
      Braak,
      CERAD,
      batch_factor
    )
}

run_gene_cognition_model_fig5 <- function(data, outcome, outcome_label) {
  model_df <- data %>%
    select(
      gene,
      abundance_z,
      all_of(outcome),
      age_death,
      sex,
      educ,
      pmi,
      Braak,
      CERAD,
      batch_factor
    ) %>%
    drop_na()

  message(
    "Running Figure 5 cognition scoring for ", outcome_label,
    ": n rows = ", nrow(model_df),
    "; n genes = ", dplyr::n_distinct(model_df$gene)
  )

  model_df %>%
    group_by(gene) %>%
    group_modify(~ {
      dat <- .x

      if (
        nrow(dat) < 40 ||
        sd(dat$abundance_z, na.rm = TRUE) == 0 ||
        sd(dat[[outcome]], na.rm = TRUE) == 0
      ) {
        return(
          tibble(
            outcome = outcome_label,
            n = nrow(dat),
            estimate = NA_real_,
            std.error = NA_real_,
            statistic = NA_real_,
            p.value = NA_real_,
            conf.low = NA_real_,
            conf.high = NA_real_,
            model_status = "skipped_low_n_or_no_variance"
          )
        )
      }

      form <- as.formula(
        paste0(
          outcome,
          " ~ abundance_z + age_death + sex + educ + pmi + Braak + CERAD + batch_factor"
        )
      )

      fit <- tryCatch(
        lm(form, data = dat),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        return(
          tibble(
            outcome = outcome_label,
            n = nrow(dat),
            estimate = NA_real_,
            std.error = NA_real_,
            statistic = NA_real_,
            p.value = NA_real_,
            conf.low = NA_real_,
            conf.high = NA_real_,
            model_status = paste0("model_error: ", fit$message)
          )
        )
      }

      broom::tidy(fit, conf.int = TRUE) %>%
        filter(term == "abundance_z") %>%
        transmute(
          outcome = outcome_label,
          n = nobs(fit),
          estimate,
          std.error,
          statistic,
          p.value,
          conf.low,
          conf.high,
          model_status = "ok"
        )
    }) %>%
    ungroup()
}

compute_gene_cognition_scores_for_fig5 <- function(gene_universe) {
  if (!exists("prot_mat_raw", envir = .GlobalEnv)) {
    stop(
      "prot_mat_raw is required to compute cognition scores for background proteins. ",
      "Run upstream proteomics setup before Figure 5.",
      call. = FALSE
    )
  }

  if (!exists("cognition_model_df", envir = .GlobalEnv)) {
    stop(
      "cognition_model_df is required to compute cognition scores for background proteins. ",
      "Run 06_build_cognition_objects.R before Figure 5.",
      call. = FALSE
    )
  }

  protein_mat <- get("prot_mat_raw", envir = .GlobalEnv)
  rownames(protein_mat) <- toupper(as.character(rownames(protein_mat)))
  colnames(protein_mat) <- as.character(colnames(protein_mat))

  gene_universe <- toupper(as.character(gene_universe))
  gene_universe <- unique(gene_universe[!is.na(gene_universe) & gene_universe != ""])

  detected_gene_universe <- intersect(gene_universe, rownames(protein_mat))
  missing_from_prot_mat <- setdiff(gene_universe, detected_gene_universe)

  message("Figure 5 cognition gene universe: ", length(gene_universe))
  message("Detected in prot_mat_raw: ", length(detected_gene_universe))
  message("Missing from prot_mat_raw: ", length(missing_from_prot_mat))

  if (length(detected_gene_universe) < 20) {
    stop(
      "Too few Figure 5 gene-universe proteins are detected in prot_mat_raw to compute cognition scores.",
      call. = FALSE
    )
  }

  gene_mat <- protein_mat[detected_gene_universe, , drop = FALSE]
  gene_mat <- apply(gene_mat, 2, as.numeric)
  rownames(gene_mat) <- detected_gene_universe

  gene_mat_z <- t(scale(t(gene_mat)))
  gene_mat_z <- gene_mat_z[
    rowSums(!is.na(gene_mat_z)) > 0,
    ,
    drop = FALSE
  ]

  abundance_long <- as.data.frame(gene_mat_z) %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "SampleID",
      values_to = "abundance_z"
    ) %>%
    mutate(
      gene = toupper(as.character(gene)),
      SampleID = as.character(SampleID),
      abundance_z = as.numeric(abundance_z)
    )

  cognition_df <- prepare_cognition_df_for_fig5(get("cognition_model_df", envir = .GlobalEnv))

  analysis_long <- abundance_long %>%
    left_join(cognition_df, by = "SampleID")

  gene_cognition_all <- bind_rows(
    run_gene_cognition_model_fig5(
      analysis_long,
      outcome = "mmse_last_valid_z",
      outcome_label = "Last-valid MMSE"
    ),
    run_gene_cognition_model_fig5(
      analysis_long,
      outcome = "dcfdx_lv_better_z",
      outcome_label = "Last-valid cognitive diagnosis"
    ),
    run_gene_cognition_model_fig5(
      analysis_long,
      outcome = "cogdx_better_z",
      outcome_label = "Final cognitive diagnosis"
    )
  ) %>%
    group_by(outcome) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      direction = case_when(
        estimate > 0 ~ "higher abundance = better cognition",
        estimate < 0 ~ "higher abundance = worse cognition",
        TRUE ~ NA_character_
      ),
      cognition_associated_nominal = model_status == "ok" & p.value < 0.05 & estimate > 0,
      cognition_associated_fdr = model_status == "ok" & q.value < 0.05 & estimate > 0
    )

  gene_cognition_summary <- gene_cognition_all %>%
    filter(model_status == "ok") %>%
    group_by(gene) %>%
    summarise(
      n_outcomes_tested = n(),

      n_nominal_better_cognition = sum(cognition_associated_nominal, na.rm = TRUE),
      n_fdr_better_cognition = sum(cognition_associated_fdr, na.rm = TRUE),

      best_p = min(p.value, na.rm = TRUE),
      best_q = min(q.value, na.rm = TRUE),

      mean_beta = mean(estimate, na.rm = TRUE),
      median_beta = median(estimate, na.rm = TRUE),

      mmse_beta = estimate[outcome == "Last-valid MMSE"][1],
      mmse_p = p.value[outcome == "Last-valid MMSE"][1],
      mmse_q = q.value[outcome == "Last-valid MMSE"][1],

      last_dx_beta = estimate[outcome == "Last-valid cognitive diagnosis"][1],
      last_dx_p = p.value[outcome == "Last-valid cognitive diagnosis"][1],
      last_dx_q = q.value[outcome == "Last-valid cognitive diagnosis"][1],

      final_dx_beta = estimate[outcome == "Final cognitive diagnosis"][1],
      final_dx_p = p.value[outcome == "Final cognitive diagnosis"][1],
      final_dx_q = q.value[outcome == "Final cognitive diagnosis"][1],

      .groups = "drop"
    ) %>%
    mutate(
      cognition_priority_score =
        n_fdr_better_cognition * 3 +
        n_nominal_better_cognition +
        if_else(!is.na(mmse_beta) & mmse_beta > 0, 0.5, 0),
      cognition_priority_tier = case_when(
        n_fdr_better_cognition >= 2 ~ "FDR-supported across multiple cognition outcomes",
        n_fdr_better_cognition == 1 ~ "FDR-supported in one cognition outcome",
        n_nominal_better_cognition >= 2 ~ "Nominal across multiple cognition outcomes",
        n_nominal_better_cognition == 1 ~ "Nominal in one cognition outcome",
        TRUE ~ "Not cognition-prioritized"
      )
    ) %>%
    arrange(desc(cognition_priority_score), best_q, desc(mean_beta))

  if (exists("cfg", envir = .GlobalEnv)) {
    fig5_tab_dir <- file.path(cfg$table_dir, "main_fig5_matched_null_specificity")
    dir.create(fig5_tab_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(gene_cognition_all, file.path(fig5_tab_dir, "Fig5_gene_universe_cognition_all_models.csv"))
    readr::write_csv(gene_cognition_summary, file.path(fig5_tab_dir, "Fig5_gene_universe_cognition_summary_by_gene.csv"))
  }

  gene_cognition_summary %>%
    select(
      gene,
      cognition_priority_score,
      n_fdr_better_cognition,
      n_nominal_better_cognition,
      mean_beta,
      median_beta,
      mmse_beta,
      mmse_q,
      last_dx_beta,
      last_dx_q,
      final_dx_beta,
      final_dx_q
    )
}

add_cognition_to_gene_table <- function(tbl, cognition_tbl, label) {
  gene_col <- require_gene_column(tbl, label)

  tbl2 <- tbl %>%
    mutate(gene = toupper(as.character(.data[[gene_col]])))

  existing_cognition_cols <- intersect(
    c(
      "cognition_priority_score",
      "n_fdr_better_cognition",
      "n_nominal_better_cognition",
      "mean_beta",
      "median_beta",
      "mmse_beta",
      "mmse_q",
      "last_dx_beta",
      "last_dx_q",
      "final_dx_beta",
      "final_dx_q"
    ),
    colnames(tbl2)
  )

  if (length(existing_cognition_cols) > 0) {
    tbl2 <- tbl2 %>% select(-all_of(existing_cognition_cols))
  }

  tbl2 %>%
    left_join(cognition_tbl, by = "gene")
}

cognition_metric_candidates <- c(
  "cognition_priority_score",
  "cognition_support_score",
  "cognition_preservation_score",
  "cognitive_preservation_score",
  "mean_cognition_priority_score"
)

null_cognition_candidates <- c(
  "null_mean_cognition_score",
  "null_mean_cognition_priority_score",
  "null_mean_cognition_support_score",
  "null_mean_cognition_preservation_score",
  "null_mean_cognitive_preservation_score"
)

existing_cognition_metric <- intersect(cognition_metric_candidates, as.character(null_summary$metric))[1]
existing_null_cognition_col <- intersect(null_cognition_candidates, colnames(null_results))[1]

if (is.na(existing_cognition_metric) || is.na(existing_null_cognition_col)) {

  hsp_gene_col <- require_gene_column(hsp_null_tbl, "Hsp60/10 null table")
  bg_gene_col <- require_gene_column(background_null_pool, "matched mitochondrial background pool")

  fig5_gene_universe <- union(
    toupper(as.character(hsp_null_tbl[[hsp_gene_col]])),
    toupper(as.character(background_null_pool[[bg_gene_col]]))
  )

  cognition_tbl <- compute_gene_cognition_scores_for_fig5(fig5_gene_universe)

  hsp_null_tbl <- add_cognition_to_gene_table(
    hsp_null_tbl,
    cognition_tbl,
    "Hsp60/10 null table"
  )

  background_null_pool <- add_cognition_to_gene_table(
    background_null_pool,
    cognition_tbl,
    "matched mitochondrial background pool"
  )

  n_hsp_with_cognition <- sum(!is.na(hsp_null_tbl$cognition_priority_score))
  n_background_with_cognition <- sum(!is.na(background_null_pool$cognition_priority_score))

  message("Hsp60/10 clients with cognition score: ", n_hsp_with_cognition, " / ", nrow(hsp_null_tbl))
  message("Background proteins with cognition score: ", n_background_with_cognition, " / ", nrow(background_null_pool))

  if (n_hsp_with_cognition < 5) {
    stop(
      "Too few Hsp60/10 clients have cognition scores after computing Figure 5 gene-universe cognition models.",
      call. = FALSE
    )
  }

  if (n_background_with_cognition < n_hsp_with_cognition) {
    stop(
      "Too few background proteins have cognition scores to build cognition matched nulls. ",
      "Need at least as many scored background proteins as scored Hsp60/10 clients.",
      call. = FALSE
    )
  }

  observed_cognition <- mean(hsp_null_tbl$cognition_priority_score, na.rm = TRUE)

  ## Prefer reusing gene-level null samples if null_results contains them.
  iter_col <- intersect(c("iteration", "iter", "null_iter", "null_iteration"), colnames(null_results))[1]
  null_gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol", "sampled_gene"), colnames(null_results))[1]

  if (!is.na(iter_col) && !is.na(null_gene_col) && nrow(null_results) > length(unique(null_results[[iter_col]]))) {
    message("Computing cognition null from gene-level null_results using columns: ", iter_col, ", ", null_gene_col)

    null_cognition_tbl <- null_results %>%
      mutate(gene = toupper(as.character(.data[[null_gene_col]]))) %>%
      left_join(cognition_tbl %>% select(gene, cognition_priority_score), by = "gene") %>%
      group_by(.data[[iter_col]]) %>%
      summarise(
        null_mean_cognition_score = mean(cognition_priority_score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(iteration_tmp = all_of(iter_col))

    null_results <- null_results %>%
      left_join(null_cognition_tbl, by = setNames("iteration_tmp", iter_col))

  } else {
    message("null_results does not contain gene-level sampled null sets; resampling cognition nulls from background pool.")

    n_iter <- nrow(null_results)
    if (!is.finite(n_iter) || n_iter < 10) n_iter <- 10000

    bin_candidates <- c(
      "abundance_bin", "mean_abundance_bin", "abundance_quantile", "abundance_decile",
      "matching_bin", "match_bin", "bin", "abundance_stratum", "abundance_match_bin"
    )

    bin_col <- intersect(bin_candidates, intersect(colnames(hsp_null_tbl), colnames(background_null_pool)))[1]

    sample_one_null <- function() {
      if (!is.na(bin_col)) {
        target_counts <- hsp_null_tbl %>%
          filter(!is.na(.data[[bin_col]]), !is.na(cognition_priority_score)) %>%
          count(across(all_of(bin_col)), name = "n")

        sampled <- lapply(seq_len(nrow(target_counts)), function(ii) {
          bin_value <- target_counts[[bin_col]][ii]
          n_needed <- target_counts$n[ii]

          pool <- background_null_pool %>%
            filter(
              .data[[bin_col]] == bin_value,
              !is.na(cognition_priority_score)
            )

          if (nrow(pool) == 0) {
            stop("No scored background proteins available for abundance bin: ", bin_value, call. = FALSE)
          }

          pool[sample(seq_len(nrow(pool)), size = n_needed, replace = nrow(pool) < n_needed), , drop = FALSE]
        }) %>%
          bind_rows()

        mean(sampled$cognition_priority_score, na.rm = TRUE)

      } else {
        pool <- background_null_pool %>%
          filter(!is.na(cognition_priority_score))

        pool <- pool[sample(seq_len(nrow(pool)), size = n_hsp_with_cognition, replace = nrow(pool) < n_hsp_with_cognition), , drop = FALSE]
        mean(pool$cognition_priority_score, na.rm = TRUE)
      }
    }

    if (!is.na(bin_col)) {
      message("Cognition null resampling is stratified by: ", bin_col)
    } else {
      warning(
        "Could not find an abundance-bin column shared by hsp_null_tbl and background_null_pool. ",
        "Cognition null will be sampled without abundance-bin stratification.",
        call. = FALSE
      )
    }

    set.seed(1405)
    null_cognition_vec <- replicate(n_iter, sample_one_null())

    if (nrow(null_results) == length(null_cognition_vec)) {
      null_results$null_mean_cognition_score <- as.numeric(null_cognition_vec)
    } else {
      null_results <- tibble(null_mean_cognition_score = as.numeric(null_cognition_vec)) %>%
        bind_cols(null_results[seq_len(min(nrow(null_results), length(null_cognition_vec))), , drop = FALSE])
    }
  }

  if (!"null_mean_cognition_score" %in% colnames(null_results)) {
    stop("Failed to create null_mean_cognition_score in null_results.", call. = FALSE)
  }

  cognition_summary_row <- tibble(
    metric = "cognition_priority_score",
    observed = observed_cognition,
    null_mean = mean(null_results$null_mean_cognition_score, na.rm = TRUE),
    null_sd = sd(null_results$null_mean_cognition_score, na.rm = TRUE),
    empirical_p_greater = (
      sum(
        null_results$null_mean_cognition_score >= observed_cognition,
        na.rm = TRUE
      ) + 1
    ) / (
      sum(is.finite(null_results$null_mean_cognition_score)) + 1
    )
  )

  null_summary <- null_summary %>%
    filter(!metric %in% cognition_metric_candidates) %>%
    bind_rows(cognition_summary_row)

  message(
    "Added cognition_priority_score to null_summary: observed = ",
    signif(observed_cognition, 4),
    "; null mean = ", signif(cognition_summary_row$null_mean, 4),
    "; empirical p = ", signif(cognition_summary_row$empirical_p_greater, 4)
  )
}


require_columns(
  null_summary,
  c("metric", "observed", "null_mean", "null_sd", "empirical_p_greater"),
  "Figure 5 null summary"
)
require_columns(
  null_results,
  "null_agora_fraction",
  "Figure 5 null results"
)
## Figure 5 now uses cognition preservation instead of the older
## paired pathology-vulnerability composite. The upstream matched-null
## builder must provide one cognition metric in null_summary and the
## matching null mean column in null_results.
cognition_metric <- intersect(
  c(
    "cognition_priority_score",
    "cognition_support_score",
    "cognition_preservation_score",
    "cognitive_preservation_score",
    "mean_cognition_priority_score"
  ),
  as.character(null_summary$metric)
)[1]

if (is.na(cognition_metric)) {
  stop(
    "Figure 5 expected a cognition metric in null_summary but did not find one.\n",
    "Tried: cognition_priority_score, cognition_support_score, cognition_preservation_score, ",
    "cognitive_preservation_score, mean_cognition_priority_score.\n",
    "Available metrics: ", paste(unique(null_summary$metric), collapse = ", "),
    call. = FALSE
  )
}

null_cognition_col <- intersect(
  c(
    "null_mean_cognition_score",
    "null_mean_cognition_priority_score",
    "null_mean_cognition_support_score",
    "null_mean_cognition_preservation_score",
    "null_mean_cognitive_preservation_score"
  ),
  colnames(null_results)
)[1]

if (is.na(null_cognition_col)) {
  stop(
    "Figure 5 expected a cognition null-distribution column in null_results but did not find one.\n",
    "Tried: null_mean_cognition_score, null_mean_cognition_priority_score, ",
    "null_mean_cognition_support_score, null_mean_cognition_preservation_score, ",
    "null_mean_cognitive_preservation_score.\n",
    "Available columns: ", paste(colnames(null_results), collapse = ", "),
    call. = FALSE
  )
}

require_values(
  null_summary,
  "metric",
  c("protein_late_decline_magnitude", "inverse_braak_magnitude", cognition_metric, "agora_fraction"),
  "Figure 5 null summary"
)

output_dir <- file.path(cfg$plot_dir, "main_fig5_matched_null_specificity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Styling
############################################################

obs_col      <- "#0B3558"
obs_col2     <- "#174E7A"
null_fill    <- "gray72"
null_line    <- "gray55"
box_fill     <- "#EFF4F8"
box_fill2    <- "#DCEAF3"
dark_fill    <- "#0B3558"

theme_set(
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15.5, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(b = 5)),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10.8),
      plot.margin = margin(8, 8, 8, 8)
    )
)

############################################################
## 2. Summary table for plotting
############################################################

figure_tbl <- null_summary %>%
  filter(metric %in% c(
    "protein_late_decline_magnitude",
    "inverse_braak_magnitude",
    cognition_metric,
    "agora_fraction"
  )) %>%
  mutate(
    metric_clean = case_when(
      metric == "protein_late_decline_magnitude" ~ "Late-stage\nprotein decline",
      metric == "inverse_braak_magnitude" ~ "Inverse Braak\nassociation",
      metric == cognition_metric ~ "Cognition\npreservation score",
      metric == "agora_fraction" ~ "Agora target\nfraction",
      TRUE ~ as.character(metric)
    ),
    metric_clean = factor(
      metric_clean,
      levels = rev(c(
        "Late-stage\nprotein decline",
        "Inverse Braak\nassociation",
        "Cognition\npreservation score",
        "Agora target\nfraction"
      ))
    ),
    fold_enrichment = observed / null_mean,
    fold_label = case_when(
      metric == "agora_fraction" ~ paste0(round(fold_enrichment, 2), "x enriched"),
      TRUE ~ paste0(round(fold_enrichment, 2), "x higher")
    ),
    value_label = paste0(
      "Obs ", number(observed, accuracy = 0.001),
      " | Null ", number(null_mean, accuracy = 0.001)
    ),
    p_label = ifelse(
      empirical_p_greater <= 1 / nrow(null_results),
      paste0("empirical p <= ", format(1 / nrow(null_results), scientific = TRUE)),
      paste0("empirical p = ", signif(empirical_p_greater, 2))
    )
  )

print(figure_tbl)

############################################################
## 3. Panel A
############################################################

panel_A <- ggplot() +
  xlim(-1.8, 14.6) +
  ylim(0, 10.4) +

  annotate(
    "label",
    x = 1.2, y = 7.9,
    label = paste0("Hsp60/10\nclients\nn = ", nrow(hsp_null_tbl)),
    size = 2.70,
    fontface = "bold",
    lineheight = 0.90,
    label.padding = unit(0.36, "lines"),
    label.r = unit(0.18, "lines"),
    fill = box_fill,
    color = "black"
  ) +
  annotate(
    "label",
    x = 11.6, y = 7.9,
    label = paste0("Non-client\nmitochondrial\nproteins\nn = ", nrow(background_null_pool)),
    size = 2.52,
    fontface = "bold",
    lineheight = 0.90,
    label.padding = unit(0.36, "lines"),
    label.r = unit(0.18, "lines"),
    fill = box_fill,
    color = "black"
  ) +

  annotate(
    "segment",
    x = 3.6, xend = 6.8,
    y = 7.9, yend = 7.9,
    linewidth = 0.95,
    arrow = arrow(length = unit(0.16, "inches"))
  ) +
  annotate(
    "label",
    x = 6.15, y = 9.45,
    label = "abundance-matched",
    size = 3.0,
    fontface = "bold",
    label.padding = unit(0.18, "lines"),
    label.r = unit(0.12, "lines"),
    fill = "white",
    color = "black"
   ) +

  annotate(
    "label",
    x = 6.15, y = 5.55,
    label = "10,000 matched\nmitochondrial null sets",
    size = 3.35,
    fontface = "bold",
    label.padding = unit(0.38, "lines"),
    label.r = unit(0.18, "lines"),
    fill = box_fill2,
    color = obs_col
  ) +

  annotate(
    "segment",
    x = 6.15, xend = 6.15,
    y = 4.75, yend = 3.95,
    linewidth = 1.0,
    arrow = arrow(length = unit(0.16, "inches"))
  ) +

  annotate(
    "label",
    x = 6.15, y = 2.25,
    label = "Hsp60/10 clients exceed\nmatched-null expectation\nacross pathology, cognition,\nand AD target metrics",
    size = 2.70,
    fontface = "bold",
    label.padding = unit(0.34, "lines"),
    label.r = unit(0.18, "lines"),
    fill = dark_fill,
    color = "white"
  ) +

  labs(
    title = "A. Matched mitochondrial specificity test",
    subtitle = "Comparing Hsp60/10 clients against abundance-matched non-client mitochondrial background"
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(10, 18, 10, 18)
  )

############################################################
## 4. Panel B
############################################################

require_columns(
  figure_tbl,
  c("metric", "metric_clean", "fold_enrichment", "fold_label", "value_label", "p_label"),
  "Figure 5 Panel B summary table"
)

xmax_B <- max(figure_tbl$fold_enrichment, na.rm = TRUE) + 0.64

panel_B <- ggplot(figure_tbl, aes(y = metric_clean)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.8,
    color = "gray55"
  ) +
  geom_segment(
    aes(x = 1, xend = fold_enrichment, yend = metric_clean),
    linewidth = 4.4,
    color = "gray82",
    lineend = "round"
  ) +
  geom_point(
    aes(x = fold_enrichment),
    size = 3.9,
    color = obs_col
  ) +
  geom_text(
    aes(x = fold_enrichment, label = fold_label),
    nudge_x = 0.12,
    hjust = 0,
    size = 3.35,
    fontface = "bold"
  ) +
  geom_text(
    aes(x = 1.01, label = value_label),
    hjust = 0,
    nudge_y = -0.22,
    size = 2.55,
    color = "gray25"
  ) +
  scale_x_continuous(
    limits = c(0.98, xmax_B),
    breaks = c(1.0, 1.5, 2.0),
    labels = function(x) paste0(x, "x")
  ) +
  labs(
    title = "B. Enrichment over matched-null expectation",
    subtitle = "Observed Hsp60/10 mean divided by the matched mitochondrial null mean",
    x = "Observed / matched-null mean",
    y = NULL,
    caption = "Dashed line = no enrichment (1.0x)."
  ) +
  coord_cartesian(clip = "off") +
  theme(
    axis.text.y = element_text(size = 8.7, face = "bold", lineheight = 0.92),
    axis.text.x = element_text(size = 7.6),
    axis.title.x = element_text(size = 9.0),
    plot.caption = element_text(size = 8.5, hjust = 0),
    plot.margin = margin(8, 42, 8, 8)
  )

############################################################
## 5. Panel C
############################################################

obs_cognition <- figure_tbl %>%
  filter(metric == cognition_metric) %>%
  pull(observed)

null_cognition <- figure_tbl %>%
  filter(metric == cognition_metric) %>%
  pull(null_mean)

cognition_p <- figure_tbl %>%
  filter(metric == cognition_metric) %>%
  pull(p_label)
cognition_p_short <- sub("^empirical ", "", cognition_p)

require_columns(null_results, null_cognition_col, "Figure 5 Panel C null results")

cognition_x_min_data <- min(null_results[[null_cognition_col]], na.rm = TRUE)
cognition_x_max_data <- max(c(null_results[[null_cognition_col]], obs_cognition), na.rm = TRUE)
cognition_x_span <- cognition_x_max_data - cognition_x_min_data
xmin_cognition <- cognition_x_min_data - 0.03 * cognition_x_span
xmax_cognition <- max(cognition_x_max_data + 0.12 * cognition_x_span, obs_cognition + 0.42)
ymax_cognition <- max(ggplot_build(
  ggplot(null_results, aes(x = .data[[null_cognition_col]])) +
    geom_histogram(bins = 44)
)$data[[1]]$count)
cognition_label_x <- obs_cognition + 0.025 * (xmax_cognition - xmin_cognition)
cognition_label <- paste0(
  "Observed mean",
  "\nObs ", number(obs_cognition, accuracy = 0.001),
  "\nNull ", number(null_cognition, accuracy = 0.001),
  "\n", cognition_p_short
)

panel_C <- ggplot(null_results, aes(x = .data[[null_cognition_col]])) +
  geom_histogram(
    bins = 44,
    fill = null_fill,
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = obs_cognition,
    color = obs_col,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  annotate(
    "label",
    x = cognition_label_x,
    y = ymax_cognition * 0.90,
    label = cognition_label,
    hjust = 0,
    vjust = 1,
    size = 2.16,
    lineheight = 0.92,
    fill = "white",
    label.padding = unit(0.18, "lines")
  ) +
  scale_x_continuous(
    breaks = c(2.0, 2.5, 3.0, 3.5, 4.0),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "C. Cognition preservation exceeds the matched-null distribution",
    subtitle = "Client-level cognition association ranking",
    x = "Mean cognition preservation score\nacross sampled null sets",
    y = "Null iterations"
  ) +
  coord_cartesian(xlim = c(xmin_cognition, xmax_cognition), clip = "off") +
  theme(
    axis.title.x = element_text(size = 9.0, lineheight = 0.92),
    plot.margin = margin(8, 18, 8, 8)
  )

############################################################
## 6. Panel D
############################################################

obs_agora <- figure_tbl %>%
  filter(metric == "agora_fraction") %>%
  pull(observed)

null_agora <- figure_tbl %>%
  filter(metric == "agora_fraction") %>%
  pull(null_mean)

agora_p <- figure_tbl %>%
  filter(metric == "agora_fraction") %>%
  pull(p_label)
agora_p_short <- sub("^empirical ", "", agora_p)

require_columns(null_results, "null_agora_fraction", "Figure 5 Panel D null results")

agora_x_min_data <- min(null_results$null_agora_fraction, na.rm = TRUE)
agora_x_max_data <- max(c(null_results$null_agora_fraction, obs_agora), na.rm = TRUE)
agora_x_span <- agora_x_max_data - agora_x_min_data
xmin_agora <- agora_x_min_data - 0.08 * agora_x_span
xmax_agora <- max(agora_x_max_data + 0.16 * agora_x_span, obs_agora + 0.035)
ymax_agora <- max(ggplot_build(
  ggplot(null_results, aes(x = null_agora_fraction)) +
    geom_histogram(bins = 26)
)$data[[1]]$count)
agora_label_x <- obs_agora + 0.025 * (xmax_agora - xmin_agora)
agora_label <- paste0(
  "Observed mean",
  "\nObs ", number(obs_agora, accuracy = 0.001),
  "\nNull ", number(null_agora, accuracy = 0.001),
  "\n", agora_p_short
)

panel_D <- ggplot(null_results, aes(x = null_agora_fraction)) +
  geom_histogram(
    bins = 26,
    fill = null_fill,
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = obs_agora,
    color = obs_col,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  annotate(
    "label",
    x = agora_label_x,
    y = ymax_agora * 0.90,
    label = agora_label,
    hjust = 0,
    vjust = 1,
    size = 2.16,
    lineheight = 0.92,
    fill = "white",
    label.padding = unit(0.18, "lines")
  ) +
  scale_x_continuous(
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "D. Agora target enrichment exceeds the matched-null distribution",
    subtitle = "External AD target nomination enrichment",
    x = "Fraction of Agora nominated targets\nacross sampled null sets",
    y = "Null iterations"
  ) +
  coord_cartesian(xlim = c(xmin_agora, xmax_agora), clip = "off") +
  theme(
    axis.title.x = element_text(size = 9.0, lineheight = 0.92),
    plot.margin = margin(8, 18, 8, 10)
  )

############################################################
## 7. Assemble final figure
############################################################

top_row <- panel_A + panel_B +
  plot_layout(widths = c(1.18, 1.02))

bottom_row <- panel_C + panel_D +
  plot_layout(widths = c(1.10, 1.10))

final_fig <- (top_row / bottom_row) +
  plot_layout(widths = c(1, 1), heights = c(1.00, 1.15)) +
  plot_annotation(
    title = "Hsp60/10 clients define a selectively AD- and cognition-relevant mitochondrial subnetwork",
    subtitle = paste0(
      "Against 10,000 abundance-matched non-client mitochondrial null sets, ",
      "Hsp60/10 clients show stronger late-stage protein decline, stronger Braak coupling, ",
      "higher cognition preservation support, and greater Agora target enrichment."
    ),
    caption = paste0(
      "Matched background = detected non-Hsp60/10 mitochondrial proteins. ",
      "Bottom panels: dashed blue line marks the observed Hsp60/10 client mean. ",
      "Empirical p-values are bounded by 10,000 null iterations."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(size = 9.7, hjust = 0, color = "gray25", margin = margin(t = 6)),
      plot.margin = margin(12, 14, 10, 14)
    )
  )

print(final_fig)

############################################################
## 8. Save
############################################################

ggsave(
  file.path(output_dir, "FIG_Hsp60_10_matched_mito_null_COGNITION.png"),
  final_fig,
  width = 15.5,
  height = 10.5,
  dpi = 450
)

ggsave(
  file.path(output_dir, "FIG_Hsp60_10_matched_mito_null_COGNITION.pdf"),
  final_fig,
  width = 15.5,
  height = 10.5
)

cat("\nSaved final figure to:\n")
cat(file.path(output_dir, "FIG_Hsp60_10_matched_mito_null_COGNITION.png"), "\n")
cat(file.path(output_dir, "FIG_Hsp60_10_matched_mito_null_COGNITION.pdf"), "\n")
