############################################################
## 20_complementary_direction_sensitivity.R
##
## Reviewer 2 complementary-direction sensitivity analysis.
##
## Reuses the canonical Figure 5 matched-null input produced by
## R/07_build_matched_null_objects.R so that the protein universe,
## abundance matching, abundance bins, and signed AD-versus-MCI
## effects are identical to the validated production analysis.
##
## Metrics:
##   signed effect              = mean(AD) - mean(MCI)
##   lower-abundance magnitude  = max(MCI - AD, 0)
##   higher-abundance magnitude = max(AD - MCI, 0)
##   absolute difference        = abs(AD - MCI)
##
## Additional analyses:
##   - upper- and lower-tail empirical P values
##   - Fisher exact test of lower-vs-higher direction proportions
##   - conventional AD-versus-MCI protein direction counts from Fig. S7
##   - canonical Figure 5 reproduction check
##
## Writes ONLY to:
##   outputs/reviewer_sensitivity/complementary_direction/
##
## Run from repository root:
##   Rscript R/reviewer_sensitivity/20_complementary_direction_sensitivity.R
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

source("R/00_config.R")

############################################################
## 0. Setup
############################################################

out_dir <- file.path(
  cfg$project_dir,
  "outputs",
  "reviewer_sensitivity",
  "complementary_direction"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

n_iter <- cfg$null_n_iter
seed <- 1300
reporting_floor <- 1 / n_iter

canonical_input_file <- file.path(
  cfg$table_dir,
  "figure5_null_input_COVARIATE_ADJUSTED.csv"
)

canonical_summary_file <- file.path(
  cfg$table_dir,
  "figure5_null_summary_COVARIATE_ADJUSTED.csv"
)

if (!file.exists(canonical_input_file)) {
  stop(
    "Canonical Figure 5 input not found: ",
    canonical_input_file,
    "\nRun the validated main pipeline through ",
    "R/07_build_matched_null_objects.R first.",
    call. = FALSE
  )
}

x <- read_csv(
  canonical_input_file,
  show_col_types = FALSE
)

required_cols <- c(
  "gene",
  "is_hsp60_10_client",
  "protein_late_effect",
  "protein_late_decline_magnitude",
  "matching_abundance",
  "abundance_bin"
)

missing_cols <- setdiff(
  required_cols,
  names(x)
)

if (length(missing_cols) > 0) {
  stop(
    "Canonical Figure 5 input is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

############################################################
## 1. Derive complementary metrics from the same signed effect
############################################################

x <- x |>
  mutate(
    signed_ad_minus_mci = .data$protein_late_effect,

    lower_abundance_magnitude =
      pmax(-.data$signed_ad_minus_mci, 0),

    higher_abundance_magnitude =
      pmax(.data$signed_ad_minus_mci, 0),

    absolute_difference =
      abs(.data$signed_ad_minus_mci),

    direction = case_when(
      .data$signed_ad_minus_mci < 0 ~ "lower_in_AD",
      .data$signed_ad_minus_mci > 0 ~ "higher_in_AD",
      TRUE ~ "equal"
    )
  )

############################################################
## 2. Verify exact reconstruction of primary metric
############################################################

max_lower_diff <- max(
  abs(
    x$lower_abundance_magnitude -
      x$protein_late_decline_magnitude
  ),
  na.rm = TRUE
)

if (
  !is.finite(max_lower_diff) ||
    max_lower_diff > 1e-12
) {
  stop(
    "Reconstructed lower-abundance metric does not match ",
    "the canonical Figure 5 metric. Maximum absolute difference = ",
    max_lower_diff,
    call. = FALSE
  )
}

hsp_tbl <- x |>
  filter(.data$is_hsp60_10_client)

background_tbl <- x |>
  filter(!.data$is_hsp60_10_client)

message(
  "Hsp60/10 clients: ",
  nrow(hsp_tbl)
)

message(
  "Non-client mitochondrial background: ",
  nrow(background_tbl)
)

if (
  nrow(hsp_tbl) != 306 ||
    nrow(background_tbl) != 609
) {
  warning(
    "Expected canonical Figure 5 universe of 306 clients ",
    "and 609 background proteins; found ",
    nrow(hsp_tbl),
    " and ",
    nrow(background_tbl),
    ".",
    call. = FALSE
  )
}

############################################################
## 3. Direction counts
############################################################

direction_counts <- x |>
  mutate(
    set = if_else(
      .data$is_hsp60_10_client,
      "Hsp60/10 clients",
      "Non-client mitochondrial background"
    )
  ) |>
  count(
    .data$set,
    .data$direction,
    name = "n"
  ) |>
  group_by(.data$set) |>
  mutate(
    total = sum(.data$n),
    fraction = .data$n / .data$total
  ) |>
  ungroup()

n_hsp_lower <- sum(
  hsp_tbl$signed_ad_minus_mci < 0,
  na.rm = TRUE
)

n_hsp_higher <- sum(
  hsp_tbl$signed_ad_minus_mci > 0,
  na.rm = TRUE
)

n_bg_lower <- sum(
  background_tbl$signed_ad_minus_mci < 0,
  na.rm = TRUE
)

n_bg_higher <- sum(
  background_tbl$signed_ad_minus_mci > 0,
  na.rm = TRUE
)

############################################################
## 4. Fisher exact test of directional enrichment
############################################################

fisher_matrix <- matrix(
  c(
    n_hsp_lower,
    n_hsp_higher,
    n_bg_lower,
    n_bg_higher
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    set = c(
      "Hsp60_10_clients",
      "nonclient_mito"
    ),
    direction = c(
      "lower_in_AD",
      "higher_in_AD"
    )
  )
)

fisher_result <- fisher.test(
  fisher_matrix,
  alternative = "two.sided"
)

fisher_summary <- tibble(
  hsp_lower = n_hsp_lower,
  hsp_higher = n_hsp_higher,

  hsp_lower_fraction =
    n_hsp_lower /
    (n_hsp_lower + n_hsp_higher),

  background_lower = n_bg_lower,
  background_higher = n_bg_higher,

  background_lower_fraction =
    n_bg_lower /
    (n_bg_lower + n_bg_higher),

  odds_ratio =
    unname(fisher_result$estimate),

  conf_low =
    fisher_result$conf.int[[1]],

  conf_high =
    fisher_result$conf.int[[2]],

  p_value_two_sided =
    fisher_result$p.value
)

############################################################
## 5. Observed client-set statistics
############################################################

observed <- hsp_tbl |>
  summarise(
    n_hsp = n(),

    mean_signed_effect =
      mean(
        .data$signed_ad_minus_mci,
        na.rm = TRUE
      ),

    mean_lower_abundance =
      mean(
        .data$lower_abundance_magnitude,
        na.rm = TRUE
      ),

    mean_higher_abundance =
      mean(
        .data$higher_abundance_magnitude,
        na.rm = TRUE
      ),

    mean_absolute_difference =
      mean(
        .data$absolute_difference,
        na.rm = TRUE
      ),

    n_lower_in_ad =
      sum(
        .data$signed_ad_minus_mci < 0,
        na.rm = TRUE
      ),

    n_higher_in_ad =
      sum(
        .data$signed_ad_minus_mci > 0,
        na.rm = TRUE
      ),

    n_equal =
      sum(
        .data$signed_ad_minus_mci == 0,
        na.rm = TRUE
      )
  )

############################################################
## 6. Exact canonical Figure 5 abundance-matched resampling
############################################################

sample_abundance_matched_null <- function(
  hsp_tbl,
  bg_tbl
) {
  hsp_tbl |>
    count(
      .data$abundance_bin,
      name = "n_needed"
    ) |>
    group_split(
      .data$abundance_bin
    ) |>
    map_dfr(
      function(bin_df) {

        this_bin <-
          bin_df$abundance_bin[[1]]

        n_needed <-
          bin_df$n_needed[[1]]

        bg_bin <- bg_tbl |>
          filter(
            .data$abundance_bin ==
              this_bin
          )

        if (nrow(bg_bin) == 0) {
          stop(
            "No background proteins in abundance bin ",
            this_bin,
            ".",
            call. = FALSE
          )
        }

        bg_bin |>
          slice_sample(
            n = n_needed,
            replace =
              nrow(bg_bin) <
              n_needed
          )
      }
    )
}

set.seed(seed)

message(
  "Running ",
  n_iter,
  " abundance-matched null sets with seed ",
  seed,
  "..."
)

null_results <- map_dfr(
  seq_len(n_iter),
  function(i) {

    sampled_bg <-
      sample_abundance_matched_null(
        hsp_tbl,
        background_tbl
      )

    tibble(
      iter = i,

      null_mean_signed_effect =
        mean(
          sampled_bg$signed_ad_minus_mci,
          na.rm = TRUE
        ),

      null_mean_lower_abundance =
        mean(
          sampled_bg$lower_abundance_magnitude,
          na.rm = TRUE
        ),

      null_mean_higher_abundance =
        mean(
          sampled_bg$higher_abundance_magnitude,
          na.rm = TRUE
        ),

      null_mean_absolute_difference =
        mean(
          sampled_bg$absolute_difference,
          na.rm = TRUE
        )
    )
  }
)

############################################################
## 7. Upper- and lower-tail empirical P values
############################################################

metric_spec <- tribble(
  ~metric,
  ~label,
  ~observed_col,
  ~null_col,
  ~primary_tail,

  "lower_abundance_magnitude",
  "Lower abundance in AD versus MCI",
  "mean_lower_abundance",
  "null_mean_lower_abundance",
  "upper",

  "higher_abundance_magnitude",
  "Higher abundance in AD versus MCI",
  "mean_higher_abundance",
  "null_mean_higher_abundance",
  "lower",

  "absolute_difference",
  "Absolute AD-versus-MCI difference",
  "mean_absolute_difference",
  "null_mean_absolute_difference",
  "upper"
)

summarize_one <- function(
  metric,
  label,
  observed_col,
  null_col,
  primary_tail
) {
  obs <-
    observed[[observed_col]][[1]]

  null <-
    null_results[[null_col]]

  null <-
    null[is.finite(null)]

  n_finite <-
    length(null)

  upper_count <-
    sum(null >= obs)

  lower_count <-
    sum(null <= obs)

  p_upper_raw <-
    upper_count /
    n_finite

  p_lower_raw <-
    lower_count /
    n_finite

  p_upper_reported <-
    max(
      p_upper_raw,
      reporting_floor
    )

  p_lower_reported <-
    max(
      p_lower_raw,
      reporting_floor
    )

  primary_p_raw <-
    if (primary_tail == "upper") {
      p_upper_raw
    } else {
      p_lower_raw
    }

  primary_p_reported <-
    if (primary_tail == "upper") {
      p_upper_reported
    } else {
      p_lower_reported
    }

  tibble(
    metric = metric,
    label = label,
    primary_tail = primary_tail,

    observed = obs,

    null_mean =
      mean(null),

    null_sd =
      sd(null),

    observed_to_null_ratio =
      obs / mean(null),

    z_score =
      (obs - mean(null)) /
      sd(null),

    upper_tail_extreme_count =
      upper_count,

    empirical_p_upper_raw =
      p_upper_raw,

    empirical_p_upper_reported =
      p_upper_reported,

    lower_tail_extreme_count =
      lower_count,

    empirical_p_lower_raw =
      p_lower_raw,

    empirical_p_lower_reported =
      p_lower_reported,

    primary_empirical_p_raw =
      primary_p_raw,

    primary_empirical_p_reported =
      primary_p_reported
  )
}

summary_tbl <- pmap_dfr(
  metric_spec,
  summarize_one
)

############################################################
## 8. Reproduce canonical Figure 5 lower-abundance result
############################################################

validation_tbl <- tibble()

if (file.exists(canonical_summary_file)) {

  canonical <- read_csv(
    canonical_summary_file,
    show_col_types = FALSE
  ) |>
    filter(
      .data$metric ==
        "protein_late_decline_magnitude"
    )

  if (nrow(canonical) == 1) {

    lower_row <- summary_tbl |>
      filter(
        .data$metric ==
          "lower_abundance_magnitude"
      )

    validation_tbl <- tibble(
      check = c(
        "observed lower-abundance mean",
        "matched-null lower-abundance mean"
      ),

      sensitivity_value = c(
        lower_row$observed,
        lower_row$null_mean
      ),

      canonical_value = c(
        canonical$observed,
        canonical$null_mean
      )
    ) |>
      mutate(
        abs_difference =
          abs(
            .data$sensitivity_value -
              .data$canonical_value
          ),

        status =
          if_else(
            .data$abs_difference <= 1e-10,
            "PASS",
            "FAIL"
          )
      )

    if (
      any(
        validation_tbl$status !=
          "PASS"
      )
    ) {
      stop(
        "Sensitivity analysis failed to reproduce ",
        "canonical Figure 5 lower-abundance values.",
        call. = FALSE
      )
    }
  }
}

############################################################
## 9. Conventional differential direction audit from Fig. S7
############################################################

sf7_file <- file.path(
  cfg$project_dir,
  "outputs",
  "supplemental_figures",
  "tables",
  "supplementary_figure_7_conventional_differential",
  "Supplementary_Table_conventional_RNA_protein_clients.csv"
)

sf7_direction <- tibble()

if (file.exists(sf7_file)) {

  sf7 <- read_csv(
    sf7_file,
    show_col_types = FALSE
  )

  needed <- c(
    "contrast",
    "protein_effect",
    "protein_significant_fdr05"
  )

  if (
    all(
      needed %in%
        names(sf7)
    )
  ) {

    sf7_direction <- sf7 |>
      filter(
        as.character(
          .data$contrast
        ) ==
          "AD_vs_MCI"
      ) |>
      summarise(
        n_shared_clients =
          n(),

        n_protein_fdr_significant =
          sum(
            .data$protein_significant_fdr05,
            na.rm = TRUE
          ),

        n_fdr_lower_in_ad =
          sum(
            .data$protein_significant_fdr05 &
              .data$protein_effect < 0,
            na.rm = TRUE
          ),

        n_fdr_higher_in_ad =
          sum(
            .data$protein_significant_fdr05 &
              .data$protein_effect > 0,
            na.rm = TRUE
          ),

        n_all_lower_in_ad =
          sum(
            .data$protein_effect < 0,
            na.rm = TRUE
          ),

        n_all_higher_in_ad =
          sum(
            .data$protein_effect > 0,
            na.rm = TRUE
          )
      )
  } else {
    warning(
      "Figure S7 table exists but lacks expected columns: ",
      paste(
        setdiff(
          needed,
          names(sf7)
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
} else {
  warning(
    "Figure S7 conventional differential table not found; ",
    "direction audit skipped.",
    call. = FALSE
  )
}

############################################################
## 10. Manuscript-ready consolidated summary
############################################################

lower_row <- summary_tbl |>
  filter(
    .data$metric ==
      "lower_abundance_magnitude"
  )

higher_row <- summary_tbl |>
  filter(
    .data$metric ==
      "higher_abundance_magnitude"
  )

absolute_row <- summary_tbl |>
  filter(
    .data$metric ==
      "absolute_difference"
  )

manuscript_summary <- tibble(
  statistic = c(
    "Hsp60/10 clients lower in AD",
    "Hsp60/10 clients higher in AD",
    "Hsp60/10 lower-in-AD fraction",
    "Background proteins lower in AD",
    "Background proteins higher in AD",
    "Background lower-in-AD fraction",
    "Direction Fisher odds ratio",
    "Direction Fisher 95% CI low",
    "Direction Fisher 95% CI high",
    "Direction Fisher two-sided P",

    "Lower-abundance observed mean",
    "Lower-abundance null mean",
    "Lower-abundance observed/null ratio",
    "Lower-abundance upper-tail empirical P",

    "Higher-abundance observed mean",
    "Higher-abundance null mean",
    "Higher-abundance observed/null ratio",
    "Higher-abundance lower-tail empirical P",

    "Absolute-difference observed mean",
    "Absolute-difference null mean",
    "Absolute-difference observed/null ratio",
    "Absolute-difference upper-tail empirical P"
  ),

  value = c(
    n_hsp_lower,
    n_hsp_higher,
    fisher_summary$hsp_lower_fraction,

    n_bg_lower,
    n_bg_higher,
    fisher_summary$background_lower_fraction,

    fisher_summary$odds_ratio,
    fisher_summary$conf_low,
    fisher_summary$conf_high,
    fisher_summary$p_value_two_sided,

    lower_row$observed,
    lower_row$null_mean,
    lower_row$observed_to_null_ratio,
    lower_row$empirical_p_upper_reported,

    higher_row$observed,
    higher_row$null_mean,
    higher_row$observed_to_null_ratio,
    higher_row$empirical_p_lower_reported,

    absolute_row$observed,
    absolute_row$null_mean,
    absolute_row$observed_to_null_ratio,
    absolute_row$empirical_p_upper_reported
  )
)

if (nrow(sf7_direction) == 1) {

  manuscript_summary <- bind_rows(
    manuscript_summary,

    tibble(
      statistic = c(
        "Shared clients in conventional AD-versus-MCI analysis",
        "FDR-supported protein clients",
        "FDR-supported protein clients lower in AD",
        "FDR-supported protein clients higher in AD",
        "All shared protein clients lower in AD",
        "All shared protein clients higher in AD"
      ),

      value = c(
        sf7_direction$n_shared_clients,
        sf7_direction$n_protein_fdr_significant,
        sf7_direction$n_fdr_lower_in_ad,
        sf7_direction$n_fdr_higher_in_ad,
        sf7_direction$n_all_lower_in_ad,
        sf7_direction$n_all_higher_in_ad
      )
    )
  )
}

############################################################
## 11. Save source and summary outputs
############################################################

write_csv(
  direction_counts,
  file.path(
    out_dir,
    "direction_counts.csv"
  )
)

write_csv(
  fisher_summary,
  file.path(
    out_dir,
    "direction_fisher_exact.csv"
  )
)

write_csv(
  observed,
  file.path(
    out_dir,
    "observed_client_statistics.csv"
  )
)

write_csv(
  null_results,
  file.path(
    out_dir,
    "matched_null_permutations.csv"
  )
)

write_csv(
  summary_tbl,
  file.path(
    out_dir,
    "matched_null_summary.csv"
  )
)

write_csv(
  validation_tbl,
  file.path(
    out_dir,
    "canonical_figure5_validation.csv"
  )
)

write_csv(
  manuscript_summary,
  file.path(
    out_dir,
    "manuscript_ready_summary.csv"
  )
)

if (nrow(sf7_direction) == 1) {
  write_csv(
    sf7_direction,
    file.path(
      out_dir,
      "conventional_DE_AD_vs_MCI_direction_summary.csv"
    )
  )
}

write_csv(
  x |>
    select(
      gene,
      is_hsp60_10_client,
      abundance_bin,
      matching_abundance,
      signed_ad_minus_mci,
      lower_abundance_magnitude,
      higher_abundance_magnitude,
      absolute_difference,
      direction
    ),
  file.path(
    out_dir,
    "gene_level_direction_metrics.csv"
  )
)

############################################################
## 12. Diagnostic plot
############################################################

plot_df <- null_results |>
  select(
    iter,
    lower_abundance_magnitude =
      null_mean_lower_abundance,
    higher_abundance_magnitude =
      null_mean_higher_abundance,
    absolute_difference =
      null_mean_absolute_difference
  ) |>
  pivot_longer(
    cols = -iter,
    names_to = "metric",
    values_to = "null_mean"
  ) |>
  left_join(
    summary_tbl |>
      select(
        metric,
        label,
        observed
      ),
    by = "metric"
  )

p <- ggplot(
  plot_df,
  aes(
    x = .data$null_mean
  )
) +
  geom_histogram(
    bins = 50,
    fill = "grey80",
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    aes(
      xintercept =
        .data$observed
    ),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(
    ~label,
    scales = "free_x",
    ncol = 1
  ) +
  labs(
    title =
      "Complementary-direction abundance-matched sensitivity",

    subtitle =
      paste0(
        "Hsp60/10 clients versus abundance-matched ",
        "non-client mitochondrial sets; ",
        format(
          n_iter,
          big.mark = ","
        ),
        " null sets"
      ),

    x = "Null-set mean",
    y = "Count"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    strip.text =
      element_text(
        face = "bold"
      ),

    plot.title =
      element_text(
        face = "bold"
      )
  )

ggsave(
  file.path(
    out_dir,
    "complementary_direction_sensitivity.pdf"
  ),
  p,
  width = 7,
  height = 8,
  units = "in",
  device = grDevices::pdf,
  bg = "white",
  useDingbats = FALSE
)

ggsave(
  file.path(
    out_dir,
    "complementary_direction_sensitivity.png"
  ),
  p,
  width = 7,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

############################################################
## 13. Provenance
############################################################

git_head <- tryCatch(
  system2(
    "git",
    c(
      "rev-parse",
      "HEAD"
    ),
    stdout = TRUE,
    stderr = FALSE
  )[1],
  error =
    function(e) NA_character_
)

input_md5 <-
  unname(
    tools::md5sum(
      canonical_input_file
    )
  )

writeLines(
  c(
    "analysis=Reviewer 2 complementary-direction sensitivity",
    paste0(
      "git_head=",
      git_head
    ),
    paste0(
      "canonical_input=",
      canonical_input_file
    ),
    paste0(
      "canonical_input_md5=",
      input_md5
    ),
    paste0(
      "seed=",
      seed
    ),
    paste0(
      "n_iter=",
      n_iter
    ),
    paste0(
      "reporting_floor=",
      reporting_floor
    ),
    paste0(
      "n_hsp=",
      nrow(hsp_tbl)
    ),
    paste0(
      "n_background=",
      nrow(background_tbl)
    ),
    paste0(
      "max_lower_metric_reconstruction_error=",
      max_lower_diff
    ),
    paste0(
      "timestamp=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S %z"
      )
    )
  ),
  file.path(
    out_dir,
    "provenance.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "session_info.txt"
  )
)

############################################################
## 14. Console summary
############################################################

cat(
  "\n============================================================\n"
)

cat(
  "COMPLEMENTARY-DIRECTION SENSITIVITY COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Direction counts:\n"
)

print(
  direction_counts
)

cat(
  "\nFisher exact direction test:\n"
)

print(
  fisher_summary
)

cat(
  "\nObserved client statistics:\n"
)

print(
  observed
)

cat(
  "\nMatched-null sensitivity summary:\n"
)

print(
  summary_tbl
)

if (nrow(validation_tbl) > 0) {
  cat(
    "\nCanonical Figure 5 reproduction:\n"
  )

  print(
    validation_tbl
  )
}

if (nrow(sf7_direction) > 0) {
  cat(
    "\nConventional AD-versus-MCI DE direction summary:\n"
  )

  print(
    sf7_direction
  )
}

cat(
  "\nManuscript-ready summary:\n"
)

print(
  manuscript_summary,
  n = Inf
)

cat(
  "\nOutputs written to:\n"
)

cat(
  out_dir,
  "\n"
)

cat(
  "\nCanonical Figure 5 and frozen supplemental outputs were not modified.\n"
)
