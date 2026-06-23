
############################################################
## 07_build_matched_null_objects.R
## Build upstream objects required by exact Figure 5 plot block.
## Null loop is migrated from ROSMAP_eQTL_snRNAseq.Rmd; source tables
## are rebuilt from adjusted all-client/protein objects.
############################################################

require_objects(
  c("all_hsp60_10_client_tbl", "prot_mat_raw", "prot_meta_adj", "protein_covars", "mitocarta_genes"),
  context = "07_build_matched_null_objects.R"
)

if (!file.exists(cfg$agora_target_file)) {
  stop("Missing AGORA nominated target file required for Figure 5: ", cfg$agora_target_file,
       "\nUpdate cfg$agora_target_file in 00_config.R, or provide agora_targets in the environment.", call. = FALSE)
}

agora_targets <- readr::read_csv(cfg$agora_target_file, show_col_types = FALSE) |>
  janitor::clean_names()

gene_col_agora <- intersect(
  c("gene", "gene_symbol", "symbol", "hgnc_symbol", "target", "target_gene", "gene_name"),
  colnames(agora_targets)
)[1]

if (is.na(gene_col_agora)) {
  stop(
    "Could not find gene column in AGORA target file. Available columns: ",
    available_cols_msg(agora_targets),
    call. = FALSE
  )
}

agora_target_symbols <- agora_targets |>
  transmute(gene = canonical_gene_symbol(.data[[gene_col_agora]])) |>
  filter(!is.na(gene), gene != "") |>
  distinct() |>
  pull(gene)

hsp_symbols <- clean_gene_symbols(all_hsp60_10_client_tbl$gene)
background_symbols <- setdiff(intersect(clean_gene_symbols(mitocarta_genes), rownames(prot_mat_raw)), hsp_symbols)
hsp_symbols_detected <- intersect(hsp_symbols, rownames(prot_mat_raw))
all_null_genes <- union(hsp_symbols_detected, background_symbols)

## Stage metadata for adjusted collapse metric.
prot_stage_meta_null <- prot_meta_adj |>
  mutate(
    stage = dplyr::recode(as.character(EmoryStrictDx.2019), Control = "Control", AsymAD = "Early AD", AD = "AD", .default = as.character(EmoryStrictDx.2019)),
    stage = factor(stage, levels = c("Control", "Early AD", "AD"))
  ) |>
  filter(SampleID %in% colnames(prot_mat_raw), !is.na(stage))

calc_background_metric_one <- function(g) {
  y <- as.numeric(prot_mat_raw[g, prot_stage_meta_null$SampleID])
  tmp <- tibble(stage = prot_stage_meta_null$stage, y = y)
  early <- tmp$y[tmp$stage == "Early AD"]
  ad <- tmp$y[tmp$stage == "AD"]
  ctrl <- tmp$y[tmp$stage == "Control"]
  protein_late_effect <- mean(ad, na.rm = TRUE) - mean(early, na.rm = TRUE)
  tibble(
    gene = g,
    protein_late_effect = protein_late_effect,
    protein_collapse_magnitude = ifelse(is.finite(protein_late_effect) && protein_late_effect < 0, abs(protein_late_effect), 0),
    matching_abundance = mean(c(early, ctrl, ad), na.rm = TRUE)
  )
}

background_metrics <- purrr::map_dfr(all_null_genes, calc_background_metric_one)

background_braak <- fit_adjusted_braak_beta(
  mat = prot_mat_raw,
  meta_df = prot_meta_adj,
  sample_col = "SampleID",
  genes = all_null_genes,
  covars = protein_covars,
  adjust_for_cerad = cfg$adjust_braak_for_cerad,
  min_n = cfg$min_n_gene_model
) |>
  transmute(gene, inverse_braak_magnitude = adjusted_inverse_braak_beta)

null_input_tbl <- tibble(gene = all_null_genes) |>
  mutate(is_hsp60_10_client = gene %in% hsp_symbols_detected) |>
  checked_left_join(background_metrics, by = "gene", label = "null genes to background collapse metrics") |>
  checked_left_join(background_braak, by = "gene", label = "null genes to adjusted Braak metrics") |>
  mutate(
    agora_nominated_target = gene %in% agora_target_symbols,
    collapse_rank = percentile01(protein_collapse_magnitude) / 100,
    inverse_braak_rank = percentile01(inverse_braak_magnitude) / 100,
    pathology_vulnerability_score = rowMeans(cbind(collapse_rank, inverse_braak_rank), na.rm = TRUE)
  )

message("Null input dimensions before complete-metric filter: ", paste(dim(null_input_tbl), collapse = " x "))
null_input_tbl <- null_input_tbl |>
  filter(
    !is.na(protein_collapse_magnitude),
    !is.na(inverse_braak_magnitude),
    !is.na(pathology_vulnerability_score),
    !is.na(matching_abundance)
  )
message("Null input dimensions after complete-metric filter: ", paste(dim(null_input_tbl), collapse = " x "))

hsp_null_tbl <- null_input_tbl %>% filter(is_hsp60_10_client)
background_null_pool <- null_input_tbl %>% filter(!is_hsp60_10_client)

cat("\nNull input summary:\n")
cat("Hsp60/10 clients:", nrow(hsp_null_tbl), "\n")
cat("Background non-Hsp mitochondrial proteins:", nrow(background_null_pool), "\n")
cat("AGORA targets in Hsp60/10:", sum(hsp_null_tbl$agora_nominated_target), "\n")
cat("AGORA targets in background:", sum(background_null_pool$agora_nominated_target), "\n")

if (nrow(hsp_null_tbl) < 50) stop("Too few Hsp60/10 clients for null.")
if (nrow(background_null_pool) < 50) stop("Too few non-Hsp mitochondrial background proteins for null.")

############################################################
## 11. Abundance-matched mitochondrial null test
############################################################

set.seed(1300)

n_iter <- cfg$null_n_iter
n_hsp <- nrow(hsp_null_tbl)

cat("\nRunning abundance-matched mitochondrial null...\n")
cat("Iterations:", n_iter, "\n")
cat("Observed Hsp60/10 n:", n_hsp, "\n")
cat("Background pool n:", nrow(background_null_pool), "\n")

null_input_tbl <- null_input_tbl %>%
  mutate(abundance_bin = ntile(matching_abundance, 10))

hsp_null_tbl <- null_input_tbl %>% filter(is_hsp60_10_client)
background_null_pool <- null_input_tbl %>% filter(!is_hsp60_10_client)

observed_stats <- hsp_null_tbl %>%
  summarise(
    observed_mean_collapse = mean(protein_collapse_magnitude, na.rm = TRUE),
    observed_mean_inverse_braak = mean(inverse_braak_magnitude, na.rm = TRUE),
    observed_mean_pathology_score = mean(pathology_vulnerability_score, na.rm = TRUE),
    observed_agora_fraction = mean(agora_nominated_target, na.rm = TRUE),
    observed_agora_n = sum(agora_nominated_target, na.rm = TRUE),
    n_hsp = n()
  )

sample_abundance_matched_null <- function(hsp_tbl, bg_tbl) {
  sampled <- hsp_tbl %>%
    count(abundance_bin, name = "n_needed") %>%
    group_split(abundance_bin) %>%
    purrr::map_dfr(function(bin_df) {
      this_bin <- bin_df$abundance_bin[1]
      n_needed <- bin_df$n_needed[1]
      bg_bin <- bg_tbl %>% filter(abundance_bin == this_bin)
      if (nrow(bg_bin) == 0) return(tibble())
      bg_bin %>% slice_sample(n = n_needed, replace = nrow(bg_bin) < n_needed)
    })
  sampled
}

null_results <- purrr::map_dfr(seq_len(n_iter), function(i) {
  sampled_bg <- sample_abundance_matched_null(hsp_tbl = hsp_null_tbl, bg_tbl = background_null_pool)
  tibble(
    iter = i,
    null_mean_collapse = mean(sampled_bg$protein_collapse_magnitude, na.rm = TRUE),
    null_mean_inverse_braak = mean(sampled_bg$inverse_braak_magnitude, na.rm = TRUE),
    null_mean_pathology_score = mean(sampled_bg$pathology_vulnerability_score, na.rm = TRUE),
    null_agora_fraction = mean(sampled_bg$agora_nominated_target, na.rm = TRUE),
    null_agora_n = sum(sampled_bg$agora_nominated_target, na.rm = TRUE),
    null_n = nrow(sampled_bg)
  )
})

obs <- observed_stats

null_summary <- tibble(
  metric = c("protein_collapse_magnitude", "inverse_braak_magnitude", "pathology_vulnerability_score", "agora_fraction"),
  observed = c(obs$observed_mean_collapse, obs$observed_mean_inverse_braak, obs$observed_mean_pathology_score, obs$observed_agora_fraction),
  null_mean = c(mean(null_results$null_mean_collapse, na.rm = TRUE), mean(null_results$null_mean_inverse_braak, na.rm = TRUE), mean(null_results$null_mean_pathology_score, na.rm = TRUE), mean(null_results$null_agora_fraction, na.rm = TRUE)),
  null_sd = c(sd(null_results$null_mean_collapse, na.rm = TRUE), sd(null_results$null_mean_inverse_braak, na.rm = TRUE), sd(null_results$null_mean_pathology_score, na.rm = TRUE), sd(null_results$null_agora_fraction, na.rm = TRUE)),
  empirical_p_greater = c(mean(null_results$null_mean_collapse >= obs$observed_mean_collapse, na.rm = TRUE), mean(null_results$null_mean_inverse_braak >= obs$observed_mean_inverse_braak, na.rm = TRUE), mean(null_results$null_mean_pathology_score >= obs$observed_mean_pathology_score, na.rm = TRUE), mean(null_results$null_agora_fraction >= obs$observed_agora_fraction, na.rm = TRUE))
) %>%
  mutate(
    z_score = (observed - null_mean) / null_sd,
    empirical_p_greater = pmax(empirical_p_greater, 1 / n_iter)
  )

output_dir <- file.path(cfg$plot_dir, "main_fig5_matched_null_specificity_exact")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tbl(null_summary, "figure5_null_summary_COVARIATE_ADJUSTED")
write_tbl(null_results, "figure5_null_results_COVARIATE_ADJUSTED")
write_tbl(null_input_tbl, "figure5_null_input_COVARIATE_ADJUSTED")
save_obj(null_results, "figure5_null_results")
save_obj(null_summary, "figure5_null_summary")
save_obj(observed_stats, "figure5_observed_stats")
message("Loaded 07_build_matched_null_objects.R")

