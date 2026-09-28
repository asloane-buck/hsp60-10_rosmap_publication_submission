############################################################
## 45_validate_Hsp_pathway_PC1_sensitivity.R
##
## Hard validation of Reviewer 2.3 Hsp60/10 RNA
## pathway-score PC1/eigengene sensitivity.
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "Hsp_pathway_PC1_sensitivity"
)

required_files <- c(
  "Hsp60_10_RNA_client_manifest.csv",
  "participant_mean_z_and_PC1_scores.csv",
  "mean_z_vs_PC1_concordance.csv",
  "PCA_variance_explained.csv",
  "PC1_gene_loadings.csv",
  "PC1_loading_summary.csv",
  "stage_score_summary.csv",
  "stage_Wilcoxon_mean_z_vs_PC1.csv",
  "stage_direction_concordance.csv",
  "analysis_provenance.csv",
  "sessionInfo.txt"
)

missing_files <- required_files[
  !file.exists(
    file.path(
      OUTDIR,
      required_files
    )
  )
]

if (length(missing_files) > 0) {
  stop(
    "Missing required PC1 sensitivity outputs: ",
    paste(missing_files, collapse = ", ")
  )
}

client_manifest <- read.csv(
  file.path(
    OUTDIR,
    "Hsp60_10_RNA_client_manifest.csv"
  )
)

scores <- read.csv(
  file.path(
    OUTDIR,
    "participant_mean_z_and_PC1_scores.csv"
  )
)

concordance <- read.csv(
  file.path(
    OUTDIR,
    "mean_z_vs_PC1_concordance.csv"
  )
)

variance_tbl <- read.csv(
  file.path(
    OUTDIR,
    "PCA_variance_explained.csv"
  )
)

loadings <- read.csv(
  file.path(
    OUTDIR,
    "PC1_gene_loadings.csv"
  )
)

loading_summary <- read.csv(
  file.path(
    OUTDIR,
    "PC1_loading_summary.csv"
  )
)

stage_summary <- read.csv(
  file.path(
    OUTDIR,
    "stage_score_summary.csv"
  )
)

stage_tests <- read.csv(
  file.path(
    OUTDIR,
    "stage_Wilcoxon_mean_z_vs_PC1.csv"
  )
)

provenance <- read.csv(
  file.path(
    OUTDIR,
    "analysis_provenance.csv"
  )
)

checks <- data.frame(
  check = character(),
  passed = logical(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_check <- function(name, pass, detail) {
  checks <<- rbind(
    checks,
    data.frame(
      check = name,
      passed = isTRUE(pass),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  )
}


## ---------------------------------------------------------
## Cohort / feature universe
## ---------------------------------------------------------

add_check(
  "participant_n_577",
  nrow(scores) == 577,
  paste("observed =", nrow(scores))
)

add_check(
  "participant_ids_unique",
  anyDuplicated(scores$sample_id) == 0,
  paste("duplicates =", anyDuplicated(scores$sample_id))
)

add_check(
  "client_n_297",
  nrow(client_manifest) == 297,
  paste("observed =", nrow(client_manifest))
)

add_check(
  "client_ids_unique",
  anyDuplicated(client_manifest$gene) == 0,
  paste("duplicates =", anyDuplicated(client_manifest$gene))
)

stage_counts <- table(scores$clinical_stage)

add_check(
  "stage_counts",
  identical(
    as.integer(stage_counts[c("AD", "MCI", "NCI")]),
    c(219L, 158L, 200L)
  ),
  paste(
    names(stage_counts),
    stage_counts,
    collapse = "; "
  )
)


## ---------------------------------------------------------
## Exact reproduction of production score
## ---------------------------------------------------------

mean_z_diff <- concordance$value[
  concordance$metric ==
    "mean_z_max_abs_diff_vs_production"
]

add_check(
  "production_mean_z_exact",
  length(mean_z_diff) == 1 &&
    is.finite(mean_z_diff) &&
    mean_z_diff <= 1e-12,
  paste("max abs diff =", mean_z_diff)
)


## ---------------------------------------------------------
## PCA construction
## ---------------------------------------------------------

pearson_r <- concordance$value[
  concordance$metric ==
    "pearson_r_mean_z_vs_PC1"
]

spearman_rho <- concordance$value[
  concordance$metric ==
    "spearman_rho_mean_z_vs_PC1"
]

pc1_var <- concordance$value[
  concordance$metric ==
    "PC1_variance_explained"
]

add_check(
  "PC1_oriented_positive",
  length(pearson_r) == 1 &&
    is.finite(pearson_r) &&
    pearson_r > 0,
  paste("Pearson r =", pearson_r)
)

add_check(
  "Pearson_expected",
  abs(pearson_r - 0.7911645) < 1e-6,
  paste("Pearson r =", pearson_r)
)

add_check(
  "Spearman_expected",
  abs(spearman_rho - 0.7825094) < 1e-6,
  paste("Spearman rho =", spearman_rho)
)

add_check(
  "PC1_variance_expected",
  abs(pc1_var - 0.2182521) < 1e-6,
  paste("PC1 variance explained =", pc1_var)
)

add_check(
  "PC1_score_centered",
  abs(mean(scores$PC1_eigengene)) < 1e-10,
  paste("mean =", mean(scores$PC1_eigengene))
)

add_check(
  "PC1_score_sd_one",
  abs(sd(scores$PC1_eigengene) - 1) < 1e-10,
  paste("SD =", sd(scores$PC1_eigengene))
)

add_check(
  "loading_n_297",
  nrow(loadings) == 297,
  paste("observed =", nrow(loadings))
)

add_check(
  "loading_summary_204_positive",
  loading_summary$n_positive == 204,
  paste("positive =", loading_summary$n_positive)
)

add_check(
  "loading_summary_93_negative",
  loading_summary$n_negative == 93,
  paste("negative =", loading_summary$n_negative)
)


## ---------------------------------------------------------
## Stage sensitivity
## ---------------------------------------------------------

add_check(
  "four_stage_tests",
  nrow(stage_tests) == 4,
  paste("observed =", nrow(stage_tests))
)

mean_mci_ad <- subset(
  stage_tests,
  score_type == "mean_z" &
    comparison == "MCI_vs_AD"
)

pc1_mci_ad <- subset(
  stage_tests,
  score_type == "PC1_eigengene" &
    comparison == "MCI_vs_AD"
)

mean_nci_mci <- subset(
  stage_tests,
  score_type == "mean_z" &
    comparison == "NCI_vs_MCI"
)

pc1_nci_mci <- subset(
  stage_tests,
  score_type == "PC1_eigengene" &
    comparison == "NCI_vs_MCI"
)

add_check(
  "MCI_to_AD_direction_concordant",
  nrow(mean_mci_ad) == 1 &&
    nrow(pc1_mci_ad) == 1 &&
    mean_mci_ad$delta_mean < 0 &&
    pc1_mci_ad$delta_mean < 0,
  paste(
    "mean-z delta =", mean_mci_ad$delta_mean,
    "; PC1 delta =", pc1_mci_ad$delta_mean
  )
)

add_check(
  "NCI_to_MCI_direction_concordant",
  nrow(mean_nci_mci) == 1 &&
    nrow(pc1_nci_mci) == 1 &&
    mean_nci_mci$delta_mean > 0 &&
    pc1_nci_mci$delta_mean > 0,
  paste(
    "mean-z delta =", mean_nci_mci$delta_mean,
    "; PC1 delta =", pc1_nci_mci$delta_mean
  )
)

add_check(
  "mean_z_MCI_AD_p_expected",
  abs(mean_mci_ad$p_value - 0.001248987) < 1e-8,
  paste("p =", mean_mci_ad$p_value)
)

add_check(
  "PC1_MCI_AD_p_expected",
  abs(pc1_mci_ad$p_value - 0.050758026) < 1e-8,
  paste("p =", pc1_mci_ad$p_value)
)


## ---------------------------------------------------------
## Provenance
## ---------------------------------------------------------

prov <- setNames(
  provenance$value,
  provenance$item
)

add_check(
  "provenance_RNA_samples",
  identical(as.character(prov["RNA_samples"]), "577"),
  paste("value =", prov["RNA_samples"])
)

add_check(
  "provenance_client_n",
  identical(
    as.character(
      prov["RNA_detected_Hsp60_10_clients"]
    ),
    "297"
  ),
  paste(
    "value =",
    prov["RNA_detected_Hsp60_10_clients"]
  )
)

add_check(
  "provenance_PCA_input",
  identical(
    as.character(prov["PCA_input"]),
    "same gene-wise z-scored client matrix"
  ),
  paste("value =", prov["PCA_input"])
)


## ---------------------------------------------------------
## Warning audit
## ---------------------------------------------------------

log_file <- file.path(
  "outputs",
  "reviewer_revisions",
  "Hsp_pathway_PC1_sensitivity_run.log"
)

if (file.exists(log_file)) {
  log_text <- readLines(
    log_file,
    warn = FALSE
  )

  warning_hits <- grep(
    "^Warning|Warning in|warning:",
    log_text,
    value = TRUE,
    ignore.case = TRUE
  )

  add_check(
    "run_log_no_warnings",
    length(warning_hits) == 0,
    paste("warning lines =", length(warning_hits))
  )
} else {
  add_check(
    "run_log_exists",
    FALSE,
    paste("Missing:", log_file)
  )
}


## ---------------------------------------------------------
## Report
## ---------------------------------------------------------

print(
  checks,
  row.names = FALSE
)

n_pass <- sum(checks$passed)
n_total <- nrow(checks)

cat(
  "\nValidation:",
  n_pass,
  "/",
  n_total,
  "checks passed\n"
)

write.csv(
  checks,
  file.path(
    OUTDIR,
    "PC1_sensitivity_validation_checks.csv"
  ),
  row.names = FALSE
)

if (!all(checks$passed)) {
  failed <- checks$check[
    !checks$passed
  ]

  stop(
    "PC1 sensitivity validation failed: ",
    paste(failed, collapse = ", ")
  )
}

cat(
  "\n============================================================\n"
)
cat(
  "Hsp60/10 PATHWAY PC1 SENSITIVITY VALIDATED\n"
)
cat(
  "============================================================\n"
)
