############################################################
## 54_audit_reporting_requirements.R
##
## Master reviewer-facing reporting audit.
##
## Audits existing validated outputs for:
##   effect estimate, SE, 95% CI, n, P, FDR,
##   and participant-level/distribution source data where needed.
##
## This script DOES NOT recompute scientific analyses.
############################################################

options(stringsAsFactors = FALSE)

OUTDIR <- file.path(
  "outputs", "reviewer_revisions", "reporting_requirements_audit"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------
## 1. Canonical reviewer-facing analyses
## ---------------------------------------------------------

R <- data.frame(
  analysis_id = c(
    "RNA_stage_DE",
    "Protein_stage_DE",
    "Hsp_pathway_stage_Fig1",
    "Protein_Braak_clients",
    "Protein_CERAD_clients",
    "Cognition_network",
    "Cognition_clients",
    "Regional_AD_gene_interaction",
    "Regional_Braak_continuous_gene_interaction",
    "Regional_Braak_categorical_gene_interaction",
    "Regional_Hsp_pathway_interaction",
    "APOE_late_stage_clients",
    "APOE_Braak_clients",
    "APOE_CERAD_clients",
    "Alternative_mechanism_markers",
    "RNA_pathway_PC1_sensitivity",
    "Cognition_matched_null"
  ),
  role = c(
    "Conventional RNA stage differential expression",
    "Conventional protein stage differential abundance",
    "Primary Hsp60/10 pathway stage trajectory",
    "Protein association with continuous Braak",
    "Protein association with CERAD",
    "Network-level cognition association",
    "Client-level cognition association",
    "Formal DLPFC-vs-STG AD interaction",
    "Formal DLPFC-vs-STG continuous Braak interaction",
    "Formal DLPFC-vs-STG categorical Braak interaction",
    "Direct Hsp60/10 pathway regional interaction",
    "APOE sensitivity: late AD-minus-MCI decline",
    "APOE sensitivity: Braak",
    "APOE sensitivity: CERAD",
    "Alternative-mechanism marker panel",
    "Mean-z versus PC1 pathway-score sensitivity",
    "Matched null / specificity analysis"
  ),
  primary = c(
    "outputs/reviewer_revisions/conventional_stage_differential_final/04_RNA_DESeq2_all_contrasts.csv",
    "outputs/reviewer_revisions/conventional_stage_differential_final/07_protein_PRIMARY_complete400_limma.csv",
    "outputs/main_figures/tables/main_fig1_all_clients/Fig1B_all_clients_side_by_side_modality_effects.csv",
    "outputs/main_figures/tables/main_fig3_all_clients/Fig3_all_clients_pathology_table.csv",
    "outputs/main_figures/tables/main_fig3_all_clients/Fig3_all_clients_pathology_table.csv",
    "outputs/reviewer_revisions/corrected_cognition_production_validation/corrected_network_cognition_models.csv",
    "outputs/main_figures/tables/main_fig4_cognition/Fig4_client_level_cognition_all_models.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/validated_final/FINAL_region_x_AD_interaction_all_proteins.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/validated_final/FINAL_region_x_Braak_continuous_interaction_all_proteins.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/validated_final/FINAL_region_x_Braak_bin3_interaction_all_proteins.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/hsp_conclusion_audit/Hsp_pathway_score_primary_interactions.csv",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_sensitivity_late_stage_clients.csv",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_sensitivity_Braak_clients.csv",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_sensitivity_CERAD_clients.csv",
    "outputs/reviewer_revisions/alternative_mechanism_marker_panel/marker_stage_limma_results.csv",
    "outputs/reviewer_revisions/Hsp_pathway_PC1_sensitivity/stage_Wilcoxon_mean_z_vs_PC1.csv",
    "outputs/reviewer_revisions/corrected_cognition_production_validation/corrected_figure5_null_summary.csv"
  ),
  support = c(
    "outputs/reviewer_revisions/conventional_stage_differential_final/20_sample_stage_audit.csv",
    "outputs/reviewer_revisions/conventional_stage_differential_final/20_sample_stage_audit.csv",
    "outputs/main_figures/tables/clinical_stage_sample_counts_before_abundance_missingness.csv",
    "outputs/main_figures/tables/main_fig3_all_clients/Fig3_all_clients_stats_combined.csv",
    "outputs/main_figures/tables/main_fig3_all_clients/Fig3_all_clients_stats_combined.csv",
    "outputs/main_figures/tables/main_fig4_cognition/Fig4_network_level_cognition_models_source_data.csv",
    "outputs/reviewer_revisions/corrected_cognition_production_validation/corrected_client_cognition_summary.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/paired_participant_analysis_manifest.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/paired_participant_analysis_manifest.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/paired_participant_analysis_manifest.csv",
    "outputs/reviewer_revisions/regional_formal_interaction/hsp_conclusion_audit/Hsp_pathway_score_participant_region_table.csv",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_sensitivity_effect_summary.csv",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_sensitivity_effect_summary.csv",
    "outputs/reviewer_revisions/APOE_sensitivity/APOE_sensitivity_effect_summary.csv",
    "outputs/reviewer_revisions/alternative_mechanism_marker_panel/marker_stage_observation_counts.csv",
    "outputs/reviewer_revisions/Hsp_pathway_PC1_sensitivity/stage_score_summary.csv",
    "outputs/main_figures/tables/main_fig5_matched_null_specificity/Fig5_gene_universe_cognition_summary_by_gene.csv"
  ),
  distribution = c(
    NA, NA,
    "outputs/main_figures/tables/main_fig1_all_clients/Fig1A_all_clients_overlay_long.csv",
    NA, NA,
    "outputs/main_figures/tables/main_fig4_cognition/Fig4A_network_diagnosis_panel_dataframe.csv",
    NA, NA, NA, NA,
    "outputs/reviewer_revisions/regional_formal_interaction/hsp_conclusion_audit/Hsp_pathway_score_participant_region_table.csv",
    NA, NA, NA, NA,
    "outputs/reviewer_revisions/Hsp_pathway_PC1_sensitivity/participant_mean_z_and_PC1_scores.csv",
    NA
  ),
  stringsAsFactors = FALSE
)

## required / not_applicable / policy_review
R$effect_req <- "required"
R$se_req <- c(
  "required","required","not_applicable","required","required",
  "required","required","required","required","required","required",
  "required","required","required","required",
  "not_applicable","not_applicable"
)
R$ci_req <- c(
  "required","required","not_applicable","required","required",
  "required","required","required","required","required","required",
  "required","required","required","required",
  "not_applicable","not_applicable"
)
R$n_req <- "required"
R$p_req <- c(
  "required","required","not_applicable","required","required",
  "required","required","required","required","required","required",
  "required","required","required","required",
  "required","required"
)
R$fdr_req <- c(
  "required","required","not_applicable","required","required",
  "policy_review","required","required","required","required",
  "policy_review","required","required","required","required",
  "policy_review","policy_review"
)
R$distribution_req <- c(
  "not_applicable","not_applicable","required","not_applicable",
  "not_applicable","required","not_applicable","not_applicable",
  "not_applicable","not_applicable","required","not_applicable",
  "not_applicable","not_applicable","not_applicable","required",
  "not_applicable"
)

## ---------------------------------------------------------
## 2. Helpers
## ---------------------------------------------------------

split_paths <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

header <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(character())
  }
  tryCatch(
    names(read.csv(
      path, nrows = 2, check.names = FALSE,
      stringsAsFactors = FALSE
    )),
    error = function(e) character()
  )
}

norm <- function(x) tolower(gsub("[^a-z0-9]+", "_", x))

detect <- function(cols, metric) {
  if (!length(cols)) return(character())
  x <- norm(cols)

  pat <- switch(
    metric,
    effect = paste(c(
      "^effect$","^estimate$","^beta$","^coefficient$",
      "^log2foldchange$","^logfc$","^delta$",
      "_effect$","^effect_","_beta$","^beta_",
      "estimate","log2_fold","log2fold","log_fc",
      "correlation","^rho$","region_x_predictor","^stage_change$"
    ), collapse = "|"),
    se = paste(c(
      "^se$","^std_error$","^std_err$","^stderr$","^lfcse$",
      "std_error","std_err","stderr","standard_error",
      "_se$","^se_"
    ), collapse = "|"),
    ci_low = "ci_low|ci_lower|conf_low|conf_lower|lower_95|lwr",
    ci_high = "ci_high|ci_upper|conf_high|conf_upper|upper_95|upr",
    n = paste(c(
      "^n$","^n_","_n$","n_complete","n_model","n_used",
      "n_obs","n_participant","sample_count","participant_count",
      "n_total","n_analy"
    ), collapse = "|"),
    p = paste(c(
      "^p$","^p_value$","^pvalue$","^p_val$",
      "^pvalue_","_pvalue$","^p_","_p$",
      "p_value","pvalue","p_val","p_region_x"
    ), collapse = "|"),
    fdr = paste(c(
      "^fdr$","^padj$","adj_p","adjusted_p",
      "q_value","qvalue","fdr"
    ), collapse = "|")
  )

  cols[grepl(pat, x, perl = TRUE)]
}

resolve <- function(primary, support, metric) {
  paths <- c(primary, split_paths(support))

  for (j in seq_along(paths)) {
    cols <- header(paths[j])

    if (metric == "ci") {
      lo <- detect(cols, "ci_low")
      hi <- detect(cols, "ci_high")
      hit <- length(lo) > 0 && length(hi) > 0
      used <- c(lo, hi)
    } else {
      used <- detect(cols, metric)
      hit <- length(used) > 0
    }

    if (hit) {
      return(list(
        yes = TRUE,
        location = ifelse(j == 1, "primary", "support"),
        file = paths[j],
        columns = paste(used, collapse = ";")
      ))
    }
  }

  list(
    yes = FALSE,
    location = "missing",
    file = NA_character_,
    columns = NA_character_
  )
}

status <- function(req, yes) {
  if (req == "not_applicable") return("NOT_APPLICABLE")
  if (req == "policy_review") {
    return(ifelse(yes, "POLICY_REVIEW_PRESENT", "POLICY_REVIEW_ABSENT"))
  }
  ifelse(yes, "PRESENT", "GAP")
}

## ---------------------------------------------------------
## 3. Build master audit
## ---------------------------------------------------------

metrics <- c("effect","se","ci","n","p","fdr")
rows <- vector("list", nrow(R))

for (i in seq_len(nrow(R))) {
  z <- lapply(metrics, function(m) {
    resolve(R$primary[i], R$support[i], m)
  })
  names(z) <- metrics

  dist_yes <- (
    !is.na(R$distribution[i]) &&
    nzchar(R$distribution[i]) &&
    file.exists(R$distribution[i])
  )

  row <- R[i, c("analysis_id","role","primary","support","distribution")]
  row$primary_exists <- file.exists(R$primary[i])

  for (m in metrics) {
    req <- R[[paste0(m, "_req")]][i]
    row[[paste0(m, "_requirement")]] <- req
    row[[paste0(m, "_status")]] <- status(req, z[[m]]$yes)
    row[[paste0(m, "_location")]] <- z[[m]]$location
    row[[paste0(m, "_file")]] <- z[[m]]$file
    row[[paste0(m, "_columns")]] <- z[[m]]$columns
  }

  row$distribution_requirement <- R$distribution_req[i]
  row$distribution_status <- status(
    R$distribution_req[i], dist_yes
  )

  rows[[i]] <- row
}

audit <- do.call(rbind, rows)

## ---------------------------------------------------------
## 4. File/schema inventory
## ---------------------------------------------------------

paths <- unique(unlist(lapply(seq_len(nrow(R)), function(i) {
  c(R$primary[i], split_paths(R$support[i]), R$distribution[i])
})))
paths <- paths[!is.na(paths) & nzchar(paths)]

schema <- do.call(rbind, lapply(paths, function(f) {
  cols <- header(f)
  data.frame(
    file = f,
    exists = file.exists(f),
    n_columns = length(cols),
    columns = paste(cols, collapse = ";"),
    effect_columns = paste(detect(cols, "effect"), collapse = ";"),
    se_columns = paste(detect(cols, "se"), collapse = ";"),
    ci_low_columns = paste(detect(cols, "ci_low"), collapse = ";"),
    ci_high_columns = paste(detect(cols, "ci_high"), collapse = ";"),
    n_columns_detected = paste(detect(cols, "n"), collapse = ";"),
    p_columns = paste(detect(cols, "p"), collapse = ";"),
    fdr_columns = paste(detect(cols, "fdr"), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

## ---------------------------------------------------------
## 5. Gap and policy tables
## ---------------------------------------------------------

all_metrics <- c(metrics, "distribution")
gaps <- list()
policy <- list()

for (i in seq_len(nrow(audit))) {
  for (m in all_metrics) {
    s <- audit[[paste0(m, "_status")]][i]

    if (s == "GAP") {
      gaps[[length(gaps) + 1]] <- data.frame(
        analysis_id = audit$analysis_id[i],
        role = audit$role[i],
        metric = m,
        primary_file = audit$primary[i],
        stringsAsFactors = FALSE
      )
    }

    if (grepl("^POLICY_REVIEW", s)) {
      policy[[length(policy) + 1]] <- data.frame(
        analysis_id = audit$analysis_id[i],
        role = audit$role[i],
        metric = m,
        status = s,
        stringsAsFactors = FALSE
      )
    }
  }
}

empty_gaps <- data.frame(
  analysis_id=character(), role=character(),
  metric=character(), primary_file=character()
)
empty_policy <- data.frame(
  analysis_id=character(), role=character(),
  metric=character(), status=character()
)

gaps <- if (length(gaps)) do.call(rbind, gaps) else empty_gaps
policy <- if (length(policy)) do.call(rbind, policy) else empty_policy

missing_files <- schema[!schema$exists, c("file","exists"), drop=FALSE]

analysis_state <- data.frame(
  analysis_id = audit$analysis_id,
  primary_exists = audit$primary_exists,
  n_hard_gaps = vapply(seq_len(nrow(audit)), function(i) {
    sum(vapply(all_metrics, function(m) {
      audit[[paste0(m, "_status")]][i] == "GAP"
    }, logical(1)))
  }, integer(1)),
  n_policy_review = vapply(seq_len(nrow(audit)), function(i) {
    sum(vapply(all_metrics, function(m) {
      grepl("^POLICY_REVIEW", audit[[paste0(m, "_status")]][i])
    }, logical(1)))
  }, integer(1)),
  stringsAsFactors = FALSE
)

analysis_state$state <- ifelse(
  !analysis_state$primary_exists,
  "SOURCE_FILE_MISSING",
  ifelse(
    analysis_state$n_hard_gaps > 0,
    "REQUIRES_REPORTING_PATCH",
    ifelse(
      analysis_state$n_policy_review > 0,
      "POLICY_DECISION_PENDING",
      "REPORTING_FIELDS_PRESENT"
    )
  )
)

metric_summary <- do.call(rbind, lapply(all_metrics, function(m) {
  s <- audit[[paste0(m, "_status")]]
  data.frame(
    metric = m,
    present = sum(s == "PRESENT"),
    hard_gaps = sum(s == "GAP"),
    policy_present = sum(s == "POLICY_REVIEW_PRESENT"),
    policy_absent = sum(s == "POLICY_REVIEW_ABSENT"),
    not_applicable = sum(s == "NOT_APPLICABLE")
  )
}))

## ---------------------------------------------------------
## 6. Statistical-policy notes
## ---------------------------------------------------------

notes <- data.frame(
  topic = c(
    "Conventional stage DE",
    "Regional gene-level testing",
    "Pathology client testing",
    "Cognition client testing",
    "Network/pathway multiplicity",
    "Distributions",
    "Cross-modal subtraction"
  ),
  policy = c(
    "Benjamini-Hochberg FDR within each tested contrast/universe.",
    "Genome-wide BH FDR retained for Hsp subsets; no Hsp-only re-FDR.",
    "BH adjustment retained from production pathology outputs.",
    "Client-level multiple-testing adjustment should remain explicit.",
    "Flagged for explicit reviewer/manuscript policy decision; not silently assigned.",
    "Participant-level source data audited separately from model SE/CI.",
    "Conventional cross-modal DE is side-by-side; no RNA-minus-protein/protein-minus-RNA subtraction."
  )
)

## ---------------------------------------------------------
## 7. Internal validation of the audit itself
## ---------------------------------------------------------

stopifnot(
  nrow(R) == 17,
  anyDuplicated(R$analysis_id) == 0,
  nrow(audit) == 17,
  nrow(analysis_state) == 17,
  nrow(metric_summary) == 7,
  nrow(gaps) == sum(vapply(all_metrics, function(m) {
    sum(audit[[paste0(m, "_status")]] == "GAP")
  }, integer(1))),
  nrow(policy) == sum(vapply(all_metrics, function(m) {
    sum(grepl("^POLICY_REVIEW", audit[[paste0(m, "_status")]]))
  }, integer(1))),
  nrow(missing_files) == sum(!schema$exists)
)

## ---------------------------------------------------------
## 8. Write
## ---------------------------------------------------------

write.csv(R, file.path(OUTDIR, "01_registry.csv"), row.names=FALSE)
write.csv(schema, file.path(OUTDIR, "02_schema_inventory.csv"), row.names=FALSE)
write.csv(audit, file.path(OUTDIR, "03_master_reporting_audit.csv"), row.names=FALSE)
write.csv(gaps, file.path(OUTDIR, "04_hard_reporting_gaps.csv"), row.names=FALSE)
write.csv(policy, file.path(OUTDIR, "05_policy_review_items.csv"), row.names=FALSE)
write.csv(metric_summary, file.path(OUTDIR, "06_metric_summary.csv"), row.names=FALSE)
write.csv(analysis_state, file.path(OUTDIR, "07_analysis_state.csv"), row.names=FALSE)
write.csv(missing_files, file.path(OUTDIR, "08_missing_source_files.csv"), row.names=FALSE)
write.csv(notes, file.path(OUTDIR, "09_statistical_policy_notes.csv"), row.names=FALSE)

git_head <- tryCatch(
  system2("git", c("rev-parse","HEAD"), stdout=TRUE, stderr=FALSE),
  error=function(e) NA_character_
)

prov <- data.frame(
  item=c("analysis_n","hard_gap_n","policy_review_n",
         "missing_source_n","git_HEAD"),
  value=c(nrow(R),nrow(gaps),nrow(policy),nrow(missing_files),
          paste(git_head, collapse=";"))
)
write.csv(prov, file.path(OUTDIR, "10_provenance.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "11_sessionInfo.txt"))

## ---------------------------------------------------------
## 9. Console report
## ---------------------------------------------------------

cat("\n============================================================\n")
cat("MASTER REVIEWER REPORTING REQUIREMENTS AUDIT\n")
cat("============================================================\n\n")

cat("Analyses registered: ", nrow(R), "\n", sep="")
cat("Missing registered source/support files: ",
    nrow(missing_files), "\n", sep="")
cat("Hard reporting gaps: ", nrow(gaps), "\n", sep="")
cat("Policy-review items: ", nrow(policy), "\n\n", sep="")

cat("METRIC SUMMARY\n")
print(metric_summary, row.names=FALSE)

cat("\nANALYSIS STATE\n")
print(analysis_state, row.names=FALSE)

cat("\nHARD GAPS\n")
if (nrow(gaps)) {
  print(gaps[,c("analysis_id","metric")], row.names=FALSE)
} else {
  cat("None detected.\n")
}

cat("\nPOLICY REVIEW\n")
if (nrow(policy)) {
  print(policy[,c("analysis_id","metric","status")], row.names=FALSE)
} else {
  cat("None.\n")
}

cat("\nOutputs: ", OUTDIR, "\n", sep="")

cat("\n============================================================\n")
cat("REPORTING REQUIREMENTS AUDIT COMPLETE\n")
cat("INTERNAL AUDIT VALIDATION PASSED\n")
cat("============================================================\n")
