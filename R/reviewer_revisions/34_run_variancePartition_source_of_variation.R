options(stringsAsFactors = FALSE)

# =============================================================================
# Reviewer 2: source-of-variation / covariate-selection analysis
#
# PRIMARY MODELS
#
# RNA:
#   expression ~ age + sex + PMI + RIN + (1 | RNA_batch)
#
# Protein:
#   abundance ~ age + sex + PMI + (1 | TMT_batch)
#
# IMPORTANT:
#   - non-residualized normalized expression/abundance only
#   - one nonredundant batch term per modality
#   - disease/pathology variables are NOT nuisance covariates in this analysis
# =============================================================================


# =============================================================================
# Packages
# =============================================================================

required_packages <- c(
  "variancePartition",
  "BiocParallel"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}


# =============================================================================
# Production objects
# =============================================================================

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {

  if (!file.exists(f)) {
    stop("Missing production script: ", f)
  }

  source(f)
}


# =============================================================================
# Output
# =============================================================================

OUTDIR <- "outputs/reviewer_revisions/variancePartition_source_of_variation"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# Helpers
# =============================================================================

zscore <- function(x) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  s <- sd(
    x,
    na.rm = TRUE
  )

  if (!is.finite(s) || s == 0) {
    stop("Cannot z-score zero-variance/nonfinite variable.")
  }

  as.numeric(scale(x))
}


encode_sex <- function(x) {

  y <- tolower(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_real_,
    length(y)
  )

  out[y %in% c("female", "f", "0")] <- 0
  out[y %in% c("male", "m", "1")] <- 1

  bad <- unique(
    y[
      !is.na(y) &
        is.na(out)
    ]
  )

  if (length(bad) > 0) {
    stop(
      "Unexpected sex values: ",
      paste(bad, collapse = ", ")
    )
  }

  out
}


summarize_varpart <- function(
  vp,
  modality,
  universe
) {

  x <- as.data.frame(vp)

  rows <- lapply(
    names(x),
    function(term) {

      v <- x[[term]]

      data.frame(
        modality = modality,
        universe = universe,
        source = term,

        n_features =
          sum(is.finite(v)),

        mean_fraction =
          mean(
            v,
            na.rm = TRUE
          ),

        median_fraction =
          median(
            v,
            na.rm = TRUE
          ),

        q25_fraction =
          unname(
            quantile(
              v,
              0.25,
              na.rm = TRUE
            )
          ),

        q75_fraction =
          unname(
            quantile(
              v,
              0.75,
              na.rm = TRUE
            )
          ),

        q90_fraction =
          unname(
            quantile(
              v,
              0.90,
              na.rm = TRUE
            )
          ),

        fraction_features_gt_0.01 =
          mean(
            v > 0.01,
            na.rm = TRUE
          ),

        fraction_features_gt_0.05 =
          mean(
            v > 0.05,
            na.rm = TRUE
          ),

        fraction_features_gt_0.10 =
          mean(
            v > 0.10,
            na.rm = TRUE
          ),

        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    rows
  )
}


top_source_table <- function(
  vp,
  modality
) {

  x <- as.data.frame(vp)

  nonresid <- setdiff(
    names(x),
    grep(
      "Residual",
      names(x),
      value = TRUE,
      ignore.case = TRUE
    )
  )

  if (length(nonresid) == 0) {
    stop(
      "No non-residual variance components found for ",
      modality
    )
  }

  m <- as.matrix(
    x[, nonresid, drop = FALSE]
  )

  winner <- apply(
    m,
    1,
    function(v) {

      if (all(!is.finite(v))) {
        return(NA_character_)
      }

      names(v)[
        which.max(v)
      ]
    }
  )

  tab <- sort(
    table(winner),
    decreasing = TRUE
  )

  data.frame(
    modality = modality,
    source = names(tab),
    n_features_top_source =
      as.integer(tab),
    fraction_features_top_source =
      as.integer(tab) /
      sum(tab),
    stringsAsFactors = FALSE
  )
}


extract_error_attribute <- function(
  vp,
  modality
) {

  possible <- c(
    "errors",
    "error",
    "error.initial"
  )

  out <- list()

  k <- 1

  for (a in possible) {

    value <- attr(
      vp,
      a
    )

    if (is.null(value)) {
      next
    }

    txt <- capture.output(
      print(value)
    )

    out[[k]] <- data.frame(
      modality = modality,
      attribute = a,
      contents = paste(
        txt,
        collapse = "\n"
      ),
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }

  if (length(out) == 0) {

    return(
      data.frame(
        modality = character(),
        attribute = character(),
        contents = character()
      )
    )
  }

  do.call(
    rbind,
    out
  )
}


# =============================================================================
# Hard object validation
# =============================================================================

required_objects <- c(
  "rna_mat_raw",
  "rna_mat",
  "rna_mat_adj",
  "prot_mat_raw",
  "prot_mat",
  "prot_mat_adj",
  "rna_meta",
  "rna_meta_adj",
  "prot_meta_adj"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing expected production objects: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}

if (
  !identical(
    colnames(rna_mat_raw),
    as.character(
      rna_meta_adj$sample_id
    )
  )
) {
  stop(
    "RNA matrix columns are not exactly aligned to rna_meta_adj."
  )
}

if (
  !identical(
    colnames(prot_mat_raw),
    as.character(
      prot_meta_adj$sample_id
    )
  )
) {
  stop(
    "Protein matrix columns are not exactly aligned to prot_meta_adj."
  )
}

if (
  !identical(
    as.character(
      rna_meta$library_batch
    ),
    as.character(
      rna_meta$sequencing_batch
    )
  )
) {
  stop(
    "RNA library_batch and sequencing_batch are no longer identical."
  )
}

if (
  !identical(
    as.character(
      prot_meta_adj$tmt_batch
    ),
    as.character(
      prot_meta_adj$batch
    )
  )
) {
  stop(
    "Protein tmt_batch and batch are no longer identical."
  )
}

if (isTRUE(all.equal(rna_mat_raw, rna_mat_adj))) {
  stop(
    "RNA raw and adjusted matrices unexpectedly identical."
  )
}

if (isTRUE(all.equal(prot_mat_raw, prot_mat_adj))) {
  stop(
    "Protein raw and adjusted matrices unexpectedly identical."
  )
}


# =============================================================================
# RNA design
# =============================================================================

rna_design <- data.frame(
  sample_id =
    as.character(
      rna_meta_adj$sample_id
    ),

  age_z =
    zscore(
      rna_meta_adj$age_num
    ),

  sex_binary =
    encode_sex(
      rna_meta_adj$sex_label
    ),

  pmi_z =
    zscore(
      rna_meta_adj$pmi_num
    ),

  rin_z =
    zscore(
      rna_meta_adj$rin_num
    ),

  sequencing_batch =
    factor(
      as.character(
        rna_meta_adj$sequencing_batch
      )
    ),

  stringsAsFactors = FALSE
)

rna_complete <- complete.cases(
  rna_design
)

if (sum(rna_complete) != 577) {
  stop(
    "Expected 577 complete RNA samples; found ",
    sum(rna_complete)
  )
}

rna_design <- rna_design[
  rna_complete,
  ,
  drop = FALSE
]

if (anyDuplicated(rna_design$sample_id)) {
  stop("Duplicate RNA sample IDs in variancePartition design.")
}

rownames(rna_design) <- rna_design$sample_id

if (nlevels(rna_design$sequencing_batch) != 9) {
  stop(
    "Expected 9 canonical RNA sequencing-batch levels; found ",
    nlevels(rna_design$sequencing_batch)
  )
}

RNA <- rna_mat_raw[
  ,
  rna_design$sample_id,
  drop = FALSE
]

if (
  nrow(RNA) != 21433 ||
  ncol(RNA) != 577
) {
  stop(
    "Unexpected RNA dimensions: ",
    nrow(RNA),
    " x ",
    ncol(RNA)
  )
}

if (anyNA(RNA)) {
  stop(
    "RNA non-residualized matrix unexpectedly contains missing values."
  )
}


# =============================================================================
# Protein design
# =============================================================================

prot_design_all <- data.frame(
  sample_id =
    as.character(
      prot_meta_adj$sample_id
    ),

  age_z =
    zscore(
      prot_meta_adj$age_num
    ),

  sex_binary =
    encode_sex(
      prot_meta_adj$sex_label
    ),

  pmi_z =
    zscore(
      prot_meta_adj$pmi_num
    ),

  tmt_batch =
    factor(
      as.character(
        prot_meta_adj$tmt_batch
      )
    ),

  stringsAsFactors = FALSE
)

prot_complete <- complete.cases(
  prot_design_all
)

if (sum(prot_complete) != 398) {
  stop(
    "Expected 398 complete protein nuisance-covariate samples; found ",
    sum(prot_complete)
  )
}

prot_design <- prot_design_all[
  prot_complete,
  ,
  drop = FALSE
]

if (anyDuplicated(prot_design$sample_id)) {
  stop("Duplicate protein sample IDs in variancePartition design.")
}

rownames(prot_design) <- prot_design$sample_id

PROT_all <- prot_mat_raw[
  ,
  prot_design$sample_id,
  drop = FALSE
]

if (ncol(PROT_all) != 398) {
  stop(
    "Expected 398 protein samples after nuisance-covariate filtering."
  )
}

protein_n_observed <- rowSums(
  !is.na(PROT_all)
)

protein_fraction_observed <- (
  protein_n_observed /
  ncol(PROT_all)
)

# Descriptive missingness audit:
# retain the >=80% flag, but do not model incomplete protein rows.
protein_keep_80 <- (
  protein_fraction_observed >= 0.80
)

# Primary variancePartition protein universe:
# complete across all 398 nuisance-complete participants.
protein_keep_complete <- (
  protein_n_observed ==
    ncol(PROT_all)
)

PROT <- PROT_all[
  protein_keep_complete,
  ,
  drop = FALSE
]

if (nrow(PROT) != 5088) {
  stop(
    "Expected 5088 proteins complete across all 398 samples; found ",
    nrow(PROT)
  )
}

if (anyNA(PROT)) {
  stop("Primary protein variancePartition matrix contains missing values.")
}


# =============================================================================
# Feature manifests
# =============================================================================

rna_manifest <- data.frame(
  feature_id = rownames(RNA),
  n_observed = ncol(RNA),
  fraction_observed = 1,
  included_primary = TRUE,
  stringsAsFactors = FALSE
)

protein_manifest <- data.frame(
  feature_id = rownames(PROT_all),
  n_observed = protein_n_observed,
  fraction_observed = protein_fraction_observed,
  observed_ge_80pct = protein_keep_80,
  complete_all_398 = protein_keep_complete,
  included_primary = protein_keep_complete,
  stringsAsFactors = FALSE
)

write.csv(
  rna_manifest,
  file.path(
    OUTDIR,
    "RNA_variancePartition_feature_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  protein_manifest,
  file.path(
    OUTDIR,
    "protein_variancePartition_feature_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  rna_design,
  file.path(
    OUTDIR,
    "RNA_variancePartition_sample_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  prot_design,
  file.path(
    OUTDIR,
    "protein_variancePartition_sample_manifest.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Report pre-model dimensions
# =============================================================================

cat("======================================================================\n")
cat("VARIANCEPARTITION SOURCE-OF-VARIATION ANALYSIS\n")
cat("======================================================================\n")

cat("\nRNA:\n")
cat("  features:", nrow(RNA), "\n")
cat("  samples:", ncol(RNA), "\n")
cat("  batch levels:", nlevels(rna_design$sequencing_batch), "\n")

cat("\nProtein:\n")
cat("  original features:", nrow(PROT_all), "\n")
cat("  nuisance-complete samples:", ncol(PROT_all), "\n")
cat(
  "  >=80% observed features:",
  sum(protein_keep_80),
  "\n"
)
cat(
  "  complete across all 398:",
  sum(protein_keep_complete),
  "\n"
)
cat(
  "  TMT batch levels:",
  nlevels(prot_design$tmt_batch),
  "\n"
)


# =============================================================================
# Formulas
# =============================================================================

rna_formula <- ~
  age_z +
  sex_binary +
  pmi_z +
  rin_z +
  (1 | sequencing_batch)

protein_formula <- ~
  age_z +
  sex_binary +
  pmi_z +
  (1 | tmt_batch)

cat("\nRNA model:\n")
print(rna_formula)

cat("\nProtein model:\n")
print(protein_formula)


# =============================================================================
# Collinearity audit
# =============================================================================

rna_cor <- variancePartition::canCorPairs(
  rna_formula,
  rna_design
)

protein_cor <- variancePartition::canCorPairs(
  protein_formula,
  prot_design
)

write.csv(
  as.data.frame(rna_cor),
  file.path(
    OUTDIR,
    "RNA_canCorPairs.csv"
  ),
  row.names = TRUE
)

write.csv(
  as.data.frame(protein_cor),
  file.path(
    OUTDIR,
    "protein_canCorPairs.csv"
  ),
  row.names = TRUE
)

cat("\nRNA canonical-correlation structure:\n")
print(rna_cor)

cat("\nProtein canonical-correlation structure:\n")
print(protein_cor)


# =============================================================================
# Parallel backend
# =============================================================================

workers <- min(
  4L,
  parallel::detectCores()
)

param <- BiocParallel::SnowParam(
  workers = workers,
  type = "SOCK",
  progressbar = TRUE
)

cat(
  "\nUsing",
  workers,
  "parallel workers.\n"
)


# =============================================================================
# RNA variance partition
# =============================================================================

cat("\n======================================================================\n")
cat("RUNNING RNA variancePartition\n")
cat("======================================================================\n")

rna_vp <- variancePartition::fitExtractVarPartModel(
  exprObj = RNA,
  formula = rna_formula,
  data = rna_design,
  REML = FALSE,
  useWeights = FALSE,
  showWarnings = TRUE,
  BPPARAM = param
)

rna_vp_df <- as.data.frame(
  rna_vp
)

rna_vp_df$feature_id <- rownames(
  rna_vp_df
)

rna_vp_df <- rna_vp_df[
  ,
  c(
    "feature_id",
    setdiff(
      names(rna_vp_df),
      "feature_id"
    )
  ),
  drop = FALSE
]

write.csv(
  rna_vp_df,
  file.path(
    OUTDIR,
    "RNA_variance_fractions_all_features.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Protein variance partition
# =============================================================================

cat("\n======================================================================\n")
cat("RUNNING PROTEIN variancePartition\n")
cat("======================================================================\n")

protein_vp <- variancePartition::fitExtractVarPartModel(
  exprObj = PROT,
  formula = protein_formula,
  data = prot_design,
  REML = FALSE,
  useWeights = FALSE,
  showWarnings = TRUE,
  BPPARAM = param
)

protein_vp_df <- as.data.frame(
  protein_vp
)

protein_vp_df$feature_id <- rownames(
  protein_vp_df
)

protein_vp_df <- protein_vp_df[
  ,
  c(
    "feature_id",
    setdiff(
      names(protein_vp_df),
      "feature_id"
    )
  ),
  drop = FALSE
]

write.csv(
  protein_vp_df,
  file.path(
    OUTDIR,
    "protein_variance_fractions_complete_all_398.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Fraction-sum validation
# =============================================================================

rna_components <- setdiff(
  names(rna_vp_df),
  "feature_id"
)

protein_components <- setdiff(
  names(protein_vp_df),
  "feature_id"
)

rna_sums <- rowSums(
  rna_vp_df[
    ,
    rna_components,
    drop = FALSE
  ],
  na.rm = FALSE
)

protein_sums <- rowSums(
  protein_vp_df[
    ,
    protein_components,
    drop = FALSE
  ],
  na.rm = FALSE
)

cat("\nVariance-fraction sum audit:\n")

cat(
  "  RNA max |sum - 1|:",
  max(
    abs(
      rna_sums - 1
    ),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "  Protein max |sum - 1|:",
  max(
    abs(
      protein_sums - 1
    ),
    na.rm = TRUE
  ),
  "\n"
)


# =============================================================================
# Summaries
# =============================================================================

rna_summary <- summarize_varpart(
  rna_vp,
  "RNA",
  "all_21433_features"
)

protein_summary <- summarize_varpart(
  protein_vp,
  "protein",
  "complete_all_398"
)

master_summary <- rbind(
  rna_summary,
  protein_summary
)

write.csv(
  master_summary,
  file.path(
    OUTDIR,
    "MASTER_variance_fraction_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Which modeled source is largest for each feature?
# =============================================================================

top_sources <- rbind(
  top_source_table(
    rna_vp,
    "RNA"
  ),
  top_source_table(
    protein_vp,
    "protein"
  )
)

write.csv(
  top_sources,
  file.path(
    OUTDIR,
    "top_nonresidual_source_counts.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Model error/warning audit
# =============================================================================

model_attributes <- rbind(
  extract_error_attribute(
    rna_vp,
    "RNA"
  ),
  extract_error_attribute(
    protein_vp,
    "protein"
  )
)

write.csv(
  model_attributes,
  file.path(
    OUTDIR,
    "variancePartition_model_error_attributes.csv"
  ),
  row.names = FALSE
)

rna_bad <- sum(
  !complete.cases(
    rna_vp_df[
      ,
      rna_components,
      drop = FALSE
    ]
  )
)

protein_bad <- sum(
  !complete.cases(
    protein_vp_df[
      ,
      protein_components,
      drop = FALSE
    ]
  )
)

cat("\nModel-result completeness:\n")

cat(
  "  RNA feature rows with any missing variance fraction:",
  rna_bad,
  "\n"
)

cat(
  "  Protein feature rows with any missing variance fraction:",
  protein_bad,
  "\n"
)


# =============================================================================
# Console summary
# =============================================================================

cat("\n======================================================================\n")
cat("VARIANCE FRACTION SUMMARY\n")
cat("======================================================================\n")

print(
  master_summary[
    ,
    c(
      "modality",
      "universe",
      "source",
      "n_features",
      "median_fraction",
      "q25_fraction",
      "q75_fraction",
      "q90_fraction",
      "fraction_features_gt_0.05"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("TOP NON-RESIDUAL SOURCE PER FEATURE\n")
cat("======================================================================\n")

print(
  top_sources,
  row.names = FALSE
)


# =============================================================================
# Save R objects for later figures / sensitivity analyses
# =============================================================================

saveRDS(
  rna_vp,
  file.path(
    OUTDIR,
    "RNA_variancePartition_results.rds"
  )
)

saveRDS(
  protein_vp,
  file.path(
    OUTDIR,
    "protein_variancePartition_results.rds"
  )
)


# =============================================================================
# Provenance
# =============================================================================

provenance <- data.frame(
  item = c(
    "RNA_expression_object",
    "RNA_features",
    "RNA_samples",
    "RNA_batch_variable",
    "RNA_batch_levels",

    "protein_expression_object",
    "protein_original_features",
    "protein_primary_features",
    "protein_samples",
    "protein_batch_variable",
    "protein_batch_levels",

    "REML",
    "protein_primary_universe",
    "parallel_workers"
  ),

  value = c(
    "rna_mat_raw",
    nrow(RNA),
    ncol(RNA),
    "sequencing_batch",
    nlevels(rna_design$sequencing_batch),

    "prot_mat_raw",
    nrow(PROT_all),
    nrow(PROT),
    ncol(PROT),
    "tmt_batch",
    nlevels(prot_design$tmt_batch),

    "FALSE",
    "complete_all_398",
    workers
  ),

  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(
    OUTDIR,
    "analysis_provenance.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)


# =============================================================================
# Hard validation
# =============================================================================

if (nrow(rna_vp_df) != nrow(RNA)) {
  stop(
    "RNA variancePartition returned ",
    nrow(rna_vp_df),
    " rows for ",
    nrow(RNA),
    " input features."
  )
}

if (!setequal(rna_vp_df$feature_id, rownames(RNA))) {
  stop("RNA variancePartition feature IDs do not match the input universe.")
}

if (nrow(protein_vp_df) != 5088) {
  stop(
    "Protein variancePartition returned ",
    nrow(protein_vp_df),
    " rows; expected exactly 5088."
  )
}

if (!setequal(protein_vp_df$feature_id, rownames(PROT))) {
  stop(
    paste0(
      "Protein variancePartition feature IDs do not match ",
      "the complete-all-398 input universe."
    )
  )
}

if (rna_bad > 0) {
  stop(
    "RNA variancePartition contains incomplete result rows: ",
    rna_bad
  )
}

if (protein_bad > 0) {
  stop(
    "Protein variancePartition contains incomplete result rows: ",
    protein_bad
  )
}

if (
  max(
    abs(
      rna_sums - 1
    ),
    na.rm = TRUE
  ) > 1e-5
) {
  stop(
    "RNA variance fractions do not sum to ~1."
  )
}

if (
  max(
    abs(
      protein_sums - 1
    ),
    na.rm = TRUE
  ) > 1e-5
) {
  stop(
    "Protein variance fractions do not sum to ~1."
  )
}

cat("\n======================================================================\n")
cat("SOURCE-OF-VARIATION ANALYSIS COMPLETED AND VALIDATED\n")
cat("======================================================================\n")

cat(
  "Outputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)
