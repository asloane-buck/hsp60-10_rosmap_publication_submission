options(stringsAsFactors = FALSE)

OUTDIR <- "outputs/reviewer_revisions/RNA_batch_adjustment_sensitivity"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R",
  "R/04_build_pathway_sets_all_clients.R"
)) {

  if (!file.exists(f)) {
    stop("Missing production script: ", f)
  }

  source(f)
}


# =============================================================================
# Helpers
# =============================================================================

residualize_preserve_mean <- function(
  Y,
  X
) {

  # Y = features x samples
  # X = samples x covariates

  if (ncol(Y) != nrow(X)) {
    stop("Y/X dimensions do not align.")
  }

  fit <- lm.fit(
    x = X,
    y = t(Y)
  )

  residuals <- fit$residuals

  # Return to features x samples and preserve each feature's original mean.
  out <- t(residuals)

  out <- sweep(
    out,
    1,
    rowMeans(
      Y,
      na.rm = TRUE
    ),
    "+"
  )

  out
}


contrast_means <- function(
  mat,
  stage
) {

  required <- c(
    "NCI",
    "MCI",
    "AD"
  )

  if (
    !all(
      required %in%
        unique(stage)
    )
  ) {
    stop("Missing one or more canonical stages.")
  }

  means <- lapply(
    required,
    function(s) {

      rowMeans(
        mat[
          ,
          stage == s,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    }
  )

  names(means) <- required

  data.frame(
    feature_id =
      rownames(mat),

    MCI_vs_NCI =
      means$MCI -
      means$NCI,

    AD_vs_MCI =
      means$AD -
      means$MCI,

    AD_vs_NCI =
      means$AD -
      means$NCI,

    stringsAsFactors = FALSE
  )
}


compare_effects <- function(
  old,
  new,
  contrast
) {

  x <- old[[contrast]]
  y <- new[[contrast]]

  data.frame(
    contrast = contrast,

    n_features =
      sum(
        is.finite(x) &
          is.finite(y)
      ),

    spearman =
      cor(
        x,
        y,
        method = "spearman",
        use = "complete.obs"
      ),

    pearson =
      cor(
        x,
        y,
        method = "pearson",
        use = "complete.obs"
      ),

    sign_agreement =
      mean(
        sign(x) ==
          sign(y),
        na.rm = TRUE
      ),

    median_abs_effect_old =
      median(
        abs(x),
        na.rm = TRUE
      ),

    median_abs_effect_batch_adjusted =
      median(
        abs(y),
        na.rm = TRUE
      ),

    median_abs_change =
      median(
        abs(y - x),
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Validate canonical RNA objects
# =============================================================================

cat("======================================================================\n")
cat("RNA BATCH-ADJUSTMENT SENSITIVITY\n")
cat("======================================================================\n")

if (
  nrow(rna_mat_raw) != 21433 ||
  ncol(rna_mat_raw) != 578
) {
  stop(
    "Unexpected rna_mat_raw dimensions."
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
    "RNA matrix and metadata are not exactly aligned."
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
    "library_batch and sequencing_batch are no longer identical."
  )
}


# =============================================================================
# Design
# =============================================================================

meta <- data.frame(
  sample_id =
    as.character(
      rna_meta_adj$sample_id
    ),

  age_num =
    as.numeric(
      rna_meta_adj$age_num
    ),

  sex_factor =
    factor(
      rna_meta_adj$sex_factor
    ),

  pmi_num =
    as.numeric(
      rna_meta_adj$pmi_num
    ),

  rin_num =
    as.numeric(
      rna_meta_adj$rin_num
    ),

  batch_factor =
    factor(
      as.character(
        rna_meta$library_batch
      )
    ),

  clinical_stage =
    factor(
      as.character(
        rna_meta_adj$clinical_stage
      ),
      levels = c(
        "NCI",
        "MCI",
        "AD"
      )
    ),

  stringsAsFactors = FALSE
)

if (!all(complete.cases(meta))) {
  stop(
    "Unexpected missing RNA design values."
  )
}

cat(
  "Samples:",
  nrow(meta),
  "\n"
)

cat(
  "Batch levels:",
  nlevels(meta$batch_factor),
  "\n"
)

cat(
  "Stage counts:\n"
)

print(
  table(
    meta$clinical_stage
  )
)


# =============================================================================
# Existing and batch-augmented model matrices
# =============================================================================

X_current <- model.matrix(
  ~ age_num +
    sex_factor +
    pmi_num +
    rin_num,
  data = meta
)

X_batch <- model.matrix(
  ~ age_num +
    sex_factor +
    pmi_num +
    rin_num +
    batch_factor,
  data = meta
)

cat("\n======================================================================\n")
cat("DESIGN MATRIX AUDIT\n")
cat("======================================================================\n")

cat(
  "Current columns:",
  ncol(X_current),
  "\n"
)

cat(
  "Current rank:",
  qr(X_current)$rank,
  "\n"
)

cat(
  "Batch-adjusted columns:",
  ncol(X_batch),
  "\n"
)

cat(
  "Batch-adjusted rank:",
  qr(X_batch)$rank,
  "\n"
)

if (
  qr(X_current)$rank !=
    ncol(X_current)
) {
  stop(
    "Current nuisance model is rank deficient."
  )
}

if (
  qr(X_batch)$rank !=
    ncol(X_batch)
) {
  stop(
    "Batch-augmented nuisance model is rank deficient."
  )
}

cat(
  "Current condition number:",
  kappa(X_current),
  "\n"
)

cat(
  "Batch-adjusted condition number:",
  kappa(X_batch),
  "\n"
)


# =============================================================================
# Reconstruct current production adjustment
# =============================================================================

cat("\nReconstructing current production adjustment...\n")

rna_current_reconstructed <-
  residualize_preserve_mean(
    rna_mat_raw,
    X_current
  )

if (
  !identical(
    dim(rna_current_reconstructed),
    dim(rna_mat_adj)
  )
) {
  stop(
    "Reconstructed/current RNA dimensions differ."
  )
}

reconstruction_difference <- (
  rna_current_reconstructed -
    rna_mat_adj
)

max_abs_reconstruction_difference <-
  max(
    abs(
      reconstruction_difference
    ),
    na.rm = TRUE
  )

median_abs_reconstruction_difference <-
  median(
    abs(
      reconstruction_difference
    ),
    na.rm = TRUE
  )

cat(
  "\nCurrent-adjustment reconstruction:\n"
)

cat(
  "  Max absolute difference:",
  max_abs_reconstruction_difference,
  "\n"
)

cat(
  "  Median absolute difference:",
  median_abs_reconstruction_difference,
  "\n"
)

# Do not proceed unless we have reproduced production closely.
if (
  max_abs_reconstruction_difference >
    1e-6
) {

  stop(
    "Current production residualization was not reproduced closely enough.\n",
    "Do not interpret the batch sensitivity until residualization method ",
    "is reconciled."
  )
}


# =============================================================================
# Create batch-adjusted sensitivity matrix
# =============================================================================

cat(
  "\nCurrent production adjustment reproduced.\n",
  "Building batch-adjusted sensitivity matrix...\n",
  sep = ""
)

rna_batch_adj <-
  residualize_preserve_mean(
    rna_mat_raw,
    X_batch
  )

rownames(rna_batch_adj) <-
  rownames(rna_mat_raw)

colnames(rna_batch_adj) <-
  colnames(rna_mat_raw)


# =============================================================================
# Incremental batch variance beyond existing nuisance model
# =============================================================================

cat("\n======================================================================\n")
cat("INCREMENTAL BATCH CONTRIBUTION\n")
cat("======================================================================\n")

Y <- t(
  rna_mat_raw
)

fit_current <- lm.fit(
  X_current,
  Y
)

fit_batch <- lm.fit(
  X_batch,
  Y
)

rss_current <- colSums(
  fit_current$residuals^2
)

rss_batch <- colSums(
  fit_batch$residuals^2
)

partial_batch_r2 <- (
  rss_current -
    rss_batch
) /
  rss_current

partial_batch_r2[
  partial_batch_r2 < 0 &
    partial_batch_r2 > -1e-12
] <- 0

incremental <- data.frame(
  feature_id =
    rownames(rna_mat_raw),

  partial_batch_R2 =
    partial_batch_r2,

  stringsAsFactors = FALSE
)

write.csv(
  incremental,
  file.path(
    OUTDIR,
    "RNA_incremental_batch_partial_R2_all_genes.csv"
  ),
  row.names = FALSE
)

cat(
  "Batch partial R2 after age/sex/PMI/RIN:\n"
)

print(
  summary(
    incremental$partial_batch_R2
  )
)

cat(
  "Fraction genes partial R2 > 0.01:",
  mean(
    incremental$partial_batch_R2 > 0.01
  ),
  "\n"
)

cat(
  "Fraction genes partial R2 > 0.05:",
  mean(
    incremental$partial_batch_R2 > 0.05
  ),
  "\n"
)

cat(
  "Fraction genes partial R2 > 0.10:",
  mean(
    incremental$partial_batch_R2 > 0.10
  ),
  "\n"
)


# =============================================================================
# Does the added adjustment actually remove batch structure?
# =============================================================================

batch_only_r2 <- function(
  y,
  batch
) {

  summary(
    lm(
      y ~ batch
    )
  )$r.squared
}

audit_genes <- unique(
  round(
    seq(
      1,
      nrow(rna_mat_raw),
      length.out = min(
        2000,
        nrow(rna_mat_raw)
      )
    )
  )
)

r2_current <- vapply(
  audit_genes,
  function(i) {

    batch_only_r2(
      rna_mat_adj[i, ],
      meta$batch_factor
    )
  },
  numeric(1)
)

r2_batch_adjusted <- vapply(
  audit_genes,
  function(i) {

    batch_only_r2(
      rna_batch_adj[i, ],
      meta$batch_factor
    )
  },
  numeric(1)
)

batch_signal <- data.frame(
  feature_id =
    rownames(rna_mat_raw)[
      audit_genes
    ],

  current_adjusted_batch_R2 =
    r2_current,

  batch_adjusted_batch_R2 =
    r2_batch_adjusted,

  stringsAsFactors = FALSE
)

write.csv(
  batch_signal,
  file.path(
    OUTDIR,
    "RNA_batch_signal_current_vs_batch_adjusted_2000genes.csv"
  ),
  row.names = FALSE
)

cat(
  "\nMedian batch R2, current adjusted:",
  median(
    r2_current
  ),
  "\n"
)

cat(
  "Median batch R2, batch-adjusted:",
  median(
    r2_batch_adjusted
  ),
  "\n"
)


# =============================================================================
# Genome-wide stage-effect robustness
# =============================================================================

stage <- as.character(
  meta$clinical_stage
)

effects_current <- contrast_means(
  rna_mat_adj,
  stage
)

effects_batch <- contrast_means(
  rna_batch_adj,
  stage
)

names(effects_current)[
  -1
] <- paste0(
  names(effects_current)[-1],
  "_current"
)

names(effects_batch)[
  -1
] <- paste0(
  names(effects_batch)[-1],
  "_batch_adjusted"
)

effect_table <- merge(
  effects_current,
  effects_batch,
  by = "feature_id",
  all = FALSE
)

write.csv(
  effect_table,
  file.path(
    OUTDIR,
    "RNA_stage_effects_current_vs_batch_adjusted_all_genes.csv"
  ),
  row.names = FALSE
)

# Create un-suffixed copies for generic comparison helper.
current_for_compare <- data.frame(
  feature_id =
    effect_table$feature_id,

  MCI_vs_NCI =
    effect_table$MCI_vs_NCI_current,

  AD_vs_MCI =
    effect_table$AD_vs_MCI_current,

  AD_vs_NCI =
    effect_table$AD_vs_NCI_current
)

batch_for_compare <- data.frame(
  feature_id =
    effect_table$feature_id,

  MCI_vs_NCI =
    effect_table$MCI_vs_NCI_batch_adjusted,

  AD_vs_MCI =
    effect_table$AD_vs_MCI_batch_adjusted,

  AD_vs_NCI =
    effect_table$AD_vs_NCI_batch_adjusted
)

stage_robustness <- do.call(
  rbind,
  lapply(
    c(
      "MCI_vs_NCI",
      "AD_vs_MCI",
      "AD_vs_NCI"
    ),
    function(contrast) {

      compare_effects(
        current_for_compare,
        batch_for_compare,
        contrast
      )
    }
  )
)

write.csv(
  stage_robustness,
  file.path(
    OUTDIR,
    "RNA_stage_effect_robustness_summary.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("GENOME-WIDE STAGE-EFFECT ROBUSTNESS\n")
cat("======================================================================\n")

print(
  stage_robustness,
  row.names = FALSE
)


# =============================================================================
# Identify canonical RNA Hsp60/10 client set from production objects
# =============================================================================

cat("\n======================================================================\n")
cat("HSP60/10 CLIENT-SET DETECTION\n")
cat("======================================================================\n")

object_names <- ls(
  envir = .GlobalEnv
)

client_object_names <- object_names[
  grepl(
    "client",
    object_names,
    ignore.case = TRUE
  )
]

extract_candidate_genes <- function(obj) {

  if (is.character(obj)) {
    return(
      unique(obj)
    )
  }

  if (is.factor(obj)) {
    return(
      unique(
        as.character(obj)
      )
    )
  }

  if (is.data.frame(obj)) {

    likely <- names(obj)[
      grepl(
        "gene|symbol",
        names(obj),
        ignore.case = TRUE
      )
    ]

    if (length(likely) == 0) {
      return(character())
    }

    vals <- unique(
      unlist(
        lapply(
          obj[
            ,
            likely,
            drop = FALSE
          ],
          as.character
        )
      )
    )

    return(vals)
  }

  character()
}

candidate_sets <- list()

for (nm in client_object_names) {

  obj <- get(
    nm,
    envir = .GlobalEnv
  )

  genes <- extract_candidate_genes(
    obj
  )

  overlap <- intersect(
    genes,
    rownames(rna_mat_raw)
  )

  if (
    length(overlap) >= 250 &&
    length(overlap) <= 350
  ) {

    candidate_sets[[nm]] <-
      sort(
        unique(overlap)
      )
  }
}

if (length(candidate_sets) == 0) {

  cat(
    "No candidate RNA client set of expected size detected automatically.\n"
  )

} else {

  candidate_summary <- data.frame(
    object =
      names(candidate_sets),

    n_RNA_genes =
      vapply(
        candidate_sets,
        length,
        integer(1)
      ),

    stringsAsFactors = FALSE
  )

  print(
    candidate_summary,
    row.names = FALSE
  )
}


# =============================================================================
# Use a unique 297-gene RNA client set if available
# =============================================================================

sets_297 <- candidate_sets[
  vapply(
    candidate_sets,
    length,
    integer(1)
  ) == 297
]

if (length(sets_297) > 0) {

  signatures <- vapply(
    sets_297,
    function(x) {
      paste(
        sort(x),
        collapse = "|"
      )
    },
    character(1)
  )

  unique_signatures <- unique(
    signatures
  )

  if (length(unique_signatures) != 1) {

    cat(
      "\nMultiple distinct 297-gene client sets were detected.\n",
      "Stopping Hsp-specific portion rather than guessing.\n",
      sep = ""
    )

  } else {

    hsp_rna <- sets_297[[1]]

    cat(
      "\nCanonical detected RNA Hsp client set:",
      length(hsp_rna),
      "genes\n"
    )

    hsp_idx <- match(
      hsp_rna,
      rownames(rna_mat_raw)
    )

    current_hsp_score <- colMeans(
      t(
        scale(
          t(
            rna_mat_adj[
              hsp_idx,
              ,
              drop = FALSE
            ]
          )
        )
      ),
      na.rm = TRUE
    )

    batch_hsp_score <- colMeans(
      t(
        scale(
          t(
            rna_batch_adj[
              hsp_idx,
              ,
              drop = FALSE
            ]
          )
        )
      ),
      na.rm = TRUE
    )

    hsp_stage <- data.frame(
      stage =
        c(
          "NCI",
          "MCI",
          "AD"
        ),

      current_score =
        vapply(
          c(
            "NCI",
            "MCI",
            "AD"
          ),
          function(s) {
            mean(
              current_hsp_score[
                stage == s
              ]
            )
          },
          numeric(1)
        ),

      batch_adjusted_score =
        vapply(
          c(
            "NCI",
            "MCI",
            "AD"
          ),
          function(s) {
            mean(
              batch_hsp_score[
                stage == s
              ]
            )
          },
          numeric(1)
        ),

      stringsAsFactors = FALSE
    )

    write.csv(
      hsp_stage,
      file.path(
        OUTDIR,
        "RNA_Hsp_pathway_stage_means_current_vs_batch_adjusted.csv"
      ),
      row.names = FALSE
    )

    cat(
      "\nHsp pathway stage means:\n"
    )

    print(
      hsp_stage,
      row.names = FALSE
    )
  }

} else {

  cat(
    "\nNo unique 297-gene RNA Hsp client object was automatically detected.\n"
  )
}


# =============================================================================
# Save sensitivity matrix and provenance
# =============================================================================

saveRDS(
  rna_batch_adj,
  file.path(
    OUTDIR,
    "rna_mat_batch_adjusted_sensitivity.rds"
  )
)

decision_metrics <- data.frame(
  metric = c(
    "max_abs_current_adjustment_reconstruction_difference",
    "median_partial_batch_R2_all_genes",
    "fraction_genes_partial_batch_R2_gt_0.01",
    "fraction_genes_partial_batch_R2_gt_0.05",
    "fraction_genes_partial_batch_R2_gt_0.10",
    "median_batch_R2_current_adjusted_2000genes",
    "median_batch_R2_batch_adjusted_2000genes"
  ),

  value = c(
    max_abs_reconstruction_difference,

    median(
      incremental$partial_batch_R2
    ),

    mean(
      incremental$partial_batch_R2 > 0.01
    ),

    mean(
      incremental$partial_batch_R2 > 0.05
    ),

    mean(
      incremental$partial_batch_R2 > 0.10
    ),

    median(
      r2_current
    ),

    median(
      r2_batch_adjusted
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  decision_metrics,
  file.path(
    OUTDIR,
    "RNA_batch_adjustment_decision_metrics.csv"
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

cat("\n======================================================================\n")
cat("DECISION METRICS\n")
cat("======================================================================\n")

print(
  decision_metrics,
  row.names = FALSE
)

cat(
  "\nRNA BATCH-ADJUSTMENT SENSITIVITY COMPLETED.\n"
)
