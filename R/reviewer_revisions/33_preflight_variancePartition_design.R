options(stringsAsFactors = FALSE)

OUTDIR <- "outputs/reviewer_revisions/variancePartition_input_audit"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  source(f)
}

cat("======================================================================\n")
cat("VARIANCEPARTITION DESIGN PREFLIGHT\n")
cat("======================================================================\n")


# =============================================================================
# Helpers
# =============================================================================

matrix_summary <- function(x, name) {

  # Use up to first 1,000 features for distribution diagnostics.
  nr <- min(1000, nrow(x))

  probe <- as.numeric(
    x[seq_len(nr), , drop = FALSE]
  )

  probe <- probe[
    is.finite(probe)
  ]

  data.frame(
    object = name,
    n_features = nrow(x),
    n_samples = ncol(x),

    n_missing = sum(is.na(x)),
    fraction_missing = mean(is.na(x)),

    min = min(x, na.rm = TRUE),
    q01 = unname(
      quantile(probe, 0.01, na.rm = TRUE)
    ),
    median = median(probe, na.rm = TRUE),
    q99 = unname(
      quantile(probe, 0.99, na.rm = TRUE)
    ),
    max = max(x, na.rm = TRUE),

    fraction_zero_probe =
      mean(probe == 0),

    fraction_negative_probe =
      mean(probe < 0),

    fraction_integer_like_probe =
      mean(
        abs(
          probe - round(probe)
        ) < 1e-8
      ),

    stringsAsFactors = FALSE
  )
}


one_to_one <- function(a, b) {

  d <- unique(
    data.frame(
      a = as.character(a),
      b = as.character(b),
      stringsAsFactors = FALSE
    )
  )

  a_to_b <- tapply(
    d$b,
    d$a,
    function(x) {
      length(unique(x))
    }
  )

  b_to_a <- tapply(
    d$a,
    d$b,
    function(x) {
      length(unique(x))
    }
  )

  all(a_to_b == 1) &&
    all(b_to_a == 1)
}


# =============================================================================
# 1. Matrix scale
# =============================================================================

scale_summary <- rbind(
  matrix_summary(
    rna_mat_raw,
    "rna_mat_raw"
  ),
  matrix_summary(
    rna_mat,
    "rna_mat"
  ),
  matrix_summary(
    rna_mat_adj,
    "rna_mat_adj"
  ),
  matrix_summary(
    prot_mat_raw,
    "prot_mat_raw"
  ),
  matrix_summary(
    prot_mat,
    "prot_mat"
  ),
  matrix_summary(
    prot_mat_adj,
    "prot_mat_adj"
  )
)

cat("\n======================================================================\n")
cat("MATRIX SCALE / MISSINGNESS\n")
cat("======================================================================\n")

print(
  scale_summary,
  row.names = FALSE
)

write.csv(
  scale_summary,
  file.path(
    OUTDIR,
    "matrix_scale_preflight.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 2. Are similarly named matrices actually identical?
# =============================================================================

object_comparison <- data.frame(
  comparison = c(
    "rna_mat_raw vs rna_mat",
    "rna_mat_raw vs rna_mat_adj",
    "rna_mat vs rna_mat_adj",
    "prot_mat_raw vs prot_mat",
    "prot_mat_raw vs prot_mat_adj",
    "prot_mat vs prot_mat_adj"
  ),

  identical = c(
    identical(rna_mat_raw, rna_mat),
    identical(rna_mat_raw, rna_mat_adj),
    identical(rna_mat, rna_mat_adj),
    identical(prot_mat_raw, prot_mat),
    identical(prot_mat_raw, prot_mat_adj),
    identical(prot_mat, prot_mat_adj)
  ),

  all_equal = c(
    isTRUE(all.equal(rna_mat_raw, rna_mat)),
    isTRUE(all.equal(rna_mat_raw, rna_mat_adj)),
    isTRUE(all.equal(rna_mat, rna_mat_adj)),
    isTRUE(all.equal(prot_mat_raw, prot_mat)),
    isTRUE(all.equal(prot_mat_raw, prot_mat_adj)),
    isTRUE(all.equal(prot_mat, prot_mat_adj))
  ),

  stringsAsFactors = FALSE
)

cat("\n======================================================================\n")
cat("OBJECT IDENTITY\n")
cat("======================================================================\n")

print(
  object_comparison,
  row.names = FALSE
)

write.csv(
  object_comparison,
  file.path(
    OUTDIR,
    "matrix_object_identity.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 3. Exact matrix / metadata ordering
# =============================================================================

alignment <- data.frame(
  modality = c(
    "RNA",
    "protein"
  ),

  matrix_columns = c(
    ncol(rna_mat_raw),
    ncol(prot_mat_raw)
  ),

  metadata_rows = c(
    nrow(rna_meta),
    nrow(prot_meta_adj)
  ),

  exact_same_order = c(
    identical(
      as.character(
        colnames(rna_mat_raw)
      ),
      as.character(
        rna_meta$sample_id
      )
    ),

    identical(
      as.character(
        colnames(prot_mat_raw)
      ),
      as.character(
        prot_meta_adj$sample_id
      )
    )
  ),

  set_equal = c(
    setequal(
      colnames(rna_mat_raw),
      rna_meta$sample_id
    ),

    setequal(
      colnames(prot_mat_raw),
      prot_meta_adj$sample_id
    )
  ),

  stringsAsFactors = FALSE
)

cat("\n======================================================================\n")
cat("MATRIX / METADATA ORDERING\n")
cat("======================================================================\n")

print(
  alignment,
  row.names = FALSE
)


# =============================================================================
# 4. RNA batch structure
# =============================================================================

cat("\n======================================================================\n")
cat("RNA BATCH STRUCTURE\n")
cat("======================================================================\n")

cat(
  "rna_batch nonmissing:",
  sum(!is.na(rna_meta$rna_batch)),
  "/",
  nrow(rna_meta),
  "\n"
)

cat(
  "library_batch levels:",
  length(unique(rna_meta$library_batch)),
  "\n"
)

cat(
  "sequencing_batch levels:",
  length(unique(rna_meta$sequencing_batch)),
  "\n"
)

cat(
  "library_batch == sequencing_batch rowwise:",
  all(
    as.character(rna_meta$library_batch) ==
      as.character(rna_meta$sequencing_batch)
  ),
  "\n"
)

cat(
  "library_batch <-> sequencing_batch one-to-one:",
  one_to_one(
    rna_meta$library_batch,
    rna_meta$sequencing_batch
  ),
  "\n"
)

cat("\nLibrary batch counts:\n")
print(
  sort(
    table(rna_meta$library_batch),
    decreasing = TRUE
  )
)

cat("\nSequencing batch counts:\n")
print(
  sort(
    table(rna_meta$sequencing_batch),
    decreasing = TRUE
  )
)

cat("\nLibrary × sequencing mapping:\n")

batch_map <- unique(
  data.frame(
    library_batch =
      as.character(
        rna_meta$library_batch
      ),
    sequencing_batch =
      as.character(
        rna_meta$sequencing_batch
      ),
    stringsAsFactors = FALSE
  )
)

batch_map <- batch_map[
  order(
    batch_map$library_batch,
    batch_map$sequencing_batch
  ),
  ,
  drop = FALSE
]

print(
  batch_map,
  row.names = FALSE
)

write.csv(
  batch_map,
  file.path(
    OUTDIR,
    "RNA_library_vs_sequencing_batch_mapping.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 5. Protein batch structure
# =============================================================================

cat("\n======================================================================\n")
cat("PROTEIN BATCH STRUCTURE\n")
cat("======================================================================\n")

cat(
  "tmt_batch levels:",
  length(
    unique(
      prot_meta_adj$tmt_batch
    )
  ),
  "\n"
)

cat(
  "batch levels:",
  length(
    unique(
      prot_meta_adj$batch
    )
  ),
  "\n"
)

cat(
  "tmt_batch == batch rowwise:",
  all(
    as.character(
      prot_meta_adj$tmt_batch
    ) ==
      as.character(
        prot_meta_adj$batch
      )
  ),
  "\n"
)

cat(
  "tmt_batch <-> batch one-to-one:",
  one_to_one(
    prot_meta_adj$tmt_batch,
    prot_meta_adj$batch
  ),
  "\n"
)

cat("\nTMT batch counts:\n")

print(
  sort(
    table(
      prot_meta_adj$tmt_batch
    ),
    decreasing = TRUE
  )
)


# =============================================================================
# 6. Candidate variance-partition cohorts
# =============================================================================

rna_design <- data.frame(
  sample_id =
    rna_meta_adj$sample_id,

  age =
    rna_meta_adj$age_num,

  sex =
    rna_meta_adj$sex_factor,

  pmi =
    rna_meta_adj$pmi_num,

  rin =
    rna_meta_adj$rin_num,

  library_batch =
    rna_meta$library_batch,

  sequencing_batch =
    rna_meta$sequencing_batch,

  clinical_stage =
    rna_meta_adj$clinical_stage,

  braak =
    rna_meta_adj$braak_num,

  cerad =
    rna_meta_adj$cerad_num,

  stringsAsFactors = FALSE
)

prot_design <- data.frame(
  sample_id =
    prot_meta_adj$sample_id,

  age =
    prot_meta_adj$age_num,

  sex =
    prot_meta_adj$sex_factor,

  pmi =
    prot_meta_adj$pmi_num,

  tmt_batch =
    prot_meta_adj$tmt_batch,

  clinical_stage =
    prot_meta_adj$clinical_stage,

  braak =
    prot_meta_adj$braak_num,

  cerad =
    prot_meta_adj$cerad_num,

  stringsAsFactors = FALSE
)


cohort_summary <- data.frame(
  design = c(
    "RNA_nuisance",
    "RNA_nuisance_plus_stage",
    "RNA_nuisance_plus_Braak",
    "RNA_nuisance_plus_CERAD",

    "protein_nuisance",
    "protein_nuisance_plus_stage",
    "protein_nuisance_plus_Braak",
    "protein_nuisance_plus_CERAD"
  ),

  n_complete = c(
    sum(
      complete.cases(
        rna_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "rin",
            "library_batch",
            "sequencing_batch"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        rna_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "rin",
            "library_batch",
            "sequencing_batch",
            "clinical_stage"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        rna_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "rin",
            "library_batch",
            "sequencing_batch",
            "braak"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        rna_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "rin",
            "library_batch",
            "sequencing_batch",
            "cerad"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        prot_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "tmt_batch"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        prot_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "tmt_batch",
            "clinical_stage"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        prot_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "tmt_batch",
            "braak"
          )
        ]
      )
    ),

    sum(
      complete.cases(
        prot_design[
          ,
          c(
            "age",
            "sex",
            "pmi",
            "tmt_batch",
            "cerad"
          )
        ]
      )
    )
  ),

  stringsAsFactors = FALSE
)

cat("\n======================================================================\n")
cat("COMPLETE-CASE DESIGN COUNTS\n")
cat("======================================================================\n")

print(
  cohort_summary,
  row.names = FALSE
)

write.csv(
  cohort_summary,
  file.path(
    OUTDIR,
    "variancePartition_candidate_cohort_counts.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 7. Existing production covariate definitions
# =============================================================================

cat("\n======================================================================\n")
cat("CURRENT PRODUCTION COVARIATE DEFINITIONS\n")
cat("======================================================================\n")

cat(
  "rna_covars:\n  ",
  paste(
    rna_covars,
    collapse = ", "
  ),
  "\n"
)

cat(
  "protein_covars:\n  ",
  paste(
    protein_covars,
    collapse = ", "
  ),
  "\n"
)

cat(
  "required_protein_covars:\n  ",
  paste(
    required_protein_covars,
    collapse = ", "
  ),
  "\n"
)


# =============================================================================
# 8. Installation environment
# =============================================================================

cat("\n======================================================================\n")
cat("INSTALLATION ENVIRONMENT\n")
cat("======================================================================\n")

cat(
  "R version:",
  R.version.string,
  "\n"
)

cat(
  "BiocManager installed:",
  requireNamespace(
    "BiocManager",
    quietly = TRUE
  ),
  "\n"
)

if (
  requireNamespace(
    "BiocManager",
    quietly = TRUE
  )
) {
  cat(
    "Bioconductor version:",
    as.character(
      BiocManager::version()
    ),
    "\n"
  )
}

cat(
  "variancePartition installed:",
  requireNamespace(
    "variancePartition",
    quietly = TRUE
  ),
  "\n"
)


cat("\n======================================================================\n")
cat("PREFLIGHT COMPLETE\n")
cat("======================================================================\n")

cat(
  "Outputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)

cat(
  "\nDo not run variancePartition yet; use this output to select the\n",
  "RNA normalization path and nonredundant batch term.\n",
  sep = ""
)
