options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  source(f)
}

cat("======================================================================\n")
cat("CANONICAL RNA BATCH-ADJUSTMENT VALIDATION\n")
cat("======================================================================\n")

checks <- data.frame(
  check = character(),
  passed = logical(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_check <- function(name, passed, detail = "") {
  checks <<- rbind(
    checks,
    data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  )
}

# ---------------------------------------------------------------------
# Cohort
# ---------------------------------------------------------------------

add_check(
  "RNA matrix n = 577",
  ncol(rna_mat_raw) == 577L,
  ncol(rna_mat_raw)
)

add_check(
  "RNA metadata n = 577",
  nrow(rna_meta_adj) == 577L,
  nrow(rna_meta_adj)
)

add_check(
  "492_120515 excluded from matrix",
  !"492_120515" %in% colnames(rna_mat_raw)
)

add_check(
  "492_120515 excluded from metadata",
  !"492_120515" %in% rna_meta_adj$sample_id
)

stage_counts <- vapply(
  c("NCI", "MCI", "AD"),
  function(x) {
    sum(as.character(rna_meta_adj$clinical_stage) == x)
  },
  integer(1)
)

add_check(
  "Stage counts = 200 / 158 / 219",
  identical(
    unname(stage_counts),
    c(200L, 158L, 219L)
  ),
  paste(
    names(stage_counts),
    stage_counts,
    sep = "=",
    collapse = ", "
  )
)

# ---------------------------------------------------------------------
# Batch / covariates
# ---------------------------------------------------------------------

add_check(
  "9 sequencing-batch levels",
  nlevels(rna_meta_adj$sequencing_batch_factor) == 9L,
  nlevels(rna_meta_adj$sequencing_batch_factor)
)

expected_covars <- c(
  "age_num",
  "sex_factor",
  "pmi_num",
  "rin_num",
  "sequencing_batch_factor"
)

add_check(
  "Canonical RNA covariates",
  identical(rna_covars, expected_covars),
  paste(rna_covars, collapse = ", ")
)

# ---------------------------------------------------------------------
# Adjusted matrix
# ---------------------------------------------------------------------

add_check(
  "Adjusted RNA dimensions = 21433 x 577",
  identical(
    dim(rna_mat_adj),
    c(21433L, 577L)
  ),
  paste(dim(rna_mat_adj), collapse = " x ")
)

add_check(
  "Adjusted RNA has no missing values",
  !anyNA(rna_mat_adj),
  paste("NA =", sum(is.na(rna_mat_adj)))
)

add_check(
  "Adjusted RNA differs from raw RNA",
  !isTRUE(all.equal(rna_mat_adj, rna_mat_raw))
)

# ---------------------------------------------------------------------
# Verify residual batch signal is removed
# ---------------------------------------------------------------------

batch <- rna_meta_adj$sequencing_batch_factor

audit_idx <- unique(
  round(
    seq(
      1,
      nrow(rna_mat_adj),
      length.out = 2000
    )
  )
)

batch_r2 <- vapply(
  audit_idx,
  function(i) {
    summary(
      lm(
        as.numeric(rna_mat_adj[i, ]) ~ batch
      )
    )$r.squared
  },
  numeric(1)
)

median_batch_r2 <- median(batch_r2)

add_check(
  "Residual batch signal essentially removed",
  median_batch_r2 < 1e-10,
  sprintf("median batch R2 = %.3e", median_batch_r2)
)

# ---------------------------------------------------------------------
# Old vs new RNA stage-effect comparison
# ---------------------------------------------------------------------

old_path <- paste0(
  "outputs/reviewer_revisions/RNA_batch_canonicalization/",
  "pre_patch_snapshot/rna_mat_adj.rds"
)

if (file.exists(old_path)) {

  old <- readRDS(old_path)

  shared_samples <- intersect(
    colnames(old),
    colnames(rna_mat_adj)
  )

  shared_genes <- intersect(
    rownames(old),
    rownames(rna_mat_adj)
  )

  old <- old[
    shared_genes,
    shared_samples,
    drop = FALSE
  ]

  new <- rna_mat_adj[
    shared_genes,
    shared_samples,
    drop = FALSE
  ]

  meta_shared <- rna_meta_adj[
    match(
      shared_samples,
      rna_meta_adj$sample_id
    ),
    ,
    drop = FALSE
  ]

  stage <- as.character(
    meta_shared$clinical_stage
  )

  get_effect <- function(mat, a, b) {
    rowMeans(
      mat[, stage == a, drop = FALSE]
    ) -
      rowMeans(
        mat[, stage == b, drop = FALSE]
      )
  }

  contrasts <- list(
    MCI_vs_NCI = c("MCI", "NCI"),
    AD_vs_MCI = c("AD", "MCI"),
    AD_vs_NCI = c("AD", "NCI")
  )

  comparison <- do.call(
    rbind,
    lapply(
      names(contrasts),
      function(nm) {

        a <- contrasts[[nm]][1]
        b <- contrasts[[nm]][2]

        old_eff <- get_effect(old, a, b)
        new_eff <- get_effect(new, a, b)

        data.frame(
          contrast = nm,
          n_genes = length(old_eff),
          spearman = cor(
            old_eff,
            new_eff,
            method = "spearman"
          ),
          pearson = cor(
            old_eff,
            new_eff
          ),
          sign_agreement = mean(
            sign(old_eff) == sign(new_eff)
          ),
          median_abs_change = median(
            abs(new_eff - old_eff)
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  cat("\n======================================================================\n")
  cat("OLD VS NEW RNA STAGE EFFECTS\n")
  cat("======================================================================\n")

  print(
    comparison,
    row.names = FALSE
  )

  dir.create(
    "outputs/reviewer_revisions/RNA_batch_canonicalization",
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    comparison,
    paste0(
      "outputs/reviewer_revisions/RNA_batch_canonicalization/",
      "old_vs_new_RNA_stage_effect_summary.csv"
    ),
    row.names = FALSE
  )

} else {

  cat(
    "\nNOTE: pre-patch rna_mat_adj.rds snapshot not found; ",
    "skipping old/new comparison.\n",
    sep = ""
  )
}

# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

cat("\n======================================================================\n")
cat("VALIDATION CHECKS\n")
cat("======================================================================\n")

print(
  checks,
  row.names = FALSE
)

dir.create(
  "outputs/reviewer_revisions/RNA_batch_canonicalization",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  checks,
  paste0(
    "outputs/reviewer_revisions/RNA_batch_canonicalization/",
    "canonical_RNA_batch_validation.csv"
  ),
  row.names = FALSE
)

if (any(!checks$passed)) {
  stop(
    "One or more canonical RNA batch-adjustment checks failed.",
    call. = FALSE
  )
}

cat("\nALL CANONICAL RNA BATCH-ADJUSTMENT CHECKS PASSED.\n")
