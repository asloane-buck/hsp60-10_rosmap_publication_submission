options(stringsAsFactors = FALSE)

old_dir <- paste0(
  "outputs/reviewer_revisions/RNA_batch_canonicalization/",
  "pre_batch_conventional_DE"
)

new_dir <- paste0(
  "outputs/reviewer_revisions/conventional_stage_differential_final"
)

find_one <- function(dir, pattern) {
  hits <- list.files(
    dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(hits) != 1L) {
    stop(
      "Expected exactly one match for ",
      pattern,
      " in ",
      dir,
      "; found ",
      length(hits),
      call. = FALSE
    )
  }

  hits
}

old_file <- find_one(
  old_dir,
  "RNA_DESeq2_all_contrasts\\.csv$"
)

new_file <- find_one(
  new_dir,
  "RNA_DESeq2_all_contrasts\\.csv$"
)

old <- read.csv(
  old_file,
  stringsAsFactors = FALSE
)

new <- read.csv(
  new_file,
  stringsAsFactors = FALSE
)

required <- c(
  "contrast",
  "feature_id",
  "gene_symbol",
  "effect",
  "p_value",
  "fdr"
)

for (nm in required) {
  if (!nm %in% names(old)) {
    stop("Old table missing: ", nm)
  }

  if (!nm %in% names(new)) {
    stop("New table missing: ", nm)
  }
}

contrasts <- c(
  "MCI_vs_NCI",
  "AD_vs_MCI",
  "AD_vs_NCI"
)

out <- list()

cat("======================================================================\n")
cat("PRE- VS POST-BATCH RNA DIFFERENTIAL ANALYSIS\n")
cat("======================================================================\n")

for (ct in contrasts) {

  o <- old[
    old$contrast == ct,
    ,
    drop = FALSE
  ]

  n <- new[
    new$contrast == ct,
    ,
    drop = FALSE
  ]

  m <- merge(
    o[, c(
      "feature_id",
      "gene_symbol",
      "effect",
      "fdr"
    )],
    n[, c(
      "feature_id",
      "gene_symbol",
      "effect",
      "fdr"
    )],
    by = "feature_id",
    suffixes = c("_old", "_new")
  )

  old_sig <- o$feature_id[
    !is.na(o$fdr) &
      o$fdr < 0.05
  ]

  new_sig <- n$feature_id[
    !is.na(n$fdr) &
      n$fdr < 0.05
  ]

  overlap <- intersect(
    old_sig,
    new_sig
  )

  row <- data.frame(
    contrast = ct,
    n_old_features = nrow(o),
    n_new_features = nrow(n),
    n_shared_features = nrow(m),

    effect_spearman = cor(
      m$effect_old,
      m$effect_new,
      method = "spearman",
      use = "complete.obs"
    ),

    effect_pearson = cor(
      m$effect_old,
      m$effect_new,
      method = "pearson",
      use = "complete.obs"
    ),

    effect_sign_agreement = mean(
      sign(m$effect_old) ==
        sign(m$effect_new),
      na.rm = TRUE
    ),

    old_fdr_sig = length(old_sig),
    new_fdr_sig = length(new_sig),
    both_fdr_sig = length(overlap),
    old_only_fdr_sig = length(
      setdiff(old_sig, new_sig)
    ),
    new_only_fdr_sig = length(
      setdiff(new_sig, old_sig)
    ),

    pct_old_sig_retained = if (
      length(old_sig) > 0
    ) {
      100 * length(overlap) /
        length(old_sig)
    } else {
      NA_real_
    },

    stringsAsFactors = FALSE
  )

  out[[ct]] <- row

  cat("\n", ct, "\n", sep = "")
  print(row, row.names = FALSE)

  old_only <- setdiff(
    old_sig,
    new_sig
  )

  new_only <- setdiff(
    new_sig,
    old_sig
  )

  make_changed_rows <- function(ids, status_label) {
    if (length(ids) == 0L) {
      return(
        data.frame(
          contrast = character(),
          status = character(),
          feature_id = character(),
          stringsAsFactors = FALSE
        )
      )
    }

    data.frame(
      contrast = rep(ct, length(ids)),
      status = rep(status_label, length(ids)),
      feature_id = ids,
      stringsAsFactors = FALSE
    )
  }

  changed <- rbind(
    make_changed_rows(old_only, "old_only"),
    make_changed_rows(new_only, "new_only")
  )

  if (nrow(changed) > 0) {

    changed <- merge(
      changed,
      m,
      by = "feature_id",
      all.x = TRUE
    )

    write.csv(
      changed,
      file.path(
        "outputs/reviewer_revisions/RNA_batch_canonicalization",
        paste0(
          "DE_changed_significance_",
          ct,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
}

summary_tbl <- do.call(
  rbind,
  out
)

write.csv(
  summary_tbl,
  paste0(
    "outputs/reviewer_revisions/RNA_batch_canonicalization/",
    "pre_vs_post_batch_RNA_DE_summary.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("SUMMARY\n")
cat("======================================================================\n")

print(
  summary_tbl,
  row.names = FALSE
)

cat("\nCOMPARISON COMPLETE.\n")
