options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  source(f)
}

OUTDIR <- "outputs/reviewer_revisions/RNA_batch_adjustment_sensitivity"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("======================================================================\n")
cat("RNA BATCH PROVENANCE AUDIT\n")
cat("======================================================================\n")

batch_cols <- names(rna_meta)[
  grepl(
    "batch|library|sequenc",
    names(rna_meta),
    ignore.case = TRUE
  )
]

id_cols <- names(rna_meta)[
  grepl(
    "sample|individual|projid|subject|participant",
    names(rna_meta),
    ignore.case = TRUE
  )
]

show_cols <- unique(
  c(
    id_cols,
    batch_cols,
    "clinical_stage",
    "rin",
    "pmi",
    "age_death",
    "sex"
  )
)

show_cols <- show_cols[
  show_cols %in% names(rna_meta)
]

cat("\nBatch-related columns:\n")
print(batch_cols)

cat("\nUnique library_batch values:\n")
print(sort(unique(as.character(rna_meta$library_batch))))

cat("\nUnique sequencing_batch values:\n")
print(sort(unique(as.character(rna_meta$sequencing_batch))))

cat("\nRows with unusual batch labels:\n")

weird <- (
  as.character(rna_meta$library_batch) %in%
    c("0", "0, 6, 7")
)

print(
  rna_meta[
    weird,
    show_cols,
    drop = FALSE
  ],
  row.names = FALSE
)

write.csv(
  rna_meta[
    weird,
    show_cols,
    drop = FALSE
  ],
  file.path(
    OUTDIR,
    "RNA_unusual_batch_samples.csv"
  ),
  row.names = FALSE
)

cat("\nCross-tab of every batch-like field against library_batch:\n")

for (col in batch_cols) {

  if (col == "library_batch") {
    next
  }

  cat("\n---", col, "---\n")

  print(
    table(
      library_batch =
        as.character(rna_meta$library_batch),
      other =
        as.character(rna_meta[[col]]),
      useNA = "ifany"
    )
  )
}

cat("\n======================================================================\n")
cat("SOURCE OBJECT SEARCH\n")
cat("======================================================================\n")

# Inspect similarly named columns in the broader original metadata object.
if (exists("analysis_meta")) {

  cols <- names(analysis_meta)[
    grepl(
      "batch|library|sequenc|sample",
      names(analysis_meta),
      ignore.case = TRUE
    )
  ]

  cat(
    "analysis_meta relevant columns:\n"
  )

  print(cols)

  if (
    "sample_id" %in% names(analysis_meta)
  ) {

    target_ids <- as.character(
      rna_meta$sample_id[weird]
    )

    hit <- analysis_meta[
      as.character(
        analysis_meta$sample_id
      ) %in% target_ids,
      cols,
      drop = FALSE
    ]

    cat(
      "\nOriginal analysis_meta rows for unusual samples:\n"
    )

    print(
      hit,
      row.names = FALSE
    )

    write.csv(
      hit,
      file.path(
        OUTDIR,
        "RNA_unusual_batch_samples_analysis_meta.csv"
      ),
      row.names = FALSE
    )
  }
}

cat("\nRNA BATCH PROVENANCE AUDIT COMPLETED.\n")
