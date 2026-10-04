#!/usr/bin/env Rscript

path <- "R/reviewer_revisions/21_run_conventional_stage_differential_analysis.R"

if (!file.exists(path)) {
  stop("Missing ", path, call. = FALSE)
}

txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

checks <- data.frame(
  check = c(
    "Uses basename() on raw count sample IDs",
    "Strips terminal .bam from raw count sample IDs",
    "Checks duplicate IDs after normalization",
    "Reports post-harmonization overlap"
  ),
  passed = c(
    grepl("basename(raw_count_sample_ids_original)", txt, fixed = TRUE),
    grepl('sub("\\\\.bam$", "", raw_count_sample_ids_normalized)', txt, fixed = TRUE),
    grepl("anyDuplicated(raw_count_sample_ids_normalized)", txt, fixed = TRUE),
    grepl("RNA production sample IDs matched after harmonization", txt, fixed = TRUE)
  )
)

print(checks, row.names = FALSE)

if (any(!checks$passed)) {
  stop("Raw-count ID harmonization patch validation FAILED.", call. = FALSE)
}

cat("\nRAW-COUNT ID HARMONIZATION PATCH VALIDATED.\n")
