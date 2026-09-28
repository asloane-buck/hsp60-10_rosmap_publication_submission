############################################################
## 61_validate_final_reproducibility_gate.R
############################################################

options(stringsAsFactors = FALSE)

OUT <- file.path(
  "outputs",
  "reviewer_revisions",
  "final_reproducibility_run"
)

required <- c(
  "SUCCESS",
  "start_utc.txt",
  "end_utc.txt",
  "git_HEAD_before.txt",
  "git_HEAD_after.txt",
  "git_status_before.txt",
  "git_status_after.txt",
  "R_code_diff_unchanged.txt",
  "phase_status.tsv",
  "sessionInfo_final.txt",
  "installed_packages.csv",
  "scientific_csv_sha256.txt",
  "warning_line_count.txt",
  "renv_lock_present.txt",
  "personal_path_hits_present.txt"
)

missing <- required[
  !file.exists(file.path(OUT, required))
]

if (length(missing)) {
  stop(
    "Missing final reproducibility artifact(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

read_scalar <- function(name) {
  trimws(readLines(file.path(OUT, name), warn = FALSE)[1])
}

phase <- read.delim(
  file.path(OUT, "phase_status.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

checks <- list()

add <- function(name, ok, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(ok),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

add(
  "all_phases_exit_zero",
  nrow(phase) >= 20L &&
    all(phase$status == 0L),
  paste0(
    "phases=", nrow(phase),
    "; nonzero=", sum(phase$status != 0L)
  )
)

add(
  "git_HEAD_unchanged",
  identical(
    read_scalar("git_HEAD_before.txt"),
    read_scalar("git_HEAD_after.txt")
  ),
  paste0(
    "before=", read_scalar("git_HEAD_before.txt"),
    "; after=", read_scalar("git_HEAD_after.txt")
  )
)

add(
  "R_code_diff_unchanged_during_execution",
  identical(
    read_scalar("R_code_diff_unchanged.txt"),
    "TRUE"
  ),
  paste0(
    "value=",
    read_scalar("R_code_diff_unchanged.txt")
  )
)

add(
  "no_tracked_personal_path_hits",
  identical(
    read_scalar("personal_path_hits_present.txt"),
    "FALSE"
  ),
  paste0(
    "value=",
    read_scalar("personal_path_hits_present.txt")
  )
)

required_outputs <- c(
  "outputs/main_figures/tables/all_hsp60_10_client_stage_long_COVARIATE_ADJUSTED.csv",
  "outputs/main_figures/tables/main_fig1_all_clients/Fig1B_all_clients_clean_effect_summary.csv",
  "outputs/main_figures/tables/main_fig3_all_clients/Fig3_all_clients_pathology_table.csv",
  "outputs/main_figures/tables/main_fig4_cognition/Fig4_client_level_cognition_all_models.csv",
  "outputs/main_figures/tables/main_fig5_matched_null_specificity/Fig5_gene_universe_cognition_summary_by_gene.csv",
  "outputs/supplemental_figures/audits/SuppFig2_matched_stage_crosswalk.csv",
  "outputs/reviewer_revisions/reporting_requirements_resolution/10_final_reporting_status.csv",
  "outputs/reviewer_revisions/remaining_analysis_gates/15_remaining_gate_status.csv"
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]

add(
  "key_outputs_present",
  length(missing_outputs) == 0L,
  if (length(missing_outputs)) {
    paste(missing_outputs, collapse = ";")
  } else {
    paste0("n=", length(required_outputs))
  }
)

rna_counts_file <- "outputs/main_figures/tables/rna_primary_cogdx_stage_counts.csv"
protein_counts_file <- "outputs/main_figures/tables/protein_primary_stage_analysis_counts.csv"

if (file.exists(rna_counts_file)) {

  rna <- read.csv(
    rna_counts_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stage_candidates <- intersect(
    c(
      "stage",
      "clinical_stage"
    ),
    names(rna)
  )

  count_candidates <- intersect(
    c(
      "n",
      "count",
      "n_samples",
      "n_participants"
    ),
    names(rna)
  )

  if (
    length(stage_candidates) == 0L ||
    length(count_candidates) == 0L
  ) {

    add(
      "RNA_primary_counts_200_158_219",
      FALSE,
      paste0(
        "Could not identify RNA stage/count columns. Available: ",
        paste(names(rna), collapse = ";")
      )
    )

  } else {

    stage_col <- stage_candidates[1]
    n_col <- count_candidates[1]

    rna_map <- setNames(
      as.integer(rna[[n_col]]),
      toupper(
        trimws(
          as.character(
            rna[[stage_col]]
          )
        )
      )
    )

    observed <- unname(
      rna_map[
        c(
          "NCI",
          "MCI",
          "AD"
        )
      ]
    )

    add(
      "RNA_primary_counts_200_158_219",
      length(observed) == 3L &&
        !anyNA(observed) &&
        identical(
          as.integer(observed),
          c(
            200L,
            158L,
            219L
          )
        ),
      paste0(
        "NCI=", observed[1],
        "; MCI=", observed[2],
        "; AD=", observed[3]
      )
    )
  }

} else {

  add(
    "RNA_primary_counts_200_158_219",
    FALSE,
    "RNA stage-count file missing"
  )
}

if (file.exists(protein_counts_file)) {

  prot <- read.csv(
    protein_counts_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_protein_count_cols <- c(
    "stage",
    "n_analyzable"
  )

  missing_protein_count_cols <- setdiff(
    required_protein_count_cols,
    names(prot)
  )

  if (length(missing_protein_count_cols) > 0L) {

    add(
      "protein_primary_complete_counts_167_96_109",
      FALSE,
      paste0(
        "Missing protein count column(s): ",
        paste(
          missing_protein_count_cols,
          collapse = ";"
        ),
        ". Available: ",
        paste(names(prot), collapse = ";")
      )
    )

  } else {

    prot_map <- setNames(
      as.integer(
        prot$n_analyzable
      ),
      toupper(
        trimws(
          as.character(
            prot$stage
          )
        )
      )
    )

    observed <- unname(
      prot_map[
        c(
          "NCI",
          "MCI",
          "AD"
        )
      ]
    )

    add(
      "protein_primary_complete_counts_167_96_109",
      length(observed) == 3L &&
        !anyNA(observed) &&
        identical(
          as.integer(observed),
          c(
            167L,
            96L,
            109L
          )
        ),
      paste0(
        "NCI=", observed[1],
        "; MCI=", observed[2],
        "; AD=", observed[3]
      )
    )
  }

} else {

  add(
    "protein_primary_complete_counts_167_96_109",
    FALSE,
    "Protein stage-count file missing"
  )
}

matched_file <- "outputs/supplemental_figures/audits/SuppFig2_matched_stage_crosswalk.csv"

if (file.exists(matched_file)) {
  matched <- read.csv(
    matched_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  protein_stage <- toupper(
    trimws(as.character(matched$protein_stage))
  )
  rna_stage <- toupper(
    trimws(as.character(matched$rna_stage))
  )

  counts <- table(
    factor(
      protein_stage,
      levels = c("NCI","MCI","AD")
    )
  )

  add(
    "matched_cohort_198_94_56_48",
    nrow(matched) == 198L &&
      length(unique(matched$individual_id)) == 198L &&
      all(protein_stage == rna_stage) &&
      identical(
        as.integer(counts),
        c(94L,56L,48L)
      ),
    paste0(
      "n=", nrow(matched),
      "; NCI=", counts["NCI"],
      "; MCI=", counts["MCI"],
      "; AD=", counts["AD"]
    )
  )
} else {
  add(
    "matched_cohort_198_94_56_48",
    FALSE,
    "SuppFig2 matched crosswalk missing"
  )
}

report_file <- "outputs/reviewer_revisions/reporting_requirements_resolution/10_final_reporting_status.csv"

if (file.exists(report_file)) {

  report <- read.csv(
    report_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  allowed_reporting_status <- c(
    "resolved",
    "not_applicable"
  )

  status_candidates <- names(report)[
    vapply(
      report,
      function(x) {

        z <- tolower(
          trimws(
            as.character(x)
          )
        )

        z <- z[
          !is.na(z) &
            nzchar(z)
        ]

        length(z) > 0L &&
          all(
            z %in%
              allowed_reporting_status
          )
      },
      logical(1)
    )
  ]

  if (length(status_candidates) == 0L) {

    add(
      "reporting_requirements_zero_unresolved",
      FALSE,
      paste0(
        "Could not identify final reporting-status column. Available: ",
        paste(
          names(report),
          collapse = ";"
        )
      )
    )

  } else {

    status_col <- status_candidates[1]

    status_values <- tolower(
      trimws(
        as.character(
          report[[status_col]]
        )
      )
    )

    add(
      "reporting_requirements_zero_unresolved",
      nrow(report) == 17L &&
        all(
          status_values %in%
            allowed_reporting_status
        ),
      paste0(
        "rows=", nrow(report),
        "; status_col=", status_col,
        "; values=",
        paste(
          sort(
            unique(
              status_values
            )
          ),
          collapse = ";"
        )
      )
    )
  }

} else {

  add(
    "reporting_requirements_zero_unresolved",
    FALSE,
    "reporting final-status file missing"
  )
}

gate_file <- "outputs/reviewer_revisions/remaining_analysis_gates/15_remaining_gate_status.csv"

if (file.exists(gate_file)) {
  gates <- read.csv(
    gate_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  add(
    "remaining_analysis_gates_no_failures",
    nrow(gates) == 4L &&
      !any(gates$analytical_status == "FAILED") &&
      sum(
        gates$analytical_status ==
          "ANALYTICALLY_REVALIDATED"
      ) == 3L &&
      sum(
        gates$analytical_status ==
          "ANALYTICALLY_RESOLVED_MANUSCRIPT_EDIT_PENDING"
      ) == 1L,
    paste(
      gates$gate,
      gates$analytical_status,
      sep = "=",
      collapse = ";"
    )
  )
} else {
  add(
    "remaining_analysis_gates_no_failures",
    FALSE,
    "remaining-gate status file missing"
  )
}

checksum_lines <- readLines(
  file.path(OUT, "scientific_csv_sha256.txt"),
  warn = FALSE
)

add(
  "scientific_checksum_manifest_nonempty",
  length(checksum_lines) >= 50L,
  paste0("entries=", length(checksum_lines))
)

warning_n <- suppressWarnings(
  as.integer(read_scalar("warning_line_count.txt"))
)

add(
  "warning_inventory_created",
  is.finite(warning_n) && warning_n >= 0L,
  paste0("warning/deprecation log lines=", warning_n)
)

renv_present <- read_scalar("renv_lock_present.txt")

add(
  "renv_lock_status_recorded",
  renv_present %in% c("TRUE","FALSE"),
  paste0("renv.lock present=", renv_present)
)

validation <- do.call(rbind, checks)
print(validation, row.names = FALSE)

npass <- sum(validation$passed)
ntotal <- nrow(validation)

cat(
  "\nValidation: ",
  npass,
  " / ",
  ntotal,
  " checks passed\n",
  sep = ""
)

failed <- validation[
  !validation$passed,
  ,
  drop = FALSE
]

write.csv(
  validation,
  file.path(OUT, "final_reproducibility_validation.csv"),
  row.names = FALSE
)

if (nrow(failed)) {
  cat("\nFAILED CHECKS:\n")
  print(failed, row.names = FALSE)

  stop(
    "Final reproducibility validation failed: ",
    nrow(failed),
    " check(s).",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("FINAL ANALYSIS REPRODUCIBILITY GATE PASSED\n")
cat("============================================================\n")
cat(
  "All analytical/reproducibility checks passed.\n",
  "The Fig1B subtraction manuscript/figure edit remains intentionally pending.\n",
  sep = ""
)

if (identical(renv_present, "FALSE")) {
  cat(
    "\nNOTE: renv.lock is not present. This does not invalidate the\n",
    "analysis run, but environment freezing should be completed before\n",
    "the public code archive is finalized.\n",
    sep = ""
  )
}
