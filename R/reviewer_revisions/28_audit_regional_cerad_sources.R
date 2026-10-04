############################################################
## 28_audit_regional_cerad_sources.R
##
## READ-ONLY diagnostic.
##
## Goal:
##   Identify the correct CERAD source for the 1,669-person
##   regional cohort without writing or modifying production
##   metadata.
##
## A candidate is useful only if it:
##   1. has a valid participant ID,
##   2. overlaps the regional cohort,
##   3. contains a CERAD field,
##   4. provides high CERAD coverage for the regional IDs,
##   5. does not have conflicting duplicate IDs.
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

regional_file <- paste0(
  "outputs/reviewer_revisions/",
  "regional_braak_fix/",
  "regional_participant_metadata_corrected_braak.csv"
)

regional <- read_csv(
  regional_file,
  show_col_types = FALSE
)

regional_ids <- unique(
  trimws(as.character(regional$individual_id))
)

regional_ids <- regional_ids[
  !is.na(regional_ids) &
    nzchar(regional_ids)
]

cat("\n============================================================\n")
cat("REGIONAL CERAD SOURCE AUDIT — READ ONLY\n")
cat("============================================================\n\n")

cat(
  "Regional cohort individuals:",
  length(regional_ids),
  "\n\n"
)

############################################################
## Candidate roots
############################################################

roots <- unique(
  c(
    ".",
    "..",
    "../Brain Region Specificity",
    "data",
    "outputs"
  )
)

roots <- roots[
  dir.exists(roots)
]

############################################################
## Candidate files
############################################################

all_files <- unique(
  unlist(
    lapply(
      roots,
      function(root) {
        list.files(
          root,
          recursive = TRUE,
          full.names = TRUE
        )
      }
    ),
    use.names = FALSE
  )
)

all_files <- all_files[
  grepl(
    "\\.(csv|rds)$",
    all_files,
    ignore.case = TRUE
  )
]

############################################################
## Restrict to plausible metadata files
############################################################

candidate_files <- all_files[
  grepl(
    paste(
      c(
        "metadata",
        "clinical",
        "analysis_meta",
        "cerad",
        "phenotype",
        "pathology",
        "individual"
      ),
      collapse = "|"
    ),
    basename(all_files),
    ignore.case = TRUE
  )
]

candidate_files <- unique(candidate_files)

cat(
  "Candidate metadata files discovered:",
  length(candidate_files),
  "\n\n"
)

############################################################
## Column aliases
############################################################

id_aliases <- c(
  "individual_id",
  "individualID",
  "IndividualID",
  "projid",
  "ProjID",
  "participant_id",
  "ParticipantID",
  "participant",
  "subject_id",
  "SubjectID"
)

cerad_aliases <- c(
  "cerad",
  "CERAD",
  "cerad_num",
  "CERAD_num",
  "ceradsc",
  "CERADSC",
  "cerad_score",
  "CERAD.Score",
  "cerad_meta"
)

############################################################
## Analyze one data frame
############################################################

audit_dataframe <- function(
  df,
  source,
  object_name = NA_character_
) {

  n_rows <- nrow(df)
  cols <- names(df)

  id_hits <- intersect(id_aliases, cols)
  cerad_hits <- intersect(cerad_aliases, cols)

  if (length(id_hits) == 0L || length(cerad_hits) == 0L) {
    return(NULL)
  }

  results <- list()

  for (id_col in id_hits) {
    ids <- trimws(as.character(df[[id_col]]))

    for (cerad_col in cerad_hits) {

      cerad <- suppressWarnings(
        as.numeric(df[[cerad_col]])
      )

      tmp <- tibble(
        .id = ids,
        .cerad = cerad
      ) |>
        filter(
          !is.na(.id),
          nzchar(.id)
        )

      duplicate_ids <- tmp |>
        count(.id) |>
        filter(n > 1)

      conflicting_ids <- tmp |>
        filter(
          .id %in% duplicate_ids$.id
        ) |>
        group_by(.id) |>
        summarise(
          n_valid_cerad = n_distinct(
            .cerad[
              !is.na(.cerad) &
                .cerad %in% 1:4
            ]
          ),
          .groups = "drop"
        ) |>
        filter(n_valid_cerad > 1)

      ## Collapse to one row per ID for coverage calculation.
      tmp_unique <- tmp |>
        group_by(.id) |>
        summarise(
          cerad_num = {
            x <- .cerad[
              !is.na(.cerad) &
                .cerad %in% 1:4
            ]

            if (length(x) == 0L) {
              NA_real_
            } else {
              unique(x)[1]
            }
          },
          .groups = "drop"
        )

      regional_match <- tmp_unique |>
        filter(.id %in% regional_ids)

      overlap <- nrow(regional_match)

      cerad_valid <- sum(
        !is.na(regional_match$cerad_num)
      )

      coverage <- cerad_valid / length(regional_ids)

      results[[length(results) + 1L]] <- tibble(
        source = source,
        object = object_name,
        rows = n_rows,
        id_column = id_col,
        cerad_column = cerad_col,
        unique_ids = nrow(tmp_unique),
        regional_overlap = overlap,
        regional_overlap_pct =
          100 * overlap / length(regional_ids),
        valid_cerad_in_regional = cerad_valid,
        cerad_coverage_pct =
          100 * coverage,
        duplicate_id_count =
          nrow(duplicate_ids),
        conflicting_duplicate_id_count =
          nrow(conflicting_ids)
      )
    }
  }

  bind_rows(results)
}

############################################################
## CSV audit
############################################################

results <- list()

cat("===== CSV SOURCES =====\n")

for (f in candidate_files[grepl(
  "\\.csv$",
  candidate_files,
  ignore.case = TRUE
)]) {

  cat(
    "Inspecting:",
    f,
    "\n"
  )

  dat <- tryCatch(
    read_csv(
      f,
      n_max = 0,
      show_col_types = FALSE,
      name_repair = "minimal"
    ),
    error = function(e) NULL
  )

  if (is.null(dat)) {
    next
  }

  cols <- names(dat)

  id_hits <- intersect(
    id_aliases,
    cols
  )

  cerad_hits <- intersect(
    cerad_aliases,
    cols
  )

  if (
    length(id_hits) == 0L ||
    length(cerad_hits) == 0L
  ) {
    next
  }

  needed_cols <- unique(
    c(id_hits, cerad_hits)
  )

  dat <- tryCatch(
    read_csv(
      f,
      col_select = all_of(needed_cols),
      show_col_types = FALSE,
      name_repair = "minimal"
    ),
    error = function(e) NULL
  )

  if (is.null(dat)) {
    next
  }

  x <- audit_dataframe(
    dat,
    source = f
  )

  if (!is.null(x)) {
    results[[length(results) + 1L]] <- x
  }
}

############################################################
## RDS audit
############################################################

cat("\n===== RDS SOURCES =====\n")

rds_candidates <- candidate_files[
  grepl(
    "\\.rds$",
    candidate_files,
    ignore.case = TRUE
  )
]

for (f in rds_candidates) {

  cat(
    "Inspecting:",
    f,
    "\n"
  )

  obj <- tryCatch(
    readRDS(f),
    error = function(e) NULL
  )

  if (is.null(obj)) {
    next
  }

  ## Direct data frame / tibble.
  if (is.data.frame(obj)) {

    x <- audit_dataframe(
      obj,
      source = f,
      object_name = "RDS"
    )

    if (!is.null(x)) {
      results[[length(results) + 1L]] <- x
    }

  } else if (is.list(obj)) {

    ## Inspect named data-frame components.
    for (nm in names(obj)) {

      component <- obj[[nm]]

      if (!is.data.frame(component)) {
        next
      }

      x <- audit_dataframe(
        component,
        source = f,
        object_name = nm
      )

      if (!is.null(x)) {
        results[[length(results) + 1L]] <- x
      }
    }
  }
}

############################################################
## Final report
############################################################

if (length(results) == 0L) {

  cat(
    "\nNO CANDIDATE SOURCE CONTAINED BOTH A RECOGNIZED ",
    "PARTICIPANT ID AND CERAD FIELD.\n\n",
    sep = ""
  )

  cat(
    "Next step would be to inspect the metadata directory ",
    "manually rather than guessing.\n"
  )

  quit(
    save = "no",
    status = 0
  )
}

audit <- bind_rows(results) |>
  arrange(
    desc(.data$cerad_coverage_pct),
    desc(.data$regional_overlap_pct),
    .data$conflicting_duplicate_id_count
  )

cat("\n============================================================\n")
cat("CERAD SOURCE AUDIT RESULTS\n")
cat("============================================================\n\n")

print(
  audit,
  n = nrow(audit),
  width = Inf
)

cat("\n============================================================\n")
cat("HIGH-COVERAGE CANDIDATES\n")
cat("============================================================\n\n")

high <- audit |>
  filter(
    regional_overlap_pct >= 90,
    cerad_coverage_pct >= 90,
    conflicting_duplicate_id_count == 0
  )

if (nrow(high) == 0L) {

  cat(
    "NO SAFE >=90% CERAD SOURCE FOUND.\n\n"
  )

  cat(
    "Do NOT write the production regional joint-pathology ",
    "metadata yet.\n",
    sep = ""
  )

} else {

  print(
    high,
    n = nrow(high),
    width = Inf
  )

  cat(
    "\nPASS: at least one candidate provides >=90% CERAD ",
    "coverage without conflicting duplicate IDs.\n",
    sep = ""
  )
}

cat(
  "\nNo production files were written.\n"
)
