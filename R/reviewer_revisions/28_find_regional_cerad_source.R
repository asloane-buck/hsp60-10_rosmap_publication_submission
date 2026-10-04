############################################################
# 28_find_regional_cerad_source.R
#
# Diagnostic only.
#
# Goal:
#   Find the local AMP-AD Diverse Cohorts metadata source that
#   corresponds to the 1,669-person regional cohort and contains
#   CERAD neuropathology.
#
# NO FILES ARE WRITTEN.
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

source("R/00_config.R")
source("R/01_utils.R")

regional_file <- file.path(
  "outputs",
  "reviewer_revisions",
  "regional_braak_fix",
  "regional_participant_metadata_corrected_braak.csv"
)

if (!file.exists(regional_file)) {
  stop(
    "Missing regional metadata:\n",
    regional_file
  )
}

regional <- readr::read_csv(
  regional_file,
  show_col_types = FALSE
)

regional_ids <- unique(
  as.character(regional$individual_id)
)

cat("\n============================================================\n")
cat("REGIONAL CERAD SOURCE DISCOVERY\n")
cat("============================================================\n")

cat(
  "Regional rows:",
  nrow(regional),
  "\n"
)

cat(
  "Regional individuals:",
  length(regional_ids),
  "\n"
)

cat(
  "Example IDs:",
  paste(head(regional_ids, 10), collapse = ", "),
  "\n"
)


############################################################
# Candidate search locations
############################################################

search_dirs <- unique(
  c(
    file.path("data", "raw"),
    file.path("data", "supplemental_inputs"),
    file.path("data", "external"),
    file.path("outputs", "reviewer_revisions"),
    file.path("outputs", "supplemental_figures")
  )
)

search_dirs <- search_dirs[
  dir.exists(search_dirs)
]

cat("\n===== SEARCH DIRECTORIES =====\n")
print(search_dirs)


############################################################
# Find candidate tabular files
############################################################

all_files <- unlist(
  lapply(
    search_dirs,
    function(x) {
      list.files(
        x,
        recursive = TRUE,
        full.names = TRUE,
        include.dirs = FALSE
      )
    }
  ),
  use.names = FALSE
)

all_files <- unique(all_files)

candidate_files <- all_files[
  grepl(
    "\\.(csv|tsv|txt|rds)$",
    all_files,
    ignore.case = TRUE
  )
]

candidate_files <- candidate_files[
  grepl(
    paste(
      c(
        "amp",
        "diverse",
        "nybb",
        "metadata",
        "clinical",
        "phenotype",
        "neuro",
        "path",
        "harmon",
        "braak",
        "cerad",
        "biospecimen"
      ),
      collapse = "|"
    ),
    basename(candidate_files),
    ignore.case = TRUE
  )
]

cat(
  "\nCandidate files:",
  length(candidate_files),
  "\n"
)


############################################################
# Helpers
############################################################

read_candidate <- function(path) {

  ext <- tolower(
    tools::file_ext(path)
  )

  tryCatch({

    if (ext %in% c("tsv", "txt")) {
      readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE
      )
    } else if (ext == "csv") {
      readr::read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE
      )
    } else if (ext == "rds") {
      obj <- readRDS(path)

      if (is.data.frame(obj)) {
        tibble::as_tibble(obj)
      } else if (
        is.list(obj) &&
        !is.null(names(obj))
      ) {
        # Try to find a data-frame element.
        df_names <- names(obj)[
          vapply(
            obj,
            is.data.frame,
            logical(1)
          )
        ]

        if (length(df_names) > 0) {
          tibble::as_tibble(
            obj[[df_names[[1]]]]
          )
        } else {
          NULL
        }
      } else {
        NULL
      }
    } else {
      NULL
    }

  }, error = function(e) {
    NULL
  })
}


id_like_columns <- function(df) {

  nm <- names(df)

  nm[
    grepl(
      "individual|participant|subject|donor|sample|specimen|patient|person|id$|_id$",
      nm,
      ignore.case = TRUE
    )
  ]
}


cerad_columns <- function(df) {

  nm <- names(df)

  nm[
    grepl(
      "cerad",
      nm,
      ignore.case = TRUE
    )
  ]
}


braak_columns <- function(df) {

  nm <- names(df)

  nm[
    grepl(
      "braak",
      nm,
      ignore.case = TRUE
    )
  ]
}


valid_cerad <- function(x) {

  y <- suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )

  !is.na(y) & y %in% 1:4
}


############################################################
# Scan files
############################################################

results <- list()

cat("\n===== SCANNING CANDIDATE FILES =====\n")

for (path in candidate_files) {

  message(
    "Checking: ",
    path
  )

  df <- read_candidate(path)

  if (is.null(df)) {
    next
  }

  if (nrow(df) == 0 || ncol(df) == 0) {
    next
  }

  id_cols <- id_like_columns(df)

  if (length(id_cols) == 0) {
    next
  }

  cerad_cols <- cerad_columns(df)
  braak_cols <- braak_columns(df)

  for (id_col in id_cols) {

    ids <- as.character(
      df[[id_col]]
    )

    ids <- ids[
      !is.na(ids) &
        nzchar(ids)
    ]

    overlap <- intersect(
      regional_ids,
      unique(ids)
    )

    if (length(overlap) == 0) {
      next
    }

    overlap_df <- df[
      as.character(df[[id_col]]) %in% regional_ids,
      ,
      drop = FALSE
    ]

    if (length(cerad_cols) > 0) {

      cerad_counts <- vapply(
        cerad_cols,
        function(x) {
          sum(
            valid_cerad(
              overlap_df[[x]]
            )
          )
        },
        numeric(1)
      )

      best_cerad_n <- max(
        cerad_counts
      )

      best_cerad_col <- cerad_cols[
        which.max(cerad_counts)
      ]

    } else {

      best_cerad_n <- 0
      best_cerad_col <- NA_character_

    }

    results[[length(results) + 1L]] <- tibble(
      path = path,
      n_rows = nrow(df),
      n_columns = ncol(df),
      id_column = id_col,
      regional_overlap = length(overlap),
      regional_overlap_pct =
        100 * length(overlap) / length(regional_ids),
      cerad_column = best_cerad_col,
      cerad_available_in_overlap = best_cerad_n,
      cerad_coverage_pct =
        ifelse(
          length(overlap) > 0,
          100 * best_cerad_n / length(overlap),
          NA_real_
        ),
      n_cerad_columns = length(cerad_cols),
      n_braak_columns = length(braak_cols),
      example_overlap =
        paste(
          head(overlap, 5),
          collapse = "; "
        )
    )
  }
}


############################################################
# Results
############################################################

if (length(results) == 0) {

  cat(
    "\nNO MATCHING METADATA TABLE FOUND.\n\n",
    "This means the AMP-AD Diverse Cohorts metadata needed ",
    "for CERAD is probably not currently present in the ",
    "searched local directories.\n",
    sep = ""
  )

  quit(
    status = 0
  )
}

results_tbl <- bind_rows(
  results
) |>
  arrange(
    desc(regional_overlap),
    desc(cerad_available_in_overlap),
    desc(cerad_coverage_pct)
  )


############################################################
# High-overlap candidates
############################################################

cat("\n============================================================\n")
cat("TOP METADATA CANDIDATES\n")
cat("============================================================\n")

print(
  results_tbl |>
    select(
      path,
      n_rows,
      id_column,
      regional_overlap,
      regional_overlap_pct,
      cerad_column,
      cerad_available_in_overlap,
      cerad_coverage_pct,
      n_braak_columns
    ) |>
    slice_head(n = 30),
  n = 30
)


############################################################
# Strong candidates
############################################################

strong <- results_tbl |>
  filter(
    regional_overlap_pct >= 90
  )

cat("\n============================================================\n")
cat("HIGH-OVERLAP CANDIDATES (>=90% REGIONAL ID OVERLAP)\n")
cat("============================================================\n")

if (nrow(strong) == 0) {

  cat(
    "NONE FOUND.\n",
    "Do not modify the regional metadata yet.\n",
    sep = ""
  )

} else {

  print(
    strong |>
      select(
        path,
        n_rows,
        id_column,
        regional_overlap,
        regional_overlap_pct,
        cerad_column,
        cerad_available_in_overlap,
        cerad_coverage_pct,
        n_braak_columns
      ),
    n = Inf
  )
}


############################################################
# Best CERAD-bearing candidate
############################################################

cerad_best <- results_tbl |>
  filter(
    !is.na(cerad_column),
    cerad_available_in_overlap > 0
  ) |>
  arrange(
    desc(regional_overlap),
    desc(cerad_available_in_overlap),
    desc(cerad_coverage_pct)
  ) |>
  slice_head(n = 10)

cat("\n============================================================\n")
cat("BEST CERAD-BEARING CANDIDATES\n")
cat("============================================================\n")

if (nrow(cerad_best) == 0) {

  cat(
    "NONE FOUND.\n",
    "The correct AMP-AD Diverse Cohorts neuropathology metadata ",
    "does not appear to be present locally.\n",
    sep = ""
  )

} else {

  print(
    cerad_best |>
      select(
        path,
        id_column,
        regional_overlap,
        regional_overlap_pct,
        cerad_column,
        cerad_available_in_overlap,
        cerad_coverage_pct
      ),
    n = 10
  )
}


############################################################
# Explicit NYBB check
############################################################

cat("\n============================================================\n")
cat("EXPLICIT NYBB ID CHECK\n")
cat("============================================================\n")

nybb_id <- "NYBB_17"

nybb_hits <- results_tbl |>
  filter(
    path != regional_file
  )

if (nrow(nybb_hits) > 0) {

  found_nybb <- logical(nrow(nybb_hits))

  for (i in seq_len(nrow(nybb_hits))) {

    df <- read_candidate(
      nybb_hits$path[[i]]
    )

    if (
      is.null(df) ||
      !(nybb_hits$id_column[[i]] %in% names(df))
    ) {
      found_nybb[[i]] <- FALSE
      next
    }

    found_nybb[[i]] <- nybb_id %in%
      as.character(
        df[[nybb_hits$id_column[[i]]]]
      )
  }

  nybb_tbl <- nybb_hits[
    found_nybb,
    ,
    drop = FALSE
  ]

  if (nrow(nybb_tbl) == 0) {

    cat(
      "No candidate file contained ",
      nybb_id,
      " in a recognized ID column.\n",
      sep = ""
    )

  } else {

    print(
      nybb_tbl |>
        select(
          path,
          id_column,
          regional_overlap,
          regional_overlap_pct,
          cerad_column,
          cerad_available_in_overlap,
          cerad_coverage_pct
        ),
      n = Inf
    )
  }
}


cat("\n============================================================\n")
cat("DISCOVERY COMPLETE — NO FILES WRITTEN\n")
cat("============================================================\n")

cat(
  "Do NOT rerun the production join yet.\n",
  "We need to select the actual AMP-AD Diverse Cohorts metadata\n",
  "source from the candidates above.\n",
  sep = ""
)
