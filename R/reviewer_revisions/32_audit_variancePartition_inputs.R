options(stringsAsFactors = FALSE)

OUTDIR <- "outputs/reviewer_revisions/variancePartition_input_audit"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

required_scripts <- c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)

missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "Missing production scripts:\n",
    paste(missing_scripts, collapse = "\n")
  )
}

cat("======================================================================\n")
cat("VARIANCEPARTITION INPUT AUDIT\n")
cat("======================================================================\n")

cat("\nLoading finalized production objects...\n")

for (f in required_scripts) {
  cat("  sourcing:", f, "\n")
  source(f)
}

cat("\nProduction objects loaded.\n")


# =============================================================================
# Helpers
# =============================================================================

is_expression_object <- function(x) {

  if (is.matrix(x)) {

    return(
      is.numeric(x) &&
        nrow(x) >= 500 &&
        ncol(x) >= 50
    )
  }

  if (is.data.frame(x)) {

    if (
      nrow(x) < 500 ||
      ncol(x) < 50
    ) {
      return(FALSE)
    }

    numeric_fraction <- mean(
      vapply(
        x,
        is.numeric,
        logical(1)
      )
    )

    return(
      numeric_fraction >= 0.80
    )
  }

  FALSE
}


is_metadata_object <- function(x) {

  is.data.frame(x) &&
    nrow(x) >= 50 &&
    nrow(x) <= 5000 &&
    ncol(x) >= 2 &&
    ncol(x) <= 500
}


clean_values <- function(x) {

  y <- as.character(x)

  y[
    trimws(y) %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "<NA>"
    )
  ] <- NA

  y
}


candidate_covariate <- function(name) {

  grepl(
    paste(
      c(
        "age",
        "sex",
        "gender",
        "pmi",
        "rin",
        "batch",
        "tmt",
        "cogdx",
        "diagn",
        "stage",
        "braak",
        "cerad",
        "apoe",
        "education",
        "race",
        "ethnicity",
        "cohort",
        "study",
        "dataset"
      ),
      collapse = "|"
    ),
    name,
    ignore.case = TRUE
  )
}


candidate_id <- function(name) {

  grepl(
    "sample|individual|participant|subject|projid|specimen|assay|(^|_)id$",
    name,
    ignore.case = TRUE
  )
}


summarize_column <- function(
  object_name,
  column_name,
  x
) {

  y <- clean_values(x)

  nonmissing <- !is.na(y)

  n_unique <- length(
    unique(y[nonmissing])
  )

  is_num <- is.numeric(x)

  if (is_num) {

    xn <- suppressWarnings(
      as.numeric(x)
    )

    num_summary <- c(
      min = suppressWarnings(
        min(xn, na.rm = TRUE)
      ),
      median = suppressWarnings(
        median(xn, na.rm = TRUE)
      ),
      max = suppressWarnings(
        max(xn, na.rm = TRUE)
      )
    )

    if (!any(is.finite(num_summary))) {
      num_summary[] <- NA_real_
    }

  } else {

    num_summary <- c(
      min = NA_real_,
      median = NA_real_,
      max = NA_real_
    )
  }

  tab <- sort(
    table(
      y,
      useNA = "no"
    ),
    decreasing = TRUE
  )

  top_values <- if (
    length(tab) == 0
  ) {
    NA_character_
  } else {

    paste(
      paste0(
        names(head(tab, 8)),
        "=",
        as.integer(head(tab, 8))
      ),
      collapse = "; "
    )
  }

  data.frame(
    object = object_name,
    column = column_name,
    n_rows = length(x),
    n_nonmissing = sum(nonmissing),
    n_missing = sum(!nonmissing),
    fraction_nonmissing = mean(nonmissing),
    n_unique = n_unique,
    numeric = is_num,
    min = num_summary["min"],
    median = num_summary["median"],
    max = num_summary["max"],
    top_values = top_values,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Object inventory
# =============================================================================

object_names <- ls(
  envir = .GlobalEnv
)

inventory <- lapply(
  object_names,
  function(nm) {

    x <- get(
      nm,
      envir = .GlobalEnv
    )

    dims <- dim(x)

    data.frame(
      object = nm,
      class = paste(
        class(x),
        collapse = ";"
      ),
      nrow = if (
        is.null(dims)
      ) {
        NA_integer_
      } else {
        dims[1]
      },
      ncol = if (
        is.null(dims)
      ) {
        NA_integer_
      } else {
        dims[2]
      },
      expression_candidate =
        is_expression_object(x),
      metadata_candidate =
        is_metadata_object(x),
      name_suggests_raw =
        grepl(
          "raw|norm|normalized|count|abundance",
          nm,
          ignore.case = TRUE
        ),
      name_suggests_adjusted =
        grepl(
          "adj|adjust|resid",
          nm,
          ignore.case = TRUE
        ),
      stringsAsFactors = FALSE
    )
  }
)

inventory <- do.call(
  rbind,
  inventory
)

write.csv(
  inventory,
  file.path(
    OUTDIR,
    "production_object_inventory.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Expression candidates
# =============================================================================

expr_inventory <- inventory[
  inventory$expression_candidate,
  ,
  drop = FALSE
]

cat("\n======================================================================\n")
cat("EXPRESSION MATRIX CANDIDATES\n")
cat("======================================================================\n")

if (nrow(expr_inventory) == 0) {

  cat("No expression-matrix candidates detected.\n")

} else {

  print(
    expr_inventory[
      ,
      c(
        "object",
        "class",
        "nrow",
        "ncol",
        "name_suggests_raw",
        "name_suggests_adjusted"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# =============================================================================
# Metadata candidates
# =============================================================================

meta_inventory <- inventory[
  inventory$metadata_candidate,
  ,
  drop = FALSE
]

cat("\n======================================================================\n")
cat("METADATA TABLE CANDIDATES\n")
cat("======================================================================\n")

if (nrow(meta_inventory) == 0) {

  cat("No metadata candidates detected.\n")

} else {

  print(
    meta_inventory[
      ,
      c(
        "object",
        "nrow",
        "ncol"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# =============================================================================
# Covariate audit
# =============================================================================

covariate_rows <- list()
k <- 1

for (nm in meta_inventory$object) {

  x <- get(
    nm,
    envir = .GlobalEnv
  )

  cols <- names(x)[
    vapply(
      names(x),
      candidate_covariate,
      logical(1)
    )
  ]

  if (length(cols) == 0) {
    next
  }

  for (col in cols) {

    covariate_rows[[k]] <-
      summarize_column(
        nm,
        col,
        x[[col]]
      )

    k <- k + 1
  }
}

if (length(covariate_rows) > 0) {

  covariate_audit <- do.call(
    rbind,
    covariate_rows
  )

} else {

  covariate_audit <- data.frame()
}

write.csv(
  covariate_audit,
  file.path(
    OUTDIR,
    "metadata_covariate_audit.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("CANDIDATE COVARIATES\n")
cat("======================================================================\n")

if (nrow(covariate_audit) == 0) {

  cat("No candidate covariate fields detected.\n")

} else {

  print(
    covariate_audit[
      ,
      c(
        "object",
        "column",
        "n_rows",
        "n_nonmissing",
        "n_missing",
        "n_unique",
        "min",
        "median",
        "max",
        "top_values"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# =============================================================================
# Expression-metadata alignment audit
# =============================================================================

alignment_rows <- list()
k <- 1

for (expr_name in expr_inventory$object) {

  expr <- get(
    expr_name,
    envir = .GlobalEnv
  )

  if (is.null(colnames(expr))) {
    next
  }

  expr_ids <- clean_values(
    colnames(expr)
  )

  expr_ids <- unique(
    expr_ids[
      !is.na(expr_ids)
    ]
  )

  for (meta_name in meta_inventory$object) {

    meta <- get(
      meta_name,
      envir = .GlobalEnv
    )

    id_cols <- names(meta)[
      vapply(
        names(meta),
        candidate_id,
        logical(1)
      )
    ]

    if (length(id_cols) == 0) {
      next
    }

    for (id_col in id_cols) {

      vals <- clean_values(
        meta[[id_col]]
      )

      vals <- unique(
        vals[
          !is.na(vals)
        ]
      )

      overlap <- length(
        intersect(
          expr_ids,
          vals
        )
      )

      if (overlap == 0) {
        next
      }

      alignment_rows[[k]] <- data.frame(
        expression_object = expr_name,
        expression_n_columns =
          length(expr_ids),
        metadata_object = meta_name,
        metadata_id_column = id_col,
        metadata_n_unique_ids =
          length(vals),
        overlap = overlap,
        fraction_expression_matched =
          overlap /
          length(expr_ids),
        fraction_metadata_matched =
          overlap /
          length(vals),
        stringsAsFactors = FALSE
      )

      k <- k + 1
    }
  }
}

if (length(alignment_rows) > 0) {

  alignment <- do.call(
    rbind,
    alignment_rows
  )

  alignment <- alignment[
    order(
      alignment$expression_object,
      -alignment$overlap,
      -alignment$fraction_expression_matched
    ),
    ,
    drop = FALSE
  ]

} else {

  alignment <- data.frame()
}

write.csv(
  alignment,
  file.path(
    OUTDIR,
    "expression_metadata_alignment.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("BEST EXPRESSION-METADATA ALIGNMENTS\n")
cat("======================================================================\n")

if (nrow(alignment) == 0) {

  cat("No direct expression-column / metadata-ID matches detected.\n")

} else {

  best <- do.call(
    rbind,
    lapply(
      split(
        alignment,
        alignment$expression_object
      ),
      function(x) {

        x[
          order(
            -x$overlap,
            -x$fraction_expression_matched
          ),
          ,
          drop = FALSE
        ][
          seq_len(
            min(
              5,
              nrow(x)
            )
          ),
          ,
          drop = FALSE
        ]
      }
    )
  )

  rownames(best) <- NULL

  print(
    best,
    row.names = FALSE
  )
}


# =============================================================================
# Explicit raw-vs-adjusted object audit
# =============================================================================

cat("\n======================================================================\n")
cat("RAW / NORMALIZED VS ADJUSTED OBJECTS\n")
cat("======================================================================\n")

relevant_objects <- inventory[
  inventory$expression_candidate &
    (
      inventory$name_suggests_raw |
      inventory$name_suggests_adjusted
    ),
  ,
  drop = FALSE
]

if (nrow(relevant_objects) == 0) {

  cat(
    "No expression candidates had obvious raw/normalized/adjusted naming.\n"
  )

} else {

  print(
    relevant_objects[
      ,
      c(
        "object",
        "nrow",
        "ncol",
        "name_suggests_raw",
        "name_suggests_adjusted"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# =============================================================================
# Specifically inspect likely canonical production objects
# =============================================================================

canonical_patterns <- c(
  "rna",
  "prot",
  "protein",
  "meta",
  "clinical",
  "adjust",
  "raw",
  "batch"
)

canonical_hits <- inventory[
  Reduce(
    `|`,
    lapply(
      canonical_patterns,
      function(p) {
        grepl(
          p,
          inventory$object,
          ignore.case = TRUE
        )
      }
    )
  ),
  ,
  drop = FALSE
]

write.csv(
  canonical_hits,
  file.path(
    OUTDIR,
    "likely_canonical_objects.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("LIKELY CANONICAL RNA / PROTEIN OBJECTS\n")
cat("======================================================================\n")

print(
  canonical_hits[
    ,
    c(
      "object",
      "class",
      "nrow",
      "ncol",
      "expression_candidate",
      "metadata_candidate"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


# =============================================================================
# Package availability
# =============================================================================

packages <- c(
  "variancePartition",
  "limma",
  "edgeR",
  "BiocParallel"
)

package_status <- data.frame(
  package = packages,
  installed = vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  ),
  version = vapply(
    packages,
    function(pkg) {

      if (
        requireNamespace(
          pkg,
          quietly = TRUE
        )
      ) {
        as.character(
          packageVersion(pkg)
        )
      } else {
        NA_character_
      }
    },
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_status,
  file.path(
    OUTDIR,
    "variancePartition_package_status.csv"
  ),
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("PACKAGE STATUS\n")
cat("======================================================================\n")

print(
  package_status,
  row.names = FALSE
)


# =============================================================================
# Session / provenance
# =============================================================================

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
cat("AUDIT COMPLETE\n")
cat("======================================================================\n")

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)

cat(
  "\nIMPORTANT:\n",
  "variancePartition should be run on the appropriate non-residualized\n",
  "expression/abundance object. Adjusted/residualized matrices are being\n",
  "inventoried here only so that we can explicitly avoid using them.\n",
  sep = ""
)

cat("\nVARIANCEPARTITION INPUT AUDIT COMPLETED.\n")
