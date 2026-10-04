#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

read_script <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

s4 <- read_script(
  "R/supplemental/13_make_supplementary_figure_4_pathology_model_robustness.R"
)
s6 <- read_script(
  "R/supplemental/15_make_supplementary_figure_6_regional_proteomics_validation.R"
)

checks <- data.frame(
  check = c(
    "Supp Fig 4 Model 1 drop_na cleaned",
    "Supp Fig 4 Model 2 drop_na cleaned",
    "Supp Fig 6 scalar ad-associated-decline case_when removed",
    "Supp Fig 6 Wilcoxon group select cleaned"
  ),
  passed = c(
    !grepl(
      "drop_na(.data$y_z, .data$braak_z)",
      s4,
      fixed = TRUE
    ),
    !grepl(
      "drop_na(.data$y_z, .data$braak_z, .data$cerad_z)",
      s4,
      fixed = TRUE
    ),
    !grepl(
      'ad_col == "adjusted_collapse_magnitude" ~',
      s6,
      fixed = TRUE
    ),
    !grepl(
      "select(.data$group, value = dplyr::all_of(value_col))",
      s6,
      fixed = TRUE
    )
  )
)

print(checks, row.names = FALSE)

if (any(!checks$passed)) {
  stop("Residual warning-cleanup static validation FAILED.", call. = FALSE)
}

message("RESIDUAL WARNING-CLEANUP STATIC CHECKS PASSED.")
