#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

source("R/00_config.R")
source("R/01_utils.R")

read_script <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

r13 <- read_script("R/13_make_main_figure_4_cognition.R")
r14 <- read_script("R/14_make_main_figure_5_matched_null_specificity.R")
s3 <- read_script("R/supplemental/12_make_supplementary_figure_3_mitochondrial_specificity.R")

checks <- tibble::tribble(
  ~check, ~passed,
  "Figure 4 MitoCarta import forced to text",
  grepl('col_types = "text"', r13, fixed = TRUE),
  "Figure 5 obsolete label.size removed",
  !grepl("label.size =", r14, fixed = TRUE),
  "Figure 5 cognition hard scale limit removed",
  !grepl("limits = c(xmin_cognition, xmax_cognition)", r14, fixed = TRUE),
  "Figure 5 AGORA hard scale limit removed",
  !grepl("limits = c(xmin_agora, xmax_agora)", r14, fixed = TRUE),
  "Supp Fig 3 scalar case_when removed",
  !grepl("is.logical(.data$agora_nominated_target) ~", s3, fixed = TRUE)
)

print(checks, n = Inf)
if (any(!checks$passed)) {
  stop("Warning-cleanup static validation FAILED.", call. = FALSE)
}
message("WARNING CLEANUP STATIC CHECKS PASSED.")
