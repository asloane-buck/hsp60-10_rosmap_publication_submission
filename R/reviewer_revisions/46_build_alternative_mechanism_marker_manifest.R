############################################################
## 46_build_alternative_mechanism_marker_manifest.R
##
## Reviewer 3.1
## Prespecified alternative-mechanism marker detectability.
##
## IMPORTANT:
## Marker list is defined before inspection of TMT results.
############################################################

options(stringsAsFactors = FALSE)

for (f in c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R"
)) {
  if (!file.exists(f)) {
    stop("Missing required production script: ", f)
  }
  source(f)
}

require_objects(
  c("prot_mat_raw", "prot_mat", "prot_meta_adj"),
  context = "46_build_alternative_mechanism_marker_manifest.R"
)

OUTDIR <- file.path(
  "outputs",
  "reviewer_revisions",
  "alternative_mechanism_marker_panel"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

marker_sets <- list(
  mitochondrial_mass_import = c(
    "TOMM20", "TOMM22", "TOMM40",
    "TIMM23", "VDAC1", "VDAC2", "CS"
  ),

  mitochondrial_biogenesis = c(
    "PPARGC1A", "PPARGC1B",
    "NRF1", "GABPA", "GABPB1",
    "ESRRA", "TFAM"
  ),

  mitophagy = c(
    "PINK1", "PRKN",
    "BNIP3", "BNIP3L",
    "FUNDC1", "OPTN",
    "CALCOCO2", "PHB2"
  ),

  autophagy_lysosome = c(
    "BECN1", "MAP1LC3B",
    "SQSTM1", "ATG5",
    "ATG7", "ULK1",
    "LAMP1", "LAMP2"
  ),

  proteasome_20S_core = c(
    paste0("PSMA", 1:7),
    paste0("PSMB", 1:7)
  )
)

## Canonicalize each marker set BEFORE constructing the
## category/gene table. clean_gene_symbols() sorts/uniques its
## input, so applying it to manifest$gene afterward would break
## category-to-gene row correspondence.
marker_sets <- lapply(
  marker_sets,
  clean_gene_symbols
)

manifest <- do.call(
  rbind,
  lapply(
    names(marker_sets),
    function(category) {
      data.frame(
        category = category,
        gene = marker_sets[[category]],
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(manifest) <- NULL

if (anyDuplicated(manifest[c("category", "gene")])) {
  stop("Duplicate category/gene entries in prespecified marker panel.")
}

expected_marker_categories <- c(
  TOMM20 = "mitochondrial_mass_import",
  TFAM = "mitochondrial_biogenesis",
  PINK1 = "mitophagy",
  ATG5 = "autophagy_lysosome",
  PSMA1 = "proteasome_20S_core"
)

for (g in names(expected_marker_categories)) {
  observed_category <- manifest$category[
    manifest$gene == g
  ]

  if (
    length(observed_category) != 1 ||
    observed_category != expected_marker_categories[[g]]
  ) {
    stop(
      "Marker/category mapping validation failed for ",
      g,
      ". Expected ",
      expected_marker_categories[[g]],
      "; observed ",
      paste(observed_category, collapse = ", ")
    )
  }
}

manifest$detected_raw_TMT <- (
  manifest$gene %in%
    rownames(prot_mat_raw)
)

manifest$detected_adjusted_TMT <- (
  manifest$gene %in%
    rownames(prot_mat)
)

manifest$n_observed_raw <- vapply(
  manifest$gene,
  function(g) {
    if (!g %in% rownames(prot_mat_raw)) {
      return(0L)
    }

    sum(
      is.finite(
        prot_mat_raw[g, ]
      )
    )
  },
  integer(1)
)

manifest$fraction_observed_raw <- (
  manifest$n_observed_raw /
    ncol(prot_mat_raw)
)

write.csv(
  manifest,
  file.path(
    OUTDIR,
    "prespecified_marker_detectability_manifest.csv"
  ),
  row.names = FALSE
)

category_summary <- manifest |>
  dplyr::group_by(category) |>
  dplyr::summarise(
    n_prespecified = dplyr::n(),
    n_detected_raw_TMT = sum(detected_raw_TMT),
    n_detected_adjusted_TMT = sum(detected_adjusted_TMT),
    pct_detected_raw_TMT =
      100 * n_detected_raw_TMT / n_prespecified,
    .groups = "drop"
  )

write.csv(
  category_summary,
  file.path(
    OUTDIR,
    "marker_detectability_by_category.csv"
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("PRESPECIFIED R3.1 MARKER DETECTABILITY\n")
cat("============================================================\n\n")

print(
  category_summary,
  row.names = FALSE
)

cat("\nMarker-level detectability:\n\n")

print(
  manifest[
    ,
    c(
      "category",
      "gene",
      "detected_raw_TMT",
      "n_observed_raw",
      "fraction_observed_raw"
    )
  ],
  row.names = FALSE
)

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)
