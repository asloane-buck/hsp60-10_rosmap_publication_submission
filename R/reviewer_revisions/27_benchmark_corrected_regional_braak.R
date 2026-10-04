options(stringsAsFactors = FALSE)

base <- "/Users/ashlynsloane/Buck Institute 25-26/ROSMAP Analysis March 26/Brain Region Specificity"

outdir <- "outputs/reviewer_revisions/regional_braak_fix/benchmark"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

corrected_file <- file.path(
  "outputs/reviewer_revisions/regional_braak_fix",
  "regional_participant_metadata_corrected_braak.csv"
)

old_results_file <- file.path(
  base,
  "REGIONAL_COVARIATE_ADJUSTED_outputs",
  "AMP_covariate_adjusted_regional_screen_labeled_Hsp60_mito.csv"
)

regions <- list(
  DLPFC = list(
    map = file.path(
      base,
      "REGIONAL_FIRST_PASS_outputs",
      "AMP_DLPFC_sample_map_WITH_PHENO.csv"
    ),
    matrix = file.path(
      base,
      "Raw Data",
      "n1086_residual_log2_batch.csv"
    )
  ),
  STG = list(
    map = file.path(
      base,
      "REGIONAL_FIRST_PASS_outputs",
      "AMP_STG_sample_map_WITH_PHENO.csv"
    ),
    matrix = file.path(
      base,
      "Raw Data",
      "n278_residual_log2_batch.TCX.csv"
    )
  )
)

required_files <- c(
  corrected_file,
  old_results_file,
  unlist(lapply(regions, unlist))
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}


clean_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.0$", "", x)
}


get_accession <- function(x) {
  x <- trimws(as.character(x))

  ifelse(
    grepl("\\|", x),
    sub("^.*\\|", "", x),
    x
  )
}


matrix_accessions <- function(mat) {
  ids <- trimws(as.character(mat[[1]]))

  data.frame(
    matrix_row_id = ids,
    accession = get_accession(ids),
    stringsAsFactors = FALSE
  )
}


extract_protein <- function(mat, accession) {

  ids <- matrix_accessions(mat)

  hit <- which(ids$accession == accession)

  if (length(hit) != 1) {
    return(NULL)
  }

  values <- suppressWarnings(
    as.numeric(mat[hit, -1])
  )

  names(values) <- names(mat)[-1]

  values
}


fit_numeric <- function(d, predictor) {

  vars <- c(
    "abundance",
    predictor,
    "age_death_num",
    "sex",
    "pmi_num"
  )

  x <- d[
    complete.cases(d[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  if (nrow(x) < 20) {
    return(NULL)
  }

  if (length(unique(x[[predictor]])) < 2) {
    return(NULL)
  }

  formula <- as.formula(
    paste(
      "abundance ~",
      predictor,
      "+ age_death_num + sex + pmi_num"
    )
  )

  fit <- lm(formula, data = x)

  co <- summary(fit)$coefficients

  if (!(predictor %in% rownames(co))) {
    return(NULL)
  }

  ci <- confint(
    fit,
    predictor,
    level = 0.95
  )

  data.frame(
    n_rows = nrow(x),
    n_people = length(unique(x$individual_id)),
    beta = unname(co[predictor, "Estimate"]),
    se = unname(co[predictor, "Std. Error"]),
    statistic = unname(co[predictor, "t value"]),
    p = unname(co[predictor, "Pr(>|t|)"]),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2])
  )
}


fit_bin3 <- function(d) {

  # Convert to the intended categorical representation BEFORE
  # defining the complete-case analysis set.
  d$braak_bin3 <- trimws(as.character(d$braak_bin3))

  d$braak_bin3[
    d$braak_bin3 %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "<NA>"
    )
  ] <- NA_character_

  d$braak_bin3 <- factor(
    d$braak_bin3,
    levels = c(
      "I-II",
      "III-IV",
      "V-VI"
    )
  )

  vars <- c(
    "abundance",
    "braak_bin3",
    "age_death_num",
    "sex",
    "pmi_num"
  )

  x <- d[
    complete.cases(d[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  x$braak_bin3 <- droplevels(x$braak_bin3)

  if (
    nrow(x) < 20 ||
    nlevels(x$braak_bin3) < 2
  ) {
    return(NULL)
  }

  # Hard check: both nested models must use this exact same dataset.
  n_expected <- nrow(x)

  fit0 <- lm(
    abundance ~
      age_death_num +
      sex +
      pmi_num,
    data = x,
    na.action = na.fail
  )

  fit1 <- lm(
    abundance ~
      braak_bin3 +
      age_death_num +
      sex +
      pmi_num,
    data = x,
    na.action = na.fail
  )

  if (
    nobs(fit0) != n_expected ||
    nobs(fit1) != n_expected ||
    nobs(fit0) != nobs(fit1)
  ) {
    stop(
      "Three-bin models were not fitted to the same rows: ",
      "expected=", n_expected,
      ", reduced=", nobs(fit0),
      ", full=", nobs(fit1)
    )
  }

  comparison <- anova(fit0, fit1)
  co <- summary(fit1)$coefficients

  get_term <- function(term) {

    if (!(term %in% rownames(co))) {
      return(
        c(
          beta = NA_real_,
          se = NA_real_,
          p = NA_real_
        )
      )
    }

    c(
      beta = co[term, "Estimate"],
      se = co[term, "Std. Error"],
      p = co[term, "Pr(>|t|)"]
    )
  }

  mid <- get_term("braak_bin3III-IV")
  high <- get_term("braak_bin3V-VI")

  data.frame(
    n_rows = nrow(x),
    n_people = length(unique(x$individual_id)),
    global_p = comparison$`Pr(>F)`[2],
    beta_III_IV_vs_I_II = unname(mid["beta"]),
    se_III_IV_vs_I_II = unname(mid["se"]),
    p_III_IV_vs_I_II = unname(mid["p"]),
    beta_V_VI_vs_I_II = unname(high["beta"]),
    se_V_VI_vs_I_II = unname(high["se"]),
    p_V_VI_vs_I_II = unname(high["p"])
  )
}


cat("============================================================\n")
cat("CORRECTED REGIONAL BRAAK BENCHMARK\n")
cat("============================================================\n")

meta <- read.csv(
  corrected_file,
  check.names = FALSE
)

meta$individual_id <- clean_id(meta$individual_id)

old <- read.csv(
  old_results_file,
  check.names = FALSE
)


# ============================================================================
# VALIDATE OLD RESULT COLUMNS
# ============================================================================

required_old_columns <- c(
  "region_short",
  "gene_symbol",
  "gene_raw",
  "beta_braak",
  "p_braak",
  "is_hsp60_client"
)

missing_old <- setdiff(
  required_old_columns,
  names(old)
)

if (length(missing_old) > 0) {
  stop(
    "Old regional result table missing columns: ",
    paste(missing_old, collapse = ", ")
  )
}


# ============================================================================
# LOAD MATRICES ONCE
# ============================================================================

matrices <- list()

for (region in names(regions)) {

  cat("\nLoading ", region, " matrix...\n", sep = "")

  matrices[[region]] <- read.csv(
    regions[[region]]$matrix,
    check.names = FALSE
  )

  cat(
    region,
    " proteins:",
    nrow(matrices[[region]]),
    "\n"
  )
}


# ============================================================================
# IDENTIFY HSP CLIENT PROTEIN IDS PRESENT IN BOTH REGIONS AND BOTH MATRICES
# ============================================================================

client_rows <- old[
  old$is_hsp60_client %in% c(TRUE, "TRUE"),
  c(
    "region_short",
    "gene_symbol",
    "gene_raw"
  ),
  drop = FALSE
]

client_rows$accession <- get_accession(
  client_rows$gene_raw
)

dl_acc <- unique(
  matrix_accessions(
    matrices$DLPFC
  )$accession
)

st_acc <- unique(
  matrix_accessions(
    matrices$STG
  )$accession
)

candidate <- unique(
  client_rows[
    client_rows$accession %in% dl_acc &
      client_rows$accession %in% st_acc,
    c(
      "gene_symbol",
      "gene_raw",
      "accession"
    )
  ]
)

# Require exact old result row for this protein identifier in each region.
keep <- logical(nrow(candidate))

for (i in seq_len(nrow(candidate))) {

  gr <- candidate$gene_raw[i]

  n_dl <- sum(
    old$region_short == "DLPFC" &
      old$gene_raw == gr,
    na.rm = TRUE
  )

  n_st <- sum(
    old$region_short == "STG" &
      old$gene_raw == gr,
    na.rm = TRUE
  )

  keep[i] <- (
    n_dl == 1 &&
      n_st == 1
  )
}

candidate <- candidate[keep, , drop = FALSE]

candidate <- candidate[
  order(
    candidate$gene_symbol,
    candidate$accession
  ),
  ,
  drop = FALSE
]

benchmark <- head(
  candidate,
  8
)

if (nrow(benchmark) < 5) {
  stop(
    "Too few uniquely identifiable Hsp60/10 client proteins ",
    "are present in both regional matrices."
  )
}

cat("\nBenchmark proteins:\n")
print(
  benchmark,
  row.names = FALSE
)


# ============================================================================
# RUN BENCHMARK
# ============================================================================

all_results <- list()
k <- 1

for (region in names(regions)) {

  cat("\n============================================================\n")
  cat(region, "\n")
  cat("============================================================\n")

  smap <- read.csv(
    regions[[region]]$map,
    check.names = FALSE
  )

  smap$individual_id <- clean_id(
    smap$individual_id
  )

  meta_keep <- meta[
    ,
    c(
      "individual_id",
      "braak_num_legacy_2v6",
      "braak_stage_num",
      "braak_bin3"
    ),
    drop = FALSE
  ]

  dmap <- merge(
    smap,
    meta_keep,
    by = "individual_id",
    all.x = TRUE,
    sort = FALSE
  )

  mat <- matrices[[region]]

  eligible_legacy <- with(
    dmap,
    !is.na(braak_num_legacy_2v6) &
      !is.na(age_death_num) &
      !is.na(sex) &
      !is.na(pmi_num)
  )

  eligible_corrected <- with(
    dmap,
    !is.na(braak_stage_num) &
      !is.na(age_death_num) &
      !is.na(sex) &
      !is.na(pmi_num)
  )

  cat(
    "Legacy eligible rows:",
    sum(eligible_legacy),
    "\n"
  )

  cat(
    "Corrected eligible rows:",
    sum(eligible_corrected),
    "\n"
  )

  if (
    sum(eligible_legacy) !=
      sum(eligible_corrected)
  ) {
    stop(
      region,
      ": corrected Braak unexpectedly changed ",
      "phenotype eligibility."
    )
  }


  for (i in seq_len(nrow(benchmark))) {

    gene <- benchmark$gene_symbol[i]
    gene_raw <- benchmark$gene_raw[i]
    accession <- benchmark$accession[i]

    cat(
      "\nTesting ",
      region,
      " | ",
      gene,
      " | ",
      accession,
      "\n",
      sep = ""
    )

    old_row <- old[
      old$region_short == region &
        old$gene_raw == gene_raw,
      ,
      drop = FALSE
    ]

    if (nrow(old_row) != 1) {
      cat(
        "  SKIP: old result rows = ",
        nrow(old_row),
        "\n",
        sep = ""
      )
      next
    }

    y <- extract_protein(
      mat,
      accession
    )

    if (is.null(y)) {
      cat(
        "  SKIP: accession not uniquely found ",
        "in abundance matrix\n"
      )
      next
    }

    matching_rows <- (
      dmap$matrix_col %in%
        names(y)
    )

    d <- dmap[
      matching_rows,
      ,
      drop = FALSE
    ]

    d$abundance <- y[
      d$matrix_col
    ]

    cat(
      "  mapped assay rows: ",
      nrow(d),
      "\n",
      sep = ""
    )

    cat(
      "  nonmissing abundance: ",
      sum(!is.na(d$abundance)),
      "\n",
      sep = ""
    )

    legacy <- fit_numeric(
      d,
      "braak_num_legacy_2v6"
    )

    continuous <- fit_numeric(
      d,
      "braak_stage_num"
    )

    bin3 <- fit_bin3(d)

    if (is.null(legacy)) {
      cat("  SKIP: legacy model failed eligibility\n")
      next
    }

    if (is.null(continuous)) {
      cat("  SKIP: corrected continuous model failed eligibility\n")
      next
    }

    if (is.null(bin3)) {
      cat("  SKIP: three-bin model failed eligibility\n")
      next
    }

    cat("  Models completed.\n")

    all_results[[k]] <- data.frame(
      region = region,
      gene_symbol = gene,
      gene_raw = gene_raw,
      accession = accession,

      old_stored_beta =
        old_row$beta_braak[1],

      old_stored_p =
        old_row$p_braak[1],

      legacy_n_rows =
        legacy$n_rows,

      legacy_n_people =
        legacy$n_people,

      legacy_beta =
        legacy$beta,

      legacy_se =
        legacy$se,

      legacy_p =
        legacy$p,

      corrected_n_rows =
        continuous$n_rows,

      corrected_n_people =
        continuous$n_people,

      corrected_beta_per_stage =
        continuous$beta,

      corrected_se =
        continuous$se,

      corrected_p =
        continuous$p,

      corrected_ci_low =
        continuous$ci_low,

      corrected_ci_high =
        continuous$ci_high,

      bin3_global_p =
        bin3$global_p,

      bin3_beta_III_IV_vs_I_II =
        bin3$beta_III_IV_vs_I_II,

      bin3_p_III_IV_vs_I_II =
        bin3$p_III_IV_vs_I_II,

      bin3_beta_V_VI_vs_I_II =
        bin3$beta_V_VI_vs_I_II,

      bin3_p_V_VI_vs_I_II =
        bin3$p_V_VI_vs_I_II,

      stringsAsFactors = FALSE
    )

    k <- k + 1
  }
}


if (length(all_results) == 0) {
  stop(
    "No benchmark models completed. ",
    "See per-protein diagnostics above."
  )
}

results <- do.call(
  rbind,
  all_results
)

results$legacy_vs_stored_difference <- (
  results$legacy_beta -
    results$old_stored_beta
)

results$sign_agreement <- (
  sign(results$legacy_beta) ==
    sign(results$corrected_beta_per_stage)
)

outfile <- file.path(
  outdir,
  "corrected_braak_benchmark_results.csv"
)

write.csv(
  results,
  outfile,
  row.names = FALSE
)


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("BENCHMARK RESULTS\n")
cat("============================================================\n")

print(
  results[
    ,
    c(
      "region",
      "gene_symbol",
      "accession",
      "old_stored_beta",
      "legacy_beta",
      "corrected_beta_per_stage",
      "legacy_p",
      "corrected_p",
      "bin3_global_p",
      "sign_agreement"
    )
  ],
  row.names = FALSE
)


cat("\n============================================================\n")
cat("VALIDATION SUMMARY\n")
cat("============================================================\n")

cat(
  "Models completed:",
  nrow(results),
  "\n"
)

cat(
  "Sign agreement legacy vs corrected:",
  sum(results$sign_agreement),
  "/",
  nrow(results),
  "\n"
)

cat(
  "Max abs legacy-vs-stored beta difference:",
  max(
    abs(
      results$legacy_vs_stored_difference
    ),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Spearman legacy vs corrected beta:",
  cor(
    results$legacy_beta,
    results$corrected_beta_per_stage,
    method = "spearman",
    use = "complete.obs"
  ),
  "\n"
)

cat("\nWROTE:\n")
cat(outfile, "\n")

cat("\nCORRECTED BRAAK BENCHMARK COMPLETED.\n")
