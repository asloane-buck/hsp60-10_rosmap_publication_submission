options(stringsAsFactors = FALSE)

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop("Package 'nlme' is required.")
}

# =============================================================================
# Paths
# =============================================================================

BASE <- Sys.getenv(
  "ROSMAP_REGIONAL_BASE",
  unset = ""
)

if (!nzchar(BASE)) {
  stop(
    "Environment variable ROSMAP_REGIONAL_BASE is required. ",
    "Set it to the local Brain Region Specificity data directory."
  )
}

BASE <- normalizePath(
  BASE,
  mustWork = TRUE
)

FINAL_DIR <- paste0(
  "outputs/reviewer_revisions/regional_formal_interaction/",
  "validated_final"
)

OUTDIR <- paste0(
  "outputs/reviewer_revisions/regional_formal_interaction/",
  "hsp_conclusion_audit"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

DLPFC_MATRIX <- file.path(
  BASE,
  "Raw Data",
  "n1086_residual_log2_batch.csv"
)

STG_MATRIX <- file.path(
  BASE,
  "Raw Data",
  "n278_residual_log2_batch.TCX.csv"
)

DLPFC_MAP <- file.path(
  BASE,
  "REGIONAL_FIRST_PASS_outputs",
  "AMP_DLPFC_sample_map_WITH_PHENO.csv"
)

STG_MAP <- file.path(
  BASE,
  "REGIONAL_FIRST_PASS_outputs",
  "AMP_STG_sample_map_WITH_PHENO.csv"
)

SOURCE_META <- file.path(
  BASE,
  "Metadata",
  "AMP-AD_DiverseCohorts_individual_metadata.csv"
)

CORRECTED_META <- file.path(
  "outputs/reviewer_revisions/regional_braak_fix",
  "regional_participant_metadata_corrected_braak.csv"
)

OLD_HSP <- file.path(
  BASE,
  "REGIONAL_COVARIATE_ADJUSTED_outputs",
  "AMP_covariate_adjusted_Hsp60_clients_DLPFC_vs_STG.csv"
)

NEW_AD <- file.path(
  FINAL_DIR,
  "FINAL_region_x_AD_interaction_Hsp60_clients.csv"
)

NEW_BRAAK <- file.path(
  FINAL_DIR,
  "FINAL_region_x_Braak_continuous_interaction_Hsp60_clients.csv"
)

NEW_BIN3 <- file.path(
  FINAL_DIR,
  "FINAL_region_x_Braak_bin3_interaction_Hsp60_clients.csv"
)

required <- c(
  DLPFC_MATRIX,
  STG_MATRIX,
  DLPFC_MAP,
  STG_MAP,
  SOURCE_META,
  CORRECTED_META,
  OLD_HSP,
  NEW_AD,
  NEW_BRAAK,
  NEW_BIN3
)

missing <- required[
  !file.exists(required)
]

if (length(missing) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing, collapse = "\n")
  )
}


# =============================================================================
# Helpers
# =============================================================================

clean_id <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "<NA>"
    )
  ] <- NA

  sub(
    "\\.0$",
    "",
    x
  )
}


clean_missing <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "<NA>",
      "missing or unknown"
    )
  ] <- NA

  x
}


get_accession <- function(x) {

  x <- trimws(
    as.character(x)
  )

  ifelse(
    grepl("\\|", x),
    sub("^.*\\|", "", x),
    x
  )
}


zscore <- function(x) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  s <- sd(
    x,
    na.rm = TRUE
  )

  if (!is.finite(s) || s == 0) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  as.numeric(
    scale(x)
  )
}


read_expression <- function(file) {

  x <- read.csv(
    file,
    check.names = FALSE
  )

  accession <- get_accession(
    x[[1]]
  )

  if (anyDuplicated(accession)) {
    stop(
      "Duplicate protein accessions in: ",
      file
    )
  }

  mat <- as.matrix(
    x[, -1, drop = FALSE]
  )

  storage.mode(mat) <- "numeric"

  rownames(mat) <- accession

  mat
}


collapse_region <- function(
  mat,
  smap,
  paired_ids
) {

  smap$individual_id <- clean_id(
    smap$individual_id
  )

  sm <- smap[
    !is.na(smap$individual_id) &
      smap$individual_id %in% paired_ids &
      smap$matrix_col %in% colnames(mat),
    ,
    drop = FALSE
  ]

  if (anyDuplicated(sm$matrix_col)) {
    stop(
      "Duplicate matrix_col values in sample map."
    )
  }

  by_person <- split(
    sm$matrix_col,
    sm$individual_id
  )

  out <- vapply(
    by_person,
    function(cols) {

      x <- rowMeans(
        mat[, cols, drop = FALSE],
        na.rm = TRUE
      )

      x[is.nan(x)] <- NA_real_

      x
    },
    FUN.VALUE = numeric(
      nrow(mat)
    )
  )

  rownames(out) <- rownames(mat)

  out
}


# =============================================================================
# Mixed-model fitting with optimizer fallback
# =============================================================================

fit_lme_fallback <- function(
  fixed,
  data,
  method = "REML"
) {

  attempt <- function(
    optimizer,
    max_iter,
    max_eval
  ) {

    warnings_seen <- character()

    fit <- withCallingHandlers(

      tryCatch(
        nlme::lme(
          fixed = fixed,
          random = ~1 | individual_id,
          data = data,
          method = method,
          na.action = na.fail,
          control = nlme::lmeControl(
            opt = optimizer,
            msMaxIter = max_iter,
            msMaxEval = max_eval,
            returnObject = TRUE
          )
        ),
        error = function(e) e
      ),

      warning = function(w) {

        warnings_seen <<- c(
          warnings_seen,
          conditionMessage(w)
        )

        invokeRestart(
          "muffleWarning"
        )
      }
    )

    list(
      fit = fit,
      warnings = unique(
        warnings_seen
      )
    )
  }


  first <- attempt(
    "optim",
    200,
    400
  )

  if (
    !inherits(first$fit, "error") &&
    length(first$warnings) == 0
  ) {

    return(
      list(
        fit = first$fit,
        optimizer = "optim",
        fallback = FALSE
      )
    )
  }


  second <- attempt(
    "nlminb",
    1000,
    2000
  )

  if (
    !inherits(second$fit, "error") &&
    length(second$warnings) == 0
  ) {

    return(
      list(
        fit = second$fit,
        optimizer = "nlminb",
        fallback = TRUE
      )
    )
  }


  first_msg <- if (
    inherits(first$fit, "error")
  ) {
    conditionMessage(first$fit)
  } else {
    paste(
      first$warnings,
      collapse = " | "
    )
  }

  second_msg <- if (
    inherits(second$fit, "error")
  ) {
    conditionMessage(second$fit)
  } else {
    paste(
      second$warnings,
      collapse = " | "
    )
  }

  stop(
    "Both optimizers failed.\n",
    "optim: ",
    first_msg,
    "\n",
    "nlminb: ",
    second_msg
  )
}


find_interaction <- function(
  terms,
  predictor
) {

  candidates <- c(
    paste0(
      "regionSTG:",
      predictor
    ),
    paste0(
      predictor,
      ":regionSTG"
    )
  )

  hit <- candidates[
    candidates %in% terms
  ]

  if (length(hit) != 1) {
    stop(
      "Could not uniquely identify interaction for ",
      predictor
    )
  }

  hit[1]
}


# =============================================================================
# Fit one pathway-score numeric interaction
# =============================================================================

fit_score_numeric <- function(
  score_long,
  predictor,
  label
) {

  needed <- c(
    "score",
    predictor,
    "region",
    "cohort",
    "age_z",
    "sex",
    "pmi_z",
    "individual_id"
  )

  d <- score_long[
    complete.cases(
      score_long[
        ,
        needed,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  counts <- table(
    d$individual_id
  )

  paired <- names(
    counts[counts == 2]
  )

  d <- d[
    d$individual_id %in% paired,
    ,
    drop = FALSE
  ]

  d$individual_id <- factor(
    d$individual_id
  )

  d$region <- factor(
    d$region,
    levels = c(
      "DLPFC",
      "STG"
    )
  )

  d$cohort <- droplevels(
    factor(
      d$cohort,
      levels = c(
        "Mayo Clinic",
        "Emory"
      )
    )
  )

  d$sex <- droplevels(
    factor(d$sex)
  )

  f <- as.formula(
    paste0(
      "score ~ region * ",
      predictor,
      " + region * cohort",
      " + region * age_z",
      " + region * sex",
      " + region * pmi_z"
    )
  )

  ans <- fit_lme_fallback(
    fixed = f,
    data = d,
    method = "REML"
  )

  tt <- summary(
    ans$fit
  )$tTable

  term <- find_interaction(
    rownames(tt),
    predictor
  )

  beta <- unname(
    tt[
      term,
      "Value"
    ]
  )

  se <- unname(
    tt[
      term,
      "Std.Error"
    ]
  )

  df <- unname(
    tt[
      term,
      "DF"
    ]
  )

  p <- unname(
    tt[
      term,
      "p-value"
    ]
  )

  crit <- qt(
    0.975,
    df = df
  )

  beta_dlpfc <- unname(
    tt[
      predictor,
      "Value"
    ]
  )

  beta_stg <- (
    beta_dlpfc +
    beta
  )

  data.frame(
    analysis = label,

    n_people = length(
      unique(
        d$individual_id
      )
    ),

    n_observations = nrow(d),

    optimizer = ans$optimizer,
    optimizer_fallback = ans$fallback,

    beta_region_x_predictor = beta,
    se_region_x_predictor = se,
    df_region_x_predictor = df,

    ci_low_region_x_predictor =
      beta - crit * se,

    ci_high_region_x_predictor =
      beta + crit * se,

    p_region_x_predictor = p,

    beta_predictor_DLPFC =
      beta_dlpfc,

    beta_predictor_STG =
      beta_stg,

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Categorical Braak score sensitivity
# =============================================================================

fit_score_bin3 <- function(
  score_long,
  label
) {

  d <- score_long

  d$braak_bin3 <- factor(
    d$braak_bin3,
    levels = c(
      "I-II",
      "III-IV",
      "V-VI"
    )
  )

  needed <- c(
    "score",
    "braak_bin3",
    "region",
    "cohort",
    "age_z",
    "sex",
    "pmi_z",
    "individual_id"
  )

  d <- d[
    complete.cases(
      d[
        ,
        needed,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  counts <- table(
    d$individual_id
  )

  paired <- names(
    counts[counts == 2]
  )

  d <- d[
    d$individual_id %in% paired,
    ,
    drop = FALSE
  ]

  d$individual_id <- factor(
    d$individual_id
  )

  d$region <- factor(
    d$region,
    levels = c(
      "DLPFC",
      "STG"
    )
  )

  d$cohort <- droplevels(
    factor(
      d$cohort,
      levels = c(
        "Mayo Clinic",
        "Emory"
      )
    )
  )

  d$sex <- droplevels(
    factor(d$sex)
  )

  d$braak_bin3 <- droplevels(
    d$braak_bin3
  )

  full_formula <- score ~
    region * braak_bin3 +
    region * cohort +
    region * age_z +
    region * sex +
    region * pmi_z

  reduced_formula <- score ~
    region +
    braak_bin3 +
    region * cohort +
    region * age_z +
    region * sex +
    region * pmi_z

  full <- fit_lme_fallback(
    fixed = full_formula,
    data = d,
    method = "ML"
  )

  reduced <- fit_lme_fallback(
    fixed = reduced_formula,
    data = d,
    method = "ML"
  )

  cmp <- anova(
    reduced$fit,
    full$fit
  )

  data.frame(
    analysis = label,

    n_people = length(
      unique(
        d$individual_id
      )
    ),

    n_observations = nrow(d),

    full_optimizer =
      full$optimizer,

    reduced_optimizer =
      reduced$optimizer,

    global_region_x_braak_bin3_p =
      cmp$`p-value`[2],

    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Existing formal interaction results
# =============================================================================

cat(
  "======================================================================\n"
)

cat(
  "HSP60/10 REGIONAL CONCLUSION AUDIT\n"
)

cat(
  "======================================================================\n"
)

ad_new <- read.csv(
  NEW_AD,
  check.names = FALSE
)

braak_new <- read.csv(
  NEW_BRAAK,
  check.names = FALSE
)

bin3_new <- read.csv(
  NEW_BIN3,
  check.names = FALSE
)

if (
  nrow(ad_new) != 300 ||
  nrow(braak_new) != 300 ||
  nrow(bin3_new) != 300
) {
  stop(
    "Expected exactly 300 Hsp clients in each validated table."
  )
}


# =============================================================================
# Descriptive Hsp interaction audit
# =============================================================================

describe_interactions <- function(
  x,
  label
) {

  data.frame(
    analysis = label,

    n = nrow(x),

    median_interaction =
      median(
        x$beta_region_x_predictor,
        na.rm = TRUE
      ),

    mean_interaction =
      mean(
        x$beta_region_x_predictor,
        na.rm = TRUE
      ),

    n_negative =
      sum(
        x$beta_region_x_predictor < 0,
        na.rm = TRUE
      ),

    n_positive =
      sum(
        x$beta_region_x_predictor > 0,
        na.rm = TRUE
      ),

    fraction_negative =
      mean(
        x$beta_region_x_predictor < 0,
        na.rm = TRUE
      ),

    n_genomewide_FDR_lt_0.05 =
      sum(
        x$p_region_x_predictor_fdr < 0.05,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )
}


interaction_summary <- rbind(
  describe_interactions(
    ad_new,
    "AD"
  ),

  describe_interactions(
    braak_new,
    "Braak_continuous"
  )
)

write.csv(
  interaction_summary,
  file.path(
    OUTDIR,
    "Hsp_formal_interaction_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nFormal Hsp interaction summary:\n"
)

print(
  interaction_summary,
  row.names = FALSE
)


# =============================================================================
# Significant Hsp proteins
# =============================================================================

ad_sig <- ad_new[
  !is.na(
    ad_new$p_region_x_predictor_fdr
  ) &
    ad_new$p_region_x_predictor_fdr < 0.05,
  ,
  drop = FALSE
]

braak_sig <- braak_new[
  !is.na(
    braak_new$p_region_x_predictor_fdr
  ) &
    braak_new$p_region_x_predictor_fdr < 0.05,
  ,
  drop = FALSE
]

write.csv(
  ad_sig,
  file.path(
    OUTDIR,
    "Hsp_AD_formal_interaction_FDR_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  braak_sig,
  file.path(
    OUTDIR,
    "Hsp_Braak_formal_interaction_FDR_significant.csv"
  ),
  row.names = FALSE
)

show_cols <- c(
  "accession",
  "gene_symbol",
  "n_people",
  "beta_region_x_predictor",
  "ci_low_region_x_predictor",
  "ci_high_region_x_predictor",
  "p_region_x_predictor",
  "p_region_x_predictor_fdr",
  "beta_predictor_DLPFC",
  "beta_predictor_STG"
)

show_cols <- show_cols[
  show_cols %in% names(ad_new)
]

cat(
  "\nGenome-wide FDR-significant Hsp region x AD interactions:\n"
)

print(
  ad_sig[
    ,
    show_cols,
    drop = FALSE
  ],
  row.names = FALSE
)

cat(
  "\nGenome-wide FDR-significant Hsp region x Braak interactions:\n"
)

print(
  braak_sig[
    ,
    show_cols,
    drop = FALSE
  ],
  row.names = FALSE
)


# =============================================================================
# Compare old regional subtraction metrics to new formal interactions
# =============================================================================

old <- read.csv(
  OLD_HSP,
  check.names = FALSE
)

if (!"gene_raw" %in% names(old)) {
  stop(
    "Old regional Hsp table lacks gene_raw."
  )
}

old$accession <- get_accession(
  old$gene_raw
)

required_old <- c(
  "accession",
  "stg_minus_dlpfc_adjusted_collapse",
  "stg_minus_dlpfc_adjusted_inverse_braak"
)

missing_old <- setdiff(
  required_old,
  names(old)
)

if (length(missing_old) > 0) {
  stop(
    "Old Hsp table lacks expected columns: ",
    paste(
      missing_old,
      collapse = ", "
    )
  )
}

old_compare <- old[
  ,
  required_old,
  drop = FALSE
]

ad_compare <- merge(
  ad_new[
    ,
    c(
      "accession",
      "gene_symbol",
      "beta_region_x_predictor"
    ),
    drop = FALSE
  ],
  old_compare[
    ,
    c(
      "accession",
      "stg_minus_dlpfc_adjusted_collapse"
    ),
    drop = FALSE
  ],
  by = "accession",
  all = FALSE
)

names(ad_compare)[
  names(ad_compare) ==
    "beta_region_x_predictor"
] <- "new_formal_interaction"

braak_compare <- merge(
  braak_new[
    ,
    c(
      "accession",
      "gene_symbol",
      "beta_region_x_predictor"
    ),
    drop = FALSE
  ],
  old_compare[
    ,
    c(
      "accession",
      "stg_minus_dlpfc_adjusted_inverse_braak"
    ),
    drop = FALSE
  ],
  by = "accession",
  all = FALSE
)

names(braak_compare)[
  names(braak_compare) ==
    "beta_region_x_predictor"
] <- "new_formal_interaction"


compare_metrics <- function(
  new,
  old
) {

  ok <- (
    is.finite(new) &
    is.finite(old)
  )

  new <- new[ok]
  old <- old[ok]

  data.frame(
    n = length(new),

    spearman_old_vs_new =
      cor(
        old,
        new,
        method = "spearman"
      ),

    spearman_old_vs_negative_new =
      cor(
        old,
        -new,
        method = "spearman"
      ),

    sign_agreement_old_vs_new =
      mean(
        sign(old) ==
          sign(new)
      ),

    sign_agreement_old_vs_negative_new =
      mean(
        sign(old) ==
          sign(-new)
      ),

    stringsAsFactors = FALSE
  )
}


old_new_summary <- rbind(

  cbind(
    analysis = "old_collapse_vs_formal_AD_interaction",

    compare_metrics(
      ad_compare$new_formal_interaction,
      ad_compare$stg_minus_dlpfc_adjusted_collapse
    )
  ),

  cbind(
    analysis = "old_inverse_Braak_vs_formal_Braak_interaction",

    compare_metrics(
      braak_compare$new_formal_interaction,
      braak_compare$stg_minus_dlpfc_adjusted_inverse_braak
    )
  )
)

write.csv(
  old_new_summary,
  file.path(
    OUTDIR,
    "old_vs_formal_regional_Hsp_comparison_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  ad_compare,
  file.path(
    OUTDIR,
    "old_vs_formal_AD_Hsp_gene_level.csv"
  ),
  row.names = FALSE
)

write.csv(
  braak_compare,
  file.path(
    OUTDIR,
    "old_vs_formal_Braak_Hsp_gene_level.csv"
  ),
  row.names = FALSE
)

cat(
  "\nOld subtraction metrics vs new formal interactions:\n"
)

print(
  old_new_summary,
  row.names = FALSE
)


# =============================================================================
# Rebuild paired regional expression matrices
# =============================================================================

cat(
  "\nRebuilding paired regional matrices for pathway-score analysis...\n"
)

dl_map <- read.csv(
  DLPFC_MAP,
  check.names = FALSE
)

st_map <- read.csv(
  STG_MAP,
  check.names = FALSE
)

dl_map$individual_id <- clean_id(
  dl_map$individual_id
)

st_map$individual_id <- clean_id(
  st_map$individual_id
)

paired_ids <- sort(
  intersect(
    unique(
      dl_map$individual_id[
        !is.na(
          dl_map$individual_id
        )
      ]
    ),
    unique(
      st_map$individual_id[
        !is.na(
          st_map$individual_id
        )
      ]
    )
  )
)

if (length(paired_ids) != 215) {
  stop(
    "Expected 215 paired participants; found ",
    length(paired_ids)
  )
}

DL_raw <- read_expression(
  DLPFC_MATRIX
)

ST_raw <- read_expression(
  STG_MATRIX
)

DL <- collapse_region(
  DL_raw,
  dl_map,
  paired_ids
)

ST <- collapse_region(
  ST_raw,
  st_map,
  paired_ids
)

hsp_accessions <- sort(
  unique(
    ad_new$accession
  )
)

hsp_accessions <- intersect(
  hsp_accessions,
  intersect(
    rownames(DL),
    rownames(ST)
  )
)

if (length(hsp_accessions) != 300) {
  stop(
    "Expected 300 shared regional Hsp clients; found ",
    length(hsp_accessions)
  )
}

DL <- DL[
  hsp_accessions,
  paired_ids,
  drop = FALSE
]

ST <- ST[
  hsp_accessions,
  paired_ids,
  drop = FALSE
]


# =============================================================================
# Hsp client coverage
# =============================================================================

coverage <- data.frame(
  accession = hsp_accessions,

  DLPFC_nonmissing =
    rowSums(
      !is.na(DL)
    ),

  STG_nonmissing =
    rowSums(
      !is.na(ST)
    ),

  DLPFC_fraction =
    rowMeans(
      !is.na(DL)
    ),

  STG_fraction =
    rowMeans(
      !is.na(ST)
    ),

  stringsAsFactors = FALSE
)

coverage$min_region_fraction <- pmin(
  coverage$DLPFC_fraction,
  coverage$STG_fraction
)

write.csv(
  coverage,
  file.path(
    OUTDIR,
    "Hsp_pathway_score_client_coverage.csv"
  ),
  row.names = FALSE
)

cat(
  "\nHsp client coverage across 215 paired participants:\n"
)

cat(
  "  Hsp clients:",
  nrow(coverage),
  "\n"
)

cat(
  "  Complete in both regions for all 215:",
  sum(
    coverage$DLPFC_nonmissing == 215 &
      coverage$STG_nonmissing == 215
  ),
  "\n"
)

cat(
  "  >=90% coverage in each region:",
  sum(
    coverage$min_region_fraction >= 0.90
  ),
  "\n"
)

cat(
  "  >=80% coverage in each region:",
  sum(
    coverage$min_region_fraction >= 0.80
  ),
  "\n"
)


# =============================================================================
# Primary pathway score
#
# Include clients with >=90% participant coverage in BOTH regions.
# Each protein is z-scored jointly across DLPFC + STG observations so that:
#   - each protein contributes equally;
#   - both regions remain on the same protein-specific scale;
#   - regional predictor interactions are preserved.
#
# No imputation.
# Each participant-region score requires >=90% of selected clients observed.
# =============================================================================

score_clients <- coverage$accession[
  coverage$min_region_fraction >= 0.90
]

if (length(score_clients) < 100) {
  stop(
    "Fewer than 100 Hsp clients meet >=90% coverage; inspect before scoring."
  )
}

DL_score_z <- matrix(
  NA_real_,
  nrow = length(score_clients),
  ncol = length(paired_ids),
  dimnames = list(
    score_clients,
    paired_ids
  )
)

ST_score_z <- DL_score_z

for (acc in score_clients) {

  combined <- c(
    DL[acc, paired_ids],
    ST[acc, paired_ids]
  )

  mu <- mean(
    combined,
    na.rm = TRUE
  )

  sig <- sd(
    combined,
    na.rm = TRUE
  )

  if (
    !is.finite(sig) ||
    sig == 0
  ) {
    next
  }

  DL_score_z[acc, ] <- (
    DL[acc, paired_ids] -
      mu
  ) / sig

  ST_score_z[acc, ] <- (
    ST[acc, paired_ids] -
      mu
  ) / sig
}

min_required <- ceiling(
  0.90 *
    length(score_clients)
)

dl_n <- colSums(
  !is.na(
    DL_score_z
  )
)

st_n <- colSums(
  !is.na(
    ST_score_z
  )
)

dl_score <- colMeans(
  DL_score_z,
  na.rm = TRUE
)

st_score <- colMeans(
  ST_score_z,
  na.rm = TRUE
)

dl_score[
  dl_n < min_required
] <- NA_real_

st_score[
  st_n < min_required
] <- NA_real_

cat(
  "\nPrimary pathway-score construction:\n"
)

cat(
  "  Selected clients:",
  length(score_clients),
  "\n"
)

cat(
  "  Required observed clients/sample:",
  min_required,
  "\n"
)

cat(
  "  DLPFC samples passing coverage:",
  sum(
    !is.na(dl_score)
  ),
  "/ 215\n"
)

cat(
  "  STG samples passing coverage:",
  sum(
    !is.na(st_score)
  ),
  "/ 215\n"
)


# =============================================================================
# Phenotypes
# =============================================================================

source <- read.csv(
  SOURCE_META,
  check.names = FALSE
)

source$individual_id <- clean_id(
  source$individualID
)

corrected <- read.csv(
  CORRECTED_META,
  check.names = FALSE
)

corrected$individual_id <- clean_id(
  corrected$individual_id
)

cohort_lookup <- source[
  ,
  c(
    "individual_id",
    "cohort"
  ),
  drop = FALSE
]

pheno <- merge(
  corrected,
  cohort_lookup,
  by = "individual_id",
  all.x = TRUE,
  validate = "one_to_one"
)

pheno <- pheno[
  pheno$individual_id %in% paired_ids,
  ,
  drop = FALSE
]

pheno$age_death_num <- suppressWarnings(
  as.numeric(
    pheno$age_death_num
  )
)

pheno$pmi_num <- suppressWarnings(
  as.numeric(
    pheno$pmi_num
  )
)

pheno$braak_stage_num <- suppressWarnings(
  as.numeric(
    pheno$braak_stage_num
  )
)

pheno$age_z <- zscore(
  pheno$age_death_num
)

pheno$pmi_z <- zscore(
  pheno$pmi_num
)

pheno$cohort <- clean_missing(
  pheno$cohort
)

pheno$sex <- clean_missing(
  pheno$sex
)

pheno$ad_outcome <- clean_missing(
  pheno$ad_outcome
)

pheno$braak_bin3 <- clean_missing(
  pheno$braak_bin3
)

pheno$ad_status <- NA_real_

pheno$ad_status[
  pheno$ad_outcome == "Control"
] <- 0

pheno$ad_status[
  pheno$ad_outcome == "AD"
] <- 1


# =============================================================================
# Long score table
# =============================================================================

make_score_rows <- function(
  region_name,
  scores,
  n_clients
) {

  d <- pheno

  d$region <- region_name

  d$score <- scores[
    d$individual_id
  ]

  d$n_Hsp_clients_observed <- n_clients[
    d$individual_id
  ]

  d
}

score_long <- rbind(
  make_score_rows(
    "DLPFC",
    dl_score,
    dl_n
  ),
  make_score_rows(
    "STG",
    st_score,
    st_n
  )
)

score_long$region <- factor(
  score_long$region,
  levels = c(
    "DLPFC",
    "STG"
  )
)

write.csv(
  score_long,
  file.path(
    OUTDIR,
    "Hsp_pathway_score_participant_region_table.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Primary pathway models
# =============================================================================

ad_long <- score_long[
  score_long$ad_outcome %in%
    c(
      "Control",
      "AD"
    ),
  ,
  drop = FALSE
]

braak_long <- score_long[
  score_long$braak_stage_num %in%
    1:6,
  ,
  drop = FALSE
]

pathway_ad <- fit_score_numeric(
  ad_long,
  predictor = "ad_status",
  label = "Hsp_score_region_x_AD"
)

pathway_braak <- fit_score_numeric(
  braak_long,
  predictor = "braak_stage_num",
  label = "Hsp_score_region_x_Braak_continuous"
)

pathway_bin3 <- fit_score_bin3(
  braak_long,
  label = "Hsp_score_region_x_Braak_bin3"
)

pathway_primary <- rbind(
  pathway_ad,
  pathway_braak
)

pathway_primary$BH_FDR_two_primary_tests <- p.adjust(
  pathway_primary$p_region_x_predictor,
  method = "BH"
)

write.csv(
  pathway_primary,
  file.path(
    OUTDIR,
    "Hsp_pathway_score_primary_interactions.csv"
  ),
  row.names = FALSE
)

write.csv(
  pathway_bin3,
  file.path(
    OUTDIR,
    "Hsp_pathway_score_Braak_bin3_sensitivity.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Complete-client sensitivity
# =============================================================================

complete_clients <- coverage$accession[
  coverage$DLPFC_nonmissing == 215 &
    coverage$STG_nonmissing == 215
]

complete_sensitivity <- NULL

if (length(complete_clients) >= 50) {

  DL_complete_z <- matrix(
    NA_real_,
    nrow = length(complete_clients),
    ncol = length(paired_ids),
    dimnames = list(
      complete_clients,
      paired_ids
    )
  )

  ST_complete_z <- DL_complete_z

  for (acc in complete_clients) {

    combined <- c(
      DL[acc, paired_ids],
      ST[acc, paired_ids]
    )

    mu <- mean(combined)

    sig <- sd(combined)

    DL_complete_z[acc, ] <- (
      DL[acc, paired_ids] -
        mu
    ) / sig

    ST_complete_z[acc, ] <- (
      ST[acc, paired_ids] -
        mu
    ) / sig
  }

  dl_complete_score <- colMeans(
    DL_complete_z
  )

  st_complete_score <- colMeans(
    ST_complete_z
  )

  complete_long <- rbind(
    make_score_rows(
      "DLPFC",
      dl_complete_score,
      rep(
        length(complete_clients),
        length(paired_ids)
      )
    ),
    make_score_rows(
      "STG",
      st_complete_score,
      rep(
        length(complete_clients),
        length(paired_ids)
      )
    )
  )

  complete_long$region <- factor(
    complete_long$region,
    levels = c(
      "DLPFC",
      "STG"
    )
  )

  complete_ad <- complete_long[
    complete_long$ad_outcome %in%
      c(
        "Control",
        "AD"
      ),
    ,
    drop = FALSE
  ]

  complete_braak <- complete_long[
    complete_long$braak_stage_num %in%
      1:6,
    ,
    drop = FALSE
  ]

  complete_sensitivity <- rbind(

    fit_score_numeric(
      complete_ad,
      predictor = "ad_status",
      label = "complete_client_Hsp_score_region_x_AD"
    ),

    fit_score_numeric(
      complete_braak,
      predictor = "braak_stage_num",
      label = "complete_client_Hsp_score_region_x_Braak"
    )
  )

  write.csv(
    complete_sensitivity,
    file.path(
      OUTDIR,
      "Hsp_pathway_score_complete_client_sensitivity.csv"
    ),
    row.names = FALSE
  )
}


# =============================================================================
# Report
# =============================================================================

cat(
  "\n======================================================================\n"
)

cat(
  "HSP PATHWAY-LEVEL FORMAL INTERACTION RESULTS\n"
)

cat(
  "======================================================================\n"
)

print(
  pathway_primary,
  row.names = FALSE
)

cat(
  "\nThree-bin Braak pathway sensitivity:\n"
)

print(
  pathway_bin3,
  row.names = FALSE
)

cat(
  "\nComplete-client sensitivity proteins:",
  length(complete_clients),
  "\n"
)

if (!is.null(complete_sensitivity)) {

  print(
    complete_sensitivity,
    row.names = FALSE
  )

} else {

  cat(
    "Complete-client sensitivity skipped because <50 clients were complete.\n"
  )
}


cat(
  "\n======================================================================\n"
)

cat(
  "INTERPRETATION KEY\n"
)

cat(
  "======================================================================\n"
)

cat(
  "For AD:\n"
)

cat(
  "  interaction beta < 0 = AD-associated Hsp score change is more negative\n"
)

cat(
  "                         in STG than DLPFC.\n"
)

cat(
  "  interaction beta > 0 = AD-associated Hsp score change is more positive\n"
)

cat(
  "                         in STG than DLPFC.\n"
)

cat(
  "\nFor continuous Braak:\n"
)

cat(
  "  interaction beta < 0 = Hsp score declines more strongly with increasing\n"
)

cat(
  "                         Braak stage in STG than DLPFC.\n"
)

cat(
  "  interaction beta > 0 = Braak association is more positive / less negative\n"
)

cat(
  "                         in STG than DLPFC.\n"
)

cat(
  "\nThese are direct region x predictor interaction tests from one combined\n"
)

cat(
  "paired mixed model, not subtraction of separately fitted regional effects.\n"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)

cat(
  "\nHSP REGIONAL CONCLUSION AUDIT COMPLETED.\n"
)
