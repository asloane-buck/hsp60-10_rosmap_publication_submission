options(stringsAsFactors = FALSE)

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop("Package 'nlme' is required.")
}

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

OUTDIR <- "outputs/reviewer_revisions/regional_formal_interaction"

DLPFC_MATRIX <- file.path(
  BASE, "Raw Data", "n1086_residual_log2_batch.csv"
)

STG_MATRIX <- file.path(
  BASE, "Raw Data", "n278_residual_log2_batch.TCX.csv"
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

EXISTING_BRAAK <- file.path(
  OUTDIR,
  "region_x_Braak_continuous_interaction_all_proteins.csv"
)


# =============================================================================
# Helpers
# =============================================================================

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "nan", "<NA>")] <- NA
  sub("\\.0$", "", x)
}


clean_missing <- function(x) {
  x <- trimws(as.character(x))

  x[x %in% c(
    "",
    "NA",
    "NaN",
    "nan",
    "<NA>",
    "missing or unknown"
  )] <- NA

  x
}


get_accession <- function(x) {
  x <- trimws(as.character(x))

  ifelse(
    grepl("\\|", x),
    sub("^.*\\|", "", x),
    x
  )
}


zscore <- function(x) {
  x <- as.numeric(x)

  s <- sd(x, na.rm = TRUE)

  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }

  as.numeric(scale(x))
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
    stop("Duplicate accessions in expression matrix.")
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

  by_person <- split(
    sm$matrix_col,
    sm$individual_id
  )

  ans <- vapply(
    by_person,
    function(cols) {

      y <- rowMeans(
        mat[, cols, drop = FALSE],
        na.rm = TRUE
      )

      y[is.nan(y)] <- NA_real_

      y
    },
    FUN.VALUE = numeric(nrow(mat))
  )

  rownames(ans) <- rownames(mat)

  ans
}


make_long <- function(
  dl,
  st,
  phenotype,
  predictor
) {

  d <- phenotype

  d$DLPFC <- dl[d$individual_id]
  d$STG <- st[d$individual_id]

  long <- rbind(
    transform(
      d,
      region = "DLPFC",
      abundance = DLPFC
    ),
    transform(
      d,
      region = "STG",
      abundance = STG
    )
  )

  long$region <- factor(
    long$region,
    levels = c("DLPFC", "STG")
  )

  required <- c(
    "abundance",
    predictor,
    "region",
    "cohort",
    "age_z",
    "sex",
    "pmi_z",
    "individual_id"
  )

  long <- long[
    complete.cases(
      long[, required, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]

  counts <- table(
    long$individual_id
  )

  complete_ids <- names(
    counts[counts == 2]
  )

  long <- long[
    long$individual_id %in% complete_ids,
    ,
    drop = FALSE
  ]

  long$individual_id <- factor(
    long$individual_id
  )

  long$cohort <- droplevels(
    factor(
      long$cohort,
      levels = c(
        "Mayo Clinic",
        "Emory"
      )
    )
  )

  long$sex <- droplevels(
    factor(long$sex)
  )

  long
}


fit_with_warning_capture <- function(
  long,
  predictor,
  optimizer = "optim",
  max_iter = 200,
  max_eval = 400
) {

  warnings_seen <- character()

  f <- as.formula(
    paste0(
      "abundance ~ region * ",
      predictor,
      " + region * cohort",
      " + region * age_z",
      " + region * sex",
      " + region * pmi_z"
    )
  )

  fit <- withCallingHandlers(

    tryCatch(
      nlme::lme(
        fixed = f,
        random = ~1 | individual_id,
        data = long,
        method = "REML",
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

      invokeRestart("muffleWarning")
    }
  )

  if (inherits(fit, "error")) {

    return(
      list(
        status = "error",
        warning = NA_character_,
        error = conditionMessage(fit),
        fit = NULL
      )
    )
  }

  status <- if (
    length(warnings_seen) > 0
  ) {
    "warning"
  } else {
    "ok"
  }

  list(
    status = status,
    warning = if (
      length(warnings_seen) > 0
    ) {
      paste(
        unique(warnings_seen),
        collapse = " | "
      )
    } else {
      NA_character_
    },
    error = NA_character_,
    fit = fit
  )
}


extract_interaction <- function(
  fit,
  predictor
) {

  tt <- summary(fit)$tTable

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

  term <- candidates[
    candidates %in% rownames(tt)
  ]

  if (length(term) != 1) {
    stop(
      "Could not uniquely identify interaction term."
    )
  }

  term <- term[1]

  c(
    beta = unname(
      tt[term, "Value"]
    ),
    se = unname(
      tt[term, "Std.Error"]
    ),
    p = unname(
      tt[term, "p-value"]
    )
  )
}


# =============================================================================
# Load metadata
# =============================================================================

cat("======================================================================\n")
cat("REGIONAL MODEL CONVERGENCE AUDIT\n")
cat("======================================================================\n")

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

pheno$braak_stage_num <- as.numeric(
  pheno$braak_stage_num
)

pheno$age_death_num <- as.numeric(
  pheno$age_death_num
)

pheno$pmi_num <- as.numeric(
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

pheno$age_z <- zscore(
  pheno$age_death_num
)

pheno$pmi_z <- zscore(
  pheno$pmi_num
)


# =============================================================================
# Load maps / identify paired participants
# =============================================================================

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

paired_ids <- intersect(
  unique(
    dl_map$individual_id[
      !is.na(dl_map$individual_id)
    ]
  ),
  unique(
    st_map$individual_id[
      !is.na(st_map$individual_id)
    ]
  )
)

if (length(paired_ids) != 215) {
  stop(
    "Expected 215 paired participants; found ",
    length(paired_ids)
  )
}


# =============================================================================
# Load and collapse matrices
# =============================================================================

cat("Loading and collapsing DLPFC/STG matrices...\n")

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

common <- sort(
  intersect(
    rownames(DL),
    rownames(ST)
  )
)

if (length(common) != 9015) {
  stop(
    "Expected 9015 shared proteins; found ",
    length(common)
  )
}

cat("Shared proteins:", length(common), "\n")


# =============================================================================
# Continuous Braak phenotype
# =============================================================================

braak <- pheno[
  pheno$individual_id %in% paired_ids &
    pheno$braak_stage_num %in% 1:6 &
    !is.na(pheno$age_z) &
    !is.na(pheno$sex) &
    !is.na(pheno$pmi_z) &
    !is.na(pheno$cohort),
  ,
  drop = FALSE
]

braak <- braak[
  order(braak$individual_id),
  ,
  drop = FALSE
]

if (nrow(braak) != 134) {
  stop(
    "Expected 134 Braak-eligible participants; found ",
    nrow(braak)
  )
}


# =============================================================================
# Rerun CONTINUOUS BRAAK ONLY with warning capture
# =============================================================================

existing <- read.csv(
  EXISTING_BRAAK,
  check.names = FALSE
)

audit <- vector(
  "list",
  length(common)
)

cat("\nAuditing continuous Braak convergence...\n")

for (i in seq_along(common)) {

  if (i %% 500 == 0) {
    cat(
      "  ",
      i,
      " / ",
      length(common),
      "\n",
      sep = ""
    )
  }

  acc <- common[i]

  long <- make_long(
    DL[acc, ],
    ST[acc, ],
    braak,
    "braak_stage_num"
  )

  n_people <- length(
    unique(long$individual_id)
  )

  if (n_people < 25) {

    audit[[i]] <- data.frame(
      accession = acc,
      n_people = n_people,
      status = "insufficient_n",
      warning = NA_character_,
      error = NA_character_,
      beta = NA_real_,
      se = NA_real_,
      p = NA_real_,
      stringsAsFactors = FALSE
    )

    next
  }

  ans <- fit_with_warning_capture(
    long,
    "braak_stage_num",
    optimizer = "optim"
  )

  est <- c(
    beta = NA_real_,
    se = NA_real_,
    p = NA_real_
  )

  if (!is.null(ans$fit)) {
    est <- extract_interaction(
      ans$fit,
      "braak_stage_num"
    )
  }

  audit[[i]] <- data.frame(
    accession = acc,
    n_people = n_people,
    status = ans$status,
    warning = ans$warning,
    error = ans$error,
    beta = unname(est["beta"]),
    se = unname(est["se"]),
    p = unname(est["p"]),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(
  rbind,
  audit
)

write.csv(
  audit,
  file.path(
    OUTDIR,
    "continuous_Braak_convergence_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Compare clean reruns with original output
# =============================================================================

comparison <- merge(
  audit[
    audit$status == "ok",
    c(
      "accession",
      "beta",
      "se",
      "p"
    )
  ],
  existing[
    ,
    c(
      "accession",
      "beta_region_x_predictor",
      "se_region_x_predictor",
      "p_region_x_predictor"
    )
  ],
  by = "accession",
  all = FALSE
)

comparison$beta_diff <- (
  comparison$beta -
  comparison$beta_region_x_predictor
)

comparison$se_diff <- (
  comparison$se -
  comparison$se_region_x_predictor
)

comparison$p_diff <- (
  comparison$p -
  comparison$p_region_x_predictor
)

write.csv(
  comparison,
  file.path(
    OUTDIR,
    "continuous_Braak_rerun_comparison.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Retry warning/error proteins using nlminb
# =============================================================================

flagged <- audit[
  audit$status %in% c(
    "warning",
    "error"
  ),
  ,
  drop = FALSE
]

retry_results <- list()

if (nrow(flagged) > 0) {

  cat(
    "\nRetrying flagged Braak models with nlminb...\n"
  )

  for (i in seq_len(nrow(flagged))) {

    acc <- flagged$accession[i]

    cat(
      "  ",
      acc,
      "\n",
      sep = ""
    )

    long <- make_long(
      DL[acc, ],
      ST[acc, ],
      braak,
      "braak_stage_num"
    )

    ans <- fit_with_warning_capture(
      long,
      "braak_stage_num",
      optimizer = "nlminb",
      max_iter = 1000,
      max_eval = 2000
    )

    est <- c(
      beta = NA_real_,
      se = NA_real_,
      p = NA_real_
    )

    if (!is.null(ans$fit)) {
      est <- extract_interaction(
        ans$fit,
        "braak_stage_num"
      )
    }

    retry_results[[i]] <- data.frame(
      accession = acc,
      original_status = flagged$status[i],
      original_warning = flagged$warning[i],
      retry_status = ans$status,
      retry_warning = ans$warning,
      retry_error = ans$error,
      retry_beta = unname(est["beta"]),
      retry_se = unname(est["se"]),
      retry_p = unname(est["p"]),
      stringsAsFactors = FALSE
    )
  }

  retry_results <- do.call(
    rbind,
    retry_results
  )

} else {

  retry_results <- data.frame(
    accession = character(),
    original_status = character(),
    original_warning = character(),
    retry_status = character(),
    retry_warning = character(),
    retry_error = character(),
    retry_beta = numeric(),
    retry_se = numeric(),
    retry_p = numeric()
  )
}

write.csv(
  retry_results,
  file.path(
    OUTDIR,
    "continuous_Braak_flagged_model_retries.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Q96B70 AD diagnostic
# =============================================================================

cat("\n======================================================================\n")
cat("Q96B70 / LENG9 AD FIT DIAGNOSTIC\n")
cat("======================================================================\n")

ad <- pheno[
  pheno$individual_id %in% paired_ids &
    pheno$ad_outcome %in% c(
      "Control",
      "AD"
    ) &
    !is.na(pheno$age_z) &
    !is.na(pheno$sex) &
    !is.na(pheno$pmi_z) &
    !is.na(pheno$cohort),
  ,
  drop = FALSE
]

ad$ad_status <- ifelse(
  ad$ad_outcome == "AD",
  1,
  0
)

q96 <- make_long(
  DL["Q96B70", ],
  ST["Q96B70", ],
  ad,
  "ad_status"
)

cat(
  "Participants:",
  length(unique(q96$individual_id)),
  "\n"
)

cat("\nAD status by participant:\n")
print(
  table(
    q96[
      !duplicated(q96$individual_id),
      "ad_status"
    ]
  )
)

cat("\nCohort x AD status:\n")
q96_people <- q96[
  !duplicated(q96$individual_id),
  ,
  drop = FALSE
]

print(
  table(
    q96_people$cohort,
    q96_people$ad_status
  )
)

cat("\nSex x AD status:\n")

print(
  table(
    q96_people$sex,
    q96_people$ad_status
  )
)

fixed_formula <- abundance ~
  region * ad_status +
  region * cohort +
  region * age_z +
  region * sex +
  region * pmi_z

mm <- model.matrix(
  fixed_formula,
  data = q96
)

cat(
  "\nFixed-effect design columns:",
  ncol(mm),
  "\n"
)

cat(
  "Fixed-effect design rank:",
  qr(mm)$rank,
  "\n"
)

cat(
  "Rank deficient:",
  qr(mm)$rank < ncol(mm),
  "\n"
)

if (qr(mm)$rank < ncol(mm)) {

  lm_alias <- alias(
    lm(
      fixed_formula,
      data = q96
    )
  )

  cat("\nAliased fixed-effect terms:\n")
  print(
    lm_alias$Complete
  )
}

cat(
  "\nRetrying Q96B70 AD with nlminb...\n"
)

q96_retry <- fit_with_warning_capture(
  q96,
  "ad_status",
  optimizer = "nlminb",
  max_iter = 1000,
  max_eval = 2000
)

cat(
  "Retry status:",
  q96_retry$status,
  "\n"
)

cat(
  "Retry warning:",
  q96_retry$warning,
  "\n"
)

cat(
  "Retry error:",
  q96_retry$error,
  "\n"
)

if (!is.null(q96_retry$fit)) {

  q96_est <- extract_interaction(
    q96_retry$fit,
    "ad_status"
  )

  cat("\nQ96B70 AD interaction after retry:\n")
  print(q96_est)
}


# =============================================================================
# Final summary
# =============================================================================

cat("\n======================================================================\n")
cat("CONVERGENCE AUDIT SUMMARY\n")
cat("======================================================================\n")

print(
  table(
    audit$status
  )
)

cat(
  "\nFlagged continuous-Braak models:",
  nrow(flagged),
  "\n"
)

if (nrow(flagged) > 0) {

  print(
    flagged[
      ,
      c(
        "accession",
        "n_people",
        "status",
        "warning",
        "error",
        "beta",
        "p"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )

  cat("\nFlagged-model retries:\n")

  print(
    retry_results,
    row.names = FALSE
  )
}

cat(
  "\nClean-rerun models compared:",
  nrow(comparison),
  "\n"
)

cat(
  "Max absolute beta difference vs original:",
  max(
    abs(comparison$beta_diff),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Max absolute SE difference vs original:",
  max(
    abs(comparison$se_diff),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Max absolute p-value difference vs original:",
  max(
    abs(comparison$p_diff),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nCONVERGENCE AUDIT COMPLETED.\n"
)
